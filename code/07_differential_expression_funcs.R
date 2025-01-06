#' @title create_results_list
#' @rdname create_results_list
#' @method create_results_list limma
#' @importFrom magrittr use_series
#' @importFrom edgeR estimateDisp glmQLFit glmQLFTest getCounts
#' @importFrom rlang set_names
#' @importFrom purrr map map_chr chuck
#' @importFrom tibble as_tibble
#' @importFrom matrixStats rowMeans2
#' @importFrom rstatix adjust_pvalue
#' @importFrom dplyr mutate inner_join filter pull rename relocate
#' @return list
create_results_list <- function(
    processed_data,
    object,
    comparison_grouping_variable,
    BPPARAM       = bpparam(),
    ...
){


  if (!"deg_results" %in% names(processed_data)){
    comparison_list = processed_data[["comparisons"]]
    design          = processed_data[['design_matrix']]

    if (!"voom_exprs" %in% names(processed_data)){
      message("Running voom...")
      voom_exprs <-
        limma::voomWithQualityWeights(
          counts    = object,
          design    = design,
          plot      = FALSE,
          save.plot = FALSE
        )
    } else {
      message("Using previous voom results")
      voom_exprs <- processed_data[["voom_exprs"]]
    }

    if (!"efit" %in% names(processed_data)){
      if (!"fit" %in% names(processed_data)){
        message("Fitting data...")
        fit <-
          limma::lmFit(
            object = voom_exprs,
            design = design,
            method = "robust"
          )
      } else {
        fit <- processed_data[["fit"]]
      }
      message("Calculating empirical Bayes stats...")
      efit <- limma::eBayes(fit)
    } else {
      message("Using previous Bayes stats...")
      efit <- processed_data[["efit"]]
    }

    contra_matrix <-
      limma::makeContrasts(
        contrasts = comparison_list,
        levels =
          c(
            levels(object[["samples"]][[comparison_grouping_variable]]),
            colnames(dplyr::select(object[["samples"]], starts_with("SV")))
          )
      )

    res <- purrr::map(
      .x = colnames(contra_matrix),
      .f = \(i) {
        message(glue::glue("Performing DEG for {i}..."))
        limma::topTreat(
          fit = efit,
          coef = i,
          number = Inf
        ) |>
          tibble::as_tibble(rownames = "gene") |>
          dplyr::arrange(dplyr::desc(logFC)) |>
          dplyr::rename(
            baseMean       = AveExpr,
            log2FoldChange = logFC,
            pvalue         = P.Value,
            padj           = adj.P.Val
          )
      }
    ) |>
      rlang::set_names(
        purrr::map_chr(
          comparison_list,
          stringr::str_remove,
          pattern = paste0(comparison_grouping_variable, "_")
        )
      )
  } else {
    message("Using previously calculated results...")
    res <- processed_data[["deg_results"]]
  }

  res
}


create_deg_tables <- function(
  deg_res,
  comparison_list,
  grouping_variable,
  direction=c("up","down"),
  logfc_cutoff = 0,
  padj_cutoff = 0.05
){
  direction <- match.arg(arg = direction, choices=c("up", "down"))

  purrr::map(deg_res, \(i){
    degs <-
      if (!tibble::is_tibble(i)){
        tibble::as_tibble(
          x = i,
          rownames = "gene"
        )
      } else {
        degs <- i
      }

    degs <-
      switch(
        EXPR = direction,
        up   =
          dplyr::filter(
            .data = degs,
            !is.na(padj) & padj <= padj_cutoff,
            log2FoldChange > logfc_cutoff
            ),
        down =
          dplyr::filter(
            .data = degs,
            !is.na(padj) & padj <= padj_cutoff,
            log2FoldChange < logfc_cutoff
          )
      )

    dplyr::mutate(
      .data = degs,
      log2FoldChange = abs(log2FoldChange),
      dplyr::across(
        .cols = where(is.numeric),
        .fns = \(x) signif(x = x, digits =  4)
        )
      ) |>
    # dplyr::top_n(
    #   n = 25,
    #   wt = log2FoldChange
    #   ) |>
    dplyr::arrange(
      dplyr::desc(
        log2FoldChange
        )
      )
    }) |>
    rlang::set_names(
      nm = names(deg_res)
    ) |>
    purrr::keep(~ nrow(.x) > 0)
}


# create_upregulation_tables <- function(
    #   results,
#   comparison_list,
#   grouping_variable
# ){
#   map(seq_along(results), function(i){
#     results[[i]] |>
#     as_tibble(
#       rownames = "gene"
#     ) |>
#     filter(
#       !is.na(padj) & padj <= 0.05,
#       log2FoldChange > 0
#     ) |>
#     mutate_at(
#       vars(-gene),
#       list(~signif(., 2)
#       )
#     ) |>
#     top_n(
#       n = 25,
#       wt = log2FoldChange
#     ) |>
#     arrange(
#       desc(
#         log2FoldChange
#       )
#     )
# }) |>
#   set_names(
#     nm = map_chr(
#       .x = comparison_list,
#       .f = str_remove,
#       pattern = str_glue("{grouping_variable}_")
#     )
#   ) |>
#   keep(~ nrow(.x) > 0)
# }

extract_de_genes <- function(
    results,
    comparison_list,
    grouping_variable
){
  purrr::map(
    .x = results,
    .f = \(i){
      if (!tibble::is_tibble(i)){
        i <- tibble::as_tibble(rownames = "gene")
      }
      i |>
        dplyr::filter(padj < 0.05) |>
        dplyr::filter(abs(log2FoldChange) >= 0.5) |>
        dplyr::pull(gene)
      }) |>
    rlang::set_names(
      nm = comparison_list
    )
}


group_degs <- function(degs, comparison_vars){
  comparison_vars_regex = paste0(comparison_vars, "_", collapse="|")

  tibble::enframe(degs) |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::transmute(
      group =
        stringr::str_extract(
          string = name,
          pattern = comparison_vars_regex
          ) |>
        stringr::str_remove(pattern = "_$"),
      comparison =
        stringr::str_remove(
          string = name,
          pattern = comparison_vars_regex
        ),
      gene_symbol = value
      ) |>
    dplyr::group_by(gene_symbol) |>
    dplyr::summarise(
      gene_symbol = gene_symbol,
      count = n(),
      comparison =
        dplyr::case_when(
          count == 1 ~ comparison,
          count > 1 ~ "multiple"
        )
    ) |>
    dplyr::select(-count) |>
    dplyr::distinct() |>
    tibble::column_to_rownames(var = "gene_symbol")
}


calc_deg_means <- function(
    exprs,
    deg_class,
    metadata,
    grouping_variable
){
  if(is.character(grouping_variable)){
    diffused_grouping_variable <- sym(grouping_variable)
  }

  diffused_grouping_variable <- enquo(diffused_grouping_variable)

  deg_means <-
    exprs |>
    t() |>
    as_tibble(rownames = "name") |>
    select(name, one_of(deg_class[["gene_symbol"]])) |>
    pivot_longer(
      -name,
      names_to = "gene",
      values_to = "expr"
    ) |>
    left_join(
      metadata
    ) |>
    group_by(
      gene,
      {{diffused_grouping_variable}}
    ) |>
    summarise(avg = mean(expr)) |>
    pivot_wider(
      names_from = gene,
      values_from = avg
      ) |>
    column_to_rownames(grouping_variable)
}


extract_top_degs <- function(up_tables, down_tables){
  up_degs <-
    purrr::map(
      .x = names(up_tables),
      .f = function(i){
        up_tables[[i]] |>
          dplyr::filter(padj < 0.05) |>
          dplyr::top_n(25, log2FoldChange) |>
          dplyr::pull(gene)
        }
      ) |>
    rlang::set_names(names(up_tables)) |>
    unlist() # unlisting a map - should we just use `walk` here instead?

  down_degs <-
    purrr::map(
      .x = names(up_tables),
      .f = function(i){
        up_tables[[i]] |>
          dplyr::filter(padj < 0.05) |>
          dplyr::top_n(25, log2FoldChange) |>
          dplyr::pull(gene)
        }
      ) |>
    rlang::set_names(names(up_tables)) |>
    unlist()

  unique(c(up_degs, down_degs))
}


#' Title
#'
#' @param object DESeqResults object
#'
#' @return
#' @export
#'
#' @examples
alt_summary <-
  function(
    object,
    lfcThreshold = 0.25
  ){
    notallzero <- sum(object[["baseMean"]] > 0)
    up <- sum(
      object[["padj"]] < 0.05 &
        object[["log2FoldChange"]] > lfcThreshold,
      na.rm = TRUE
    )
    down <- sum(
      object[["padj"]] < 0.05 &
        object$log2FoldChange < lfcThreshold,
      na.rm = TRUE
    )
    outlier <-
      sum(
        object$baseMean > 0 &
          is.na(object$pvalue)
      )

    filterThresh <-
      filterThreshold(
        object = object,
        alfa = 0.1,
        padj_method = "fdr"
      )

    ft <-
      ifelse(
        test = is.null(filterThresh),
        yes  = 0,
        no   = round(filterThresh)
      )

    filt <- sum(!is.na(object[["pvalue"]]) & is.na(object[["padj"]]))

    total <- nrow(object)

    tibble::tibble(
      up = up,
      down = down,
      outlier = outlier,
      ft = ft,
      lowcounts = filt,
      total = total
    )
  }


#' filterThreshold
#'
#' @param object tibble from \code{create_results_list}
#' @param filter the vector of filter statistics over which the independent
#' filtering will be optimized. By default the mean of normalized counts is
#' used.
#' @param alfa the significance cutoff used for optimizing the independent
#' filtering (by default 0.1). If the adjusted p-value cutoff (FDR) will be a
#' value other than 0.1, alpha should be set to that value.
#' @param padj_method the method to use for adjusting p-values, see \code{?p.adjust}
#'
#' @return
#' @export
filterThreshold <- function(
    object,
    filter      = NULL,
    alfa        = 0.1,
    padj_method = "fdr"
){
  # Borrowed from {DESeq2} results.R `pvalueAdjustment()`
  filter        <- filter %||% object[["baseMean"]]

  lowerQuantile <- mean(filter == 0)
  upperQuantile <-
    ifelse(
      test = lowerQuantile < .95,
      yes  = 0.95,
      no   = 1
    )

  theta <-
    seq(
      lowerQuantile,
      upperQuantile,
      length = 50
    )


  # do filtering using genefilter
  filtPadj <- genefilter::filtered_p(
    filter = filter,
    test   = object[["pvalue"]],
    theta  = theta,
    method = "fdr"
  )

  numRej  <-
    matrixStats::colSums2(
      x = filtPadj < alfa,
      na.rm = TRUE
    )

  lo_fit <-
    stats::lowess(
      x = numRej ~ theta,
      f = 0.2
    )
  if (max(numRej) <= 10) {
    j <- 1
  } else {
    residual <- if (all(numRej==0)) {
      0
    } else {
      numRej[numRej > 0] - lo_fit[["y"]][numRej > 0]
    }

    thresh <- max(lo_fit[["y"]]) - sqrt(mean(residual^2))

    j <- if (any(numRej > thresh)) {
      which(numRej > thresh)[1]
    } else {
      1
    }
  }

  cutoffs <- stats::quantile(filter, theta)
  cutoffs[j]
}


deg_pathway_enrichment <-
  function(
    results,
    fcThreshold,
    fcPvalue,
    species = "human",
    direction = c("up", "down"),
    ontology  = c("ALL", "BP", "CC", "MF"),
    min_filtered_genes = 6
    ){
    direction <- match.arg(arg = direction, choices=c("up", "down"))
    ontology  <- match.arg(arg = ontology, choices=c("ALL", "BP", "CC", "MF"))

    filtered_results <-
      switch(
        EXPR = direction,
        up   = dplyr::filter(.data = results, padj <= fcPvalue, log2FoldChange >= fcThreshold),
        down = dplyr::filter(.data = results, padj <= fcPvalue, log2FoldChange >= -(fcThreshold))
      )

    # Technically, the enrich-series of functions work with gene symbols,
    # but using bitr ensures we eliminate those that do not map to a Entrez ID
    # and prevents a "Error in names(x) <- value :'names' attribute [2] must be the same length as the vector [1]"

    orgdb <- findOrgDb(target_species = species)

    filtered_gene_list <- dplyr::pull(.data = filtered_results, gene)
    if (length(filtered_gene_list) > min_filtered_genes){
      degIDs <-
        clusterProfiler::bitr(
          geneID   = filtered_gene_list,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = orgdb,
          drop     = TRUE
          ) |>
        dplyr::pull("ENTREZID")

      go_pathways <-
        clusterProfiler::enrichGO(
          gene          = degIDs,
          OrgDb         = orgdb,
          keyType       = "ENTREZID",
          ont           = ontology,
          pAdjustMethod = "fdr",
          readable      = TRUE
          )

      species <- switch(
        EXPR = species,
        `Homo sapiens` = "human",
        human          = "human",
        `H sapiens`    = "human",
        `Mus musculus` = "mouse",
        mouse          = "mouse",
        `M musculus`   = "mouse"
      )
      react_pathways <-
        ReactomePA::enrichPathway(
          gene          = degIDs,
          organism      = species,
          pAdjustMethod = "fdr",
          readable      = TRUE
        )
    } else {
      go_pathways <- NULL
      react_pathways <- NULL
    }

    list(
      gene_ontology = go_pathways,
      reactome      = react_pathways
    )
  }


get_enrichment_fcs <- function(enrichResult, degResult){

  if (!is.null(enrichResult[["reactome"]])){
    genes_in_reactome <- slot(enrichResult[["reactome"]], "gene2Symbol")
    reactome_lfc <-
      degResult |>
      dplyr::filter(gene %in% genes_in_reactome) |>
      dplyr::pull(log2FoldChange, gene)
  } else {
    reactome_lfc <- NULL
  }

  if (!is.null(enrichResult[["gene_ontology"]])){
    genes_in_go       <- slot(enrichResult[["gene_ontology"]], "gene2Symbol")
    go_lfc <-
      degResult |>
      dplyr::filter(gene %in% genes_in_go) |>
      dplyr::pull(log2FoldChange, gene)
  } else {
    go_lfc <- NULL
  }

  list(
    gene_ontology = go_lfc,
    reactome      = reactome_lfc
  )
}

