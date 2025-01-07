import_metadata <- function(
    metadata_file,
    comparison_grouping_variable,
    sample_name_column,
    grouping_column,
    project_column,
    regression_columns,
    filter_column,
    filter_value                  = NULL,
    extra_columns                 = NULL,
    metadata_sheet                = "main",
    groups_to_include             = NULL,
    groups_to_exclude             = NULL,
    samples_to_exclude            = NULL,
    skip_lines                    = NULL,
    control_group                 = "control"
){
  `%||%` <- rlang::`%||%`

  # Setup variables for non-standard evaluation
  diffused_comparison_sym <- rlang::sym(comparison_grouping_variable)
  diffused_comparison_sym <- rlang::enquo(diffused_comparison_sym)

  diffused_sample_name <- rlang::sym(sample_name_column)
  diffused_sample_name <- rlang::enquo(diffused_sample_name)

  diffused_project <- rlang::sym(project_column)
  diffused_project <- rlang::enquo(diffused_project)

  study_metadata <-
    read_md_file(
      path = metadata_file,
      skip = skip_lines,
      sheet = metadata_sheet
    ) |>
    dplyr::select(
      tidyselect::all_of(
        c(
          sample_name_column,
          grouping_column,
          project_column,
          regression_columns,
          filter_column,
          extra_columns
        )
      )
    ) |>
    dplyr::filter(
      !is.na(sample_name_column),
      {{diffused_comparison_sym}} %in% (groups_to_include %||% unique(.data[[comparison_grouping_variable]])),
      !{{diffused_comparison_sym}} %in% groups_to_exclude,
      !{{diffused_sample_name}} %in% samples_to_exclude
    ) |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::where(is.character) & -tidyselect::all_of(sample_name_column),
        forcats::as_factor
      ),
      {{diffused_sample_name}} :=
        janitor::make_clean_names(
          string = {{diffused_sample_name}},
          case = "all_caps"
        ),
      {{diffused_comparison_sym}} :=
        tolower({{diffused_comparison_sym}}) |>
        forcats::fct_relevel(control_group)
    )
  dplyr::distinct(study_metadata)
}


import_counts <- function(
    directory,
    metadata
){

  tx_files <-
    dir(
      path = directory,
      pattern = "quant.sf.gz",
      # pattern = "fixed_quant.gz",
      recursive = TRUE,
      full.name = TRUE
    ) |>
    grep(
      pattern = "Undetermined|NONE",
      invert = TRUE,
      value = TRUE
    )

  tx_sample_names <-
    tx_files |>
    stringr::str_split(pattern = "/") |>
    purrr::map_chr(\(x) purrr::pluck(.x = x, length(x)-1)) |>
    stringr::str_remove(pattern = '(_[L|S][[:digit:]]+)+') |>
    janitor::make_clean_names(case = "all_caps")

  tx_files <-
    rlang::set_names(
      x = tx_files,
      nm = tx_sample_names
    )

  tx_files
}

prep_data_import <- function(
    count_files,
    sample_metadata,
    annotations,
    aligner            = "salmon",
    # minimum_gene_count = 1,
    removal_pattern    = "^RNA5",
    only_hugo          = TRUE
){
  if (tibble::is_tibble(sample_metadata)){
    sample_metadata <- tibble::olumn_to_rownames(sample_metadata, var = "sample_name")
  }

  common_names <-
    intersect(
      x = names(count_files),
      y = rownames(sample_metadata)
    )

  count_files <- magrittr::extract(count_files, common_names)
  filtered_metadata <- magrittr::extract(sample_metadata, common_names,)

  if (!is.data.frame(filtered_metadata)){
    filtered_metadata <- as.data.frame(filtered_metadata)
    rownames(filtered_metadata) <- common_names
    colnames(filtered_metadata) <- colnames(sample_metadata)
  }
  message("Importing count files")
  counts <-
    tximport::tximport(
      files    = count_files,
      type     = aligner,
      txIn     = TRUE,
      txOut    = FALSE,
      tx2gene  = annotations,
      importer = data.table::fread
    )

  filtered_counts <-
    purrr::map(
      .x = c("abundance", "counts", "length"),
      .f = \(x){
        cleaned_counts <-
          counts[[x]] |>
          tibble::as_tibble(rownames = "gene_symbol") |>
          dplyr::filter(
            # gene_symbol %in% genes_with_passing_counts,
            stringr::str_detect(
              string = gene_symbol,
              pattern = removal_pattern,
              negate = TRUE
            )
          )

        if (isTRUE(only_hugo)){
          cleaned_counts <-
            dplyr::mutate(
              hugo = HGNChelper::checkGeneSymbols(gene_symbol)[["Suggested.Symbol"]]
            ) |>
            dplyr::filter(
              .data = cleaned_counts,
              !is.na(hugo)
            ) |>
            dplyr::group_by(hugo) |>
            dplyr::mutate(across(where(is.numeric), sum)) |>
            dplyr::slice(1) |>
            dplyr::ungroup() |>
            dplyr::select(
              -gene_symbol,
              gene_symbol = hugo
            ) |>
            tibble::column_to_rownames("gene_symbol") |>
            as.matrix()
        } else {
          cleaned_counts <-
            cleaned_counts |>
          dplyr::mutate(across(-matches("gene_symbol"), as.numeric)) |>
            tibble::column_to_rownames("gene_symbol") |>
            as.matrix()
        }
      })
  filtered_counts[["countsFromAbundance"]] <- counts[["countsFromAbundance"]]
  filtered_counts <-
    rlang::set_names(
      x = filtered_counts,
      nm = names(counts)
    )

  list(
    metadata = filtered_metadata,
    counts   = filtered_counts
  )
}

create_final_md <- function(
    md,
    tx_files,
    study_design,
    comparison_group,
    control_group,
    sample_name
){

  if(is.character(comparison_group)){
    diffused_comparison_group <- rlang::sym(comparison_group)
    diffused_comparison_group <- rlang::enquo(diffused_comparison_group)
  } else {
    diffused_comparison_group <- rlang::enquo(comparison_group)
  }

  diffused_sample_name <- rlang::sym(sample_name)
  diffused_sample_name <- rlang::enquo(diffused_sample_name)

  final_md <-
    dplyr::filter(
      .data = md,
      {{diffused_sample_name}} %in% names(tx_files),
      stringr::str_detect(
        string = {{diffused_sample_name}},
        pattern = "_2$",
        negate = TRUE
      )
    ) |>
    dplyr::mutate(
      {{diffused_comparison_group}} :=
        forcats::as_factor({{diffused_comparison_group}}) |>
        forcats::fct_relevel(control_group),
      sample_name =
        janitor::make_clean_names(
          {{diffused_sample_name}},
          case="all_caps",
          allow_dupes = FALSE
        )
    ) |>
    dplyr::group_by(sample_name) |>
    # dplyr::slice_max(n = 1, order_by = initial_concentration_ng_ul, with_ties = FALSE) |>
    tibble::column_to_rownames(var="sample_name")

  final_md
}


#' @rdname remove_outliers
#' @method remove_outliers DGEList
#' @importFrom irlba prcomp_irlba
#' @importFrom limma voom
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate inner_join filter pull
#' @return list
remove_outliers <-
  function(
    object,
    design,
    pc1_zscore_cutoff,
    pc2_zscore_cutoff = NULL
  ){

    v <- limma::voom(object, design)
    pca_res =
      irlba::prcomp_irlba(v$E)[['rotation']] |>
      tibble::as_tibble() |>
      dplyr::mutate(
        sample_name = colnames(v$E),
        pc1_zscore = abs((PC1 - mean(PC1))/sd(PC1)),
        pc2_zscore = abs((PC2 - mean(PC2))/sd(PC2)),
      ) |>
      dplyr::inner_join(tibble::as_tibble(v$targets, rownames = "sample_name"))

    pc1_outliers <-
      dplyr::filter(
        .data = pca_res,
        pc1_zscore >= pc1_zscore_cutoff
      ) |>
      dplyr::pull(sample_name)

    if (!is.null(pc2_zscore_cutoff)){
      pc2_outliers <-
        dplyr::filter(
          .data = pca_res,
          pc2_zscore >= pc2_zscore_cutoff
        ) |>
        dplyr::pull(sample_name)
    } else {
      pc2_outliers <- NULL
    }

    outliers <- unique(c(pc1_outliers, pc2_outliers))

    if (length(outliers > 0)){
      object <- object[,colnames(object) %nin% outliers]
    }

    return(
      list(
        count_object = object,
        pca          = pca_res,
        removed      = outliers
      )
    )
  }


#' @title process_counts
#'
#' @description Read in RNAseq counts and perform QC,
#' normalization, modeling, and DEG analysis
#'
#' @param ... Parameters to pass along to the actual functions
#' that will process the RNAseq data
#' @param method Process the data using limma, edgeR, or DESeq2?
#'
#' @return
#' @export
#'
#' @examples

process_counts <-
  function(
    imported_counts,
    comparison_grouping_variable,
    batch_variable               = NULL,
    study_design                 = NULL,
    pc1_zscore_threshold         = 2,
    pc2_zscore_threshold         = 2,
    BPPARAM                      = BPPARAM,
    use_combat                   = FALSE,
    num_sva                      = 2,
    sva_control_genes            = NULL,
    minimum_gene_count           = 1,
    control_group                = "control",
    ...
  ){
    diffused_grouping_variable <- rlang::sym(comparison_grouping_variable)
    diffused_grouping_variable <- rlang::enquo(diffused_grouping_variable)

    if (is.null(study_design)){
      study_design = as.formula(paste("~", comparison_grouping_variable))
    }

    message("Correcting for effective library sizes")
    # Obtaining per-observation scaling factors for length, adjusted to avoid
    # changing the magnitude of the counts.
    normCts <- imported_counts[["counts"]][["counts"]]/imported_counts[["counts"]][["length"]]

    # Computing effective library sizes from scaled counts, to account for
    # composition biases between samples.
    eff.lib <- edgeR::calcNormFactors(normCts) * matrixStats::colSums2(normCts)

    # Multiply each gene by the library size for each sample and take the log
    # (sweep applies a function either rowwise or column wise; STATS is a vector
    # equal in length to the chosen axis to use as an argument to the applied function)
    normMat <-
      sweep(
        x      = imported_counts[["counts"]][["length"]]/exp(matrixStats::rowMeans2(log(imported_counts[["counts"]][["length"]]))),
        MARGIN = 2,
        STATS  = eff.lib,
        FUN    = "*"
      ) |>
      log()

    sample_grouping <-
      tidyr::unite(
        data = imported_counts[["metadata"]],
        col = "grouping",
        {{diffused_grouping_variable}}
      ) |>
      dplyr::pull(grouping)

    # Combining effective library sizes with the length factors, and calculating
    # offsets for a log-link GLM.
    message("Creating initial data object...")
    pre_qc_dge <-
      edgeR::DGEList(
        counts       = imported_counts[["counts"]][["counts"]],
        samples      = imported_counts[["metadata"]],
        group        = sample_grouping,
        lib.size     = eff.lib,
        remove.zeros = TRUE
      ) %>%
      edgeR::scaleOffset(
        # y = .,
        offset = normMat[rownames(.[["counts"]]),]
      ) %>%
      magrittr::extract(
        edgeR::filterByExpr(
          y               = .,
          group           = sample_grouping,
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )


    message("Creating preliminary study design...")
    preliminary_design <-
      model.matrix(
        object = study_design,
        data   = magrittr::extract(
          imported_counts[["metadata"]],
          colnames(pre_qc_dge),
        )
      ) %>%
      magrittr::set_colnames(
        c("Intercept", colnames(.)[2:ncol(.)])
      )

    message("Removing outliers...")
    outlier_qc <-
      remove_outliers(
        object            = pre_qc_dge,
        pc1_zscore_cutoff = pc1_zscore_threshold,
        pc2_zscore_cutoff = pc2_zscore_threshold,
        design            = preliminary_design
      )

    # Creating a DGEList object for use in edgeR.
    pre_sva_dge <-
      edgeR::scaleOffset(
        y                = outlier_qc[["count_object"]],
        offset           =
          normMat[
            rownames(outlier_qc[["count_object"]][["counts"]]),
            colnames(outlier_qc[["count_object"]][["counts"]])
          ]
      )

    pre_sva_dge <-
      magrittr::extract(
        pre_sva_dge,
        edgeR::filterByExpr(
          y               = pre_sva_dge,
          group           = sample_grouping,
          # keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    message("Running surrogate variable analysis...")
    sva_res <-
      calc_sva(
        object        = pre_sva_dge,
        model_design  = study_design,
        n.sva         = num_sva,
        control_genes = sva_control_genes
      )

    post_qc_dge <-
      edgeR::calcNormFactors(
        object = sva_res[["data_object"]]
      )

    # groups <- paste(c("0", comparison_grouping_variable), collapse = " + ")

    mm <-
      model.matrix(
        object = sva_res[["design"]],
        data = sva_res[["data_object"]][["samples"]]
      )

    message("Running voom...")
    voom_exprs <-
      limma::voomWithQualityWeights(
        counts    = post_qc_dge,
        design    = mm,
        plot      = TRUE,
        save.plot = TRUE
      )

    message("Fitting data...")
    fit <-
      limma::lmFit(
        object = voom_exprs,
        design = mm,
        method = "robust"
      )

    comparisons <-
      tibble::as_tibble(
        gtools::combinations(
          n = length(unique(imported_counts[["metadata"]][[comparison_grouping_variable]])),
          r = 2,
          v = unique(as.character(imported_counts[["metadata"]][[comparison_grouping_variable]])),
          repeats.allowed = FALSE
        )
      ) %>%
      rlang::set_names(
        nm = c("V1","V2")
      ) %>%
      dplyr::transmute(
        name = paste0(V2, "_vs_", V1),
        compare = paste(V2, "-", V1)
      ) %>%
      tibble::deframe()

    coeff <-
      stringr::str_remove(
        string = comparisons,
        pattern = "\\s-\\s[:graph:]+"
      ) %>%
      stringr::str_trim() %>%
      paste0(comparison_grouping_variable, .)

    contra_matrix <-
      limma::makeContrasts(
        contrasts = comparisons,
        levels =
          c(
            levels(post_qc_dge[["samples"]][[comparison_grouping_variable]]),
            str_extract_all(
              string  = as.character(sva_res[["design"]])[[2]],
              pattern = "SV[0-9]")[[1]]
          )
      )

    efit <- limma::eBayes(limma::contrasts.fit(fit, contra_matrix))

    res = purrr::map(
      .x = colnames(contra_matrix),
      .f = \(i) {
        message(stringr::str_glue("Performing DEG for {i}..."))
        limma::topTreat(
          fit = efit,
          coef = i,
          number = Inf
        ) %>%
          tibble::as_tibble(rownames = "gene") %>%
          dplyr::arrange(dplyr::desc(logFC)) %>%
          dplyr::rename(
            baseMean       = AveExpr,
            log2FoldChange = logFC,
            pvalue         = P.Value,
            padj           = adj.P.Val
          )
      }
    ) %>%
      rlang::set_names(
        nm = names(comparisons)
      )

    list(
      raw_counts                 = edgeR::getCounts(post_qc_dge),
      normalized_counts          = edgeR::cpm(post_qc_dge, log = TRUE, shrunk = TRUE),
      variance_stabilized_counts = voom_exprs[["E"]],
      outlier_samples            = outlier_qc[["removed"]],
      qc_pca                     = outlier_qc[["pca"]],
      sva_graph_data             = sva_res[["sva"]],
      design                     = sva_res[["design"]],
      dataset                    = post_qc_dge,
      comparisons                = comparisons,
      design_matrix              = mm,
      deg_results                = res,
      voom_exprs                 = voom_exprs,
      efit                       = efit,
      fit                        = fit
    )
  }

calc_sva <-
  function(
    object,
    model_design = NULL,
    batch_var = NULL,
    n.sva = 2,
    control_genes = NULL
  ){

    batch_var <- batch_var %||% 1

    svseq <-
      sva::svaseq(
        dat = object[["counts"]],
        mod = model.matrix(as.formula(model_design), object[["samples"]]),
        mod0 = model.matrix(as.formula(paste("~", 1)), object[["samples"]]),
        n.sv = n.sva
      )

    if (is.null(control_genes)){
      svseq <-
        sva::svaseq(
          dat = object[["counts"]],
          mod = model.matrix(as.formula(model_design), object[["samples"]]),
          mod0 = model.matrix(as.formula(paste("~", 1)), object[["samples"]]),
          n.sv = n.sva
        )
    } else {
      if (!all(control_genes %in% rownames(filtered_dat))){
        control_genes <- control_genes[which(control_genes %in% rownames(filtered_dat))]
        message(paste("Not all provided control genes are present in the data. Using those that were found:", paste(control_genes, collapse = ", ")))
      }
      control_prob <- as.integer(rownames(filtered_dat) %in% control_genes)
      svseq <-
        sva::svaseq(
          dat = object[["counts"]],
          mod = model.matrix(as.formula(model_design), object[["samples"]]),
          mod0 = model.matrix(as.formula(paste("~", 1)), object[["samples"]]),
          n.sv = n.sva,
          method = "supervised",
          controls = control_prob
        )
    }

    svseq[["sv"]] <-
      rlang::set_names(
        x = tibble::as_tibble(svseq[["sv"]], .name_repair = "unique"),
        nm = paste0("SV", seq(ncol(svseq[["sv"]])))
      ) |>
      dplyr::mutate(sample_name = rownames(object[["samples"]]))

    object[["samples"]] <-
      dplyr::left_join(
        tibble::as_tibble(object[["samples"]], rownames = "sample_name"),
        svseq[["sv"]]
      ) |>
      tibble::column_to_rownames(var = "sample_name")

    svseq$sv <- tibble::column_to_rownames(svseq[["sv"]],"sample_name")

    design_formula <-
      as.formula(
        paste("~",
              paste(as.character(model_design)[[2]],
                    paste(colnames(svseq[["sv"]]),
                          collapse = " + "),
                    sep = " + "),
              collapse = " "))

    ret_vals = list(
      data_object = object,
      sva = svseq,
      design = design_formula
    )
  }


plot_sva <- function(sva_graph_data){
  sva_graph_data |>
    purrr::pluck("sv") |>
    tibble::as_tibble(rownames = "sample_name") |>
    dplyr::select(
      sample_name,
      tidyselect::starts_with("SV")
    ) |>
    tidyr::pivot_longer(
      -sample_name,
      names_to = "covar"
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = sample_name,
        y = value
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "red"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(covar)
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          size = 9,
          hjust = 1,
          vjust = 1
        )
    )
}


read_md_file <- function(path, ...){
  if (!is.na(readxl::excel_format(path))){
    imported_file <-
      readxl::read_excel(
        path = path,
        na = c("n/a", "0", "unk", "NA", "NaN"),
        ...
      )
  } else {
    ext <- stringr::str_split(path, "\\.", simplify=TRUE)

    if (ext[[length(ext)]] == "csv"){
      delimiter = ","
    } else if (ext[[length(ext)]] == "tsv"){
      delimiter = "\t"
    }

    imported_file <-
      readr::read_delim(
        file = path,
        delim = delimiter,
        trim_ws = TRUE,
        na      = c("n/a", "0", "unk", "NA", "NaN")
      )
  }

  imported_file
}


extract_transformed_data <- function(data_obj){
  data_obj |>
    tibble::as_tibble(rownames = "gene") |>
    dplyr::mutate(hugo = HGNChelper::checkGeneSymbols(gene)[["Suggested.Symbol"]]) |>
    dplyr::filter(!is.na(hugo)) |>
    dplyr::group_by(hugo) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(
      -gene,
      gene = hugo
    ) |>
    tibble::column_to_rownames("gene") |>
    as.matrix()
}


filter_imported_data <- function(obj, ...){
  md <-
    obj$metadata |>
    dplyr::filter(...)

  counts <-
    purrr::map(
      .x = c("abundance", "counts", "length"),
      .f = \(x){
        obj[["counts"]][[x]] |>
          tibble::as_tibble(rownames = "gene_name") |>
          dplyr::select(gene_name, tidyselect::all_of(rownames(md))) |>
          tibble::column_to_rownames(var = "gene_name") |>
          as.matrix()
      }
    )
  names(counts) <- c("abundance", "counts", "length")
  counts[["countsFromAbundance"]] <- obj[["countsFromAbundance"]]
  return_obj <- list(
    counts = counts,
    metadata = md
  )
}

filter_processed_data <- function(
    previous_object,
    samples_to_keep,
    comparison_grouping_variable,
    ...
  ){

  samples_to_keep <-
    dplyr::intersect(
      x = samples_to_keep,
      y = colnames(previous_object[["variance_stabilized_counts"]])
      )

  filtered_dataset <-
    with(
      data = previous_object[["dataset"]],
      expr = edgeR::DGEList(
        counts = counts[,samples_to_keep],
        samples = samples[samples_to_keep,]
      )
    )
  filtered_dataset[["offset"]] <- previous_object[["offset"]][,samples_to_keep]
  filtered_dataset <- edgeR::calcNormFactors(filtered_dataset)

  filtered_sva_graph_data <-
    with(
      data = previous_object[["sva_graph_data"]],
      expr =
        list(
          sv = sv[samples_to_keep, ],
          pprob.gam = pprob.gam,
          pprob.b = pprob.b,
          n.sv = n.sv
        )
    )

  filtered_qc_pca <-
    previous_object[["qc_pca"]] %>%
    dplyr::filter(sample_name %in% samples_to_keep)

  filtered_outlier_samples <-
    previous_object[["outlier_samples"]][previous_object[["outlier_samples"]] %in% samples_to_keep]

  filtered_design_matrix <-
    previous_object[["design_matrix"]][samples_to_keep,]
  attributes(filtered_design_matrix)[["assign"]] <-
    attributes(previous_object[["design_matrix"]])[["assign"]]
  attributes(filtered_design_matrix)[["contrasts"]] <-
    attributes(previous_object[["design_matrix"]])[["contrasts"]]

  message("Re-running voom...")
  filtered_voom_exprs <-
    limma::voomWithQualityWeights(
      counts    = filtered_dataset,
      design    = filtered_design_matrix,
      plot      = FALSE,
      save.plot = FALSE
    )

  message("Re-fitting data...")
  filtered_fit <-
    limma::lmFit(
      object = filtered_voom_exprs,
      design = filtered_design_matrix,
      method = "robust"
    )

  coeff <-
    stringr::str_remove(
      string = previous_object[["comparisons"]],
      pattern = "\\s-\\s[:graph:]+"
    ) %>%
    stringr::str_trim() %>%
    paste0(comparison_grouping_variable, .)

  contra_matrix <-
    limma::makeContrasts(
      contrasts = previous_object[["comparisons"]],
      levels =
        c(
          levels(filtered_dataset[["samples"]][[comparison_grouping_variable]]),
          str_extract_all(
            string  = as.character(previous_object[["design"]])[[2]],
            pattern = "SV[0-9]")[[1]]
        )
    )

  filtered_efit <-
    limma::eBayes(
      fit =
        limma::contrasts.fit(
          fit = filtered_fit,
          contrasts = contra_matrix
          )
      )

  filtered_res = purrr::map(
    .x = colnames(contra_matrix),
    .f = \(i) {
      message(stringr::str_glue("Performing DEG for {i}..."))
      limma::topTreat(
        fit = filtered_efit,
        coef = i,
        number = Inf
      ) %>%
        tibble::as_tibble(rownames = "gene") %>%
        dplyr::arrange(dplyr::desc(logFC)) %>%
        dplyr::rename(
          baseMean       = AveExpr,
          log2FoldChange = logFC,
          pvalue         = P.Value,
          padj           = adj.P.Val
        )
    }
  ) %>%
    rlang::set_names(
      nm = names(previous_object[["comparisons"]])
    )

  list(
    raw_counts                 = edgeR::getCounts(filtered_dataset),
    normalized_counts          = edgeR::cpm(filtered_dataset, log = TRUE, shrunk = TRUE),
    variance_stabilized_counts = filtered_voom_exprs[["E"]],
    outlier_samples            = filtered_outlier_samples,
    qc_pca                     = filtered_qc_pca,
    sva_graph_data             = filtered_sva_graph_data,
    dataset                    = filtered_dataset,
    comparisons                = previous_object[["comparisons"]],
    design_matrix              = filtered_design_matrix,
    deg_results                = filtered_res,
    voom_exprs                 = filtered_voom_exprs,
    efit                       = filtered_efit,
    fit                        = filtered_fit
  )
}
