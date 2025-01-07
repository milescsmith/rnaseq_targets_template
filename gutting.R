pre_qc_process <-
  function(
    imported_counts,
    comparison_grouping_variable,
    study_design                 = NULL,
    minimum_gene_count           = 1,
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

    list(
      dge      = pre_qc_dge,
      grouping = sample_grouping
    )
  }

create_pre_sva_object <-
  function(
    outlier_qc,
    sample_grouping,
    minimum_gene_count
  ){
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
            min.count       = minimum_gene_count,
            min.total.count = 10
          ),
        )
      pre_sva_dge
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



    message("Running surrogate variable analysis...")

    

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