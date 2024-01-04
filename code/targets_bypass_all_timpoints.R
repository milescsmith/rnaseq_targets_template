source("code/project_parameters.R")

source("code/01_import_funcs.R")
source("code/02_filtering_funcs.R")
source("code/04_module_funcs.R")
source("code/05_cluster_funcs.R")
source("code/06_dimensional_reduction_funcs.R")
source("code/07_differential_expression_funcs.R")
source("code/08_WGCNA_funcs.R")
source("code/09_ml_funcs.R")
source("code/10_viral_transcript_funcs.R")
source("code/11_stats_testing_funcs.R")
# source("code/12_table_outputs.R")
source("code/13_pathways.R")
source("code/14_report_funcs.R")
source("code/97_misc_functions.R")
source("code/98_palettes_funcs.R")
source("code/99_output_funcs.R")

options(tidyverse.quiet = TRUE)
future::plan(strategy = future::multisession)
utils::globalVariables("where")
project_params[["comparison_grouping_variable"]] <- "responder_timepoint"
raw_metadata = "metadata/all_blast.csv"
md = import_metadata(
  metadata_file                 = raw_metadata,
  metadata_sheet                = project_params[["metadata_sheet"]],
  comparison_grouping_variable  = project_params[["comparison_grouping_variable"]],
  projects_to_include           = project_params[["projects_to_include"]],
  projects_to_exclude           = project_params[["projects_to_exclude"]],
  groups_to_include             = project_params[["groups_to_include"]],
  groups_to_exclude             = project_params[["groups_to_exclude"]],
  sample_name_column            = project_params[["sample_name_column"]],
  grouping_column               = project_params[["grouping_column"]],
  project_column                = project_params[["project_column"]],
  regression_columns            = project_params[["regression_columns"]],
  filter_column                 = project_params[["filter_column"]],
  filter_value                  = project_params[["filter_value"]],
  extra_columns                 = project_params[["extra_columns"]],
  samples_to_exclude            = project_params[["manual_sample_removal"]],
  extra_controls_metadata_file  = project_params[["extra_controls_sample_sheet"]],
  extra_controls_metadata_sheet = project_params[["main_sample_sheet"]],
  extra_controls_metadata_skip  = project_params[["main_sample_sheet_skip"]],
  extra_controls_ident_col      = project_params[["extra_controls_ident_col"]],
  extra_controls_ident          = "Control",
  skip_lines                    = project_params[["skip_lines"]],
  control_group                 = project_params[["control_group"]]
)

md <- mutate(
  .data = md,
  timepoint =
    case_when(
      responder == "control" ~ "control",
      str_detect(string = Sample_Alias, pattern = "BSL$") ~ "baseline",
      str_detect(string = Sample_Alias, pattern = "03",) ~ "month_3",
      TRUE ~ "remove"
    ),
    responder_timepoint =
      if_else(
        condition = responder == "control",
        true = "control",
        false = paste(responder, timepoint, sep="_")
      ) |>
      fct_relevel(
        "control",
        "non_responder_baseline",
        "non_responder_month_3",
        "responder_baseline",
        "responder_month_3"
      )
) |>
  filter(timepoint != "remove")

seq_file_directory = project_params[["sequencing_file_directory"]]

tx_files =
  import_counts(
    directory = seq_file_directory,
    metadata  = md
  )

annotation_file = project_params[["annotation_file"]]

annot = read_csv(annotation_file)

#### final_md ####
final_md    =
  create_final_md(
    md               = md,
    tx_files         = tx_files,
    comparison_group = project_params[["comparison_grouping_variable"]],
    control_group    = project_params[["control_group"]],
    sample_name      = project_params[["sample_name_column"]]
  ) |>
  mutate(across(where(is.factor), fct_drop))

#### imported_data ####
imported_data =
  prep_data_import(
    count_files        = tx_files,
    sample_metadata    = final_md,
    aligner            = project_params[["aligner"]],
    annotations        = annot,
    # minimum_gene_count = 1,
    removal_pattern    = "^RNA5",
    only_hugo          = project_params[["only_hugo_named_genes"]]
  )

#### processed_data ####
processed_data = process_counts(
  imported_counts              = imported_data,
  comparison_grouping_variable = "responder_timepoint",
  batch_variable               = project_params[["batch_variable"]],
  study_design                 = ~responder_timepoint,
  pc1_zscore_threshold         = project_params[["pc1_zscore_threshold"]],
  pc2_zscore_threshold         = project_params[["pc2_zscore_threshold"]],
  BPPARAM                      = BPPARAM,
  use_combat                   = project_params[["use_combat"]],
  minimum_gene_count           = project_params[["minimum_gene_count"]],
  sva_control_genes            = project_params[["sva_control_genes"]],
  num_sva                      = project_params[["sva_num"]],
  method                       = project_params[["process_method"]],
  control_group                = project_params[["control_group"]]
)

imported_counts              = imported_data
comparison_grouping_variable = project_params[["grouping_column"]]
batch_variable               = project_params[["batch_variable"]]
study_design                 = ~responder_timepoint
pc1_zscore_threshold         = project_params[["pc1_zscore_threshold"]]
pc2_zscore_threshold         = project_params[["pc2_zscore_threshold"]]
BPPARAM                      = BPPARAM
use_combat                   = project_params[["use_combat"]]
minimum_gene_count           = project_params[["minimum_gene_count"]]
sva_control_genes            = project_params[["sva_control_genes"]]
num_sva                      = project_params[["sva_num"]]
method                       = project_params[["process_method"]]
control_group                = project_params[["control_group"]]


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

    concentration_col <-
      grep(
        x = colnames(imported_counts[["metadata"]]),
        pattern="concentration",
        value = TRUE
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
          keep.lib.sizes  = FALSE,
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
        message(glue::glue("Performing DEG for {i}..."))
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

    processing_results <- list(
      raw_counts                 = edgeR::getCounts(post_qc_dge),
      normalized_counts          = edgeR::cpm(post_qc_dge, log = TRUE, shrunk = TRUE),
      variance_stabilized_counts = voom_exprs[["E"]],
      outlier_samples            = outlier_qc[["removed"]],
      qc_pca                     = outlier_qc[["pca"]],
      sva_graph_data             = sva_res[["sva"]],
      dataset                    = post_qc_dge,
      comparisons                = comparisons,
      design_matrix              = mm,
      deg_results                = res,
      voom_exprs                 = voom_exprs,
      efit                       = efit,
      fit                        = fit
    )

