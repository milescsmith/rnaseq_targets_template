# TODO: why are we processing baseline and controls twice?
# TODO: change to processing once, control for subject_id
# and then slice out the BL and controls.

library(targets)
library(tarchetypes)

source(here::here("code/project_parameters.R"))

source(here::here("code/01_import_funcs.R"))
source(here::here("code/02_filtering_funcs.R"))
source(here::here("code/04_module_funcs.R"))
source(here::here("code/05_cluster_funcs.R"))
source(here::here("code/06_dimensional_reduction_funcs.R"))
source(here::here("code/07_differential_expression_funcs.R"))
source(here::here("code/08_WGCNA_funcs.R"))
source(here::here("code/09_ml_funcs.R"))
source(here::here("code/10_viral_transcript_funcs.R"))
source(here::here("code/11_stats_testing_funcs.R"))
# source("code/12_table_outputs.R")
source(here::here("code/13_pathways.R"))
source(here::here("code/14_report_funcs.R"))
source(here::here("code/97_misc_functions.R"))
source(here::here("code/98_palettes_funcs.R"))
source(here::here("code/99_output_funcs.R"))
source(here::here("code/column_and_gene_lists.R"))
options(tidyverse.quiet = TRUE)
# future::plan(
#   strategy = future::multisession,
#   workers = parallelly::availableCores(which = "max")
# )
utils::globalVariables("where")
targets::tar_config_set(
  workers = parallelly::availableCores(which = "max")
)
# reticulate::use_condaenv(condaenv = "reticulate")
# targets::tar_option_set(debug = "hsic_down_genes_use_f60d376b", cue = tar_cue(mode = "never"))

list(

  tar_option_set(workspace_on_error = TRUE),
  targets::tar_target(
    name = sample_species_org,
    command = findOrgDb(project_params[["sample_species"]])
  ),

  targets::tar_target(
    name = samples_to_remove,
    command    = project_params[["manual_sample_removal"]],
    cue = targets::tar_cue(mode = "never"),
    deployment = "main"
  ),
  tarchetypes::tar_file_read(
    name = md,
    command = project_params[["metadata_file"]],
    read =
      import_metadata(
        metadata_file                 = !!.x,
        comparison_grouping_variable  = project_params[["comparison_grouping_variable"]],
        sample_name_column            = project_params[["sample_name_column"]],
        grouping_column               = project_params[["grouping_column"]],
        project_column                = project_params[["project_column"]],
        regression_columns            = project_params[["regression_columns"]],
        filter_column                 = project_params[["filter_column"]],
        filter_value                  = project_params[["filter_value"]],
        extra_columns                 = project_params[["extra_columns"]],
        metadata_sheet                = project_params[["metadata_sheet"]],
        groups_to_include             = project_params[["groups_to_include"]],
        groups_to_exclude             = project_params[["groups_to_exclude"]],
        samples_to_exclude            = project_params[["manual_sample_removal"]],
        skip_lines                    = project_params[["skip_lines"]],
        control_group                 = project_params[["control_group"]]
      ),
    packages =
      c(
        "janitor",
        "readxl",
        "dplyr",
        "tidyselect",
        "tidyr",
        "stringr"
      ),
    deployment = "main"
  ),

  tarchetypes::tar_file(
    name    = seq_file_directory,
    command = project_params[["sequencing_file_directory"]],
    cue     = targets::tar_cue(mode = "never"),
    deployment = "main"
  ),

  targets::tar_target(
    name = tx_files,
    command    =
      import_counts(
        directory = seq_file_directory,
        metadata  = md
      ),
    packages =
      c(
        "purrr",
        "magrittr",
        "stringr"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  #### annot ####
  tarchetypes::tar_file_read(
    name = annot,
    command = project_params[["annotation_file"]],
    read = readr::read_csv(!!.x),
    packages   = "readr"
  ),

  #### final_md ####
  targets::tar_target(
    name = final_md,
    command    =
      create_final_md(
        md               = md,
        tx_files         = tx_files,
        comparison_group = project_params[["comparison_grouping_variable"]],
        control_group    = project_params[["control_group"]],
        sample_name      = project_params[["sample_name_column"]]
      ),
    packages = c(
      "forcats",
      "stringr",
      "dplyr",
      "rlang",
      "tibble"
    ),
    cue = targets::tar_cue(mode = "never")
  ),

  #### imported_data ####
  targets::tar_target(
    name = imported_data,
    command =
      prep_data_import(
        count_files        = tx_files,
        sample_metadata    = final_md,
        aligner            = project_params[["aligner"]],
        annotations        = annot,
        # minimum_gene_count = 1,
        removal_pattern    = "^RNA5",
        only_hugo          = project_params[["only_hugo_named_genes"]]
      ),
    packages =
      c(
        "magrittr",
        "tximport",
        "data.table",
        "glue",
        "tibble",
        "dplyr",
        "purrr",
        "HGNChelper",
        "stringr"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = processed_data,
    command = process_counts(
      imported_counts              = imported_data,
      comparison_grouping_variable = project_params[["grouping_column"]],
      batch_variable               = project_params[["batch_variable"]],
      study_design                 = project_params[["study_design"]],
      pc1_zscore_threshold         = project_params[["pc1_zscore_threshold"]],
      pc2_zscore_threshold         = project_params[["pc2_zscore_threshold"]],
      BPPARAM                      = BPPARAM,
      use_combat                   = project_params[["use_combat"]],
      minimum_gene_count           = project_params[["minimum_gene_count"]],
      prune_majority_zero          = TRUE,
      sva_control_genes            = project_params[["sva_control_genes"]],
      num_sva                      = project_params[["sva_num"]],
      method                       = project_params[["process_method"]],
      control_group                = project_params[["control_group"]]
    ),
    packages =
      c(
        "edgeR",
        "DESeq2",
        "sva",
        "Rfast",
        "limma",
        "dplyr",
        "purrr",
        "tibble",
        "stringr",
        "rlang",
        "tidyr",
        "magrittr",
        "gtools"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = qc_pca,
    command = processed_data[["qc_pca"]]
  ),

  targets::tar_target(
    name = outlier_samples,
    command = processed_data[["outlier_samples"]]
  ),

  targets::tar_target(
    name = sva_graph_data,
    command = plot_sva(processed_data[["sva_graph_data"]]),
    packages =
      c(
        "ggplot2",
        "cowplot",
        "tidyselect",
        "tidyr",
        "tibble",
        "purrr"
      ),
    deployment = "main"
  ),

  targets::tar_target(
    name = pbmc_normalization_genes,
    command = c("ACTB", "HMBS", "HPRT1"),
    deployment = "main"
  ),
  targets::tar_target(
    name = vsc_exprs,
    command = processed_data[["variance_stabilized_counts"]],
    packages =
      c(
        "tibble",
        "dplyr",
        "HGNChelper"
      )
  ),

  # This should be changed into a list that we can walk through
  targets::tar_target(
    name = banchereau_module_file,
    command = project_params[["banchereau_modules"]],
    format = "file"
  ),

  targets::tar_target(
    name = banchereau_module_annotations_file,
    command = project_params[["banchereau_module_annotations"]],
    format = "file"
  ),

  targets::tar_target(
    name = banchereau_modules,
    command = create_module_list(banchereau_module_file),
    packages =
      c(
        "purrr",
        "tibble",
        "dplyr",
        "tidyselect",
        "rlang",
        "stringr"
      ),
    deployment = "main"
  ),

  targets::tar_target(
    name = responder_up_genes,
    command = up_tables |>
      with(responder_vs_non_responder) |>
      with(gene),
    deployment = "main"
  ),
  targets::tar_target(
    name = responder_down_genes,
    command = down_tables |>
      with(responder_vs_non_responder) |>
      with(gene),
    deployment = "main"
  ),

  #### top5_up ####
  targets::tar_target(
    name = top5_up_in_responders,
    command = res |>
      with(responder_vs_non_responder) |>
      dplyr::filter(
        padj < 0.05,
        gene %in% rownames(all_vsc_exprs)
      ) |>
      dplyr::slice_max(n = 5, order_by = log2FoldChange, with_ties = TRUE) |>
      # if there are ties, take from those the ones with the lowest padj
      dplyr::group_by(gene) |>
      dplyr::slice_min(n = 5, order_by = padj, with_ties = FALSE) |>
      dplyr::pull(gene),
    packages = "dplyr",
    deployment = "main"
  ),
  targets::tar_target(
    name = top5_down_in_responders,
    command = res |>
      with(responder_vs_non_responder) |>
      dplyr::filter(
        padj < 0.05,
        gene %in% rownames(all_vsc_exprs)
      ) |>
      dplyr::slice_min(n = 5, order_by = log2FoldChange, with_ties = TRUE) |>
      # if there are ties, take from those the ones with the lowest padj
      dplyr::group_by(gene) |>
      dplyr::slice_min(n = 5, order_by = padj, with_ties = FALSE) |>
      dplyr::pull(gene),
    packages = "dplyr",
    deployment = "main"
  ),
  targets::tar_target(
    name = deg_exprs,
    command =
      vsc_exprs |>
      tibble::as_tibble(rownames="gene") |>
      dplyr::filter(
        gene %in% c(responder_up_genes, responder_down_genes)
      ) |>
      tidyr::pivot_longer(
        -gene,
        names_to = "sample_name",
        values_to = "exprs"
      ) |>
      dplyr::left_join(
        dplyr::select(
          .data = final_md,
          sample_name,
          responder_status
        )
      ),
    packages = c(
      "tibble",
      "dplyr",
      "tidyr"
    ),
    deployment = "main"
  ),

  targets::tar_target(
    name = deg_exprs_stats,
    command =
      deg_exprs |>
      dplyr::group_by(gene) |>
      rstatix::wilcox_test(
        formula = exprs ~ responder_status,
        p.adjust.method = "fdr"
      ) |>
      rstatix::remove_ns() |>
      rstatix::add_xy_position(
        group = "responder_status",
        scales = "free_y"
      ) |>
      rstatix::add_significance(),
    packages = c(
      "dplyr",
      "rstatix"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = significant_degs,
    command =
      deg_exprs_stats |>
      dplyr::select(
        gene,
        group1,
        group2,
        p.adj.signif
      ) |>
      tidyr::pivot_wider(
        names_from = c("group1", "group2"),
        values_from = "p.adj.signif",
        names_glue = "{group1}_vs_{group2}",
        values_fn = forcats::as_factor
      ) |>
      dplyr::filter(
        if_any(
          .cols = tidyselect::where(is.factor),
          .fns = ~ .x != "ns"
        )
      ) |>
      dplyr::pull(gene),
    packages = c(
      "dplyr",
      "tidyr",
      "forcats",
      "tidyselect"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = deg_exprs_graph,
    command =
      deg_exprs |>
      dplyr::filter(gene %in% significant_degs) |>
      ggpubr::ggboxplot(
        x = "responder_status",
        fill = "responder_status",
        y = "exprs",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm() +
      ggpubr::stat_pvalue_manual(
        data = dplyr::filter(
          .data = deg_exprs_stats,
          gene %in% significant_degs
        )
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle=45,
            hjust=1,
            vjust=1
          )
      ) +
      ggforce::facet_wrap_paginate(
        ggplot2::vars(gene),
        ncol = 3,
        nrow = 3,
        scales = "free_y"
      ),
    packages = c(
      "dplyr",
      "ggpubr",
      "ggbeeswarm",
      "ggplot2",
      "ggforce"
    ),
    deployment = "main"
  ),

  #### top degs ####
  # perform PCA using only the top up- and down-regulated genes
  # only need the first PC

  # revised:
  # 1) need gene loadings so we can project all samples into metagene space
  # 2) if we calculate the score on all samples, but then subset on just the
  # baseline and control samples, we get a numerically stable score; if we
  # calculate on just baseline and controls, we get differences that are stable,
  # but the absolute values are not
  tarchetypes::tar_file_read(
    name = module_annotation,
    command = project_params[["banchereau_module_annotations"]],
    read =
      readr::read_csv(!!.x) |>
      dplyr::mutate(type = forcats::as_factor(type)),
    packages =
      c(
        "readr",
        "dplyr",
        "forcats",
        "magrittr"
      )
  ),

  tarchetypes::tar_file_read(
    name = ldg_modules,
    command = project_params[["ldg_modules"]],
    read = create_module_list(!!.x),
    packages =
      c(
        "purrr",
        "tibble",
        "tidyr",
        "readr",
        "dplyr"
      )
  ),

  tarchetypes::tar_file_read(
    name = metasignature_module,
    command = project_params[["metasignature_modules"]],
    read = create_module_list(!!.x),
    packages =
      c(
        "purrr",
        "tibble",
        "tidyr",
        "readr",
        "dplyr"
      )
  ),

  targets::tar_target(
    name = down_genes_use,
    command =
      scoreEigengenes(
        object = processed_data[["dataset"]],
        module_list = banchereau_modules,
        score_func = 'rsvd'
      ) |>
      scoreEigengenes(
        module_list = ldg_modules,
        score_func = 'rsvd'
      ) |>
      scoreEigengenes(
        module_list = metasignature_module,
        score_func = 'rsvd'
      ),
    packages =
      c(
        "moduleScoreR",
        "magrittr"
      )
  ),

  targets::tar_target(
    name = variable_down_deg_pca,
    command =
      calculate_module_score(
        expr_mat = all_vsc_exprs,
        gene_list = down_genes_use,
        score_func = project_params[["module_scoring_func"]]
      ),
    pattern = map(down_genes_use, down_number_of_genes),
    iteration = "list",
    packages = "magrittr",
    deployment = "worker"
  ),
  targets::tar_target(
    name = variable_down_deg_tbl,
    command =
      purrr::map2(
        .x = variable_down_deg_pca,
        .y = seq_along(variable_down_deg_pca),
        .f = \(x, idx) {
          with(x, scores) |>
            tibble::as_tibble(rownames = "sample_name") |>
            magrittr::set_colnames(c("sample_name", paste0("PC", idx + 1)))
        }
      ) |>
      purrr::reduce(left_join),
    packages = c(
      "purrr",
      "dplyr",
      "magrittr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = variable_up_deg_scores,
    command =
      variable_up_deg_tbl |>
      dplyr::filter(sample_name %in% bl_ctrl_md[["sample_name"]]) |>
      dplyr::left_join(bl_ctrl_md),
    packages =
      c(
        "dplyr",
        "purrr",
        "rlang"
      ),
    deployment = "main"
  ),

  targets::tar_target(
    name = annotated_mod_list,
    command = {
      dplyr::mutate(
        .data = annotated_modules,
        module_type =
          paste(
            module,
            type,
            sep = " - "
          )
      ) %>%
        dplyr::select(-type) %>%
        tibble::deframe()
    },
    packages =
      c(
        "dplyr",
        "tibble",
        "magrittr"
      )
  ),

  targets::tar_target(
    name = annotated_module_scores,
    command =
      dplyr::select(
        .data = module_scores,
        sample_name,
        tidyselect::all_of(annotated_modules[["module"]])
      ),
    packages =
      c(
        "dplyr",
        "tidyselect"
      )
  ),

  targets::tar_target(
    name = hsic_down_genes_use,
    command =
      hsic_lasso_select_features(
        expr_mat = vsc_exprs,
        metadata = dplyr::filter(bl_ctrl_md, responder_status != "control"), # improves selection
        gene_list = down_genes_use[[length(down_genes_use)]],
        sample_var = "sample_name",
        classification_var = "responder_status",
        num_features = down_number_of_genes
      ),
    pattern = map(down_number_of_genes),
    iteration = "list",
    packages = c(
      "reticulate",
      "rlang",
      "dplyr",
      "tibble",
      "tidyselect"
    ),
    deployment = "worker"
  ),

  targets::tar_target(
    name = variable_hsic_down_deg_pca,
    command =
      calculate_module_score(
        expr_mat = all_vsc_exprs,
        gene_list = hsic_down_genes_use,
        score_func = project_params[["module_scoring_func"]]
      ),
    pattern = map(hsic_down_genes_use, down_number_of_genes),
    iteration = "list",
    packages = "magrittr",
    deployment = "worker"
  ),
  targets::tar_target(
    name = variable_hsic_down_deg_tbl,
    command =
      purrr::map2(
        .x = variable_hsic_down_deg_pca,
        .y = seq_along(variable_hsic_down_deg_pca),
        .f = \(x, idx) {
          with(x, scores) |>
            tibble::as_tibble(rownames = "sample_name") |>
            magrittr::set_colnames(c("sample_name", paste0("PC", idx + 1)))
        }
      ) |>
      purrr::reduce(left_join),
    packages = c(
      "purrr",
      "dplyr",
      "magrittr"
    ),
    deployment = "main"
  ),

  targets::tar_target(
    name = variable_hsic_down_deg_scores,
    command =
      variable_hsic_down_deg_tbl |>
      dplyr::filter(sample_name %in% bl_ctrl_md[["sample_name"]]) |>
      dplyr::left_join(bl_ctrl_md),
    packages = c(
      "dplyr",
      "purrr",
      "rlang"
    ),
    deployment = "main"
  ),

  targets::tar_target(
    name = minimum_hsic_down_degs,
    command =
      variable_hsic_down_deg_scores |>
      dplyr::mutate(
        sample_name = colnames(vsc_exprs)
      ) |>
      dplyr::left_join(
        dplyr::select(
          .data = bl_ctrl_md,
          sample_name,
          responder_status
        )
      ) |>
      tidyr::pivot_longer(
        cols = starts_with("PC"),
        names_to = "number_degs",
        values_to = "module_score"
      ) |>
      dplyr::mutate(
        number_degs =
          stringr::str_remove(
            string = number_degs,
            pattern = "^PC"
          ) |>
          as.integer()
      ),
    packages = c(
      "dplyr",
      "tidyr",
      "stringr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_hsic_down_degs_stats,
    command =
      minimum_hsic_down_degs |>
      dplyr::group_by(number_degs) |>
      rstatix::wilcox_test(
        formula = module_score ~ responder_status,
        p.adjust.method = "fdr"
      ) |>
      rstatix::remove_ns() |>
      rstatix::add_xy_position(
        group = "responder_status",
        scales = "free_y"
      ) |>
      rstatix::add_significance(),
    packages = c(
      "dplyr",
      "rstatix"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_hsic_down_degs_graph,
    command =
      ggpubr::ggboxplot(
        data = minimum_hsic_down_degs,
        x = "responder_status",
        fill = "responder_status",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm() +
      ggpubr::stat_pvalue_manual(
        data = minimum_hsic_down_degs_stats
      ) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(number_degs),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(title = "Metagene scores using indicated number of most down-regulated in responders genes"),
    packages = c(
      "ggpubr",
      "ggbeeswarm",
      "targets",
      "ggplot2"
    ),
    deployment = "main"
  ),

  targets::tar_target(
    name = disease_class_tbl,
    command =
      final_md |>
      tibble::as_tibble(rownames = "sample_name") |>
      dplyr::select(sample_name, disease_class),
    deployment = "main",
    packages = c("tibble", "dplyr")
  ),

  targets::tar_target(
    name = hsic_up_genes_use,
    command =
      ident_clusters(
        expr_mat = vsc_exprs,
        max_k = 10,
        sig_pc_method = "horn"
      ),
    packages =
      c(
        "randomForest",
        "fpc",
        "PCAtools",
        "rlang",
        "tibble",
        "forcats",
        "dplyr",
        "magrittr"
      )
  ),

  targets::tar_target(
    name = study_md,
    command =
      dplyr::left_join(
        getMetaData(dataset_with_scores),
        sample_cluster_info[["clusters"]]
      ),
    packages = c(
      "dplyr",
      "purrr",
      "rlang"
    ),
    deployment = "main"
  ),

  targets::tar_target(
    name = minimum_hsic_up_degs,
    command =
      variable_hsic_up_deg_scores |>
      dplyr::mutate(
        sample_name = colnames(vsc_exprs)
      ) |>
      dplyr::left_join(
        dplyr::select(
          .data = bl_ctrl_md,
          sample_name,
          responder_status
        )
      ) |>
      tidyr::pivot_longer(
        cols = starts_with("PC"),
        names_to = "number_degs",
        values_to = "module_score"
      ) |>
      dplyr::mutate(
        number_degs =
          stringr::str_remove(
            string = number_degs,
            pattern = "^PC"
          ) |>
          as.integer()
      ),
    packages = c(
      "dplyr",
      "tidyr",
      "stringr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_hsic_up_degs_stats,
    command =
      minimum_hsic_up_degs |>
      dplyr::group_by(number_degs) |>
      rstatix::wilcox_test(
        formula = module_score ~ responder_status,
        p.adjust.method = "fdr"
      ) |>
      rstatix::remove_ns() |>
      rstatix::add_xy_position(
        group = "responder_status",
        scales = "free_y"
      ) |>
      rstatix::add_significance(),
    packages = c(
      "dplyr",
      "rstatix"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_hsic_up_degs_graph,
    command =
      ggpubr::ggboxplot(
        data = minimum_hsic_up_degs,
        x = "responder_status",
        fill = "responder_status",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm() +
      ggpubr::stat_pvalue_manual(
        data = minimum_hsic_up_degs_stats
      ) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(number_degs),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(title = "Metagene scores using indicated number of genes up-regulated in responders as selected by HSIC-Lasso"),
    packages = c(
      "ggpubr",
      "ggbeeswarm",
      "targets",
      "ggplot2"
    ),
    deployment = "main"
  ),

  #### HSIC combination of genes ####

  targets::tar_target(
    name = hsic_down_combos,
    command = gtools::combinations(
      n = length(hsic_down_genes_use[[length(hsic_down_genes_use)]]),
      r = project_params[["hsic_gene_set_length"]],
      v = hsic_down_genes_use[[length(hsic_down_genes_use)]]
    ),
    packages = "gtools"
  ),

  targets::tar_target(
    umap_results,
    run_umap(
      expr_data    = vsc_exprs,
      metadata     = study_md
    ),
    packages =
      c(
        "uwot",
        "tibble",
        "parallel",
        "rlang",
        "dplyr"
      )
  ),

  # TODO: {targets} should be able to handle mapping
  # the comparison list to the function instead of us mapping it

  #### res ####
  targets::tar_target(
    name = res,
    command =
      create_results_list(
        comparison_list              = processed_data[["comparisons"]],
        method                       = project_params[['process_method']],
        object                       = dataset_with_scores,
        comparison_grouping_variable = project_params[["comparison_grouping_variable"]],
        shrink_lfc                   = project_params[["shrink_lfc"]],
        BPPARAM                      = project_params[["BPPARAM"]],
        design                       = processed_data[['design_matrix']]
      ),
    packages =
      c(
        "purrr",
        "DESeq2",
        "tibble",
        "rlang",
        "magrittr",
        "edgeR",
        "matrixStats",
        "rstatix",
        "dplyr",
        "stringr"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  # TODO: generalize create_deg_tables?
  targets::tar_target(
    name = minimum_hsic_down_combo_sets,
    command = variable_hsic_down_combo_sets_scores |>
      dplyr::mutate(
        sample_name = colnames(vsc_exprs)
      ) |>
      dplyr::left_join(
        dplyr::select(
          .data = bl_ctrl_md,
          sample_name,
          responder_status
        )
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::contains("/"),
        names_to = "gene_set",
        values_to = "module_score"
      ),
    deployment = "main",
    packages = c(
      "dplyr",
      "tidyr",
      "tidyselect"
    )
  ),

  targets::tar_target(
    name = minimum_hsic_down_combo_sets_stats,
    command = minimum_hsic_down_combo_sets |>
      dplyr::group_by(gene_set) |>
      rstatix::wilcox_test(
        formula = module_score ~ responder_status,
        p.adjust.method = "fdr"
      ) |>
      rstatix::remove_ns() |>
      rstatix::add_xy_position(
        group = "responder_status",
        scales = "free_y"
      ) |>
      rstatix::add_significance(),
    deployment = "main",
    packages =
      c(
        "dplyr",
        "rstatix"
      )
  ),

  targets::tar_target(
    name = minimum_hsic_down_combo_sets_graph,
    command = ggpubr::ggboxplot(
      data = minimum_hsic_down_combo_sets,
      x = "responder_status",
      fill = "responder_status",
      y = "module_score",
      palette = "Set1"
    ) +
      ggbeeswarm::geom_beeswarm() +
      ggpubr::stat_pvalue_manual(
        data = minimum_hsic_down_combo_sets_stats
      ) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(gene_set),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(title = "Metagene scores using sets of HSIC-Lasso identified down-regulated in responders genes"),
    deployment = "main",
    packages =
      c(
        "ggpubr",
        "ggbeeswarm",
        "ggplot2",
        "ggforce"
      )
  ),

  targets::tar_target(
    name = up_enrichment,
    command =
      rlang::set_names(
        x = unnamed_up_enrichment,
        nm = names(res)
      ) |>
      purrr::map(
        .f = purrr::discard,
        .p = is.null
      ) |>
      purrr::map(
        .f = purrr::discard,
        .p = empty_enrichment
      )
  ),

  targets::tar_target(
    name = distinguishing_hsic_down_combo_sets_graph,
    command = minimum_hsic_down_combo_sets |>
      dplyr::filter(gene_set %in% distinguishing_hsic_down_gene_sets) |>
      ggpubr::ggboxplot(
        x = "responder_status",
        fill = "responder_status",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm() +
      ggpubr::stat_pvalue_manual(
        data = minimum_hsic_down_combo_sets_stats |> dplyr::filter(gene_set %in% distinguishing_hsic_down_gene_sets)
      ) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(gene_set),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(title = "Metagene scores using sets of HSIC-Lasso identified down-regulated genes that distinguish responders from non-responders"),
    deployment = "main",
    packages = c(
      "dplyr",
      "ggpubr",
      "ggbeeswarm",
      "ggplot2",
      "ggforce"
    )
  ),

  targets::tar_target(
    name = hsic_up_combos,
    command = gtools::combinations(
      n = length(hsic_up_genes_use[[length(hsic_up_genes_use)]]),
      r = project_params[["hsic_gene_set_length"]],
      v = hsic_up_genes_use[[length(hsic_up_genes_use)]]
    ),
    packages = "gtools"
  ),

  targets::tar_target(
    name = number_of_up_hsic_sets,
    command = seq(nrow(hsic_up_combos))
  ),

  targets::tar_target(
    name = hsic_up_gene_combo_sets,
    command = hsic_up_combos[number_of_up_hsic_sets,],
    pattern = map(number_of_up_hsic_sets),
    iteration = "list",
    deployment = "main",
  ),

  targets::tar_target(
    name = variable_hsic_up_combo_sets_pca,
    command =
      rlang::set_names(
        x = unnamed_down_enrichment,
        nm = names(res)
      ) |>
      purrr::map(
        .f = purrr::discard,
        .p = is.null
      ) |>
      purrr::map(
        .f = purrr::discard,
        .p = empty_enrichment
      )
  ),

  targets::tar_target(
    name = unnamed_down_enrichment_degs,
    command =
      if(!is.null(unlist(down_enrichment))){
        get_enrichment_fcs(
          enrichResult = down_enrichment[[res_length]],
          degResult = res[[res_length]]
        )
      } else {
        list()
      },
    pattern = map(res_length),
    iteration = "list",
    packages = "dplyr"
  ),

  targets::tar_target(
    name = down_enrichment_degs,
    command =
      if(
        length(unnamed_down_enrichment_degs[[1]][["gene_ontology"]] > 0) |
        length(unnamed_down_enrichment_degs[[1]][["reactome"]] > 0)
      ){
        rlang::set_names(x = unnamed_down_enrichment_degs, nm = names(res))
      } else {
        list()
      }
  ),

  targets::tar_target(
    name = minimum_hsic_up_combo_sets,
    command = variable_hsic_up_combo_sets_scores |>
      dplyr::mutate(
        sample_name = colnames(vsc_exprs)
      ) |>
      dplyr::left_join(
        dplyr::select(
          .data = bl_ctrl_md,
          sample_name,
          responder_status
        )
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::contains("/"),
        names_to = "gene_set",
        values_to = "module_score"
      ),
    deployment = "main",
    packages = c(
      "dplyr",
      "tidyr",
      "tidyselect"
    )
  ),

  targets::tar_target(
    name = minimum_hsic_up_combo_sets_stats,
    command = minimum_hsic_up_combo_sets |>
      dplyr::group_by(gene_set) |>
      rstatix::wilcox_test(
        formula = module_score ~ responder_status,
        p.adjust.method = "fdr"
      ) |>
      rstatix::remove_ns() |>
      rstatix::add_xy_position(
        group = "responder_status",
        scales = "free_y"
      ) |>
      rstatix::add_significance(),
    deployment = "main",
    packages =
      c(
        "dplyr",
        "rstatix"
      )
  ),

  targets::tar_target(
    name = minimum_hsic_up_combo_sets_graph,
    command = ggpubr::ggboxplot(
      data = minimum_hsic_up_combo_sets,
      x = "responder_status",
      fill = "responder_status",
      y = "module_score",
      palette = "Set1"
    ) +
      ggbeeswarm::geom_beeswarm() +
      ggpubr::stat_pvalue_manual(
        data = minimum_hsic_up_combo_sets_stats
      ) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(gene_set),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(
        title =
          paste(
            "Metagene scores for down-regulated,",
            "HSIC-Lasso-identified gene modules that",
            "distinguish responders from non-responders"
          )
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)),
    deployment = "main",
    packages =
      c(
        "ggpubr",
        "ggbeeswarm",
        "ggplot2",
        "ggforce"
      )
  ),

  targets::tar_target(
    name = distinguishing_hsic_up_gene_sets,
    command = minimum_hsic_up_combo_sets |>
      dplyr::group_by(
        gene_set,
        responder_status
      ) |>
      dplyr::summarise(
        min = min(module_score),
        max = max(module_score)
      ) |>
      tidyr::pivot_wider(
        names_from = responder_status,
        values_from = c(min, max)
      ) |>
      dplyr::mutate(
        overlap =
          DescTools::Overlap(
            c(min_non_responder, max_non_responder),
            c(min_responder, max_responder)
          )
      ) |>
      dplyr::ungroup() |>
      dplyr::filter(overlap <= 0) |>
      # dplyr::slice_min(
      #   order_by = overlap,
      #   with_ties = FALSE
      # ) |>
      dplyr::pull(gene_set),
    deployment = "main",
    packages = c(
      "dplyr",
      "tidyr",
      "DescTools"
    )
  ),

  targets::tar_target(
    name = distinguishing_hsic_up_combo_sets_graph,
    command = minimum_hsic_up_combo_sets |>
      dplyr::filter(gene_set %in% distinguishing_hsic_up_gene_sets) |>
      ggpubr::ggboxplot(
        x = "responder_status",
        fill = "responder_status",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm() +
      ggpubr::stat_pvalue_manual(
        data = minimum_hsic_up_combo_sets_stats |> dplyr::filter(gene_set %in% distinguishing_hsic_up_gene_sets)
      ) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(gene_set),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(title = "Metagene scores for up-regulated, HSIC-Lasso-identified gene modules that distinguish responders from non-responders"),
    deployment = "main",
    packages = c(
      "dplyr",
      "ggpubr",
      "ggbeeswarm",
      "ggplot2",
      "ggforce"
    )
  ),

  #### HSIC with baseline and month 3 ####

  targets::tar_target(
    name = minimum_hsic_down_degs_all_samples,
    command =
      variable_hsic_down_deg_tbl |>
      dplyr::left_join(
        dplyr::select(
          .data = all_passing_md,
          sample_name,
          subject_ref,
          visit,
          responder_status
        )
      ) |>
      tidyr::pivot_longer(
        cols = starts_with("PC"),
        names_to = "number_degs",
        values_to = "module_score"
      ) |>
      dplyr::mutate(
        number_degs =
          stringr::str_remove(
            string = number_degs,
            pattern = "^PC"
          ) |>
          as.integer(),
        responder_timepoint =
          dplyr::case_when(
            visit == "control" ~ "control",
            visit == "BL" ~ paste(responder_status, "baseline"),
            visit == "V3" ~ paste(responder_status, "month 3")
          )|>
          forcats::fct_relevel(
            "control",
            "non_responder baseline",
            "non_responder month 3",
            "responder baseline",
            "responder month 3"
          )
      ) |>
      dplyr::select(
        responder_timepoint,
        module_score,
        subject_ref,
        number_degs
      ),
    packages = c(
      "dplyr",
      "tidyr",
      "stringr",
      "forcats"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_hsic_down_degs_stats_all_samples,
    command =
      minimum_hsic_down_degs_all_samples |>
      dplyr::group_by(number_degs) |>
      rstatix::wilcox_test(
        formula = module_score ~ responder_timepoint,
        p.adjust.method = "fdr"
      ) |>
      rstatix::remove_ns() |>
      rstatix::add_xy_position(
        group = "responder_timepoint",
        scales = "free_y"
      ) |>
      rstatix::add_significance(),
    packages = c(
      "dplyr",
      "rstatix"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_hsic_down_degs_all_samples_graph,
    command =
      ggpubr::ggboxplot(
        data = minimum_hsic_down_degs_all_samples,
        x = "responder_timepoint",
        fill = "responder_timepoint",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::stat_pvalue_manual(
        data = minimum_hsic_down_degs_stats_all_samples
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(number_degs),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(title = "Metagene scores using indicated number of HSIC-Lasso-chosen genes in responders genes"),
    deployment = "main"
  ),

  targets::tar_target(
    name = sft,
    command = WGCNA::pickSoftThreshold(
      data = vsc_top,
      powerVector =
        c(
          seq(10),
          seq(
            from = 12,
            to = 30,
            by = 1
          )
        ),
      verbose = 5
    ),
    packages = "WGCNA"
  ),

  targets::tar_target(
    name = minimum_hsic_up_degs_stats_all_samples,
    command =
      minimum_hsic_up_degs_all_samples |>
      dplyr::group_by(number_degs) |>
      rstatix::wilcox_test(
        formula = module_score ~ responder_timepoint,
        p.adjust.method = "fdr"
      ) |>
      rstatix::remove_ns() |>
      rstatix::add_xy_position(
        group = "responder_timepoint",
        scales = "free_y"
      ) |>
      rstatix::add_significance(),
    packages = c(
      "dplyr",
      "rstatix"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_hsic_up_degs_all_samples_graph,
    command =
      ggpubr::ggboxplot(
        data = minimum_hsic_up_degs_all_samples,
        x = "responder_timepoint",
        fill = "responder_timepoint",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::stat_pvalue_manual(
        data = minimum_hsic_up_degs_stats_all_samples
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(number_degs),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(title = "Metagene scores using indicated number of HSIC-Lasso-chosen genes in responders genes"),
    deployment = "main"
  ),

  #### minimum number necessary genes####

  targets::tar_target(
    name = minimum_down_degs_no_overlap,
    command = 25, # Joan wants this set at 25, doesn't believe 10
    # minimum_down_degs |>
    # dplyr::group_by(
    #   number_degs,
    #   responder_status
    # ) |>
    # dplyr::summarise(
    #   min = min(module_score),
    #   max = max(module_score)
    # ) |>
    # tidyr::pivot_wider(
    #   names_from = responder_status,
    #   values_from = c(min, max)
    # ) |>
    # dplyr::mutate(
    #   overlap =
    #     DescTools::Overlap(
    #       c(min_non_responder, max_non_responder),
    #       c(min_responder, max_responder)
    #     )
    # ) |>
    # dplyr::ungroup() |>
    # dplyr::slice_min(
    #   order_by = overlap,
    #   with_ties = FALSE
    # ) |>
    # dplyr::pull(number_degs),
    packages = c(
      "dplyr",
      "DescTools"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = distinguishing_down_degs,
    command = down_genes_use[[which(purrr::map_int(.x = down_genes_use, .f = length) == minimum_down_degs_no_overlap)]],
    deployment = "main",
    packages = "purrr"
  ),
  targets::tar_target(
    name = distinguishing_up_degs,
    command = up_genes_use[[which(purrr::map_int(.x = up_genes_use, .f = length) == minimum_up_degs_no_overlap)]],
    deployment = "main",
    packages = "purrr"
  ),
  targets::tar_target(
    name = minimum_down_degs_significant,
    command =
      minimum_down_degs_stats |>
      dplyr::filter(
        group1 == "non_responder",
        group2 == "responder"
      ) |>
      dplyr::slice_min(
        order_by = number_degs,
        n = 1,
        with_ties = FALSE
      ) |>
      pull(number_degs),
    packages = "dplyr",
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_up_degs_significant,
    command =
      minimum_up_degs_stats |>
      dplyr::filter(
        group1 == "non_responder",
        group2 == "responder"
      ) |>
      dplyr::slice_min(
        order_by = number_degs,
        n = 1,
        with_ties = FALSE
      ) |>
      pull(number_degs),
    packages = "dplyr",
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_up_degs_no_overlap,
    command = 25, # Joan wants this set at 25, doesn't believe 10
    # minimum_up_degs |>
    # dplyr::group_by(
    #   number_degs,
    #   responder_status
    # ) |>
    # dplyr::summarise(
    #   min = min(module_score),
    #   max = max(module_score)
    # ) |>
    # tidyr::pivot_wider(
    #   names_from = responder_status,
    #   values_from = c(min, max)
    # ) |>
    # dplyr::mutate(
    #   overlap =
    #     DescTools::Overlap(
    #       c(min_non_responder, max_non_responder),
    #       c(min_responder, max_responder)
    #     )
    # ) |>
    # dplyr::ungroup() |>
    # dplyr::slice_min(
    #   order_by = overlap,
    #   with_ties = FALSE
    # ) |>
    # dplyr::pull(number_degs),
    packages = c(
      "dplyr",
      "DescTools"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_down_degs_no_overlap_table,
    command =
      down_genes_use[[minimum_down_degs_no_overlap - 1]] |>
      knitr::kable(
        caption = "Maximum expression levels in controls of responder down-regulated genes"
      ) |>
      kableExtra::collapse_rows(
        columns = 1,
        valign="middle"
      ) |>
      kableExtra::kable_styling(
        bootstrap_options =
          c(
            "striped",
            "hover",
            "condensed",
            "responsive"
          ),
        latex_options =
          c(
            c(
              "striped",
              "hold_position",
              "repeat_header"
            )
          ),
        position = "center",
        full_width = TRUE,
        font_size = 12
      ),
    packages = c(
      "knitr",
      "kableExtra"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = down_degs_table,
    command =
      responder_down_genes |>
      tibble::tibble() |>
      dplyr::rename(`Down-regulated genes` = "responder_down_genes") |>
      knitr::kable() |>
      kableExtra::collapse_rows(
        columns = 1,
        valign="middle"
      ) |>
      kableExtra::kable_styling(
        bootstrap_options =
          c(
            "striped",
            "hover",
            "condensed",
            "responsive"
          ),
        latex_options =
          c(
            c(
              "striped",
              "hold_position",
              "repeat_header"
            )
          ),
        position = "center",
        full_width = TRUE,
        font_size = 12
      ),
    packages = c(
      "knitr",
      "kableExtra"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = responder_down_genes_functions,
    command =
      biomaRt::getBM(
        attributes = c(
          "hgnc_symbol",
          "entrezgene_description",
          "goslim_goa_description"
        ),
        filters = "hgnc_symbol",
        values = responder_down_genes,
        mart = biomaRt::useMart(
          dataset = "hsapiens_gene_ensembl",
          biomart = "ENSEMBL_MART_ENSEMBL",
          host = "https://useast.ensembl.org"
        )
      ) |>
      dplyr::group_by(hgnc_symbol) |>
      tidyr::nest(
        goslim_goa_description = goslim_goa_description
      ) |>
      dplyr::mutate(
        goslim_goa_description = purrr::map_chr(
          .x = goslim_goa_description,
          .f = \(x) stringr::str_c(
            x[[1]],
            collapse = ", "
          )
        )
      ) |>
      dplyr::group_by(hgnc_symbol) |>
      dplyr::slice_head(n = 1),
    packages = c("biomaRt", "tidyr","dplyr", "stringr","purrr")
  ),
  targets::tar_target(
    name = unnamed_wgcna_rf_models,
    command   =
      rf_classifier(
        dataset            = wgcna_scores_with_md,
        classification_var = comparison_groups,
        starts_with("ME"),
        train_proportion   = 0.75,
        print              = FALSE
      ),
    pattern   = map(comparison_groups),
    iteration = "list",
    packages  =
      c(
        "rlang",
        "ranger",
        "randomForest",
        "dplyr",
        "forcats",
        "rsample",
        "recipes",
        "parsnip",
        "dials",
        "workflows",
        "tune"
      ),
    cue = targets::tar_cue("never")
  ),

  targets::tar_target(
    name = wgcna_rf_models,
    command  =
      rlang::set_names(
        x    = unnamed_wgcna_rf_models,
        nm   = comparison_groups
      ),
    packages = "rlang"
  ),

  targets::tar_target(
    name = responder_up_genes_functions,
    command =
      biomaRt::getBM(
        attributes = c(
          "hgnc_symbol",
          "entrezgene_description",
          "goslim_goa_description"
        ),
        filters = "hgnc_symbol",
        values = responder_up_genes,
        mart = biomaRt::useMart(
          dataset = "hsapiens_gene_ensembl",
          biomart = "ENSEMBL_MART_ENSEMBL",
          host = "https://useast.ensembl.org"
        )
      ) |>
      dplyr::group_by(hgnc_symbol) |>
      tidyr::nest(
        goslim_goa_description = goslim_goa_description
      ) |>
      dplyr::mutate(
        goslim_goa_description = purrr::map_chr(
          .x = goslim_goa_description,
          .f = \(x) stringr::str_c(
            x[[1]],
            collapse = ", "
          )
        )
      ) |>
      dplyr::group_by(hgnc_symbol) |>
      dplyr::slice_head(n = 1),
    packages = c("biomaRt", "tidyr","dplyr", "stringr","purrr")
  ),
  targets::tar_target(
    name = only_hugo_approved_genes,
    command =
      dplyr::left_join(
        with(res, responder_vs_non_responder),
        magrittr::set_colnames(
          HGNChelper::checkGeneSymbols(
            x = with(res, responder_vs_non_responder) |> with(gene),
            unmapped.as.na = TRUE
          ),
          c("gene", "approved", "suggested_symbol")
        )
      ) |>
      dplyr::filter(approved == TRUE) |>
      dplyr::pull(gene),
    packages = c("dplyr", "magrittr", "HGNChelper"),
    deployment = "main"
  ),
  targets::tar_target(
    name = top5_up,
    command = res |>
      with(responder_vs_non_responder) |>
      #filter(gene %in% only_hugo_approved_genes) |>
      dplyr::filter(gene %in% rownames(all_vsc_exprs)) |>
      dplyr::slice_max(
        n = 5,
        order_by = padj,
        with_ties = FALSE
      ) |>
      dplyr::pull(gene),
    packages = "dplyr"
  ),
  targets::tar_target(
    name = top5_down,
    command = res |>
      with(responder_vs_non_responder) |>
      #filter(gene %in% only_hugo_approved_genes) |>
      dplyr::filter(gene %in% rownames(all_vsc_exprs)) |>
      dplyr::slice_min(
        n = 5,
        order_by = padj,
        with_ties = FALSE
      ) |>
      dplyr::pull(gene),
    packages = "dplyr",
    deployment = "main"
  ),
  targets::tar_target(
    name = top5_up_exprs,
    command = all_vsc_exprs[top5_up, ],
    deployment = "main"
  ),
  targets::tar_target(
    name = top5_down_exprs,
    command = all_vsc_exprs[top5_down, ],
    deployment = "main"
  ),
  targets::tar_target(
    name = top5_up_exprs_pivot,
    command =
      top5_up_exprs |>
      tibble::as_tibble(rownames = "gene") |>
      tidyr::pivot_longer(
        -gene,
        names_to = "sample_name",
        values_to = "baseline_exprs"
      ) |>
      dplyr::left_join(
        dplyr::select(
          final_md,
          sample_name,
          responder_status,
          subject_ref
        )
      ),
    packages = c(
      "tibble",
      "tidyr",
      "dplyr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = top5_down_exprs_pivot,
    command =
      top5_down_exprs |>
      tibble::as_tibble(rownames = "gene") |>
      tidyr::pivot_longer(
        -gene,
        names_to = "sample_name",
        values_to = "baseline_exprs"
      ) |>
      dplyr::left_join(
        x   = module_scores,
        y   = tibble::as_tibble(
          x        = annotation_info,
          rownames = "sample_name"
        )
      ),
    packages =
      c(
        "dplyr",
        "tibble"
      )
  ),

  targets::tar_target(
    name = unnamed_module_rf_models,
    command   =
      {
        message(comparison_groups)
        rf_classifier(
          dataset            = md_with_module_scores,
          classification_var = comparison_groups,
          tidyselect::matches("^M[[:digit:]]+"),
          train_proportion   = 0.75,
          print              = FALSE
        )
      },
    pattern   = map(comparison_groups),
    iteration = "list",
    packages  =
      c(
        "rlang",
        "ranger",
        "randomForest",
        "dplyr",
        "forcats",
        "rsample",
        "recipes",
        "parsnip",
        "dials",
        "workflows",
        "tune"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = max_top5_down_values,
    command =
      top5_down_exprs_pivot |>
      dplyr::filter(responder_status == "control") |>
      dplyr::group_by(gene) |>
      dplyr::summarise(
        max_expr = max(baseline_exprs),
        mean_expr = mean(baseline_exprs),
        std_expr = sd(baseline_exprs)
      ) |>
      dplyr::select(
        gene,
        max_expr
      ),
    packages = "dplyr",
    deployment = "main"
  ),

  targets::tar_target(
    name = control_down_exprs_pivot,
    command =
      all_vsc_exprs |>
      tibble::as_tibble(rownames = "gene") |>
      dplyr::select(
        gene,
        tidyselect::all_of(all_passing_control_samples)
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(all_passing_control_samples),
        names_to = "sample_name",
        values_to = "baseline_exprs"
      ) |>
      dplyr::left_join(
        dplyr::select(
          .data = all_passing_md,
          sample_name,
          subject_ref,
          responder_status
        )
      ) |>
      dplyr::mutate(
        month3_exprs = baseline_exprs,
        timepoint = "control"
      ),
    packages = c(
      "tibble",
      "tidyr",
      "tidyselect",
      "dplyr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = baseline_down_exprs_pivot,
    command =
      all_vsc_exprs |>
      tibble::as_tibble(rownames = "gene") |>
      dplyr::select(
        gene,
        tidyselect::all_of(all_passing_baseline_samples)
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(all_passing_baseline_samples),
        names_to = "sample_name",
        values_to = "baseline_exprs"
      ) |>
      dplyr::left_join(
        dplyr::select(
          .data = all_passing_md,
          sample_name,
          subject_ref,
          responder_status
        )
      ) |>
      dplyr::left_join(max_top5_down_values) |>
      dplyr::mutate(
        baseline_diff = baseline_exprs - max_expr,
        timepoint = "baseline"
      ),
    packages = c(
      "tibble",
      "tidyr",
      "tidyselect",
      "dplyr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = month3_down_exprs_pivot,
    command =
      all_vsc_exprs |>
      tibble::as_tibble(rownames = "gene") |>
      dplyr::select(
        gene,
        tidyselect::all_of(all_passing_month3_samples)
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(all_passing_month3_samples),
        names_to = "sample_name",
        values_to = "month3_exprs"
      ) |>
      dplyr::left_join(
        dplyr::select(
          .data = all_passing_md,
          sample_name,
          subject_ref,
          responder_status
        )
      ) |>
      dplyr::left_join(max_top5_down_values) |>
      dplyr::mutate(
        month3_diff = month3_exprs - max_expr,
        timepoint = "month 3"
      ),
    packages = c(
      "tibble",
      "tidyr",
      "tidyselect",
      "dplyr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = baseline_vs_month3_top5_down_graph,
    command =
      dplyr::inner_join(
        x = baseline_down_exprs_pivot,
        y = month3_down_exprs_pivot,
        by = c("subject_ref", "gene", "responder_status", "max_expr"),
        keep = FALSE
      ) |>
      tidyr::pivot_longer(
        cols = c(baseline_diff, month3_diff),
        names_pattern = "(\\w+)_",
        names_to = "timepoint",
        values_to = "difference"
      ) |>
      dplyr::mutate(
        responder_timepoint =
          paste(
            responder_status,
            timepoint
          ) |>
          forcats::fct_relevel(
            "non_responder baseline",
            "non_responder month3",
            "responder baseline",
            "responder month3"
          )
      ) |>
      dplyr::filter(gene %in% top5_down) |>
      ggpubr::ggviolin(
        x = "responder_timepoint",
        y = "difference",
        fill = "responder_timepoint",
        palette = "Paired"
      ) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(2,1,6,5)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::facet_wrap(
        facets = ggplot2::vars(gene),
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = "Difference in expression vs control upper limit",
        subtitle = "Down-regulated in responders vs non-responders, baseline vs month 3"
      ),
    packages = c(
      "dplyr",
      "tidyr",
      "forcats",
      "ggpubr",
      "ggplot2",
      "RColorBrewer",
      "ggbeeswarm"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = max_top5_up_values,
    command =
      top5_up_exprs_pivot |>
      dplyr::filter(responder_status == "control") |>
      dplyr::group_by(gene) |>
      dplyr::summarise(
        max_expr = max(baseline_exprs),
        mean_expr = mean(baseline_exprs),
        std_expr = sd(baseline_exprs)
      ) |>
      dplyr::select(
        gene,
        max_expr
      ),
    packages = "dplyr",
    deployment = "main"
  ),

  targets::tar_target(
    name = control_up_exprs_pivot,
    command =
      all_vsc_exprs |>
      tibble::as_tibble(rownames = "gene") |>
      dplyr::select(
        gene,
        tidyselect::all_of(all_passing_control_samples)
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(all_passing_control_samples),
        names_to = "sample_name",
        values_to = "baseline_exprs"
      ) |>
      dplyr::left_join(
        study_md,
        all_module_scores
      ),
    packages = c(
      "tibble",
      "tidyr",
      "tidyselect",
      "dplyr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = baseline_up_exprs_pivot,
    command =
      rlang::set_names(
        nm =
          targets_recode(
            target_list = names(all_module_scores_with_md),
            thing_to_unquote_splice = annotated_mod_list
          ),
        x  = all_module_scores_with_md
      ),
    packages = c(
      "tibble",
      "tidyr",
      "tidyselect",
      "dplyr"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = baseline_vs_month3_top5_up_graph,
    command =
      dplyr::inner_join(
        x = baseline_up_exprs_pivot,
        y = month3_up_exprs_pivot,
        by = c("subject_ref", "gene", "responder_status", "max_expr"),
        keep = FALSE
      ) |>
      tidyr::pivot_longer(
        cols = c(baseline_diff, month3_diff),
        names_pattern = "(\\w+)_",
        names_to = "timepoint",
        values_to = "difference"
      ) |>
      dplyr::mutate(
        responder_timepoint =
          paste(
            responder_status,
            timepoint
          ) |>
          forcats::fct_relevel(
            "non_responder baseline",
            "non_responder month3",
            "responder baseline",
            "responder month3"
          )
      ) |>
      dplyr::filter(gene %in% top5_up) |>
      ggpubr::ggviolin(
        x = "responder_timepoint",
        y = "difference",
        fill = "responder_timepoint",
        palette = "Paired"
      ) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(2,1,6,5)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::facet_wrap(
        facets = ggplot2::vars(gene),
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = "Difference in expression vs control upper limit",
        subtitle = "Up-regulated in responders vs non-responders, baseline vs month 3"
      ),
    packages = c(
      "dplyr",
      "tidyr",
      "forcats",
      "ggpubr",
      "ggplot2",
      "RColorBrewer",
      "ggbeeswarm"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = baseline_vs_month3_with_ctrls_down_graph,
    command =
      module_scores_mat <-
      all_module_scores |>
      dplyr::select(
        sample_name,
        tidyselect::starts_with("ME"),
        tidyselect::one_of(
          c(
            names(ldg_modules),
            names(banchereau_modules),
            names(metasignature_module)
          )
        )
      ) |>
      ggpubr::ggviolin(
        x = "responder_timepoint",
        y = "exprs",
        fill = "responder_timepoint",
        palette = "Paired"
      ) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::labs(
        title = "Down-regulated in responders vs non-responders, baseline vs month 3"
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(gene),
        scales = "free_y",
        nrow = 3,
        ncol = 3
      ),
    packages = c(
      "dplyr",
      "tidyr",
      "forcats",
      "ggpubr",
      "ggplot2",
      "RColorBrewer",
      "ggbeeswarm"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = baseline_vs_month3_with_ctrls_up_graph,
    command =
      dplyr::bind_rows(
        dplyr::rename(
          .data = baseline_up_exprs_pivot,
          exprs = baseline_exprs
        ) |>
          dplyr::select(
            sample_name,
            gene,
            exprs,
            subject_ref,
            responder_status,
            timepoint
          ),
        dplyr::rename(
          .data = month3_up_exprs_pivot,
          exprs = month3_exprs
        ) |>
          dplyr::select(
            sample_name,
            gene,
            exprs,
            subject_ref,
            responder_status,
            timepoint
          ),
        dplyr::rename(
          .data = control_up_exprs_pivot,
          exprs = month3_exprs
        ) |>
          dplyr::select(
            sample_name,
            gene,
            exprs,
            subject_ref,
            responder_status,
            timepoint
          )
      ) |>
      dplyr::filter(gene %in% responder_up_genes) |>
      dplyr::mutate(
        responder_timepoint =
          dplyr::if_else(
            condition = responder_status == "control",
            true = "control",
            false = paste(responder_status, timepoint)
          ) |>
          forcats::fct_relevel(
            "control",
            "non_responder baseline",
            "non_responder month 3",
            "responder baseline",
            "responder month 3"
          )
      ) |>
      ggpubr::ggviolin(
        x = "responder_timepoint",
        y = "exprs",
        fill = "responder_timepoint",
        palette = "Paired"
      ) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::facet_wrap(
        facets = ggplot2::vars(gene),
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = "Up-regulated in responders vs non-responders, baseline vs month 3"
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(gene),
        scales = "free_y",
        nrow = 3,
        ncol = 3
      ),
    packages = c(
      "dplyr",
      "tidyr",
      "forcats",
      "ggpubr",
      "ggplot2",
      "RColorBrewer",
      "ggbeeswarm"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = individual_minimum_down_degs_baseline_vs_month3_graph,
    command =
      dplyr::select(
        .data = all_module_scores_with_md,
        sample_name,
        tidyselect::matches("^M[[:digit:]]+"),
        tidyselect::matches("mg"),
        tidyselect::starts_with("ldg"),
        tidyselect::one_of(comparison_groups)
      ) %>%
      tidyr::pivot_longer(
        cols          = c(
          tidyselect::matches("^M[[:digit:]]+"),
          tidyselect::matches("mg"),
          tidyselect::starts_with("ldg")
        ),
        names_to      = "module",
        values_to     = "score"
      ),
    packages =
      c(
        "dplyr",
        "tidyselect",
        "tidyr"
      )
  ),
  targets::tar_target(
    name = individual_minimum_up_degs_baseline_vs_month3_graph,
    command =
      dplyr::bind_rows(
        dplyr::rename(
          .data = baseline_up_exprs_pivot,
          exprs = baseline_exprs
        ) |>
          dplyr::select(
            sample_name,
            gene,
            exprs,
            subject_ref,
            responder_status,
            timepoint
          ),
        dplyr::rename(
          .data = month3_up_exprs_pivot,
          exprs = month3_exprs
        ) |>
          dplyr::select(
            sample_name,
            gene,
            exprs,
            subject_ref,
            responder_status,
            timepoint
          ),
        dplyr::rename(
          .data = control_up_exprs_pivot,
          exprs = month3_exprs
        ) |>
          dplyr::select(
            sample_name,
            gene,
            exprs,
            subject_ref,
            responder_status,
            timepoint
          )
      ) |>
      dplyr::filter(
        .data = module_scores_pivot,
        module %in% annotated_modules[["module"]]
      ),
    packages = "dplyr"
  ),

  targets::tar_target(
    name = minimum_down_degs_all_samples,
    command =
      dplyr::select(
        .data = all_module_scores_with_md,
        sample_name,
        tidyselect::starts_with("ME"),
        tidyselect::one_of(comparison_groups)
      ) %>%
      tidyr::pivot_longer(
        cols          = c(
          tidyselect::starts_with("ME")
        ),
        names_to      = "module",
        values_to     = "score"
      ),
    packages =
      c(
        "dplyr",
        "tidyselect",
        "tidyr"
      )
  ),

  targets::tar_target(
    name = distinguishing_down_metagene_timepoints_graph,
    command =
      minimum_down_degs_all_samples |>
      dplyr::filter(number_degs == minimum_down_degs_no_overlap) |>
      ggpubr::ggboxplot(
        x = "responder_timepoint",
        fill = "responder_timepoint",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::stat_pvalue_manual(
        data = dplyr::filter(
          .data = minimum_down_degs_stats_all_samples,
          number_degs == minimum_down_degs_no_overlap)
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggplot2::labs(title = "Metagene from down-regulated in responders genes"),
    deployment = "main"
  ),
  targets::tar_target(
    name = significant_down_metagene_timepoints_graph,
    command =
      minimum_down_degs_all_samples |>
      dplyr::filter(number_degs == length(responder_down_genes)) |>
      ggpubr::ggboxplot(
        x = "responder_timepoint",
        fill = "responder_timepoint",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::stat_pvalue_manual(
        data = dplyr::filter(
          .data = minimum_down_degs_stats_all_samples,
          number_degs == length(responder_down_genes))
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggplot2::labs(title = "Metagene from down-regulated in responders genes"),
    deployment = "main"
  ),
  targets::tar_target(
    name = up_genes_use_all,
    command =
      magrittr::use_series(up_tables, "responder_vs_non_responder") |>
      dplyr::filter(gene %in% rownames(all_vsc_exprs)) |>
      dplyr::slice_min(
        order_by = log2FoldChange,
        n = up_number_of_genes,
        with_ties = FALSE
      ) |>
      dplyr::slice_min(
        order_by = padj,
        n = up_number_of_genes,
        with_ties = FALSE
      ) |>
      dplyr::pull(gene),
    pattern = map(up_number_of_genes),
    iteration = "list",
    packages = c("magrittr", "dplyr"),
    deployment = "worker"
  ),
  targets::tar_target(
    name = variable_up_deg_pca_all_samples,
    command =
      calculate_module_score(
        expr_mat = all_vsc_exprs,
        gene_list = up_genes_use_all,
        score_func = project_params[["module_scoring_func"]],
        absolute = project_params[["absolute_module_scores"]]
      ) |>
      tibble::as_tibble(
        .name_repair = \(y) {
          y = paste(
            "PC",
            up_number_of_genes,
            sep="_"
          )
        }
      ),
    pattern = map(up_genes_use_all, up_number_of_genes),
    iteration = "list",
    packages = "magrittr",
    deployment = "worker"
  ),
  targets::tar_target(
    name = variable_up_deg_scores_all_samples,
    command =
      dplyr::bind_cols(
        variable_up_deg_pca_all_samples,
        .name_repair = "unique"
      ),
    packages = "dplyr",
    deployment = "main"
  ),

  targets::tar_target(
    name = minimum_up_degs_all_samples,
    command =
      variable_up_deg_tbl |>
      dplyr::left_join(
        dplyr::select(
          .data = all_passing_md,
          sample_name,
          subject_ref,
          visit,
          responder_status
        )
      ) |>
      tidyr::pivot_longer(
        cols = starts_with("PC"),
        names_to = "number_degs",
        values_to = "module_score"
      ) |>
      dplyr::mutate(
        number_degs =
          stringr::str_remove(
            string = number_degs,
            pattern = "^PC"
          ) |>
          as.integer(),
        responder_timepoint =
          dplyr::case_when(
            visit == "control" ~ "control",
            visit == "BL" ~ paste(responder_status, "baseline"),
            visit == "V3" ~ paste(responder_status, "month 3")
          )|>
          forcats::fct_relevel(
            "control",
            "non_responder baseline",
            "non_responder month 3",
            "responder baseline",
            "responder month 3"
          )
      ) |>
      dplyr::select(
        responder_timepoint,
        module_score,
        subject_ref,
        number_degs
      ),
    packages = c(
      "dplyr",
      "tidyr",
      "stringr",
      "forcats"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_up_degs_stats_all_samples,
    command =
      minimum_up_degs_all_samples |>
      dplyr::group_by(number_degs) |>
      rstatix::wilcox_test(
        formula = module_score ~ responder_timepoint,
        p.adjust.method = "fdr"
      ) |>
      rstatix::remove_ns() |>
      rstatix::add_xy_position(
        group = "responder_timepoint",
        scales = "free_y"
      ) |>
      rstatix::add_significance(),
    packages = c(
      "dplyr",
      "rstatix"
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = minimum_up_degs_all_samples_graph,
    command =
      ggpubr::ggboxplot(
        data = minimum_up_degs_all_samples,
        x = "responder_timepoint",
        fill = "responder_timepoint",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm() +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::stat_pvalue_manual(
        data = minimum_up_degs_stats_all_samples
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggforce::facet_wrap_paginate(
        facets = ggplot2::vars(number_degs),
        shrink = TRUE,
        ncol=3,
        nrow=3,
        scales="free_y"
      ) +
      ggplot2::labs(
        title = stringr::str_wrap(
          string =
            paste(
              "Metagene scores using indicated number of most",
              "up-regulated in responders genes"
            )
        ),
        width = 40,
      ),
    deployment = "main"
  ),
  targets::tar_target(
    name = distinguishing_up_metagene_timepoints_graph,
    command =
      minimum_up_degs_all_samples |>
      dplyr::filter(number_degs == minimum_up_degs_no_overlap) |>
      ggpubr::ggboxplot(
        x = "responder_timepoint",
        fill = "responder_timepoint",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::stat_pvalue_manual(
        data =dplyr::filter(
          .data = minimum_up_degs_stats_all_samples,
          number_degs == minimum_up_degs_no_overlap)
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggplot2::labs(title = "Metagene from up-regulated in responders genes"),
    deployment = "main"
  ),
  targets::tar_target(
    name = significant_up_metagene_timepoints_graph,
    command =
      minimum_up_degs_all_samples |>
      dplyr::filter(number_degs == length(responder_up_genes)) |>
      ggpubr::ggboxplot(
        x = "responder_timepoint",
        fill = "responder_timepoint",
        y = "module_score",
        palette = "Set1"
      ) +
      ggbeeswarm::geom_beeswarm(dodge.width = 1) +
      ggplot2::scale_fill_manual(
        values = RColorBrewer::brewer.pal(
          n = 6,
          name = "Paired"
        )[c(6,2,1,4,3)],
        guide = ggplot2::guide_legend(title = NULL)
      ) +
      ggplot2::geom_line(ggplot2::aes(group = subject_ref)) +
      ggpubr::stat_pvalue_manual(
        data = dplyr::filter(
          .data = minimum_up_degs_stats_all_samples,
          number_degs == length(responder_up_genes)
        )
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
      ) +
      ggplot2::labs(title = "Metagene from up-regulated in responders genes"),
    deployment = "main"
  ),

  #### clinical metadata ####
  tarchetypes::tar_file_read(
    name = clinical_md,
    command = project_params[["clinical_file"]],
    read =
      readxl::read_excel(
        path = !!.x,
        sheet = "clinical",
        .name_repair = janitor::make_clean_names
      ),
    packages = "readxl",
    deployment = "main"
  ),

  targets::tar_target(
    name = clinical_md_tbl,
    command = clinical_md |>
      dplyr::filter(visit_ref %in% bl_ctrl_md[["visit_ref"]]) |>
      dplyr::select(
        visit_ref,
        malar_rash,
        discoid_rash,
        photosensitivity,
        oral_ulcers,
        arthritis,
        serositis,
        renal_disorder,
        neurologic_disorder,
        hematologic_disorder,
        immunologic_disorder,
        apls,
        anti_phospholipid_history,
        thrombosis,
        pregnancy,
        thrombocytopenia_itp,
        nephritis,
        ra,
        sjogrens,
        uctd,
        ana,
        selena = selena_sledai_pga_0_3,
        sledai = sledai_total_score,
        bilag = bilag_total,
        acr = acr_total
      ) |>
      dplyr::full_join(
        dplyr::select(
          .data = bl_ctrl_md,
          responder_status,
          age,
          sex,
          race_code,
          visit_ref
        )
      ) |>
      dplyr::mutate(
        dplyr::across(
          .cols = tidyselect::one_of(factor_columns),
          .fns = as.factor
        ),
        dplyr::across(
          .cols = tidyselect::all_of(bool_columns),
          .fns = as.logical
        ),
        race_code = dplyr::case_when(
          race_code == "[EA]" ~ "EA",
          race_code == "[AA]" ~ "AA",
          race_code == "[AI]" ~ "AI",
          race_code == "[A]" ~ "A",
          # race_code == "[H]" ~ "H",
          TRUE ~ "mixed"
        ) |> as.factor(),
        dplyr::across(
          .cols = tidyselect::all_of(score_columns),
          .fns = as.numeric
        )
      ),
    packages = c("dplyr", "tidyselect"),
    deployment = "main"
  ),

  targets::tar_target(
    name = logical_characteristics,
    command =
      purrr::map(
        .x = module_comparisons_stats,
        .f = dplyr::filter,
        module %in% annotated_modules$module
      ),
    packages =
      c(
        "purrr",
        "rlang",
        "dplyr",
        "tidyr",
        "tibble",
        "rstatix"
      )
  ),

  targets::tar_target(
    name = logical_stats_tbl,
    command =
      c(
        generatePalettes(
          .data = study_md,
          c(
            project_params[["comparison_grouping_variable"]],
            one_of(project_params[["heatmap_row_annotations"]])
          )
        ),
        generatePalettes(
          .data = module_annotation,
          .cols = "type",
          use_palettes = "ggsci::default_ucscgb"
        ),
        list(
          "chr" = paletteer::paletteer_d(
            palette = getRandomPalette(2),
            n = 2
          ) |>
            as.character() |>
            rlang::set_names(nm = c("X", "Y"))
        )
      ),
    packages =
      c(
        "purrr",
        "rlang",
        "dplyr",
        "tidyr",
        "tibble",
        "rstatix"
      )
  ),

  targets::tar_target(
    name = continuous_characteristics,
    command =
      clinical_md_tbl %>%
      dplyr::select(
        where(is.numeric)
      ) |>
      colnames(),
    deployment = "main",
    packages = "dplyr"
  ),

  targets::tar_target(
    name = continuous_stats_tbl,
    command =
      clinical_md_tbl |>
      dplyr::select(
        responder_status,
        tidyselect::all_of(continuous_characteristics)
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(continuous_characteristics),
        names_to = "characteristic"
      ) |>
      dplyr::group_by(characteristic) |>
      rstatix::kruskal_test(formula = value ~ responder_status) |>
      dplyr::mutate(
        `P-value` =
          dplyr::case_when(
            p == 1 ~ ">0.99",
            # p >= 0.05 ~ "ns",
            TRUE ~ format(x = p, digits = 2, drop0trailing = TRUE)
          )
      ) |>
      dplyr::select(
        characteristic,
        `P-value`
      ),
    deployment = "main",
    packages =
      c(
        "dplyr",
        "tidyselect",
        "tidyr",
        "rstatix"
      )
  ),

  targets::tar_target(
    name = output_metadata,
    command  =
      writeMetaData(
        object = dataset_with_scores,
        output_name = "processed_data/sample_metadata.csv.gz"
      ),
    format   = "file",
    packages =
      c(
        "data.table",
        "pluck"
      ),
  ),

  targets::tar_target(
    name = symptoms_columns,
    command = c(
      "malar_rash",
      "discoid_rash",
      "photosensitivity",
      "oral_ulcers",
      "arthritis",
      "serositis",
      "renal_disorder",
      "neurologic_disorder",
      "hematologic_disorder",
      "immunologic_disorder",
      "apls",
      "anti_phospholipid_history",
      "thrombosis",
      "pregnancy",
      "thrombocytopenia_itp",
      "nephritis",
      "ra",
      "sjogrens",
      "uctd"
    ),
    deployment = "main"
  ),

  targets::tar_target(
    name = score_columns,
    command = c(
      "selena",
      "sledai",
      "bilag",
      "acr"
    ),
    deployment = "main"
  ),

  targets::tar_target(
    name = pval_tbl,
    command =
      dplyr::bind_rows(
        logical_stats_tbl,
        race_stats_tbl,
        continuous_stats_tbl
      ) |>
      dplyr::mutate(
        Characteristic =
          targets_recode(
            characteristic,
            list(
              "race_code"                 = "Race/Ethnicity, n (%)",
              "ana"                       = " Antinuclear Antibody",
              "anti_phospholipid_history" = " Antiphospholipid Syndrome",
              "discoid_rash"              = " Discoid Rash",
              "hematologic_disorder"      = " Hematologic Disorder",
              "immunologic_disorder"      = " Immunologic Disorder",
              "malar_rash"                = " Malar Rash",
              "nephritis"                 = " Nephritis",
              "oral_ulcers"               = " Oral Ulcers",
              "photosensitivity"          = " Photosensitivity",
              "pregnancy"                 = " Pregnancy",
              "renal_disorder"            = " Renal Disorder",
              "thrombocytopenia_itp"      = " Thrombocytopenia",
              "thrombosis"                = " Thrombosis",
              "acr"                       = "Total ACR score",
              "age"                       = "Age at visit, median (IQR)",
              "bilag"                     = "BILAG index, mean (IQR)",
              "selena"                    = "Physician Global Assessment, mean (IQR)",
              "sledai"                    = "SLEDAI, mean (IQR)"
            )
          )
      ) |>
      select(-characteristic)
  ),

  #### targets for manuscript tables and figures ####
  targets::tar_target(
    name = output_wgcna_module_genes,
    command  =
      writeData(
        object = wgcna_module_genes,
        output_name = "processed_data/wgcna_module_genes.csv.gz"
      ),
    format   = "file",
    packages = "data.table"
  ),

  tar_quarto(
    name = primary_report,
    path          = "analysis/report.qmd",
    packages      =
      c(
        "rmarkdown",
        "knitr",
        "targets",
        "here",
        "ggplot2",
        "ggpubr",
        "rlang",
        "magrittr",
        "janitor",
        "kableExtra",
        "flextable",
        "genefilter",
        "tidyselect",
        "purrr",
        "formattable",
        "grid",
        "stringr"
      ),
    quiet = FALSE
  )
)
