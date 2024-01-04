#### Set options ####
options(future.globals.maxSize    = +Inf)
Sys.setenv('RSTUDIO_PANDOC'       = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()
BPPARAM =
  BiocParallel::SnowParam(
    workers       = parallel::detectCores()-2,
    exportglobals = FALSE,
    progressbar   = TRUE,
    type = "SOCK"
    )

BiocParallel::register(BPPARAM)

`%nin%` <- purrr::negate(`%in%`)
`%||%` <- rlang::`%||%`

#### Setup project variables ####
project_params = list(

  sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/novaseq/",
  metadata_file                   = "metadata/metadata.csv",
  metadata_sheet                  = NULL,
  skip_lines                      = 0,

  main_sample_sheet               = NULL,
  main_sample_sheet_skip          = 0,
  extra_controls_sample_sheet     = "metadata/NovaSeq_Sample_List.xlsx",
  extra_controls_ident_col        = "disease_class",
  extra_controls_iden             = "control",
  annotation_file                 = "references/gencode_v32_virus_tx2gene_v1.2.csv",
  clinical_file                   = "metadata/BLAST_200312.xlsx",

  #### Metadata columns ####
  sample_species                  = "Homo sapiens",
  sample_name_column              = "sample_name",
  grouping_column                 = "responder_status",
  project_column                  = "study",
  regression_columns              = c("age", "sex", "race_code", "final_concentration_ng_ul"),
  filter_column                   = "initial_concentration_ng_ul",
  filter_value                    = 1.5,
  extra_columns                   = c("visit_ref", "subject_ref", "sample_alias", "visit"),
  heatmap_row_annotations         = c("responder", "sex",  "race_code", "cluster"),

  #### Setup project variables ####
  projects_to_include             = "BLAST",
  projects_to_exclude             = c("ALE06", "Xencor"),

  groups_to_include               = c("responder", "non_responder", "control"),
  groups_to_exclude               = NULL,
  study_design                    = ~ group,

  comparison_grouping_variable    = "group",
  comparison_groups               = "group",
  batch_variable                  = NULL,
  control_group                   = "case",
  experimental_group              = "control",
  manual_sample_removal           = NULL,

  aligner                         = "salmon",
  only_hugo_named_genes           = FALSE,
  minimum_gene_count              = 1,
  initial_concentration_threshold = 1.5,
  pc1_zscore_threshold            = 2,
  pc2_zscore_threshold            = 2.5,
  sva_num                         = 3,
  use_combat                      = FALSE,
  process_method                  = "limma",
  lfcThreshold                    = 0.01,
  deg_substantial_threshold       = 0.01, # which is about a 2.5-fold change
  deg_padj_threshold              = 0.1,
  module_scoring_func             = "nmf",
  absolute_module_scores          = FALSE,
  hsic_gene_set_length            = 2,

  BPPARAM                         = BPPARAM,

  # number of top variable genes to use for WGCNA
  n_var_genes                     = 20000,

  banchereau_modules              = "references/banchereau_modules.csv",
  banchereau_module_annotations   = "references/banchereau_module_annotations.csv",
  ldg_modules                     = "references/ldg_modules.csv",
  metasignature_modules           = "references/metasignature_module.csv",

  cytokine_data                   = "metadata/cyto_ana_data.csv",
  antibody_data                   = "metadata/cyto_ana_data.csv"
  )
