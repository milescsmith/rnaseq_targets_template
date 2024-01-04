readr::read_csv("vsc_exprs.csv/vsc_exprs.csv") #use this as example data
library(tidyverse)
library(reticulate)
library(targets)

tar_load(vsc_exprs)
tar_load(bl_ctrl_md)
tar_load(up_genes_use)
tar_load(down_genes_use)
use_condaenv("reticulate")



hsic_up_module <- calculate_module_score(
    expr_mat = vsc_exprs,
    gene_list = hsic_genes,
    score_func = "nmf"
  )

  hsic_up_scores <-
    hsic_up_module$scores |>
    as_tibble(rownames = "sample_name") |>
    magrittr::set_colnames(c("sample_name", "module_score")) |>
    dplyr::filter(sample_name %in% bl_ctrl_md[["sample_name"]]) |>
    dplyr::left_join(bl_ctrl_md) |>
    dplyr::select(sample_name, module_score, responder_status)

hsic_stats <- hsic_up_scores |>
  rstatix::wilcox_test(
    formula = module_score ~ responder_status,
    p.adjust.method = "fdr"
  ) |>
  rstatix::remove_ns() |>
  rstatix::add_xy_position(
    group = "responder_status",
    scales = "free_y"
  ) |>
  rstatix::add_significance()

ggpubr::ggboxplot(
  data = hsic_up_scores,
  x = "responder_status",
  fill = "responder_status",
  y = "module_score",
  palette = "Set1"
) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(
    data = hsic_stats
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
  ggplot2::labs(title = "Metagene scores using indicated number of most up-regulated in responders genes")
