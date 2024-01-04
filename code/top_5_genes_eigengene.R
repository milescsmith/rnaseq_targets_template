vsc_tbl <-
  full_vsc_exprs %>%
  t() %>%
  as_tibble(rownames = "sample_name")

filtered_res <-
  magrittr::use_series(res, "responder_vs_non_responder") |>
  filter(
    padj < 0.01,
    gene %in% colnames(vsc_tbl)
    ) %>%
  left_join(
    HGNChelper::checkGeneSymbols(.[["gene"]]) %>% rename(gene = x)
  )

x <- 5

top_5_down <-
  magrittr::use_series(res, "responder_vs_non_responder") |>
  filter(
    padj < 0.01,
    gene %in% colnames(vsc_tbl)
  ) |>
  slice_min(
    order_by = log2FoldChange,
    n = x
  ) %>%
  slice_min(
    order_by = padj,
    n = x
  ) |>
  pull(gene)

approved_top_5_down <-
  filtered_res |>
  filter(
    Approved == TRUE
  ) |>
  slice_min(
    order_by = log2FoldChange,
    n = x
  ) %>%
  slice_min(
    order_by = padj,
    n = x
  ) |>
  pull(gene)

top_5_up <-
  magrittr::use_series(res, "responder_vs_non_responder") |>
  filter(
    padj < 0.01,
    gene %in% colnames(vsc_tbl)
    ) |>
  slice_max(
    order_by = log2FoldChange,
    n = x
  ) %>%
  slice_min(
    order_by = padj,
    n = x
  ) |>
  pull(gene)

approved_top_5_up <-
  filtered_res |>
  filter(
    Approved == TRUE
  ) |>
  slice_max(
    order_by = log2FoldChange,
    n = x
  ) %>%
  slice_min(
    order_by = padj,
    n = x
  ) |>
  pull(gene)

plot_metagene <- function(genes_use){

  all_pca <-
    rsvd::rsvd(
      A = full_vsc_exprs[
        genes_use,
      ],
      k = 1
    )

  gene_loadings <-
    magrittr::set_rownames(
      x = all_pca$u,
      value = genes_use
    ) |>
    magrittr::set_colnames(value="PC1") |>
    tibble::as_tibble(rownames="gene")

  pc1_score <-
    magrittr::set_rownames(
      x = all_pca$v,
      value = colnames(full_vsc_exprs)
    ) |>
    magrittr::set_colnames(value="PC1") |>
    tibble::as_tibble(rownames="sample_name") |>
    dplyr::left_join(full_md) |>
    filter(
      group == "control" | str_detect(Sample_Alias, "03$") | str_detect(Sample_Alias, "BSL$")
    ) |>
    mutate(
      timepoint =
        case_when(
          group == "control" ~ "none",
          str_detect(Sample_Alias, "03$") ~ "month 3",
          str_detect(Sample_Alias, "BSL$") ~ "baseline"
        ),
      responder_timepoint =
        if_else(
          condition = responder == "control",
          true = "control",
          false = paste(responder, timepoint)
        ) |>
        fct_relevel(
          "control",
          "non_responder baseline",
          "non_responder month 3",
          "responder baseline",
          "responder month 3"
        )
    )

  pc1_stats <-
    pc1_score %>%
    rstatix::wilcox_test(
      formula = PC1 ~ responder_timepoint,
      p.adjust.method = "fdr"
    ) |>
    rstatix::add_xy_position(
      group = "responder_timepoint",
      scales = "free_y",
    ) |>
    rstatix::add_significance()

  ggpubr::ggviolin(
    data = pc1_score,
    x="responder_timepoint",
    y="PC1",
    fill="responder_timepoint"
  ) +
    ggbeeswarm::geom_beeswarm(dodge.width=1) +
    geom_line(aes(group=SubjectRef))+
    ggpubr::stat_pvalue_manual(data = pc1_stats) +
    scale_fill_manual(
      values = RColorBrewer::brewer.pal(n=6, name="Paired")[c(4,2,1,6,5)],
      guide = guide_legend(nrow = 5, byrow = TRUE, title = NULL)
    ) +
    ggpubr::theme_pubr(base_size = 9) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "right"
    )
}

cowplot::plot_grid(
  plot_metagene(top_5_down) +
    labs(
      title = str_wrap("Metagene scores using the 5 most down-regulated in genes responders", width = 50),
      x = "",
      y = "Eigengene score"
      ) +
    theme(
      title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
      ),
  plot_metagene(top_5_up) +
    labs(
      title = str_wrap("Metagene scores using the 5 most up-regulated in genes responders", width = 50),
      x = "",
      y = "Eigengene score"
    ) +
    theme(
      title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
    ),
  ncol=2
  )


full_vsc_exprs[top_5_down,] %>%
  t() %>%
  as_tibble(rownames = "sample_name") %>%
  dplyr::left_join(full_md) |>
  filter(
    group == "control" | str_detect(Sample_Alias, "03$") | str_detect(Sample_Alias, "BSL$")
  ) |>
  mutate(
    timepoint =
      case_when(
        group == "control" ~ "none",
        str_detect(Sample_Alias, "03$") ~ "month 3",
        str_detect(Sample_Alias, "BSL$") ~ "baseline"
      ),
    responder_timepoint =
      if_else(
        condition = responder == "control",
        true = "control",
        false = paste(responder, timepoint)
      ) |>
      fct_relevel(
        "control",
        "non_responder baseline",
        "non_responder month 3",
        "responder baseline",
        "responder month 3"
      )
  ) %>%
  write_csv(file = "top_5_down_normalized_expression.csv")

full_vsc_exprs[top_5_up,] %>%
  t() %>%
  as_tibble(rownames = "sample_name") %>%
  dplyr::left_join(full_md) |>
  filter(
    group == "control" | str_detect(Sample_Alias, "03$") | str_detect(Sample_Alias, "BSL$")
  ) |>
  mutate(
    timepoint =
      case_when(
        group == "control" ~ "none",
        str_detect(Sample_Alias, "03$") ~ "month 3",
        str_detect(Sample_Alias, "BSL$") ~ "baseline"
      ),
    responder_timepoint =
      if_else(
        condition = responder == "control",
        true = "control",
        false = paste(responder, timepoint)
      ) |>
      fct_relevel(
        "control",
        "non_responder baseline",
        "non_responder month 3",
        "responder baseline",
        "responder month 3"
      )
  ) %>%
  write_csv(file = "top_5_up_normalized_expression.csv")
