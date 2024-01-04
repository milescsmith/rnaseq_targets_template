vsc_exprs <- full_data$variance_stabilized_counts

full_md <-
  full_data$dataset$samples |>
  as_tibble(rownames="sample_name") |>
  select(
    -VisitRef,
    -SubjectRef
    ) |>
  left_join(novaseq_ref_list)

control_md <- full_md |> filter(group == "control")
baseline_md <- full_md |> filter(str_detect(Sample_Alias, "BSL$"))
month3_md <- full_md |> filter(str_detect(Sample_Alias, "03$"))

baseline_top5_up_exprs <- vsc_exprs[top_5_up, c(control_md[["sample_name"]], baseline_md[["sample_name"]])]
baseline_top5_down_exprs <- vsc_exprs[top_5_down,c(control_md[["sample_name"]], baseline_md[["sample_name"]])]


baseline_top5_down_exprs_pivot <-
  baseline_top5_down_exprs |>
  as_tibble(rownames="gene") |>
  pivot_longer(
    -gene,
    names_to = "sample_name",
    values_to = "exprs"
  ) |>
  left_join(
    select(
      full_md,
      sample_name,
      responder
    )
  ) |>
  mutate(
    timepoint = if_else(
      condition = responder == "control",
      true = "none",
      false = "baseline"
    )
  ) |>
  left_join(select(full_md, sample_name, SubjectRef))

baseline_top5_up_exprs_pivot <-
  baseline_top5_up_exprs |>
  as_tibble(rownames="gene") |>
  pivot_longer(
    -gene,
    names_to = "sample_name",
    values_to = "exprs"
  ) |>
  left_join(
    select(
      full_md,
      sample_name,
      responder
    )
  ) |>
  mutate(
    timepoint = if_else(
      condition = responder == "control",
      true = "none",
      false = "baseline"
    )
  ) |>
  left_join(select(full_md, sample_name, SubjectRef))

baseline_top5_down_exprs_stats <-
  baseline_top5_down_exprs_pivot |>
  group_by(gene) |>
  rstatix::wilcox_test(formula = exprs ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

bsl_5_down <- baseline_top5_down_exprs_pivot |>
  ggpubr::ggboxplot(
    x="responder",
    y="exprs",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = baseline_top5_down_exprs_stats) +
  facet_wrap(vars(gene), scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  ggtitle("Top 5 down-regulated in responders vs non-responders at baseline")

baseline_top5_up_exprs_stats <-
  baseline_top5_up_exprs_pivot |>
  group_by(gene) |>
  rstatix::wilcox_test(formula = exprs ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

bsl_5_up <- baseline_top5_up_exprs_pivot |>
  ggpubr::ggboxplot(
    x="responder",
    y="exprs",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = baseline_top5_up_exprs_stats) +
  facet_wrap(vars(gene), scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  ggtitle("Top 5 up-regulated in responders vs non-responders at baseline")

month3_top5_up_exprs <- vsc_exprs[top_5_up, c(control_md[["sample_name"]], month3_md[["sample_name"]])]

month3_top5_down_exprs <- vsc_exprs[top_5_down,c(control_md[["sample_name"]], month3_md[["sample_name"]])]

month3_top5_down_exprs_pivot <-
  month3_top5_down_exprs |>
  as_tibble(rownames="gene") |>
  pivot_longer(
    -gene,
    names_to = "sample_name",
    values_to = "exprs"
  ) |>
  left_join(
    select(
      full_md,
      sample_name,
      responder
    )
  ) |>
  mutate(
    timepoint = if_else(
      condition = responder == "control",
      true = "none",
      false = "month 3"
    )
  ) |>
  left_join(select(full_md, sample_name, SubjectRef))

month3_top5_up_exprs_pivot <-
  month3_top5_up_exprs |>
  as_tibble(rownames="gene") |>
  pivot_longer(
    -gene,
    names_to = "sample_name",
    values_to = "exprs"
  ) |>
  left_join(
    select(
      full_md,
      sample_name,
      responder
    )
  ) |>
  mutate(
    timepoint = if_else(
      condition = responder == "control",
      true = "none",
      false = "month 3"
    )
  ) |>
  left_join(select(full_md, sample_name, SubjectRef))


month3_top5_down_exprs_stats <-
  month3_top5_down_exprs_pivot |>
  group_by(gene) |>
  rstatix::wilcox_test(formula = exprs ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

month_3_down <- month3_top5_down_exprs_pivot |>
  ggpubr::ggboxplot(
    x="responder",
    y="exprs",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = month3_top5_down_exprs_stats) +
  facet_wrap(vars(gene), scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  ggtitle("Top 5 down-regulated in responders vs non-responders at month 3")

month3_top5_up_exprs_stats <-
  month3_top5_up_exprs_pivot |>
  group_by(gene) |>
  rstatix::wilcox_test(formula = exprs ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

month_3_up <- month3_top5_up_exprs_pivot |>
  ggpubr::ggboxplot(
    x="responder",
    y="exprs",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = month3_top5_up_exprs_stats) +
  facet_wrap(vars(gene), scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  ggtitle("Top 5 up-regulated in responders vs non-responders at month 3")

cowplot::plot_grid(bsl_5_down, month_3_down)
cowplot::plot_grid(bsl_5_up, month_3_up)

top5_down_exprs_pivot <-
  bind_rows(
    baseline_top5_down_exprs_pivot,
    month3_top5_down_exprs_pivot
    ) |>
  distinct()

top5_up_exprs_pivot <-
  bind_rows(
    baseline_top5_up_exprs_pivot,
    month3_top5_up_exprs_pivot
  ) |>
  distinct()

top5_down_exprs_pivot |>
  mutate(
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
  ) |>
  ggpubr::ggviolin(
    x="responder_timepoint",
    y="exprs",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm(dodge.width=1) +
  geom_line(aes(group=SubjectRef))+
  # ggpubr::stat_pvalue_manual(data = month3_top5_up_exprs_stats) +
  facet_wrap(vars(gene), scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  ggtitle("Top 5 up-regulated in responders vs non-responders at month 3")
