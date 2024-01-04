# get maximum values for the controls
max_top5_down_values <-
  controls_top5_down_exprs_pivot |>
  filter(responder == "control") |>
  group_by(gene) |>
  summarise(
    max_expr = max(baseline_exprs),
    mean_expr = mean(baseline_exprs),
    std_expr = sd(baseline_exprs)
  ) |>
  select(gene, max_expr)

baseline_down_exprs <- full_vsc_exprs[top_5_down, baseline_md[["sample_name"]]]
baseline_down_exprs_pivot <-
  baseline_exprs |>
  as_tibble(rownames = "gene") |>
  pivot_longer(
    cols = one_of(baseline_md[["sample_name"]]),
    names_to = "sample_name",
    values_to = "baseline_exprs"
  ) |>
  left_join(
    select(
      .data = baseline_md,
      sample_name,
      SubjectRef,
      responder
    )
  ) |>
  left_join(max_top5_down_values) |>
  mutate(
    baseline_diff = baseline_exprs - max_expr,
    timepoint = "baseline"
    )

month3_down_exprs <- full_vsc_exprs[top_5_down, month3_md[["sample_name"]]]
month3_down_exprs_pivot <-
  month3_down_exprs |>
  as_tibble(rownames = "gene") |>
  pivot_longer(
    cols = one_of(month3_md[["sample_name"]]),
    names_to = "sample_name",
    values_to = "month3_exprs"
  ) |>
  left_join(
    select(
      .data = month3_md,
      sample_name,
      SubjectRef,
      responder
    )
  ) |>
  left_join(max_top5_down_values) |>
  mutate(
    month3_diff = month3_exprs - max_expr,
    timepoint = "month 3"
    )

baseline_vs_month3 <-
  inner_join(
    x = baseline_down_exprs_pivot,
    y = month3_down_exprs_pivot,
    by = c("SubjectRef", "gene", "responder", "max_expr"),
    keep = FALSE
    )

baseline_vs_month3 |>
  pivot_longer(
    cols = c(baseline_diff, month3_diff),
    names_pattern = "(\\w+)_",
    names_to = "timepoint",
    values_to = "difference"
  ) |>
  mutate(
    responder_timepoint = paste(responder, timepoint) |>
      fct_relevel(
        "non_responder baseline",
        "non_responder month3",
        "responder baseline",
        "responder month3"
      )
  ) |>
  ggpubr::ggviolin(
    x="responder_timepoint",
    y="difference",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm(dodge.width=1) +
  geom_line(aes(group=SubjectRef)) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  facet_wrap(facets = vars(gene), scales = "free_y") +
  ggtitle("Down-regulated in responders vs non-responders, baseline vs month 3")
