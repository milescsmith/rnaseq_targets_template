group_pal2$direction <- unlist(list(up = RColorBrewer::brewer.pal(4, "Set1")[[3]], down = "#984EA3"))
group_pal2$responder <- unlist(list(responder = "#377EB8", non_responder = "#E41A1C"))

vsc_exprs[up_tables$`responder - non_responder`$gene,] |>
  t() |>
  ComplexHeatmap::Heatmap() +
  ComplexHeatmap::HeatmapAnnotation(df = tibble::column_to_rownames(dplyr::select(study_md, group, sample_name), "sample_name"), which = "row")

tirosh_score_modules(
  expr_obj = vsc_exprs,
  module_list = list(up_tables$`responder - non_responder`$gene)
) |>
  tibble::as_tibble(rownames = "sample_name") |>
  dplyr::left_join(
    dplyr::select(
      study_md,
      sample_name,
      responder
    )
  )

tirosh_score_modules(
  expr_obj = vsc_exprs,
  module_list = list(up_tables$`responder - non_responder`$gene)
) |>
  tibble::as_tibble(rownames = "sample_name") |>
  dplyr::left_join(
    dplyr::select(
      study_md,
      sample_name,
      responder
    )
  ) |>
  ggpubr::ggboxplot(
    x = "responder",
    y = "V1",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm()

tirosh_score_modules(
  expr_obj = vsc_exprs,
  module_list = list(down_tables$`responder - non_responder`$gene)
) |>
  tibble::as_tibble(rownames = "sample_name") |>
  dplyr::left_join(
    dplyr::select(
      study_md,
      sample_name,
      responder
    )
  ) |>
  ggpubr::ggboxplot(
    x = "responder",
    y = "V1",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm()

DEG_modules <-
  tibble::tibble(
    sample_name = colnames(vsc_exprs),
    DEGUM = tirosh_score_modules(
      expr_obj = vsc_exprs,
      module_list = list(up_tables$`responder - non_responder`$gene)
    ) |> purrr::pluck(1),
    DEGDM = tirosh_score_modules(
      expr_obj = vsc_exprs,
      module_list = list(down_tables$`responder - non_responder`$gene)
    ) |> purrr::pluck(1)
  ) |>
  dplyr::left_join(
    dplyr::select(study_md, sample_name, responder)
  )


ggplot(
  data = DEG_modules,
  mapping =
    aes(
      x = DEGUM,
      y = DEGDM,
      color = responder
    )
  ) +
  geom_point()

groupedComplexHeatmap(
  expr_mat            = vsc_exprs,
  gene_list           = c("BATF2", "CXCL10", "DDX60", "EPSTI1", "HERC5", "HES4", "IFI44", "IFI44L", "IFIT1", "IFIT3", "IFITM3", "ISG15", "LAMP3", "LY6E", "MX1", "OAS1", "OAS2", "OAS3", "OASL", "OTOF", "RSAD2", "RTP4", "SERPING1", "TRIM6", "XAF1", "AIM2", "APOL6", "ATF3", "CARD17", "CCL8", "CEACAM1", "DDX58", "DHX58", "EIF2AK2", "FBXO6", "GALM", "GBP1", "GBP3", "GBP4", "GBP5", "GBP6", "HELZ2", "HERC6", "IDO1", "IFI35", "IFIH1", "IFIT2", "IFIT5", "IFITM1", "IRF7", "LAP3", "LGALS3BP", "MOV10", "MT2A", "PARP10", "PARP12", "PARP14",
                          "PARP9", "PLSCR1", "PML", "SAMD9L", "SCO2", "SOCS1", "STAT1", "STAT2", "TIMM10", "TNFAIP6", "TNFSF10", "TRIM22", "UBE2L6", "WARS", "ZBP1", "ZNF684", "ABCA1", "ACTA2", "ADAR", "BST2", "BTN3A1", "C1QA", "CASP1", "CHMP5", "CPT1B", "DHRS9", "DRAP1", "DYNLT1", "ETV7", "GADD45B", "GBP2", "HSH2D", "IFI16", "IRF9", "ISG20", "LGALS9", "LHFPL2", "MDK", "NBN", "NCOA7", "NMI", "NT5C3A", "NTNG2", "PHF11", "PSMB9", "RBCK1", "REC8", "RHBDF2", "SAMD9", "SP100", "SP110", "SP140", "SRBD1", "TAP1", "TAP2", "TCN2",
                          "TDRD7", "TMEM140", "TRAFD1", "TRANK1", "TRIM21", "TRIM25", "TRIM38", "TRIM5", "TRIM56", "TYMP", "UNC93B1", "ZC3HAV1", "ZNFX1"),
  md                  = study_md,
  annotation_palettes = group_pal2,
  row_grouping        = "responder",
  row_annotation      = c("responder", "RaceCode"),
  col_grouping        = "module",
  col_annotation      = tibble::as_tibble(module_ISGs, rownames="obsv"),
  scale_exprs         = TRUE,
  col_fontsize        = 6
)

groupedComplexHeatmap(
  expr_mat            = vsc_exprs,
  gene_list           = c("AGAP13P","AC040169.1","PRR27","AC016542.1","ARLNC1","AC026150.2",
                          "AC109992.2","AC108125.1","ANKRD18DP","RN7SL801P","AP005482.3","AC090950.1",
                          "AC100814.1","OR13K1P","AC016168.4","MAP6D1","ZC2HC1B","AC007923.1",
                          "TAF1A-AS1","VN2R10P","MAST1","AC010267.1","SDC1","SHOX2",
                          "AL139280.2","AC138965.2","RPL22P1","AL121845.2","AC027644.4","AC069287.1",
                          "NUTM2E","GIMAP3P","NBPF13P","FAM230H","RGPD3","AC011497.2",
                          "AC021148.1","CYP27C1","AC009879.3","PMS2P11","AL158207.1","AC234775.4",
                          "NF1P8","RCC2P6","AL354863.1","C1QTNF4","OSCP1","SNORA12",
                          "AC133644.1","TMEM249","AC020658.2","FAM170B","AC007192.2"),
  md                  = study_md,
  annotation_palettes = group_pal2,
  row_grouping        = "responder",
  row_annotation      = c("responder", "RaceCode"),
  scale_exprs         = TRUE,
  col_annotation      = deg_list,
  col_grouping        = "direction"
)

deg_module_stats <-
  deg_modules |>
  pivot_longer(
    -sample,
    names_to = "direction",
    values_to = "expr"
  ) |>
  rename(sample = "sample_name") |>
  left_join(
    as_tibble(select(study_md, sample_name, responder))
  ) |>
  group_by(direction) |>
  wilcox_test(formula = expr ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |> rstatix::add_significance()

deg_modules |>
  pivot_longer(
    -sample,
    names_to = "direction",
    values_to = "expr"
  ) |>
  rename(sample = "sample_name") |>
  left_join(
    as_tibble(select(study_md, sample_name, responder))
  ) |>
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "expr", palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = deg_module_stats) +
  facet_wrap(facets = vars(direction))

tibble::enframe(
  x     = wgcna_modules[["colors"]],
  name = "gene",
  value = "module"
)

merge$newMEs |> as_tibble(rownames = "sample_name") |> pivot_longer(-sample_name, names_to = "module", values_to = "eigen") |> left_join(select(study_md, sample_name, responder)) |> group_by(module) |>  wilcox_test(formula = eigen ~ responder) |> add_significance() |> filter(p.signif != "ns")

merge =
  mergeCloseModules(
    vsc_top,
    wgcna_modules$colors,
    cutHeight = 0.25,
    verbose = 3
  )

responder_df =
  column_to_rownames(
    select(
      study_md,
      sample_name,
      responder
    ),
    var = "sample_name"
  )

ComplexHeatmap::Heatmap(
  matrix = t(plier_c2_react$B),
  col = viridisLite::cividis(n = 3),
  border = TRUE,
  column_names_rot = 45,
  right_annotation =
    ComplexHeatmap::rowAnnotation(
      df = responder_df
    ),
  row_split = responder_df
)


plier_c2_react_tbl <-
  plier_c2_react$B |>
  as_tibble(rownames = "pathway") |>
  mutate(
    pathway =
      str_replace_all(
        string = pathway,
        pattern = "_",
        replacement = " "
      )
  ) |>
  pivot_longer(
    -pathway,
    names_to = "sample_name",
    values_to = "expr"
  ) |>
  left_join(as_tibble(responder_df, rownames = "sample_name"))

plier_c2_react_stats <-
  plier_c2_react_tbl |>
  group_by(pathway) |>
  wilcox_test(formula = expr ~ responder) |>
  add_xy_position()

plier_c2_react_tbl |>
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "expr",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(plier_c2_react_stats) +
  facet_wrap(
    facets = vars(pathway),
    scales = "free_y",
    labeller =
      label_wrap_gen(
        width = 25,
        multi_line = TRUE
      ),
    nrow = 3
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 9
      ),
    axis.text.y =
      element_text(
        size = 9
      ),
    strip.text.x =
      element_text(
        size = 9
      )
  ) +
  scale_y_continuous(expand = expansion(mult = 0.2))

plier_c7_immune

ComplexHeatmap::Heatmap(
  matrix = t(plier_c7_immune$B),
  col = viridisLite::cividis(n = 3),
  border = TRUE,
  column_names_rot = 45,
  right_annotation =
    ComplexHeatmap::rowAnnotation(
      df = responder_df
    ),
  row_split = responder_df
)

plier_c7_immune_tbl <-
  plier_c7_immune$B |>
  as_tibble(rownames = "pathway") |>
  mutate(
    pathway =
      str_replace_all(
        string = pathway,
        pattern = "_",
        replacement = " "
      )
  ) |>
  pivot_longer(
    -pathway,
    names_to = "sample_name",
    values_to = "expr"
  ) |>
  left_join(as_tibble(responder_df, rownames = "sample_name"))

plier_c7_immune_stats <-
  plier_c7_immune_tbl |>
  group_by(pathway) |>
  wilcox_test(formula = expr ~ responder) |>
  add_xy_position()

walk(seq(2), \(x) {
  boxplt <-
    ggpubr::ggboxplot(
      data = plier_c7_immune_tbl,
      x = "responder",
      fill = "responder",
      y = "expr",
      palette = "Set1"
    ) +
    ggbeeswarm::geom_beeswarm() +
    ggpubr::stat_pvalue_manual(data = plier_c7_immune_stats) +
    facet_wrap(
      facets = vars(pathway),
      scales = "free_y",
      ,
      nrow = 3
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          vjust = 1,
          hjust = 1,
          size = 9
        ),
      axis.text.y =
        element_text(
          size = 9
        ),
      strip.text.x =
        element_text(
          size = 9
        )
    ) +
    scale_y_continuous(expand = expansion(mult = 0.2))
  print(boxplt)
}


c5_mf <- read.gmt("references/c5.go.mf.v7.4.symbols.gmt")
c5_mf_mat <-
  c5_mf %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = "term",
    values_from = "present",
    values_fill = 0
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

plier_c5_mf <- PLIER(vsc_exprs, c5_mf_mat)

plier_c5_mf_tbl <-
  plier_c5_mf$B |>
  as_tibble(rownames = "pathway") |>
  mutate(
    pathway =
      str_replace_all(
        string = pathway,
        pattern = "_",
        replacement = " "
      )
  ) |>
  pivot_longer(
    -pathway,
    names_to = "sample_name",
    values_to = "expr"
  ) |>
  left_join(as_tibble(responder_df, rownames = "sample_name"))

plier_c5_mf_stats <-
  plier_c5_mf_tbl |>
  group_by(pathway) |>
  wilcox_test(formula = expr ~ responder) |>
  add_xy_position()


ggpubr::ggboxplot(
    data = plier_c5_mf_tbl,
    x = "responder",
    fill = "responder",
    y = "expr",
    palette = "Set1"
  ) +
    ggbeeswarm::geom_beeswarm() +
    ggpubr::stat_pvalue_manual(data = plier_c5_mf_stats) +
    facet_wrap(
      facets = vars(pathway),
      scales = "free_y",
      labeller =
        label_wrap_gen(
          width = 25,
          multi_line = TRUE
        ),
      nrow = 3
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          vjust = 1,
          hjust = 1,
          size = 9
        ),
      axis.text.y =
        element_text(
          size = 9
        ),
      strip.text.x =
        element_text(
          size = 9
        )
    ) +
    scale_y_continuous(expand = expansion(mult = 0.2))

c5_bp <- read.gmt("references/c5.go.bp.v7.4.symbols.gmt")
c5_bp_mat <-
  c5_bp %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = "term",
    values_from = "present",
    values_fill = 0
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

plier_c5_bp <- PLIER(vsc_exprs, c5_bp_mat)

plier_c5_bp_tbl <-
  plier_c5_bp$B |>
  as_tibble(rownames = "pathway") |>
  mutate(
    pathway =
      str_replace_all(
        string = pathway,
        pattern = "_",
        replacement = " "
      )
  ) |>
  pivot_longer(
    -pathway,
    names_to = "sample_name",
    values_to = "expr"
  ) |>
  left_join(as_tibble(responder_df, rownames = "sample_name"))

plier_c5_bp_stats <-
  plier_c5_bp_tbl |>
  group_by(pathway) |>
  wilcox_test(formula = expr ~ responder) |>
  add_xy_position()

ggpubr::ggboxplot(
  data = plier_c5_bp_tbl,
  x = "responder",
  fill = "responder",
  y = "expr",
  palette = "Set1"
) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = plier_c5_bp_stats) +
  facet_wrap(
    facets = vars(pathway),
    scales = "free_y",
    labeller =
      label_wrap_gen(
        width = 25,
        multi_line = TRUE
      ),
    nrow = 3
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 9
      ),
    axis.text.y =
      element_text(
        size = 9
      ),
    strip.text.x =
      element_text(
        size = 9
      )
  ) +
  scale_y_continuous(expand = expansion(mult = 0.2))

c3_reg <- read.gmt("references/c3.all.v7.4.symbols.gmt")
c3_reg_mat <-
  c3_reg %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = "term",
    values_from = "present",
    values_fill = 0
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

plier_c3_reg <- PLIER(vsc_exprs, c3_reg_mat)

plier_c3_reg_tbl <-
  plier_c3_reg$B |>
  as_tibble(rownames = "pathway") |>
  mutate(
    pathway =
      str_replace_all(
        string = pathway,
        pattern = "_",
        replacement = " "
      )
  ) |>
  pivot_longer(
    -pathway,
    names_to = "sample_name",
    values_to = "expr"
  ) |>
  left_join(as_tibble(responder_df, rownames = "sample_name"))

plier_c3_reg_stats <-
  plier_c3_reg_tbl |>
  group_by(pathway) |>
  wilcox_test(formula = expr ~ responder) |>
  add_xy_position()


ggpubr::ggboxplot(
  data = plier_c3_reg_tbl,
  x = "responder",
  fill = "responder",
  y = "expr",
  palette = "Set1"
) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = plier_c3_reg_stats) +
  facet_wrap(
    facets = vars(pathway),
    scales = "free_y",
    labeller =
      label_wrap_gen(
        width = 25,
        multi_line = TRUE
      ),
    nrow = 3
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 9
      ),
    axis.text.y =
      element_text(
        size = 9
      ),
    strip.text.x =
      element_text(
        size = 9
      )
  ) +
  scale_y_continuous(expand = expansion(mult = 0.2))



wgcna_stats <-
  wgcna_tbl |>
  group_by(pathway) |>
  wilcox_test(formula = expr ~ responder) |>
  add_xy_position()

ggpubr::ggboxplot(
  data = wgcna_tbl,
  x = "responder",
  fill = "responder",
  y = "expr",
  palette = "Set1"
) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = wgcna_stats) +
  facet_wrap(
    facets = vars(pathway),
    scales = "free_y",
    labeller =
      label_wrap_gen(
        width = 25,
        multi_line = TRUE
      ),
    nrow = 3
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 9
      ),
    axis.text.y =
      element_text(
        size = 9
      ),
    strip.text.x =
      element_text(
        size = 9
      )
  ) +
  scale_y_continuous(expand = expansion(mult = 0.2))


filter(
  .data = MEplotting,
  module %in% unique(MEplotting[["module"]])[1:9]
  ) |>
  ggplot(
    mapping =
      aes(
        x = GeneRatio,
        y = Description,
        color = p.adjust
        )
    ) +
  geom_point() +
  scale_color_viridis_c(guide = guide_colourbar(barwidth = 10)) +
  scale_y_discrete(
    labels = stringr::str_wrap(
      paste(
        MEplotting[["ID"]],
        MEplotting[["Description"]],
        sep = ":"
      ),
      width = 40
    )
  ) +
  facet_wrap(
    facets = vars(module),
    scales = "free",
    ncol = 3,
    nrow = 3
    ) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x =
      element_text(
        size = 9
      ),
    axis.text.y =
      element_text(
        size = 9
      ),
    strip.text =
      element_text(
        size = 9
      )
  )

filter(
  .data = MEplotting,
  module %in% unique(MEplotting[["module"]])[10:16]
) |>
  ggplot(
    mapping =
      aes(
        x = GeneRatio,
        y = Description,
        color = p.adjust
      )
  ) +
  geom_point() +
  scale_color_viridis_c(guide = guide_colourbar(barwidth = 10)) +
  scale_y_discrete(
    labels = stringr::str_wrap(
      paste(
        MEplotting[["ID"]],
        MEplotting[["Description"]],
        sep = ":"
      ),
      width = 40
    )
  ) +
  facet_wrap(
    facets = vars(module),
    scales = "free",
    ncol = 3,
    nrow = 3
  ) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x =
      element_text(
        size = 9
      ),
    axis.text.y =
      element_text(
        size = 9
      ),
    strip.text =
      element_text(
        size = 9
      )
  )

ht_opt$TITLE_PADDING = unit(c(6, 60), "points")

groupedComplexHeatmap(
  expr_mat            = module_scores_mat,
  gene_list           = c("M1.1", "M1.2", "M2.2", "M2.3", "M3.1", "M3.2", "M3.3", "M3.4", "M3.5", "M3.6", "M4.1", "M4.2", "M4.3", "M4.5", "M4.6", "M4.7", "M4.10", "M4.11", "M4.13", "M4.14", "M4.15", "M5.1", "M5.6", "M5.7", "M5.9", "M5.10", "M5.12", "M5.15", "M6.2", "M6.6", "M6.11", "M6.12", "M6.13", "M6.16", "M6.18", "M7.1", "M8.46", "M8.83", "M9.42", "ldg_a", "ldg_b", "ldg_1.1", "ldg_2.1", "mg"),
  md                  = study_md,
  annotation_palettes = group_pal2,
  row_grouping        = "responder",
  row_annotation      = c("responder", "RaceCode"),
  col_grouping        = "type",
  col_annotation      = dplyr::rename(module_annotation, obsv = 'module'),
  scale_exprs         = TRUE
)

groupedComplexHeatmap(
  expr_mat            = module_scores_mat,
  gene_list           = c("M1.1", "M1.2", "M2.1", "M2.2", "M2.3", "M3.1", "M3.2", "M3.3", "M3.4", "M3.5", "M3.6", "M4.1", "M4.2", "M4.3", "M4.4", "M4.5", "M4.6", "M4.7", "M4.8", "M4.9", "M4.10", "M4.11", "M4.12", "M4.13", "M4.14", "M4.15", "M4.16", "M5.1", "M5.2", "M5.3", "M5.4", "M5.5", "M5.6", "M5.7", "M5.8", "M5.9", "M5.10", "M5.11", "M5.12", "M5.13", "M5.14", "M5.15", "M6.1", "M6.2", "M6.3", "M6.4", "M6.5", "M6.6", "M6.7", "M6.8", "M6.9", "M6.10", "M6.11", "M6.12", "M6.13", "M6.14", "M6.15", "M6.16", "M6.17", "M6.18",
                          "M6.19", "M6.20", "M7.1", "M7.2", "M7.3", "M7.4", "M7.5", "M7.6", "M7.7", "M7.8", "M7.9", "M7.10", "M7.11", "M7.12", "M7.13", "M7.14", "M7.15", "M7.16", "M7.17", "M7.18", "M7.19", "M7.20", "M7.21", "M7.22", "M7.23", "M7.24", "M7.25", "M7.26", "M7.27", "M7.28", "M7.29", "M7.30", "M7.31", "M7.32", "M7.33", "M7.34", "M7.35", "M8.1", "M8.2", "M8.3", "M8.4", "M8.5", "M8.6", "M8.7", "M8.8", "M8.9", "M8.10", "M8.11", "M8.12", "M8.13", "M8.14", "M8.15", "M8.16", "M8.17", "M8.18", "M8.19", "M8.20", "M8.21",
                          "M8.22", "M8.23", "M8.24", "M8.25", "M8.26", "M8.27", "M8.28", "M8.29", "M8.30", "M8.31", "M8.32", "M8.33", "M8.34", "M8.35", "M8.36", "M8.37", "M8.38", "M8.39", "M8.40", "M8.41", "M8.42", "M8.43", "M8.44", "M8.45", "M8.46", "M8.47", "M8.48", "M8.49", "M8.50", "M8.51", "M8.52", "M8.53", "M8.54", "M8.55", "M8.56", "M8.57", "M8.58", "M8.59", "M8.60", "M8.61", "M8.62", "M8.63", "M8.64", "M8.65", "M8.66", "M8.67", "M8.68", "M8.69", "M8.70", "M8.71", "M8.72", "M8.73", "M8.74", "M8.75", "M8.76", "M8.77",
                          "M8.78", "M8.79", "M8.80", "M8.81", "M8.82", "M8.83", "M8.84", "M8.85", "M8.86", "M8.87", "M8.88", "M8.89", "M8.90", "M8.91", "M8.92", "M8.93", "M8.94", "M8.95", "M8.96", "M8.97", "M8.98", "M8.99", "M8.100", "M8.101", "M8.102", "M8.103", "M8.104", "M8.105", "M8.106", "M8.107", "M8.108", "M8.109", "M8.110", "M8.111", "M9.1", "M9.2", "M9.3", "M9.4", "M9.5", "M9.6", "M9.7", "M9.8", "M9.9", "M9.10", "M9.11", "M9.12", "M9.13", "M9.14", "M9.15", "M9.16", "M9.17", "M9.18", "M9.19", "M9.20", "M9.21", "M9.22",
                          "M9.23", "M9.24", "M9.25", "M9.26", "M9.27", "M9.28", "M9.29", "M9.30", "M9.31", "M9.32", "M9.33", "M9.34", "M9.35", "M9.36", "M9.37", "M9.38", "M9.39", "M9.40", "M9.41", "M9.42", "M9.43", "M9.44", "M9.45", "M9.46", "M9.47", "M9.48", "M9.49", "M9.50", "M9.51", "M9.52"),
  md                  = study_md,
  annotation_palettes = group_pal,
  row_grouping        = "responder",
  row_annotation      = c("responder", "RaceCode"),
  col_grouping        = "type",
  col_annotation      = dplyr::rename(module_annotation, obsv = 'module'),
  scale_exprs         = TRUE, col_fontsize = 6
)
