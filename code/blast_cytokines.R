#### import cytokines ####
cytokines = c(
  "April_TNFSF13", "BAFF_BLyS_TNFSF13B",
  "betaNGF", "BMP2", "CCL2_JE_MCP1", "CCL3_MIP1alpha", "CCL4_MIP1beta",
  "CXCL9_MIG", "CXCL10_IP10", "CCL11_Eotaxin", "CD25_IL2Ralpha", "CD40Ligand_TNFSF5",
  "ESelectin_CD62E", "Fas_TNFRSF6_CD95", "FasLigand_TNFSF6", "GCSF", "ICAM_CD54",
  "IFNalpha", "IFNbeta", "IFNgamma", "IL1ra_IL1F3", "IL1alpha_IL1F1", "IL1beta_IL1F2",
  "IL2", "IL4", "IL5", "IL6", "IL7", "IL8_CXCL8", "IL10", "IL12p70", "IL17_IL17A", "IL21",
  "IL15", "IL13", "IL23", "IL33", "Leptin_OB", "LIF", "PDGFBB", "SCF_ckitLigand",
  "Syndecan1_CD138", "Resistin", "TNFalpha", "TNFRII_TNFRSF1B", "TNFRI_TNFRSF1A",
  "TRAIL_TNFSF10", "VCAM1_CD106", "VEGF"
)

cytokine_dataset <-
  readxl::read_excel(
    path = "metadata/BLAST analysis dataset.xlsx",
    sheet="main",
    skip=1
    ) %>%
  select(
    NovaSeq_Sample_ID, Cytokine_Sample_ID, Responder, one_of(cytokines)
  ) %>%
  mutate(
    NovaSeq_Sample_ID = janitor::make_clean_names(NovaSeq_Sample_ID, case="screaming_snake"),
    Responder =
      dplyr::recode(
        .x = Responder,
        `1` = "responder",
        `0` = "non-responder"
        ) %>%
      forcats::as_factor()
    ) %>%
  pivot_longer(
    cols = one_of(cytokines),
    names_to="cytokine",
    values_to = "exprs"
    )

cytokine_dataset_stats <-
  cytokine_dataset %>%
  filter(
    str_detect(string=Cytokine_Sample_ID, pattern="BSL$"),
    !is.na(Responder),
    !is.na(exprs),
    cytokine != "LIF"
    ) %>%
  group_by(cytokine) %>%
  rstatix::wilcox_test(formula = exprs ~ Responder) |>
  rstatix::add_xy_position(group = "Responder", scales = "free_y") |>
  rstatix::add_significance()

cytokine_dataset %>%
  filter(
    str_detect(string=Cytokine_Sample_ID, pattern="BSL$"),
    !is.na(Responder),
    !is.na(exprs),
    cytokine != "LIF",
    cytokine %in% pull(filter(cytokine_dataset_stats, p.signif != "ns"), cytokine)
  ) %>%
  ggpubr::ggboxplot(
    x="Responder",
    y="exprs",
    fill="Responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm(groupOnX = TRUE) +
  ggpubr::stat_pvalue_manual(data = filter(cytokine_dataset_stats, p.signif != "ns") ) +
  facet_wrap(vars(cytokine), scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
