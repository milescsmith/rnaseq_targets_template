cytokine_names <- c("baff" = "BAFF",
                    "blc_cxcl13" = "BCL-CXCL13",
                    "gm_csf" = "GM-CSF",
                    "gro_alpha" = "GROalpha",
                    "ifn_alpha" = "IFNalpha",
                    "ifn_gamma" = "IFNgamma",
                    "il1alpha" = "IL1alpha",
                    "il1beta" = "IL1beta",
                    "il_1ra" = "IL1RA",
                    "il_10" = "IL10",
                    "il_12p70" = "IL12p70",
                    "il_13" = "IL13",
                    "il_15" = "IL15",
                    "il_17a" = "IL17a",
                    "il_18" = "IL18",
                    "il_2" = "IL2",
                    "il_21" = "IL21",
                    "il22" = "IL22",
                    "il_23" = "IL23",
                    "il_27" = "IL27",
                    "il_31" = "IL31",
                    "il_4" = "IL4",
                    "il5" = "IL5",
                    "il_7" = "IL7",
                    "il_8" = "IL8",
                    "il_9" = "IL9",
                    "ip_10" = "IP10",
                    "mcp1" = "MCP1",
                    "mip1alpha" = "MIP1alpha",
                    "mip_1beta" = "MIP1beta",
                    "rantes" = "RANTES",
                    "s_cd40l" = "sCD40L",
                    "scf" = "SCF",
                    "sdf_1alpha" = "SDF1alpha",
                    "s_icam_1" = "sICAM1",
                    "tnf_alpha" = "TNFalpha",
                    "tnf_beta" = "TNFbeta")

antibody_columns <-
  c(
    "ana_flag",
    "ds_dna", "ds_dna_flag",
    "chromatin", "chromatin_flag",
    "ribosomal_p", "ribosomal_p_flag",
    "ssa_ro_composite","ssa_ro_composite_flag",
    "ssb_la", "ssb_la_flag",
    "sm","sm_flag",
    "sm_rnp", "sm_rnp_flag",
    "rnp_composite", "rnp_composite_flag",
    "centromere_b", "centromere_b_flag",
    "scl_70", "scl_70_flag",
    "jo_1", "jo_1_flag"
  )

acr_sle_criteria <-
  c(
    "cluster",
    "discoid_rash",
    "malar_rash",
    "ana",
    "photo_sensitivity",
    "blood_disorder",
    "oral_ulcers",
    "immune_disorder"
  )

categorical_characteristics <-
  c(
    "female",
    "male",
    "african_american",
    "other",
    "smoking_now",
    "smoking_then",
    "non_smoker",
    "chronic_cle",
    "scle",
    "dle",
    "lp",
    "tle",
    "no_meds",
    "topicals",
    "antimalarials"
  )

continuous_characteristics <-
  c(
    "age_at_visit",
    "clasi_activity",
    "clasi_damage"
  )
