readxl::excel_sheets(here::here("metadata/BLAST_200312.xlsx"))
clinical_md <-
  readxl::read_excel(
    path = here::here("metadata/BLAST_200312.xlsx"),
    sheet = "clinical",
    .name_repair = janitor::make_clean_names
    )

factor_columns <-
  c(
    "sex",
    "race_code"
  )

bool_columns <-
  c(
    "apls",
    "anti_phospholipid_history",
    "thrombosis",
    "pregnancy",
    "thrombocytopenia_itp",
    "nephritis",
    "ra",
    "sjogrens",
    "uctd",
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
    "ana"
  )

symptoms_columns = c(
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
)

score_columns = c(
  "selena",
  "sledai",
  "bilag",
  "acr"
)

targets::tar_load(bl_ctrl_md)

clinical_md_tbl <-
  clinical_md |>
  dplyr::filter(
    visit_ref %in% bl_ctrl_md[["visit_ref"]]
  ) |>
  dplyr::select(
    # -responder_1_non_responder_0,
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
      # race_code == "[AI]" ~ "AI",
      # race_code == "[A]" ~ "A",
      # race_code == "[H]" ~ "H",
      TRUE ~ "other"
    ),
    dplyr::across(
      .cols = tidyselect::all_of(score_columns),
      .fns = as.numeric
    )
  )

clinical_md_tbl |>
  dplyr::group_by(responder_status) |>
  dplyr::mutate(
    responder_status =
      stringr::str_replace(
        string = responder_status,
        pattern = "_",
        replacement = "-"
        ),
    responder_status = stringr::str_glue("{snakecase::to_any_case(string = as.character(responder_status), 'title', sep_out='-')} (N={n()})")
  ) |>
  dplyr::summarise(
    across(
      one_of(bool_columns),
      .fns = sum,
      .names = "{.col}_count"
    ),
    across(
      one_of(bool_columns),
      .fns = \(x) sum(x)/n(),
      .names = "{.col}_pct"
    ),
    median_age = median(age),
    age_iqr = matrixStats::iqr(age),
    aa_count = sum(race_code == "AA"),
    ea_count = sum(race_code == "EA"),
    non_aa_ea_count = sum(race_code == "other"),
    aa_pct = aa_count/n(),
    ea_pct = ea_count/n(),
    non_aa_ea_pct = non_aa_ea_count/n(),
    female_count = sum(sex == "Female"),
    male_count = sum(sex == "Male"),
    female_pct = female_count/n(),
    male_pct = male_count/n(),
    dplyr::across(
      .cols = tidyselect::one_of(score_columns),
      .fns = mean,
      .names = "mean_{.col}"
    ),
    dplyr::across(
      .cols = tidyselect::one_of(score_columns),
      .fns = max,
      .names = "max_{.col}"
    ),
    dplyr::across(
      .cols = tidyselect::one_of(score_columns),
      .fns = min,
      .names = "min_{.col}"
    ),
    dplyr::across(
      .cols = tidyselect::one_of(score_columns),
      .fns = matrixStats::iqr,
      .names = "iqr_{.col}"
    )
  ) |>
  dplyr::transmute(
    responder_status                        = responder_status,
    `Sex, n (%)`                         = c("","",""),
    ` Female`                            = stringr::str_glue("{female_count} ({round(female_pct, digits = 3)*100})"),
    ` Male`                              = stringr::str_glue("{male_count} ({round(male_pct, digits = 3)*100})"),
    `Race/Ethnicity, n (%)`              = c("","", ""),
    ` White`                             = stringr::str_glue("{ea_count} ({round(ea_pct, digits = 3)*100})"),
    ` Black`                             = stringr::str_glue("{aa_count} ({round(aa_pct, digits = 3)*100})"),
    ` Other`                             = stringr::str_glue("{non_aa_ea_count} ({round(non_aa_ea_pct, digits = 3)*100})"),
    `Age at visit, median (IQR)`         = stringr::str_glue("{round(median_age, digits = 0)} ({round(median_age-age_iqr, 0)}-{round(median_age+age_iqr,0)})"),
    `Total ACR score` = stringr::str_glue("{round(mean_acr, digits = 0)} ({max_acr}-{min_acr})"),
    `ACR SLE criterion, n (%)`           = c("", "", ""),
    ` Malar Rash`                        = stringr::str_glue("{malar_rash_count} ({round(malar_rash_pct, 3)*100})"),
    ` Discoid Rash`                      = stringr::str_glue("{discoid_rash_count} ({round(discoid_rash_pct, 3)*100})"),
    ` Photosensitivity`                  = stringr::str_glue("{photosensitivity_count} ({round(photosensitivity_pct, 3)*100})"),
    ` Oral Ulcers`                       = stringr::str_glue("{oral_ulcers_count} ({round(oral_ulcers_pct, 3)*100})"),
    ` Renal Disorder`                    = stringr::str_glue("{renal_disorder_count} ({round(renal_disorder_pct, 3)*100})"),
    ` Hematologic Disorder`              = stringr::str_glue("{hematologic_disorder_count} ({round(hematologic_disorder_pct, 3)*100})"),
    ` Immunologic Disorder`              = stringr::str_glue("{immunologic_disorder_count} ({round(immunologic_disorder_pct, 3)*100})"),
    ` Antinuclear Antibody`              = stringr::str_glue("{ana_count} ({round(ana_pct, 3)*100})"),
    ` Antiphospholipid Syndrome`         = stringr::str_glue("{anti_phospholipid_history_count} ({round(anti_phospholipid_history_pct, 3)*100})"),
    ` Thrombosis` = stringr::str_glue("{thrombosis_count} ({round(thrombosis_pct, 3)*100})"),
    ` Pregnancy` = stringr::str_glue("{pregnancy_count} ({round(pregnancy_pct, 3)*100})"),
    ` Thrombocytopenia` = stringr::str_glue("{thrombocytopenia_itp_count} ({round(thrombocytopenia_itp_pct, 3)*100})"),
    ` Nephritis` = stringr::str_glue("{nephritis_count} ({round(nephritis_pct, 3)*100})"),
    `Physician Global Assessment, mean (IQR)` = stringr::str_glue("{round(mean_selena, digits = 0)} ({round(mean_selena-iqr_selena, 0)}-{round(mean_selena+iqr_selena,0)})"),
    `SLEDAI, mean (IQR)` = stringr::str_glue("{round(mean_sledai, digits = 0)} ({round(mean_sledai-iqr_sledai, 0)}-{round(mean_sledai+iqr_sledai,0)})"),
    `BILAG index, mean (IQR)` = stringr::str_glue("{round(mean_bilag, digits = 0)} ({round(mean_bilag-iqr_bilag, 0)}-{round(mean_bilag+iqr_bilag,0)})")
  ) |>
  dplyr::mutate(dplyr::across(tidyselect::everything(), as.character)) |>
  tidyr::pivot_longer(-responder_status, names_to = "Characteristic", values_transform = as.character) |>
  tidyr::pivot_wider(names_from = "responder_status") |>
  dplyr::mutate(
    across(
      .cols = everything(),
      .fns = \(x) if_else(condition = stringr::str_detect(string = x, pattern = "NA"), true = "", false = x)
    )
  ) |>
  flextable::qflextable() |>
  flextable::align(align = "left", part = "body", j = 1) |>
  flextable::align(align = "left", part = "header", j = 1) |>
  flextable::bold(i = 1, part = "header") |>
  flextable::footnote(
    i = 7,
    j = 1,
    value =
      flextable::as_paragraph(
        "Other includes 1 Hispanic and 1 American Indian Non-Responder; 1 Hispanic Responder; and 2 American Indian, 1 Asian, 2 Hispanic, 1 Mixed ethnicity healthy controls"
      ),
    ref_symbols = c("a"),
    part = "body"
  ) |>
  flextable::add_footer_lines(
    value =
      flextable::as_paragraph(
        "Abbreviations: ACR – American College of Rheumatology, SLEDAI – Systemic Lupus Erythematosus Disease Activity Index, BILAG - British Isles Lupus Assessment Group, IQR - interquartile range"
      )
  ) |>
  flextable::padding(
    j = 1,
    i = c(
      2,3,
      5,6,7,
      11:23
      ),
    padding.left = 15
  )

