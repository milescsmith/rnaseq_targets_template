#### BLAST50

# if the pipeline haven't already been run, run the following
targets::tar_make(names = all_of(!!c("vsc_exprs", "up_tables","down_tables", "study_md", "res", "group_pal", "ISGs", "module_ISGs")))

library(targets)
library(tidyverse)
library(magrittr)

# else
tar_load(vsc_exprs)
tar_load(up_tables)
tar_load(down_tables)
tar_load(group_pal)
tar_load(ISGs)
tar_load(module_ISGs)
tar_load(study_md)
tar_load(res)
source("~/workspace/analysis/rnaseq/blast_2021_10_14/code/98_palettes_funcs.R")

responder_up_genes <- up_tables$`responder - non_responder`$gene
top_5_up <-
  res$`responder_vs_non_responder` %>%
  top_n(5, log2FoldChange) %>%
  pull(gene)

responder_down_genes <- down_tables$`responder - non_responder`$gene
top_5_down <-
  res$`responder_vs_non_responder` %>%
  top_n(-5, log2FoldChange) %>%
  pull(gene)

only_approved_genes <-
  left_join(
    res$`responder_vs_non_responder`,
    set_colnames(
      HGNChelper::checkGeneSymbols(
        x = res$`responder_vs_non_responder`$gene,
        unmapped.as.na = TRUE
        ),
      c("gene", "approved", "suggested_symbol")
      )
    ) %>%
  filter(approved == TRUE)

top_5_approved_up <-
  only_approved_genes %>%
  top_n(5, log2FoldChange) %>%
  pull(gene)

top_5_approved_down <-
  only_approved_genes %>%
  top_n(-5, log2FoldChange) %>%
  pull(gene)


#### top degs ####
# perform PCA using only the top up- and down-regulated genes
# only need the first PC
all_vsc_exprs <- vsc_exprs
all_study_md <- study_md

study_md <- filter(study_md, responder != "control") |> mutate(across(where(is.factor), fct_drop))

vsc_exprs <- vsc_exprs[, study_md %>% pull(sample_name)]

deg_pca <-
  irlba::irlba(
    A = vsc_exprs[
      c(
        responder_up_genes,
        responder_down_genes
        ),
      ],
    nv=1
    )

# left singular vector should have gene loadings
gene_loadings <-
  magrittr::set_rownames(
    x = deg_pca$u,
    value = c(
      responder_up_genes,
      responder_down_genes
      )
    ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
pc1_score <-
  magrittr::set_rownames(
    x = deg_pca$v,
    value = colnames(vsc_exprs)
    ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(study_md)

pc1_stats <-
  pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |> rstatix::add_significance()

pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = pc1_stats) +
  ggpubr::theme_pubr()

#### upregulated ####

up_deg_pca <-
  irlba::irlba(
    A = vsc_exprs[
      responder_up_genes,
    ],
    nv=1
  )

# left singular vector should have gene loadings
up_gene_loadings <-
  magrittr::set_rownames(
    x = up_deg_pca$u,
    value = responder_up_genes
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
up_pc1_score <-
  magrittr::set_rownames(
    x = up_deg_pca$v,
    value = colnames(vsc_exprs)
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(study_md)

up_pc1_stats <-
  up_pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

up_pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = up_pc1_stats) +
  ggpubr::theme_pubr()

#### downregulated ####

down_deg_pca <-
  irlba::irlba(
    A = vsc_exprs[
      responder_down_genes,
    ],
    nv=1
  )

# left singular vector should have gene loadings
down_gene_loadings <-
  magrittr::set_rownames(
    x = down_deg_pca$u,
    value = responder_down_genes
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
down_pc1_score <-
  magrittr::set_rownames(
    x = down_deg_pca$v,
    value = colnames(vsc_exprs)
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(study_md)

down_pc1_stats <-
  down_pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

down_pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = down_pc1_stats) +
  ggpubr::theme_pubr()

#### top 5 down ####

top_5_down_deg_pca <-
  irlba::irlba(
    A = vsc_exprs[
      top_5_down,
    ],
    nv=1
  )

# left singular vector should have gene loadings
top_5_down_gene_loadings <-
  magrittr::set_rownames(
    x =top_5_down_deg_pca$u,
    value = top_5_down
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
top_5_down_pc1_score <-
  magrittr::set_rownames(
    x = top_5_down_deg_pca$v,
    value = colnames(vsc_exprs)
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(study_md)

top_5_down_pc1_stats <-
  top_5_down_pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

top_5_down_pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = top_5_down_pc1_stats) +
  ggpubr::theme_pubr()

#### top 5 up ####

top_5_up_deg_pca <-
  irlba::irlba(
    A = vsc_exprs[
      top_5_up,
    ],
    nv=1
  )

# left singular vector should have gene loadings
top_5_up_gene_loadings <-
  magrittr::set_rownames(
    x =top_5_up_deg_pca$u,
    value = top_5_up
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
top_5_up_pc1_score <-
  magrittr::set_rownames(
    x = top_5_up_deg_pca$v,
    value = colnames(vsc_exprs)
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(study_md)

top_5_up_pc1_stats <-
  top_5_up_pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

top_5_up_pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = top_5_up_pc1_stats) +
  ggpubr::theme_pubr()

#### top 5 approved down ####

top_5_approved_down_deg_pca <-
  irlba::irlba(
    A = vsc_exprs[
      top_5_approved_down,
    ],
    nv=1
  )

# left singular vector should have gene loadings
top_5_approved_down_gene_loadings <-
  magrittr::set_rownames(
    x =top_5_approved_down_deg_pca$u,
    value = top_5_approved_down
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
top_5_approved_down_pc1_score <-
  magrittr::set_rownames(
    x = top_5_approved_down_deg_pca$v,
    value = colnames(vsc_exprs)
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(study_md)

top_5_approved_down_pc1_stats <-
  top_5_approved_down_pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

top_5_approved_down_pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = top_5_approved_down_pc1_stats) +
  ggpubr::theme_pubr() +
  ggtitle("PC1 of top 5 HUGO approved responder downregulated genes")

#### top 5 approved up ####

top_5_approved_up_deg_pca <-
  irlba::irlba(
    A = vsc_exprs[
      top_5_approved_up,
    ],
    nv=1
  )

# left singular vector should have gene loadings
top_5_approved_up_gene_loadings <-
  magrittr::set_rownames(
    x =top_5_approved_up_deg_pca$u,
    value = top_5_approved_up
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
top_5_approved_up_pc1_score <-
  magrittr::set_rownames(
    x = top_5_approved_up_deg_pca$v,
    value = colnames(vsc_exprs)
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(study_md)

top_5_approved_up_pc1_stats <-
  top_5_approved_up_pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

top_5_approved_up_pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = top_5_approved_up_pc1_stats) +
  ggpubr::theme_pubr() +
  ggtitle("PC1 of top 5 HUGO approved responder upregulated genes")

top_5_approved_up
top_5_approved_down

#### with controls ####

#### top 5 approved down with controls ####

top_5_approved_down_deg_pca <-
  irlba::irlba(
    A = all_vsc_exprs[
      top_5_approved_down,
    ],
    nv=1
  )

# left singular vector should have gene loadings
top_5_approved_down_gene_loadings <-
  magrittr::set_rownames(
    x = top_5_approved_down_deg_pca$u,
    value = top_5_approved_down
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
top_5_approved_down_pc1_score <-
  magrittr::set_rownames(
    x = top_5_approved_down_deg_pca$v,
    value = colnames(all_vsc_exprs)
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(all_study_md)

top_5_approved_down_pc1_stats <-
  top_5_approved_down_pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

top_5_approved_down_pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = top_5_approved_down_pc1_stats) +
  ggpubr::theme_pubr() +
  ggtitle("PC1 of top 5 HUGO approved responder downregulated genes")

#### top 5 approved up with controls ####

top_5_approved_up_deg_pca <-
  irlba::irlba(
    A = all_vsc_exprs[
      top_5_approved_up,
    ],
    nv=1
  )

# left singular vector should have gene loadings
top_5_approved_up_gene_loadings <-
  magrittr::set_rownames(
    x = top_5_approved_up_deg_pca$u,
    value = top_5_approved_up
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="gene")

# right singular vector has module score
top_5_approved_up_pc1_score <-
  magrittr::set_rownames(
    x = top_5_approved_up_deg_pca$v,
    value = colnames(all_vsc_exprs)
  ) |>
  magrittr::set_colnames(value="PC1") |>
  tibble::as_tibble(rownames="sample_name") |>
  dplyr::left_join(all_study_md)

top_5_approved_up_pc1_stats <-
  top_5_approved_up_pc1_score %>%
  rstatix::wilcox_test(formula = PC1 ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

top_5_approved_up_pc1_score %>%
  ggpubr::ggboxplot(
    x = "responder",
    fill = "responder",
    y = "PC1",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = top_5_approved_up_pc1_stats) +
  ggpubr::theme_pubr() +
  ggtitle("PC1 of top 5 HUGO approved responder upregulated genes")

vsc_controls <- all_vsc_exprs[,filter(all_study_md, responder == "control") |> pull(sample_name)]

controls_top5_up_exprs <- vsc_controls[top_5_approved_up,]
controls_top5_down_exprs <- vsc_controls[top_5_approved_down,]

controls_top5_up_exprs_pivot <-
  controls_top5_up_exprs |>
  as_tibble(rownames="gene") |>
  pivot_longer(
    -gene,
    names_to = "sample_name",
    values_to = "exprs"
    ) |>
  left_join(
    select(
      all_study_md,
      sample_name,
      responder
    )
  )

controls_top5_up_exprs_pivot |>
  group_by(gene, responder) |>
  summarise(max = max(exprs), mean = mean(exprs), std = sd(exprs))

all_top5_up_exprs <- all_vsc_exprs[top_5_approved_up,]
all_top5_down_exprs <- all_vsc_exprs[top_5_approved_down,]

top5_up_exprs_pivot <-
  all_top5_up_exprs |>
  as_tibble(rownames="gene") |>
  pivot_longer(
    -gene,
    names_to = "sample_name",
    values_to = "exprs"
  ) |>
  left_join(
    select(
      all_study_md,
      sample_name,
      responder
    )
  )

top5_up_exprs_pivot |>
  group_by(
    gene,
    responder
    ) |>
  summarise(
    max = max(exprs),
    mean = mean(exprs),
    std = sd(exprs)
    ) |>
  ggplot(
    aes(
      x=responder,
      y=mean,
      fill=responder
      ), color = "black"
    ) +
  geom_col() +
  facet_wrap(vars(gene)) +
  scale_fill_brewer(palette = "Set1") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

top5_down_exprs_pivot <-
  all_top5_down_exprs |>
  as_tibble(rownames="gene") |>
  pivot_longer(
    -gene,
    names_to = "sample_name",
    values_to = "exprs"
  ) |>
  left_join(
    select(
      all_study_md,
      sample_name,
      responder
    )
  )

top5_down_exprs_stats <-
  top5_down_exprs_pivot %>%
  group_by(gene) %>%
  rstatix::wilcox_test(formula = exprs ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

top5_down_exprs_pivot |>
  # group_by(
  #   gene,
  #   responder
  # ) |>
  # summarise(
  #   max = max(exprs),
  #   mean = mean(exprs),
  #   std = sd(exprs)
  # ) |>
  ggpubr::ggboxplot(
      x="responder",
      y="exprs",
      fill="responder",
      palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = top5_down_exprs_stats) +
  facet_wrap(vars(gene)) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  ggtitle("Top 5 down-regulated in responders vs non-responders")

top5_up_exprs_stats <-
  top5_up_exprs_pivot %>%
  group_by(gene) %>%
  rstatix::wilcox_test(formula = exprs ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_xy_position() |>
  rstatix::add_significance()

top5_up_exprs_pivot |>
  # group_by(
  #   gene,
  #   responder
  # ) |>
  # summarise(
  #   max = max(exprs),
  #   mean = mean(exprs),
  #   std = sd(exprs)
  # ) |>
  ggpubr::ggboxplot(
    x="responder",
    y="exprs",
    fill="responder",
    palette = "Set1"
  ) +
  ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(data = top5_up_exprs_stats) +
  facet_wrap(vars(gene)) +
  scale_fill_brewer(palette = "Set1") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  ggtitle("Top 5 up-regulated in responders vs non-responders")
