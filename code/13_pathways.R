plot_gex_violins <- function(
  dataset,
  goi,
  stats_data,
  x_var,
  y_var,
  palette = "ggsci::uniform_startrek"
){

  ggpubr::ggviolin(
    data =
      dplyr::filter(
        .data = dataset,
        gene == goi
        ),
    x = x_var,
    y = y_var,
    fill = x_var,
    draw_quantiles = c(0.25, 0.5, 0.75),
    trim = TRUE,
    add.params = list(
      size = 0.75,
      alpha = 0.5
      )
    ) +
    ggbeeswarm::geom_quasirandom(size = 1) +
    paletteer::scale_fill_paletteer_d(palette) +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(gene),
      scales = "free_y"
    ) +
    ggpubr::stat_pvalue_manual(
      data =
        dplyr::filter(
          .data = stats_data,
          gene == goi
          ),
      label = "p.adj.signif",
      hide.ns = TRUE
    ) +
    ggplot2::labs(
      title = goi,
      caption =
        stringr::str_glue(
          "Mann - Whitney, Benjamini - Hochberg adjusted value.",
          "<br> * <i>p</i> <0.05, ",
          "** <i>p</i> <0.01, ",
          "*** <i>p</i> <0.001, ",
          "**** <i>p</i> <0.0001"
        )
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3)) +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
          ),
      plot.caption = ggplot2::element_markdown()
      )
  }


extract_pathway_exprs <- function(
    exprs_mat,
    gene_list,
    metadata
    ){
  tibble::as_tibble(
    x = exprs_mat,
    rownames = "gene"
    ) |>
  dplyr::filter(gene %in% gene_list) |>
  tidyr::pivot_longer(
      -gene,
      names_to = "sample_name",
      values_to = "expression"
      ) |>
  dplyr::left_join(
    y =
      dplyr::select(
        .data = metadata,
        sample_name,
        study_group,
        grant_defined_severity,
        project_group
        )
    )
}


calc_gex_stats <- function(expr_tbl){
  expr_tbl |>
  dplyr::group_by(gene) |>
  rstatix::wilcox_test(
    formula         = expression ~ study_group,
    p.adjust.method = "BH"
  ) |>
  rstatix::add_significance() |>
  grouped_add_xy_positions(
    stats_tbl          = .,
    data_tbl           = expr_tbl,
    group_var          = gene,
    compare_value      = expression,
    percent_shift_down = 0.99
  )
}

grouped_complex_heatmap <- function(
  exprs_mat,
  exprs_stats,
  sample_info_tbl,
  pathway,
  gene_list,
  group_var,
  group1,
  group2,
  group1_color = "#440154FF",
  group2_color = "#FDE725FF",
  row_name_fontsize = 12
){

  if(is.character(group_var)){
    diffused_group_var <- rlang::sym(group_var)
  }

  diffused_group_var <- rlang::enquo(diffused_group_var)

  ha_vec <-
    exprs_stats |>
    purrr::pluck(pathway) |>
    dplyr::select(
      gene,
      p.adj.signif
    )

  ha_stats <-
    gene_list |>
    purrr::pluck(pathway) |>
    tibble::enframe() |>
    dplyr::rename(gene = value) |>
    dplyr::filter(
      gene %in% exprs_mat[[pathway]][["gene"]],
      !gene %in% ha_vec[["gene"]]
    ) |>
    dplyr::select(-name) |>
    dplyr::mutate(p.adj.signif = "ns") |>
    dplyr::bind_rows(ha_vec) |>
    tibble::column_to_rownames("gene")

  group1_pathway_exprs_mat <-
    exprs_mat |>
    purrr::pluck(pathway) |>
    dplyr::select(
      gene,
      sample_name,
      expression,
    ) |>
    tidyr::pivot_wider(
      names_from = "gene",
      values_from = "expression"
    ) |>
    tibble::column_to_rownames("sample_name") |>
    as.matrix() |>
    scale() |>
    magrittr::extract(
      filter(
        .data =
          purrr::pluck(
            .x = exprs_mat,
            pathway
          ),
        {{diffused_group_var}} == group1
      ) |>
        dplyr::pull(sample_name) |>
        unique(),
    ) |>
    t()

  group2_pathway_exprs_mat <-
    exprs_mat |>
    purrr::pluck(pathway) |>
    dplyr::select(
      gene,
      sample_name,
      expression,
    ) |>
    tidyr::pivot_wider(
      names_from = "gene",
      values_from = "expression"
    ) |>
    tibble::column_to_rownames("sample_name") |>
    as.matrix() |>
    scale() |>
    magrittr::extract(
      dplyr::filter(
        .data =
          purrr::pluck(
            .x = exprs_mat,
            pathway
          ),
        {{diffused_group_var}} == group2
        # study_group == group2
      ) |>
        dplyr::pull(sample_name) |>
        unique(),
    ) |>
    t()

  group2_heatmap <-
    ComplexHeatmap::Heatmap(
      matrix            = group2_pathway_exprs_mat,
      right_annotation  =
        ComplexHeatmap::rowAnnotation(
          padj = ComplexHeatmap::anno_text(
            ha_stats |>
              tibble::as_tibble(rownames = "gene") |>
              dplyr::arrange(gene) |>
              tibble::deframe(),
            gp = grid::gpar(fontsize = row_name_fontsize)
            )
          ),
      top_annotation    =
        ComplexHeatmap::columnAnnotation(
          group_var = ComplexHeatmap::anno_block(
            gp = grid::gpar(fill=group2_color)
            )
          ),
      cluster_rows      = FALSE,
      cluster_columns   = TRUE,
      show_row_names    = TRUE,
      show_column_names = FALSE,
      column_split      =
        dplyr::filter(
          .data = sample_info_tbl,
          study_group == group2
        ) |>
        tibble::column_to_rownames("sample_name"),
      col               =
        circlize::colorRamp2(
          seq(-4,4),
          viridis::viridis(length(seq(-4,4)))
          ),
      name              = pathway,
      row_names_gp      = grid::gpar(fontsize = row_name_fontsize)
    )

  group1_heatmap <-
    ComplexHeatmap::Heatmap(
      matrix            = group1_pathway_exprs_mat,
      top_annotation    =
        ComplexHeatmap::columnAnnotation(
          group_var = ComplexHeatmap::anno_block(
            gp = grid::gpar(fill=group1_color)
            )
          ),
      cluster_rows      = FALSE,
      cluster_columns   = TRUE,
      show_row_names    = TRUE,
      show_column_names = FALSE,
      column_split      =
        dplyr::filter(
          .data = sample_info_tbl,
          {{diffused_group_var}} == group1
          # study_group == group1
        ) |>
        tibble::column_to_rownames("sample_name"),
      col               =
        circlize::colorRamp2(
          seq(-4,4),
          viridis::viridis(length(seq(-4,4)))
          ),
      name              = pathway,
      row_names_gp = gpar(fontsize = row_name_fontsize)
    )

  ht_list = group1_heatmap + group2_heatmap
  m1      = group1_pathway_exprs_mat
  m2      = group2_pathway_exprs_mat
  rg      = range(c(group1_pathway_exprs_mat, group2_pathway_exprs_mat))
  rg[1]   = rg[1] - (rg[2] - rg[1])* 0.02
  rg[2]   = rg[2] + (rg[2] - rg[1])* 0.02

  anno_multiple_boxplot = function(index) {
    nr = length(index)

    pushViewport(
      grid::viewport(
        xscale = rg,
        yscale = c(0.5, nr + 0.5)
      )
    )

    for(i in seq_along(index)) {
      grid::grid.rect(
        y = nr-i+1,
        height = 1,
        default.units = "native"
      )

      grid::grid.boxplot(
        m1[ index[i], ],
        pos = nr-i+1 + 0.2,
        box_width = 0.3,
        gp = grid::gpar(fill = group1_color),
        direction = "horizontal"
      )
      grid::grid.boxplot(
        m2[ index[i], ],
        pos = nr-i+1 - 0.2,
        box_width = 0.3,
        gp = grid::gpar(fill = group2_color),
        direction = "horizontal"
      )
    }
    grid::grid.xaxis()
    grid::popViewport()
  }

  ht_list =
    ComplexHeatmap:rowAnnotation(
      boxplot              = anno_multiple_boxplot,
      width                = unit(2, "cm"),
      show_annotation_name = FALSE
    ) +
    ht_list

  lgd <-
    ComplexHeatmap::Legend(
      labels    = c(group1, group2),
      title     = "Expression",
      legend_gp = grid::gpar(fill = c(group1_color, group2_color))
    )

  ComplexHeatmap::draw(
    object              = ht_list,
    padding             = unit(c(20, 2, 2, 2), "mm"),
    heatmap_legend_list = list(lgd)
  )
}

empty_enrichment <- function(er){
  sig_er <- dplyr::filter(er@result, p.adjust <= er@pvalueCutoff)
  if (nrow(sig_er) > 0){
    FALSE
  } else {
    TRUE
  }
}
