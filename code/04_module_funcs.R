#' @title extract_module_scores
#'
#' @description Extract the module scores from the
#' the metadata associated with a data object
#'
#' @param object A data object (either DESeqDataSet or DGEList)
#' @param ... a list of module names
#'
#' @return A tibble with the sample names and module values
#' @export
#'
extract_module_scores <- function(object, ...){
  UseMethod("extract_module_scores", object)
}

#' @rdname extract_module_scores
#' @method extract_module_scores DESeqDataSet
#' @importFrom purrr chuck
#' @importFrom rlang dots_list
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom tidyselect one_of
#' @return tibble
#'
extract_module_scores.DESeqDataSet <-
  function(
    object,
    ...
  ){
    module_names <- rlang::dots_list(...)
    object %>%
      purrr::chuck("colData") %>%
      tibble::as_tibble(
        rownames = "sample_name"
      ) %>%
      dplyr::select(
        sample_name,
        tidyselect::one_of(!!!module_names)
      )
  }

#' @rdname extract_module_scores
#' @method extract_module_scores DGEList
#' @importFrom purrr chuck
#' @importFrom rlang dots_list
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom tidyselect one_of
#' @return tibble
#'
extract_module_scores.DGEList <-
  function(
    object,
    ...
  ){
    module_names <- rlang::dots_list(...)
    object %>%
      purrr::chuck("samples") %>%
      tibble::as_tibble(
        rownames = "sample_name"
      ) %>%
      dplyr::select(
        sample_name,
        tidyselect::one_of(!!!module_names)
      )
  }

create_module_list <- function(module_file){
  readr::read_csv(module_file) %>%
    tidyr::nest(data = value) %>%
    tibble::deframe() %>%
    purrr::map(.f = dplyr::pull, value)
}


create_module_table <- function(
    ...,
    module_annotation
){
  module_list <- rlang::list2(...)

  module_tbl <-
    purrr::map_dfr(
      module_list,
      \(x) tibble::enframe(x) %>% tidyr::unnest(cols = "value")
    ) %>%
    dplyr::rename(
      module = name,
      gene   = value
    ) %>%
    dplyr::inner_join(module_annotation)

  module_tbl
}

extract_module_genes <- function(
    module_table,
    exprs_mat,
    module_annotation = "Interferon"
){
  dplyr::group_by(
    .data = module_table,
    gene
  ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      type == module_annotation,
      gene %in% rownames(exprs_mat)
    ) %>%
    dplyr::select(module, gene) %>%
    dplyr::arrange(module) %>%
    tibble::column_to_rownames("gene")
}

prior_calculate_module_score <- function(
    expr_mat,
    gene_list,
    score_func = c("tirosh", "rsvd", "irlba", "pca", "svd", "nmf"),
    absolute = FALSE
){
  score_func = match.arg(score_func)

  calc_score <-
    switch(
      score_func,
      tirosh = tirosh_score_module(expr_obj = expr_mat, features = gene_list),
      rsvd = rsvd::rsvd(A = expr_mat[gene_list,], k = 1)
        |> magrittr::use_series("v"),
      irlba = irlba::irlba(A = expr_mat[gene_list,], right_only = TRUE, nv = 1)
        |> magrittr::use_series("v"),
      pca = irlba::prcomp_irlba(x = t(expr_mat[gene_list,]), n = 1)
        |> magrittr::use_series("x"),
      svd = svd(x = expr_mat[gene_list,], nv = 1)
        |> magrittr::use_series("v"),
      nmf = {
        library(NMF)
        NMF::nmf(x = expr_mat[gene_list,], rank = 1)@fit@H |> t()
      }
    ) %>%
    magrittr::set_rownames(colnames(expr_mat))


  if (isTRUE(absolute)) {
    abs(calc_score)
  } else {
    calc_score
  }
}

calculate_module_score <- function(
    expr_mat,
    gene_list,
    score_func = c("rsvd", "irlba", "pca", "svd", "rrpca", "tirosh", "nmf"),
    ...
    #absolute = FALSE
){
  score_func <- match.arg(score_func)
  calc_score <-
    switch(
      score_func,
      rsvd = {
        result <- rsvd::rsvd(A = expr_mat[gene_list,], k = 1)
        list(
          scores   = result[['v']] |>
            scale(center = TRUE, scale = TRUE) |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(colnames(expr_mat)),
          loadings = result[['u']] |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(gene_list)
        )
      },
      irlba = {
        result <- irlba::irlba(A = expr_mat[gene_list,], right_only = FALSE, nv = 1)
        list(
          scores   = result[['v']] |>
            scale(center = TRUE, scale = TRUE) |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(colnames(expr_mat)),
          loadings = result[['u']] |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(gene_list)
        )
      },
      pca = {
        result <- irlba::prcomp_irlba(x = t(expr_mat)[,gene_list], n = 1)
        list(
          scores = result[['x']] |>
            scale(center = TRUE, scale = TRUE) |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(colnames(expr_mat)),
          loadings = result[['rotation']] |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(gene_list)
        )
      },
      svd = {
        result <- svd(x = expr_mat[gene_list,], nv = 1, nu = 1)
        list(
          scores   = result[['v']] |>
            scale(center = TRUE, scale = TRUE) |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(colnames(expr_mat)),
          loadings = result[['u']] |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(gene_list)
        )
      },
      rrpca = {
        result <- rsvd::rrpca(A = t(expr_mat[gene_list,]))
        list(
          scores = result[["L"]][,1] |>
            as.matrix() |>
            scale(center = TRUE, scale = TRUE) |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(colnames(expr_mat)),
          loadings = result[["S"]][1,] |>
            as.matrix() |>
            magrittr::set_colnames("PC1") |>
            magrittr::set_rownames(gene_list)
        )
      },
      tirosh = {
        result <- tirosh_score_module(expr_obj = expr_mat, features = gene_list) |> rename(PC1 = "value")
          list(
            scores = result,
            loadings = list()
          )
      },
      nmf = {
        library(NMF)
        expr_mat <- expr_mat[gene_list,]
        if (any(expr_mat < 0)){
          expr_mat <- (expr_mat + abs(min(expr_mat)))
        }
        result <- NMF::nmf(x = expr_mat[gene_list,], rank = 1)@fit
        list(
          scores = result@H |> t() |> magrittr::set_colnames("PC1"),
          loadings = result@W |> magrittr::set_colnames("PC1")
        )
      }
    )

  # if (isTRUE(absolute)) {
  #   abs(calc_score)
  # } else {
  #   calc_score
  # }
  calc_score
}

#' @title geometric_normalize_vsc
#'
#' @description given a set of reference genes an a gene expression matrix,
#' perform a per-sample geometric normalization. See
#' doi:10.1186/gb-2002-3-7-research0034 for more info
#'
#' @param exprs a gene-by-sample numeric matrix
#' @param normalization_genes a character vector with two or more genes to be
#' used for calculating the normalization factor.
#'
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom dplyr filter group_by summarise left_join mutate select
#' @importFrom tidyr pivot_longer pivot_wider
geometric_normalize_vsc <- function(exprs, normalization_genes){
  normalization_factors <-
    exprs %>%
    tibble::as_tibble(rownames = "gene") %>%
    dplyr::filter(gene %in% normalization_genes) %>%
    tidyr::pivot_longer(
      -gene,
      names_to = "sample_name"
    ) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::summarise(
      geo_mean = exp(mean(log(value)))
    )

  exprs_long <-
    exprs %>%
    tibble::as_tibble(rownames = "gene") %>%
    tidyr::pivot_longer(
      -gene,
      names_to = "sample_name"
    )

  norm_exprs <-
    dplyr::left_join(
      exprs_long,
      normalization_factors
    ) %>%
    dplyr::mutate(norm_value = value / geo_mean) %>%
    dplyr::select(
      gene,
      sample_name,
      norm_value
    ) %>%
    tidyr::pivot_wider(
      id_cols = "gene",
      names_from = "sample_name",
      values_from = "norm_value"
    ) %>%
    tibble::column_to_rownames("gene")

  norm_exprs
}

# Adapted from Seurat::AddModuleScore, which in turn took it from Tirosh (2006)
# TODO: this should probably be moved to the {moduleScoreR} package
# this version differs from others you may find elsewhere in that it
# functions on one module at a time
tirosh_score_module <- function(
    expr_obj,
    features,
    breaks = 25,
    num_ctrls = 100
) {

  data_avg <- Matrix::rowMeans(x = expr_obj)
  data_avg <- data_avg[order(data_avg)]
  data_cut <-
    ggplot2::cut_number(
      x = data_avg + rnorm(n = length(data_avg)) / 1e30,
      n = num_ctrls,
      labels = FALSE,
      right = FALSE
    ) %>%
    rlang::set_names(names(data_avg))

  # use only the module genes that are present in our dataset
  features_use <- features[which(features %in% rownames(expr_obj))]

  # find controls for each module gene
  ctrl_use <-
    purrr::map(
      .x = features_use,
      .f = \(x) {
        names(
          x = sample(
            x = data_cut[which(data_cut == data_cut[x])],
            size = num_ctrls,
            replace = FALSE
          )
        )
      }
    ) %>%
    unlist() %>%
    unique()

  ctrl_scores <- Matrix::colMeans(x = expr_obj[ctrl_use, ])
  features_scores <-  Matrix::colMeans(expr_obj[features_use, , drop = FALSE])

  (features_scores - ctrl_scores) %>%
    tibble::enframe() %>%
    tibble::column_to_rownames(var = "name")
}
