# This is a newer version (better?), but it isn't what was used in the CLE paper
# ident_clusters <- function(
#   expr_mat,
#   sig_pc_method = c("elbow", "horn"),
#   bootmethod = "bojit",
#   max_k = 10,
#   ...
# ){
#   sig_pc_method = match.arg(sig_pc_method)
#
#   expr_mat <- expr_mat[which(apply(expr_mat, 1, var) != 0),]
#
#   message("Performing PCA on expression matrix")
#   pca_res <-
#     PCAtools::pca(
#       mat    = expr_mat,
#       center = TRUE,
#       scale  = FALSE,
#       removeVar = 0.1
#     )
#
#   pcs_use <-
#     switch(
#       sig_pc_method,
#       elbow = PCAtools::findElbowPoint(pca_res$variance),
#       horn = PCAtools::parallelPCA(mat = expr_mat)[["n"]]
#     )
#
#   if (max_k >= nrow(pca_res[["rotated"]])){
#     message("max_k was set to a value higher than the ",
#             "number of samples, while k *must* be less ",
#             "than the number of samples minus 1. ",
#             "Adjusting max_k...")
#     max_k <- nrow(pca_res[["rotated"]])-1
#   }
#
#   cbs <- future.apply::future_lapply(
#     X = seq(2,max_k),
#     FUN = \(y) {
#       fpc::clusterboot(
#         data          = pca_res[["rotated"]][,1:pcs_use],
#         clustermethod = fpc::claraCBI,
#         k             = y,
#         bootmethod    = "bojit"
#       )
#     }
#   )
#
#   # TODO: handle situations where there are no
#   # identifiable stable clusters
#   optimal_k <-
#     which(
#       purrr::imap_lgl(
#         .x = cbs,
#         .f = \(x, y) {
#           all(cbs[[y]][["bojitmean"]] > 0.6)
#           }
#         )
#       )
#
#     if (length(optimal_k) == 0){
#       sample_clusters <-
#         tibble::tibble(
#           sample_name = colnames(expr_mat),
#           cluster     = 1
#         )
#     } else {
#       sample_clusters <-
#         cbs[[max(optimal_k)]][["partition"]] |>
#         rlang::set_names(colnames(expr_mat)) |>
#         tibble::enframe(
#           name  = "sample_name",
#           value = "cluster"
#         ) |>
#         dplyr::mutate(
#           cluster = forcats::as_factor(cluster)
#         )
#     }
#
#   ret_values <-
#     list(
#       kmeans_res     = cbs,
#       k              = max(optimal_k) + 1,
#       clusters       = sample_clusters
#     )
#
#   ret_values
# }

ident_clusters <- function(
    expr_mat,
    optimal_k_method = "Tibs2001SEmax",
    nstart = 25,
    max_k = 50,
    B = 100,
    d_power = 2
    ){
  module_rf <-
    randomForest::randomForest(
      x = expr_mat,
      y = NULL,
      prox = T
      )

  rf_distance_mat <-
    stats::dist(1 - module_rf[["proximity"]]) %>%
    as.matrix()

  library(future)
  library(matrixStats)
  library(FasterMatrixMath)
  kmeans_gap_stat <-
    cluster::clusGap(
      x = rf_distance_mat,
      FUNcluster = kmeans,
      nstart = nstart,
      K.max = max_k,
      B = B,
      d.power = d_power,
      future_plan = "multisession",
      parallel = TRUE
    )

  new_optimal_k <-
    with(
      data = kmeans_gap_stat,
      expr = cluster::maxSE(
        Tab[,"gap"],
        Tab[,"SE.sim"],
        method=optimal_k_method
      )
    )

  k_clusters <-
    stats::kmeans(
      x = rf_distance_mat,
      centers = new_optimal_k,
      nstart = 25
    )

  sample_clusters <-
    tibble::enframe(x = k_clusters[["cluster"]],
            name = "sample_name",
            value = "cluster")

  list(
    kmeans_res = k_clusters,
    rf_distance = rf_distance_mat,
    clusters = sample_clusters,
    gap_stat = kmeans_gap_stat
  )
}

rf_based_ident_clusters <- function(
    expr_mat,
    optimal_k_method = "Tibs2001SEmax",
    nstart = 25,
    max_k = 50,
    B = 100,
    d_power = 2
    ){
  module_rf <-
    randomForest::randomForest(
      x = expr_mat,
      y = NULL,
      prox = T
      )

  rf_distance_mat <-
    stats::dist(1 - module_rf[["proximity"]]) %>%
    as.matrix()

  library(future)
  library(matrixStats)
  library(FasterMatrixMath)
  kmeans_gap_stat <-
    cluster::clusGap(
      x = rf_distance_mat,
      FUNcluster = kmeans,
      nstart = nstart,
      K.max = max_k,
      B = B,
      d.power = d_power,
      future_plan = "multisession",
      parallel = TRUE
    )

  new_optimal_k <-
    with(
      data = kmeans_gap_stat,
      expr = cluster::maxSE(
        Tab[,"gap"],
        Tab[,"SE.sim"],
        method=optimal_k_method
      )
    )

  k_clusters <-
    stats::kmeans(
      x = rf_distance_mat,
      centers = new_optimal_k,
      nstart = 25
    )

  sample_clusters <-
    tibble::enframe(x = k_clusters[["cluster"]],
            name = "sample_name",
            value = "cluster")

  list(
    kmeans_res = k_clusters,
    rf_distance = rf_distance_mat,
    clusters = sample_clusters,
    gap_stat = kmeans_gap_stat
  )
}
