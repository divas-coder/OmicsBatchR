#' Benchmark batch correction methods using LISI and kBET scores
#'
#' @description
#' `benchmark_methods()` computes integration quality metrics for each
#' corrected Seurat object in an `OmicsBatchR` result list.
#'
#' Two metrics are computed:
#'
#' - **iLISI** (Integration LISI): measures how well batches are mixed.
#'   Higher = better mixing (range: 1 to n_batches).
#' - **cLISI** (Cell-type LISI): measures how well cell types are preserved.
#'   Lower = better separation (range: 1 to n_celltypes).
#' - **kBET acceptance rate**: measures local batch mixing via k-nearest
#'   neighbours. Higher = better (range: 0–1).
#'
#' @param results An `OmicsBatchR` list from `run_batch_correction()`.
#' @param batch_key Character. Batch column name in metadata.
#' @param celltype_key Character or NULL. Cell type column for cLISI.
#'   If NULL, cLISI is skipped. Default: `NULL`.
#' @param reduction_map Named character vector mapping method names to their
#'   reduction names. If NULL, auto-detected. Default: `NULL`.
#' @param k Integer. Number of neighbours for kBET. Default: `30`.
#' @param verbose Logical. Default: `TRUE`.
#'
#' @return A data frame with one row per method and columns:
#'   `method`, `iLISI`, `cLISI` (if celltype_key given), `kBET_acceptance`,
#'   `overall_score` (composite rank score).
#'
#' @examples
#' \dontrun{
#' results <- run_batch_correction(obj, batch_key = "batch")
#' scores  <- benchmark_methods(results, batch_key = "batch",
#'                              celltype_key = "cell_type")
#' print(scores)
#' }
#'
#' @export
benchmark_methods <- function(results,
                               batch_key     = "batch",
                               celltype_key  = NULL,
                               reduction_map = NULL,
                               k             = 30,
                               verbose       = TRUE) {

  if (!inherits(results, "OmicsBatchR")) {
    stop("'results' must be an OmicsBatchR object from run_batch_correction().")
  }

  # Default reduction names per method
  default_reductions <- c(
    harmony    = "harmony",
    seurat_cca = "pca_cca",
    combat     = "pca_combat",
    fastmnn    = "mnn",
    scanorama  = "scanorama",
    bbknn      = "bbknn_pca",
    scvi       = "scvi"
  )

  if (!is.null(reduction_map)) {
    default_reductions[names(reduction_map)] <- reduction_map
  }

  scores_list <- list()

  for (method in names(results)) {
    if (verbose) .msg(sprintf("Benchmarking: %s", toupper(method)))

    obj       <- results[[method]]
    reduction <- default_reductions[method]

    # Check reduction exists
    if (!reduction %in% names(obj@reductions)) {
      warning(sprintf("Reduction '%s' not found for method '%s'. Skipping.",
                      reduction, method))
      next
    }

    embedding <- Seurat::Embeddings(obj, reduction = reduction)
    batch_vec <- obj@meta.data[[batch_key]]

    # ── iLISI score ──────────────────────────────────────────
    ilisi <- tryCatch({
      lisi_scores <- .compute_lisi(embedding, batch_vec, k = k)
      mean(lisi_scores, na.rm = TRUE)
    }, error = function(e) NA_real_)

    # ── cLISI score ──────────────────────────────────────────
    clisi <- NA_real_
    if (!is.null(celltype_key) && celltype_key %in% colnames(obj@meta.data)) {
      celltype_vec <- obj@meta.data[[celltype_key]]
      clisi <- tryCatch({
        lisi_scores <- .compute_lisi(embedding, celltype_vec, k = k)
        mean(lisi_scores, na.rm = TRUE)
      }, error = function(e) NA_real_)
    }

    # ── kBET acceptance rate ──────────────────────────────────
    kbet_score <- tryCatch({
      .compute_kbet(embedding, batch_vec, k = k)
    }, error = function(e) NA_real_)

    scores_list[[method]] <- data.frame(
      method          = method,
      iLISI           = round(ilisi, 4),
      cLISI           = round(clisi, 4),
      kBET_acceptance = round(kbet_score, 4),
      stringsAsFactors = FALSE
    )
  }

  scores_df <- do.call(rbind, scores_list)
  rownames(scores_df) <- NULL

  # ── Compute composite overall score ──────────────────────────
  scores_df <- .compute_overall_score(scores_df)

  # Sort by overall score descending
  scores_df <- scores_df[order(scores_df$overall_score, decreasing = TRUE), ]

  class(scores_df) <- c("OmicsBatchR_benchmark", "data.frame")

  if (verbose) {
    .success("Benchmarking complete!")
    cat("\n")
    print(scores_df, digits = 3, row.names = FALSE)
    cat(sprintf("\nBest method: %s (score = %.3f)\n",
                scores_df$method[1], scores_df$overall_score[1]))
  }

  return(scores_df)
}


#' Automatically select the best batch correction method
#'
#' @description
#' `select_best_method()` runs benchmarking and returns the corrected Seurat
#' object from the best-performing method based on the composite score.
#'
#' @param results An `OmicsBatchR` list from `run_batch_correction()`.
#' @param batch_key Character. Batch column name.
#' @param celltype_key Character or NULL. Cell type column for cLISI.
#' @param verbose Logical. Default: `TRUE`.
#'
#' @return A list with:
#'   - `best_method`: name of the best method
#'   - `seurat_obj`: corrected Seurat object
#'   - `benchmark_scores`: full benchmark table
#'
#' @examples
#' \dontrun{
#' results <- run_batch_correction(obj, batch_key = "batch")
#' best    <- select_best_method(results, batch_key = "batch",
#'                               celltype_key = "cell_type")
#' best$best_method
#' # "harmony"
#' }
#'
#' @export
select_best_method <- function(results,
                                batch_key    = "batch",
                                celltype_key = NULL,
                                verbose      = TRUE) {

  scores <- benchmark_methods(results,
                               batch_key    = batch_key,
                               celltype_key = celltype_key,
                               verbose      = verbose)

  best_method <- scores$method[1]

  if (verbose) {
    cat("\n")
    .success(sprintf("Selected best method: %s", toupper(best_method)))
  }

  return(list(
    best_method      = best_method,
    seurat_obj       = results[[best_method]],
    benchmark_scores = scores
  ))
}


# ── Internal: LISI computation ─────────────────────────────────────────────────
# Simplified LISI (Local Inverse Simpson's Index)
# For full implementation see: Korsunsky et al. 2019 (Nature Methods)
.compute_lisi <- function(embedding, labels, k = 30) {
  n_cells   <- nrow(embedding)
  labels    <- as.factor(labels)
  n_labels  <- nlevels(labels)

  # Find k nearest neighbours
  knn_result <- BiocNeighbors::findKNN(embedding, k = k, BNPARAM = BiocNeighbors::AnnoyParam())
  knn_idx    <- knn_result$index

  lisi_scores <- numeric(n_cells)
  for (i in seq_len(n_cells)) {
    neighbour_labels <- labels[knn_idx[i, ]]
    label_counts     <- table(neighbour_labels)
    label_props      <- label_counts / sum(label_counts)
    # Inverse Simpson's index
    lisi_scores[i]   <- 1 / sum(label_props^2)
  }

  return(lisi_scores)
}


# ── Internal: kBET approximation ───────────────────────────────────────────────
# Simplified kBET: fraction of neighbourhoods where batch composition
# matches global batch composition (chi-squared test acceptance rate)
.compute_kbet <- function(embedding, batch_labels, k = 30) {
  n_cells      <- nrow(embedding)
  batch_labels <- as.factor(batch_labels)
  global_props <- prop.table(table(batch_labels))

  knn_result <- BiocNeighbors::findKNN(embedding, k = k, BNPARAM = BiocNeighbors::AnnoyParam())
  knn_idx    <- knn_result$index

  # Subsample 1000 cells for speed
  test_cells  <- sample(seq_len(n_cells), min(1000, n_cells))
  accept_vec  <- logical(length(test_cells))

  for (j in seq_along(test_cells)) {
    i                  <- test_cells[j]
    neighbour_labels   <- batch_labels[knn_idx[i, ]]
    observed_counts    <- table(factor(neighbour_labels, levels = levels(batch_labels)))
    expected_counts    <- global_props * k

    # Chi-squared test
    chi_result         <- suppressWarnings(
      chisq.test(observed_counts, p = global_props)
    )
    accept_vec[j]      <- chi_result$p.value > 0.05
  }

  return(mean(accept_vec))
}


# ── Internal: composite score ──────────────────────────────────────────────────
.compute_overall_score <- function(scores_df) {
  # Normalise each metric to 0-1, then combine
  norm01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
  }

  # iLISI: higher = better mixing (good)
  ilisi_norm <- norm01(scores_df$iLISI)

  # cLISI: lower = better cell-type preservation (invert)
  if (all(is.na(scores_df$cLISI))) {
    clisi_norm <- rep(0.5, nrow(scores_df))
  } else {
    clisi_norm <- 1 - norm01(scores_df$cLISI)
  }

  # kBET: higher = better (good)
  kbet_norm <- norm01(scores_df$kBET_acceptance)

  # Weighted composite: iLISI (40%) + cLISI (30%) + kBET (30%)
  scores_df$overall_score <- round(
    0.40 * ilisi_norm + 0.30 * clisi_norm + 0.30 * kbet_norm,
    4
  )

  return(scores_df)
}
