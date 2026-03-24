#' Run multiple batch correction methods on a Seurat object
#'
#' @description
#' `run_batch_correction()` is the main entry point of OmicsBatchR. It accepts
#' a Seurat object and runs one or more batch correction methods, returning a
#' list of corrected Seurat objects ready for downstream analysis.
#'
#' Supported methods:
#' - **harmony** — Harmony integration (fast, widely used)
#' - **seurat_cca** — Seurat CCA (canonical correlation analysis)
#' - **combat** — ComBat-seq (count-level batch correction)
#' - **fastmnn** — Fast mutual nearest neighbours (batchelor)
#' - **scanorama** — Scanorama (Python-based, via reticulate)
#' - **bbknn** — Batch balanced KNN (Python-based, via reticulate)
#' - **scvi** — scVI deep learning model (Python-based, via reticulate)
#'
#' @param seurat_obj A Seurat object with a batch variable in metadata.
#' @param batch_key Character. Column name in `seurat_obj@meta.data` that
#'   identifies the batch. Default: `"batch"`.
#' @param methods Character vector. One or more of: `"harmony"`,
#'   `"seurat_cca"`, `"combat"`, `"fastmnn"`, `"scanorama"`, `"bbknn"`,
#'   `"scvi"`, or `"all"`. Default: `"all"`.
#' @param n_dims Integer. Number of PCA dimensions to use. Default: `30`.
#' @param n_features Integer. Number of variable features to use. Default: `2000`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#' @param python_env Character. Name of the conda environment for
#'   Python-based methods (scVI, Scanorama, BBKNN). Default: `"r-reticulate"`.
#'
#' @return A named list of Seurat objects (class `OmicsBatchR`), one per
#'   method. Each object has a corrected dimensionality reduction and UMAP
#'   computed on the corrected embedding.
#'
#' @examples
#' \dontrun{
#' library(OmicsBatchR)
#' library(Seurat)
#'
#' results <- run_batch_correction(
#'   seurat_obj = your_seurat,
#'   batch_key  = "batch",
#'   methods    = "all"
#' )
#'
#' # Run specific methods only
#' results <- run_batch_correction(
#'   seurat_obj = your_seurat,
#'   batch_key  = "batch",
#'   methods    = c("harmony", "seurat_cca", "fastmnn")
#' )
#' }
#'
#' @export
run_batch_correction <- function(seurat_obj,
                                 batch_key   = "batch",
                                 methods     = "all",
                                 n_dims      = 30,
                                 n_features  = 2000,
                                 verbose     = TRUE,
                                 python_env  = "r-reticulate") {

  # ── Validate inputs ─────────────────────────────────────────────────────────
  .validate_inputs(seurat_obj, batch_key, methods)

  all_methods <- c("harmony", "seurat_cca", "combat", "fastmnn",
                   "scanorama", "bbknn", "scvi")
  if ("all" %in% methods) methods <- all_methods

  python_methods <- c("scanorama", "bbknn", "scvi")
  needs_python   <- any(methods %in% python_methods)

  if (needs_python) .setup_python(python_env, verbose)

  # ── Preprocessing ────────────────────────────────────────────────────────────
  if (verbose) .msg("Preprocessing Seurat object...")
  seurat_obj <- .preprocess_seurat(seurat_obj, n_features, n_dims)

  # ── Run each method ──────────────────────────────────────────────────────────
  results <- list()

  for (method in methods) {
    if (verbose) .msg(sprintf("Running: %s", toupper(method)))

    tryCatch({
      result <- switch(method,
        harmony    = .run_harmony(seurat_obj, batch_key, n_dims),
        seurat_cca = .run_seurat_cca(seurat_obj, batch_key, n_dims, n_features),
        combat     = .run_combat(seurat_obj, batch_key, n_dims),
        fastmnn    = .run_fastmnn(seurat_obj, batch_key, n_dims),
        scanorama  = .run_scanorama(seurat_obj, batch_key, n_dims),
        bbknn      = .run_bbknn(seurat_obj, batch_key, n_dims),
        scvi       = .run_scvi(seurat_obj, batch_key, n_dims)
      )
      results[[method]] <- result
      if (verbose) .success(sprintf("%s complete", toupper(method)))
    }, error = function(e) {
      warning(sprintf("Method '%s' failed: %s", method, e$message))
      if (verbose) .warn(sprintf("Skipping %s due to error: %s", method, e$message))
    })
  }

  if (length(results) == 0) stop("All batch correction methods failed.")

  class(results)              <- c("OmicsBatchR", "list")
  attr(results, "batch_key") <- batch_key
  attr(results, "n_dims")    <- n_dims
  attr(results, "methods_run") <- names(results)

  if (verbose) .success(sprintf("Done! %d method(s) completed successfully.", length(results)))
  return(results)
}


# ── Internal: preprocessing ───────────────────────────────────────────────────

.preprocess_seurat <- function(obj, n_features, n_dims) {
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_features, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = n_dims, verbose = FALSE)
  return(obj)
}


# ── Internal: Harmony ─────────────────────────────────────────────────────────

.run_harmony <- function(obj, batch_key, n_dims) {
  obj <- harmony::RunHarmony(
    obj,
    group.by.vars = batch_key,
    reduction     = "pca",
    dims.use      = seq_len(n_dims),
    verbose       = FALSE
  )
  obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = seq_len(n_dims),
                         reduction.name = "umap_harmony", verbose = FALSE)
  obj@misc$correction_method <- "harmony"
  return(obj)
}


# ── Internal: Seurat CCA ──────────────────────────────────────────────────────

.run_seurat_cca <- function(obj, batch_key, n_dims, n_features) {
  obj_list <- Seurat::SplitObject(obj, split.by = batch_key)

  obj_list <- lapply(obj_list, function(x) {
    x <- Seurat::NormalizeData(x, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(x, nfeatures = n_features, verbose = FALSE)
    return(x)
  })

  anchors <- Seurat::FindIntegrationAnchors(
    object.list     = obj_list,
    dims            = seq_len(n_dims),
    anchor.features = n_features,
    verbose         = FALSE
  )

  integrated <- Seurat::IntegrateData(anchorset = anchors, dims = seq_len(n_dims),
                                      verbose = FALSE)
  Seurat::DefaultAssay(integrated) <- "integrated"
  integrated <- Seurat::ScaleData(integrated, verbose = FALSE)
  integrated <- Seurat::RunPCA(integrated, npcs = n_dims,
                               reduction.name = "pca_cca", verbose = FALSE)
  integrated <- Seurat::RunUMAP(integrated, reduction = "pca_cca",
                                dims = seq_len(n_dims),
                                reduction.name = "umap_seurat_cca",
                                verbose = FALSE)
  integrated@misc$correction_method <- "seurat_cca"
  return(integrated)
}


# ── Internal: ComBat-seq ──────────────────────────────────────────────────────

.run_combat <- function(obj, batch_key, n_dims) {
  count_mat <- as.matrix(Seurat::GetAssayData(obj, slot = "counts"))
  batch_vec <- obj@meta.data[[batch_key]]

  corrected_counts <- sva::ComBat_seq(counts = count_mat, batch = batch_vec)

  obj_combat <- Seurat::SetAssayData(obj, slot = "counts", new.data = corrected_counts)
  obj_combat <- Seurat::NormalizeData(obj_combat, verbose = FALSE)
  obj_combat <- Seurat::FindVariableFeatures(obj_combat, verbose = FALSE)
  obj_combat <- Seurat::ScaleData(obj_combat, verbose = FALSE)
  obj_combat <- Seurat::RunPCA(obj_combat, npcs = n_dims,
                               reduction.name = "pca_combat", verbose = FALSE)
  obj_combat <- Seurat::RunUMAP(obj_combat, reduction = "pca_combat",
                                dims = seq_len(n_dims),
                                reduction.name = "umap_combat",
                                verbose = FALSE)
  obj_combat@misc$correction_method <- "combat"
  return(obj_combat)
}


# ── Internal: FastMNN ─────────────────────────────────────────────────────────

.run_fastmnn <- function(obj, batch_key, n_dims) {
  sce     <- Seurat::as.SingleCellExperiment(obj)
  mnn_out <- batchelor::fastMNN(
    sce,
    batch   = sce[[batch_key]],
    d       = n_dims,
    BSPARAM = BiocSingular::RandomParam()
  )

  mnn_embedding <- SingleCellExperiment::reducedDim(mnn_out, "corrected")
  rownames(mnn_embedding) <- colnames(obj)

  obj[["mnn"]] <- Seurat::CreateDimReducObject(
    embeddings = mnn_embedding,
    key        = "MNN_",
    assay      = Seurat::DefaultAssay(obj)
  )
  obj <- Seurat::RunUMAP(obj, reduction = "mnn", dims = seq_len(n_dims),
                         reduction.name = "umap_fastmnn", verbose = FALSE)
  obj@misc$correction_method <- "fastmnn"
  return(obj)
}


# ── Internal: Scanorama (Python) ──────────────────────────────────────────────

.run_scanorama <- function(obj, batch_key, n_dims) {
  scanorama    <- reticulate::import("scanorama")
  batch_levels <- unique(obj@meta.data[[batch_key]])

  datasets   <- list()
  gene_lists <- list()

  for (b in batch_levels) {
    cells         <- rownames(obj@meta.data[obj@meta.data[[batch_key]] == b, , drop = FALSE])
    mat           <- t(as.matrix(Seurat::GetAssayData(obj, slot = "data")[, cells]))
    datasets[[b]]   <- mat
    gene_lists[[b]] <- colnames(mat)
  }

  integrated_data <- scanorama$integrate(datasets, gene_lists)
  corrected       <- do.call(rbind, integrated_data[[1]])

  cell_order <- unlist(lapply(batch_levels, function(b) {
    rownames(obj@meta.data[obj@meta.data[[batch_key]] == b, , drop = FALSE])
  }))
  rownames(corrected) <- cell_order
  corrected           <- corrected[colnames(obj), ]

  obj[["scanorama"]] <- Seurat::CreateDimReducObject(
    embeddings = corrected[, seq_len(min(n_dims, ncol(corrected))), drop = FALSE],
    key        = "SCANORAMA_",
    assay      = Seurat::DefaultAssay(obj)
  )
  obj <- Seurat::RunUMAP(obj, reduction = "scanorama",
                         dims = seq_len(min(n_dims, ncol(corrected))),
                         reduction.name = "umap_scanorama", verbose = FALSE)
  obj@misc$correction_method <- "scanorama"
  return(obj)
}


# ── Internal: BBKNN (Python) ──────────────────────────────────────────────────

.run_bbknn <- function(obj, batch_key, n_dims) {
  anndata <- reticulate::import("anndata")
  bbknn   <- reticulate::import("bbknn")
  sc      <- reticulate::import("scanpy")

  mat    <- t(as.matrix(Seurat::GetAssayData(obj, slot = "data")))
  obs_df <- obj@meta.data
  adata  <- anndata$AnnData(X = mat, obs = obs_df)

  sc$tl$pca(adata, n_comps = as.integer(n_dims))
  bbknn$bbknn(adata, batch_key = batch_key)
  sc$tl$umap(adata)

  umap_coords <- adata$obsm[["X_umap"]]
  pca_coords  <- adata$obsm[["X_pca"]]
  rownames(umap_coords) <- colnames(obj)
  rownames(pca_coords)  <- colnames(obj)

  obj[["bbknn_pca"]] <- Seurat::CreateDimReducObject(
    embeddings = pca_coords, key = "BBKNN_",
    assay = Seurat::DefaultAssay(obj)
  )
  obj[["umap_bbknn"]] <- Seurat::CreateDimReducObject(
    embeddings = umap_coords, key = "UMAP_BBKNN_",
    assay = Seurat::DefaultAssay(obj)
  )
  obj@misc$correction_method <- "bbknn"
  return(obj)
}


# ── Internal: scVI (Python) ───────────────────────────────────────────────────

.run_scvi <- function(obj, batch_key, n_dims) {
  scvi    <- reticulate::import("scvi")
  anndata <- reticulate::import("anndata")

  counts_mat <- t(as.matrix(Seurat::GetAssayData(obj, slot = "counts")))
  obs_df     <- obj@meta.data
  adata      <- anndata$AnnData(X = counts_mat, obs = obs_df)

  scvi$model$SCVI$setup_anndata(adata, batch_key = batch_key)
  model <- scvi$model$SCVI(adata, n_latent = as.integer(n_dims))
  model$train(max_epochs = as.integer(400), early_stopping = TRUE)

  latent              <- model$get_latent_representation()
  rownames(latent)    <- colnames(obj)

  obj[["scvi"]] <- Seurat::CreateDimReducObject(
    embeddings = latent, key = "SCVI_",
    assay = Seurat::DefaultAssay(obj)
  )
  obj <- Seurat::RunUMAP(obj, reduction = "scvi", dims = seq_len(n_dims),
                         reduction.name = "umap_scvi", verbose = FALSE)
  obj@misc$correction_method <- "scvi"
  return(obj)
}


# ── Internal: input validation ────────────────────────────────────────────────

.validate_inputs <- function(obj, batch_key, methods) {
  if (!inherits(obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.", call. = FALSE)
  }
  if (!batch_key %in% colnames(obj@meta.data)) {
    stop(sprintf(
      "batch_key '%s' not found in seurat_obj@meta.data.\nAvailable columns: %s",
      batch_key, paste(colnames(obj@meta.data), collapse = ", ")
    ), call. = FALSE)
  }
  valid_methods <- c("all", "harmony", "seurat_cca", "combat",
                     "fastmnn", "scanorama", "bbknn", "scvi")
  bad <- setdiff(methods, valid_methods)
  if (length(bad) > 0) {
    stop(sprintf("Unknown method(s): %s\nValid: %s",
                 paste(bad, collapse = ", "),
                 paste(valid_methods, collapse = ", ")), call. = FALSE)
  }
  n_batches <- length(unique(obj@meta.data[[batch_key]]))
  if (n_batches < 2) {
    stop("At least 2 batches are required for batch correction.", call. = FALSE)
  }
}


# ── Internal: Python setup ────────────────────────────────────────────────────

.setup_python <- function(python_env, verbose) {
  if (verbose) .msg("Setting up Python environment...")
  tryCatch(
    reticulate::use_condaenv(python_env, required = FALSE),
    error = function(e) {
      warning("Could not activate conda env '", python_env,
              "'. Trying default Python.")
    }
  )
}


# ── Internal: message helpers (cli only — no crayon) ─────────────────────────

.msg     <- function(x) cli::cli_alert_info(x)
.success <- function(x) cli::cli_alert_success(x)
.warn    <- function(x) cli::cli_alert_warning(x)
