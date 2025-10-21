
library(Seurat)
library(Matrix)
library(cluster)       
library(clusterCrit)   

#tsv_path   <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/smart/astro_928/expmat.tsv"
#tsv_path   <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/smart/GSE151334_counts.human.tsv"
#outdir     <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/smart/benchmarking"
#tsv_outdir <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/smart/benchmarking"

tsv_path   <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/outputvisium/gene_by_xy_matrix.tsv"  # 行=gene, 列=spot
outdir     <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/outputvisium/"
tsv_outdir <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/outputvisium/"
metrics_csv <- file.path(tsv_outdir, "spaceranger_metrics_summary1012.csv")
#tsv_path   <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/finalexpmat_filtered.tsv" 
#outdir     <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/"
#tsv_outdir <- "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/"
#metrics_csv <- file.path(tsv_outdir, "astro_metrics_summary1012.csv")
n_subsamples     <- 20
pcs              <- 50
seed_base        <- 1
target_clusters  <- 15                  
downsample_probs <- c(0.5)
max_cells_for_sil <- 10000              

findResolutionForK <- function(seurat_obj, target_clusters) {
  low <- 0; high <- 4
  for (i in seq_len(50)) {
    mid <- (low + high) / 2
    message(sprintf("Testing resolution: %.6f", mid))
    seurat_obj <- suppressMessages(
      Seurat::FindClusters(seurat_obj, resolution = mid, verbose = FALSE)
    )
    n_clusters <- length(unique(seurat_obj@meta.data$seurat_clusters))
    if (n_clusters == target_clusters) {
      message(sprintf("Found resolution: %.6f", mid))
      return(list(resolution = mid))
    }
    if (n_clusters > target_clusters) high <- mid else low <- mid
  }
  message("Could not reach target number of clusters")
  return(NULL)
}


subsample_counts_dgC <- function(dgc, fraction) {
  if (!is(dgc, "dgCMatrix")) dgc <- as(dgc, "dgCMatrix")
  out <- dgc
  out@x <- as.numeric(stats::rbinom(length(dgc@x), size = dgc@x, prob = fraction))
  nz <- out@x != 0
  out@x <- out@x[nz]
  new_i <- vector("list", ncol(out))
  p <- out@p
  k <- 1L
  for (j in seq_len(ncol(out))) {
    idx <- seq.int(p[j] + 1L, p[j+1L])
    keep <- nz[idx]
    new_i[[j]] <- out@i[idx][keep]
  }
  out@i <- as.integer(unlist(new_i, use.names = FALSE))
  lens <- vapply(new_i, length, integer(1))
  out@p <- c(0L, cumsum(as.integer(lens)))
  out
}


#read_tsv_to_seurat <- function(path) {
#  message("Reading TSV: ", path)
#  mat <- read.table(path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
#  mat <- as.matrix(mat)
#  storage.mode(mat) <- "numeric"
#  dgc <- as(mat, "dgCMatrix")
#  CreateSeuratObject(counts = dgc)
#}
read_tsv_to_seurat <- function(path, dedup = c("sum", "unique", "first")) {
  dedup <- match.arg(dedup)
  message("Reading TSV: ", path)
  
  tab <- read.table(path, header = TRUE, sep = "\t",
                    check.names = FALSE, quote = "", comment.char = "")
  genes <- as.character(tab[[1]])
  mat   <- as.matrix(tab[ , -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  
  if (any(duplicated(colnames(mat)))) {
    warning("Duplicated column (spot) names detected; making them unique with make.unique().")
    colnames(mat) <- make.unique(colnames(mat))
  }
  
  if (any(duplicated(genes))) {
    ndup <- sum(duplicated(genes))
    message("Found duplicated gene names: ", ndup, " duplicates.")
    if (dedup == "sum") {
      mat <- rowsum(mat, group = genes, reorder = FALSE)
      genes <- rownames(mat)  
    } else if (dedup == "first") {
      keep <- !duplicated(genes)
      mat  <- mat[keep, , drop = FALSE]
      genes <- genes[keep]
    } else if (dedup == "unique") {
      genes <- make.unique(genes)
    }
  }
  
  rownames(mat) <- genes
  dgc <- as(mat, "dgCMatrix")
  CreateSeuratObject(counts = dgc)
}


compute_metrics <- function(obj, pcs = 50, max_cells_for_sil = 10000) {
  emb <- as.matrix(obj@reductions$pca@cell.embeddings[, seq_len(min(pcs, ncol(obj@reductions$pca@cell.embeddings))), drop=FALSE])
  cl  <- obj@meta.data$seurat_clusters
  lbl <- as.integer(as.factor(cl))  # 转为 1..K
  
  set.seed(1)
  idx <- seq_len(nrow(emb))
  if (length(idx) > max_cells_for_sil) idx <- sample(idx, max_cells_for_sil)
  sil_mean <- NA_real_
  if (length(unique(lbl[idx])) > 1) {
    d <- dist(emb[idx, , drop=FALSE])
    sil <- silhouette(lbl[idx], d)
    sil_mean <- mean(sil[, "sil_width"])
  }
  

  ch <- NA_real_; db <- NA_real_
  crit <- tryCatch({
    intCriteria(traj = emb, part = lbl, c("Calinski_Harabasz", "Davies_Bouldin"))
  }, error = function(e) NULL)
  if (!is.null(crit)) {
    if (!is.null(crit$calinski_harabasz)) ch <- as.numeric(crit$calinski_harabasz)
    if (!is.null(crit$davies_bouldin))    db <- as.numeric(crit$davies_bouldin)
  }
  
  data.frame(
    n_cells      = nrow(emb),
    n_clusters   = length(unique(lbl)),
    silhouette   = sil_mean,
    CH           = ch,
    DB           = db,
    stringsAsFactors = FALSE
  )
}

process_tsv <- function() {
  # 读 TSV → Seurat
  base_name <- tools::file_path_sans_ext(basename(tsv_path))
  obj0 <- read_tsv_to_seurat(tsv_path)
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tsv_outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (!file.exists(metrics_csv)) {
    header <- data.frame(file=character(), down_prob=numeric(), replicate=integer(),
                         resolution=numeric(), n_cells=integer(), n_clusters=integer(),
                         silhouette=numeric(), CH=numeric(), DB=numeric())
    write.table(header, file = metrics_csv, sep = ",", row.names = FALSE, quote = TRUE, col.names = TRUE)
  }
  
  for (down_prob in downsample_probs) {
    for (i in seq_len(n_subsamples)) {
      set.seed(seed_base + i)
      
      counts0 <- obj0@assays$RNA@counts
      counts_sub <- subsample_counts_dgC(counts0, fraction = down_prob)
      obj <- CreateSeuratObject(counts = counts_sub)
      
      obj <- subset(obj, subset = nCount_RNA > 0)
      if (ncol(obj) < target_clusters) {
        message(sprintf("[skip] prob=%.2f rep=%d: cells=%d < target_clusters=%d",
                        down_prob, i, ncol(obj), target_clusters))
        next
      }
      
      obj <- SCTransform(obj, verbose = FALSE)
      obj <- RunPCA(obj, npcs = pcs, verbose = FALSE)
      obj <- FindNeighbors(obj, dims = 1:pcs, verbose = FALSE)
      
      res_res <- findResolutionForK(obj, target_clusters = target_clusters)
      if (is.null(res_res)) {
        message(sprintf("[skip] prob=%.2f rep=%d: resolution not found", down_prob, i))
        next
      }
      best_res <- res_res$resolution
      obj <- FindClusters(obj, resolution = best_res, verbose = FALSE)
      
      m <- compute_metrics(obj, pcs = pcs, max_cells_for_sil = max_cells_for_sil)
      
      line <- data.frame(
        file       = base_name,
        down_prob  = down_prob,
        replicate  = i,
        resolution = best_res,
        m,
        check.names = FALSE
      )
      write.table(line, file = metrics_csv, sep = ",", row.names = FALSE, quote = TRUE, col.names = FALSE, append = TRUE)
      
   
      pca_df <- as.data.frame(obj@reductions$pca@cell.embeddings)
      comb   <- cbind(Cell = rownames(pca_df), pca_df, Clusters = obj@meta.data$seurat_clusters)
      out_tsv <- file.path(tsv_outdir, sprintf("%s_down%.2f_%d_pca_clustering.tsv", base_name, down_prob, i))
      #write.table(comb, file = out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
      
      
      out_rds <- file.path(outdir, sprintf("%s_down%.2f_%d.rds", base_name, down_prob, i))
      #saveRDS(obj, out_rds)
      
      message(sprintf("[ok] %s prob=%.2f rep=%d | res=%.4f | sil=%.3f | CH=%.1f | DB=%.3f",
                      base_name, down_prob, i, best_res, m$silhouette, m$CH, m$DB))
    }
  }
  
  message("All done. Metrics at: ", metrics_csv)
}

process_tsv()

