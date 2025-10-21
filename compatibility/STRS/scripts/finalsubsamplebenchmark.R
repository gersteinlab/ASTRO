library(data.table)
library(ggplot2)
library(Seurat)
library(stringr)
library(cluster)        
library(clusterCrit)
library(future)
plan(multicore, workers = 4)
options(future.globals.maxSize = 5 * 1024^3)

quickload <- function(filein) {
  matin <- fread(filein, header = TRUE);rown <- unlist(matin[,1]);matin <- as.sparse(matin[, -1, with=FALSE]); rownames(matin) <- rown
  return(matin)
}

findResolutionForK <- function(seurat_obj, target_clusters) {
  low <- 0
  high <- 4
  best_res <- NULL
  best_diff <- Inf
  
  for (i in seq_len(50)) {
    mid <- (low + high) / 2
    print(paste("Testing resolution:", mid))
    seurat_obj <- suppressMessages(Seurat::FindClusters(seurat_obj, resolution = mid, verbose = FALSE))
    temp <- seurat_obj@meta.data$seurat_clusters
    n_clusters <- length(unique(temp))
    
    if (n_clusters == target_clusters) {
      best_res <- mid
      print(paste("Found resolution:", best_res))
      return(list(resolution = best_res, clusters = temp))
    }
    if (n_clusters > target_clusters) {
      high <- mid
    } else {
      low <- mid
    }
  }
  print("Could not reach target number of clusters")
  return(NULL)
}

subsample_counts_dgC <- function(dgc, fraction) {
  dgc <- as(dgc, "dgCMatrix")
  dgc@x <- as.numeric(stats::rbinom(length(dgc@x), size = as.integer(dgc@x), prob = fraction))  
  Matrix::drop0(dgc)
}

getafile = function(expr_mat,dims_to_use,outputfile){
  metrics_list <- list()
  for (ii in 1:40) {
    set.seed(ii)
    subdatain = subsample_counts_dgC(expr_mat,0.5)
    
    keep_g <- Matrix::rowSums(subdatain) > 0
    keep_c <- Matrix::colSums(subdatain) > 0
    subdatain <- subdatain[keep_g, keep_c, drop = FALSE]
    
    subdatain <- CreateSeuratObject(subdatain)
    subdatain <- SCTransform(subdatain,  verbose = FALSE)
    subdatain <- RunPCA(subdatain, assay = "SCT", verbose = FALSE)
    subdatain <- FindNeighbors(subdatain, reduction = "pca", dims = dims_to_use, verbose = FALSE)
    #ElbowPlot(subdatain, ndims = 50)
    res_obj <- findResolutionForK(subdatain, 15)  
    if (is.null(res_obj)) {
      metrics_list[[length(metrics_list) + 1]] <- data.frame(rep = NA, resolution = NA, silhouette = NA, CH = NA, DB = NA)
    }else{
      res_val   <- res_obj$resolution
      clusters  <- res_obj$clusters                    
      clusters  <- as.integer(factor(clusters))        
      pca_embed <- subdatain@reductions$pca@cell.embeddings[, dims_to_use]
      dmat      <- dist(pca_embed)
      sil_avg   <- mean(cluster::silhouette(clusters, dmat)[, 3])
      
      int_res   <- clusterCrit::intCriteria( as.matrix(pca_embed), clusters, c("Calinski_Harabasz", "Davies_Bouldin"))
      ch_idx    <- int_res$calinski_harabasz
      db_idx    <- int_res$davies_bouldin
      
      metrics_list[[length(metrics_list) + 1]] <- data.frame(rep = ii, resolution = res_val, silhouette = sil_avg, CH = ch_idx, DB = db_idx)
    }            
  }
  metrics_df <- do.call(rbind, metrics_list)
  write.table(metrics_df, outputfile, sep = '\t', row.names = FALSE)
}

summarylist = list()


# subdatain <- CreateSeuratObject(expr_mat)
# subdatain <- SCTransform(subdatain,  verbose = FALSE)
# subdatain <- RunPCA(subdatain, assay = "SCT", verbose = FALSE)
# ElbowPlot(subdatain, ndims = 50)

pixelin <- read.table('cache/pixelin.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
dims_to_use <- 1:15
expr_mat = Read10X('data/GSM6034864/STAR/'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,1]]
subdatain <- CreateSeuratObject(expr_mat)
subdatain <- SCTransform(subdatain,  verbose = FALSE)
subdatain <- RunPCA(subdatain, assay = "SCT", verbose = FALSE)
ElbowPlot(subdatain, ndims = 50)+ggtitle('their')
outputfile = 'result/theirdownsample.txt'
getafile(expr_mat,dims_to_use,outputfile)
expr_mat = quickload('result/output/filteredout/addlowq_expmat.tsv')
expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
subdatain <- CreateSeuratObject(expr_mat)
subdatain <- SCTransform(subdatain,  verbose = FALSE)
subdatain <- RunPCA(subdatain, assay = "SCT", verbose = FALSE)
ElbowPlot(subdatain, ndims = 50)+ggtitle('ASTRO')
outputfile = 'result/astrodownsample.txt'
getafile(expr_mat,dims_to_use,outputfile)


