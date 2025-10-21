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

getoneperformance <- function(expr_mat, thisname, target_ks = 15, dims_to_use = 1:50, seedi = 9){
  set.seed(seedi)
  datain <- CreateSeuratObject(expr_mat)
  datain <- SCTransform(datain, verbose = FALSE)
  datain <- RunPCA(datain, assay = "SCT", verbose = FALSE)
  datain <- FindNeighbors(datain, reduction = "pca", dims = dims_to_use, verbose = FALSE)
  
  res_obj <- findResolutionForK(datain, target_ks)  
  if (is.null(res_obj)) {
    return(c(name = thisname, k = target_ks, resolution = NA, silhouette = NA, CH = NA, DB = NA))
  }else{                
    
    res_val   <- res_obj$resolution
    clusters  <- res_obj$clusters                    
    clusters  <- as.integer(factor(clusters))        
    
    pca_embed <- datain@reductions$pca@cell.embeddings[, dims_to_use]
    dmat      <- dist(pca_embed)
    sil_avg   <- mean(cluster::silhouette(clusters, dmat)[, 3])
    
    int_res   <- clusterCrit::intCriteria( as.matrix(pca_embed), clusters, c("Calinski_Harabasz", "Davies_Bouldin"))
    ch_idx    <- int_res$calinski_harabasz
    db_idx    <- int_res$davies_bouldin
    
    return(c(name = thisname, k = target_ks, resolution = res_val, silhouette = sil_avg, CH = ch_idx, DB = db_idx))
  }
}

summarylist = list()
#  summarydf <- do.call(rbind, summarylist);summarydf
#  write.table(summarydf, outputfile, sep = '\t', row.names = FALSE)

#{outputfile=  'benchbypixelin100.txt'
#pixelin <- read.table('pixelin.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)}

{outputfile=  'benchbypixelin10.txt'
pixelin <- read.table('pixelin10.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)}

{outputfile=  'benchbypixelin200.txt'
pixelin <- read.table('pixelin200.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)}

{outputfile=  'benchbypixelin10_1_50.txt'
  pixelin <- read.table('pixelin10.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)}

expr_mat = Read10X('data/STARREO/'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,1]]
summarylist[['starreo_highexp_points']] <- getoneperformance(expr_mat, 'starreo_highexp_points')
expr_mat = Read10X('data/STAR/'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,1]]
summarylist[['star_highexp_points']] <- getoneperformance(expr_mat, 'star_highexp_points')
##
expr_mat = quickload('output/filteredout/addlowq_expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
summarylist[['addlowq_expmat']] <- getoneperformance(expr_mat, 'addlowq_expmat')
expr_mat = quickload('output/expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
summarylist[['expmat_tsv_new']] <- getoneperformance(expr_mat, 'expmat_tsv_new')
expr_mat = quickload('/gpfs/gibbs/pi/gerstein/cz357/25fall/ASTRO/STRS_control_output/finalexpmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
summarylist[['finalexpmat_old']] <- getoneperformance(expr_mat, 'finalexpmat_old')
expr_mat = quickload('/gpfs/gibbs/pi/gerstein/cz357/25fall/ASTRO/STRS_control_output/expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
summarylist[['expmat_tsv_old']] <- getoneperformance(expr_mat, 'expmat_tsv_old')
##


expr_mat = quickload('output/filteredout/addlowq_expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
names_collapsed = gsub('(__exon|__intron)', '', rownames(expr_mat));expr_mat <- rowsum(expr_mat, group = names_collapsed, reorder = FALSE)
summarylist[['addlowq_expmat_collapsed']] <- getoneperformance(expr_mat, 'addlowq_expmat_collapsed')
expr_mat = quickload('output/expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
names_collapsed = gsub('(__exon|__intron)', '', rownames(expr_mat));expr_mat <- rowsum(expr_mat, group = names_collapsed, reorder = FALSE)
summarylist[['expmat_tsv_new_collapsed']] <- getoneperformance(expr_mat, 'expmat_tsv_new_collapsed')
expr_mat = quickload('/gpfs/gibbs/pi/gerstein/cz357/25fall/ASTRO/STRS_control_output/finalexpmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
names_collapsed = gsub('(__exon|__intron)', '', rownames(expr_mat));expr_mat <- rowsum(expr_mat, group = names_collapsed, reorder = FALSE)
summarylist[['finalexpmat_old_collapsed']] <- getoneperformance(expr_mat, 'finalexpmat_old_collapsed')
expr_mat = quickload('/gpfs/gibbs/pi/gerstein/cz357/25fall/ASTRO/STRS_control_output/expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
names_collapsed = gsub('(__exon|__intron)', '', rownames(expr_mat));expr_mat <- rowsum(expr_mat, group = names_collapsed, reorder = FALSE)
summarylist[['expmat_tsv_old_collapsed']] <- getoneperformance(expr_mat, 'expmat_tsv_old_collapsed')

##
expr_mat = quickload('/gpfs/gibbs/pi/gerstein/cz357/25fall/ASTRO/STRS_control_output/finalexpmat.tsv'); expr_mat = expr_mat[!grepl('\\-\\-', rownames(expr_mat)),];expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
summarylist[['finalexpmat_old_nodoubledash']] <- getoneperformance(expr_mat, 'finalexpmat_old_nodoubledash')
expr_mat = quickload('/gpfs/gibbs/pi/gerstein/cz357/25fall/ASTRO/STRS_control_output/expmat.tsv'); expr_mat = expr_mat[!grepl('\\-\\-', rownames(expr_mat)),];expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
summarylist[['expmat_old_nodoubledash']] <- getoneperformance(expr_mat, 'expmat_old_nodoubledash')
expr_mat = quickload('output/filteredout/addlowq_expmat.tsv'); expr_mat = expr_mat[!grepl('\\+', rownames(expr_mat)),];expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
summarylist[['addlowq_expmat_noplus']] <- getoneperformance(expr_mat, 'addlowq_expmat_noplus')
expr_mat = quickload('output/expmat.tsv'); expr_mat = expr_mat[!grepl('\\+', rownames(expr_mat)),];expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
summarylist[['expmat_tsv_new_noplus']] <- getoneperformance(expr_mat, 'expmat_tsv_new_noplus')

## 
expr_mat = quickload('output/filteredout/addlowq_expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
filter1 = !grepl('^(mt|ENSMUS|Rps|Rpl|Gm1|Gm2|Gm3|Gm4|Gm5|Gm6|Gm7|Gm8|Gm9|0|1|2|3|4|5|6|7|8|9|A03|A13|A23|A33|A43|A53|A53|A63|A73|A83|A93|AC0|AC1|AL|AU|AW|B02|B13|B23|B43|B93|BC0|C03|C13|C23|C33|C43|C53|C63|C92|D03|D13|D23|D33|D43|D53|D63|D73|D83|D93|E03|E13|E23|E33|E43|E53|F63|F73|F83|F93|G43|G53|G63|G73|I83|n-R5|RP23|RP24)',rownames(expr_mat));expr_mat = expr_mat[filter1,]
summarylist[['addlowq_expmat_geneselected']] <- getoneperformance(expr_mat, 'addlowq_expmat_geneselected')
expr_mat = quickload('output/expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
filter1 = !grepl('^(mt|ENSMUS|Rps|Rpl|Gm1|Gm2|Gm3|Gm4|Gm5|Gm6|Gm7|Gm8|Gm9|0|1|2|3|4|5|6|7|8|9|A03|A13|A23|A33|A43|A53|A53|A63|A73|A83|A93|AC0|AC1|AL|AU|AW|B02|B13|B23|B43|B93|BC0|C03|C13|C23|C33|C43|C53|C63|C92|D03|D13|D23|D33|D43|D53|D63|D73|D83|D93|E03|E13|E23|E33|E43|E53|F63|F73|F83|F93|G43|G53|G63|G73|I83|n-R5|RP23|RP24)',rownames(expr_mat));expr_mat = expr_mat[filter1,]
summarylist[['expmat_tsv_new_geneselected']] <- getoneperformance(expr_mat, 'expmat_tsv_new_geneselected')
expr_mat = quickload('/gpfs/gibbs/pi/gerstein/cz357/25fall/ASTRO/STRS_control_output/finalexpmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
filter1 = !grepl('^(mt|ENSMUS|Rps|Rpl|Gm1|Gm2|Gm3|Gm4|Gm5|Gm6|Gm7|Gm8|Gm9|0|1|2|3|4|5|6|7|8|9|A03|A13|A23|A33|A43|A53|A53|A63|A73|A83|A93|AC0|AC1|AL|AU|AW|B02|B13|B23|B43|B93|BC0|C03|C13|C23|C33|C43|C53|C63|C92|D03|D13|D23|D33|D43|D53|D63|D73|D83|D93|E03|E13|E23|E33|E43|E53|F63|F73|F83|F93|G43|G53|G63|G73|I83|n-R5|RP23|RP24)',rownames(expr_mat));expr_mat = expr_mat[filter1,]
summarylist[['finalexpmat_old_geneselected']] <- getoneperformance(expr_mat, 'finalexpmat_old_geneselected')
expr_mat = quickload('/gpfs/gibbs/pi/gerstein/cz357/25fall/ASTRO/STRS_control_output/expmat.tsv'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,2]]
filter1 = !grepl('^(mt|ENSMUS|Rps|Rpl|Gm1|Gm2|Gm3|Gm4|Gm5|Gm6|Gm7|Gm8|Gm9|0|1|2|3|4|5|6|7|8|9|A03|A13|A23|A33|A43|A53|A53|A63|A73|A83|A93|AC0|AC1|AL|AU|AW|B02|B13|B23|B43|B93|BC0|C03|C13|C23|C33|C43|C53|C63|C92|D03|D13|D23|D33|D43|D53|D63|D73|D83|D93|E03|E13|E23|E33|E43|E53|F63|F73|F83|F93|G43|G53|G63|G73|I83|n-R5|RP23|RP24)',rownames(expr_mat));expr_mat = expr_mat[filter1,]
summarylist[['expmat_tsv_old_geneselected']] <- getoneperformance(expr_mat, 'expmat_tsv_old_geneselected')


#expr_mat = Read10X('data/kallisto/'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,1]]
#summarylist[['kallisto_highexp_points']] <- getoneperformance(expr_mat, 'kallisto_highexp_points')
#expr_mat = Read10X('data/kallistoREO/'); expr_mat <- expr_mat[, colnames(expr_mat) %in% pixelin[,1]]
#summarylist[['kallistoreo_highexp_points']] <- getoneperformance(expr_mat, 'kallistoreo_highexp_points')

summarydf <- do.call(rbind, summarylist)
write.table(summarydf, outputfile, sep = '\t', row.names = FALSE)

