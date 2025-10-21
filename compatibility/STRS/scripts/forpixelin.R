library(data.table)
library(ggplot2)
library(Seurat)
library(stringr)
library(cluster)        
library(clusterCrit)
library(future)
plan(multicore, workers = 4)

matin <- fread('cache/finalexpmat.tsv ', header = TRUE)
allgenename = unlist(matin[,1])
matin = matin[unlist(matin[,1]) %in% allgenename,]
gene_ids <- matin[[1]]                      
expr_mat <- as.matrix(matin[, -1, with=FALSE])
rownames(expr_mat) <- gene_ids
cell_sum <- Matrix::colSums(expr_mat) 

expr_mat = Read10X('data/GSM6034864/STAR/')
oldcell_sum <- Matrix::colSums(expr_mat) 

pixelin <- read.table('data/visium-v1_coordinates.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
pixelin = data.frame(seq = pixelin[,1], name = paste0(pixelin[,2],'x',pixelin[,3]))
seq2name = pixelin[,2]
names(seq2name) = pixelin[,1]

plot(density(cell_sum))
plot(density(oldcell_sum))
names(oldcell_sum) = seq2name[names(oldcell_sum)]

sum(is.na(seq2name[names(oldcell_sum)]))
selectpixels = intersect(names(cell_sum[oldcell_sum > 100]), names(cell_sum[cell_sum > 100]))
outputframe = pixelin[pixelin[,2] %in% selectpixels,]

write.table(outputframe, file = "cache/pixelin100.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)

selectpixels = intersect(names(cell_sum[oldcell_sum > 10]), names(cell_sum[cell_sum > 10]))
outputframe = pixelin[pixelin[,2] %in% selectpixels,]
write.table(outputframe, file = "cache/pixelin10.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)
selectpixels = intersect(names(cell_sum[oldcell_sum > 200]), names(cell_sum[cell_sum > 200]))
outputframe = pixelin[pixelin[,2] %in% selectpixels,]
write.table(outputframe, file = "cache/pixelin.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)
