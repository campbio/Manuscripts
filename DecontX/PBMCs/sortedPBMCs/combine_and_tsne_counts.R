library(Matrix)
library(Rtsne)

setwd("/restricted/projectnb/camplab/projects/celda/Datasets/10X")

makeMatrix = function(dataset.name, 
                      path.suffix="/filtered_matrices_mex/hg19/") {
    path = paste0(dataset.name, path.suffix)
    mm <- as.matrix(readMM(paste0(path, "matrix.mtx")))
    barcodes <- read.table(paste0(path, "barcodes.tsv"), sep = "\t")
    genes <- read.table(paste0(path, "genes.tsv"), sep = "\t")
    rownames(mm) <- genes$V2
    colnames(mm) <- barcodes$V1
    return(mm)
}

dataset.names = c("CD4_Helper_Tcells", 
                  "CD34_cells", 
                  "CD56_Natural_Killer_Cells", 
                  "CD4_CD25_Regulatory_Tcells", 
                  "CD8_CD45RA_Naive_Cytotoxic_Tcells",
                  "CD4_CD45RA_CD25neg_Naive_Tcells",
                  "CD8_Cytotoxic_Tcells", 
                  #"CD4_CD45RO_Memory_ Tcells", 
                  "CD14_Monocytes", 
                  "CD19_Bcells")

all.counts = lapply(dataset.names, makeMatrix)
dataset.labels = lapply(1:length(all.counts),
                        function(idx) {
                            num.reps = ncol(all.counts[[idx]])
                            return(rep(dataset.names[idx], times=num.reps))
                        })
dataset.labels = unlist(dataset.labels)

all.immune.counts = do.call(cbind, all.counts)
saveRDS(all.immune.counts, "Combined_Immune_Datasets/all_immune_counts.Rds")

filtered.immune.counts =  all.immune.counts[rowSums(all.immune.counts > 2) > 2, 
                                            colSums(all.immune.counts > 2) > 2] 
filtered.immune.tsne = Rtsne(t(filtered.immune.counts))
filtered.immune.tsne$tsne.df = as.data.frame(filtered.immune.tsne$Y)
colnames(filtered.immune.tsne$tsne.df) = c("tsne1", "tsne2")
filtered.immune.tsne$filtered.immune.tsne$dataset = dataset.labels
saveRDS(filtered.immune.tsne, file="Combined_Immune_Datasets/all_immune_counts_tsne.Rds")
