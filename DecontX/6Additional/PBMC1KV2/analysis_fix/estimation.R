library(Matrix)
library(SingleCellExperiment)
library(devtools)
load_all("~/gitProjects/celda", recompile=T)
source("/restricted/projectnb/camplab/home/syyang/contamination/automaticClustering/decontx.R")

path = "/restricted/projectnb/camplab/home/syyang/contamination/data/"
dataname = "PBMC1KV2"

data.path = file.path(path, dataname, "filtered_feature_bc_matrix")
files = list.files(data.path)

counts.path = file.path(data.path, "matrix.mtx.gz")
features.path = file.path(data.path, "features.tsv.gz") 
barcodes.path = file.path(data.path, "barcodes.tsv.gz") 

counts = readMM( counts.path) 
features = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcodes = read.delim(barcodes.path, header = FALSE, stringsAsFactors = FALSE)      


rownames(counts) = features[, 2]
colnames(counts) = barcodes[, 1]


counts2 = counts[ rowSums(counts) != 0, ]
counts2 <- counts2[rowSums(counts2>2) >2, ]


# build sce object
sce <- SingleCellExperiment(assays = list(counts = counts2))


dbscan.eps = 1
L = 100
sce_dcon = autoDecontX(sce, varGenes=nrow(sce), seed=12345, L=L, dbscan.eps=dbscan.eps, decontxIter=200)

saveRDS(sce_dcon, paste0("filterG_sce_dcon_eps", dbscan.eps, "L", L, ".rds")) 
