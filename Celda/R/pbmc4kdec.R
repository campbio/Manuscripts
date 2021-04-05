# Load and preprocess PBMC4K data

library(TENxPBMCData)
library(celda)
library(Seurat)

pbmc4k <- TENxPBMCData(dataset = "pbmc4k")
counts(pbmc4k) <- as.matrix(counts(pbmc4k))
colnames(pbmc4k) <- colData(pbmc4k)$Sequence

pbmc4kdec <- decontX(pbmc4k)
decontXcounts(pbmc4kdec) <- as.matrix(decontXcounts(pbmc4kdec))

pbmc4kseurat <- CreateSeuratObject(counts = decontXcounts(pbmc4kdec),
    project = "pbmc4kdec", min.cells = 3, min.features = 200)

pbmc4kseurat <- NormalizeData(pbmc4kseurat)
pbmc4kseurat <- FindVariableFeatures(pbmc4kseurat)
vf1 <- VariableFeatures(pbmc4kseurat)

# pbmc4kseurat2 <- CreateSeuratObject(counts = decontXcounts(pbmc4kdec),
#     project = "pbmc4kdec", min.cells = 0, min.features = 0)
# pbmc4kseurat2 <- FindVariableFeatures(pbmc4kseurat2)
# vf2 <- VariableFeatures(pbmc4kseurat)
# all(vf1 == vf2)
#[1] TRUE

# unique gene names as rownames
rd <- rowData(pbmc4kdec)
rd$runique <- rd$Symbol_TENx
rd[duplicated(rd$Symbol_TENx), "runique"] <-
    paste0(rd[duplicated(rd$Symbol_TENx), "runique"], ".1")
rowData(pbmc4kdec) <- rd
rownames(pbmc4kdec) <- rowData(pbmc4kdec)$runique

pbmc4kf <- pbmc4kdec[rowData(pbmc4kdec)$ENSEMBL_ID %in% vf1, ]
rownames(pbmc4kf) <- rowData(pbmc4kf)$Symbol_TENx

altExp(pbmc4kdec, "featureSubset") <- pbmc4kf
saveRDS(pbmc4kdec, file = "../data/pbmc4kdec.rds")
