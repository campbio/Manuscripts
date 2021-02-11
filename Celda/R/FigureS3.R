# Figure S3: Seurat and scran clusterings

library(celda)
library(Seurat)
library(scran)
library(scater)

sce <- readRDS("../Data/sce.rds")

K <- 20
L <- 80

# Seurat
pbmc4kseurat <- CreateSeuratObject(counts = decontXcounts(sce),
    project = "pbmc4kdec", min.cells = 3, min.features = 200)
dim(pbmc4kseurat)

pbmc4kseurat <- NormalizeData(pbmc4kseurat)
pbmc4kseurat <- FindVariableFeatures(pbmc4kseurat)
pbmc4kseurat <- ScaleData(pbmc4kseurat)
pbmc4kseurat <- RunPCA(pbmc4kseurat)
#ElbowPlot(pbmc4kseurat, ndims = 50)
pbmc4kseurat <- FindNeighbors(pbmc4kseurat, dims = seq(22))
pbmc4kseurat <- FindClusters(pbmc4kseurat)
pbmc4kseurat <- RunUMAP(pbmc4kseurat, dims = seq(22),
    n.neighbors = 10, min.dist = 0.5)


g1 <- DimPlot(pbmc4kseurat, reduction = "umap", label = TRUE, pt.size = 0.7) +
    theme(legend.position = "none")

markers <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMK", #CD8+ T
    "MS4A1", "CD19", #B
    "CD4", "CCR7", "CD27", "SELL", # Naive CD4+ T
    "IL7R", "S100A4", # Memory CD4+ T
    "NKG7", "GZMB", "GNLY", # NK
    "FCGR3A", # FCGR3A+ mono
    "CD14", "LYZ", "S100A9", # CD14+ mono
    "MKI67", "CENPF", "CENPM", # Activated T
    "FCER1A", "HLA-DQA1", "HLA-DPB1", "HLA-DRB1", "CST3", # DC
    "CLEC4C", "PLAC8", "IRF7", "IRF8", # pDC
    "PPBP", "ITGA2B", "CXCR4", # Mk
    "CD34", # CD34+
    "MT-CO1", "MT-CO2", "MT-CO3", "IGKC",
    "IGHG1", "IGHG2", "IGHG3", # pc
    "IGLC2", "IGLC3")

len <- length(markers)
grids <- len / 9
if (len %% 9 != 0) {
    grids <- as.integer(grids + 1)
}
for (i in seq(grids)) {
    nums <- seq(((i - 1) * 9) + 1, min(i * 9, len))
    g <- FeaturePlot(pbmc4kseurat, features = markers[nums],
        reduction = "umap", order = TRUE, pt.size = 0.1)
    print(g)
}

# scran
set.seed(12345)
ae <- altExp(sce)
assay(ae, "counts") <- NULL

ae <- logNormCounts(ae, assay.type = "decontXcounts")
ae <- runPCA(ae, ntop = 2000) # no scaling

output <- getClusteredPCs(reducedDim(ae, "PCA"))
npcs <- metadata(output)$chosen
reducedDim(ae, "PCAsub") <- reducedDim(ae, "PCA")[, seq(npcs), drop=FALSE]
# npcs # 30

graph <- buildSNNGraph(ae, use.dimred = "PCAsub")
cluster <- igraph::cluster_walktrap(graph)$membership

colLabels(ae) <- factor(cluster)
table(colLabels(ae))
length(unique(colLabels(ae)))
# 24

ae <- runUMAP(ae, dimred = "PCAsub", n_neighbors = 5)

g2 <- plotUMAP(ae, colour_by = "label", text_by = "label", point_alpha = 1) +
    scale_color_manual(values = distinctColors(length(unique(colLabels(ae))))) +
    theme(legend.position = "none")

pdf("../Figures/FigureS3.pdf")
print(g1)
print(g2)
dev.off()


markers <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMK", #CD8+ T
    "MS4A1", "CD19", #B
    "CD4", "CCR7", "CD27", "SELL", # Naive CD4+ T
    "IL7R", "S100A4", # Memory CD4+ T
    "NKG7", "GZMB", "GNLY", # NK
    "FCGR3A", # FCGR3A+ mono
    "CD14", "LYZ", "S100A9", # CD14+ mono
    "MKI67", "CENPF", "CENPM", # Activated T
    "FCER1A", "HLA-DQA1", "HLA-DPB1", "HLA-DRB1", "CST3", # DC
    "CLEC4C", "PLAC8", "IRF7", "IRF8", # pDC
    "PPBP", "ITGA2B", "CXCR4", # Mk
    "CD34", # CD34+
    "MT-CO1", "MT-CO2", "MT-CO3", "IGKC",
    "IGHG1", "IGHG2", "IGHG3", # pc
    "IGLC2", "IGLC3")


len <- length(markers)
grids <- len / 9
if (len %% 9 != 0) {
    grids <- as.integer(grids + 1)
}

altExp(sce) <- ae
for (i in seq(grids)) {
    nums <- seq(((i - 1) * 9) + 1, min(i * 9, len))
    g <- plotDimReduceFeature(sce,
        altExpName = "featureSubset",
        useAssay = "decontXcounts",
        reducedDimName = "UMAP",
        features = markers[nums],
        decreasing = FALSE,
        size = 0.1)
    print(g)
}





