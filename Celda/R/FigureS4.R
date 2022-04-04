
library(Matrix)
library(TENxPBMCData)
library(SingleCellExperiment)
library(celda)
library(singleCellTK)
library(Seurat)
library(ggplot2)
library(data.table)

pbmc33k <- TENxPBMCData("pbmc33k")

colnames(pbmc33k) <- colData(pbmc33k)$Barcode
counts(pbmc33k) <- as(counts(pbmc33k), "dgCMatrix")

pbmc33k <- decontX(pbmc33k)

decUMAP <- reducedDim(pbmc33k, "decontX_UMAP")
altExp(pbmc33k, "allFeatures") <- pbmc33k

pbmc33k <- seuratFindHVG(
    pbmc33k,
    useAssay = "decontXcounts",
    hvgMethod = "vst")

o <- head(order(
    rowData(pbmc33k)$seurat_variableFeatures_vst_varianceStandardized,
    decreasing = TRUE),
    n = 2000)
vst2000 <- pbmc33k[o, ]
altExp(vst2000, "allFeatures") <- NULL
altExp(pbmc33k, "VST2000") <- vst2000

pbmc33k <- recursiveSplitModule(pbmc33k,
    useAssay = "decontXcounts",
    altExpName = "VST2000",
    initialL = 3,
    maxL = 150)

sce1 <- pbmc33k

g1 <- plotGridSearchPerplexity(sce1, altExpName = "VST2000", sep = 10)
g1 <- g1 + scale_x_discrete(breaks = seq(0, 150, 10),
    limits = factor(seq(0, 152, 1))) +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

g2 <- plotRPC(sce1, altExpName = "VST2000", sep = 10, n = 20)
g2 <- g2 + scale_x_discrete(breaks = seq(0, 150, 10),
    limits = factor(seq(0, 152, 1))) +
    ylab("RPC of L") +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

L <- 80

pbmc33kL80 <- subsetCeldaList(pbmc33k, list(L = L), altExpName = "VST2000")

# ~ 15 minutes
pbmc33k <- recursiveSplitCell(pbmc33kL80,
    useAssay = "decontXcounts",
    altExpName = "VST2000",
    initialK = 3,
    maxK = 30,
    yInit = celdaModules(pbmc33kL80, altExpName = "VST2000"))

sce2 <- pbmc33k

g3 <- plotGridSearchPerplexity(sce2, altExpName = "VST2000", sep = 5)
g3 <- g3 + scale_x_discrete(breaks = seq(0, 30, 5),
    limits = factor(seq(0, 31, 1))) +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")

g4 <- plotRPC(sce2, altExpName = "VST2000", sep = 5, n = 3)
g4 <- g4 + scale_x_discrete(breaks = seq(0, 30, 5),
    limits = factor(seq(0, 31, 1))) +
    ylab("RPC of K") +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")

K <- 20
pbmc33k <- subsetCeldaList(pbmc33k, list(K = K), altExpName = "VST2000")

assay(pbmc33k, "normcounts") <- NormalizeData(decontXcounts(pbmc33k),
    normalization.method = "RC",
    scale.factor = 1)

pbmc33k <- celdaUmap(
    pbmc33k,
    useAssay = "decontXcounts",
    altExpName = "VST2000")

#saveRDS(pbmc33k, file = "../Data/pbmc33k_umap.rds")

selectedMarkers <- c("CD3D", # T cells
    "CD4",
    "CD8A",
    "CD19", # B
    "NKG7",
    #"KLRD1",
    "FCGR3A", # FCGR3A+ monocytes
    "CLEC4C",
    "CD14", # CD14+ monocytes
    #"MKI67", # proliferating T
    "IGJ",
    #"IGLL5",
    "FCER1A", # DC
    "CD34",
    "ITGA2B") # Mk

g5 <- plotDimReduceFeature(x = pbmc33k,
    features = selectedMarkers,
    displayName = "Symbol_TENx",
    reducedDimName = "celda_UMAP",
    useAssay = "normcounts",
    altExpName = "VST2000",
    size = 0.1,
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    normalize = FALSE)
g5 <- g5 + theme(strip.text = element_text(size = 18),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18))


labels <- c("1: CD14+ Mono",
    "2: CD14+ Mono",
    "3: CD14+ Mono",
    "4: CD14+ Mono",
    "5: Mk",
    "6: FCGR3A+ Mono",
    "7: FCGR3A+ Mono",
    "8: pDC & CD34+",
    "9: DC",
    "10: B",
    "11: B",
    "12: Mk",
    "13: NK",
    "14: NK",
    "15: NKT",
    "16: CD8+ T",
    "17: T",
    "18: CD4+ T",
    "19: Plasma cell",
    "20: Plasma cell")

g6 <- plotDimReduceCluster(
    pbmc33k,
    reducedDimName = "celda_UMAP",
    altExpName = "VST2000",
    labelClusters = TRUE,
    size = 0.1,
    xlab = "UMAP_1",
    ylab = "UMAP_2")
g6 <- g6 + scale_color_manual(labels = labels,
    values = distinctColors(length(labels))) +
    theme(axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)) +
    guides(colour = guide_legend(override.aes = list(size = 4)))

pdf("../results/FigureS4.pdf")
print(g6)
print(g5)
print(g1)
print(g2)
print(g3)
print(g4)
dev.off()

