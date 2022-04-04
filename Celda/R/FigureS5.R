
library(TENxPBMCData)
library(Matrix)
library(SingleCellExperiment)
library(celda)
library(singleCellTK)
library(Seurat)
library(ggplot2)
library(data.table)

pbmc68k <- TENxPBMCData("pbmc68k")

colnames(pbmc68k) <- colData(pbmc68k)$Barcode
counts(pbmc68k) <- as(counts(pbmc68k), "dgCMatrix")

pbmc68k <- decontX(pbmc68k)

decUMAP <- reducedDim(pbmc68k, "decontX_UMAP")
altExp(pbmc68k, "allFeatures") <- pbmc68k

pbmc68k <- seuratFindHVG(
    pbmc68k,
    useAssay = "decontXcounts",
    hvgMethod = "vst")

o <- head(order(
    rowData(pbmc68k)$seurat_variableFeatures_vst_varianceStandardized,
    decreasing = TRUE),
    n = 2000)
vst2000 <- pbmc68k[o, ]
altExp(vst2000, "allFeatures") <- NULL
altExp(pbmc68k, "VST2000") <- vst2000

pbmc68k <- recursiveSplitModule(pbmc68k,
    useAssay = "decontXcounts",
    altExpName = "VST2000",
    initialL = 3,
    maxL = 100)

sce1 <- pbmc68k
g1 <- plotGridSearchPerplexity(sce1, altExpName = "VST2000", sep = 10)
g1 <- g1 + scale_x_discrete(breaks = seq(0, 100, 10),
    limits = factor(seq(0, 102, 1))) +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")
g2 <- plotRPC(sce1, altExpName = "VST2000", sep = 10, n = 10)
g2 <- g2 + scale_x_discrete(breaks = seq(0, 100, 10),
    limits = factor(seq(0, 102, 1))) +
    ylab("RPC of L") +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")

L <- 80
pbmc68kL80 <- subsetCeldaList(pbmc68k, list(L = L), altExpName = "VST2000")

pbmc68k <- recursiveSplitCell(pbmc68kL80,
    useAssay = "decontXcounts",
    altExpName = "VST2000",
    initialK = 3,
    maxK = 30,
    yInit = celdaModules(pbmc68kL80, altExpName = "VST2000"))

sce2 <- pbmc68k

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
pbmc68k <- subsetCeldaList(sce2, list(K = K), altExpName = "VST2000")

assay(pbmc68k, "normcounts") <- NormalizeData(decontXcounts(pbmc68k),
    normalization.method = "RC",
    scale.factor = 1)

pbmc68k <- celdaUmap(
    pbmc68k,
    useAssay = "decontXcounts",
    altExpName = "VST2000")

#saveRDS(pbmc68k, file = "../Data/pbmc68k_umap.rds")

selectedMarkers <- c("CD3D", # T cells
    "CD4",
    "CD8A",
    "CD19", # B
    "NKG7",
    #"KLRD1",
    "GNLY",
    "FCGR3A", # FCGR3A+ monocytes
    "CLEC4C",
    "CD14", # CD14+ monocytes
    #"MKI67", # proliferating T
    "IGJ",
    #"IGLL5",
    "FCER1A", # DC
    #"CD34",
    "ITGA2B") # Mk

g5 <- plotDimReduceFeature(x = pbmc68k,
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

labels <- c("1: Mk",
    "2: FCGR3A+ Mono",
    "3: FCGR3A+ Mono",
    "4: CD14+ Mono",
    "5: CD14+ Mono",
    "6: B",
    "7: B",
    "8: B",
    "9: DC",
    "10: pDC",
    "11: NK",
    "12: Plasma cell",
    "13: NKT",
    "14: T",
    "15: CD4+ T",
    "16: CD4+ T",
    "17: CD8+ T",
    "18: CD4+ T",
    "19: NKT",
    "20: CD4+ T")

g6 <- plotDimReduceCluster(
    pbmc68k,
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

pdf("../results/FigureS5.pdf")
print(g6)
print(g5)
print(g1)
print(g2)
print(g3)
print(g4)
dev.off()
