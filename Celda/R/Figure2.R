# Figure 2
library(celda)
library(ggplot2)

sce <- readRDS("../data/sce.rds")

g3.1 <- plotDimReduceCluster(sce, reducedDimName = "celda_UMAP",
    labelClusters = FALSE, xlab = "UMAP_1", ylab = "UMAP_2") +
    theme(axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)) +
    guides(colour = guide_legend(override.aes = list(size = 4)))

labels <- c("1: Plasma",
    paste0(seq(2, 4), rep(": B", 3)),
    paste0(seq(5, 6), rep(": DC", 2)),
    "7: pDC",
    "8: CD34+",
    "9: NK",
    "10: Mk",
    paste0(seq(11, 13), ": ", rep("CD14+ mono", 3)),
    "14: FCGR3A+ mono",
    "15: Proliferating T",
    "16: T memory",
    "17: T cytotoxic",
    "18: NKT",
    "19: NKT",
    "20: T")

g3.2 <- g3.1 +
    scale_color_manual(labels = labels, values = distinctColors(20))
print(g3.2)

selectedMarkers <- c("CD3D", # leukocytes
    "CD19", # B
    "KLRD1", # NK "NCAM1", "KLRD1", "TBX21", "EOMES", "KLRC1", "NKG7", "GNLY"
    "FCGR3A", # FCGR3A+ monocytes
    "CD14", # CD14+ monocytes
    #"MKI67", # proliferating T
    "FCER1A", # DC
    "CLEC4C", # pDC
    "ITGA2B", # Mk
    "CD34") # CD34+ HSC/progenitors

g4.1 <- plotDimReduceFeature(x = sce,
    reducedDimName = "celda_UMAP",
    useAssay = "decontXcounts",
    features = selectedMarkers, size = 0.3, trim = c(-2, 2),
    xlab = "UMAP_1", ylab = "UMAP_2", #colorMid = "floralwhite",
    normalize = TRUE)
g4.2 <- g4.1 + theme(strip.text = element_text(size = 18),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18))
print(g4.2)


pdf("../results/Figure2.pdf", width = 10)
print(g3.2)
print(g4.2)
dev.off()
