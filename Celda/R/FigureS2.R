# Figure S2: T-cell marker gene expressions

library(celda)
library(ggplot2)

sce <- readRDS("../Data/sce.rds")

markers <- c("CD3E", "CD3G",
    "CD8A", "CD8B",
    "CCR7", "SELL", "CD27", "LEF1",
    "GNLY", "KLRG1",  # NK
    "GZMA", "GZMH",
    "CD4", "IL7R",
    "MKI67", "IL2RA", "CENPF", "CENPM")

len <- length(markers)
grids <- len / 9
if (len %% 9 != 0) {
    grids <- as.integer(grids + 1)
}
glist <- vector("list", length = grids)

for (i in seq(grids)) {
    nums <- seq(((i - 1) * 9) + 1, min(i * 9, len))
    gtcell <- plotDimReduceFeature(x = sce,
        reducedDimName = "celda_UMAP",
        useAssay = "decontXcounts",
        features = markers[nums],
        size = 0.1,
        trim = c(-2, 2),
        xlab = "UMAP_1", ylab = "UMAP_2",
        normalize = TRUE)
    gtcell <- gtcell + theme(strip.text = element_text(size = 18),
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 18)) +
        theme(legend.position = "none")
    glist[[i]] <- gtcell
}

pdf("../Figures/FigureS2.pdf")
print(glist)
dev.off()
