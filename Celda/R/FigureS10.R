# Figure S10: Heatmaps of first 9 PCs

library(celda)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(data.table)
library(patchwork)
library(ComplexHeatmap)
library(gridExtra)
library(scater)

sce <- readRDS("../data/sce.rds")

L <- 80

colorLow = "blue4"
colorMid = "grey90"
colorHigh = "firebrick1"
limits = c(-2, 2)
midpoint = 0

seuratobj <- CreateSeuratObject(counts = decontXcounts(sce),
    project = "pca", min.cells = 3, min.features = 200)

seuratobj <- NormalizeData(seuratobj, verbose = TRUE)
seuratobj <- FindVariableFeatures(seuratobj, nfeatures = 2000,
    verbose = TRUE)

seuratobj <- ScaleData(seuratobj, verbose = TRUE)
seuratobj <- RunPCA(seuratobj, npcs = L, verbose = TRUE)


pchm <- function(seuratobj,
    dim,
    ncells = 100,
    nfeatures = 30,
    balanced = TRUE,
    annotation_height = unit(1.5, "mm"),
    fontsize = 12,
    showHeatmapLegend = FALSE,
    showAnnoLegend = FALSE,
    showColumnTitle = FALSE,
    col = circlize::colorRamp2(c(-2, 0, 2),
        c("#1E90FF", "#FFFFFF", "#CD2626")),
    ...) {


    cells <- TopCells(
        object = seuratobj[["pca"]],
        dim = dim,
        ncells = ncells,
        balanced = balanced)

    if (balanced) {
        cells$negative <- rev(x = cells$negative)
    }
    cells <- unlist(x = unname(obj = cells))

    tfea <- TopFeatures(seuratobj[["pca"]],
        dim = dim,
        nfeatures = nfeatures,
        balanced = balanced)
    tfea <- unique(c(tfea$negative, rev(tfea$positive)))

    datascale <- seuratobj@assays$RNA@scale.data

    dts <- datascale[tfea, cells]

    cols <- distinctColors(20)
    dtannot <- data.table(cellid = cells,
        celda_cluster = celdaClusters(sce[, cells]))
    dtannot[, color := cols[celda_cluster]]
    colmap <- c(dtannot$color)
    names(colmap) <- dtannot$celda_cluster

    hm1 <- Heatmap(matrix = dts,
        col = col,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = showHeatmapLegend,
        row_names_gp = gpar(fontsize = fontsize),
        top_annotation = HeatmapAnnotation(cellcluster = dtannot$celda_cluster,
            col = list(cellcluster = colmap),
            show_legend = showAnnoLegend,
            show_annotation_name = FALSE,
            simple_anno_size = annotation_height),
        column_title = ifelse(showColumnTitle,
            paste0("PC ", dim), ""),
        ...)
    return(hm1)
}


pchmdims <- function(seuratobj,
    dims,
    nrow,
    ncol,
    ncells = 100,
    nfeatures = 30,
    balanced = TRUE,
    annotation_height = unit(1.5, "mm"),
    fontsize = 12,
    gap = unit(2, "mm"),
    padding = unit(rep(5.5, 4), "points"),
    showColumnTitle = FALSE,
    ...) {

    pchm <- lapply(X = dims,
        FUN = pchm,
        seuratobj = seuratobj,
        #dim = dims,
        ncells = ncells,
        nfeatures = nfeatures,
        balanced = balanced,
        annotation_height = annotation_height,
        fontsize = fontsize,
        showColumnTitle = showColumnTitle,
        ... = ...)

    glist <- lapply(pchm, function(x) {
        grid.grabExpr(draw(x, gap = gap, padding = padding))
    })

    g <- arrangeGrob(grobs = glist, nrow = nrow,
        ncol = ncol)
    return(g)
}

fontsize <- 5
padlr <- 17


cols <- Seurat::PurpleAndYellow()
col2 <- circlize::colorRamp2(seq(-2, 2, length.out = 50), cols)

g1 <- pchmdims(seuratobj, seq(9), nrow = 3, ncol = 3,
    #width = 10,
    row_names_max_width = unit(7, "mm"),
    fontsize = fontsize,
    padding = unit(c(6, padlr, 6, padlr), "points"),
    col = col2)

pdf("../results/FigureS4.pdf")
grid.draw(g1)
dev.off()


# Percent variance explained
sce2 <- SingleCellExperiment(assays = list(
    counts = seuratobj@assays$RNA@scale.data))
scaterpca <- scater::runPCA(sce2, ntop = 2000, exprs_values = "counts",
    ncomponents = 80)
red <- reducedDim(scaterpca, "PCA")
str(red)
attr(red, "percentVar")[1:9]
# [1] 8.8633114 3.2198463 1.9634761 1.4881054 0.9732545 0.8353006 0.7893381
# [8] 0.5367381 0.4450517

