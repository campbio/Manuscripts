# Figure 5: Comparison of PCA and modules

library(celda)
library(Seurat)
library(SingleCellExperiment)
library(circlize)
library(data.table)
library(ComplexHeatmap)
library(ggplot2)
library(scales)
library(gridExtra)
library(cowplot)
library(gtable)

sce <- readRDS("../data/sce.rds")
L <- 80

seuratobj <- CreateSeuratObject(counts = decontXcounts(sce),
    project = "pca", min.cells = 3, min.features = 200)

seuratobj <- NormalizeData(seuratobj, verbose = TRUE)
seuratobj <- FindVariableFeatures(seuratobj, nfeatures = 2000,
    verbose = TRUE)

seuratobj <- ScaleData(seuratobj, verbose = TRUE)
seuratobj <- RunPCA(seuratobj, npcs = L, verbose = TRUE)
seuratobj <- RunUMAP(seuratobj, dims = 1:15)

ae <- altExp(sce)
seuratobj@reductions$umap@cell.embeddings[, "UMAP_1"] <-
    reducedDim(ae, "celda_UMAP")[, "celda_UMAP1"]
seuratobj@reductions$umap@cell.embeddings[, "UMAP_2"] <-
    reducedDim(ae, "celda_UMAP")[, "celda_UMAP2"]

# PC heatmap
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

cols <- Seurat::PurpleAndYellow()
col2 <- circlize::colorRamp2(seq(-2, 2, length.out = 50), cols)

pchg2 <- pchm(seuratobj, dim = 2, col = col2,
    annotation_height = unit(3, "mm"),
    showHeatmapLegend = TRUE,
    heatmap_legend_param = list(grid_width = unit(6, "mm"),
        legend_height = unit(3, "cm"),
        title = ""))

# pc umap colored by PC2 score
pc2score <- seuratobj@reductions$pca@cell.embeddings[, "PC_2"]
# all(names(pc2score) == colnames(ae))
# [1] TRUE
seuratobj@meta.data$pc2score <- pc2score

pcumapscore <- FeaturePlot(seuratobj, "pc2score", order = TRUE) +
    scale_color_gradientn(colors = c(muted("green"),
        "lightgrey",
        muted("blue")),
        values = rescale(c(min(pc2score), 0, max(pc2score))))
pcumapscorenl <- pcumapscore + NoLegend() + ggtitle(NULL) +
    theme(axis.line = element_line(color = "grey90"),
        panel.border = element_rect(color = "grey90"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title = element_blank())

# all(colnames(seuratobj) == colnames(ae))
# [1] TRUE
pc2ld <- seuratobj@reductions$pca@feature.loadings[, "PC_2"]
names(pc2ld) <- rownames(seuratobj@reductions$pca@feature.loadings)
pc2ld <- pc2ld[order(pc2ld)]
poscorgenes <- names(pc2ld[(length(pc2ld) - 14):length(pc2ld)])
negcorgenes <- names(pc2ld[1:15])

scaleddata <- seuratobj@assays$RNA@scale.data

psd <- scaleddata[poscorgenes, ]
avgpsd <- apply(psd, 2, mean)

nsd <- scaleddata[negcorgenes, ]
avgnsd <- apply(nsd, 2, mean)

seuratobj@meta.data$avgpsd <- avgpsd
seuratobj@meta.data$avgnsd <- avgnsd

# Ppositively correlated genes
pcumappos <- FeaturePlot(seuratobj, "avgpsd", order = TRUE) +
    scale_color_gradientn(colors = c(muted("green"),
        "lightgrey",
        muted("blue")),
        values = rescale(c(min(avgpsd), 0, max(avgpsd))))
pcumapposnl <- pcumappos + NoLegend() + ggtitle(NULL) +
    theme(axis.line = element_line(color = "grey90"),
        panel.border = element_rect(color = "grey90"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title = element_blank())

# Negatively correlated genes
pcumapneg <- FeaturePlot(seuratobj, "avgnsd", order = TRUE) +
    scale_color_gradientn(colors = c(muted("green"),
        "lightgrey",
        muted("blue")),
        values = rescale(c(min(avgnsd), 0, max(avgnsd))))
pcumapnegnl <- pcumapneg + NoLegend() + ggtitle(NULL) +
    theme(axis.line = element_line(color = "grey90"),
        panel.border = element_rect(color = "grey90"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title = element_blank())

# celda module heatmaps and umaps
mhf <- function(sce,
    featureModule,
    useAssay = "decontXcounts",
    topCells = 50,
    rowFontSize = 5,
    showFeaturenames = FALSE,
    showModuleLabel = FALSE,
    topAnnotationHeight = 3,
    row_names_max_width = unit(7, "mm"),
    ...) {

    return(moduleHeatmap(sce,
        useAssay = useAssay,
        featureModule = featureModule,
        topCells = topCells,
        rowFontSize = rowFontSize,
        showFeaturenames = showFeaturenames,
        showModuleLabel = showModuleLabel,
        topAnnotationHeight = topAnnotationHeight,
        row_names_max_width = row_names_max_width,
        ...))
}

# modules of positively correlated genes
# S100A6 CST7 SRGN GZMA CCL5
# IFITM1 ANXA1 NKG7 CTSW CD7
# CD3D TRAC TMSB4X IL32 S100A4

posmodules <- c(49, 40, 34, 68, 37, 32, 39, 74, 78, 72, 51)
moduleumapposlist <- vector("list", length = length(posmodules))
moduleheatmapposlist <- vector("list", length = length(posmodules))

for (i in seq_along(posmodules)) {
    moduleumapposlist[[i]] <- plotDimReduceModule(sce,
        useAssay = "decontXcounts",
        reducedDimName = "celda_UMAP",
        modules = posmodules[i],
        size = 0.1) +
        theme(axis.line = element_line(color = "grey90"),
            panel.border = element_rect(color = "grey90"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            strip.text = element_blank(),
            axis.title = element_blank())
}

for (i in seq_along(posmodules)) {
    moduleheatmapposlist[[i]] <-  mhf(sce = sce,
        featureModule = posmodules[[i]])
}

rawidth <- 0.2
padding <- 7
unit <- "mm"

# heatmap annotations for positive 15 genes
hmposannot1 = rowAnnotation(foo = anno_mark(at = 1,
    padding = padding,
    labels = c("S100A6")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot2 = rowAnnotation(foo = anno_mark(at = c(1, 2, 3),
    padding = padding,
    labels = c("NKG7", "GZMA", "CST7")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot3 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("SRGN")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot4 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("CCL5")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot5 = rowAnnotation(foo = anno_mark(at = c(1, 2),
    padding = padding,
    labels = c("CD7", "IFITM1")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot6 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("ANXA1")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot7 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("CTSW")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot8 = rowAnnotation(foo = anno_mark(at = c(1, 2),
    padding = padding,
    labels = c("TRAC", "CD3D")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot9 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("TMSB4X")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot10 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("IL32")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmposannot11 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("S100A4")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))

hmposannotlist <- list(hmposannot1, hmposannot2, hmposannot3, hmposannot4,
    hmposannot5, hmposannot6, hmposannot7, hmposannot8,
    hmposannot9, hmposannot10, hmposannot11)

# modules of negatively correlated genes
negmodules <- c(10, 5, 12, 6, 14, 7)
moduleumapneglist <- vector("list", length = length(negmodules))
moduleheatmapneglist <- vector("list", length = length(negmodules))


for (i in seq_along(negmodules)) {
    moduleumapneglist[[i]] <- plotDimReduceModule(sce,
        useAssay = "decontXcounts",
        reducedDimName = "celda_UMAP",
        modules = negmodules[i],
        size = 0.1) +
        theme(axis.line = element_line(color = "grey90"),
            panel.border = element_rect(color = "grey90"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            strip.text = element_blank(),
            axis.title = element_blank())
}

for (i in seq_along(negmodules)) {
    moduleheatmapneglist[[i]] <- mhf(sce = sce,
        featureModule = negmodules[[i]])
}

hmnegannot1 = rowAnnotation(foo = anno_mark(at = seq(7),
    padding = padding,
    labels = c("CD79A", "CD79B", "IGHM", "MS4A1",
        "IGHD", "LINC00926", "CD22")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmnegannot2 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("IGKC")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmnegannot3 = rowAnnotation(foo = anno_mark(at = c(2, 1),
    padding = padding,
    labels = c("BANK1", "TCL1A")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmnegannot4 = rowAnnotation(foo = anno_mark(at = c(1),
    padding = padding,
    labels = c("CD74")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmnegannot5 = rowAnnotation(foo = anno_mark(at = c(2, 1),
    padding = padding,
    labels = c("HLA-DPA1", "HLA-DPB1")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))
hmnegannot6 = rowAnnotation(foo = anno_mark(at = c(3, 1),
    padding = padding,
    labels = c("HLA-DQB1", "HLA-DRA")),
    #padding = unit(1, unit),
    width = unit(rawidth, unit))

hmnegannotlist <- list(hmnegannot1, hmnegannot2, hmnegannot3, hmnegannot4,
    hmnegannot5, hmnegannot6)

# combine plots

# PC heatmap and UMAPs
pcgrob <- arrangeGrob(grobs = list(grid.grabExpr(draw(pchg2)),
    pcumapscorenl,
    pcumapnegnl,
    pcumapposnl),
    nrow = 2,
    ncol = 2,
    newpage = FALSE)

pcgrobpad <- arrangeGrob(grobs = lapply(
    list(as_gtable(grid.grabExpr(draw(pchg2))),
        as_gtable(pcumapscorenl),
        as_gtable(pcumapnegnl),
        as_gtable(pcumapposnl)),
    gtable_add_padding, unit(15, "mm")),
    nrow = 2,
    ncol = 2,
    #widths = c(2, 1.8),
    newpage = FALSE)

# module heatmaps and UMAPs

modulegrobposlist <- vector("list", length = 6)
modulegrobposlist2 <- vector("list", length = 5)

# first 6 of moduleheatmapposlist
for (i in seq(6)) {
    modulegrobposlist[[i]] <- arrangeGrob(grobs = list(grid.grabExpr(
        draw(moduleheatmapposlist[[i]] + hmposannotlist[[i]])),
        moduleumapposlist[[i]]),
        nrow = 2,
        ncol = 1,
        newpage = FALSE)
}

# last 6 of moduleheatmapposlist
for (i in seq(7, 11)) {
    modulegrobposlist2[[i - 6]] <- arrangeGrob(grobs = list(grid.grabExpr(
        draw(moduleheatmapposlist[[i]] + hmposannotlist[[i]])),
        moduleumapposlist[[i]]),
        nrow = 2,
        ncol = 1,
        newpage = FALSE)
}

modulegrobneglist <- vector("list", length = length(moduleheatmapneglist))

# moduleheatmapneglist
for (i in seq_along(moduleheatmapneglist)) {
    modulegrobneglist[[i]] <- arrangeGrob(grobs = list(grid.grabExpr(
        draw(moduleheatmapneglist[[i]] + hmnegannotlist[[i]])),
        moduleumapneglist[[i]]),
        nrow = 2,
        ncol = 1,
        newpage = FALSE)
}

margin <- unit(10, "mm")

garr1 <- arrangeGrob(grobs = lapply(modulegrobneglist,
    gtable_add_padding, margin),
    nrow = 1,
    ncol = 6,
    #widths = c(2, 1.8),
    newpage = FALSE)
garr2 <- arrangeGrob(grobs = lapply(modulegrobposlist,
    gtable_add_padding, margin),
    nrow = 1,
    ncol = 6,
    newpage = FALSE)
garr3 <- arrangeGrob(grobs = lapply(modulegrobposlist2,
    gtable_add_padding, margin),
    nrow = 1,
    ncol = 6,
    #widths = c(2, 1.8),
    newpage = FALSE)

moduleplots <- arrangeGrob(grobs = list(garr1, garr2, garr3),
    nrow = 3,
    ncol = 1,
    #widths = c(2, 1.8),
    newpage = FALSE)

figure5 <- arrangeGrob(grobs = list(pcgrobpad, moduleplots),
    nrow = 2,
    ncol = 1,
    heights = c(1.5, 2),
    newpage = FALSE)

pdf("../results/Figure5.pdf", width = 20, height = 20 * 4 / 3)
grid.draw(figure5)
dev.off()
