# Figure 3: Probability heatmaps, module heatmaps and UMAPs

library(celda)
library(ComplexHeatmap)
library(ggplot2)
library(gridExtra)

sce <- readRDS("../Data/sce.rds")

panel_fun = function(index, nm = NULL) {
    topFeatures <- NULL
    showModuleLabel <- FALSE
    topAnnotationHeight <- 3
    topCells <- 50

    if (index == 10) { # B
        g <- moduleHeatmap(sce,
            useAssay = "decontXcounts",
            featureModule = index,
            topFeatures = topFeatures,
            topCells = topCells,
            showFeaturenames = FALSE,
            showModuleLabel = showModuleLabel,
            topAnnotationHeight = topAnnotationHeight,
            right_annotation = rowAnnotation(mhg = anno_mark(at = c(1:2, 4, 14),
                labels = c("CD79A", "CD79B", "MS4A1", "CD19"),
                padding = unit(1, "cm")),
                annotation_width = unit(2, "cm")))
        # width = unit(5, "cm"),
        # heatmap_width = unit(8, "cm"))

        g2 <- plotDimReduceModule(sce,
            useAssay = "decontXcounts",
            reducedDimName = "celda_UMAP",
            modules = index,
            xlab = "UMAP_1",
            ylab = "UMAP_2",
            size = 0.1) +
            xlab(NULL) +
            ylab(NULL) +
            theme(axis.line = element_line(color = "grey90"),
                panel.border = element_rect(color = "grey90"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "none",
                strip.text = element_blank())
    } else if (index == 21) { # pDC
        g <- moduleHeatmap(sce,
            useAssay = "decontXcounts",
            featureModule = index,
            topFeatures = topFeatures,
            topCells = topCells,
            showFeaturenames = FALSE,
            showModuleLabel = showModuleLabel,
            topAnnotationHeight = topAnnotationHeight,
            right_annotation = rowAnnotation(mhg = anno_mark(
                at = c(1, 2, 4, 14),
                labels = c("ITM2C", "IRF7", "LILRA4", "CLEC4C"),
                padding = unit(1, "cm")),
                annotation_width = unit(2, "cm")))
        # width = unit(5, "cm"),
        # heatmap_width = unit(8, "cm"))

        g2 <- plotDimReduceModule(sce,
            useAssay = "decontXcounts",
            reducedDimName = "celda_UMAP",
            modules = index,
            xlab = "UMAP_1",
            ylab = "UMAP_2",
            size = 0.1) +
            xlab(NULL) +
            ylab(NULL) +
            theme(axis.line = element_line(color = "grey90"),
                panel.border = element_rect(color = "grey90"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "none",
                strip.text = element_blank())
    } else if (index == 43) { # CD14+
        g <- moduleHeatmap(sce,
            useAssay = "decontXcounts",
            featureModule = index,
            topFeatures = topFeatures,
            topCells = topCells,
            showFeaturenames = FALSE,
            showModuleLabel = showModuleLabel,
            topAnnotationHeight = topAnnotationHeight,
            right_annotation = rowAnnotation(mhg = anno_mark(
                at = c(1, 2, 3, 4, 5),
                labels = c("S100A9", "S100A8", "S100A12", "VCAN", "CD14"),
                padding = unit(1, "cm")),
                annotation_width = unit(2, "cm")))

        g2 <- plotDimReduceModule(sce,
            useAssay = "decontXcounts",
            reducedDimName = "celda_UMAP",
            modules = index,
            xlab = "UMAP_1",
            ylab = "UMAP_2",
            size = 0.1) +
            xlab(NULL) +
            ylab(NULL) +
            theme(axis.line = element_line(color = "grey90"),
                panel.border = element_rect(color = "grey90"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "none",
                strip.text = element_blank())
    } else if (index == 74) { # T
        g <- moduleHeatmap(sce,
            useAssay = "decontXcounts",
            featureModule = index,
            topFeatures = topFeatures,
            topCells = topCells,
            showFeaturenames = FALSE,
            showModuleLabel = showModuleLabel,
            topAnnotationHeight = topAnnotationHeight,
            right_annotation = rowAnnotation(mhg = anno_mark(
                at = c(1, 2, 3, 5),
                labels = c("TRAC", "CD3D", "TRBC1", "CD3G"),
                padding = unit(1, "cm")),
                annotation_width = unit(2, "cm")))

        g2 <- plotDimReduceModule(sce,
            useAssay = "decontXcounts",
            reducedDimName = "celda_UMAP",
            modules = index,
            xlab = "UMAP_1",
            ylab = "UMAP_2",
            size = 0.1) +
            xlab(NULL) +
            ylab(NULL) +
            theme(axis.line = element_line(color = "grey90"),
                panel.border = element_rect(color = "grey90"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "none",
                strip.text = element_blank())
    }
    g <- grid.grabExpr(draw(g))
    g <- grid.arrange(g, g2, ncol = 2, widths = c(2, 1.8), newpage = FALSE)
    grid.draw(g)
}

collabels <- c("Plasma", rep("B", 3), rep("DC", 2), "pDC",
    "CD34+", "NK", "Mk", rep("CD14+ mono", 3), "FCGR3A+ mono",
    "Proliferating T",
    "Memory T",
    "Naive CD8+ T",
    "NKT",
    "Cytotoxic T",
    "T helper")

cellTypeCol <- distinctColors(20)
names(cellTypeCol) <- as.character(seq(20))

L <- 80
gprob <- celdaProbabilityMap(sce,
    title1 = NULL,
    title2 = NULL,
    useAssay = "decontXcounts",
    showColumnNames = FALSE,
    showRowNames = FALSE,
    width = unit(11, "cm"),
    row_title = "Gene modules",
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 22),
    row_names_side = "right",
    showHeatmapLegend = FALSE,
    heatmapLegendParam = list(title = "Probability",
        legend_width = grid::unit(6, "cm"),
        legend_direction = "horizontal",
        nrow = 1,
        ncol = 2,
        by_row = TRUE,
        title_position = "lefttop"),
    right_annotation = rowAnnotation(
        rowname = anno_text(paste0("L", seq(L)), gp = gpar(fontsize = 8))),
    top_annotation = HeatmapAnnotation(
        border = FALSE,
        gap = unit(1, "mm"),
        rowname = anno_text(collabels,
            rot = 45,
            just = "left",
            location = 0),
        `Cell cluster` = anno_text(seq(20),
            rot = 0,
            just = "bottom",
            location = 0.2,
            height = unit(5, "mm"),
            gp = gpar(col = "black",
                fill = cellTypeCol,
                lwd = 1)),
        show_legend = FALSE))

g4 <- celdaProbabilityMap(sce,
    title2 = NULL,
    useAssay = "decontXcounts",
    showColumnNames = FALSE,
    showRowNames = FALSE,
    width = unit(11, "cm"),
    row_title = "Transcriptional modules",
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 22),
    row_names_side = "right",
    showHeatmapLegend = FALSE,
    right_annotation = rowAnnotation(
        rowname = anno_text(paste0("L", seq(L)), gp = gpar(fontsize = 8)),
        modulehm = anno_zoom(align_to = as.list(c(10, 21, 43, 74)),
            panel_fun = panel_fun,
            which = "row",
            side = "right",
            #size = 5,
            width = unit(15, "cm"),
            gap = unit(3, "mm"),
            link_width = unit(1.5, "cm"),
            link_gp = gpar(col = "white",
                fill = "grey",
                alpha = 0.3)),
        gap = unit(1, "mm")),
    top_annotation = HeatmapAnnotation(
        border = FALSE,
        gap = unit(1, "mm"),
        rowname = anno_text(collabels,
            rot = 45,
            just = "left",
            location = 0),
        `Cell cluster` = anno_text(seq(20),
            rot = 0,
            just = "bottom",
            location = 0.2,
            height = unit(5, "mm"),
            gp = gpar(col = "black",
                fill = cellTypeCol,
                lwd = 1)),
        show_legend = FALSE))

pdf("../Figures/Figure3.pdf", width = 18, height = 12)
draw(gprob[, 1] + g4[, 2])
dev.off()
