
library(celda)
library(SingleCellExperiment)
library(data.table)
library(ggplot2)

sce <- readRDS("../../Data/sce.rds")
cc <- readRDS("../../Data/qubic2_cell_clusters.rds")

alt <- altExp(sce)
dt <- data.table(reducedDim(alt, "celda_UMAP"), keep.rownames = TRUE)

for (i in seq(length(cc))) {
    dt[, paste0("BC", i)] <- 0
    dt[rn %in% cc[[i]], ][, paste0("BC", i)] <- 1
}

dtm <- melt(dt, id.vars = c("rn", "celda_UMAP1", "celda_UMAP2"))
dtm <- dtm[order(value), ]
dtm[, value := factor(value)]
cols <- c("grey90", "firebrick1")

l <- 4
glist <- vector("list", length = l)

for (i in seq(l)) {
    dti <- dtm[variable %in% paste0("BC", seq((i - 1) * 20 + 1, i * 20)), ]

    glist[[i]] <- ggplot(data = dti) +
        geom_point(aes(x = celda_UMAP1, y = celda_UMAP2, col = value),
            size = 0.01) +
        facet_wrap(vars(variable)) +
        ggplot2::theme_bw() +
        theme(strip.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.spacing = unit(0, "lines"),
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            strip.text = element_text(size = 18),
            legend.position = "none") +
        scale_color_manual(values = cols)
}

pdf("../../results/FigureS7.pdf")
print(glist)
dev.off()
