# Figure 4: Module UMAPs

library(celda)
library(grid)
library(gridExtra)
library(ggplot2)

sce <- readRDS("../data/sce.rds")

m2 <- c(12, 44, 40, 65)
m3 <- c(15, 45, 47, 75)
m4 <- c(7, 14, 6, 33)
m5 <- c(24, 62, 28, 80)

pf <- function(sce, module, size = 0.1) {
    g <- plotDimReduceModule(sce,
        useAssay = "decontXcounts",
        reducedDimName = "celda_UMAP",
        modules = module,
        xlab = "UMAP_1",
        ylab = "UMAP_2",
        size = size) +
        xlab(NULL) +
        ylab(NULL) +
        theme(axis.line = element_line(color = "grey90"),
            panel.border = element_rect(color = "grey90"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            #strip.text = element_blank(),
            legend.position = "none")
    return(g)
}


g21 <- pf(sce, m2[1])
g22 <- pf(sce, m2[2])
g23 <- pf(sce, m2[3])
g24 <- pf(sce, m2[4])

g31 <- pf(sce, m3[1])
g32 <- pf(sce, m3[2])
g33 <- pf(sce, m3[3])
g34 <- pf(sce, m3[4])

g41 <- pf(sce, m4[1])
g42 <- pf(sce, m4[2])
g43 <- pf(sce, m4[3])
g44 <- pf(sce, m4[4])

g51 <- pf(sce, m5[1])
g52 <- pf(sce, m5[2])
g53 <- pf(sce, m5[3])
g54 <- pf(sce, m5[4])

glist <- list(g21, g22, g23, g24,
    g31, g32, g33, g34,
    g41, g42, g43, g44,
    g51, g52, g53, g54)

margin = theme(plot.margin = unit(c(2,0,2,0), "mm"))
garr2 <- arrangeGrob(grobs = lapply(glist, "+", margin),
    nrow = 4,
    ncol = 4,
    #widths = c(2, 1.8),
    newpage = FALSE)

pdf("../results/Figure4.pdf")
grid.draw(garr2)
dev.off()

