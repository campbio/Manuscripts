# Figure S3: Marker gene expressions in plasma cell and B cells

library(celda)
library(SingleCellExperiment)
library(ggplot2)
library(data.table)

sce <- readRDS("../Data/sce.rds")

.themePublication <- function(base_size = 12,
    base_family = "sans") {
    (ggthemes::theme_foundation(base_size = base_size,
        base_family = base_family) +
            ggplot2::theme(plot.title = ggplot2::element_text(
                face = "bold",
                size = ggplot2::rel(1),
                hjust = 0.5),
                text = ggplot2::element_text(),
                panel.background = ggplot2::element_rect(color = NA),
                plot.background = ggplot2::element_rect(color = NA),
                panel.border = ggplot2::element_rect(color = NA),
                axis.title = ggplot2::element_text(
                    face = "bold",
                    size = ggplot2::rel(1)),
                axis.title.y = ggplot2::element_text(angle = 90,
                    vjust = 2),
                axis.title.x = ggplot2::element_text(vjust = -0.2),
                axis.text = ggplot2::element_text(),
                axis.line = ggplot2::element_line(color = "black"),
                axis.ticks = ggplot2::element_line(),
                panel.grid.major = ggplot2::element_line(color = "#f0f0f0"),
                panel.grid.minor = ggplot2::element_blank(),
                legend.key = ggplot2::element_rect(color = NA),
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.size = ggplot2::unit(0.5, "cm"),
                legend.margin = ggplot2::margin(0),
                legend.title = ggplot2::element_text(face = "bold"),
                plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
                strip.background = ggplot2::element_rect(
                    color = "#f0f0f0", fill = "#f0f0f0"),
                strip.text = ggplot2::element_text(face = "bold")
            ))
}

de1 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(15, 16, 17, 18, 19, 20))

de1[Gene %in% c("CD3D", "CD3E", "CD3G"), ]
#    Gene Pvalue  Log2_FC FDR
# 1: CD3D      0 6.200828   0
# 2: CD3E      0 5.565626   0
# 3: CD3G      0 4.046326   0

de2 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(17, 18, 19),
    c2 = c(15, 16, 20))

de2[Gene %in% c("CD8A", "CD8B"), ]
#    Gene        Pvalue  Log2_FC           FDR
# 1: CD8B 5.910581e-176 3.287508 4.978778e-172
# 2: CD8A 7.534401e-169 3.024783 5.077282e-165

de3 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(17),
    c2 = c(18, 19))

de3[Gene %in% c("CCR7"), ]
#    Gene       Pvalue  Log2_FC          FDR
# 1: CCR7 3.415604e-88 3.279663 3.110415e-85

de4 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(18),
    c2 = c(17, 19))

de4[Gene %in% c("GNLY", "KLRG1", "GZMA", "GZMH"), ]
#     Gene        Pvalue  Log2_FC           FDR
# 1:  GZMA 5.602638e-104 4.657639 3.775506e-100
# 2:  GZMH  5.513208e-88 3.599870  2.322025e-84
# 3: KLRG1  1.439579e-48 3.012203  1.054460e-45
# 4:  GNLY  2.251357e-34 2.259016  8.719222e-32

de5 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(1),
    c2 = c(2, 3, 4))

de5[Gene %in% c("IGHG1", "IGHG3", "IGLC2"), ]
#     Gene       Pvalue  Log2_FC          FDR
# 1: IGHG1 8.906256e-11 1.583893 2.326259e-08
# 2: IGHG3 1.030868e-10 1.483946 2.651455e-08
# 3: IGLC2 5.638171e-04 3.802954 1.860652e-02

plasmac <- which(celdaClusters(sce) == 1)
# [1] 2904

sum(counts(sce)[, plasmac])
# [1] 48443

bcells <- which(celdaClusters(sce) %in% c(2, 3, 4))

pcgenes <- c("CD79A", "CD79B", "CD19", "IGKC", "IGHG1", "IGHG3", "IGLC2",
    "IGLC3")
# plasma cell UMI counts
counts(sce)[pcgenes, plasmac]
# CD79A CD79B  CD19  IGKC IGHG1 IGHG3 IGLC2 IGLC3
#    10     8     2     2  1982  1409  8029  1480

# mean B-cell UMI counts
apply(counts(sce)[pcgenes, bcells], MARGIN = 1, mean)
#     CD79A      CD79B       CD19       IGKC      IGHG1      IGHG3      IGLC2
# 5.9934747  5.5644372  0.5791191 23.3964111  0.6084829  0.2805873  7.0864600
#     IGLC3
# 6.1011419

# median B-cell UMI counts
apply(counts(sce)[pcgenes, bcells], MARGIN = 1, median)
# CD79A CD79B  CD19  IGKC IGHG1 IGHG3 IGLC2 IGLC3
#     5     5     0    21     0     0     1     0

sort(counts(sce)[, plasmac], decreasing = TRUE)[1:20]
# IGLC2  IGHG1  IGLC3  IGHG3 MALAT1    B2M  RPL41  RPS18  RPLP1   RPS2  RPL13
#  8029   1982   1480   1409    580    399    359    339    339    269    259
# HLA-B TMSB10  RPL10  RPL21   RPS6  RPS19 RPL18A RPL13A   RPS3
#   239    234    232    220    207    205    199    193    191

# fraction of "IGHG1", "IGHG3", "IGLC2", "IGLC3" counts in plasma cell
sum(counts(sce)[c("IGHG1", "IGHG3", "IGLC2", "IGLC3"), plasmac]) /
    sum(counts(sce)[, plasmac])
#[1] 0.2662923


gdt <- data.table(fts = c(pcgenes, "Total counts"),
    pcell = c(counts(sce)[pcgenes, plasmac], sum(counts(sce)[, plasmac])),
    bcellavg = c(apply(counts(sce)[pcgenes, bcells], MARGIN = 1, mean),
        mean(colSums(counts(sce)[, bcells]))))
gdt[, fts := factor(fts, levels = c(pcgenes, "Total counts"))]
gdt[, pcelllog2 := log10(pcell + 1)]
gdt[, bcellavglog2 := log10(bcellavg + 1)]
gdt[, log2fc := log2(pcell / bcellavg)]
gdt[, ytext := max(.SD), .SDcols = c("pcelllog2", "bcellavglog2"),
    by = seq(nrow(gdt))]
gdtm <- melt(gdt, id.vars = c("fts", "log2fc"),
    measure.vars = c("pcelllog2", "bcellavglog2"))


# figure S3
g <- ggplot() +
    geom_bar(data = gdtm, aes(x = fts, y = value, fill = variable),
        stat = "identity", position = position_dodge()) +
    geom_text(data = gdt,
        aes(x = fts, y = ytext, label = round(log2fc, digits = 2)),
        vjust = -0.3, size = 5) +
    ylab("Counts + 1") +
    xlab(NULL) +
    scale_fill_discrete(labels = c("Plasma cell", "B-cell average")) +
    scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5),
        labels = c(1, 10, 100, "1,000", "10,000", "100,000"),
        limits = c(NA, 5)) +
    .themePublication() +
    theme(legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf("../Figures/FigureS3.pdf", width = 10)
print(g)
dev.off()



