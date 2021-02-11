library(celda)
library(ggplot2)

pbmc4kdec <- readRDS("../Data/pbmc4kdec.rds")

# selecting the number of L and K
rsm <- recursiveSplitModule(pbmc4kdec, useAssay = "decontXcounts", maxL = 200)
saveRDS(rsm, file = "../Data/rsm.rds")
rsm <- readRDS("../Data/rsm.rds")

rsmp <- plotGridSearchPerplexity(rsm, sep = 10)
rsmp <- rsmp + scale_x_discrete(breaks = seq(0, 200, 20),
    limits = factor(seq(0, 200, 1))) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
g1 <- plotGridSearchPerplexityDiff(rsm, sep = 10, n = 30)
g1 <- g1 + scale_x_discrete(breaks = seq(0, 200, 20),
    limits = factor(seq(0, 200, 1))) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
    ylab("RPC of L")

print(rsmp)
print(g1)
L <- 80

pbmc4kfL80 <- subsetCeldaList(rsm, list(L = L))

rsc <- recursiveSplitCell(pbmc4kfL80,
    useAssay = "decontXcounts",
    initialK = 3,
    maxK = 30,
    yInit = celdaModules(pbmc4kfL80))

saveRDS(rsc, file = "../Data/rsc_l80.rds")
rsc <- readRDS("../Data/rsc_l80.rds")

rscp <- plotGridSearchPerplexity(rsc, sep = 1) +
    scale_x_discrete(breaks = seq(0, 30, 5),
        limits = factor(seq(0, 30, 1))) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
g2 <- plotGridSearchPerplexityDiff(rsc, n = 1)
g2 <- g2 + scale_x_discrete(breaks = seq(0, 30, 5),
    limits = factor(seq(0, 30, 1))) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.key.size = unit(16, "mm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
    ylab("RPC of K")

print(rscp)
print(g2)
K <- 20

pdf("../Figures/FigureS1.pdf", width = 10)
print(rsmp)
print(g1)
print(rscp)
print(g2)
dev.off()

pbmc4kfK20L80 <- subsetCeldaList(rsc,
    params = list(K = K, L = L))
pbmc4kfK20L80 <- celdaUmap(pbmc4kfK20L80,
    useAssay = "decontXcounts",
    nNeighbors = 10,
    minDist = 0.5)

# Manual reordering of transcriptional modules and cell clusters
sce <- recodeClusterY(pbmc4kfK20L80,
    c(c(23, 74, 22), # Plasma
        c(21, seq(24, 26), seq(28, 33)), # B
        c(27, 8, 40, 67), # DC
        c(7, seq(34, 38), 66, 39, 64), # pDC
        c(41, 76, 79, 80), # CD34+
        c(78, 77, 12, 13, seq(43, 48), 68), # NK
        c(1), # Mk
        c(seq(2, 6), seq(10, 11), 16), #CD14+
        c(9, 65, 14, 15, seq(17, 19), 42, 20), # FCGR3A+
        c(71, 75), # proliferating T
        c(70, seq(55, 56), 49, 57, 59), # T cytotoxic
        c(50, seq(53, 54)), # GZMK+ T
        60, # middle T
        c(seq(51, 52), 58, seq(61, 62)), # T
        c(63, 69, seq(72, 73))), # everything
    c(seq(1, 2), # plasma
        seq(3, 12), #B
        seq(13, 16), #DC
        seq(17, 25), #pDC
        seq(26, 29), # CD34+
        seq(30, 40), # NK
        41, # Mk
        seq(42, 49), # CD14+
        seq(50, 58), # FCGR3A+
        59, # proliferating T
        seq(60, 64), # T cytotoxic
        seq(65, 67), # GZMK+ T
        68, # middle T
        seq(69, 73), # T
        seq(74, 80))) # everything

sce <- recodeClusterZ(sce,
    c(1, 7, 5, 6, 3, 4, 2, 8, 15, 17, 16, 19, 20, 18, 11, 13, 12, 9, 10, 14),
    seq(20))
saveRDS(sce, file = "../Data/sce.rds")
