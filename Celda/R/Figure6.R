# Figure 6: ARI plots and UMAP examples

library(celda)
library(ggplot2)
library(data.table)
library(grid)
library(gridExtra)
library(SingleCellExperiment)
library(Seurat)


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
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                legend.key = ggplot2::element_rect(color = NA),
                legend.position = "right",
                legend.direction = "vertical",
                # legend.key.size = ggplot2::unit(0.2, "cm"),
                # legend.margin = ggplot2::margin(0),
                # legend.title = ggplot2::element_text(face = "bold"),
                plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
                strip.background = ggplot2::element_rect(
                    color = "#f0f0f0", fill = "#f0f0f0"),
                strip.text = ggplot2::element_text(face = "bold")
            ))
}

# See scripts in Figure6_simulations for simulation details
b1 <- fread("../Data/Figure6_simulations/simulated_2000var_beta1_delta1.csv")
b5 <- fread("../Data/Figure6_simulations/simulated_2000var_beta5_delta5.csv")
b10 <- fread("../Data/Figure6_simulations/simulated_2000var_beta10_delta10.csv")
b20 <- fread("../Data/Figure6_simulations/simulated_2000var_beta20_delta20.csv")
b30 <- fread("../Data/Figure6_simulations/simulated_2000var_beta30_delta30.csv")
b40 <- fread("../Data/Figure6_simulations/simulated_2000var_beta40_delta40.csv")

dt <- rbind(b1, b5, b10, b20, b30, b40)

dt2 <- dt[, .(llist, filteredL, seuratnotieari,
    celdacgari, beta, delta)]

dtm <- melt(dt2, id = c("llist", "beta", "delta"),
    measure = c("seuratnotieari", "celdacgari"))

# median ARIs for 6 combinations
dtm[, .(beta, variable, value)][, lapply(.SD, median), by = .(beta, variable)]

# median
dtm1 <- dtm[beta == 1, ]
dtmed1 <- dtm1[, .(median = median(value),
    q75 = quantile(value)[4],
    q25 = quantile(value)[2]), by = .(llist, variable)]

dtm5 <- dtm[beta == 5, ]
dtmed5 <- dtm5[, .(median = median(value),
    q75 = quantile(value)[4],
    q25 = quantile(value)[2]), by = .(llist, variable)]

dtm10 <- dtm[beta == 10, ]
dtmed10 <- dtm10[, .(median = median(value),
    q75 = quantile(value)[4],
    q25 = quantile(value)[2]), by = .(llist, variable)]

dtm20 <- dtm[beta == 20, ]
dtmed20 <- dtm20[, .(median = median(value),
    q75 = quantile(value)[4],
    q25 = quantile(value)[2]), by = .(llist, variable)]

dtm30 <- dtm[beta == 30, ]
dtmed30 <- dtm30[, .(median = median(value),
    q75 = quantile(value)[4],
    q25 = quantile(value)[2]), by = .(llist, variable)]

dtm40 <- dtm[beta == 40, ]
dtmed40 <- dtm40[, .(median = median(value),
    q75 = quantile(value)[4],
    q25 = quantile(value)[2]), by = .(llist, variable)]


g8.0 <- ggplot() +
    geom_line(data = dtmed1, aes(x = llist, y = median, color = variable)) +
    geom_pointrange(data = dtmed1, aes(x = llist, y = median,
        ymin = q75,
        ymax = q25,
        color = variable)) +
    xlab("Number of gene modules before \nselection of variable genes") +
    ylab("ARI") +
    scale_colour_hue(name = "Methods",
        labels = c("PCA_loading_rank",
            "Celda_CG")) +
    ggtitle(expression(~beta~" = 1, "~delta~" = 1")) +
    .themePublication() +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(-0.01, 1)) +
    xlim(0, 200) +
    theme(legend.position = "none")

g8.1 <- ggplot() +
    geom_line(data = dtmed5, aes(x = llist, y = median, color = variable)) +
    geom_pointrange(data = dtmed5, aes(x = llist, y = median,
        ymin = q75,
        ymax = q25,
        color = variable)) +
    xlab("Number of gene modules before \nselection of variable genes") +
    ylab("ARI") +
    scale_colour_hue(name = "Methods",
        labels = c("PCA_loading_rank",
            "Celda_CG")) +
    ggtitle(expression(~beta~" = 5, "~delta~" = 5")) +
    .themePublication() +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(-0.01, 1)) +
    xlim(0, 200) +
    theme(legend.position = "none")

g8.2 <- ggplot() +
    geom_line(data = dtmed10, aes(x = llist, y = median, color = variable)) +
    geom_pointrange(data = dtmed10, aes(x = llist, y = median,
        ymin = q75,
        ymax = q25,
        color = variable)) +
    xlab("Number of gene modules before \nselection of variable genes") +
    ylab("ARI") +
    scale_colour_hue(name = "Methods",
        labels = c("PCA_loading_rank",
            "Celda_CG")) +
    ggtitle(expression(~beta~" = 10, "~delta~" = 10")) +
    .themePublication() +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(-0.01, 1)) +
    xlim(0, 200) +
    theme(legend.position = "none")


g8.3 <- ggplot() +
    geom_line(data = dtmed20, aes(x = llist, y = median, color = variable)) +
    geom_pointrange(data = dtmed20, aes(x = llist, y = median,
        ymin = q75,
        ymax = q25,
        color = variable)) +
    xlab("Number of gene modules before \nselection of variable genes") +
    ylab("ARI") +
    scale_colour_hue(name = "Methods",
        labels = c("PCA_loading_rank",
            "Celda_CG")) +
    ggtitle(expression(~beta~" = 20, "~delta~" = 20")) +
    .themePublication() +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(-0.01, 1)) +
    xlim(0, 200) +
    theme(legend.position = "none")

g8.4 <- ggplot() +
    geom_line(data = dtmed30, aes(x = llist, y = median, color = variable)) +
    geom_pointrange(data = dtmed30, aes(x = llist, y = median,
        ymin = q75,
        ymax = q25,
        color = variable)) +
    xlab("Number of gene modules before \nselection of variable genes") +
    ylab("ARI") +
    scale_colour_hue(name = "Methods",
        labels = c("PCA_loading_rank",
            "Celda_CG")) +
    ggtitle(expression(~beta~" = 30, "~delta~" = 30")) +
    .themePublication() +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(-0.01, 1)) +
    xlim(0, 200) +
    theme(legend.position = "none")

g8.5 <- ggplot() +
    geom_line(data = dtmed40, aes(x = llist, y = median, color = variable)) +
    geom_pointrange(data = dtmed40, aes(x = llist, y = median,
        ymin = q75,
        ymax = q25,
        color = variable)) +
    xlab("Number of gene modules before \nselection of variable genes") +
    ylab("ARI") +
    scale_colour_hue(name = "Methods",
        labels = c("PCA_loading_rank",
            "Celda_CG")) +
    ggtitle(expression(~beta~" = 40, "~delta~" = 40")) +
    .themePublication() +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(-0.01, 1)) +
    xlim(0, 200) +
    theme(legend.position = "none")

ariplots <- arrangeGrob(grobs = list(g8.0, g8.1, g8.2, g8.3, g8.4, g8.5),
    nrow = 1, newpage = FALSE)


# UMAPs of 2000 most vatiable genes
seed <- 12354

getVariableGenes <- function(simsce, nfeatures = 2000) {
    simseurat <- suppressWarnings(CreateSeuratObject(counts = counts(simsce),
        project = "simulated", min.cells = 0, min.features = 0))

    simseurat <- NormalizeData(simseurat, verbose = FALSE)
    simseurat <- FindVariableFeatures(simseurat, nfeatures = nfeatures,
        verbose = FALSE)
    vf <- VariableFeatures(simseurat)
    return(list(vf, simseurat))
}


seuratldsdt <- function(simseurat, L) {
    simseurat <- ScaleData(simseurat, verbose = FALSE)
    simseurat <- RunPCA(simseurat, npcs = L, verbose = FALSE)

    lds <- simseurat@reductions$pca@feature.loadings
    ldsdt <- as.data.table(lds, keep.rownames = TRUE)
    return(ldsdt)
}


# UMAP using true cell cluster labels
getUMAP2 <- function(seed, beta, delta,
    model = "celda_CG",
    S = 1,
    CRange = c(4000, 6000),
    NRange = c(1000, 10000),
    G = 33000,
    K = 20,
    L = 100,
    nNeighbors = 10,
    minDist = 0.5) {

    sc <- simulateCells(model = model,
        S = S,
        CRange = CRange,
        NRange = NRange,
        G = G,
        K = K,
        L = L,
        seed = seed,
        beta = beta,
        delta = delta)

    rd <- rowData(sc)
    colnames(rd)[2] <- "truelabels"
    rd$rownames <- sub("_", "-", rd$rownames)
    rownames(sc) <- rd$rownames
    rowData(sc) <- rd

    vsseuratobj <- getVariableGenes(sc, nfeatures = 2000)
    ldsdt <- seuratldsdt(vsseuratobj[[2]], L)

    scf <- sc[rd$rownames %in% ldsdt$rn, ]

    altExp(sc, "featureSubset") <- scf

    # UMAP of 2000 variable genes
    cts <- counts(scf)

    z <- colData(sc)$celda_cell_cluster
    ctsbycl <- celda:::.colSumByGroup(cts, group = z, K = K)

    ctsnm <- normalizeCounts(t(ctsbycl))
    ctsnm <- sqrt(ctsnm)

    umap <- uwot::umap(t(ctsnm),
        n_neighbors = nNeighbors,
        #pca = doPCA,
        min_dist = minDist)
    dt <- data.table(umap)
    colnames(dt) <- c("UMAP_1", "UMAP_2")
    tlb <- rowData(altExp(sc))
    tlb <- factor(tlb$truelabels)
    levels(tlb) <- seq_along(unique(tlb))
    dt[, truelbl := tlb]
    return(dt)
}

dt1 <- getUMAP2(seed = seed,
    beta = 1,
    delta = 1)

dt5 <- getUMAP2(seed = seed,
    beta = 5,
    delta = 5)

dt10 <- getUMAP2(seed = seed,
    beta = 10,
    delta = 10)

dt20 <- getUMAP2(seed = seed,
    beta = 20,
    delta = 20)

dt30 <- getUMAP2(seed = seed,
    beta = 30,
    delta = 30)

dt40 <- getUMAP2(seed = seed,
    beta = 40,
    delta = 40)


plotgeneUMAP <- function(dt, size = 0.7) {
    clusterColors <- distinctColors(nlevels(as.factor(dt$truelbl)))
    g <- ggplot(dt) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = truelbl),
            size = size) +
        xlab(NULL) +
        ylab(NULL) +
        theme(panel.border = element_rect(color = "grey90", fill = NA),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            strip.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.spacing = unit(0, "lines"),
            panel.background = ggplot2::element_blank()) +
        ggplot2::scale_color_manual(values = clusterColors)
    return(g)
}

gb1 <- plotgeneUMAP(dt1)
gb5 <- plotgeneUMAP(dt5)
gb10 <- plotgeneUMAP(dt10)
gb20 <- plotgeneUMAP(dt20)
gb30 <- plotgeneUMAP(dt30)
gb40 <- plotgeneUMAP(dt40)

umapplots <- arrangeGrob(grobs = list(gb1, gb5, gb10, gb20, gb30, gb40),
    nrow = 1, newpage = FALSE)

figure6 <- arrangeGrob(grobs = list(ariplots, umapplots),
    nrow = 2, newpage = FALSE)

pdf("../Figures/Figure6.pdf", width = 22)
grid.draw(figure6)
dev.off()
