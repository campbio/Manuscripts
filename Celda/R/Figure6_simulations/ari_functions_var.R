
library(celda)
library(Seurat)
library(SingleCellExperiment)
library(data.table)
library(mclust)
library(clue)


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


getSeuratGeneModule <- function(ldsdt, rddt, sampleTie = FALSE) {
    for (i in seq(2, ncol(ldsdt))) {
        # rank by absolute value of loadings
        #ldsdt <- ldsdt[order(abs(ldsdt[, i, with = F]), decreasing = TRUE), ]

        ldsdt <- ldsdt[order(ldsdt[, i, with = F], decreasing = TRUE), ]

        if (any(ldsdt[seq(50), ][[i]] < 0)) {
            warning("PC_", i - 1,
                ": Negative loading score found in 50 largest loadings!")
        }

        if (any(ldsdt[seq(nrow(ldsdt) - 49, nrow(ldsdt)), ][[i]] > 0)) {
            warning("PC_", i - 1,
                ": Positive loading score found in 50 smallest loadings!")
        }

        medpos50 <- median(ldsdt[seq(50), ][[i]])
        medneg50 <- median(ldsdt[seq(nrow(ldsdt) - 49, nrow(ldsdt)), ][[i]])

        if (medpos50 + medneg50 > 0) {
            # rank by positive loadings
            ldsdt[, paste0("PC", i - 1, "_rank") := seq(nrow(ldsdt))]
        } else if (medpos50 + medneg50 < 0) {
            # rank by negative loadings
            ldsdt[, paste0("PC", i - 1, "_rank") := rev(seq(nrow(ldsdt)))]
        }
    }
    dt2 <- ldsdt[, c("rn",
        grep("^PC[0-9]", colnames(ldsdt), value = TRUE)), with = FALSE]
    dt2[, seurat_clust := 0]
    isTie <- rep(FALSE, nrow(dt2))
    highrank <- rep(0, nrow(dt2))

    for (i in seq(nrow(dt2))) {
        mv <- min(dt2[i,
            grep("^PC[0-9]", colnames(ldsdt), value = TRUE), with = FALSE])
        wmv <- which(dt2[i, grep("^PC[0-9]", colnames(ldsdt),
            value = TRUE), with = FALSE] == mv)
        if (length(wmv) == 1) {
            dt2[i, seurat_clust := wmv]
        } else {
            isTie[i] <- TRUE
            if (isTRUE(sampleTie)) {
                #print(paste0("Row ", i))
                # dt2[i, clust := sample(wmv, 1)]
                dt2[i, seurat_clust := sample(wmv, 1)]
            }
        }
        highrank[i] <- mv
    }
    dt2[, isTie := isTie]
    dt2[, highrank := highrank]
    dt2 <- dt2[, .(rn, seurat_clust, isTie, highrank)]
    dt2 <- merge(dt2, rddt, by.x = "rn", by.y = "rownames", sort = FALSE)
    return(dt2)
}


getLoadingPC <- function(ldsdt, rddt) {
    dt2 <- ldsdt[, c("rn",
        grep("^PC_[0-9]", colnames(ldsdt), value = TRUE)), with = FALSE]
    dt2[, pcmaxclust := 0]
    dt2[, pcminclust := 0]

    for (i in seq(nrow(dt2))) {
        maxld <- max(dt2[i,
            grep("^PC_[0-9]", colnames(dt2), value = TRUE),
            with = FALSE])
        minld <- min(dt2[i,
            grep("^PC_[0-9]", colnames(dt2), value = TRUE),
            with = FALSE])
        dt2[i, pcmaxclust := which(dt2[i,
            grep("^PC_[0-9]", colnames(dt2), value = TRUE),
            with = FALSE] == maxld)]
        dt2[i, pcminclust := which(dt2[i,
            grep("^PC_[0-9]", colnames(dt2), value = TRUE),
            with = FALSE] == minld)]
    }
    dt2 <- dt2[, .(rn, pcmaxclust, pcminclust)]
    dt2 <- merge(dt2, rddt, by.x = "rn", by.y = "rownames", sort = FALSE)
    return(dt2)
}


getCeldaGModule <- function(sce, L, seed, nchains = 1) {
    # Celda_G
    sceg <- celda_G(sce, L = L, verbose = FALSE, seed = seed,
        nchains = nchains)

    rdg <- rowData(altExp(sceg))
    rdgdt <- as.data.table(rdg)
    rdgdt[, rownames := sub("_", "-", rownames)]
    colnames(rdgdt)[3] <- "celdaGmodule"
    return(rdgdt)
}


getCeldaCGModule <- function(sce, K, L, seed, nchains = 1) {
    # Celda_CG
    scecg <- celda_CG(sce, K = K, L = L, verbose = FALSE, seed = seed,
        nchains = nchains)

    rdcg <- rowData(altExp(scecg))
    rdcgdt <- as.data.table(rdcg)
    rdcgdt[, rownames := sub("_", "-", rownames)]
    colnames(rdcgdt)[3] <- "celdaCGmodule"
    return(rdcgdt)
}


getari <- function(moduleclusters, truelabels) {

    if (0 %in% moduleclusters) {
        stop("0 in 'moduleclusters'!")
    }

    # real cluster X Celda_G cluster
    # tb <- table(truelabels, moduleclusters)
    # tb <- tb[, clue::solve_LSAP(tb, maximum = TRUE)]
    #
    # moduleclusters2 <- plyr::mapvalues(moduleclusters,
    #     as.numeric(rownames(tb)),
    #     as.numeric(colnames(tb)))

    ari <- mclust::adjustedRandIndex(moduleclusters, truelabels)
    return(ari)
}


# G number of genes
# L number of gene modules
simulatedari <- function(model = "celda_CG",
    L,
    S = 1,
    CRange = c(4000, 6000),
    NRange = c(1000, 10000),
    G = 2000,
    K = 20,
    seed,
    beta,
    delta) {

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

    vgres <- getVariableGenes(sc)
    scf <- sc[rd$rownames %in% vgres[[1]], ]

    # L after gene filtering
    fL <- length(unique(rowData(scf)$truelabels))
    altExp(sc, "featureSubset") <- scf

    ldsdt <- seuratldsdt(vgres[[2]], fL)

    rddt <- rowData(altExp(sc))
    rddt <- as.data.table(rddt)

    pcadt <- getLoadingPC(ldsdt, rddt)

    seuratdt <- getSeuratGeneModule(ldsdt, rddt, sampleTie = FALSE)
    seuratdt2 <- getSeuratGeneModule(ldsdt, rddt, sampleTie = TRUE)

    gdt <- getCeldaGModule(sc, fL, seed = seed)
    cgdt <- getCeldaCGModule(sc, K, fL, seed = seed)

    seuratari <- getari(seuratdt[isTie == FALSE, seurat_clust],
        seuratdt[isTie == FALSE, truelabels])
    seuratari2 <- getari(seuratdt2[, seurat_clust], seuratdt2[, truelabels])

    pcaldarimax <- getari(pcadt$pcmaxclust, pcadt$truelabels)
    pcaldarimin <- getari(pcadt$pcminclust, pcadt$truelabels)

    gari <- getari(gdt[, celdaGmodule], gdt[, truelabels])
    cgari <- getari(cgdt[, celdaCGmodule], cgdt[, truelabels])

    return(c(seuratnotie = seuratari,
        seuratwithtie = seuratari2,
        pcaldarimax = pcaldarimax,
        pcaldarimin = pcaldarimin,
        celdagari = gari,
        celdacgari = cgari,
        filteredL = fL))
}

