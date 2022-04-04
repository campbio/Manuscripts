
library(SingleCellExperiment)
library(celda)
library(benchmarkme)

# print(get_ram())
# print(get_cpu())
# print(get_platform_info())

sce <- readRDS("../../Data/pbmc68k_umap.rds")

n <- 20
K <- 20
L <- 80

reslistcg <- vector("list", length = n)

seeds <- c(11, 12, 13, 14, 16, 17, 18, 19, 21, 24,
    seq(25, 29), c(31, 33, 35, 36, 37))

for (i in seq(n)) {
    reslistcg[[i]] <- system.time(
        res <- celda_CG(sce,
            useAssay = "decontXcounts",
            altExpName = "VST2000",
            K = K,
            L = L,
            seed = seeds[i],
            nchains = 1,
            zInitialize = "random",
            yInitialize = "random"),
        gcFirst = TRUE)
}

dfcg <- as.data.frame(do.call(rbind, reslistcg))

write.csv(dfcg, file = "../../Data/pbmc68k_cg_2000g.csv")
