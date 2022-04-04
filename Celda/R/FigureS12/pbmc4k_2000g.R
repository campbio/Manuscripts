
library(SingleCellExperiment)
library(celda)
library(benchmarkme)

# print(get_ram())
# print(get_cpu())
# print(get_platform_info())

sce <- readRDS("../../Data/sce.rds")
altExp(sce, "VST2000") <- altExp(sce)

n <- 20

reslistcg <- vector("list", length = n)

for (i in seq(n)) {
    reslistcg[[i]] <- system.time(
        res <- celda_CG(sce,
            useAssay = "decontXcounts",
            altExpName = "VST2000",
            K = 20,
            L = 80,
            seed = i,
            nchains = 1,
            zInitialize = "random",
            yInitialize = "random"),
        gcFirst = TRUE)
}

dfcg <- as.data.frame(do.call(rbind, reslistcg))

write.csv(dfcg, file = "../../Data/pbmc4k_cg_2000g.csv")
