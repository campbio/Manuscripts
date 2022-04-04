
library(SingleCellExperiment)
library(celda)
library(Seurat)
library(benchmarkme)

# print(get_ram())
# print(get_cpu())
# print(get_platform_info())

sce <- readRDS("../../Data/sce.rds")

pbmc4kseurat <- CreateSeuratObject(counts = decontXcounts(sce),
    project = "pbmc4kdec", min.cells = 3, min.features = 200)

nfeat <- 1000

pbmc4kseurat <- NormalizeData(pbmc4kseurat)
pbmc4kseurat <- FindVariableFeatures(pbmc4kseurat, nfeatures = nfeat)
vf1 <- VariableFeatures(pbmc4kseurat)

pbmc4kf <- sce[vf1, ]

altExp(sce, paste0("VST", nfeat)) <- pbmc4kf

n <- 20

reslistcg <- vector("list", length = n)

for (i in seq(n)) {
    reslistcg[[i]] <- system.time(
        res <- celda_CG(sce,
            useAssay = "decontXcounts",
            altExpName = paste0("VST", nfeat),
            K = 20,
            L = 80,
            seed = i,
            nchains = 1,
            zInitialize = "random",
            yInitialize = "random"),
        gcFirst = TRUE)
}

dfcg <- as.data.frame(do.call(rbind, reslistcg))

# print("dfcg")
# print(paste0("Number of genes: ", nfeat))
# print(dfcg)

write.csv(dfcg, file = paste0("../../Data/pbmc4k_cg_", nfeat, "g.csv"))
