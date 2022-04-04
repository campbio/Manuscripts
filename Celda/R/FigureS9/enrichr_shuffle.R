
library(data.table)
library(SingleCellExperiment)
library(enrichR)


parseenriched <- function(enriched, dbs, p = 0.05) {
    for (i in seq(length(dbs))) {
        dti <- enriched[[i]]
        if (nrow(dti) == 0) {
            next
        }

        if (any(dti$Adjusted.P.value < p)) {
            #return (TRUE)
            return (sum(dti$Adjusted.P.value < p))
        }
    }
    return (FALSE)
}


sce <- readRDS("../../Data/sce.rds")
alt <- altExp(sce)

cd <- colData(alt)
rd <- rowData(alt)

dbs <- c("GO_Biological_Process_2021")

mds <- rd$celda_feature_module
l <- length(unique(mds))

set.seed(123)
n <- 100
rnames <- rd$Symbol_TENx

enr100 <- vector("integer", length = n)
enrterm <- vector("list", length = n)

for (i in seq(n)) {
    enum <- 0
    enrtermi <- vector("integer", length = l)
    rni <- rnames[sample(nrow(rd))]

    for (j in seq(l)) {
        message("Module ", j)
        gnamej <- rni[which(rd$celda_feature_module == j)]

        enriched <- enrichr(gnamej, databases = dbs)

        enr <- parseenriched(enriched, dbs = dbs, p = 0.05)
        if (enr) {
            enum <- enum + 1
            enrtermi[j] <- enr
        }
        message("Number of enriched modules: ", enum)
    }

    message("There are ", enum, " enriched gene modules out of ",
        length(unique(rd$celda_feature_module)), " modules")
    enr100[i] <- enum
    enrterm[[i]] <- enrtermi
}

#print(enr100)

dterm <- as.data.table(enrterm)
colnames(dterm) <- paste0("shuffle_", seq(n))

fwrite(data.table(enr100), file = "../../Data/enr_num_bp.csv")
fwrite(dterm, file = "../../Data/shuffle_enr_term_bp.csv")
