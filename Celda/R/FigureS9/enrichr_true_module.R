
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

l <- length(unique(rd$celda_feature_module))

enum <- 0
terms <- vector("integer", length = l)

set.seed(123)
for (i in seq(l)) {
    message("Module ", i)
    gname <- rd[which(rd$celda_feature_module == i), "Symbol_TENx"]
    enriched <- enrichr(gname, databases = dbs)

    enr <- parseenriched(enriched, dbs = dbs, p = 0.05)
    if (enr) {
        enum <- enum + 1
        terms[i] <- enr
    }
    message("Number of enriched modules: ", enum)
}

message("There are ", enum, " enriched gene modules out of ",
    length(unique(rd$celda_feature_module)), " modules")

fwrite(data.table(terms), file = "../../Data/true_module_term_bp.csv")
