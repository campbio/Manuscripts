library(data.table)

file <- "../../Data/pbmc4k_dec_cpm.tab.chars.blocks"
con <- file(file, open = "r")
f <- readLines(con)

geneclstlist <- list()
cellclstlist <- list()

for (i in seq(length(f))) {
    line <- f[i]
    if (startsWith(f[i], " Genes ")) {
        geneclstlist <- c(geneclstlist, line)
    } else if (startsWith(f[i], " Conds ")) {
        cellclstlist <- c(cellclstlist, line)
    }
}
close(con)

gcl <- lapply(geneclstlist,
    function(x) {
        x2 <- gsub("^.*\\: ", "", x)
        x3 <- gsub("\ $", "", x2)
        x4 <- tstrsplit(x2, "\\ ")
        x5 <- unlist(x4)
        x6 <- sapply(x5, function(x) {
            spl <- tstrsplit(x, "\\_")
            if (length(spl) > 2) {
                stop("Gene name length longer than 2!")
            }
            return(spl[[1]])
        })
        return(x6)
    })

ccl <- lapply(cellclstlist,
    function(x) {
        x2 <- gsub("^.*\\: ", "", x)
        x3 <- gsub("\ $", "", x2)
        x4 <- unlist(tstrsplit(x2, "\\ "))
        return(x4)
    })

length(unique(unlist(gcl)))
#[1] 148

length(unique(unlist(ccl)))
#[1] 333

saveRDS(gcl, "../../Data/qubic2_gene_clusters.rds")
saveRDS(ccl, "../../Data/qubic2_cell_clusters.rds")
