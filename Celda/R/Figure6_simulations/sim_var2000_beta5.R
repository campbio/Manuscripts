
library(BiocParallel)
source("./ari_functions_var.R")

cores <- 28

seed <- 12345
betatests <- c(1, 5, 10, 20, 40, 80, 160)
deltatests <- c(1, 5, 10, 20, 40, 80, 160)
G <- 33000

rep = 10
llistunique <- seq(10, 200, 10)
llist <- rep(llistunique, times = rep)
#llist <- c(10, 20)
klist <- rep(20, length(llist))
# seuratnotieari <- vector("numeric", length = length(llist))
# seuratwithtieari <- vector("numeric", length = length(llist))
# celdagari <- vector("numeric", length = length(llist))
# celdacgari <- vector("numeric", length = length(llist))
seeds <- seq(seed, seed - 1 + length(llist))
beta <- betatests[2]
delta <- deltatests[2]


simall <- function(i, klist, llist, seeds, beta, delta, G) {
    # for windows SnowParam
    # source("./20201026_ari_functions.R")
    message(Sys.time(), " Simulating L = ", llist[i])
    res <- simulatedari(L = llist[i],
        seed = seeds[i],
        beta = beta,
        delta = delta,
        G = G)
    resdf <- data.frame(klist = klist[i],
        llist = llist[i],
        seuratnotieari = res["seuratnotie"],
        seuratwithtieari = res["seuratwithtie"],
        pcaldarimax = res["pcaldarimax"],
        pcaldarimin = res["pcaldarimin"],
        celdagari = res["celdagari"],
        celdacgari = res["celdacgari"],
        filteredL = res["filteredL"],
        seed = seeds[i],
        beta = beta,
        delta = delta)
    print(resdf)
}

resL <- BiocParallel::bplapply(X = seq(length(llist)),
    FUN = simall,
    BPPARAM = BiocParallel::MulticoreParam(workers = cores, log = TRUE,
        progressbar = TRUE),
    # BPPARAM = BiocParallel::SnowParam(workers = cores, log = TRUE,
    # progressbar = TRUE),
    klist = klist,
    llist = llist,
    seeds = seeds,
    beta = beta,
    delta = delta,
    G = G)

resdt <- as.data.table(plyr::rbind.fill(resL))

fwrite(resdt, file = paste0("simulated_2000var_beta",
    beta, "_delta",
    delta, ".csv"))

sessionInfo()
