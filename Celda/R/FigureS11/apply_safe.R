## Apply SAFE clustering
ProgramDir="SAFEclustering/gpmetis_and_shmetis_for_Linux"

#source("SAFE_2.1_Linux/individual_clustering_modified.R")
#source("SAFE_2.1_Linux/SAFE_modified.R")
suppressPackageStartupMessages({
	library("SAFEclustering")
	#library(Seurat)
	#library(cidr)
	#library(rtsne)
	#library(ADPclust)
})

apply_SAFE <- function(sce, params, k) {
	(seed <- round(1e6*runif(1)))
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      indc <- individual_clustering(inputTags = dat, 
																		mt_filter = TRUE, 
																		SC3 = TRUE, gene_filter = FALSE, CIDR = TRUE, nPC.cidr = NULL, 
																		Seurat = TRUE, nGene_filter = FALSE, nPC.seurat = NULL, resolution = 0.7, tSNE = TRUE, dimensions = 3, 
																		perplexity = 30, SEED = seed)
      safe <- SAFE(cluster_results = indc, k_min = k, k_max = k, program.dir = ProgramDir, 
                   MCLA = TRUE, HGPA = FALSE, CSPA = FALSE, 
                   cspc_cell_max = NULL)
      
      cluster <- safe$optimal_clustering
      names(cluster) <- colnames(dat)
    })
    ## Determine number of clusters automatically
    safe <- SAFE(cluster_results = indc, k_min = 2, k_max = NULL, program.dir = ProgramDir,
                 MCLA = TRUE, HGPA = FALSE, CSPA = FALSE, 
                 cspc_cell_max = NULL)
    est_k <- safe$optimal_k
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = est_k)
  }, error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), 
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
