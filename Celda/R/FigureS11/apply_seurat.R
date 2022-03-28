## Apply Seurat

suppressPackageStartupMessages({
  library(Seurat)
})

apply_Seurat <- function(sce, params, resolution) {
  (seed <- round(1e6*runif(1)))
  tryCatch({
    dat <- counts(sce)
		# Todo: Any simple filtering on the count matrix is here
		#
    st <- system.time({
      data <- CreateSeuratObject(counts = dat, min.cells = params$min.cells,
                                 min.features = params$min.genes, project = "scRNAseq", 
                                 display.progress = TRUE) 
			# log-normalize count cross features after scaling within each cell w/ scale.factor
      data <- NormalizeData(object = data, normalization.method = "LogNormalize", 
                            scale.factor = 1e4, display.progress = FALSE)
			# Find varibale features/genes
			# only if none is variable by their used method, all genes will be used then
      data <- ScaleData(object = data, display.progress = FALSE)
			# PCA on cellxgene matrix, weighted on gene/feature's variance 
			# features with 0 variance is dropped 
			# TODO: might change npcs = min(dim(sce)) ??
      data <- RunPCA(object = data, assay = "RNA", features = rownames(data@assays$RNA@counts), do.print = FALSE, 
                     npcs = params$dims.use, seed.use = seed)
			#  k.param = 20 is as default
			#  `feature` is set as null, and hence the variable genes/features are used
			data <- FindNeighbors(object = data,k.param = 20, features = NULL, compute.SNN = TRUE, reduction="pca")
			# resolution (>1.0, or <1.0) to find (larger, or smaller) #clusters
      data <- FindClusters(object = data, resolution = resolution, verbose = TRUE, 
                           random.seed = seed)
      cluster <- data$seurat_clusters
    })
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), 
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}


###############################################################################
# Argument "params" for each method:
		## copied from get_res.R

#' params_list = list("celdaC" = list(range_clusters = c(-----)),
#'                   "Seurat" = list(range_resolutions =c(-----) , min.cells = ----- , min.genes = ----, dims.use = ----) 
#'                  )  

# Seurat[dims.use] -->>  RunPCA(npcs): number oc PC components. If we use CeldaCG, ideally set this the same as L ? 


#<<<

###############################################################################
# Notes on Parameters used in each method 

#Seurat
#' @param min.cells Include features detected in at least this many cells. Will
#' subset the counts matrix as well. To reintroduce excluded features, create a
#' new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are
#' detected.

#<<<
