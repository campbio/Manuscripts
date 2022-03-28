suppressPackageStartupMessages({
	    library(celda)
})

print(packageVersion(pkg="celda"))
print("^ celda version should be: 1.9.1")

apply_celdaC = function(sce, params, k) {
	(seed <- round(1e6*runif(1)))
	# sce: dataset
	# params: other parameters (other than k) - not used for now
	# k: K in celda

	tryCatch(
		{
			#make-sure/convert the data-format for celda function 
			sce_sf = selectFeatures( sce, minCount = params$minCount, minCell = params$minCell, useAssay = "counts", altExpName = "featureSubset" )

			st_ = system.time({
				#clustering using celda and time it
				res_celdac = celda_C(sce_sf, K = k, verbose = FALSE, nchains = 3, seed=seed)
				cluster_ = celdaClusters(res_celdac)
				names(cluster_) = colnames(SingleCellExperiment::altExp(res_celdac, "featureSubset"))
			
			})

			st <- c(user.self = st_[["user.self"]], sys.self = st_[["sys.self"]], 
					user.child = st_[["user.child"]], sys.child = st_[["sys.child"]],
					elapsed = st_[["elapsed"]])
			list(st=st, cluster=cluster_, est_k=NA)

		}, error = function(e) {
			list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
									elapsed = NA), 
					 cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
					 est_k = NA)

		})


}

apply_celdaCG = function(sce, params, k) {
	(seed <- round(1e6*runif(1)))
	# sce: dataset
	# params: other parameters (other than k) - not used for now
	# k: K in celda

	tryCatch(
		{
			#make-sure/convert the data-format for celda function 
			sce_sf = selectFeatures( sce, minCount = params$minCount, minCell = params$minCell, useAssay = "counts", altExpName = "featureSubset" )
			fixedL = min(params$L, dim(SingleCellExperiment::altExp(sce_sf, "featureSubset"))[1])

			st_ = system.time({
				#clustering using celda and time it
				res_celdacg = celda_CG(sce_sf, K = k, L = fixedL, verbose = TRUE, nchains = 3, seed=seed)
				cluster_ = celdaClusters(res_celdacg)
				names(cluster_) = colnames(SingleCellExperiment::altExp(res_celdacg, "featureSubset"))
			
			})

			st <- c(user.self = st_[["user.self"]], sys.self = st_[["sys.self"]], 
					user.child = st_[["user.child"]], sys.child = st_[["sys.child"]],
					elapsed = st_[["elapsed"]])
			list(st=st, cluster=cluster_, est_k=NA)

		}, error = function(e) {
			list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
									elapsed = NA), 
					 cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
					 est_k = NA)

		})


}

