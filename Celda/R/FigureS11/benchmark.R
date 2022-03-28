source("apply_celda.R")  # load apply_celdaC function
source("apply_seurat.R")
source("apply_cidr.R")
source("apply_sc3.R")
source("apply_sc3svm.R")
source("apply_rtsnekmeans.R")
source("apply_pcareduce.R")
source("apply_pcahc.R")
source("apply_pcakmeans.R")
source("apply_tscan.R")
source("apply_flowsom.R")
source("apply_safe.R")


library(ExperimentHub)
eh <- ExperimentHub()

# cache data
data_hub = query(eh, "DuoClustering2018")
print(head(unique(mcols(data_hub)[,'title']), 3))

#datasets
## e.g., dataset_names =c("sce_filteredM3Drop10_Zhengmix8eq", "sce_filteredHVG10_Zhengmix8eq") 
allscenames_index = grep("sce_", unique(mcols(data_hub)[,'title'])) # 48 of them 
allscenames = unique(mcols(data_hub)[,'title'])[allscenames_index]


library(DuoClustering2018)
## Run clustering
record_clustering = function(method, params, sce, scename, n_rep, seed=12345) {
	set.seed(seed)
	L <- lapply(seq_len(n_rep), function(i) {  ## For each run
		cat(paste0("run = ", i, "\n"))
		if (method == "Seurat") {
			tmp <- lapply(params$range_resolutions, function(resolution) {  
				## For each resolution
				cat(paste0("resolution = ", resolution, "\n"))
				## Run clustering
				res <- get(paste0("apply_", method))(sce = sce, params = params, 
																						 resolution = resolution)
				
				## Put output in data frame
				df <- data.frame(dataset = scename, 
												 method = method, 
												 cell = names(res$cluster),
												 run = i,
												 k = length(unique(res$cluster)),
												 resolution = resolution,
												 cluster = res$cluster,
												 stringsAsFactors = FALSE, row.names = NULL)
				tm <- data.frame(dataset = scename, 
												 method = method,
												 run = i, 
												 k = length(unique(res$cluster)),
												 resolution = resolution,
												 user.self = res$st[["user.self"]],
												 sys.self = res$st[["sys.self"]],
												 user.child = res$st[["user.child"]],
												 sys.child = res$st[["sys.child"]],
												 elapsed = res$st[["elapsed"]],
												 stringsAsFactors = FALSE, row.names = NULL)
				kest <- data.frame(dataset = scename, 
													 method = method,
													 run = i, 
													 k = length(unique(res$cluster)),
													 resolution = resolution,
													 est_k = res$est_k,
													 stringsAsFactors = FALSE, row.names = NULL)
				list(clusters = df, timing = tm, kest = kest)
			})  ## End for each resolution
		} else {
      tmp <- lapply(params$range_clusters, function(k) {  ## For each k
        cat(paste0("k = ", k, "\n"))
        ## Run clustering
        res <- get(paste0("apply_", method))(sce = sce, params = params, k = k)
            
        ## Put output in data frame
        df <- data.frame(dataset = scename, 
                         method = method, 
                         cell = names(res$cluster),
                         run = i,
                         k = k,
                         resolution = NA, 
                         cluster = res$cluster,
                         stringsAsFactors = FALSE, row.names = NULL)
        tm <- data.frame(dataset = scename, 
                         method = method,
                         run = i,  
                         k = k,
                         resolution = NA, 
                         user.self = res$st[["user.self"]],
                         sys.self = res$st[["sys.self"]],
                         user.child = res$st[["user.child"]],
                         sys.child = res$st[["sys.child"]],
                         elapsed = res$st[["elapsed"]],
                         stringsAsFactors = FALSE, row.names = NULL)
        kest <- data.frame(dataset = scename, 
                           method = method,
                           run = i,  
                           k = k,
                           resolution = NA, 
                           est_k = res$est_k,
                           stringsAsFactors = FALSE, row.names = NULL)
        list(clusters = df, timing = tm, kest = kest)
      })  ## End for each k
    }
		
		## Summarize across different values of k
		assignments <- do.call(rbind, lapply(tmp, function(w) w$clusters))
		timings <- do.call(rbind, lapply(tmp, function(w) w$timing))
		k_estimates <- do.call(rbind, lapply(tmp, function(w) w$kest))
		list(assignments = assignments, timings = timings, k_estimates = k_estimates)
	})  ## End for each run


## Summarize across different runs
	assignments <- do.call(rbind, lapply(L, function(w) w$assignments))
	timings <- do.call(rbind, lapply(L, function(w) w$timings))
	k_estimates <- do.call(rbind, lapply(L, function(w) w$k_estimates))

## Add true group for each cell
	truth <- data.frame(cell = as.character(rownames(SingleCellExperiment::colData(sce))),
											trueclass = as.character(SingleCellExperiment::colData(sce)$phenoid),
											stringsAsFactors = FALSE)
	assignments$trueclass <- truth$trueclass[match(assignments$cell, truth$cell)]


## Combine results
	res <- list(assignments = assignments, timings = timings,
							k_estimates = k_estimates)


	df <- dplyr::full_join(res$assignments %>%
													 dplyr::select(dataset, method, cell, run, k, 
																				 resolution, cluster, trueclass),
												 res$k_estimates %>%
													 dplyr::select(dataset, method, run, k, 
																				 resolution, est_k)
	) %>% dplyr::full_join(res$timings %>% dplyr::select(dataset, method, run, k,
																											 resolution, elapsed))
	return(df)
}


