NFEATURES = 2000
CURRENT_DIR=normalizePath(".")
SUB_DIR = paste("sce_nfeature", NFEATURES, sep="_")
SCE_DIR = file.path(CURRENT_DIR, SUB_DIR)
source("benchmark.R")

## Set number of times to run clustering for each k
n_rep <- 5


args =  commandArgs(trailingOnly = TRUE)
saveres_dir = normalizePath(args[1])
override = args[2] #always


print(saveres_dir)
# For all available datasets -- see benchmark.R script
datasets =c("sce_filteredM3Drop10_Zhengmix8eq", "sce_filteredHVG10_Zhengmix8eq", "sce_filteredExpr10_Zhengmix8eq",
						"sce_filteredM3Drop10_Zhengmix4eq", "sce_filteredHVG10_Zhengmix4eq", "sce_filteredExpr10_Zhengmix4eq",
						"sce_filteredM3Drop10_Zhengmix4uneq", "sce_filteredHVG10_Zhengmix4uneq", "sce_filteredExpr10_Zhengmix4uneq", 
            "sce_filteredM3Drop10_Koh", "sce_filteredHVG10_Koh", "sce_filteredExpr10_Koh", 
            "sce_filteredM3Drop10_KohTCC", "sce_filteredHVG10_KohTCC", "sce_filteredExpr10_KohTCC", 
						"sce_filteredM3Drop10_Trapnell", "sce_filteredHVG10_Trapnell", "sce_filteredExpr10_Trapnell", 
						"sce_filteredM3Drop10_TrapnellTCC", "sce_filteredHVG10_TrapnellTCC", "sce_filteredExpr10_TrapnellTCC", 
						"sce_filteredM3Drop10_Kumar", "sce_filteredHVG10_Kumar", "sce_filteredExpr10_Kumar",
						"sce_filteredM3Drop10_KumarTCC", "sce_filteredHVG10_KumarTCC", "sce_filteredExpr10_KumarTCC",
						"sce_filteredM3Drop10_SimKumar8hard", "sce_filteredHVG10_SimKumar8hard", "sce_filteredExpr10_SimKumar8hard", 
						"sce_filteredM3Drop10_SimKumar4hard", "sce_filteredHVG10_SimKumar4hard", "sce_filteredExpr10_SimKumar4hard", 
						"sce_filteredM3Drop10_SimKumar4easy", "sce_filteredHVG10_SimKumar4easy", "sce_filteredExpr10_SimKumar4easy")
						
    
# Don't use simulated datasets
# Simulation process is not how celda works and you don't know that is the real machnism how cells' expression work
# simulation process is briefly mentioned in the Duoclustering2018 paper
datasets =c("sce_filteredM3Drop10_Zhengmix8eq", "sce_filteredHVG10_Zhengmix8eq", "sce_filteredExpr10_Zhengmix8eq",
            "sce_filteredM3Drop10_Zhengmix4eq", "sce_filteredHVG10_Zhengmix4eq", "sce_filteredExpr10_Zhengmix4eq",
            "sce_filteredM3Drop10_Zhengmix4uneq", "sce_filteredHVG10_Zhengmix4uneq", "sce_filteredExpr10_Zhengmix4uneq", 
            "sce_filteredM3Drop10_Koh", "sce_filteredHVG10_Koh", "sce_filteredExpr10_Koh", 
            "sce_filteredM3Drop10_KohTCC", "sce_filteredHVG10_KohTCC", "sce_filteredExpr10_KohTCC", 
            "sce_filteredM3Drop10_Trapnell", "sce_filteredHVG10_Trapnell", "sce_filteredExpr10_Trapnell", 
            "sce_filteredM3Drop10_TrapnellTCC", "sce_filteredHVG10_TrapnellTCC", "sce_filteredExpr10_TrapnellTCC", 
            "sce_filteredM3Drop10_Kumar", "sce_filteredHVG10_Kumar", "sce_filteredExpr10_Kumar",
            "sce_filteredM3Drop10_KumarTCC", "sce_filteredHVG10_KumarTCC", "sce_filteredExpr10_KumarTCC")

    
#datasets=c("sce_filteredM3Drop10_KohTCC")


# Params choices are adpted from the paper:
# A systematic performance evaluation of clustering methods for single-cell RNA-seq data 
# https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/Rscripts/parameter_settings/generate_parameter_settings.R
truencluster = list()
params_list  = list()
for (scename in datasets) {
  #sce = get(scename)()
	# Use most variable genes subsetted and saved as .rds
	sce = readRDS(file.path(SCE_DIR, paste0(scename, ".rds")))
  nc = length(unique(as.character(SingleCellExperiment::colData(sce)$phenoid)))
  truencluster[[scename]] = nc
  min_c = max(2, nc-3)
  max_c = nc+3
  range_clusters = c(min_c:max_c)
  params_list[[scename]] <- list("celdaC" = list("range_clusters" = c(min_c:max_c), "minCount"=1, "minCell"=1), 
                                 "celdaCG" = list("range_clusters" = c(min_c:max_c), "minCount"=1, "minCell"=1, "L"=60), 
                                 "Seurat" = list("range_resolutions" = seq(0.3, 1.5, by = 0.1) , min.cells = 0, min.genes =0, dims.use = 30) 
                                 )   
}

# Special datasets param-settting from pp:  A systematic performance evaluation of clustering methods for single-cell RNA-seq data. F1000Research 7:1141 (2018).
# https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/Rscripts/parameter_settings/generate_parameter_settings.R
for (scename in datasets) { 
  if (length(grep("Zhengmix8eq", scename)) > 0) {
    if (length(grep("filteredExpr10", scename)) > 0) {
      params_list[[scename]][["Seurat"]][["range_resolutions"]]  = c(0.01, 0.3, 0.35, 0.4, 0.5, 1.3, 1.4, 1.6, 1.7, 1.8)
    } else if (length(grep("filteredHVG10", scename)) > 0) {
      params_list[[scename]][["Seurat"]][["range_resolutions"]]  = c(0.1, 0.3, 1, 1.1, 1.3, 1.4, 1.6, 2)
    } else if (length(grep("filteredM3Drop10", scename)) > 0) {
      params_list[[scename]][["Seurat"]][["range_resolutions"]]  = c(0.1, 0.7, 1, 1.3, 1.6, 2.5, 3)
    }   
  } else if (length(grep("Zhengmix4eq|Zhengmix4uneq", scename)) > 0) {
      params_list[[scename]][["Seurat"]][["range_resolutions"]]  =  c(0.05, 0.1, 0.2, seq(0.3, 1.5, by = 0.1))
			if (length(grep("(filteredM3Drop10).*(Zhengmix4eq)", scename)) > 0) {
				params_list[[scename]][["Seurat"]][["range_resolutions"]]  =  c(0.05, seq(0.1, 0.2, by=0.02), seq(0.3, 1.5, by = 0.1))
			}
			if (length(grep("(filteredHVG10).*(Zhengmix4uneq)", scename)) > 0) {
				params_list[[scename]][["Seurat"]][["range_resolutions"]]  =  c(seq(0.01, 0.05, by=0.005), 0.1, 0.2, seq(0.3, 1.5, by = 0.1))
			}
			if (length(grep("(filteredM3Drop10).*(Zhengmix4uneq)", scename)) > 0) {
				params_list[[scename]][["Seurat"]][["range_resolutions"]]  =  c(seq(0.12, 0.14, by=0.002), seq(0.3, 1.5, by = 0.1))
			}
  } else if (length(grep("Koh|KohTCC", scename)) > 0 ) { 
      if (length(grep("(filteredM3Drop10).*(Koh$)", scename)) > 0) {
        params_list[[scename]][["Seurat"]][["range_resolutions"]]  = seq(0.3, 3.6, by = 0.1)
      } else if (length(grep("(filteredExpr10).*(KohTCC$)", scename)) > 0) {
        params_list[[scename]][["Seurat"]][["range_resolutions"]]  =  sort(c(1.73, 1.76, seq(0.3, 2.1, by = 0.1)))
      } else {
        params_list[[scename]][["Seurat"]][["range_resolutions"]]  = seq(0.3, 2.1, by = 0.1)
      }
  } else if (length(grep("(filteredM3Drop10).*(SimKumar8hard)", scename)) > 0) {
    params_list[[scename]][["Seurat"]][["range_resolutions"]]  = sort(c(1.13, 1.15, 1.18, seq(0.3, 1.5, by = 0.1)))
  }
}

print(params_list)


# specify clustering method to use
# and their parameters to be specified
method_list = c( "celdaCG", "Seurat", "CIDR", "RtsneKmeans", "PCAHC", "SC3svm", "SC3", "TSCAN", "FlowSOM", "pcaReduce", "PCAKmeans", "SAFE")


# record 'method' clustering performance(speed+accuracy)  
for (scename in datasets) {
	#sce = get(scename)()
	# Use most variable genes subsetted and saved as .rds
	sce = readRDS(file.path(SCE_DIR, paste0(scename, ".rds")))

	for (method in method_list) {

		cat(scename, "==============================================================================================using     ", method, "\n")
		cat("---its' dimension:", dim(sce), "\n")
		file_name = paste(method, scename, ".rds", sep="_")
		file_path = file.path(saveres_dir, file_name)
		df = record_clustering(method, params_list[[scename]][[method]], sce, scename, n_rep, seed=123) 
		saveRDS(object=df, file=file_path)

	}

}

