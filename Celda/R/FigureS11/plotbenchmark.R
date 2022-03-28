
suppressPackageStartupMessages({
  library(ExperimentHub)
  library(SingleCellExperiment)
  library(DuoClustering2018)
  library(plyr)
	library(ggplot2)
})

#colors <- c(CIDR = "#332288", FlowSOM = "#6699CC", PCAHC = "#88CCEE", 
#            PCAKmeans = "#44AA99", pcaReduce = "#117733",
#            RtsneKmeans = "#999933", Seurat = "#DDCC77", SC3svm = "#661100", 
#            SC3 = "#CC6677", TSCAN = "grey34", ascend = "orange", SAFE = "black",
#            monocle = "red", RaceID2 = "blue")

#method_colors = c(Seurat = "#DDCC77", celdaC="black", celdaCG = "#117733", 
method_colors = c(Seurat = "#DDCC77", celdaCG = "#117733", 
	CIDR = "#332288", FlowSOM = "#6699CC", PCAHC = "#88CCEE", 
  RtsneKmeans = "#999933", SC3svm = "#661100",  pcaReduce = "darkolivegreen2",
  SC3 = "#CC6677", TSCAN = "grey34", SAFE = "black", "PCAKmeans" = "orange")



ari_df <- function(x) {
  stopifnot(methods::is(x, "data.frame"))
  stopifnot(methods::is(x[, 1], "character"))

  x <- dplyr::select(x, -cell)
  columns <- utils::combn(ncol(x), 2)
  ari_nk <- array(NA, ncol(columns))
  for (i in seq_len(length(ari_nk))) {
    ari_nk[i] <- mclust::adjustedRandIndex(x[, columns[1, i]],
                                           x[, columns[2, i]])
  }
  data.frame(ari.stab = ari_nk, run1 = columns[1, ], run2 = columns[2, ],
             stringsAsFactors = FALSE)
}



res_lists = c("res_CeldaCG_L100", "res_CeldaCG_L200", "res_CeldaCG_L30", "res_CeldaCG_L70", "res_CeldaCG_L90",
"res_CeldaCG_L150", "res_CeldaCG_L250", "res_CeldaCG_L60",  "res_CeldaCG_L80")

ref_dir = "res_NoCG"
ref_path = normalizePath(ref_dir)
ref_files = list.files(ref_path)

# load celdaC & Seurat results in a separate dir if provided
ref_df = data.frame()
if (length(ref_files) > 0) {
  for (f in ref_files) {
    if (length(grep(".*\\.rds$", f)) != 0) {
      file_path = file.path(ref_dir, f)
      df = readRDS(file_path)
      ref_df = rbind(ref_df, df)
    }
  }
}


# results from benchmarking each clsutering methods are saved in a .rds file 
# in the form of a data.frame
args = commandArgs(trailingOnly = TRUE)
res_lists = args[1]
for (res_dir in res_lists) { 
    	res_path = normalizePath(res_dir)
    	res_files = list.files(res_path)
    
    # load results
    	combined_df = data.frame()
    	if (length(res_files) > 0) {
    
    		for (f in res_files) {
    			if (length(grep(".*\\.rds$", f)) != 0) {
    				file_path = file.path(res_dir, f)
    				df = readRDS(file_path)
    				combined_df = rbind(combined_df, df) 
    			}
    		}
    
    	}
    
    	combined_df = rbind(combined_df, ref_df)

		## Get ARI
			res_summary <- combined_df %>%
				dplyr::group_by(dataset, method, run, k) %>%
				dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
												 truenclust = length(unique(trueclass)),
												 estnclust = unique(est_k),
												 elapsed = stats::median(elapsed)) %>%
				tidyr::separate(dataset, sep = "_", into = c("sce", "filtering",
																										 "dataset")) %>%
				dplyr::select(-sce) %>% dplyr::ungroup()

		## Save combined_df
			write.csv(combined_df, file.path(res_dir, "combined_df.csv"))
		## Save ARI
			write.csv(res_summary, file.path(res_dir, "median_ari.csv"))

    #plot_dir =args[2]
    	plot_dir = res_dir

    ## Performance plot
    	perf <- plot_performance(combined_df, method_colors = method_colors)
    	perf$median_ari_heatmap_truek <- perf$median_ari_heatmap_truek + geom_text(aes(label = round(medianARI, 2)))
			perf$median_ari_vs_k <- perf$median_ari_vs_k + geom_point(size = 2)
    
    ### Performance: Heatmap of median ARI at true k
    	ggsave(plot = perf$median_ari_heatmap_truek, file=file.path(plot_dir, "median_ari_heatmap_truek.pdf"), width = 20, height = 12)
    	ggsave(plot = perf$median_ari_vs_k, file=file.path(plot_dir, "median_ari_vs_k.pdf"), width = 20, height = 12)  

}

##  Timing plot
#timing <- plot_timing(combined_df, method_colors = method_colors, scaleMethod = "celdaCG")
#names(timing)
#
#ggsave(plot=timing$time_boxplot, file=file.path(plot_dir, "time_boxplot.pdf"), width=14, height=9)
#ggsave(plot=timing$time_vs_k, file=file.path(plot_dir, "time_vs_k.pdf"), width=14, height=9)
#ggsave(plot=timing$time_normalized_by_ref, file=file.path(plot_dir, "time_normalized_by_ref.pdf"), width=14, height=9)
#
#
## Calculate ARI at each K for each datasets
#res_summary <- combined_df %>% 
#    dplyr::group_by(dataset, method, run, k) %>% 
#    dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
#                     truenclust = length(unique(trueclass)),
#                     estnclust = unique(est_k),
#                     elapsed = stats::median(elapsed)) %>% 
#    tidyr::separate(dataset, sep = "_", into = c("sce", "filtering",
#                                                 "dataset")) %>% 
#    dplyr::select(-sce) %>% dplyr::ungroup()
#
#saveRDS(object = res_summary, file=file.path(plot_dir, "calculated_ARIs.rds")


### Stability plot
#stab <- plot_stability(combined_df, method_colors = method_colors)
#
#ggsave(plot = stab$stability_vs_k, file=file.path(plot_dir, "stability_vs_k.pdf"))
#ggsave(plot = stab$stability_heatmap_truek, file=file.path(plot_dir, "stability_heatmap_truek.pdf"))
#
#
### Entropy plot 
#entr <- plot_entropy(combined_df, method_colors = method_colors)
#
#ggsave(plot = entr$entropy_vs_k, file=file.path(plot_dir, "entropy_vs_k.pdf"))
#ggsave(plot = entr$normentropy, file=file.path(plot_dir, "normentropy.pdf"))
#ggsave(plot = entr$deltaentropy_truek, file=file.path(plot_dir, "deltaentropy_truek.pdf"))
#

rm(ref_df)
