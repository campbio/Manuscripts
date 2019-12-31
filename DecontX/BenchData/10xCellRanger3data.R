library(ggplot2) 
library(SingleCellExperiment)
library(celda)



## Load DecontX results on datasets processed by Cell Ranger 3.0 from 10X Genomics 
dcon_PBMCV2.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/PBMC1KV2/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_PBMCV3.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/PBMC1KV3/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_BrainV2.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/BrainV2/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_BrainV3.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/BrainV3/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_HeartV2.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/HeartV2/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_HeartV3.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/HeartV3/analysis_fix/filterG_sce_dcon_eps1L100.rds")



# Plot contamination level for these 4 datasets 
extractEstCon = function( sce_dcon) {
	df = data.frame("estConp" = sce_dcon$DecontX_Contamination, "cluster" = sce_dcon$DecontX_QuickCluster, "data" = deparse(substitute(sce_dcon)))
  df$data = gsub(".*_", "", df$data)
	df$data = gsub("\\..*", "", df$data)
	df$protocols = gsub("(.*)(V.)", "\\2", df$data)
	df$name = gsub("(.*)(V.)", "\\1_1K", df$data)
	df$median = paste0(round( 100 * median(df$estConp) , 2), "%") 
	df$median[2: nrow(df)] = "" 
	return(df)
}

#estCon_conca = rbind(extractEstCon(dcon_PBMC4k), extractEstCon(dcon_Brain1KCR2),   extractEstCon(dcon_PBMC1KV2), extractEstCon(dcon_PBMC1KV3), extractEstCon(dcon_BrainV2), extractEstCon(dcon_BrainV3) )
estCon_conca = rbind(extractEstCon(dcon_PBMCV2.filterG), extractEstCon(dcon_PBMCV3.filterG), extractEstCon(dcon_BrainV2.filterG),  extractEstCon(dcon_BrainV3.filterG), extractEstCon(dcon_HeartV2.filterG), extractEstCon(dcon_HeartV3.filterG) )

plot.dataViolin = function( estConp, color = c("red4", "grey", "yellow", "green", "black")) {
	       df = estConp
       p = ggplot( df, aes( x = data, y = estConp) ) +
				       labs( color = "10X protocol") +
							       ylab("Estimated Contamination") +
										 xlab("Dataset") +
										 #geom_jitter( width = 0.2, alpha = 1, size = 0.3, aes(color = as.factor(cluster))  ) +
										 geom_jitter( width = 0.2, alpha = 1, size = 0.3, aes(color = protocols)  ) +
										 geom_violin( trim=T, scale = "width", fill = "grey", alpha = 0.5 ) +
										 scale_color_manual( values = color ) +
										 guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)) ) + 
										theme(panel.background=element_rect(fill="white", color="grey"),
													 panel.grid = element_line("grey"),
														#legend.position="none",
														legend.key=element_rect(fill="white", color="white"),
													panel.grid.minor = element_blank(),
														 panel.grid.major = element_blank(),
														 text=element_text(size=10),
														axis.text.x = element_text(size=8, angle = 45, hjust = 1),
														axis.text.y=element_text(size=8),
													 legend.key.size = grid::unit(8, "mm"),
													 legend.text = element_text(size=8),
													legend.title = element_text(size = 8),
													strip.text.x = element_text(size = 10)  ) + 
                    scale_x_discrete(labels = c("BrainV2" = "Brain_1K", "BrainV3" = "Brain_1K", "HeartV2" = "Heart_1K", "HeartV3" = "Heart_1K", "PBMCV2" = "PBMC_1K", "PBMCV3" = "PBMC_1K"))
return( p)
}


plt.10xCellRanger3data = plot.dataViolin(estCon_conca )  
#ggsave(plot=p, file = "violin_estCon3types.pdf", width = 7, height = 4)


# umap of marker genes 
DimReduce = function( sce, featureName, genelist) {
umap = sce@reducedDims$QC_UMAP 
counts = assay(sce, "counts")
	counts = as.matrix(counts) 
	geneName = featureName[featureName$V1 %in% rownames(counts), "V2"]
	rownames(counts) = geneName
  p =celda::plotDimReduceFeature(counts = counts , dim1 = umap[, 1], dim2 = umap[, 2], features = genelist)
	return(p)
}
#DimReduce( dcon_BrainV2, featureName =BrainV2genes, genelist = BrainMarker[, 2])  


# Function to plot tsne with predicted values (estimated-contamination in DecontX, or predited doublet-score in Scrublet)
plot.est = function( dim1, dim2, scaleValue, size = 1, varLabel, colorLow = "grey80", colorHigh = "blue4", colorMid = NULL) {
	  m = scaleValue
  df = data.frame( x = dim1, y = dim2, m = m )
	  p = ggplot( df, aes( x = x, y = y) ) + geom_point( stat ="identity", size = size, aes(color = m) ) +
			    xlab("Dimension_1") + ylab("Dimension_2") +
				 scale_colour_gradient2(low = colorLow,high = colorHigh,mid = colorMid, guide = "colorbar",
														midpoint = ( (max(m) + min(m))/ 2), name = varLabel ) +
				 theme_bw() +
					theme(strip.background = element_blank(),
							 panel.grid.major = element_blank(),
							panel.grid.minor = element_blank(),
							panel.spacing = unit(0, "lines"),
						 panel.background = element_blank(),
						 axis.line = element_line(colour = "black"),
						legend.title = element_text(size = 8)
							)

							        return(p )
}

# umap of contamination levels
DimReduceEst = function( sce ) {
	umap = sce@reducedDims$QC_UMAP
  p = plot.est(dim1 =  umap[, 1], dim2 = umap[, 2], scaleValue = sce$DecontX_Contamination, varLabel = "estConp")
	return(p)
}

## umap plot of clusters
getUmap = function(sce, color = color19) { 
   umap = sce@reducedDims$QC_UMAP
   qz =   sce$DecontX_QuickCluster
   df = data.frame( dim1 = umap[, 1], dim2 = umap[, 2], qz = qz ) 
	 p = celda::plotDimReduceCluster(dim1 = umap[, 1], dim2 = umap[, 2], cluster = qz, labelClusters = TRUE)
	 			return(p)
}

dataset.list = c("dcon_PBMCV2.filterG", "dcon_PBMCV3.filterG", "dcon_BrainV2.filterG", "dcon_BrainV3.filterG", 
								 "dcon_HeartV2.filterG", "dcon_HeartV3.filterG") 

dataset = "dcon_PBMCV2.filterG"


analysis_path = "/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/pbmc-ana/4k/gene2col2/experimentalEnvironmentEffect"

for (dataset in dataset.list) {
dataname = gsub("(.*_)(.*)(\\..*)", "\\2", dataset) 
p = getUmap( get(dataset) ) + ggtitle( paste0( dataname, " auto clustering w/ eps 1") ) 
#ggsave(plot = p, file = paste0(analysis_path, "/_auto_cluster_umap_", dataname, "_filterG.pdf") ) 
p = DimReduceEst( get(dataset) ) + ggtitle( paste0( dataname, " contamination level") ) 
#ggsave(plot = p, file = paste0(analysis_path, "/_estConp_", dataname, "_filterG.pdf") )
}

for (dataset in dataset.list) { 
	cat("summary statistics of estimated contamination proportion in ", dataset,  " data: \n") ; print( summary( get(dataset)$DecontX_Contamination ) )
}

for (dataset in dataset.list) {
	cat("total ", ncol(get(dataset)), " cells in ", dataset, "\n")
}
### plot marker genes on Brain cells  
