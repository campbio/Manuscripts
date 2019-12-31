RNGkind(sample.kind = "Rounding")
library(devtools)
load_all("~/gitProjects/celda_irisapo", recompile=T)   # on fix branch, celda version 1.3.1
library(ggplot2) 
library(gtable)
library(gridExtra) # arrange plot
library(SingleCellExperiment)
source("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/pbmc-ana/4k/gene2col2/functions.R")
BrainMarker = read.table("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/pbmc-ana/4k/gene2col2/experimentalEnvironmentEffect/BrainInfo/geneList.txt", sep=" ", stringsAsFactors=FALSE, header=T)

# empty plot 
plt.empt = ggplot() + theme_void()

######## New results
## Load DecontX results on datasets processed by Cell Ranger 3.0 from 10X Genomics 
dcon_PBMCV2.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/PBMC1KV2/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_PBMCV3.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/PBMC1KV3/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_BrainV2.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/BrainV2/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_BrainV3.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/BrainV3/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_HeartV2.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/HeartV2/analysis_fix/filterG_sce_dcon_eps1L100.rds")
dcon_HeartV3.filterG = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/HeartV3/analysis_fix/filterG_sce_dcon_eps1L100.rds")


dataset.list = c("dcon_PBMCV2.filterG", "dcon_PBMCV3.filterG", "dcon_BrainV2.filterG", "dcon_BrainV3.filterG", 
								 "dcon_HeartV2.filterG", "dcon_HeartV3.filterG") 


for (dataset in dataset.list) { 
	cat("summary statistics of estimated contamination proportion in ", dataset,  " data: \n") ; print( summary( get(dataset)$DecontX_Contamination ) )
}
for (dataset in dataset.list) {
	cat("total ", ncol(get(dataset)), " cells in ", dataset, "\n")
}

# functions to extract information for plotting 
extractEstCon = function( sce_dcon) {
	df = data.frame("dim_1" = sce_dcon@reducedDims$QC_UMAP[, 1], "dim_2" = sce_dcon@reducedDims$QC_UMAP[, 2], 
									"estConp" = sce_dcon$DecontX_Contamination, "cluster" = sce_dcon$DecontX_QuickCluster, 
									"data" = deparse(substitute(sce_dcon)))
  df$data = gsub(".*_", "", df$data)
	return(df)
}

estCon_conca = rbind(extractEstCon(dcon_PBMCV2.filterG), extractEstCon(dcon_PBMCV3.filterG), extractEstCon(dcon_BrainV2.filterG),  extractEstCon(dcon_BrainV3.filterG), extractEstCon(dcon_HeartV2.filterG), extractEstCon(dcon_HeartV3.filterG) )

p = plot.dataViolin(estCon_conca, color = color19 ) + ggtitle("Contamination distributions across tissues and 10X protocols")
#ggsave(plot=p, file = "violin_estCon3types.pdf", width = 7, height = 4)


dataset = "dcon_PBMCV2.filterG"

plist.dcon = list()
for (dataset in dataset.list) {
  dataname = gsub("(.*_)(.*)(\\..*)", "\\2", dataset) 
  df = extractEstCon( get(dataset) ) 
  plist.dcon[[ dataname ]] = plot.est( dim1 = df$dim_1, dim2 = df$dim_2, scaleValue = df$estConp, size = 0.3, varLabel = "Contamination (%)") + ggtitle( dataname )  
}

plist.cluster = list() 
for (dataset in dataset.list) {
  dataname = gsub("(.*_)(.*)(\\..*)", "\\2", dataset)
  df = extractEstCon( get(dataset) )
	plist.cluster[[ dataname ]] = plot.cluster(dim1 = df$dim_1, dim2 = df$dim_2, cluster = df$cluster, size = 0.3, varLabel = "Cluster") + ggtitle( dataname)  
}

#X11(); plist.dcon$PBMCV2
#X11(); plist.cluster$PBMCV2

gdecontx = ggplot_gtable(ggplot_build(plist.dcon$HeartV2 + theme( legend.position = "bottom") ))
glegend.decontx = gdecontx$grob[[which(sapply(gdecontx$grobs, function(x) x$name) == "guide-box")]]

gcluster = ggplot_gtable(ggplot_build(plist.cluster$HeartV2 + theme( legend.position = "bottom") ))
glegend.cluster = gcluster$grob[[which(sapply(gdecontx$grobs, function(x) x$name) == "guide-box")]]

pg.clusterV2= arrangeGrob(plist.cluster$PBMCV2 + theme( legend.position = "none" ), 
												 plist.cluster$BrainV2 + theme( legend.position = "none" ), 
												 plist.cluster$HeartV2 + theme( legend.position = "none" ), 
												 plt.empt,
												 ncol = 1)
pg.clusterV3= arrangeGrob(
												 plist.cluster$PBMCV3 + theme( legend.position = "none" ), 
												 plist.cluster$BrainV3 + theme( legend.position = "none" ), 
												 plist.cluster$HeartV3 + theme( legend.position = "none" ), 
												 plt.empt,
												 ncol = 1 )  
pg.dconV2 = arrangeGrob( plist.dcon$PBMCV2 + theme( legend.position = "none" ) , 
												 plist.dcon$BrainV2 + theme( legend.position = "none" ), 
												 plist.dcon$HeartV2 + theme( legend.position = "none" ), 
												 glegend.decontx, 
												 ncol = 1) 
pg.dconV3 = arrangeGrob(
												 plist.dcon$PBMCV3 + theme( legend.position = "none" ), 
												 plist.dcon$BrainV3 + theme( legend.position = "none" ), 
												 plist.dcon$HeartV3 + theme( legend.position = "none" ), 
												 glegend.cluster,
												 ncol = 1 )  

pdf("Sup_Fig_8.pdf", width = 8.5, height = 8)
grid.arrange( pg.clusterV2, pg.dconV2 , pg.clusterV3, pg.dconV3,   ncol = 4, widths = c(4, 4, 4,  4),  top = "Supplementary Figure 8")
dev.off()

