library(gridExtra)
library(ggplot2)

# empty plot 
plt.empt = ggplot() + theme_void()


# load single cell data (3 cell lines in 3 different protocol: 10x. dropseq and celseq2) 
load("/restricted/projectnb/camplab/home/syyang/contamination/data/cellBenchData/CellBench_data_08112019/data/sincell_with_class.RData")
# load single cell data (5 cell lines in 2 different protocols: 10x and celseq2, while celseq2 has 3 replicates)
load("/restricted/projectnb/camplab/home/syyang/contamination/data/cellBenchData/CellBench_data_08112019/data/sincell_with_class_5cl.RData")


# Load decontX result for sinCell_3cl
setwd("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/cellBenchData-ana/analysis_fix/sinCell3/")
dcon.10x = readRDS("dcon.10x.rds")
dcon.celseq2 = readRDS("dcon.celseq2.rds")
dcon.dropseq = readRDS("dcon.dropseq.rds")


# Load decontX result for sinCell_5cl
setwd("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/cellBenchData-ana/analysis_fix/sinCell5/")
dcon.10x_5cl = readRDS("dcon.10x.rds")
dcon.celseq2_p1 = readRDS("dcon.celseq2_p1.rds")
dcon.celseq2_p2 = readRDS("dcon.celseq2_p2.rds")
dcon.celseq2_p3 = readRDS("dcon.celseq2_p3.rds")


addDecontX = function( sce, dcontx) { 
       require(SingleCellExperiment)
       colData(sce)$DecontX_Contamination = dcontx$resList$estConp
       sce@assays$data$decontX = dcontx$resList$estNativeCounts
       return(sce) 
}

extract.dcon = function( sce ) { 
       df = data.frame( DecontX_Contamination = sce@colData$DecontX_Contamination, cell_line_demuxlet = sce@colData$cell_line_demuxlet, doublet = sce@colData$demuxlet_cls, total_cell_line = length(unique(sce@colData$cell_line_demuxlet)) ) 
       return( df ) 
}

rbind.dcon = function( dflist ) { 
    addnames = function( oneDf, oneName) { oneDf$datasetName = rep(oneName, nrow(oneDf[1]) ) 
                                           return( oneDf) } 
    namelist = as.list(names( dflist )) 
    names(namelist) = names( dflist ) 
    nameddf = mapply( addnames ,dflist, namelist, SIMPLIFY=FALSE) 
    rbinddf = do.call(rbind, nameddf) 
    return(rbinddf)
}

createSce = function(sceDcon) { 
    decontX.counts = sceDcon@assays$data$decontX
    colnames(decontX.counts) = rownames(sceDcon@colData) 
    sce = SingleCellExperiment( assays = list( counts = decontX.counts)  , 
			      colData = sceDcon@colData[, c("cell_line_demuxlet", "demuxlet_cls", "DecontX_Contamination")] 
			       ) 
    return(sce) 
}


datasets = list( sc_10x_3cl = sce_sc_10x_qc, sc_celseq2_3cl = sce_sc_CELseq2_qc, sc_dropseq_3cl = sce_sc_Dropseq_qc, 
                 sc_10x_5cl = sce_sc_10x_5cl_qc, 
                 sc_celseq2_5cl_p1 = sc_Celseq2_5cl_p1, sc_celseq2_5cl_p2 = sc_Celseq2_5cl_p2, sc_celseq2_5cl_p3 = sc_Celseq2_5cl_p3)  

dcon = list( sc_10x_3cl = dcon.10x, sc_celseq2_3cl = dcon.celseq2, sc_dropseq_3cl = dcon.dropseq, 
	     sc_10x_5cl = dcon.10x_5cl, 
	     sc_celseq2_5cl_p1 = dcon.celseq2_p1, sc_celseq2_5cl_p2 = dcon.celseq2_p2, sc_celseq2_5cl_p3 = dcon.celseq2_p3)




datasets.wDcon = mapply(addDecontX, datasets, dcon)


datasets.Dcon = lapply(datasets.wDcon, createSce)



dcon.res = lapply( datasets.wDcon, extract.dcon) 


binded.dcon = rbind.dcon( dcon.res) 

combinationnames = paste(binded.dcon$total_cell_line, binded.dcon$datasetName) 
for( i in unique( combinationnames) ) { 
	cat("total_cell_line datasetName : ", i , "\n")
	print(summary( binded.dcon$DecontX_Contamination[ combinationnames == i ] ))
}


# set working directory to save plots
setwd("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/cellBenchData-ana/analysis_fix")

# Violin plot of decontX results separated in singlets and doublets, each protocol, number of cell lines,
pViolin = function(bindedDf) { 
	bindedDf$protocol = gsub( "sc_([^_]+)_.*", "\\1", bindedDf$datasetName)
	bindedDf$part = gsub(".*_" , "", bindedDf$datasetName)
	bindedDf$part[ ! bindedDf$part  %in% c("p1", "p2", "p3") ]   = "p1"
	bindedDf$doublet = plyr::mapvalues( bindedDf$doublet, from = c("DBL", "SNG"), to = c("doublet", "singlet") ) 
	require(ggplot2) 
	bindedDf$datasetName = factor( bindedDf$datasetName, levels = c("sc_10x_3cl", "sc_10x_5cl", "sc_dropseq_3cl", "sc_celseq2_3cl", "sc_celseq2_5cl_p1", "sc_celseq2_5cl_p2", "sc_celseq2_5cl_p3") ) 
	p = ggplot( bindedDf, aes( x = datasetName, y = DecontX_Contamination) ) + 
		ylab("Estimated contamination (%)") + 
		xlab("Protocol") + 
		labs( color = "Protocol" ) + 
	#p = ggplot( bindedDf, aes( x = protocol, y = DecontX_Contamination, fill = part) ) + 
		geom_jitter(  width = 0.3, alpha = 1, size = 0.5, aes( color = doublet ) ) + 
		geom_violin( trim=T, scale = "width",  alpha = 0.5, fill = "grey", aes(color = protocol), draw_quantiles= c(0.5) ) + 
		scale_color_manual( values = c(H1975 = "#E41A1C", H2228 = "#377EB8", HCC827 = "#4DAF4A", H838 = "#984EA3", A549 = "#FFFF33", "10x" = "orangered", celseq2 = "slateblue", dropseq = "olivedrab", sortseq = "purple2", "SNG" = "grey", "DBL"="red4", "singlet" = "grey", "doublet" = "red4") )  + 
		facet_grid( .~ total_cell_line , scale = 'free_x', drop = TRUE, labeller = as_labeller(c("3" = "3 cell line mixtures", "5" = "5 cell line mixtures")) ) +
                     guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), reverse = TRUE))  + 
 	                     theme(panel.background=element_rect(fill="white", color="grey"),
	                          panel.grid = element_line("grey"), 
				  #legend.position="none",
	                          legend.key=element_rect(fill="white", color="white"),
   	                          panel.grid.minor = element_blank(),
				  panel.grid.major = element_blank(),
	                          text=element_text(size=10),
	       	                  axis.text.x = element_text(size=8, angle = 45, hjust = 1),
	       	                  #axis.text.x = element_text(size=8),
		                  axis.text.y=element_text(size=9),
		                  legend.key.size = grid::unit(8, "mm"),
			          legend.text = element_text(size=8), 
				  legend.title = element_text(size = 8), 
			          strip.text.x = element_text(size = 10)  )   + 
	        scale_x_discrete(labels = c("sc_10x_3cl" = "10x", "sc_celseq2_3cl" = "Celseq2", "sc_dropseq_3cl" = "Dropseq", "sc_10x_5cl" = "10x", "sc_celseq2_5cl_p1" = "Celseq2_p1", "sc_celseq2_5cl_p2" = "Celseq2_p2", "sc_celseq2_5cl_p3" = "Celseq2_p3") )   +  
	scale_y_continuous( limits = c(0, 1) ) 
	return( p ) 
}



plt.vio = pViolin( bindedDf = binded.dcon)  + ggtitle("Benchmark single cell datasets") 
#X11(); plt.vio 




### mixRNA analysis 

## load mixRNA data 
load("/restricted/projectnb/camplab/home/syyang/contamination/data/cellBenchData/CellBench_data_08112019/data/mRNAmix_qc.RData")

setwd("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/cellBenchData-ana/analysis_fix/mixRNA/sce2_qc") 
dcon.mixRNA.celseq2 = readRDS("dcon.sce2_qc.rds") 
get.sce2_qc = readRDS("get.sce2_qc.rds")



setwd("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/cellBenchData-ana/analysis_fix/mixRNA/sce8_qc") 
dcon.mixRNA.sortseq = readRDS("dcon.sce8_qc.rds") 
get.sce8_qc = readRDS("get.sce8_qc.rds")

mixRNA_cellFilter = function(sce) { 
    group = paste(sce$H2228_prop,sce$H1975_prop,sce$HCC827_prop) 
    sce = sce[ , group != "0.33 0.33 0.33"  ]
    colData(sce)$cellLine_prop = paste(sce$H2228_prop,sce$H1975_prop,sce$HCC827_prop)
    colData(sce)$dominantCellLine = c("H2228", "H1975", "HCC827")[apply( colData(sce)[, c("H2228_prop", "H1975_prop", "HCC827_prop")], 1, which.max)]
    colData(sce)$dominantCellLine_proportion = apply(colData(sce)[, c("H2228_prop", "H1975_prop", "HCC827_prop")], 1, max)  
    return(sce) 
}


mixRNA.datasets = list( celseq2_mixRNA = sce2_qc, sortseq_mixRNA = sce8_qc)  
lapply( mixRNA.datasets, FUN=function(i)  dim(i@assays$data$counts) )

mixRNA.subsets = lapply( mixRNA.datasets, mixRNA_cellFilter)
lapply( mixRNA.subsets, FUN=function(i)  dim(i@assays$data$counts) )

mixRNA.dcon = list( celseq2_mixRNA = dcon.mixRNA.celseq2, sortseq_mixRNA = dcon.mixRNA.sortseq )

mixRNA.wDcon = mapply(addDecontX, mixRNA.subsets, mixRNA.dcon)

df.dcon = function( sce ) { 
       df = data.frame( DecontX_Contamination = sce@colData$DecontX_Contamination, dominantCellLine= sce@colData$dominantCellLine, dominantCellLine_proportion = sce@colData$dominantCellLine_proportion ) 
       return( df ) 
}

mixRNA.df = lapply( mixRNA.wDcon, df.dcon)

binded.mixRNA = rbind.dcon( mixRNA.df)


# Violin plot of decontX results separated in singlets and doublets, each protocol, number of cell lines,
mixRNA_pViolin = function(bindedDf) { 
	bindedDf$protocol = gsub( "(.*)_(.*)", "\\1", bindedDf$datasetName)
	bindedDf$xAxis = paste0(bindedDf$protocol, bindedDf$dominantCellLine_proportion)
	require(ggplot2) 
	bindedDf$total_cell_line = "mixRNA"  # only used for facet_grid in ggplot
	bindedDf$dominantCellLine_proportion = factor( bindedDf$dominantCellLine_proportion ) 
	p = ggplot( bindedDf, aes( x = xAxis, y = DecontX_Contamination) ) + 
		ylab("Estimated contamination (%)") + 
		xlab("Protocol") + 
		labs( color = "Predominant cell line\npercentage (%)" ) +
	#p = ggplot( bindedDf, aes( x = protocol, y = DecontX_Contamination, fill = part) ) + 
		geom_jitter(  width = 0.3, alpha = 1, size = 0.5, aes( color = dominantCellLine_proportion ) ) + 
		geom_violin( trim=T, scale = "width",  alpha = 0.5, fill = "grey", aes(color = protocol), draw_quantiles= c(0.5) ) + 
		scale_color_manual( values = c(H1975 = "#E41A1C", H2228 = "#377EB8", HCC827 = "#4DAF4A", H838 = "#984EA3", A549 = "#FFFF33", "10x" = "red4", celseq2 = "slateblue", dropseq = "olivedrab3", sortseq = "plum",  "0.68" = "red4", "1" = "grey") )  + 
		facet_grid( .~ total_cell_line , scale = 'free_x', drop = TRUE, labeller = as_labeller(c("mixRNA" = "mixRNA datasets") ) ) +
                     guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))  + 
 	                     theme(panel.background=element_rect(fill="white", color="grey"),
	                          panel.grid = element_line("grey"), 
				  #legend.position="none",
	                          legend.key=element_rect(fill="white", color="white"),
   	                          panel.grid.minor = element_blank(),
				  panel.grid.major = element_blank(),
	                          text=element_text(size=10),
	       	                  axis.text.x = element_text(size=8, angle = 45, hjust = 1),
		                  axis.text.y=element_text(size=9),
		                  legend.key.size = grid::unit(8, "mm"),
			          legend.text = element_text(size=8), 
				  legend.title = element_text(size = 8), 
			          strip.text.x = element_text(size = 10)  )   + 
		scale_y_continuous( limits = c(0, 1) ) + 
	scale_x_discrete( labels = c("celseq20.68" = "Celseq2", "celseq21" = "Celseq2", "sortseq0.68" = "Sortseq", "sortseq1" = "Sortseq") ) 
	return( p ) 
}


p.mixRNA = mixRNA_pViolin(binded.mixRNA ) + ggtitle( "Benchmark mixRNA datasets") 
#X11(); p.mixRNA 

for ( i in unique(binded.mixRNA$dominantCellLine_proportion) ) { 
    for (j in unique(binded.mixRNA$datasetName) ) { 
    cat("dominantCellLine_proportion is : ", i, "in dataset -- ", j, "has contamination summary stats: \n")
    print(summary( binded.mixRNA$DecontX_Contamination[   binded.mixRNA$dominantCellLine_proportion == i & binded.mixRNA$datasetName == j ]  ) )  
    }
}






# PCA plots showing the doublets/singlets by Demuxlet estimation is consistant with decontX contamination% estimation
 ####
 #### code from benchdata paper 
 ####
library(ggplot2)

scran_norm = function(sce){
    require(scran)
    tp = system.time({
		    sce = computeSumFactors(sce)
		    sce = normalize(sce) # goes to `logcounts` by default
		      })
  
    method_name = "scran"
    method_type = "norm"
    if (!is.null(metadata(sce)$running_time)){
        metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
    } else{
        metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
    }
    return(sce)
}









# datasets.filter = lapply(datasets, gene_filter) 
 res = lapply( datasets, scran_norm)
 resPCA = lapply( res, scater::runPCA, exprs_values="logcounts")



 res.dcon = lapply(datasets.Dcon, scran_norm) 
 resPCA.dcon = lapply( res.dcon, scater::runPCA, exprs_values="logcounts") 


plot.pca = function(resPCA) { 
    require(ggplot2) 
    dr = reducedDim(resPCA)
    p = ggplot( data.frame( resPCA@reducedDims$PCA ), aes( x = PC1, y = PC2) ) + 
	    geom_point( size = 0.5,  aes( color = resPCA@colData$cell_line_demuxlet, 
			     shape = resPCA@colData$demuxlet_cls    ) ) + 
            scale_shape_manual(values=c("SNG"=19, "DBL"=0), labels = c("SNG"="singlet", "DBL"="doublet"))+
	    scale_color_manual(values=c(H1975 = "#E41A1C", H2228 = "#377EB8", HCC827 = "#4DAF4A", H838 = "#984EA3", A549 = "#cccb34")) + 
	    xlab(paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=3),"%)")) + ylab(paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=3),"%)") )  + 
	    labs( color = "Cell line", shape = "Demuxlet prediction" ) +  
	    theme(text=element_text(size=8 ),
		  axis.text.x = element_text(size=8),
		  axis.text.y = element_text(size=8), 
                  panel.background=element_blank(),
                  panel.grid = element_blank(), 
		  axis.line = element_line(colour = "black"), 
                  legend.key=element_rect(fill="white", color="white"), 
		  legend.title = element_text(size = 8)
		  ) + 
	    guides(colour = guide_legend(override.aes = list(size = 2)), 
	           shape = guide_legend(override.aes = list(size = 2) )  )
    return(p) 
}

plt.pre = lapply( resPCA, plot.pca)
plt.pos = lapply( resPCA.dcon, plot.pca)

X11();plt.pre$sc_10x_3cl
X11();plt.pre$sc_10x_5cl
X11();plt.pos$sc_10x_5cl

gpca = ggplot_gtable(ggplot_build(plt.pre$sc_10x_5cl))
glegend.pca = gpca$grob[[which(sapply(gpca$grobs, function(x) x$name) == "guide-box")]]

pg.pca = arrangeGrob( plt.pre$sc_10x_3cl + ggtitle("sc_10x (n = 902)\non original counts") + theme(legend.position = "none"),  
		     plt.pos$sc_10x_3cl + ggtitle("sc_10x (n = 902)\non decontaminated counts") + theme(legend.position = "none") , 
		     glegend.pca,
		     plt.pre$sc_10x_5cl + ggtitle("sc_10x_5cl (n = 3,918)\non original counts") + theme(legend.position = "none") , 
		     plt.pos$sc_10x_5cl + ggtitle("sc_10x_5cl (n = 3,918)\non decontaminated counts") + theme(legend.position = "none"), 
		     plt.empt,
                    ncol = 3, widths = c(3,3,3) )  

pg.pca = arrangeGrob( plt.pre$sc_10x_3cl + ggtitle("sc_10x (n = 902)\non original counts") + theme(legend.position = "none"),  
		     plt.pos$sc_10x_3cl + ggtitle("sc_10x (n = 902)\non decontaminated counts") + theme(legend.position = "none") , 
		     plt.pre$sc_10x_5cl + ggtitle("sc_10x_5cl (n = 3,918)\non original counts") + theme(legend.position = "none") , 
		     plt.pos$sc_10x_5cl + ggtitle("sc_10x_5cl (n = 3,918)\non decontaminated counts") + theme(legend.position = "none"), 
				 plt.empt,
                    ncol = 5, widths = c(3,3,3,3, 1) )  


pg.vio = arrangeGrob( p.mixRNA, plt.vio, ncol = 2) 

source("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/cellBenchData-ana/analysis_fix/10xCellRanger3data.R")
plt.3data = plt.10xCellRanger3data + ggtitle("Contamination distributions between different tissues and 10X protocols") + scale_y_continuous( limits = c(0, 1.25), breaks = seq(0,1,0.25)) + geom_text( aes(label = median, y = 1.15) , size = 4)  

pg.10xCellRanger3data = arrangeGrob(plt.3data, glegend.pca, widths = c(4, 1) )


setwd("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/cellBenchData-ana/analysis_fix")
pdf("Fig6_bench_L100.pdf", width = 8.5, height = 8)  
grid.arrange(pg.vio, pg.pca, pg.10xCellRanger3data, ncol = 1, top = "Figure 6", heights = c(4, 3.5,  4)  )
dev.off()




pg.pca2 = arrangeGrob( plt.pre$sc_celseq2_3cl + ggtitle("celseq2_3cl\non original counts") + theme(legend.position = "none"), 
		      plt.pos$sc_celseq2_3cl + ggtitle("celseq2_3cl\non decontaminated counts") + theme(legend.position = "none"), 
                      glegend.pca,
		      plt.pre$sc_dropseq_3cl + ggtitle("dropseq_3cl\non original counts") + theme(legend.position = "none"),
		      plt.pos$sc_dropseq_3cl + ggtitle("dropseq_3cl\non decontaminated counts") + theme(legend.position = "none"),
		      plt.empt,
		      plt.pre$sc_celseq2_5cl_p1 + ggtitle("celseq2_5cl_p1\non original counts") + theme(legend.position = "none"), 
		      plt.pos$sc_celseq2_5cl_p1 + ggtitle("celseq2_5cl_p1\non decontaminated counts") + theme(legend.position = "none"), 
		      plt.empt,
		      plt.pre$sc_celseq2_5cl_p2 + ggtitle("celseq2_5cl_p2\non original counts") + theme(legend.position = "none"), 
		      plt.pos$sc_celseq2_5cl_p2 + ggtitle("celseq2_5cl_p2\non decontaminated counts") + theme(legend.position = "none"), 
		      plt.empt,
		      plt.pre$sc_celseq2_5cl_p3 + ggtitle("celseq2_5cl_p3\non original counts") + theme(legend.position = "none"), 
		      plt.pos$sc_celseq2_5cl_p3 + ggtitle("celseq2_5cl_p3\non decontaminated counts") + theme(legend.position = "none"), 
		      plt.empt, 
		      ncol = 3, widths = c(3,3,4) )

setwd("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/cellBenchData-ana/analysis_fix")
pdf("Sup_Fig_7.pdf", width = 8.5, height = 8.5)
grid.arrange( pg.pca2, ncol = 1, top = "Supplementary Figure 7") 
dev.off()





 # Compare contamination% in singlets VS in doublets (SNG/DBL estimation is from benchdata-paper) 
 summary (dcon.10x$resList$estConp[sce_sc_10x_qc@colData$demuxlet_cls == "DBL"  ] )
 summary (dcon.10x$resList$estConp[sce_sc_10x_qc@colData$demuxlet_cls == "SNG"  ] )



 summary( dcon.celseq2$resList$estConp[sce_sc_CELseq2_qc@colData$demuxlet_cls == "DBL"  ] ) 
 summary( dcon.celseq2$resList$estConp[sce_sc_CELseq2_qc@colData$demuxlet_cls == "SNG"  ] ) 



 summary( dcon.dropseq$resList$estConp[sce_sc_Dropseq_qc@colData$demuxlet_cls == "DBL"  ] )
 summary( dcon.dropseq$resList$estConp[sce_sc_Dropseq_qc@colData$demuxlet_cls == "SNG"  ] )







