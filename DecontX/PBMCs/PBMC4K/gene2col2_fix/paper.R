#brary(gtable) load colSumByGroup function
library(devtools)
load_all(path="~/gitProjects/celda_irisapo", recompile=T)
library(ggplot2) 
library(gtable) 
library(reshape2)
library(ggrepel)
library(gridExtra) # arrange plot
library(RcppCNPy) # load result from Scrublet 
source("functions.R")

# Empty plot 
plt.empt = ggplot() + theme_void() 


# load PBMC-4K data 
source("loaddata.R")
K = 19 ; L = 150

cluster.pre = readRDS("rsm.L150.K.celda_yusukeCode.rds")
cluster.pre = subsetCeldaList( cluster.pre, list(L=L, K = K ) ) 
cluster.pos = readRDS( paste0("post.cluster.L", L, "K",K, ".rds") )

cluster.pre.reZ = recodeClusterZ( cluster.pre, from = 1:19, to = c(1,10,2,3,     5,7,6,4,   9,8,13,17,   14,15,16,11,   12,18,19 ))
cluster.pos.reZ = reorderZ( ori.celda = cluster.pre.reZ, decon.celda = cluster.pos ) 

dcon = readRDS(paste0("dcon.L", L, "K", K,".rds"))
pbmc4K_dcon = round(dcon$resList$estNativeCounts)
cat("summary statistics of estimated contamination proportion in PBMC4K data: \n") ; print( summary( dcon$resList$estConp ) )  


# load Scrublet predicted doublets 
doublet_scores = npyLoad("doublet_scores.npy")
predicted_doublets = npyLoad("predicted_doublets.npy", type = "integer") 
predicted_doublets[ predicted_doublets == 0 ]  = "singlet"
predicted_doublets[ predicted_doublets == 1 ]  = "doublet"


# load combined immune cells data 
source("loaddataCI.R")
#source("functionsCI.R")

cluster.sort = readRDS("../../../Immune-positivecontrol-ana/rsm.L76rsm.K.celda.rds" )
cluster.sort  = subsetCeldaList( cluster.sort, list( L = 76, K= 21 ) )

dcon.immune = readRDS("../../../Immune-positivecontrol-ana/gene2col2_filtering/dcon.K21L76.rds")
immuneCB.c = convert.sortLabel( ori.label = mat.s.label, new.label = c( "Tcells", "Natural_Killer_Cells", "Bcells", "Monocytes", "CD34_cells")  )
immuneCB.c[ immuneCB.c == "Natural_Killer_Cells"  ] = "NKcells"  # rename 
immuneCB.re= adjust.label( celda.cluster = cluster.sort@clusters$z, sorted.label = immuneCB.c )

# Show celda corrected clustering result compared with flowcytomerty sorted cell type label
table( mat.s.label, cluster.sort@clusters$z ) 

for (cell in c("Monocytes", "Tcells", "Bcells", "NKcells")) {
	counts.cell = mat.s[, immuneCB.re== cell]
	cat("in Sorted PBMCs cell type ", cell, "has ", ncol(counts.cell), "cells \n", "       summary statistics of total UMIs is \n")
	print(summary(colSums(counts.cell)))
}

#cat("summary statistics of estimated contamination proportion in sorted immune cells data: \n") ; print( summary( dcon.immune$resList$estConp ) )  

g = c("IL7R","CD3D", "CD3E", "GNLY", "IRF7", "LYZ", "S100A8" ,"S100A9" , "FCGR3A", "MS4A1", "FCER1A", "IGHG1","IGHG2",  "PPBP", "CD34", "CD14") 
markerGenes = grep(paste(g, collapse="|"), rownames(pbmc4k_select) , value=TRUE)


# Tsne of combined immune cells data 
#RNGkind(sample.kind="Rounding")  # to be consistent with R/3.5.* random seed generating function 
#set.seed(12345) 
#tsne.immune = celdaTsne( counts = mat.s, celdaMod = cluster.sort, seed = 12345 ) 
#plt.immune = plotDimReduceCluster(dim1=tsne.immune[,1], dim2=tsne.immune[,2], cluster=cluster.sort@clusters$z, size = 0.3, labelClusters = TRUE ) + theme( legend.position = "none" ) + ggtitle("Sorted PBMCs")  + theme( plot.margin=unit(c(0, 1,0, 0) , "lines")) 


# >>>  Tsne plots 
RNGkind(sample.kind="Rounding")  # to be consistent with R/3.5.* random seed generating function 
set.seed(12345) 
tsne.pre = celdaTsne(counts=pbmc4k_select, celdaMod = cluster.pre.reZ, seed = 12345)
#saveRDS(tsne.pre, "tsne.pre.rds")
plt.t1 = plotDimReduceCluster(dim1=tsne.pre[,1], dim2=tsne.pre[,2], cluster=cluster.pre.reZ@clusters$z, size = 0.3, labelClusters = TRUE ) + theme( legend.position = "none" ) + ggtitle("Original 4K PBMC")  + theme( plot.margin=unit(c(0, 1,0, 0) , "lines")) 
tsne.pre = readRDS("tsne.pre.rds")


# according to marker-tsne plot
T.cluster = c(5, 6, 7, 9, 8, 11) # cluster-16 is activated T cells (CD3D/CD8A/HMGB2/HMGN2/MKI67+) 
NK.cluster = c(10) 
B.cluster = c(2, 3)     # MS4A1 B cells  
Monocytes.cluster = c(13, 14, 15 )  # CD14+ or CD16+ monocytes 
#MonocytesCD14.cluster = c(13, 14) # CD14+ monocytes
#MonocytesCD16.cluster = c(15) # CD16+ monocytes (FCGR3A+)
Dendtritic.cluster = c(18, 19) # FCER1A+ Dendritic cells
pDendritic.cluster = c(4) # IRF7+ Plasmacytoid-Dendritic cells 
multiplets.cluster = c(17)  
Megakaryocytes.cluster = c(16) # PPBP+
Plasma.cluster = c(1) # IGHG1/2/3 
CD34Cell.cluster = c(12) # CD34


RNGkind(sample.kind="Rounding")  # to be consistent with R/3.5.* random seed generating function 
set.seed(12345)
tsne.pos = celdaTsne(pbmc4K_dcon, celdaMod = cluster.pos.reZ, seed = 12345) 
#saveRDS(tsne.pos, "tsne.pos.rds") 
plt.t2 = plotDimReduceCluster( dim1= tsne.pos[,1], dim2=tsne.pos[,2], cluster= cluster.pos.reZ@clusters$z, size=0.3, labelClusters = TRUE) + theme( legend.position = "none") + ggtitle("Decontaminated 4K PBMC") + theme(  plot.margin=unit(c(0, 1,0, 0) , "lines")) 
tsne.pos = readRDS("tsne.pos.rds")


# >>>  Silhouette width plots 
silhouette.pre = s.width( counts=pbmc4k_select, z=cluster.pre.reZ@clusters$z) 
silhouette.pos = s.width( counts=pbmc4K_dcon, z=cluster.pos.reZ@clusters$z) 
#saveRDS(silhouette.pos, "silhouette.pos.rds")
silhouette.pre  =  readRDS("silhouette.pre.rds")
silhouette.pos  = readRDS("silhouette.pos.rds") 

cat("The mean silhouette width on original PBMC 4K data after cpm normalization is : \n ", mean( silhouette.pre[, 3] ), "\n" )  
cat("The mean silhouette width on decontaminated PBMC 4K data after cpm normalization is : \n ", mean( silhouette.pos[, 3] ), "\n" )  

for ( cluster in 1:19) {
  cat("The summary statistics of silhouette width in cluster", cluster, "\n",  
			"     in original 4K PBMC counts is \n",  
			summary(silhouette.pre[ silhouette.pre[, "cluster"]== cluster, "sil_width" ]), "\n",  
			"     in decontaminated 4K PBMC counts is \n", 
			summary(silhouette.pos[ silhouette.pos[, "cluster"]== cluster, "sil_width" ]), "\n")
}

color.tsne = exc.color( plt.ggplot = plt.t2, z = cluster.pos.reZ@clusters$z) 
plt.sil = silhouette.plot(silhouette.ori = silhouette.pre, silhouette.dec = silhouette.pos, K=K, color = color.tsne) + xlab("data") + theme( plot.margin=unit(c(2, 0,2, 1) , "lines")) + guides(col = guide_legend(nrow = 10)) 
plt.silDiff = silhouetteDiff.plot(silhouette.ori = silhouette.pre, silhouette.dec = silhouette.pos, K=K, color = color.tsne, exK = c(1)  ) + xlab("4K PBMC") + theme( plot.margin=unit(c(0, 1,0, 1) , "lines")) 

plt.method = plot.methodViol(estConp = dcon$resList$estConp, doublet = predicted_doublets) + xlab("Scrublet predicted type") + ylab("DecontX estimated\ncontamination proportion") + ggtitle("Contamination in\npredicted doublets") + scale_y_continuous(limits = c(0, 1))  
	# + theme(plot.margin=unit(c(0, 4,0, 1), "lines"),   legend.margin = margin(0,0,0,0,"mm") , legend.position = c(1.25, 0.8)  ) 


g <- ggplot_gtable(ggplot_build(plt.sil))
glegend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]


pg.tsne = arrangeGrob( plt.empt,  plt.t1, plt.t2, ncol = 3, widths = c(0.5, 4, 4) ) 
pg.sil = arrangeGrob(  plt.empt, plt.silDiff + theme(legend.position = "none"), glegend, plt.method, ncol = 4 , widths = c(0.2, 2.8, 1, 3.5) ) 

pdf("Fig5.pdf", width = 8.5, height = 7)
grid.arrange( pg.tsne, pg.sil, ncol = 1 , top="Figure 5") 
dev.off()



# >>> Scatter plots 
TB.Genes = c( "IL32","CD3E","CD3D", "MS4A1","CD79A","CD79B" ) 
	     c("HLA-DRB1","HLA-DPB1","HLA-DRA","HLA-DPA1")  # no "HLA-DQA2" in pbmc4K data , "CD74" is the farest gene in B cell

z.TB = label.AB( z=cluster.pre.reZ@clusters$z, A.cluster = B.cluster, B.cluster = T.cluster, A.name="Bcells", B.name="Tcells" ) 
df.BT.gene = cluster.gene(counts=pbmc4k_select, z.AB=z.TB)

z.TB.pos = label.AB( z=cluster.pos.reZ@clusters$z, A.cluster = B.cluster, B.cluster = T.cluster , A.name="Bcells", B.name="Tcells" ) 
df.BT.decon = cluster.gene(counts=pbmc4K_dcon, z.AB=z.TB.pos)

z.TB.immune = label.AB( z = immuneCB.re, A.name = "Bcells", B.name = "Tcells" ) 
df.BT.immune = cluster.gene(counts = mat.s, z.AB = z.TB.immune ) 


p.TB.ori = plot.gene.exp( df.AB.gene = df.BT.gene, x="log2.Tcells", y="log2.Bcells", Genes = TB.Genes, colors = color30 ) + ggtitle("Original 4K PBMC (profiled\nin the same channel)") + xlab("") + ylab("")  
p.TB.decon = plot.gene.exp( df.BT.decon, x="log2.Tcells", y="log2.Bcells", Genes = TB.Genes, colors = color30 ) +ggtitle("Decontaminated 4K PBMC (profiled\nin the same channel)") + xlab("") + ylab("")
p.TB.immune = plot.gene.exp( df.BT.immune, x="log2.Tcells", y="log2.Bcells", Genes = TB.Genes, colors = color30 ) + ggtitle("Sorted PBMCs (profiled\nin different channels)") + xlab("") + ylab("")


df.TB.all = rbind.geneExp( geneExp.ori = df.BT.gene, geneExp.pos = df.BT.decon, geneExp.immune = df.BT.immune )  
p.TB.all =  plot.gene.exp( df.TB.all, x="log2.Tcells", y="log2.Bcells", Genes = TB.Genes, colors = color30 ) + facet_grid(. ~ datatype_f ) + ylab("Average expression\nof B cell (log2)") + xlab("Average expression of T cell (log2)") 


# >>> Stackbar plots 
genes =  c("MS4A1", "CD79A", "CD79B", "CD3D", "CD3E", "LYZ", "S100A8", "S100A9", "GNLY") 
#gene.group = list( "Bcell.marker" = c("MS4A1", "CD79B", "CD79A"), "Tcell.marker"=c("CD3D", "CD3E"), "NKcell.marker"=c("GNLY"), "Monocytes.marker"=c("LYZ", "S100A8", "S100A9") )
gene.group = list( "Bcells.marker" = c("MS4A1", "CD79B", "CD79A"), "Tcells.marker"=c("CD3D", "CD3E"), "NKcells.marker"=c("GNLY"), "Monocytes.marker"=c("LYZ", "S100A8", "S100A9") )

#cellType.list = list( "Bcells" = B.cluster, "Tcells" = T.cluster, "NKcells" = NK.cluster, "Monocytes" = Monocytes.cluster ) 
cellType.list = list( "Bcells" = B.cluster, "Tcells" = T.cluster, "NKcells" = NK.cluster, "Monocytes" = c(Monocytes.cluster) ) 
z2cellType = match.z.celltype( z = cluster.pre.reZ@clusters$z, cellType.list = cellType.list )  

markerG.ori = normExp.ctype(pbmc4k_select, z.cell = z2cellType, genes = genes, type="none") 
markerG.dcon = normExp.ctype(pbmc4K_dcon, z.cell = z2cellType, genes =genes,  type="none") 
markerG.immune = normExp.ctype(mat.s, z.cell = immuneCB.re, genes = genes, type ="none") 


ori.tb = gene.observed.tb( gene.exp = markerG.ori, gene.group = gene.group) 
dcon.tb = gene.observed.tb( gene.exp = markerG.dcon, gene.group = gene.group) 
immune.tb = gene.observed.tb( gene.exp = markerG.immune, gene.group = gene.group) 

a.ori = stack.tb(ori.tb=ori.tb, cellTypes = c("Tcells", "Bcells", "NKcells", "Monocytes") ) 
a.dcon = stack.tb(ori.tb = dcon.tb, cellTypes=c("Tcells", "Bcells", "NKcells", "Monocytes") ) 
a.immune = stack.tb(ori.tb = immune.tb, cellTypes=c("Tcells", "Bcells", "NKcells", "Monocytes") )


cat.stack( a.immune) 
cat.stack( a.ori) 
cat.stack( a.dcon) 


plt.stackO = plt.stackbar(a.ori, plot.status='gt0') + ggtitle( "Original 4K PBMC counts (profiled in the same channel)" ) + xlab("Gene markers") +ylab("Percentage of cells expressing\ncell-type specific markers") + theme( plot.margin=unit(c(0, 1,0, 0) , "lines")) + scale_y_continuous(limits = c(0, 110))
plt.stackP = plt.stackbar( a.dcon, plot.status='gt0' ) + ggtitle( "Decontaminated 4K PBMC counts (profiled in the same channel)")  + xlab("Gene markers") +ylab("Percentage of cells expressing\ncell-type specific markers") + theme( plot.margin=unit(c(0, 1,0, 0) , "lines")) + scale_y_continuous(limits = c(0, 110))
plt.stackI = plt.stackbar( a.immune, plot.status = "gt0" ) + ggtitle("Sorted PBMCs counts (profiled in different channels)") + xlab("Gene markers") +ylab("Percentage of cells expressing\ncell-type specific markers")+ theme( plot.margin=unit(c(0, 1,0, 0) , "lines")) + scale_y_continuous(limits = c(0, 110))

pg.stackI = arrangeGrob( plt.empt, plt.stackI, ncol=2, widths = c(0.1,8.4) ) 
pg.stackO = arrangeGrob( plt.empt, plt.stackO, ncol=2, widths = c(0.1,8.4) ) 
pg.stackP = arrangeGrob( plt.empt, plt.stackP, ncol=2, widths = c(0.1,8.4) ) 

pdf("Fig4.pdf", width = 8.5, height = 9)
grid.arrange( p.TB.all , pg.stackI, pg.stackO, pg.stackP, ncol =1 , top="Figure 4" ) 
dev.off()

pg.stackI = arrangeGrob( plt.empt, plt.stackI + xlab("") + ylab(""), ncol=2, widths = c(0.1,8.4) ) 
pg.stackO = arrangeGrob( plt.empt, plt.stackO + xlab("") + ylab(""), ncol=2, widths = c(0.1,8.4) ) 
pg.stackP = arrangeGrob( plt.empt, plt.stackP + xlab("") + ylab(""), ncol=2, widths = c(0.1,8.4) ) 
pg.stack = arrangeGrob( pg.stackI, pg.stackO, pg.stackP, left="Percentage of cells expressing cell-type specific markers", bottom="Gene markers") 

pdf("Fig4_v3.pdf", width = 8.5, height = 9)
grid.arrange( p.TB.all , pg.stack , ncol =1 , heights = c(1, 3),  top="Figure 4" ) 
dev.off()



### Supplementary Figures >>>>>  
### Supplementary Figures >>>>>

# Gene list  
BM.Genes = c("MS4A1","CD79A","CD79B","LTB",  "S100A8","S100A9","CST3","LYZ" )  
NKT.Genes = c("CD3D","CD3E",  "NKG7","GNLY","GZMB","GZMA","CST7","TYROBP" )
MT.Genes = c("CD3D","CD3E","LTB",  "S100A8","S100A9","HLA-DRA","LYZ" )
NKB.Genes = c("MS4A1","CD79A", "CD79B","LTB",  "NKG7","GNLY","GZMB","GZMA","CST7", "TYROBP" )
NKM.Genes = c("S100A8","S100A9","LYZ",  "NKG7","GNLY","GZMB","GZMA","CST7" )


# supplementaryplots of PBMC4K data  
z.BM = label.AB( z=cluster.pre.reZ@clusters$z, A.cluster = B.cluster , B.cluster = Monocytes.cluster, A.name="Bcells", B.name="Monocytes" ) 
df.BM.gene = cluster.gene(counts=pbmc4k_select, z.AB=z.BM)


z.BM.pos = label.AB( z=cluster.pos.reZ@clusters$z, A.cluster = B.cluster , B.cluster = Monocytes.cluster, A.name="Bcells", B.name="Monocytes" ) 
df.BM.dcon = cluster.gene(counts=pbmc4K_dcon, z.AB=z.BM.pos)


z.BM.immune = label.AB( z = immuneCB.re, A.name = "Bcells", B.name = "Monocytes")
df.BM.immune = cluster.gene(counts = mat.s, z.AB = z.BM.immune ) 


df.BM.all = rbind.geneExp( geneExp.ori = df.BM.gene, geneExp.pos = df.BM.dcon, geneExp.immune = df.BM.immune) 
p.BM.all = plot.gene.exp( df.BM.all,x="log2.Monocytes", y="log2.Bcells", Genes = BM.Genes, colors = color30 ) + facet_grid(. ~ datatype_f ) + ylab("Average expression\nof B cell (log2)") + xlab("Average expression of monocyte (log2)") 



z.NKT = label.AB( z=cluster.pre.reZ@clusters$z, A.cluster = NK.cluster ,  B.cluster =T.cluster , A.name="NKcells", B.name="Tcells" ) 
df.NKT.gene = cluster.gene(counts=pbmc4k_select, z.AB=z.NKT)


z.NKT.pos = label.AB( z=cluster.pos.reZ@clusters$z, A.cluster = NK.cluster ,  B.cluster =T.cluster , A.name="NKcells", B.name="Tcells" ) 
df.NKT.dcon = cluster.gene(counts=pbmc4K_dcon, z.AB=z.NKT.pos)

z.NKT.immune = label.AB(z = immuneCB.re, A.name ="NKcells", B.name = "Tcells")
df.NKT.immune = cluster.gene(counts = mat.s, z.AB = z.NKT.immune)

df.NKT.all = rbind.geneExp( geneExp.ori = df.NKT.gene, geneExp.pos = df.NKT.dcon, geneExp.immune = df.NKT.immune) 
p.NKT.all = plot.gene.exp( df.NKT.all, x="log2.Tcells", y="log2.NKcells", Genes = NKT.Genes, colors = color30 ) + facet_grid(. ~ datatype_f ) + ylab("Average expression\nof NK cell (log2)") + xlab("Average expression of T cell (log2)") 



z.MT = label.AB( z=cluster.pre.reZ@clusters$z, A.cluster = Monocytes.cluster, B.cluster = T.cluster , A.name="Monocytes", B.name="Tcells" ) 
df.MT.gene = cluster.gene(counts=pbmc4k_select, z.AB=z.MT)

z.MT.pos = label.AB( z=cluster.pos.reZ@clusters$z, A.cluster = Monocytes.cluster, B.cluster = T.cluster , A.name="Monocytes", B.name="Tcells" ) 
df.MT.dcon = cluster.gene(counts=pbmc4K_dcon, z.AB=z.MT.pos)

z.MT.immune = label.AB( z = immuneCB.re, A.name = "Monocytes", B.name = "Tcells" ) 
df.MT.immune = cluster.gene(counts=mat.s, z.AB=z.MT.immune)

df.MT.all = rbind.geneExp( geneExp.ori = df.MT.gene, geneExp.pos = df.MT.dcon, geneExp.immune = df.MT.immune) 
p.MT.all = plot.gene.exp( df.MT.all,x = "log2.Monocytes", y = "log2.Tcells", Genes = MT.Genes, colors = color30 ) + facet_grid(. ~ datatype_f ) + ylab("Average expression\nof T cell (log2)") + xlab("Average expression of monocyte (log2)") 



z.NKB = label.AB( z=cluster.pre.reZ@clusters$z , A.cluster = NK.cluster ,  B.cluster = B.cluster , A.name="NKcells", B.name="Bcells" ) 
df.NKB.gene = cluster.gene(counts=pbmc4k_select, z.AB=z.NKB)

z.NKB.pos = label.AB( z=cluster.pos.reZ@clusters$z, A.cluster = NK.cluster ,  B.cluster = B.cluster , A.name="NKcells", B.name="Bcells" )
df.NKB.dcon = cluster.gene(counts=pbmc4K_dcon, z.AB=z.NKB.pos) 

z.NKB.immune = label.AB( z = immuneCB.re, A.name = "NKcells", B.name = "Bcells" ) 
df.NKB.immune = cluster.gene(counts=mat.s, z.AB=z.NKB.immune)

df.NKB.all = rbind.geneExp( geneExp.ori = df.NKB.gene, geneExp.pos = df.NKB.dcon, geneExp.immune = df.NKB.immune) 
p.NKB.all = plot.gene.exp( df.NKB.all, x="log2.NKcells", y="log2.Bcells", Genes = NKB.Genes, colors = color30 ) + facet_grid(. ~ datatype_f ) + ylab("Average expression\nof B cell (log2)") + xlab("Average expression of NK cell (log2)") 



z.NKM = label.AB( z=cluster.pre.reZ@clusters$z, A.cluster = NK.cluster ,  B.cluster = Monocytes.cluster , A.name="NKcells", B.name="Monocytes" ) 
df.NKM.gene = cluster.gene(counts=pbmc4k_select, z.AB=z.NKM)

z.NKM.pos = label.AB(z=cluster.pos.reZ@clusters$z, A.cluster = NK.cluster ,  B.cluster = Monocytes.cluster , A.name="NKcells", B.name="Monocytes" )
df.NKM.dcon = cluster.gene(counts=pbmc4K_dcon, z.AB=z.NKM.pos)

z.NKM.immune = label.AB( z = immuneCB.re, A.name = "NKcells", B.name = "Monocytes" ) 
df.NKM.immune = cluster.gene(counts=mat.s, z.AB=z.NKM.immune)

df.NKM.all = rbind.geneExp( geneExp.ori = df.NKM.gene, geneExp.pos = df.NKM.dcon, geneExp.immune = df.NKM.immune) 
p.NKM.all = plot.gene.exp( df.NKM.all, x="log2.NKcells", y="log2.Monocytes"  , Genes = NKM.Genes, colors = color30 ) + facet_grid(. ~ datatype_f ) + ylab("Average expression\nof monocyte (log2)") + xlab("Average expression of NK cell (log2)") 

p.list = list(p.BM.all, p.NKT.all, p.MT.all, p.NKB.all, p.NKM.all )
pdf("Sup_Fig_3.pdf", width = 8.5, height = 9 ) 
marrangeGrob( p.list,  ncol=1, nrow =3 ,  top="Supplementary Figure 3" ) 
dev.off()                             



gene.group = list( "Bcells.marker" = c("MS4A1", "CD79B", "CD79A"), "Tcells.marker"=c("CD3D", "CD3E"), "NKcells.marker"=c("GNLY"), "Monocytes.marker"=c("LYZ", "S100A8", "S100A9") )
tb.gene.pre = marker.counts.tb( pbmc4k_select, z.cell = z2cellType, gene.group = gene.group, cellType = c("Tcells",  "NKcells"), onlyMarkers = FALSE)
plt.genesUMIs.pre = plot.markerUMIs(tb.gene.pre, color = color30) + ggtitle("UMI counts of marker genes in each cell type") 


tb.gene.pos = marker.counts.tb( pbmc4K_dcon, z.cell = z2cellType, gene.group = gene.group, cellType = c("Tcells",  "NKcells"), onlyMarkers = FALSE)
plt.genesUMIs.pos = plot.markerUMIs(tb.gene.pre, color = color30) + ggtitle("UMI counts of marker genes in each cell type") 


cellType = c("Tcells",  "NKcells")
for ( cType in unique(tb.gene.pre$cellType) ) {
	cat("Original counts \n")
	cat("   cell type ", cType,  ": \n")

	markers.cellTypes = do.call( base::c, gene.group) 
	names( markers.cellTypes )  = gsub("(.*)(\\..*)", "\\1", names(markers.cellTypes) ) 

	for (gene in markers.cellTypes[ names(markers.cellTypes) %in% cellType] ) { 
	cat("          has summary statistis of the expression of gene ", gene, "\n")
	print(summary( tb.gene.pre[tb.gene.pre$cellType == cType & tb.gene.pre$genes == gene, "counts"] )) 
	}
}

cellType = c("Tcells",  "NKcells")
for ( cType in unique(tb.gene.pos$cellType) ) {
	cat("Decontaminated counts \n")
	cat("   cell type ", cType,  ": \n")

	markers.cellTypes = do.call( base::c, gene.group) 
	names( markers.cellTypes )  = gsub("(.*)(\\..*)", "\\1", names(markers.cellTypes) ) 

	for (gene in markers.cellTypes[names(markers.cellTypes) %in% cellType  ] ) { 
	cat("          has summary statistis of the expression of gene ", gene, "\n")
	print(summary( tb.gene.pos[tb.gene.pos$cellType == cType & tb.gene.pos$genes == gene, "counts"] )) 
	}
}

tb.gene.pre$Data = "Original counts"
tb.gene.pos$Data = "Decontaminated counts"
tb.gene.cb = rbind(tb.gene.pre, tb.gene.pos) 
tb.gene.cb$Data = factor(tb.gene.cb$Data, levels = c("Original counts", "Decontaminated counts") ) 
plt.genesUMIs.cb = plot.genesUMIs(tb.gene.cb[ tb.gene.cb$cellType == "NKcells", ], color = color30) + ggtitle("UMI counts of marker genes in NK-cells") + facet_wrap( .~Genes, scales="free" )


NKcell.genesUMI = tb.gene.cb[ tb.gene.cb$cellType == "NKcells", ]
for (gene in gene.group$Tcells.marker) {
  NK.geneUMItb= NKcell.genesUMI[ NKcell.genesUMI$genes == gene, ]
  res.pairedTtest= t.test(x = NK.geneUMItb$counts[NK.geneUMItb$Data == "Original counts"], y = NK.geneUMItb$counts[NK.geneUMItb$Data == "Decontaminated counts"], paired = TRUE)
	cat("Paired t test on T-cell marker gene ", gene, "'s expression in NK-cells before and after decontamination \n")
  print(res.pairedTtest)

}

pdf( "Sup_Fig_3F.pdf", width = 8.5, height = 3 )
#grid.arrange( plt.genesUMIs.pre, plt.genesUMIs.pos, ncol = 2, top = "Supplementary Figure 3")
grid.arrange(plt.genesUMIs.cb, top = "Supplementary Figure 3F")
dev.off()









## tsne plot 
plt.tF = plotDimReduceFeature(dim1=tsne.pre[,1], dim2=tsne.pre[,2],counts = pbmc4k_select, features = markerGenes ) 
pdf("Sup_Fig_4.pdf", height = 8.5, width=8.5)
grid.arrange( plt.tF, ncol = 1, top="Supplementary Figure 4")
dev.off()


a = celdaProbabilityMap(counts = pbmc4k_select, celdaMod = cluster.pre.reZ) 
pdf("Sup_Fig_5.pdf", width = 8.5, height = 8)  
grid.arrange( a, top = "Supplementary Figure 5") 
dev.off()


# scatter and Tsne plots showing correlation of the results of DecontX and Scrublet 
plt.decontx = plot.est( dim1 = tsne.pos[,1], dim2 = tsne.pos[,2], size=0.5, scaleValue= dcon$resList$estConp, varLabel = "Percentage") + ggtitle("DecontX estimated contamination level") # + theme(plot.margin=unit(c(0, 1,0, 0), "lines"), legend.margin = margin(0,0,0,0, "mm"))   
plt.scrublet = plot.est( dim1 = tsne.pos[,1], dim2 = tsne.pos[,2], size=0.5, scaleValue= doublet_scores, varLabel = "Doublet\nscore")+ ggtitle("Scrublet prediction in 4K PBMC") + theme(plot.margin=unit(c(0, 0,0, 0), "lines") )
plt.doublet = plot.doublet(dim1 = tsne.pos[,1], dim2 = tsne.pos[,2], size=0.5, doubletPre=predicted_doublets, varLabel = "Predicted\ndoublets") + ggtitle("Scrublet predicted doublets") + labs(color = "Scrublet\nprediction")  #+ theme( legend.position = "none",plot.margin=unit(c(0, 1,0, 0) , "lines")   )  

welchTTest = t.test( dcon$resList$estConp[ predicted_doublets == "singlet"] , dcon$estList$estConp[ predicted_doublets == "doublet"] , var.equal = FALSE)  
cat("p-value of the contamination proportion difference in doublets and non-doublets is: \n",  format.pval( welchTTest$p.value ), "\n") 

for( predicted_type in c("singlet", "doublet") ) {
  cat("Summary statistics of estimated contamination proportion in singlets: \n")
	print(summary( dcon$resList$estConp [ predicted_doublets == predicted_type] ))
}

plt.silhouetteHistgram = histgram.silhouette( silhouette.pre, silhouette.pos, cluster = c(2:19) ) + ggtitle("Distribution of silhouette width within each cluster") 

cellType.list = list( "Bcells" = B.cluster, "Tcells" = T.cluster, "NKcells" = NK.cluster, "Monocytes" = c(Monocytes.cluster ) ) 
z2cellType = match.z.celltype( z = cluster.pre.reZ@clusters$z, cellType.list = cellType.list )  


#plt.violinUMIs = plot.UMIsViolin(counts = pbmc4k_select, cluster = cluster.pre.reZ@clusters$z) + ggtitle("Distribution of total UMIs")
plt.violinUMIs = plot.UMIsViolin(counts = pbmc4k_select[, grepl("[^0-9]", z2cellType)], cluster = z2cellType[grepl("[^0-9]", z2cellType)], xlab = "Cell type", labs = "Cell Type", color = c("#14CC14", "#8514CC", "#14CCCC", "#E6A52E") ) + ggtitle("Distribution of total UMIs") + scale_y_continuous( limits = c(0, 30000), breaks = seq(0,20000, 5000) ) + geom_text( aes(label = medianUMIs, y = 29000) , size = 4) + geom_text(aes(label = UMIs.pct, y = 24000), size = 3) 
plt.violinDcon = plot.estConpViolin( estConp = dcon$resList$estConp, cluster = cluster.pre.reZ@clusters$z) + ggtitle("Distribution of contamination levels")


gene.group = list( "Bcells.marker" = c("MS4A1", "CD79B", "CD79A"), "Tcells.marker"=c("CD3D", "CD3E"), "NKcells.marker"=c("GNLY"), "Monocytes.marker"=c("LYZ", "S100A8", "S100A9") )
tb.melt.pre = marker.counts.tb( pbmc4k_select, z.cell = z2cellType, gene.group = gene.group, cellType = c("Tcells", "Bcells", "NKcells", "Monocytes")) 
plt.markerUMIs.pre = plot.markerUMIs(tb.melt.pre, color = color30) + ggtitle("UMI counts of marker genes in each cell type") 
#ggsave(plot = plt.markerUMIs.pre, file = "marker_ori.pdf") 


for ( cType in unique(tb.melt.pre$cellType) ) {
	cat("cell type ", cType,  ": \n")

	markers.cellTypes = do.call( base::c, gene.group) 
	names( markers.cellTypes )  = gsub("(.*)(\\..*)", "\\1", names(markers.cellTypes) ) 

	for (gene in markers.cellTypes[ names(markers.cellTypes) == cType] ) { 
	cat("          has summary statistis of the expression of gene ", gene, "\n")
	print(summary( tb.melt.pre[ tb.melt.pre$cellType == cType & tb.melt.pre$genes == gene, "counts"] )) 
	}
}



gdecontx = ggplot_gtable(ggplot_build(plt.decontx + theme(legend.position = "bottom")))
glegend.decontx = gdecontx$grob[[which(sapply(gdecontx$grobs, function(x) x$name) == "guide-box")]]

gdoublet= ggplot_gtable(ggplot_build(plt.doublet + theme(legend.position = "bottom")))
glegend.doublet =  gdoublet$grob[[which(sapply(gdoublet$grobs, function(x) x$name) == "guide-box")]] 


gdmarkerUMI = ggplot_gtable(ggplot_build(plt.markerUMIs.pre ))
plegend.markerUMIs = gdmarkerUMI$grob[[which(sapply(gdmarkerUMI$grobs, function(x) x$name) == "guide-box")]]

pg.violinUMIs = arrangeGrob( plt.violinUMIs + theme(legend.position = "none"), plt.markerUMIs.pre + theme(legend.position = "none"), plegend.markerUMIs,   ncol = 3, widths =c(4, 4, 2) )  
pg.estConp = arrangeGrob(plt.violinDcon + theme(legend.position = "none") )
pg.sil = arrangeGrob( plt.silhouetteHistgram ) 
pg.cmb = arrangeGrob( plt.doublet + theme(legend.position = "none"), plt.decontx + theme(legend.position = "none") , plt.violinDcon + theme(legend.position = "none"),   ncol = 3, widths = c(3, 3, 5) ) 
pg.legend = arrangeGrob(  glegend.doublet, glegend.decontx,  ncol = 2, widths = c(1, 1 )) 

for (cell in c("Monocytes", "Tcells", "Bcells", "NKcells")) {
	counts.cell = pbmc4k_select[, z2cellType == cell]
	cat("cell type ", cell, "has ", ncol(counts.cell), "cells \n", "       summary statistics of total UMIs is \n")
	print(summary(colSums(counts.cell)))
	cat("    UMIs % for this cell type is", sum(counts.cell) / sum(pbmc4k_select), "\n" )
}


pdf("Sup_Fig_6.pdf", width = 8.5, height = 8.5) 
grid.arrange(pg.violinUMIs,  pg.sil,  pg.cmb , pg.legend,   ncol = 1, heights = c(2.5,5, 3, 0.6),  top = "Supplementary Figure 6") 
dev.off()





