wd <- "/restricted/projectnb/camplab/home/syyang/contamination/data/pbmc/4k" 
celda.wd <- "/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/pbmc-ana/4k/gene2col2"

#library(celda)
library(devtools)
load_all("~/gitProjects/celda", recompile=T)
library(Matrix)
library(SummarizedExperiment)
library(MAST)
library(ggplot2) 

data <- readMM( file.path(wd, "data/matrix.mtx")  )
barcodes <- read.table( file.path(wd, "data/barcodes.tsv")  , sep = "\t" )
genes <- read.table( file.path(wd, "data/genes.tsv")  ,sep = "\t"  )

pbmc4k <- as.matrix(data)
rownames(pbmc4k) <- paste(genes$V1, "_", genes$V2, sep = "")
colnames(pbmc4k) <- barcodes$V1

pbmc4k_select <- pbmc4k[rowSums(pbmc4k>2) > 2,]
dim(pbmc4k_select)


set.seed(12345) 
rsmL = recursiveSplitModule( counts = pbmc4k_select, maxL=200) 
saveRDS(rsmL, file.path(celda.wd, 'rsm.L.celda_yusukeCode.rds') ) 
plot.rsmL = plotGridSearchPerplexity(celdaList = rsmL ) 
ggsave(plot =  plot.rsmL, filename = file.path( celda.wd, 'rsm.L.celda_yusukeCode.pdf' ) ) 


module.split.select = subsetCeldaList(rsmL, list(L=150))
rsmK = recursiveSplitCell(pbmc4k_select, initialK = 3, maxK=40, yInit =clusters(module.split.select)$y)
saveRDS(rsmK, "rsm.L150.K.celda_yusukeCode.rds") 
rsmK = readRDS("rsm.L150.K.celda_yusukeCode.rds") 

plot.rsmK = plotGridSearchPerplexity(rsmK) 
ggsave( plot = plot.rsmK, "rsm.L150.K.celda_yusukeCode.pdf") 


res.19 = subsetCeldaList( rsmK, params = list(K = 19, L = 150))
dcon.K19 = decontX( counts = pbmc4k_select, z = res.19@clusters$z ) 
saveRDS( dcon.K19, "dcon.L150K19.rds" )

counts.dcon = round( dcon.K19$resList$estNativeCounts)
set.seed(12345) 
p.LK = celda_CG( counts = counts.dcon, zInit = res.19@clusters$z, K = 19, L=150 )
saveRDS( p.LK, "post.cluster.L150K19.rds" )


table( res.19@clusters$z, p.LK@clusters$z ) 


set.seed(12345) 
tsne.K19 = celdaTsne( counts = pbmc4k_select, celdaMod = res.19 ) 
plt.1 = plotDimReduceCluster( dim1 = tsne.K19[,1], dim2 = tsne.K19[, 2], cluster = res.19@clusters$z, labelClusters = T)
ggsave( plot = plt.1, "tsne.L150K19.pdf") 

set.seed(12345) 
tsne.P = celdaTsne( counts = counts.dcon , celdaMod = p.LK ) 
plt.2 = plotDimReduceCluster( dim1 = tsne.P[,1], dim2 = tsne.P[,2], cluster = p.LK@clusters$z, labelClusters = T ) 
ggsave( plot = plt.2, "tsne.p.L150K19.pdf" ) 
