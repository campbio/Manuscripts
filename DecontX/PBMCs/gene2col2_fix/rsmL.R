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
class(pbmc4k_select) = "integer"
dim(pbmc4k_select)


set.seed(12345) 
rsmL = recursiveSplitModule( counts = pbmc4k_select, maxL=150) 
saveRDS(rsmL, file.path(celda.wd, 'rsm.L.celda.rds') ) 
plot.rsmL = plotGridSearchPerplexity(celdaList = rsmL ) 
ggsave(plot =  plot.rsmL, filename = file.path( celda.wd, 'rsm.L.celda.pdf' ) ) 






