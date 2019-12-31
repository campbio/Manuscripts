wd <- "/restricted/projectnb/camplab/home/syyang/contamination/data/pbmc/4k_duplicate/raw_gene_bc_matrices/" 

library(Matrix)
library(SummarizedExperiment)

rawData <- readMM(paste0(wd, "GRCh38/matrix.mtx")  )
rawBarcodes <- read.table(paste0(wd, "GRCh38/barcodes.tsv")  , sep = "\t" )
rawGenes <- read.table(paste0(wd, "GRCh38/genes.tsv")  ,sep = "\t"  )

rownames(rawData) <- paste(rawGenes$V1, "_", rawGenes$V2, sep = "")
colnames(rawData) <- rawBarcodes$V1
pbmcRaw <- as.matrix(rawData)



#library(DropletUtils) 
#set.seed(12345) 
#drop.out = emptyDrops( rawData) 

#saveRDS(drop.out, "out.EmptDrops.rds") 

