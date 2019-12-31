wd <- "/restricted/projectnb/camplab/home/syyang/contamination/data/pbmc/4k/" 

library(Matrix)
library(SummarizedExperiment)

data <- readMM(paste0(wd, "data/matrix.mtx")  )
barcodes <- read.table(paste0(wd, "data/barcodes.tsv")  , sep = "\t" )
GENES <- read.table(paste0(wd, "data/genes.tsv")  ,sep = "\t"  )

pbmc4k <- as.matrix(data)
rownames(pbmc4k) <- paste(GENES$V1, "_", GENES$V2, sep = "")
colnames(pbmc4k) <- barcodes$V1

pbmc4k_select <- pbmc4k[rowSums(pbmc4k>2) >2, ] 
class(pbmc4k_select) = "integer"
dim(pbmc4k_select)


#data.duplicate = readMM("/restricted/projectnb/camplab/home/syyang/contamination/data/pbmc/4k_duplicate/filtered_gene_bc_matrices/GRCh38/matrix.mtx")
#sum( data != data.duplicate )  
print( " 4k_duplicate is the same dataset as 4k " ) 
