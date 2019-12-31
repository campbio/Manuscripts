# load data 
load("/restricted/projectnb/camplab/home/syyang/contamination/data/cellBenchData/CellBench_data_08112019/data/mRNAmix_qc.RData")

cell.lines = c("H2228_prop", "H1975_prop", "HCC827_prop")

#   sce2_qc  (dataset) is a mixture of mRNA from "H2228", "H1975" and "HCC827" 3 cell lines in different proportions. 




counts.n.z = function(sce, cellLine = cell.lines) { 
    cellLine_prop = data.frame(H2228_prop = sce@colData@listData$H2228_prop,
			       H1975_prop = sce@colData@listData$H1975_prop,
			       HCC827_prop = sce@colData@listData$HCC827_prop) 

    # filter "cell"s that have one dominant cell line ( >= 50% mRNA ) 
    filter.colIndx = apply( cellLine_prop, 1, FUN=function(i) {max(i) > 0.5} )
    counts = sce@assays$data$counts[, filter.colIndx]

    z_dominant = colnames(cellLine_prop)[apply(cellLine_prop, 1, FUN=function(i) {which.max(i)})] 
    z_dominant = z_dominant[filter.colIndx]
    z.int = plyr::mapvalues(x=z_dominant, from = cellLine, to = 1:length(cellLine) ) 
    z.int = as.integer(z.int) 
    return(list("counts" = counts, "z_dominant"= z_dominant, "z" = z.int, "cellLine_prop"=cellLine_prop, "filter.colIndx" = filter.colIndx) )   
}


get.sce2_qc= counts.n.z( sce = sce2_qc)
saveRDS(get.sce2_qc, "get.sce2_qc.rds") 

# Run DecontX
RNGkind(sample.kind = "Rounding")
library(devtools)
load_all("~/gitProjects/celda_irisapo", recompile=T)   # on fix branch, celda version 1.3.1


set.seed(12345)
dcon = decontX( counts = get.sce2_qc$counts , z = get.sce2_qc$z) 
saveRDS( dcon, "dcon.sce2_qc.rds" )  
