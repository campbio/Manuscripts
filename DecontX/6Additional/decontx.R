RNGkind(sample.kind = "Rounding")
library(devtools)
load_all("~/gitProjects/celda_irisapo", recompile=T)   # on fix branch, celda version 1.3.1



autoDecontX = function(sce, varGenes=5000, seed=12345, L=50, dbscan.eps=1, decontxIter=200) {
    require(Matrix) 
    require(dbscan)
    require(uwot)
    require(scater)
    require(scran)
    
    sce = scater::normalizeSCE(sce)  # add the log2 normalized counts into sce object 

    if( nrow(sce) <= varGenes) {
	  topVariableGenes = 1:nrow(sce)
    } else if( nrow(sce) > varGenes ) { 
	  # Use the top most variable genes to do rough clustering (celda_CG & Umap & DBSCAN graph algorithm) 
	  mvTrend = scran::trendVar(sce, use.spikes=FALSE) 
	  decomposeTrend = scran::decomposeVar(sce, mvTrend) 
	  topVariableGenes = order(decomposeTrend$bio, decreasing=TRUE)[1:varGenes]
    } 
	countsFiltered = as.matrix(counts(sce[topVariableGenes, ]))
    storage.mode(countsFiltered) = "integer"
    
    # Celda clustering using recursive module splitting
    L = min(L, nrow(countsFiltered))    
		RNGkind(sample.kind = "Rounding")
		set.seed(12345)
   	initial.module.split = recursiveSplitModule(countsFiltered, initialL=L, maxL=L, perplexity=FALSE, verbose=FALSE)
    initial.modules.model = subsetCeldaList(initial.module.split, list(L=L))
	  fun_wd = getwd()
		cat("current working directory is \n", fun_wd)
	  saveRDS(initial.module.split, file.path(fun_wd, "celda_obj.rds") )

    # UMAP dimension reduction
 		RNGkind(sample.kind = "Rounding")
		set.seed(12345)
    fm = factorizeMatrix(countsFiltered, initial.modules.model, type="counts")
    res.umap = uwot::umap(t(sqrt(fm$counts$cell)), n_neighbors=15, min_dist = 0.01, spread = 1)
    reducedDim(sce, "QC_UMAP") = res.umap
    
    # Use dbSCAN on the UMAP to identify broad cell types
    totalClusters = 1
    while(totalClusters <= 1 & dbscan.eps > 0) {
      res.dbscan = dbscan(res.umap, dbscan.eps)
      dbscan.eps = dbscan.eps - (0.25 * dbscan.eps)
      totalClusters = length(unique(res.dbscan$cluster))
    }  
    
    # If the counts is not matrix in the sce object, have to covnert the count into matrix to feed into celda_CG 
    rm(countsFiltered)
    counts = counts(sce) 
    counts.classType = class(counts)     
	counts = as.matrix(counts) 
	storage.mode(counts) = "integer"    
	
    if( totalClusters >=2) {
	  # Perform decontamination with decontX 
	  DecontXObject = decontX(counts=counts, z=res.dbscan$cluster, seed=seed, maxIter=decontxIter)
	  DecontXObject$resList$estNativeCounts <- round(DecontXObject$resList$estNativeCounts) 

	  if( counts.classType == "matrix") {
		  SummarizedExperiment::assay(sce, "decontX") <- DecontXObject$resList$estNativeCounts 
	  } else if (counts.classType == "dgCMatrix") {
	  SummarizedExperiment::assay(sce, "decontX") <- as(DecontXObject$resList$estNativeCounts, "dgCMatrix") 
	  } else {
		  SummarizedExperiment::assay(sce, "decontX") <- DelayedArray::DelayedArray(DecontXObject$resList$estNativeCounts)
	  }

	  colData(sce)$DecontX_Contamination = DecontXObject$resList$estConp
      colData(sce)$DecontX_QuickCluster = res.dbscan$cluster	
    } else {
	  colData(sce)$DecontX_Contamination = 0
	  colData(sce)$DecontX_QuickCluster = 1	    
    }
    
    return(sce) 
}
