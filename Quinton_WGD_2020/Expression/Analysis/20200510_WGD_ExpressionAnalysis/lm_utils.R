
residual.matrix <- function(the.matrix, the.model, the.coefficients=NULL) {
  # Correct for covariates using a linear model
  # Note that no "NA"s can be in the model matrix
  # Code from Marc Lenburg, 2010
  
  if(is.null(the.coefficients)) {
    # Solves the coefficients
    resid.data <- as.matrix(the.matrix) %*% (diag(dim(the.matrix)[2]) - the.model %*% solve(t(the.model) %*% the.model) %*% t(the.model))
  } else {
    # Uses user-defined coeffiencts
    resid.data <- as.matrix(the.matrix) - t(the.model %*% t(the.coefficients))
  }
  
  rownames(resid.data) <- rownames(the.matrix)
  colnames(resid.data) <- colnames(the.matrix)
  return(resid.data)
}


bwr = colorRampPalette(c("blue", "white", "red"))(100)

intersectSeveral <- function(...) { Reduce(intersect, list(...)) } 

## Custom functions
ttest = function(d, model, model2=NULL) {
  require(limma)
  require(sva)
  fit = lmFit(d, model)
  fit.t = fit$coef / fit$stdev.unscaled / fit$sigma
  fit.t.p = 2*(1-pt(abs(fit.t), fit$df.residual))
  fit.t.p.fdr = apply(fit.t.p, 2, p.adjust, method="fdr")
  
  f.p=NULL
  if(!is.null(model2)) {
    f.p = f.pvalue(d, model, model2)
  }
  
  return(list(Log2_FC=fit$coefficients, Tstat=fit.t, Pvalue=fit.t.p, FDR=fit.t.p.fdr, ANOVA.p=f.p))
}

gsea.rnk = function(genes, stat, filename) {
  na.ind = is.na(stat)
  new.rank = aggregate(stat[!na.ind], by=list(genes[!na.ind]), mean)
  write.table(new.rank, filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

is.outlier = function(v, coef=3) {
  outlier.q = as.numeric(quantile(v, c(0.25, 0.75), na.rm=TRUE))
  outlier.range = coef*abs(outlier.q[2] - outlier.q[1])
  outlier = (v > (outlier.q[2] + outlier.range)) | (v < (outlier.q[1] - outlier.range))
  return(outlier)
}

spearman = function(m, p) {
  res = c()
  for(i in 1:nrow(m)) {
    s = cor.test(as.numeric(m[i,]), as.numeric(p), method="spearman")
    res = cbind(res, c(Log2_FC=as.numeric(s$estimate), Tstat=as.numeric(s$statistic), Pvalue=as.numeric(s$p.value)))
  }
  res = t(rbind(res, FDR=p.adjust(res["Pvalue",], 'fdr')))
  return(res)
}

wilcoxon = function(m, p, sort=TRUE) {
  
  if(length(unique(p)) != 2) {
    stop("The vector p can only have 2 levels")
  }
  
  res = c()
  for(i in 1:nrow(m)) {
    s = wilcox.test(as.numeric(m[i,]) ~ p)
    a = aggregate(as.numeric(m[i,]), by=list(p), median, na.rm = TRUE)

    temp = c(a[,2], as.numeric(s$statistic), as.numeric(s$p.value))
    names(temp) = c(paste(a[1,1], "_Median", sep=""), paste(a[2,1], "_Median", sep=""), "Statistic", "Pvalue")
    
    res = rbind(res, temp)
  }
  rownames(res) = NULL
  res = data.frame(Gene=rownames(m), res, FDR=p.adjust(res[,"Pvalue"], 'fdr'), stringsAsFactors=FALSE, check.names=FALSE)
  
  if(sort == TRUE) {
    res = res[order(res$Pvalue),]
  }
 
  return(res)
}

gsea.gs.fdr.fc.cutoff = function(genes, fdr, fc, fdr.cutoff=0.05, fc.cutoff=1, prefix="gene_set") {
  na.ind = is.na(genes) | is.infinite(genes) | genes == ""
  genes = genes[!na.ind]
  fdr = fdr[!na.ind]
  fc = fc[!na.ind]
  
  ind.up = fdr < fdr.cutoff & fc > fc.cutoff
  ind.down = fdr < fdr.cutoff & fc < -fc.cutoff
  
  filename = paste(prefix, ".gmt", sep="")
  gs.up.name = paste(prefix, "_UP", sep="")
  gs.down.name = paste(prefix, "_DOWN", sep="")
   
  write.table(matrix(c(gs.up.name, "na", unique(genes[which(ind.up)])), nrow=1), filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(matrix(c(gs.down.name, "na", unique(genes[which(ind.down)])), nrow=1), filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
}


gsea.gs.fdr.fc.top = function(genes, fdr, fc, fdr.cutoff=0.05, top=200, prefix="gene_set") {
  ind = which(!(is.na(genes) | is.infinite(genes) | genes == "" | is.na(fc) | is.infinite(fc)) & fdr < fdr.cutoff & fc > 0)
  genes.up = genes[ind]
  fc.up = fc[ind]
  fc.o1 = order(fc.up, decreasing=TRUE)
  up.id = genes.up[head(fc.o1, n=top)]
  
  ind = which(!(is.na(genes) | is.infinite(genes) | genes == "" | is.na(fc) | is.infinite(fc)) & fdr < fdr.cutoff & fc < 0)
  genes.down = genes[ind]
  fc.down = fc[ind]
  fc.o1 = order(fc.down, decreasing=FALSE)
  down.id = genes.down[head(fc.o1, n=top)]
  

  filename = paste(prefix, ".gmt", sep="")
  gs.up.name = paste(prefix, "_UP", sep="")
  gs.down.name = paste(prefix, "_DOWN", sep="")

  write.table(matrix(c(gs.up.name, "na", up.id), nrow=1), filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(matrix(c(gs.down.name, "na", down.id), nrow=1), filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
}

gsea.gs = function(up, down, prefix="gene_set") {  
  
  filename = paste(prefix, ".gmt", sep="")
  gs.up.name = paste(prefix, "_UP", sep="")
  gs.down.name = paste(prefix, "_DOWN", sep="")
  
  write.table(matrix(c(gs.up.name, "na", up), nrow=1), filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(matrix(c(gs.down.name, "na", down), nrow=1), filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
}

gsea.multi.gs = function(gs.list, prefix="gene_set") {  
  filename = paste(prefix, ".gmt", sep="")
  unlink(filename)
  n = names(gs.list)
  
  for(i in 1:length(gs.list)) {
    gs.name = paste(prefix, n[i], sep="_")
    write.table(matrix(c(gs.name, "na", gs.list[[i]]), nrow=1), filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}




lme.test = function(gene.matrix, fixed, random, covar, lme.method="ML", verbose=TRUE) {

  require(nlme)
  
  ## Set up an initial formula with the first gene 
  gene.matrix = as.matrix(gene.matrix, data=covar)
  gene = gene.matrix[1,]
  form.fixed =  as.formula(paste(c("gene", as.character(fixed)), collapse=" "))
  
  ## Set up the results matrices
  m = model.matrix(form.fixed)
  lme.tstat <- matrix(NA, nrow=nrow(gene.matrix), ncol=ncol(m), dimnames=list(rownames(gene.matrix), colnames(m)))
  lme.coef <- matrix(NA, nrow=nrow(gene.matrix), ncol=ncol(m), dimnames=list(rownames(gene.matrix), colnames(m)))
  lme.pval <- matrix(NA, nrow=nrow(gene.matrix), ncol=ncol(m), dimnames=list(rownames(gene.matrix), colnames(m)))

  ## Iterate through each gene and apply the model 
  for(i in 1:nrow(gene.matrix)) {

    ## Get new gene and set up the formula
    gene = gene.matrix[i,]
    new.data = data.frame(gene, covar)
    form.fixed =  as.formula(paste(c("gene", as.character(fixed)), collapse=" "))

    ## Use try in case the model fails
    try({

      model <- lme(fixed=form.fixed, random=random, data=new.data, method=lme.method)
      s <- summary(model)$tTable
      lme.tstat[i,] = s[,4]
      lme.pval[i,] = s[,5]
      lme.coef[i,] = s[,1]
    
    }, TRUE)
    
    if(verbose == TRUE & i %% 1000 == 0) {
      cat("Iteration:", i, "\n")
    }
  }     
  
    fdr = apply(lme.pval, 2, p.adjust, method="fdr") 
    return(list(Log2_FC=lme.coef, Tstat=lme.tstat, Pvalue=lme.pval, FDR=fdr))
}

