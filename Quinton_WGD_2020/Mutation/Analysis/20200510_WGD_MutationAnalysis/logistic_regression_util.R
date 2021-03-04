## Function to test each row of a matrix with the FET
logistic.regression = function(mat, formula, covar, verbose=TRUE) {
  
  temp.model = model.matrix(formula, data=covar)
  nc = ncol(temp.model)
  pval = matrix(NA, ncol=nc, nrow=nrow(mat), dimnames=list(rownames(mat), colnames(temp.model)))
  estimate = matrix(NA, ncol=nc, nrow=nrow(mat), dimnames=list(rownames(mat), colnames(temp.model)))
  zval = matrix(NA, ncol=nc, nrow=nrow(mat), dimnames=list(rownames(mat), colnames(temp.model)))
  
  for(i in 1:nrow(mat)) {
  
    ## Set up data.frame and formula
    try({
      temp.data=data.frame(expr=mat[i,]>0, covar)
      new.formula = as.formula(paste("expr ~ ", as.character(formula)[2], sep=""))
  
      ## Run GLM
      g = glm(formula=new.formula, data=temp.data, family="binomial", maxit=100)
    
      ## Save results in matrices
      s = summary(g)$coefficients
      pval[i,] = s[,4]
      estimate[i,] = s[,1]
      zval[i,] = s[,3]
       
    }, FALSE)
    
    if(verbose == TRUE & i %% 100 == 0) {
      cat("Iteration:", i, "\n")
    }
  }
  
  fdr = apply(pval, 2, p.adjust, 'fdr')
  return(list(Pvalue=pval, Estimate=estimate, Zvalue=zval, FDR=fdr))
}  

