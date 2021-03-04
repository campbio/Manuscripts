lm_expression_with_cnv <- function(expression, cnv, covars) {
  
  # Convert to matrix if not already
  expression = as.matrix(expression)
  cnv = as.matrix(cnv)
  
  #result <- matrix(NA, ncol = (ncol(covars) + 2) * 5, nrow = nrow(expression) )
  var.names = c("Expression", "CNV", colnames(covars))
  form = as.formula(paste0("Expression ~ CNV + ", paste(colnames(covars), collapse = "+")))
  coef = c()
  sterr = c()
  pval = c()
  tstat = c()
  for(i in 1:nrow(expression)) {
    data = cbind(Expression=expression[i,], CNV=cnv[i,], covars)
    error = try({
    model = lm(formula = form, data = data)
    }, silent = TRUE)
    
    if(class(error) == "try-error") {
      coef = cbind(coef, NA)
      sterr = cbind(sterr, NA)
      tstat = cbind(tstat, NA)
      pval = cbind(pval, NA)
      
    } else {
      s = summary(model)$coefficients  
      coef = cbind(coef, s[,1])
      sterr = cbind(sterr, s[,2])
      tstat = cbind(tstat, s[,3])
      pval = cbind(pval, s[,4])
    }
    
    if(i %% 1000 == 0) {
      cat("Analysis of", i, "genes completed\n")
    }
  }
  fdr = apply(pval, 1, p.adjust, method="fdr")
  rownames(coef) = paste0(var.names, "_Estimate")
  rownames(sterr) = paste0(var.names, "_StdError")
  rownames(tstat) = paste0(var.names, "_Tstat")
  rownames(pval) = paste0(var.names, "_Pvalue")
  colnames(fdr) = paste0(var.names, "_FDR")
  
  result = data.frame(Gene = rownames(expression),
                      t(coef),t(sterr),t(tstat),t(pval),fdr)
  return(result)
}