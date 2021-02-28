source("~/GIT/utilities/R/lm_utils.R")
source("~/GIT/utilities/R/mut_utils.R")
library(stringr)

ceres.17q2 <- read.table(gzfile("../../Data/gene_effect_17Q2.csv.gz"), header = T, sep = ",", row.names=1, check.names=F)
ceres.18q3 <- t(read.table(gzfile("../../Data/gene_effect_18Q3.csv.gz"), header = T, sep = ",", row.names=1, check.names=F))
ceres.combined.18q3 <- read.table(gzfile("../../Data/D2_combined_gene_dep_scores_18Q3.csv.gz"), header = T, sep = ",", row.names=1, check.names=F)
ceres.combined.18q3 = ceres.combined.18q3[rowSums(is.na(ceres.combined.18q3)) < 5,]
Absolute_data <- read.table("../../Data/CCLE_combined.table.txt", header = T, sep = "\t", row.names=1)


enrich = function(ceres.data, absolute, cutoff = -1) {
  require(pROC)
  i = intersect(colnames(ceres.data), absolute[,1])
  abs.i = Absolute_data[i,]
  wgd = ifelse(abs.i$Genome.doublings > 0, "WGD", "Not_WGD")

  ceres.i = ceres.data[,i]

  res.fet.full = fet(ceres.i < cutoff, wgd, reorder=FALSE)
  colnames(res.fet.full) = paste0("FET_", colnames(res.fet.full))
  res.wilcox.full = wilcoxon(ceres.i, wgd, sort=FALSE)  
  colnames(res.wilcox.full) = paste0("Wilcoxon_", colnames(res.wilcox.full))
    
  ceres.i.select = rowSums(ceres.i < cutoff) > 4
  fet.select.fdr = res.fet.full$FET_Pvalue
  fet.select.fdr[!ceres.i.select] = NA
  fet.select.fdr = p.adjust(fet.select.fdr, 'fdr')
  wilcox.select.fdr = res.wilcox.full$Wilcoxon_Pvalue
  wilcox.select.fdr[!ceres.i.select] = NA
  wilcox.select.fdr = p.adjust(wilcox.select.fdr, 'fdr')
  
  auc = rep(NA, nrow(ceres.i))
  for(j in 1:nrow(ceres.i)) {
    auc[j] = auc(roc(as.factor(wgd), as.numeric(ceres.i[j,])))
  }
  res = data.frame(res.fet.full, "FET_FDR_Filter"=fet.select.fdr, AUC=auc, res.wilcox.full, "Wilcoxon_FDR_Filter"=wilcox.select.fdr)
  return(list(res=res, data=rbind(WGD=wgd, ceres.i)))
}


getTumorType = function(cn) {
  s = str_split(cn, "_", simplify=T)
  s2 = s[,-1]
  s3 = apply(s2, 1, paste, collapse="_")
  s4 = gsub("_+$", "", s3)
  return(s4)
}

ceres.17q2.tt = getTumorType(colnames(ceres.17q2))
ceres.18q3.tt = getTumorType(colnames(ceres.18q3))
ceres.combined.18q3.tt = getTumorType(colnames(ceres.combined.18q3))

tumor.type = unique(c(ceres.17q2.tt, ceres.18q3.tt, ceres.combined.18q3.tt))
min.n = 20
for(i in tumor.type) { 
  print(i)
  ix = ceres.17q2.tt == i
  cn = intersect(colnames(ceres.17q2)[ix], rownames(Absolute_data))
  if(length(cn) >= min.n) {
    print("CERES 17q2")
	res.ceres.17q2 = enrich(ceres.17q2[,cn], Absolute_data[cn,])
	fn = paste0("20181220_CERES_", i, "_17Q2_WGD_results.txt")
	write.table(res.ceres.17q2$res, fn, quote=FALSE, row.names=FALSE, sep="\t")
	fn = paste0("20181108_CERES_", i, "_17Q2_WGD_data.txt")
	write.table(data.frame(Gene=rownames(res.ceres.17q2$data), res.ceres.17q2$data), fn, quote=FALSE, row.names=FALSE, sep="\t")
  }
  ix = ceres.18q3.tt == i
  cn = intersect(colnames(ceres.18q3)[ix], rownames(Absolute_data))
  if(length(cn) >= min.n) {  
    print("CERES 18q3")  
	res.ceres.18q3 = enrich(ceres.18q3[,cn], Absolute_data[cn,])
	fn = paste0("20181220_CERES_", i, "_18Q3_WGD_results.txt")
	write.table(res.ceres.18q3$res, fn, quote=FALSE, row.names=FALSE, sep="\t")
	fn = paste0("20181220_CERES_", i, "_18Q3_WGD_data.txt")  
	write.table(data.frame(Gene=rownames(res.ceres.18q3$data), res.ceres.18q3$data), fn, quote=FALSE, row.names=FALSE, sep="\t")
  }
  ix = ceres.combined.18q3.tt == i
  cn = intersect(colnames(ceres.combined.18q3)[ix], rownames(Absolute_data))
  if(length(cn) >= min.n) {
    print("CERES Combined 18q3")
    res.ceres.combined = enrich(ceres.combined.18q3[,cn], Absolute_data[cn,])
    fn = paste0("20181108_D2combined_", i, "_18Q3_WGD_results.txt")
    write.table(res.ceres.combined$res, fn, quote=FALSE, row.names=FALSE, sep="\t")
    fn = paste0("20181220_D2combined_", i, "_18Q3_WGD_data.txt")  
    write.table(data.frame(Gene=rownames(res.ceres.combined$data), res.ceres.combined$data), fn, quote=FALSE, row.names=FALSE, sep="\t")
  }  
}


res.ceres.17q2 = enrich(ceres.17q2, Absolute_data)
res.ceres.18q3 = enrich(ceres.18q3, Absolute_data)
res.ceres.combined = enrich(ceres.combined.18q3, Absolute_data)

write.table(res.ceres.17q2$res, "20181220_CERES_17Q2_WGD_results.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(data.frame(Gene=rownames(res.ceres.17q2$data), res.ceres.17q2$data), "20181220_CERES_17Q2_WGD_data.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(res.ceres.18q3$res, "20181220_CERES_18Q3_WGD_results.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(data.frame(Gene=rownames(res.ceres.18q3$data), res.ceres.18q3$data), "20181220_CERES_18Q3_WGD_data.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(res.ceres.combined$res, "20181220_D2combined_18Q3_WGD_results.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(data.frame(Gene=rownames(res.ceres.combined$data), res.ceres.combined$data), "20181220_D2combined_18Q3_WGD_data.txt", quote=FALSE, row.names=FALSE, sep="\t")


sessionInfo()
