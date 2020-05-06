#### Compare met and local dataset
#### HS 2020
library(pwr)
library(tibble)
library(ggplot2)
library(tidyr)
library(fBasics)
library(rapportools)

file_direc="PCa/"
file_type="met"

raw_met <- read.csv("PCa/met.csv", header=TRUE)
raw_local <-read.csv("PCa/local.csv", header=TRUE)
raw_inter <- raw_local[0,]
raw_p <- raw_local[0,5:10]
colnames(raw_p)<-c('aa_freq_local','aa_freq_met','eur_freq_local','eur_freq_met','p_aa','p_eur')
colnames(raw_inter)<-c("gene",'class','aadenom_local','aadenom_met','eurdenom_local','eurdenom_met','aacount_local','aacount_met','eurcount_local','eurcount_met')
m=1
for (i in 1:nrow(raw_local)) {
  raw_met_match=raw_met[raw_met$gene==lapply(raw_local[i,1],as.character) & raw_met$class==raw_local[i,2],]
  if(!is.empty(raw_met_match)[1]) {
    raw_inter[m,1]=lapply(raw_local[i,1],as.character)
    raw_inter[m,2]=raw_local[i,2]
    
    raw_inter[m,7]=raw_local[i,7]
    raw_inter[m,9]=raw_local[i,8]
    raw_inter[m,3]=raw_local[i,5]
    raw_inter[m,5]=raw_local[i,6]
    
    raw_inter[m,8]=raw_met_match$aa_count
    raw_inter[m,10]=raw_met_match$eur_count
    raw_inter[m,4]=raw_met_match$aa_denom
    raw_inter[m,6]=raw_met_match$eur_denom

    list_aa <- data.frame(c(raw_inter[m,7],raw_inter[m,3]-raw_inter[m,7]),
                               c(raw_inter[m,8],raw_inter[m,4]-raw_inter[m,8]))
    list_eur <- data.frame(c(raw_inter[m,9],raw_inter[m,5]-raw_inter[m,9]),
                             c(raw_inter[m,10],raw_inter[m,6]-raw_inter[m,10]))
    prefish_aa <-as.matrix(sapply(list_aa, as.numeric)) 
    prefish_eur <-as.matrix(sapply(list_eur,as.numeric))
    
    raw_p[m,1]=raw_local[i,4]
    raw_p[m,3]=raw_local[i,3]
    raw_p[m,2]=raw_met_match$aa_freq
    raw_p[m,4]=raw_met_match$eur_freq
    
    raw_p[m,5]=fisher.test(prefish_aa,alternative="two.sided")[1]
    raw_p[m,6]=fisher.test(prefish_eur,alternative="two.sided")[1]
    
    m=m+1
  }
}

raw_inter$aa_freq_local=raw_p$aa_freq_local
raw_inter$aa_freq_met=raw_p$aa_freq_met
raw_inter$eur_freq_local=raw_p$eur_freq_local
raw_inter$eur_freq_met=raw_p$eur_freq_met
raw_inter$p_aa=raw_p$p_aa
raw_inter$p_eur=raw_p$p_eur
raw_inter$p_BH_aa=p.adjust(raw_p$p_aa,method="BH")
raw_inter$p_BH_eur=p.adjust(raw_p$p_eur,method="BH")
write.csv(raw_inter, file = "PCa/Inter.csv",row.names=FALSE, na="")

raw_data=raw_inter

nC=ncol(raw_data)
nR=nrow(raw_data)
data_sort_class <- raw_data[order(raw_data$class),]
head(data_sort_class,5)
file_direc="PCa/"
file_type="Inter"
unique_class <- unique(raw_data$class)

for (i in 1:length(unique_class)){
  #Break down the data w.r.t class
  name_1 <- paste("df",file_type, unique_class[i], sep = "_")
  assign(name_1, raw_data[raw_data$class==unique_class[i],])  
  data_tmp <- eval(parse(text = name_1))
  p_check=replicate(dim(data_tmp)[1], 0)
  p_check[data_tmp$p_BH_eur<=0.05]=1
  p_check <- as.factor(p_check)
  data_tmp$pcheck <- p_check
  assign(name_1, data_tmp)  

  savedirec=paste(file_direc,"EURFreq","_",file_type,"_",unique_class[i],".eps",sep="")
  cols <- c("0"="black","1"="magenta3")
  shapes <-c("0"=1,"1"=16)
  
  sp_tmp1 <- ggplot(data_tmp, aes(x=eur_freq_met, y=eur_freq_local, group=pcheck))+geom_point(aes(shape=pcheck, color=pcheck),size=2, show.legend = FALSE)+scale_shape_manual(values=shapes)+scale_color_manual(values=cols)
  sp_tmp2 <- sp_tmp1 + geom_text(aes(label=ifelse(p_BH_eur <= 0.05,as.character(gene),'')), hjust=0,vjust=0, show.legend = FALSE,size=2) + geom_abline(intercept = 0, slope = 1) + xlim(0,0.5) + ylim(0,0.5)     
  sp_tmp3 <- sp_tmp2 + labs(title=paste(file_type,unique_class[i],sep="_"), x="EUR Mutation Freq from Met", y="EUR Mutation Freq from Local")
  sp <- sp_tmp3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black")) #+ scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  ggsave(file=savedirec,width=20,height=20,units="cm")
  #adjust label position by adding "position=position_jitter(width=0,height=0.001)" to geom_text
}

