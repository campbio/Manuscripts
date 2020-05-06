
#### HS 2020
####  Scatterplot of Foundation Cohort and gene set bar chart analysis
library(pwr)
library(tibble)
library(ggplot2)
library(tidyr)
library(fBasics)
library(rapportools)
library(statmod)


#### import foundation dataset, combined, local and met csv. 
file_direc="PCa/"
file_type="Combined"
filename=paste(file_direc,file_type,".csv",sep="")
raw_data <- read.csv(filename, header=TRUE)
raw_data1 <- read.csv("PCa/local.csv", header=TRUE)
raw_data2 <-read.csv("PCa/met.csv", header=TRUE)
raw_data3=raw_data2



raw_data1$gene=gene
raw_data2$gene=gene
raw_data_3$gene=gene
m=0
#power for proportion
for (i in 1:nrow(raw_data)) {
print(power.fisher.test(raw_data$aa_freq[i], raw_data$eur_freq[i], raw_data$aa_denom[i], raw_data$eur_denom[i], alpha=0.05, alternative="two.sided"))
}
#power based on sample size to acquire a 0.05 level
for (i in 1:nrow(raw_data)) {
  print(pwr.2p2n.test(h=0.5,n1=raw_data$aa_denom[i],n2=raw_data$eur_denom[i],sig.level = 0.05)$power)
}

#Below section is to isolate from all alteration class to exclude deletion, amplification and rearrangement. 
#This way everything left would be "true" mutation. 
for (i in 1:nrow(raw_data)) {
  if(raw_data$class[i]=="all") {
    raw_data_temp=raw_data[raw_data$gene==raw_data$gene[i],]
    raw_data$aa_freq[i]=raw_data$aa_freq[i]-sum(raw_data_temp[raw_data_temp$class %in% c("rearrangement","deletion","amplification"),4])
    raw_data$eur_freq[i]=raw_data$eur_freq[i]-sum(raw_data_temp[raw_data_temp$class %in% c("rearrangement","deletion","amplification"),3])
    raw_data$aa_count[i]=raw_data$aa_count[i]-sum(raw_data_temp[raw_data_temp$class %in% c("rearrangement","deletion","amplification"),7])
    raw_data$eur_count[i]=raw_data$eur_count[i]-sum(raw_data_temp[raw_data_temp$class %in% c("rearrangement","deletion","amplification"),8])
    
  }}

raw_data_mut=raw_data[raw_data$class=="all" & raw_data$eur_freq+raw_data$aa_freq!=0,]
raw_data_mut$class="mut"
for (i in 1:nrow(raw_data_mut)) {
  
list_prefish <- data.frame(c(raw_data_mut$aa_count[i],raw_data_mut$aa_denom[i]-raw_data_mut$aa_count[i]),
                           c(raw_data_mut$eur_count[i],raw_data_mut$eur_denom[i]-raw_data_mut$eur_count[i]))
prefish_mx <-as.matrix(sapply(list_prefish, as.numeric))  
raw_data_mut[i,9]=fisher.test(prefish_mx,alternative="two.sided")[1]
raw_data_mut[i,11]=power.fisher.test(raw_data_mut$aa_freq[i], raw_data_mut$eur_freq[i], raw_data_mut$aa_denom[i], raw_data_mut$eur_denom[i], alpha=0.05, alternative="two.sided")

}
pvalues=raw_data_mut$p
raw_data_mut$p_BH=p.adjust(pvalues,method="BH")

write.csv(raw_data_mut, file = "PCa/Local_mutation.csv",row.names=FALSE, na="")

#Write mutation frequency (all mutation excluding the insertation/deletion/fusion/amplification) to csv file. 

for (i in 1:nrow(raw_data2)) {
  raw_data2_match=raw_data1[raw_data1$gene==lapply(raw_data2[i,1],as.character) & raw_data1$class==raw_data2[i,2],]
  if(!is.empty(raw_data2_match)[1]) {
    raw_data2[i,5]=raw_data2_match$aa_denom+raw_data2[i,5]
    raw_data2[i,6]=raw_data2_match$eur_denom+raw_data2[i,6]
    raw_data2[i,3]=(raw_data2[i,8]+raw_data2_match$eur_count)/raw_data2[i,6]
    raw_data2[i,4]=(raw_data2[i,7]+raw_data2_match$aa_count)/raw_data2[i,5]
    raw_data2[i,7]=raw_data2[i,7]+raw_data2_match$aa_count
    raw_data2[i,8]=raw_data2[i,8]+raw_data2_match$eur_count
    list_prefish <- data.frame(c(raw_data2[i,7],raw_data2[i,5]-raw_data2[i,7]),
                               c(raw_data2[i,8],raw_data2[i,6]-raw_data2[i,8]))
    prefish_mx <-as.matrix(sapply(list_prefish, as.numeric))  
    raw_data2[i,9]=fisher.test(prefish_mx,alternative="two.sided")[1]
    
    m=m+1
  }
}
pvalues=raw_data2$p
raw_data2$p_BH=p.adjust(pvalues,method="BH")
write.csv(raw_data2, file = "PCa/Combined.csv",row.names=FALSE, na="")
raw_data<-raw_data2

nC=ncol(raw_data)
nR=nrow(raw_data)
data_sort_class <- raw_data[order(raw_data$class),]
head(data_sort_class,5)
#df3 <- rbind(df1, df2)
#data_sort_p <-raw_data[order(raw_data$p),-1]
#head(data_sort_p,5)
#Extract how many classes are there
unique_class <- unique(raw_data$class)

for (i in 1:length(unique_class)){
#Break down the data w.r.t class
name_1 <- paste("df",file_type, unique_class[i], sep = "_")
assign(name_1, raw_data[raw_data$class==unique_class[i],])  
data_tmp <- eval(parse(text = name_1))
p_check=replicate(dim(data_tmp)[1], 0)
p_check[data_tmp$p_BH<=0.05]=1
p_check <- as.factor(p_check)
data_tmp$pcheck <- p_check
assign(name_1, data_tmp)  
#if ((length(data1_tmp)!=1)&&(length(data2_tmp)!=1)) {
#  M1=mean(data1_tmp) 
#  M2=mean(data2_tmp)
#  S1=sd(data1_tmp)
#  S2=sd(data2_tmp)
#  DD = (M1 - M2)/sqrt(((S1^2) + (S2^2))/2)  
#  ans=pwr.t.test(d=DD,power=.8,sig.level=.05,type="two.sample",alternative="two.sided")
#  print(unique_class[i])
#  print(ans$n) }
#plot(data1_tmp, data2_tmp, main=unique_class[i], xlab=name_1, ylab=name_2, pch=19)

savedirec=paste(file_direc,"Freq","_",file_type,"_",unique_class[i],".eps",sep="")
cols <- c("0"="black","1"="magenta3")
shapes <-c("0"=1,"1"=16)

sp_tmp1 <- ggplot(data_tmp, aes(x=eur_freq, y=aa_freq, group=pcheck))+geom_point(aes(shape=pcheck, color=pcheck),size=2, show.legend = FALSE)+scale_shape_manual(values=shapes)+scale_color_manual(values=cols)
sp_tmp2 <- sp_tmp1 + geom_text(aes(label=ifelse(p_BH <= 0.05,as.character(gene),'')), hjust=0,vjust=0, show.legend = FALSE,size=2) + geom_abline(intercept = 0, slope = 1) + xlim(0,0.5) + ylim(0,0.5)     
sp_tmp3 <- sp_tmp2 + labs(title=paste(file_type,unique_class[i],sep="_"), x="EUR Gene Mutation Frequency", y="AFR Gene Mutation Frequency")
sp <- sp_tmp3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) #+ scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
#ggsave(file=savedirec,width=20,height=20,units="cm")
#adjust label position by adding "position=position_jitter(width=0,height=0.001)" to geom_text
}

#p.adjust(p, method = BH, n = length(p))
raw_data5 <- raw_data[raw_data$aa_count>5 & raw_data$eur_count>5,]
p=raw_data5$p
raw_data5$p_BH=p.adjust(p, method = 'BH', n = length(p))
write.csv(raw_data5, file = "PCa/Combined_greaterthan5.csv",row.names=FALSE, na="")


# This is an example of a single fisher's test
list_prefish <- data.frame(c(37,53-37),
                           c(87,123-87))
prefish_mx <-as.matrix(sapply(list_prefish, as.numeric))  
fisher.test(prefish_mx,alternative="two.sided")
power.fisher.test(0.06,0.21, 50, 50, alpha=0.05, alternative="two.sided")
power.prop.test(n=50,p1=0.06,p2=0.21,sig.level=0.1)

