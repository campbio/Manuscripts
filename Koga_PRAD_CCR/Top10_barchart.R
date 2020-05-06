setwd("/Users/sw/Desktop/HS/")
library(pwr)
library(tibble)
library(ggplot2)
library(tidyr)
library(fBasics)
library(rapportools)
library(statmod)

#### Read in either met or primary dataset (in csv)
file_direc="PCa/"
file_type="met"
filename=paste(file_direc,file_type,".csv",sep="")
raw_data <- read.csv(filename, header=TRUE)

nC=ncol(raw_data)
nR=nrow(raw_data)
data_sort_class <- raw_data[order(raw_data$class),]

raw_data<-raw_data[!(raw_data$class=='all'),]
raw_aasort <-raw_data[order(-raw_data$aa_freq),]
raw_eursort <-raw_data[order(-raw_data$eur_freq),]


list_fin=raw_aasort[1:20,]
cols <- c("point"="chartreuse4","truncation"="azure4","deletion"="mediumpurple","amplification"="red","rearrangement"="black","others"="yellow")
barc_tmp1<- ggplot()+geom_bar(data=list_fin,aes(y=aa_freq,x=reorder(paste(gene,class,sep="_"),-aa_freq),fill=class),stat="identity",width = 0.5)+theme_bw()
barc_tmp2<- barc_tmp1 + labs(title = "Top 20 Gene Alterations for AA in Local",x="Gene Alterations", y="Alteration Frequency") + scale_fill_manual(values = cols)
  barc_tmp2+guides(fill=guide_legend(title="Alteration Class")) +theme(axis.text.x=element_text(angle =45, vjust = 0.5))+ylim(0,0.3)

  barname=paste("PCa/Top20_AA_",file_type,".eps",sep="")
  
  ggsave(file=barname,width = 20,height = 20,units = "cm")

  list_fin=raw_eursort[1:20,]
  cols <- c("point"="chartreuse4","truncation"="azure4","deletion"="mediumpurple","amplification"="red","rearrangement"="black","others"="yellow")
  barc_tmp1<- ggplot()+geom_bar(data=list_fin,aes(y=eur_freq,x=reorder(paste(gene,class,sep="_"),-eur_freq),fill=class),stat="identity",width = 0.5)+theme_bw()
  barc_tmp2<- barc_tmp1 + labs(title = "Top 20 Gene Alterations for EUR in Local",x="Gene Alterations", y="Alteration Frequency") + scale_fill_manual(values = cols)
  barc_tmp2+guides(fill=guide_legend(title="Alteration Class")) +theme(axis.text.x=element_text(angle =45, vjust = 0.5))+ylim(0,0.4)
  
  barname=paste("PCa/Top20_EUR_",file_type,".eps",sep="")
  
  ggsave(file=barname,width = 20,height = 20,units = "cm")
  

  
  