##### Gene set analysis

#### HS 2020

library(pwr)
library(tibble)
library(ggplot2)

#Enter the working directory
file_direc="PCa/"
#Enter the file name (local or met or combined)
file_type="local"
filename=paste(file_direc,file_type,".csv",sep="")
raw_data <- read.csv(filename, header=TRUE)

nC=ncol(raw_data)
nR=nrow(raw_data)

#Extract how many classes are there
unique_class <- unique(raw_data$class)


#2015 Robinson et al. + 2019 Abida et al.,, 
# Point mutation: SPOP, FOXA1, TP53
#copy number alterations: MYC, RB1, PTEN, CHD1
c_repairs=c('BRCA1','BRCA2','CHEK1','CHEK2','BAP1','BARD1','BRIP1','BACH1','FAM175A','MRE11A','RAD50','NBN','RAD51','RAD51B','XRCC2','XRCC3','RPA1','BLM','GEN1','FANCA','FANCC','FANCD2','FANCE','FANCF','FANCG','FANCI','FANCL','FANCM','FANCQ','ATM','CDK12','MLH1','MSH2','ATR','PALB2')
c_ChromatinModifier=c('KMT2C','KMT2D','KDM6A','CHD1')
c_CellCycle=c('RB1','CDKN1B','CDKN2A','CCND1','CDK4')
c_epigenetic=c('AR1D1A','ARID2','ARID4A','KDM6A','KDM6A','KMT2A','KMT2C','KMT2D','MBD1','SETD2','SETDB1','SMARCA1','SMARCAD1')
c_PI3K=c('AKT1','MTOR','PIK3CA','PIK3CB','PIK3R1','PIK3R2','PIK3R3','PTEN','TSC1','TSC2')
c_RASMAPK=c('BRAF','HRAS','KRAS','MAP2K1','MAP3K1','NF1','RASA1')
c_WNT=c('APC','AXIN1','AXIN2','CTNNB1','MF43','RNF43')
#Lysine methyltransferases
c_KMT=c('EHMT2','EHMT1','SETDB1','SETDB2','MLL','KMT2B','MLL3','MLL2','SETD2','WHSC1L1','KMT2C')
c_PIDMYC=c('AX1N1','HBP1','SUPT3H','SUPT7L','ACTL6A','MYC','PML','GSK3B','MAX','PPP2CA','TAF10','ZBTB17','PAK2','SKP2','PIN1','PPP2R5A','TAF12','TAF9','CDKN2A','KAT2A','KAT5','FBXW7','RUVBL2','RUVBL1','TRRAP')


list_repairs=raw_data[which(raw_data$gene %in% c_repairs),]  
list_ChromatinModifier=raw_data[which(raw_data$gene %in% c_ChromatinModifier),]  
list_CellCycle=raw_data[which(raw_data$gene %in% c_CellCycle),]  
list_epigenetic=raw_data[which(raw_data$gene %in% c_epigenetic),]  
list_PI3K=raw_data[which(raw_data$gene %in% c_PI3K),]  
list_RASMAPK=raw_data[which(raw_data$gene %in% c_RASMAPK),] 
list_WNT=raw_data[which(raw_data$gene %in% c_WNT),]  
list_KMT=raw_data[which(raw_data$gene %in% c_KMT),]
list_PIDMYC=raw_data[which(raw_data$gene %in% c_PIDMYC),]

path_list=c('repairs','ChromatinModifier','CellCycle','epigenetic','PI3K','RASMAPK','WNT','KMT','PIDMYC')

ylim_list=c(0.12,0.1,0.12,0.085,0.4,0.04,0.175,0.15,0.3)
for (j in 1:length(path_list)) {
  
  name_listtmp <- paste("list",path_list[j], sep = "_")
  assign(name_1, raw_data[raw_data$class==unique_class[i],])  
  list_tmp <- eval(parse(text = name_listtmp))
  list_repair=list_tmp[list_tmp$class=="all",]   
  list_gene=as.character(factor(list_repair$gene))
  #shapiro.test(list_repair$aa_freq)
  #shapiro.test(list_repair$eur_freq)
  
  
  for (i in 1:length(unique(list_tmp$class))) {
    
    list_byclass=list_tmp[list_tmp$class==unique(list_tmp$class)[i],]
    aa_byclass=list_byclass$aa_freq
    eur_byclass=list_byclass$eur_freq
    if(length(aa_byclass)>1) {
      print(unique(list_tmp$class)[i])
      print(t.test(aa_byclass,eur_byclass,paired=TRUE))
    }
  }
  
  for (i in 1:length(unique(list_gene))) {
    aa_sum_class=sum(list_tmp[list_tmp$gene==unique(list_gene)[i] & list_tmp$class!="all",4])
    eur_sum_class=sum(list_tmp[list_tmp$gene==unique(list_gene)[i] & list_tmp$class!="all",3])
    list_tmp[list_tmp$gene==unique(list_gene)[i] & list_tmp$class=="all",4]=max(c(list_tmp[list_tmp$gene==unique(list_gene)[i] & list_tmp$class=="all",4]-aa_sum_class,0))
    list_tmp[list_tmp$gene==unique(list_gene)[i] & list_tmp$class=="all",3]=max(c(list_tmp[list_tmp$gene==unique(list_gene)[i] & list_tmp$class=="all",3]-eur_sum_class,0))
    
  }
  
  
  sub_gene <-rep(array(list_tmp$gene),2)
  sub_p <-rep(list_tmp$p_BH,2)
  sub_aa <-c(list_tmp$aa_freq)
  sub_eur <-c(list_tmp$eur_freq)
  freq <- c(sub_aa,sub_eur)
  sub_class <-rep(array(list_tmp$class),2)
  sub_class[sub_class=="all"]="others"
  le=length(sub_aa)
  type <- c(rep("AA",le),rep("EUR",le))
  sub_gene[sub_gene==unique(sub_gene[sub_p<=0.05 & sub_class=="others"])]=paste(sub_gene[sub_gene==unique(sub_gene[sub_p<=0.05 & sub_class=="others"])],"*",sep="")
  list_fin=data.frame(sub_gene,freq,sub_class,type,sub_p)
  
  cols <- c("point"="green","truncation"="black","deletion"="blue","amplification"="red","reaarangement"="maganta4","others"="grey")
  barc_tmp1<- ggplot()+geom_bar(data=list_fin,aes(y=freq,x=type,fill=sub_class),stat="identity",position="stack",width = 0.5)+theme_bw()+facet_grid(~sub_gene)
  barc_tmp2<- barc_tmp1 + labs(title = paste(path_list[j],"Genes on Local",sep=" "),
                               x="Patient Groups", y="Mutation Frequency") + scale_fill_manual(values = cols)
  barc_tmp2+guides(fill=guide_legend(title="Mutation Class")) +ylim(0,ylim_list[j])
  
  barname=paste("PCa/Bar_",file_type,"_",path_list[j],".eps",sep="")
  ggsave(file=barname,width = 20,height = 20,units = "cm")
  
  
}
