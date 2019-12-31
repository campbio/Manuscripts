RNGkind(sample.kind = "Rounding")
library(devtools)
load_all("~/gitProjects/celda_irisapo", recompile=T)   # on fix branch, celda version 1.3.1
library(ggplot2)
library(gridExtra) 
source("functions.R") 


#qEmpty plots to fill spaces 
plt.empt = ggplot() + theme_void()  

# load stats data on original sequence data 
df.csv = read.csv("/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/analysis/gem_classification.csv", header=TRUE)


# load original data and decontaminated result
ori = readRDS("specificCellsWithLabelRownameEmsemble.rds")
dcon = readRDS("dcon.rds")

listPrePos  = colStat( dcon = dcon, ori = ori ) 



# Scatter plot of UMI counts  
plt.UMI= plotUMI( df.csv = df.csv[ df.csv$call != "Multiplet", ], color = c('slateblue', '#619CFF'), vline=TRUE  ) +  theme(legend.position = "right")  + xlab("hg19 UMI counts") + ylab("mm10 UMI counts")  + ggtitle("Original UMI counts\n(human-mouse mixture)")  + scale_x_continuous(breaks = c(0,20000,40000,60000,80000) , limits = c(0,90000) )  

p.UMI.hg = plotUMI( df.csv = df.csv[df.csv$call == "hg19" & df.csv$hg19 >= 0, ], color = c('slateblue') ) + theme(legend.position = "none") +  xlab("") + ylab("") +  scale_y_continuous(limits=c(0,1000), breaks = c(0, 500, 1000), labels = c(" 0", " 500", " 1000"), expand = c(0,0)) + scale_x_continuous(breaks =c(20000,40000,60000), expand=c(0,0)) 
p.UMI.mm = plotUMI( df.csv = df.csv[df.csv$call == "mm10" & df.csv$mm10 >= 0, ], color = c('#619CFF') ) + theme(legend.position = "none") +  xlab("")  + ylab("") +   scale_x_continuous(limits=c(0,1000), breaks=c(0,500,1000), expand =c(0,0))  + scale_y_continuous(breaks=c(10000,20000,30000),expand=c(0,0)) 


conp.tb = conpCal( summary.tb = df.csv ) 
print(table( conp.tb$call ))
cat("Using the original count matrix, the summary statistics of the proportion of exogenous transctrips for each droplet is: \n") ; summary(conp.tb$conp) 


p.UMIdcon = plotUMI(df.csv = df.csv[ df.csv$call != "Multiplet", ], color =  c('slateblue', '#619CFF'), vline=TRUE) + theme(legend.position = "none")  +  xlab("hg19 UMI counts") + ylab("mm10 UMI counts")  + ggtitle("Decontaminated UMI counts\n(human-mouse mixture)") + geom_point( data = listPrePos$decon.df, aes(x = hg19, y = mm10), color = mapvalues(x=listPrePos$decon.df$call, from=unique(listPrePos$decon.df$call),  to=c('darkorange1', '#00BA38') ) , size = 0.5) + scale_x_continuous(limits = c(0,90000), breaks = c(0,20000,40000,60000,80000))  
+ scale_x_continuous(expand=expand_scale(mult=c(0,0.1))) + scale_y_continuous(expand=expand_scale(mult = c(0, 0.1))) 
 

p.UMIdcon.hg = plotUMI( df.csv = df.csv[df.csv$call == "hg19" & df.csv$hg19 >= 0, ], color = c('slateblue'), legend.position = "right" ) + geom_point( data = listPrePos$decon.df[ listPrePos$decon.df$call == "hg19" & listPrePos$ori.df$hg19 >= 0,  ], aes( x=hg19, y=mm10, fill = "Decontaminated counts"), size = 0.5, color = c("#00BA38")   ) + guides( colour= guide_legend( override.aes = list( size = 2.5 ), title = ""), fill = guide_legend( override.aes = list( size = 2.5  ), title = "Human cell"   )  )    + scale_color_manual( labels = "Original counts", values = c('slateblue')  )  + xlab("") + ylab("") + ggtitle( "") +  scale_y_continuous(limits=c(0,1000), breaks = c(0, 500, 1000), labels = c(" 0", " 500", " 1000"), expand = c(0,0)) + scale_x_continuous(breaks =c(20000,40000,60000), expand=c(0,0) ) 

p.UMIdcon.mm = plotUMI( df.csv = df.csv[df.csv$call == "mm10" & df.csv$mm10 >= 0, ], color = c('#619CFF' ) , legend.position = "right") + geom_point( data = listPrePos$decon.df[ listPrePos$decon.df$call == "mm10" & listPrePos$ori.df$mm10 >= 0,  ], aes( x=hg19, y=mm10, fill = "Decontaminated counts"), size = 0.5, color = c("darkorange1") ) + guides( color = guide_legend( override.aes = list( size = 2.5) , title = ""), fill = guide_legend( override.aes = list( size = 2.5), title = "Mouse cell") ) + scale_color_manual( labels = "Original counts", values = c('#619CFF')  )  + xlab("") + ylab("")  + scale_x_continuous(limits=c(0,1000), breaks=c(0,500,1000), expand = c(0,0))  + scale_y_continuous(breaks=c(10000,20000,30000),expand=c(0,0) ) 


conp.tb.dcon = conpCal( summary.tb = listPrePos$decon.df )  
cat("after decontaminated by DecontX, the summary statistics of the proportion of exogenous transctrips for each droplet is: \n") ; summary(conp.tb.dcon$conp) 


# Violin plot of distribution of contamination level (proportion) 
plt.viol = plotViolin( conp.tb = conp.tb)  + xlab("Cell type") + theme( legend.position = "none") + ggtitle("Distribution of\ncontamination level") + ylab("Proportion of\nexogenous transcripts") 

call.list = c("hg19", "mm10")
for (call.type in call.list) {
	cell.type = ifelse(call.type == "hg19", "human", "mouse")
  cat("summary of proportions of exogenous transctrips in ", cell.type, "cells is \n")  
	print(summary( conp.tb$conp [ conp.tb$call == call.type ] ) ) 
}

# real contamination proportion vs estimated contamination proportion
df.conp = df.Conp( ori=ori, dcon=dcon)

conP.hg = plot.performance( df.conp, type="hg19" ) + scale_color_manual( values=c("slateblue" ) ) + theme( legend.position = "none") + ggtitle("Performance on\nhuman cells") + xlab("Proportion of mouse transcripts")  + ylab("Estimated contamination (%)")  

conP.mm = plot.performance( df.conp, type="mm10" ) + scale_color_manual( values=c("#619CFF") ) +  theme( legend.position = "none") + ggtitle("Performance on\nmouse cells") + xlab("Proportion of human transcripts") +  ylab("Estimated contamination (%)") 


for ( call.type in call.list ) {
	cell.type = ifelse(call.type == "hg19", "human", "mouse")
  cat("correlation of estimated contamination percentage and exogenous transcripts percentage in ", cell.type, " cells: \n" , cor( x = df.conp$tru.conp[ df.conp$call == call.type ], df.conp$est.conp [ df.conp$call == call.type ] )  )
	cat("\n")
  rmse(x = df.conp$tru.conp[ df.conp$call == call.type ], df.conp$est.conp [ df.conp$call == call.type ] )
  cat("RMSE is : \n", sqrt(mean( (df.conp$tru.conp[ df.conp$call == call.type ] - df.conp$est.conp[ df.conp$call == call.type ])^2 ) )   )  
	cat("\n")
}

cat("correlation of estimated contamination percentage and non-native transcripts percentage in each cell  \n" , cor( x = df.conp$tru.conp, df.conp$est.conp )  )
rmse(x = df.conp$tru.conp, df.conp$est.conp )
cat("RMSE is: \n", sqrt(mean( (df.conp$tru.conp - df.conp$est.conp ) ^2 ) )  ) 


# plot to show that mouse--specific gene expression in human cells is highly correlated to the distribution of the gene expression of an average mouse cell. 
counts.mouseG = listPrePos$counts.mouseGene
counts.humanG = listPrePos$counts.humanGene

mmGeneProp = geneProp( counts = counts.mouseG, z = ori$z.tb$call, curcellType = "mm10", otherType = "hg19" ) 
hgGeneProp = geneProp( counts = counts.humanG, z = ori$z.tb$call, curcellType = "hg19", otherType = "mm10" ) 

df.hm = rbind( mmGeneProp, hgGeneProp ) 


pltGene.mm = plot.geneExp( df.exp = mmGeneProp , color = "#619CFF")  + ylab("Gene expression of\nan average mouse cell") + xlab("Proportion of mouse-specific gene\ncounts in human cells") + ggtitle("Mouse-specific genes") #+ theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  

pltGene.hg = plot.geneExp( df.exp = hgGeneProp, color='slateblue') + ylab("Gene expression of\nan average human cell") + xlab("Proportion of human-specific gene\ncounts in mouse cells") + ggtitle("Human-specific genes") #+ theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  

pltGene.mm = plot.geneExp( df.exp = mmGeneProp , color = "black")  + ylab("Gene expression of\nan average mouse cell") + xlab("Proportion of mouse-specific gene\ncounts in human cells") + ggtitle("Mouse-specific genes") #+ theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  
pltGene.hg = plot.geneExp( df.exp = hgGeneProp, color='black') + ylab("Gene expression of\nan average human cell") + xlab("Proportion of human-specific gene\ncounts in mouse cells") + ggtitle("Human-specific genes") #+ theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  

for (call.type in call.list) {
	df.type = ifelse(call.type == "hg19", "hgGeneProp", "mmGeneProp") 
	cell.type = ifelse(call.type == "hg19", "human", "mouse")
  cat( "Correlation between the correlation of \n", 
			 "               1) proportions of mouse-specific genes in exogenous transcripts cross all ", cell.type, " cells \n", 
			 "           and 2) the distribution of those genes' expression in an average ", cell.type, "cell is: \n ", 
			 "               ", cor( get(df.type)$real.averageExp.cpm, get(df.type)$prop.inContaminatedTranscripts) )    
	cat("\n")
}


gb.p.UMI = ggplotGrob( plt.UMI + theme( panel.grid.major = element_blank(), legend.position="none")  )  
gb.p.UMI.hg = ggplotGrob( p.UMI.hg+ theme(plot.margin=unit(rep(0,4),"null") ))  
gb.p.UMI.mm = ggplotGrob( p.UMI.mm+ theme(plot.margin=unit(rep(0,4),"null") ))  

gtb = gtable(widths = unit(rep(1,12), "null"), height = unit(rep(1, 12) , "null"))  
gtb = gtable_add_grob(gtb,  gb.p.UMI, t=1, b=12,l=1, r=12) 
gtb = gtable_add_grob(gtb, gb.p.UMI.hg, t=8, b=10, l=7, r=11) 
gtb = gtable_add_grob(gtb, gb.p.UMI.mm, t=3, b=7, l=3, r=5)
#grid.arrange( gtb) 

pg.pltGene = arrangeGrob( pltGene.mm, pltGene.hg, ncol = 1 ) 

g.UMI = ggplot_gtable(ggplot_build(plt.UMI))
glegend.UMI = g.UMI$grob[[which(sapply(g.UMI$grobs, function(x) x$name) == "guide-box")]]
pg.UMIz = arrangeGrob( gtb,  glegend.UMI, ncol = 2, widths = c(2, 0.5)) 
pg.GnV = arrangeGrob( pltGene.mm + ggtitle(""), pltGene.hg + ggtitle(""), plt.viol, ncol = 3 , widths = c(4, 4, 3) ) 
pdf("Fig2.pdf", width = 8.5, height = 8.5)
grid.arrange( pg.UMIz, pg.GnV, ncol = 1, heights = c(3.2, 2) , top = "Figure 2") 
dev.off()



gb.p.UMIdcon = ggplotGrob(p.UMIdcon + theme(panel.grid.major = element_blank(), legend.position="none"))
gb.p.UMIdcon.hg = ggplotGrob( p.UMIdcon.hg + theme(legend.position="none", plot.margin=unit(rep(0,4),"null") ) )
gb.p.UMIdcon.mm = ggplotGrob( p.UMIdcon.mm + theme(legend.position="none", plot.margin=unit(rep(0,4),"null") ) )

gtbdcon = gtable(widths = unit(rep(1,12), "null"), height = unit(rep(1, 13) , "null")) 
gtbdcon = gtable_add_grob(gtbdcon, gb.p.UMIdcon, t=1, b=13,l=1, r=12) 
gtbdcon = gtable_add_grob(gtbdcon, gb.p.UMIdcon.hg, t=7, b=10, l=7, r=11)
gtbdcon = gtable_add_grob(gtbdcon, gb.p.UMIdcon.mm, t=3, b=7, l=3, r=5)

g.UMIdcon.hg = ggplot_gtable(ggplot_build(p.UMIdcon.hg))
glegend.UMIdcon.hg = g.UMIdcon.hg$grob[[which(sapply(g.UMIdcon.hg$grobs, function(x) x$name) == "guide-box")]] 
g.UMIdcon.mm = ggplot_gtable(ggplot_build(p.UMIdcon.mm))
glegend.UMIdcon.mm = g.UMIdcon.mm$grob[[which(sapply(g.UMIdcon.mm$grobs, function(x) x$name) == "guide-box")]] 

pg.estP = arrangeGrob(conP.hg, conP.mm, plt.empt,  ncol = 3, widths = c(2, 2, 1) ) 
pg.legend = arrangeGrob( glegend.UMIdcon.hg, glegend.UMIdcon.mm, ncol = 1 ) 
pg.UMIdcon.z = arrangeGrob( gtbdcon,  pg.legend,
			   ncol = 2, widths = c(4,  1.2)  ) 
pdf("Fig3.pdf", width = 8.5, height = 8.5)
grid.arrange( pg.UMIdcon.z, pg.estP, ncol =1, heights = c(3.2, 2),    top = "Figure 3") 
dev.off()



# Supplementary scatter plot:  the estimated gene-level contamination distribution correlated with the expression of an average huamn/mouse cell
mmGeneProp = attach.estContaminatedDist(df.geneProp =  mmGeneProp, z = ori$z$call,  dcon = dcon, curcellType = "mm10") 
hgGeneProp = attach.estContaminatedDist( df.geneProp = hgGeneProp, z = ori$z$call, dcon = dcon, curcellType = "hg19" )

plt.cpr.mm = plot.estGeneExp( df.geneProp=mmGeneProp, color = "#619CFF" ) + ylab("Gene expression of an average mouse cell") + xlab("DecontX contamination distribution\nin human cells") + ggtitle("Mouse-specific genes")  + theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  
plt.cpr.mm = plot.estGeneExp( df.geneProp=mmGeneProp, color = "black" ) + ylab("Gene expression of an average mouse cell") + xlab("DecontX contamination distribution\nin human cells") + ggtitle("Mouse-specific genes")  + theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  


plt.cpr.hg = plot.estGeneExp( df.geneProp = hgGeneProp, color = "slateblue" )+ ylab("Gene expression of an average human cell") + xlab("DecontX contamination distribution\nin mouse cells") + ggtitle("Human-specific genes") + theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  
plt.cpr.hg = plot.estGeneExp( df.geneProp = hgGeneProp, color = "black" )+ ylab("Gene expression of an average human cell") + xlab("DecontX contamination distribution\nin mouse cells") + ggtitle("Human-specific genes") + theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  

for (call.type in call.list) {
	cell.type = ifelse(call.type == "hg19", "human", "mouse")
  tb.realExp = get( ifelse(call.type == "hg19", "hgGeneProp", "mmGeneProp") ) 
  cat( "correlation of 1) proportions of ",  cell.type, " genes' expression in an average ", cell.type, "cells \n", 
			 "           and 2) probability of ", cell.type, " genes' expression in estimated ", cell.type, " expression distribution is \n")
			print( cor( tb.realExp$real.averageExp.cpm, tb.realExp$estContaminationDist.crossSameGenome) )  
}

pdf("Sup_Fig_2.pdf", width = 8, height = 4)
grid.arrange( plt.cpr.mm, plt.cpr.hg, ncol = 2, top = "Supplementary Figure 2") 
dev.off()


plt.UMIhg = plotUMI( df.csv = df.csv[ df.csv$call == "hg19", ], color = "slateblue")  + xlab("") + ylab("") + ggtitle("   ")  
plt.UMImm = plotUMI( df.csv = df.csv[ df.csv$call == "mm10", ], color = '#619CFF')  + xlab("") + ylab("") + ggtitle("  \n  ")

plt.UMIall= plotUMI( df.csv = df.csv, color = c('slateblue','#619CFF', 'grey'), legend.position = "right", label = list("hg19"="Human", "mm10"="Mouse", "Multiplet" = "Multiplet")  ) + xlab("hg19 UMI counts") + ylab("mm10 UMI counts")  + ggtitle("Original UMI counts (human-mouse mixture)")     
pdf("Sup_Fig_1.pdf", width = 5, height = 5) 
grid.arrange( plt.UMIall, top = "Supplementary Figure 1") 
dev.off()
