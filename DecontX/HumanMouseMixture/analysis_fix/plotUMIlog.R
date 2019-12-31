library(devtools)
load_all(path="~/gitProjects/celda", recompile=F)
library(ggplot2)
library(gridExtra) 
source("functions.R") 
library(gtable)

# Empty plots to fill spaces 
plt.empt = ggplot() + theme_void()  

# load stats data on original sequence data 
df.csv = read.csv("/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/analysis/gem_classification.csv", header=TRUE)


# load original data and decontaminated result
ori = readRDS("specificCellsWithLabelRownameEmsemble.rds")
dcon = readRDS("dcon.rds")

# transfer count to log2-scale
ori.log2 = ori
dcon.log2 = dcon
ori.log2$count = log2( ori$count + 1) 
dcon.log2$resList$estNativeCounts = log2( dcon$resList$estNativeCounts + 1) 

listPrePos.log2  = colStat( dcon = dcon.log2, ori = ori.log2 ) 
listPrePos = colStat( dcon = dcon, ori = ori) 


# reculculate df.csv on log2-scale 


# Scatter plot of UMI counts  
plt.UMI= plotUMI( df.csv = df.csv[ df.csv$call != "Multiplet", ], color = c('slateblue', '#619CFF')  ) +  theme(legend.position = "right")  + xlab("hg19 UMI counts") + ylab("mm10 UMI counts")  + ggtitle("Original UMI counts\n(human-mouse mixture)")  

p.UMI.log2= plotUMI( df.csv = listPrePos.log2$ori.df, color = c('slateblue', '#619CFF')  ) +  theme(legend.position = "right")  + xlab("hg19 UMI counts (log2)") + ylab("mm10 UMI counts (log2)")  + ggtitle("Original UMI counts\n(human-mouse mixture)")  



p.UMI.hg.log2 = plotUMI( df.csv = listPrePos.log2$ori.df[listPrePos.log2$ori.df$call == 'hg19' & listPrePos.log2$ori.df$hg19 > 5000,  ] , color = c('slateblue') ) + theme(legend.position = "none") +  xlab("") + ylab("") +  scale_y_continuous( breaks = c(0, 500, 1000), labels = c("  0", "  500", "  1000") ) 
p.UMI.mm.log2 = plotUMI( df.csv = listPrePos.log2$ori.df[listPrePos.log2$ori.df$call == "mm10" & listPrePos.log2$ori.df$mm10 > 5000,  ], color = c('#619CFF') ) + theme(legend.position = "none") +  xlab("")  + ylab("")  + scale_x_continuous( breaks = c(0,1000,2000) ) + scale_y_continuous( breaks = c(0 ,6000, 8000,10000) )  



p.UMIdcon.log2 = plotUMI(df.csv = listPrePos.log2$ori.df, color =  c('slateblue', '#619CFF')) + theme(legend.position = "none")  +  xlab("hg19 UMI counts (log2)") + ylab("mm10 UMI counts (log2)")  + ggtitle("Decontaminated UMI counts\n(human-mouse mixture)") + geom_point( data = listPrePos.log2$decon.df, aes(x = hg19, y = mm10), color = mapvalues(x=listPrePos$decon.df$call, from=unique(listPrePos$decon.df$call),  to=c('darkorange1', '#00BA38') ) , size = 0.5)  

p.UMIdcon.hg.log2 = plotUMI( df.csv = listPrePos.log2$ori.df[listPrePos.log2$ori.df$call == 'hg19' & listPrePos.log2$ori.df$hg19 > 5000,  ]   ,color = c('slateblue'), legend.position = "right" ) + geom_point( data = listPrePos.log2$decon.df[ listPrePos.log2$decon.df$call == "hg19" & listPrePos.log2$ori.df$hg19 >5000 ,  ], aes( x=hg19, y=mm10, fill = "Decontaminated"), size = 0.5, color = c("#00BA38")   ) + guides( colour= guide_legend( override.aes = list( size = 2.5 ), title = "hg19 call"), fill = guide_legend( override.aes = list( size = 2.5  ), title = ""   )  )    + scale_color_manual( labels = "Original", values = c('slateblue')  )  + xlab("") + ylab("") + ggtitle( "") + scale_y_continuous( breaks = c(0, 500, 1000) )   

p.UMIdcon.mm.log2 = plotUMI( df.csv = listPrePos.log2$ori.df[listPrePos.log2$ori.df$call == 'mm10' & listPrePos.log2$ori.df$mm10 > 5000, ], color = c('#619CFF' ) , legend.position = "right") + geom_point( data = listPrePos.log2$decon.df[ listPrePos.log2$decon.df$call == "mm10" & listPrePos.log2$ori.df$mm10 >5000,  ], aes( x=hg19, y=mm10, fill = "Decontaminated"), size = 0.5, color = c("darkorange1") ) + guides( color = guide_legend( override.aes = list( size = 2.5) , title = "mm10 call"), fill = guide_legend( override.aes = list( size = 2.5), title = "") ) + scale_color_manual( labels = "Original", values = c('#619CFF')  )  + xlab("") + ylab("") + ggtitle("") + scale_x_continuous( breaks = c(0,1000,2000) ) + scale_y_continuous( breaks = c(0 ,6000, 8000,10000) )  


# real contamination proportion vs estimated contamination proportion
df.conp = df.Conp( ori=ori, dcon=dcon)

conP.hg = plot.performance( df.conp, type="hg19" ) + scale_color_manual( values=c("slateblue" ) ) + theme( legend.position = "none") + ggtitle("Performance on\nhuman cells") + xlab("Proportion of mouse transcripts")  + ylab("Estimated contamination (%)")  

conP.mm = plot.performance( df.conp, type="mm10" ) + scale_color_manual( values=c("#619CFF") ) +  theme( legend.position = "none") + ggtitle("Performance on\nmouse cells") + xlab("Proportion of human transcripts") +  ylab("Estimated contamination (%)") 


cat("correlation of estimated contamination percentage and non-native transcripts percentage in human cells: \n" , cor( x = df.conp$tru.conp[ df.conp$call == "hg19" ], df.conp$est.conp [ df.conp$call == "hg19" ] )  )
rmse(x = df.conp$tru.conp[ df.conp$call == "hg19" ], df.conp$est.conp [ df.conp$call == "hg19" ] )
cat("correlation of estimated contamination percentage and non-native transcripts percentage in mouse cells  \n" , cor( x = df.conp$tru.conp[ df.conp$call == "mm10" ], df.conp$est.conp [ df.conp$call == "mm10" ] )  )
rmse(x = df.conp$tru.conp[ df.conp$call == "mm10" ], df.conp$est.conp [ df.conp$call == "mm10" ] )


# Violin plot of distribution of contamination level (proportion) 
conp.tb = conpCal( summary.tb = df.csv ) 

plt.viol = plotViolin( conp.tb = conp.tb)  + xlab("Cell type") + theme( legend.position = "none") + ggtitle("Distribution of\ncontamination level") + ylab("Proportion of\nexogenous transcripts") 
theme(legend.position=c(1.15, 0.6), plot.margin=unit(c(1, 3,1, 1) , "lines"), legend.background = element_rect( fill = NA) )


cat("summary of contamination proportion of human cells \n");  summary( conp.tb$conp [ conp.tb$call == "hg19" ] )  
cat("summary of contamination proportion of mouse cells \n"); summary( conp.tb$conp [ conp.tb$call == "mm10" ] )  


# plot to show that mouse--specific gene expression in human cells is highly correlated to the distribution of the gene expression of an average mouse cell. 
counts.mouseG = listPrePos$counts.mouseGene
counts.humanG = listPrePos$counts.humanGene

mmGeneProp = geneProp( counts = counts.mouseG, z = ori$z.tb$call, curcellType = "mm10", otherType = "hg19" ) 
hgGeneProp = geneProp( counts = counts.humanG, z = ori$z.tb$call, curcellType = "hg19", otherType = "mm10" ) 

df.hm = rbind( mmGeneProp, hgGeneProp ) 


pltGene.mm = plot.geneExp( df.exp = mmGeneProp , color = "#619CFF")  + ylab("Gene expression of\nan average mouse cell") + xlab("Proportion of mouse-specific gene\ncounts in human cells") + ggtitle("Mouse-specific genes") #+ theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  

pltGene.hg = plot.geneExp( df.exp = hgGeneProp, color='slateblue') + ylab("Gene expression of\nan average human cell") + xlab("Proportion of human-specific gene\ncounts in mouse cells") + ggtitle("Human-specific genes") #+ theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  


cat( "correlation is:", cor( mmGeneProp$real.averageExp.cpm, mmGeneProp$prop.inContaminatedTranscripts)  ) 
cat( "correlation is:", cor( hgGeneProp$real.averageExp.cpm, hgGeneProp$prop.inContaminatedTranscripts) )  

gb.p.UMI.log2 = ggplotGrob( p.UMI.log2 + theme( panel.grid.major = element_blank(), legend.position="none")  )  
gb.p.UMI.hg.log2 = ggplotGrob( p.UMI.hg.log2 + theme(plot.margin=unit(rep(0,4),"null") ))  
gb.p.UMI.mm.log2 = ggplotGrob( p.UMI.mm.log2 + theme(plot.margin=unit(rep(0,4),"null") ))  


gtb = gtable(widths = unit(rep(1,12), "null"), height = unit(rep(1, 12) , "null"))  
gtb = gtable_add_grob(gtb,  gb.p.UMI.log2, t=1, b=12,l=1, r=12) 
gtb = gtable_add_grob(gtb, gb.p.UMI.hg.log2, t=7, b=9, l=7, r=11) 
gtb = gtable_add_grob(gtb, gb.p.UMI.mm.log2, t=3, b=7, l=4, r=6)
#grid.arrange( gtb) 


pg.pltGene = arrangeGrob( pltGene.mm, pltGene.hg, ncol = 1 ) 


g.UMI = ggplot_gtable(ggplot_build(p.UMI.log2))
glegend.UMI = g.UMI$grob[[which(sapply(g.UMI$grobs, function(x) x$name) == "guide-box")]]
pg.UMIz = arrangeGrob( gtb,  glegend.UMI, ncol = 2, widths = c(2, 0.5)) 
pg.GnV = arrangeGrob( pltGene.mm, pltGene.hg, plt.viol, ncol = 3 , widths = c(4, 4, 3) ) 
pdf("Fig2_log2.pdf", width = 8.5, height = 8.5)
grid.arrange( pg.UMIz, pg.GnV, ncol = 1, heights = c(3.2, 2) , top = "Figure 2") 
dev.off()


gb.p.UMIdcon.log2 = ggplotGrob(p.UMIdcon.log2 + theme(panel.grid.major = element_blank(), legend.position="none"))
gb.p.UMIdcon.hg.log2 = ggplotGrob( p.UMIdcon.hg.log2 + theme(legend.position="none", plot.margin=unit(rep(0,4),"null") ) )
gb.p.UMIdcon.mm.log2 = ggplotGrob( p.UMIdcon.mm.log2 + theme(legend.position="none", plot.margin=unit(rep(0,4),"null") ) )

gtbdcon = gtable(widths = unit(rep(1,12), "null"), height = unit(rep(1, 13) , "null")) 
gtbdcon = gtable_add_grob(gtbdcon, gb.p.UMIdcon.log2, t=1, b=13,l=1, r=12) 
gtbdcon = gtable_add_grob(gtbdcon, gb.p.UMIdcon.hg.log2, t=7, b=10, l=7, r=11)
gtbdcon = gtable_add_grob(gtbdcon, gb.p.UMIdcon.mm.log2, t=3, b=7, l=4, r=6)

g.UMIdcon.hg = ggplot_gtable(ggplot_build(p.UMIdcon.hg.log2))
glegend.UMIdcon.hg = g.UMIdcon.hg$grob[[which(sapply(g.UMIdcon.hg$grobs, function(x) x$name) == "guide-box")]] 
g.UMIdcon.mm = ggplot_gtable(ggplot_build(p.UMIdcon.mm.log2))
glegend.UMIdcon.mm = g.UMIdcon.mm$grob[[which(sapply(g.UMIdcon.mm$grobs, function(x) x$name) == "guide-box")]] 

pg.estP = arrangeGrob(conP.hg, conP.mm, plt.empt,  ncol = 3, widths = c(2, 2, 1) ) 
pg.legend = arrangeGrob( glegend.UMIdcon.hg, glegend.UMIdcon.mm, ncol = 1 ) 
pg.UMIdcon.z = arrangeGrob( gtbdcon,  pg.legend,
			   ncol = 2, widths = c(2,  1)  ) 
pdf("Fig3_log2.pdf", width = 8.5, height = 8.5)
grid.arrange( pg.UMIdcon.z, pg.estP, ncol =1, heights = c(3.2, 2),    top = "Figure 3") 
dev.off()


# Supplementary scatter plot:  the estimated gene-level contamination distribution correlated with the expression of an average huamn/mouse cell
mmGeneProp = attach.estContaminatedDist(df.geneProp =  mmGeneProp, z = ori$z$call,  dcon = dcon, curcellType = "mm10") 
hgGeneProp = attach.estContaminatedDist( df.geneProp = hgGeneProp, z = ori$z$call, dcon = dcon, curcellType = "hg19" )

plt.cpr.mm = plot.estGeneExp( df.geneProp=mmGeneProp, color = "#619CFF" ) + ylab("Gene expression of an average mouse cell") + xlab("DecontX contamination distribution\nin human cells") + ggtitle("Mouse-specific genes")  + theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  


plt.cpr.hg = plot.estGeneExp( df.geneProp = hgGeneProp, color = "slateblue" )+ ylab("Gene expression of an average human cell") + xlab("DecontX contamination distribution\nin mouse cells") + ggtitle("Human-specific genes") + theme( plot.margin=unit(c(1, 3,1, 1) , "lines"))  


cat( "correlation is:", cor( mmGeneProp$real.averageExp.cpm, mmGeneProp$estContaminationDist.crossSameGenome) )  
cat( "correlation is:", cor( hgGeneProp$real.averageExp.cpm, hgGeneProp$estContaminationDist.crossSameGenome) )  

pdf("Sup_Fig_2.pdf", width = 8, height = 4)
grid.arrange( plt.cpr.mm, plt.cpr.hg, ncol = 2, top = "Supplementary Figure 2") 
dev.off()


plt.UMIall= plotUMI( df.csv = df.csv, legend.position = c(0.8, 0.8)) + xlab("hg19 UMI counts") + ylab("mm10 UMI counts")  + ggtitle("Original UMI counts\n(human-mouse mixture)")  
plt.UMIhg = plotUMI( df.csv = df.csv[ df.csv$call == "hg19", ], color = "slateblue")  + xlab("") + ylab("") + ggtitle("   ")  
plt.UMImm = plotUMI( df.csv = df.csv[ df.csv$call == "mm10", ], color = '#619CFF')  + xlab("") + ylab("") + ggtitle("  \n  ")

plt.UMIall= plotUMI( df.csv = df.csv, legend.position = "right"  ) + xlab("hg19 UMI counts") + ylab("mm10 UMI counts")  + ggtitle("Original UMI counts (human-mouse mixture)")     
pdf("Sup_Fig_1.pdf", width = 5, height = 5) 
grid.arrange( plt.UMIall, top = "Supplementary Figure 1") 
dev.off()
