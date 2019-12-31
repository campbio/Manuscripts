plotUMI = function( df.csv, color =  c("grey",'slateblue','#619CFF'), labels = list("hg19"="Human", "mm10"="Mouse"),   legend.position = "none" , vline=FALSE ) {
    df.csv$call = as.character( df.csv$call ) 
    uni.call = unique(df.csv$call) 
    for ( i in uni.call) { 
        df.csv$call[ df.csv$call == i ] = labels[[i]] 
    }
    if (vline == TRUE) {
	p = ggplot() + geom_hline(yintercept=0, color="grey", size = 0.5) + geom_vline(xintercept=0, color="grey", size = 0.5) 
    } else { 
	p = ggplot()
    }
    p  =  p + geom_point(data = df.csv, aes(x=hg19, y=mm10, color = call ),   size=0.5  ) + 

	    scale_color_manual(values=color, name ="Cell type" ) + 
            guides( colour= guide_legend( override.aes = list( size = 2.5) )  ) +  	

            theme(panel.background=element_rect(fill="white", color="grey"), 
		  panel.grid = element_line("grey"), legend.position=legend.position, 
		  legend.key=element_rect(fill="white", color="white"), 
		  panel.grid.minor = element_blank(),
		  text=element_text(size=10), 
		  axis.text.x = element_text(size=10), 
		  axis.text.y=element_text(size=10), 
		  legend.key.size = grid::unit(8, "mm"), 
		  legend.text = element_text(size=10)  )
    return(p) 	    
}

conpCal = function( summary.tb) { 
    ttl.UMI = summary.tb$hg19 + summary.tb$mm10 
    conp = apply( summary.tb[, c("hg19", "mm10") ] , 1, min ) 
    summary.tb$conp = conp / ttl.UMI 
    return( summary.tb) 
}


plotViolin = function( conp.tb ) {

    p = ggplot( conp.tb[ conp.tb$call != "Multiplet", ], aes(x = call, y = conp) ) + 
            geom_jitter( height = 0, size = 0.1, color = "black", alpha = 0.2) + 	
	    geom_violin( aes( fill = call) , trim = TRUE, scale = "width"  ) + 
    ylab("Proportion of exogenous transcripts") + 
    ggtitle("Distribution of contamination level") + 
    scale_fill_manual( values = c('slateblue','#619CFF') ) + 
    scale_x_discrete( labels = c("mm10" = "Mouse", "hg19" = "Human" ) ) + 
            theme(panel.background=element_rect(fill="white", color="grey"), 
		  panel.grid = element_line("grey"), legend.position="none", 
		  legend.key=element_rect(fill="white", color="white"), 
		  panel.grid.minor = element_blank(),
		  text=element_text(size=10), 
		  axis.text.x = element_text(size=10), 
		  axis.text.y=element_text(size=10), 
		  legend.key.size = grid::unit(8, "mm"), 
		  legend.text = element_text(size=10)  )
	return(p) 
	    
}

df.Conp = function( ori, dcon) {
	# no filtering on the original count matrix
	N.by.C = colSums( ori$count) 
	tru.contaminationProportion = 1- colSums( ori$real.rmat) / N.by.C
	est.contaminationProportion = dcon$resList$estConp 
	df.conp = data.frame( tru.conp = tru.contaminationProportion, 
		     	      est.conp = est.contaminationProportion,
		 	      z = ori$z.tb$z.int,
			      call = ori$z.tb$call )
	return( df.conp) 
}	


colStat = function( dcon, ori ) { 
    mm10.rowIndex = grep( "^mm10", rownames( dcon$resList$estNativeCounts) , perl=TRUE) 
    hg19.rowIndex = grep( "^hg19", rownames( dcon$resList$estNativeCounts) , perl=TRUE)
    
    mm10.totalN = colSums( dcon$resList$estNativeCounts[ mm10.rowIndex , ] )
    hg19.totalN = colSums( dcon$resList$estNativeCounts[ hg19.rowIndex , ]  ) 

    decon.df = data.frame( call = as.character(ori$z.tb$call), hg19 = hg19.totalN, mm10 = mm10.totalN ) 

    mm10.rowIndex.ori = grep( "^mm10", rownames(ori$count), perl=TRUE) 
    hg19.rowIndex.ori = grep( "^hg19", rownames(ori$count), perl=TRUE) 

    ori.df = data.frame( call = as.character(ori$z.tb$call), 
			 hg19 = colSums( ori$count[ hg19.rowIndex.ori, ] ) , 
			 mm10 = colSums( ori$count[ mm10.rowIndex.ori, ] ) ) 
    
    counts.noEmpGene = ori$count [ rowSums( ori$count ) > 0 , ] 
    counts.mouseGene = counts.noEmpGene[ mm10.rowIndex, ] 
    counts.humanGene = counts.noEmpGene[ hg19.rowIndex, ] 

    return( list("ori.df" = ori.df,  "decon.df" =decon.df, "counts.mouseGene"=counts.mouseGene, "counts.humanGene"=counts.humanGene   )  ) 
}


plot.performance = function( df.conp, type="none") {
	if ( type != "none" ) { df.conp = df.conp[  df.conp$call == type,  ] } 
	p = ggplot( df.conp , aes( x = tru.conp, y = est.conp) ) + geom_point( aes( col= call) ) +
		geom_abline( intercept = 0, slope=1, color="grey" ) +
	      	ylab("Estimated contamination proportion") + xlab("Proportion of impossible transcripts") + labs(color=" ")  + 
		theme(panel.background=element_rect(fill="white", color="grey"), panel.grid = element_line("grey"), 
		      #legend.position=c(1.15, 0.9), plot.margin=unit(c(2, 8,2, 2) , "lines"), 
		  panel.grid.minor = element_blank(),
		  text=element_text(size=10), 
		  axis.text.x = element_text(size=10), 
		  axis.text.y=element_text(size=10), 
		  legend.key.size = grid::unit(8, "mm"),  
		  legend.text = element_text(size=10)  )
	return(p)  
} 



rmse = function( x, y) { 
	rootMeanSquareError = sqrt( mean( (x-y)^2 ) ) 
	return(rootMeanSquareError) 
}


geneProp = function( counts, z, curcellType, otherType ) {
    counts.ccType = counts[ , z == curcellType ] 
    average.exp = rowSums( counts.ccType ) /  sum( counts.ccType ) 

    counts.cpmNor = celda::normalizeCounts( counts = counts, "cpm") 
    counts.cpmNor.ccType = counts.cpmNor[, z == curcellType ] 
    average.exp.cpm = rowSums( counts.cpmNor.ccType ) / ncol( counts.cpmNor.ccType ) 

    counts.otrType = counts[, z == otherType ] 
    prop.inConTranscripts = rowSums( counts.otrType ) / sum( counts.otrType) 


    df = data.frame(real.averageExp = average.exp,  
		    real.averageExp.cpm = average.exp.cpm, 
		    prop.inContaminatedTranscripts = prop.inConTranscripts ) 
	return(df) 
}

attach.estContaminatedDist = function( df.geneProp, z, dcon, curcellType )  { 
    genome.rowIndex = grep( paste0("^", curcellType), rownames(dcon$resList$estNativeCounts) , perl=TRUE) 

    estContaminatedCounts = dcon$resList$estNativeCounts[, z == curcellType]
    estContaminatedCounts.crossSameGenome = estContaminatedCounts[ genome.rowIndex, ]
    estContaminationDist.crossSameGenome = rowSums(estContaminatedCounts.crossSameGenome )/ 
	    sum( estContaminatedCounts.crossSameGenome ) 
    df.geneProp$estContaminationDist.crossSameGenome = estContaminationDist.crossSameGenome
    return(df.geneProp) 
}



plot.geneExp = function( df.exp, color ='red'   ) { 
	p = ggplot( df.exp , aes( y= real.averageExp.cpm,  x=  prop.inContaminatedTranscripts  ) ) + geom_point( size = 1, color = color) + 
		scale_x_continuous( breaks = c(0, 0.005, 0.010, 0.015) ) + 
		theme(panel.background=element_rect(fill="white", color="grey"), panel.grid = element_line("grey"), 
		      #legend.position=c(1.15, 0.9), plot.margin=unit(c(2, 8,2, 2) , "lines"), 
		  panel.grid.minor = element_blank(),
		  text=element_text(size=10), 
		  axis.text.x = element_text(size=10), 
		  axis.text.y=element_text(size=10), 
		  legend.key.size = grid::unit(8, "mm"),  
		  legend.text = element_text(size=10)  )
       return( p) 
}

plot.estGeneExp = function( df.geneProp, color = "red" ) { 
    p = ggplot( df.geneProp , aes( y = real.averageExp.cpm, x = estContaminationDist.crossSameGenome ) ) + 
	   geom_point( size = 1, color = color) + 
           theme(panel.background=element_rect(fill="white", color="grey"), panel.grid = element_line("grey"), 
		      #legend.position=c(1.15, 0.9), plot.margin=unit(c(2, 8,2, 2) , "lines"), 
		  panel.grid.minor = element_blank(),
		  text=element_text(size=10), 
		  axis.text.x = element_text(size=10), 
		  axis.text.y=element_text(size=10), 
		  legend.key.size = grid::unit(8, "mm"),  
		  legend.text = element_text(size=10)  )
       return( p) 
}
