load("/restricted/projectnb/camplab/home/syyang/contamination/data/cellBenchData/CellBench_data_08112019/data/sincell_with_class.RData")




# Compare contamination% in different protocols
dcon.10x = readRDS("dcon.10x.rds")
dcon.celseq2 = readRDS("dcon.celseq2.rds")
dcon.dropseq = readRDS("dcon.dropseq.rds")



rbind.res = function( ...) { 
	List <- function(...) {
		 names <- as.list(substitute(list(...)))[-1L]
                  setNames(list(...), names)
	}
	l = List(...) 
	ldf = lapply( l, FUN=function(i) { data.frame( dcontxEstP = i$resList$estConp ) } ) 
	sumstat = lapply(l, FUN=function(i) summary(i$resList$estConp) )
	cdf = do.call( rbind, ldf) 
	cdf$dataset = gsub("\\.[^\\.]*$", "",  rownames(cdf) ) 
	return(list("cdf" = cdf, "summaryStat" = sumstat))
}


res = rbind.res( dcon10x =dcon.10x, celseq2 = dcon.celseq2, dropseq = dcon.dropseq) 

res$summaryStat

library(ggplot2)
p.violin = ggplot( res$cdf, aes( x = dataset, y = dcontxEstP )  ) + 
	geom_jitter( width = 0.2, alpha = 1, size = 0.3, aes(color = dataset)  ) +
	 geom_violin( trim=T, scale = "width", fill = "grey", alpha = 0.5 )  


 ggsave( plot = p.violin, filename = "decontEstP.png")




 # Compare contamination% in singlets VS in doublets (SNG/DBL estimation is from benchdata-paper) 
 summary (dcon.10x$resList$estConp[sce_sc_10x_qc@colData$demuxlet_cls == "DBL"  ] )
 summary (dcon.10x$resList$estConp[sce_sc_10x_qc@colData$demuxlet_cls == "SNG"  ] )



 summary( dcon.celseq2$resList$estConp[sce_sc_CELseq2_qc@colData$demuxlet_cls == "DBL"  ] ) 
 summary( dcon.celseq2$resList$estConp[sce_sc_CELseq2_qc@colData$demuxlet_cls == "SNG"  ] ) 



 summary( dcon.dropseq$resList$estConp[sce_sc_Dropseq_qc@colData$demuxlet_cls == "DBL"  ] )
 summary( dcon.dropseq$resList$estConp[sce_sc_Dropseq_qc@colData$demuxlet_cls == "SNG"  ] )
