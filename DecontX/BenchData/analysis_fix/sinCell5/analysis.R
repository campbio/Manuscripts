dcon.10x = readRDS("dcon.10x.rds")
dcon.celseq2_p1 = readRDS("dcon.celseq2_p1.rds")
dcon.celseq2_p2 = readRDS("dcon.celseq2_p2.rds")
dcon.celseq2_p3 = readRDS("dcon.celseq2_p3.rds")



cbres = function( ...) { 
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


res = cbres( dcon10x =dcon.10x, celseq2_p1 = dcon.celseq2_p1, celseq2_p2 = dcon.celseq2_p2, celseq2_p3 = dcon.celseq2_p3) 

res$summaryStat

library(ggplot2)
p.violin = ggplot( res$cdf, aes( x = dataset, y = dcontxEstP )  ) + 
	geom_jitter( width = 0.2, alpha = 1, size = 0.3, aes(color = dataset)  ) +
	 geom_violin( trim=T, scale = "width", fill = "grey", alpha = 0.5 )  


 ggsave( plot = p.violin, filename ="decontEstP.png" ) 

