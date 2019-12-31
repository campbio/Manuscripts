library(Matrix)


mat.mm = readMM("/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/raw_gene_bc_matrices/mm10/matrix.mtx")
 dim(mat.mm) 
##  [output] 27998 737280
mat.hg = readMM('/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/raw_gene_bc_matrices/hg19/matrix.mtx')
 dim(mat.hg) 
##  [output] 32738 737280
barcode = read.table("/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/raw_gene_bc_matrices/mm10/barcodes.tsv", sep="\t")
## length(barcode) == ncol(mat.mm) == ncol(mat.hg)

# mat.mm.filtered = readMM('/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/filtered_gene_bc_matrices/mm10/matrix.mtx')
## dim(mat.mm.filtered) 
##   [output]: 27998  6657 
# mat.hg.filtered = readMM('/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/filtered_gene_bc_matrices/hg19/matrix.mtx') #
## dim(mat.hg.filtered)
## [output]:  32738  6893 

gene.mm = read.table("/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/raw_gene_bc_matrices/mm10/genes.tsv", sep='\t')
gene.hg = read.table("/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/raw_gene_bc_matrices/hg19/genes.tsv", sep='\t')

gem_classification = read.csv("/restricted/projectnb/camplab/home/syyang/contamination/data/10xWeb-12k/analysis/gem_classification.csv", header=TRUE)
allcells.index =  barcode[,1] %in% gem_classification$barcode
sum(allcells.index)

mat.allcells.mm = mat.mm[, allcells.index]
rmat.allcells.mm = as.matrix(mat.allcells.mm)
colnames(rmat.allcells.mm) = barcode[allcells.index,1]
rownames(rmat.allcells.mm) = gene.mm$V2

matconvert = function(sparse.mat, cell.filter.index ,  gene.name, barcode.name) {
	if ( class(cell.filter.index) == 'logical' & length(cell.filter.index) != ncol(sparse.mat)  ) { stop("wrong cell.filter.index") }
	mat.filtered = sparse.mat[, allcells.index ]
	barcode.filtered = barcode.name[cell.filter.index]
	rmat.filtered = as.matrix(mat.filtered)
	
	colnames(rmat.filtered) = barcode.filtered
	rownames(rmat.filtered) = gene.name
	return(rmat.filtered)
}
rmat.allcells.mm = matconvert(sparse.mat=mat.mm, cell.filter.index=allcells.index, gene.name=gene.mm$V2, barcode.name=barcode[,1])
rmat.allcells.hg = matconvert(sparse.mat=mat.hg, cell.filter.index=allcells.index, gene.name=gene.hg$V2, barcode.name=barcode[,1])

rmat.allcells = rbind(rmat.allcells.mm, rmat.allcells.hg)
saveRDS(rmat.allcells, '/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/10xWeb-12k-ana/allcells/allcells.rds')

specificCells.name = gem_classification$barcode[gem_classification$call != "Multiplet"]
rmat.specificCells = rmat.allcells[ , colnames(rmat.allcells) %in% specificCells.name  ]
saveRDS(rmat.specificCells, "/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/10xWeb-12k-ana/specific_cells/specificCells.rds")

# prepare cell label using gem_classification file from 10X Genomics 
library(plyr)
df.1 = data.frame( rmat.colname = colnames(rmat.allcells) )
df = merge( df.1, gem_classification, by.x="rmat.colname", by.y="barcode" ) 
df$z = factor(df$call)
levels(df$z)
df$z = revalue(df$z, c("Multiplet" =3, "hg19"=1, "mm10"=2 ) )
saveRDS(df, "/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/10xWeb-12k-ana/allcells/zTable.rds")

df.specificCells = data.frame( barcode.specificCells = colnames(rmat.specificCells) ) 
df.specificCells = merge( df.specificCells, df,  by.x="barcode.specificCells", by.y="rmat.colname") 
saveRDS(df.specificCells, "/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/10xWeb-12k-ana/specific_cells/zTable.rds")

rmat.specificCells.list = list("count" = rmat.specificCells, "z.tb" = df.specificCells)
saveRDS(rmat.specificCells.list, "/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/10xWeb-12k-ana/specific_cells/specificCellsWithLabel.rds")



# use species genome to get native expressed count matrix real.rmat
rmat.specificCells.list$z.tb$z.int = as.integer(as.character(rmat.specificCells.list$z.tb$z))
real.rmat = matrix(0, ncol=ncol(rmat.specificCells.list$count), nrow=nrow(rmat.specificCells.list$count) )
colnames(real.rmat) = colnames(rmat.specificCells.list$count)   
rownames(real.rmat) = rownames(rmat.specificCells.list$count) 

hg.index = grep("^hg19", rownames(real.rmat), perl=TRUE) 
real.rmat[ hg.index, rmat.specificCells.list$z.tb$call == "hg19"   ] = rmat.specificCells.list$count[ hg.index, rmat.specificCells.list$z.tb$call == "hg19"  ]
mm.index = grep("^mm10", rownames(real.rmat), perl=TRUE) 
real.rmat[ mm.index, rmat.specificCells.list$z.tb$call == "mm10"   ] = rmat.specificCells.list$count[ mm.index, rmat.specificCells.list$z.tb$call == "mm10"  ]

  
rmat.specificCells.list = append(rmat.specificCells.list, list("real.rmat"=real.rmat) ) 
saveRDS(rmat.specificCells.list, "/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/10xWeb-12k-ana/specific_cells/specificCellsWithLabel.rds")





