# The combined immune data is sorted PBMCs data
print("load combined immune data") 
# load original count matrix 
mat <- readRDS("/restricted/projectnb/camplab/projects/celda/Datasets/10X/Combined_Immune_Datasets/all_immune_counts.Rds")
mat.s <- mat[rowSums(mat > 2) > 2, colSums(mat > 2) >2 ] 

# load cell sequncing labels ( total 9 different types of cells)
mat.wlabel = readRDS("/restricted/projectnb/camplab/projects/celda/Datasets/10X/Combined_Immune_Datasets/all_immune_counts_tsne.Rds")
mat.s.label = mat.wlabel$filtered.immune.tsne$dataset


