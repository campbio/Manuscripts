import scrublet as rc 
import matplotlib.pyplot as plt 
import scipy.io 
from scipy.sparse import csc_matrix 
import numpy as np 

wd = "/restricted/projectnb/camplab/home/syyang/contamination/data/pbmc/4k/" 

counts = scipy.io.mmread( wd + 'data/matrix.mtx' ) 

geneIndex = ((counts > 2).sum( axis = 1 ) > 2 )  
counts_filter = counts.toarray()[ np.array(geneIndex).reshape(-1), : ] 

print( counts_filter.shape) 

counts_csc = csc_matrix( counts_filter.T ) 
genes = np.array( rc.load_genes( wd + 'data/genes.tsv', delimiter='\t', column=1)) [ np.array(geneIndex).reshape(-1)  ]


scrub = rc.Scrublet(counts_csc, expected_doublet_rate=0.06)


doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)


np.save( "doublet_scores.npy", doublet_scores  ) 
np.save( "predicted_doublets.npy", predicted_doublets * 1 )
