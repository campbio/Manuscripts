specific.cells = readRDS("/restricted/projectnb/camplab/home/syyang/contamination/contamination-output/10xWeb-12k-ana/specific_cells/specificCellsWithLabelRownameEmsemble.rds")
count = specific.cells$count[ rowSums(specific.cells$count) > 0 , ]
#count = specific.cells$count[ rowSums(specific.cells$count > 2) > 2, ]
z = specific.cells$z$z
z = as.integer(as.character(z)) 

max.iter = 60

RNGkind(sample.kind = "Rounding")
library(devtools)
load_all("~/gitProjects/celda_irisapo", recompile=T)   # on fix branch, celda version 1.3.1


set.seed(12345) 
dcon <- decontX(count , z = z,  maxIter = max.iter)


saveRDS(dcon, "dcon.rds" )

