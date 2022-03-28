NFEATURE=2000
CurrentDir=normalizePath(".")
SubSce_Dir = file.path(CurrentDir, paste("sce_nfeature", NFEATURE, sep="_") )
library(Seurat)

library(ExperimentHub)
eh <- ExperimentHub()

# cache data
data_hub = query(eh, "DuoClustering2018")
print(head(unique(mcols(data_hub)[,'title']), 3)) 

#datasets
## e.g., dataset_names =c("sce_filteredM3Drop10_Zhengmix8eq", "sce_filteredHVG10_Zhengmix8eq") 
allscenames_index = grep("sce_", unique(mcols(data_hub)[,'title'])) # 48 of them 
allscenames = unique(mcols(data_hub)[,'title'])[allscenames_index]


library(DuoClustering2018)

datasets =c("sce_filteredM3Drop10_Zhengmix8eq", "sce_filteredHVG10_Zhengmix8eq", "sce_filteredExpr10_Zhengmix8eq",
            "sce_filteredM3Drop10_Zhengmix4eq", "sce_filteredHVG10_Zhengmix4eq", "sce_filteredExpr10_Zhengmix4eq",
            "sce_filteredM3Drop10_Zhengmix4uneq", "sce_filteredHVG10_Zhengmix4uneq", "sce_filteredExpr10_Zhengmix4uneq", 
            "sce_filteredM3Drop10_Koh", "sce_filteredHVG10_Koh", "sce_filteredExpr10_Koh", 
            "sce_filteredM3Drop10_KohTCC", "sce_filteredHVG10_KohTCC", "sce_filteredExpr10_KohTCC", 
            "sce_filteredM3Drop10_Trapnell", "sce_filteredHVG10_Trapnell", "sce_filteredExpr10_Trapnell", 
            "sce_filteredM3Drop10_TrapnellTCC", "sce_filteredHVG10_TrapnellTCC", "sce_filteredExpr10_TrapnellTCC", 
            "sce_filteredM3Drop10_Kumar", "sce_filteredHVG10_Kumar", "sce_filteredExpr10_Kumar",
            "sce_filteredM3Drop10_KumarTCC", "sce_filteredHVG10_KumarTCC", "sce_filteredExpr10_KumarTCC")



for (scename in datasets) {
	  sce = get(scename)()
		dat <- counts(sce)

		data <- CreateSeuratObject(counts = dat, min.cells = 0,
															 min.features = 0, project = "scRNAseq", 
															 display.progress = TRUE) 	
		data <- NormalizeData(object = data, normalization.method = "LogNormalize",
													scale.factor = 1e4, display.progress = FALSE)

		sub_dat <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
		sub_genes <- VariableFeatures(sub_dat)

		sub_sce <- sce[ sub_genes, ]

		file_name <- file.path(SubSce_Dir, paste0(scename, ".rds"))
		saveRDS(sub_sce, file_name)

}
