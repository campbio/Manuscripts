## replace output with the absolute path of the folders on your machine
input="/path/to/directory/Hong_SCTK-QC"
output="/path/to/directory/Hong_SCTK-QC/Results"

# Install SCTK-QC
Rscript ${input}/Scripts/SCTK_install.R

# SCTK-QC pipeline

## The following code will execute the SCTK_runQC pipeline script. 
## Please refer to https://github.com/compbiomed/singleCellTK/blob/master/exec/SCTK_runQC_Documentation.md for details on each parameter. 

## PBMC1K, Gencode version 27, V2 Chemistry
Rscript ${input}/Scripts/SCTK_runQC.R \
-P SceRDS \
-s g27_pbmc1k_v2 \
-o ${output}/Pipeline_Output \
-F SCE \
-c ${input}/Data/gencode27_pbmc1k_v2.rds \
-d Cell \
-M FALSE

## PBMC1K, Gencode version 27, V3 Chemistry
Rscript ${input}/Scripts/SCTK_runQC.R \
-P SceRDS \
-s g27_pbmc1k_v3 \
-o ${output}/Pipeline_Output \
-F SCE \
-c ${input}/Data/gencode27_pbmc1k_v3.rds \
-d Cell \
-M FALSE

## PBMC1K, Gencode version 34, V2 Chemistry
Rscript ${input}/Scripts/SCTK_runQC.R \
-P SceRDS \
-s g34_pbmc1k_v2 \
-o ${output}/Pipeline_Output \
-F SCE \
-c ${input}/Data/gencode34_pbmc1k_v2.rds \
-d Cell \
-M FALSE

## PBMC1K, Gencode version 34, V3 Chemistry
Rscript ${input}/Scripts/SCTK_runQC.R \
-P SceRDS \
-s g34_pbmc1k_v3 \
-o ${output}/Pipeline_Output \
-F SCE \
-c ${input}/Data/gencode34_pbmc1k_v3.rds \
-d Cell \
-M FALSE

## SMART-Seq2, Replicate 1
Rscript ${input}/Scripts/SCTK_runQC.R \
-P SceRDS \
-s HCA_SCP424_PBMC1 \
-o ${output}/Pipeline_Output \
-F SCE \
-r NULL \
-c ${input}/Data/SM2_HCA_PBMC1.rds \
-d Cell \
-M FALSE

## SMART-Seq2, Replicate 2
Rscript ${input}/Scripts/SCTK_runQC.R \
-P SceRDS \
-s HCA_SCP424_PBMC2 \
-o ${output}/Pipeline_Output \
-F SCE \
-r NULL \
-c ${input}/Data/SM2_HCA_PBMC2.rds \
-d Cell \
-M FALSE

# Generate figures used in manuscript

R -e "rmarkdown::render('${input}/Scripts/SCTK_Generate_Figures.Rmd', output_dir = '${output}/Figures')" --args ${input} ${output}
