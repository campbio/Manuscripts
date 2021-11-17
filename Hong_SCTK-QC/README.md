# Comprehensive generation, visualization, and reporting of quality control metrics for single-cell RNA sequencing data

The SCTK-QC pipeline within the [singleCellTK](https://github.com/compbiomed/singleCellTK) R package can be used tostreamline and standardize QC analysis for scRNA-seq data across a variety of platforms. Features in this pipeline include the ability to import data from 11 different preprocessing tools or file formats, perform empty droplet detection with 2 different algorithms, generate standard quality control metrics such as number of UMIs per cell or the percentage of mitochondrial counts, predict doublets using 6 different algorithms, and estimate ambient RNA. QC data can be exported to R and/or Python objects used in popular down-stream workflows.

## Repository Structure

* Scripts folder contains the shell and R scripts which will be used to generate the results and figures.
* Data folder contains R objects in RDS formats.
* Results folder contains a html file containing the figures presented in the study. Additionally, upon running the pipeline, the output will be stored in this directory as well by default.

## How to regenerate the results

To regenerate the results, please follow these steps:
1. Set the file path of this repository as the `input` in the `execute_SCTK-QC.sh` script.
2. Set the desired file path for the output of the pipeline as the `output` in the `execute_SCTK-QC.sh` script.
3. Execute the `execute_SCTK-QC.sh` script with the following command in your terminal:

 `bash execute_SCTK-QC.sh`

The outputs of the pipeline will be stored in the specified `output` directory.
