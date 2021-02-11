# Package installations

install.packages("renv")
library(renv)

# If you are attempting to reproduce this analysis, you can
# 1) install/load `renv`,
# 2) copy/paste the `renv.lock` file into your current working directory,
# 3) run the `renv::restore()` command to automatically install all of the
#    packages (except Celda) with the same versions used for the analyses, and
# 4) install Celda version v1.7.3 from GitHub.

# Make sure the 'renv.lock' file is in your current working directory.
renv::restore()

# run this if needed
devtools::install_github("campbio/celda@v1.7.3")

# Load R packages
library(celda)
library(Seurat)
library(TENxPBMCData)
library(ggplot2)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(data.table)
library(scales)
library(cowplot)
library(gtable)
library(uwot)
library(scran)
library(scuttle)
library(scater)

