#!/usr/bin/env Rscript --vanilla

if (!require("devtools", character.only = TRUE)) {
    install.packages("devtools", dependencies = TRUE)
}

if (!require("singleCellTK", character.only = TRUE)) {
    devtools::install_github("compbiomed/singleCellTK", dependencies = TRUE)
}
