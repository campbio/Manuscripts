#!/usr/bin/env Rscript --vanilla

if (!require("devtools", character.only = TRUE)) {
    install.packages("devtools", dependencies = TRUE)
}

if (!require("optparse", character.only = TRUE)) {
    install.packages("optparse", dependencies = TRUE, repos = "http://cran.us.r-project.org")
}

if (!require("rmarkdown", character.only = TRUE)) {
    install.packages("rmarkdown", dependencies = TRUE,repos = "http://cran.us.r-project.org")
}

if (!require("singleCellTK", character.only = TRUE)) {
    devtools::install_github("compbiomed/singleCellTK", dependencies = TRUE)
}
