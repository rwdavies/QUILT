#!/usr/bin/env Rscript

required_packages <- c("optparse")
for(package in required_packages) {
    if (!suppressPackageStartupMessages(require(package, character.only = TRUE))) {
        out <- install.packages(package, repos="http://cran.rstudio.com/")
        out <- require(package, character.only = TRUE)
        if (!out) {
            stop(paste0("Failed to install package:", package))
        }
    }
}
