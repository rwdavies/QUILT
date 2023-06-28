#!/usr/bin/env Rscript

required_packages <- c("proftools", "Rcpp", "RcppArmadillo", "optparse", "devtools", "testthat", "roxygen2", "data.table", "RcppEigen")
for(package in required_packages) {
    if (!suppressPackageStartupMessages(require(package, character.only = TRUE))) {
        out <- install.packages(package, repos="http://cran.rstudio.com/")
        out <- require(package, character.only = TRUE)
        if (!out) {
            stop(paste0("Failed to install package:", package))
        }
    }
}
if (!suppressPackageStartupMessages(require("rrbgen")))
    install.packages("https://github.com/rwdavies/rrbgen/raw/master/releases/rrbgen_0.0.4.tar.gz", repos=NULL)

    check <- as.logical(Sys.getenv("DEV_STITCH")) == TRUE
    if (is.na(check)) {
        check <- FALSE
    }
    if (check) {
        install_github("rwdavies/STITCH", subdir = "STITCH", upgrade = "never")
    }

if (!suppressPackageStartupMessages(require("STITCH"))) {
    check <- as.logical(Sys.getenv("DEV_STITCH")) == TRUE
    if (is.na(check)) {
        check <- FALSE
    }
    if (check) {
        install_github("rwdavies/STITCH", subdir = "STITCH", upgrade = "never")
    } else {
        install.packages("https://github.com/rwdavies/STITCH/releases/download/1.6.10/STITCH_1.6.10.tar.gz", repos=NULL)
    }
}



