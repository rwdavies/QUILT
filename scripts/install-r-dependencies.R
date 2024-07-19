#!/usr/bin/env Rscript

required_packages <- c("proftools", "Rcpp", "RcppArmadillo", "optparse", "devtools", "testthat", "roxygen2", "data.table", "RcppEigen", "microbenchmark")
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
    install.packages("https://github.com/rwdavies/rrbgen/releases/download/0.0.6/rrbgen_0.0.6.tar.gz", repos=NULL)

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
        install.packages("https://github.com/rwdavies/STITCH/releases/download/1.7.0/STITCH_1.7.0.tar.gz", repos=NULL)
    }
}


if (!suppressPackageStartupMessages(require("mspbwt"))) {
    install.packages("https://github.com/rwdavies/mspbwt/releases/download/0.1.0/mspbwt_0.1.0.tar.gz", repos=NULL)
}
