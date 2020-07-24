#!/usr/bin/env Rscript

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}

## specify package
pkg <- "QUILT"
root_dir <- getwd()

## if does not exist make releases folder
dir.create("releases", showWarnings = FALSE)


## documentation
devtools::document(pkg = pkg)

## make the tarball
print("Build")
package_tarball <- devtools::build(pkg = pkg, manual = TRUE)
version <- unlist(
    strsplit(unlist(strsplit(basename(package_tarball), "QUILT_")), ".tar.gz")
)

## move tarball to releases
release_package_tarball <- file.path("releases", paste0("QUILT_", version, ".tar.gz"))
system(paste0("mv '", package_tarball, "' ", release_package_tarball))
package_tarball <- release_package_tarball

## install from tarball
result <- tryCatch({
    install.packages(package_tarball)
}, warning = function(w) {
    stop(paste0(pkg, " installation failed"))
}, error = function(e) {
    stop(paste0(pkg, " installation failed"))       
})

## build PDF - disabled for now, pdflatex broken on rescomp
## pdf_name <- file.path("releases", paste0("QUILT_", version, ".pdf"))
## args = c(
##     "CMD", "Rd2pdf", pkg, "--batch", "--force", "--no-preview",
##     paste0("--output=", pdf_name)
## )
## system2("R", args)
