#!/usr/bin/env Rscript

library("proftools")
library("QUILT")
options(scipen = 999)

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
        stitch_dir <- getwd()
    }
}

profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)

profile_start <- Sys.time()
##################
prepared_reference_filename <- Sys.getenv("prepared_reference_filename")
##     phasefile = phasefile,
setwd("/data/smew1/rdavies/quilt_data/hrc_2021_04_20/")
dir <- "/data/smew1/rdavies/quilt_data/hrc_2021_04_20/2021_04_20_bams"
f <- dir(dir)[grep(".1.0.bam", dir(dir))]
f <- f[-grep("bai", f)]
bamlist <- tempfile()
write.table(matrix(paste0(dir, "/", f), ncol = 1), file = bamlist, row.names = FALSE, col.names = FALSE, quote = FALSE)
QUILT(
    outputdir = "/data/smew1/rdavies/quilt_data/hrc_2021_04_20/",
    chr = "chr20",
    regionStart = 2000000,
    regionEnd = 4000000,
    buffer = 500000,
    bamlist = bamlist,
    bqFilter = 10,
    nCores = 1,
    prepared_reference_filename = prepared_reference_filename,
    nGibbsSamples = 3,
    n_seek_its = 2
)
##################

profile_end <- Sys.time()

Rprof(NULL)
pd <- readProfileData(profout)
title <- Sys.getenv("TITLE")

output_plot <- Sys.getenv("OUTPUT_PLOT")
if (output_plot == "") {
    setwd(stitch_dir)
    output_plot <- file.path("profile.pdf")
}
pdf(output_plot, height = 24, width = 24)
par(mfrow = c(3, 1), oma=c(0, 0, 3, 0))
flameGraph(pd, order = "time", main = "Time")
flameGraph(pd, order = "alpha", main = "Alphabetically")
flameGraph(pd, order = "hot", main = "Hot path")
title(title, outer=TRUE)
dev.off()

print(profile_end - profile_start)
