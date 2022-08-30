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
dir <- "/well/davies/users/dcc832/QUILT_nicola_testing_2022_08_26/"
f <- function(x) file.path(dir, x)
CHR="chr20"
REGIONSTART=10000001
REGIONEND=12000000
BUFFER=500000
QUILT(
    outputdir=f("quilt_mspbwt_test"),
    chr=CHR,
    regionStart=REGIONSTART,
    regionEnd=REGIONEND,
    buffer=BUFFER,
    bamlist=f("bamlist.txt"),
    posfile=f("posfile.txt"),
    use_mspbwt=TRUE,
    make_plots=FALSE,
    block_gibbs_iterations=3,
    n_gibbs_burn_in_its=6,
    n_seek_its=5,
    nGibbsSamples=2,
    Ksubset = 100,
    Knew = 100
)
##################

profile_end <- Sys.time()

Rprof(NULL)
pd <- readProfileData(profout)
title <- Sys.getenv("TITLE")

output_plot <- Sys.getenv("OUTPUT_PLOT")
if (output_plot == "") {
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
