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
outputdir <- tempdir()
setwd("~/Google\ Drive/Papers/2020\ -\ STITCH\ Haplotagging\ lcWGS/reviewer_package/QUILT_review_package_2020_07_27/")
system("head -n1 bamlist.1.0.txt > bamlist.1.0.onesample.txt")
if (1 == 0) {
    load(file = "quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData")
        nMaxDH <- 2 ** 8 - 1
        ref_error <- 1e-3
        out <- make_rhb_t_equality(
            rhb_t = rhb_t,
            nMaxDH = nMaxDH,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
        distinctHapsB <- out[["distinctHapsB"]]
        distinctHapsIE <- out[["distinctHapsIE"]]            
        hapMatcher <- out[["hapMatcher"]]
        eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
        eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]
    
    save(
        eMatDH_special_grid_which,
        eMatDH_special_values_list,
        af, cM, cM_grid, distinctHapsB, distinctHapsIE, dl, grid, hapMatcher, inRegion2, L, L_grid, nGrids, nSNPs, pos, ref_alleleCount, ref_error, reference_samples, rh_in_L, rhb_t, sigmaCurrent_m,
        file = "quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData"
    )
}
##    
## ARGH have_truth_haplotypes <- FALSE
if (Sys.getenv("USE_PHASEFILE") == "TRUE") {
    phasefile <- "phase.chr20.2000001.4000000.txt"
} else {
    phasefile <- ""
}
QUILT(
    outputdir = "quilt_output",
    chr = "chr20",
    regionStart = 2000001,
    regionEnd = 4000000,
    buffer = 500000,
    bamlist = "bamlist.1.0.txt",
    posfile = "pos.chr20.2000001.4000000.txt",
    phasefile = phasefile,
    bqFilter = 10,
    nCores = 1,
    prepared_reference_filename = "quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData",
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
