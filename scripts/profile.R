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
    library("QUILT")
setwd("/data/smew1/rdavies/riyan_debug_2021_03_15")
i_chr <- 20
set.seed(i_chr)
regionStart <- regionEnd <- buffer <- NA
chr <- paste0("chr", i_chr)
bams <- dir("bams")
bams <- bams[-grep(".bai", bams)][1]
write.table(
    matrix(paste0("bams/", bams), ncol = 1),
    file = "bamlist.txt",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = ""
)
## 
QUILT(
    outputdir = "/data/smew1/rdavies/riyan_debug_2021_03_15",
    chr = chr,
    regionStart = regionStart,
    regionEnd = regionEnd,
    buffer = buffer,
    bamlist = "bamlist.txt",
    override_default_params_for_small_ref_panel = TRUE,
    nCores = 1,
    nGibbsSamples = 3,
    print_extra_timing_information = TRUE
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
