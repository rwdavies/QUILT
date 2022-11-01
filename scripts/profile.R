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

setwd(stitch_dir)
profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)





dir <- "/well/davies/users/dcc832/QUILT_testing_2022_10_14/"
setwd(dir)
f <- function(x) file.path(dir, x)
CHR="chr20"
REGIONSTART=10000001
REGIONEND=12000000
BUFFER=500000


if (1 == 0) {
    
    setwd("/well/davies/users/dcc832/QUILT_testing_2022_10_14/")
    REF_DIR="/well/davies/shared/hrc_ega/sinan_2020_06_20"
    RECOMB_DIR="/well/davies/shared/recomb/CEU/"
    QUILT_prepare_reference(
        outputdir = f("quilt_mspbwt_test2"),
        reference_exclude_samplelist_file = "excludelist.txt",
        nGen=100,
        reference_haplotype_file="temp3.hap.gz",
        reference_legend_file="temp3.legend.gz",
        reference_sample_file=file.path(REF_DIR, "hg38_liftover_chr20.samples"),
        genetic_map_file=file.path(RECOMB_DIR, paste0("CEU-", CHR, "-final.b38.txt.gz")),
        chr=CHR,
        regionStart=REGIONSTART,
        regionEnd=REGIONEND,
        buffer=BUFFER,
        use_mspbwt = TRUE,
        mspbwt_nindices = 4L
    )


}

profile_start <- Sys.time()

##################
QUILT(
    outputdir=f("quilt_mspbwt_test2"),
    chr=CHR,
    regionStart=REGIONSTART,
    regionEnd=REGIONEND,
    buffer=BUFFER,
    bamlist=f("bamlist.txt"),
    posfile=f("posfile.txt"),
    phasefile=f("phasefile.txt"),    
    use_mspbwt=TRUE,
    make_plots=FALSE,
    use_sample_is_diploid = TRUE
)

##################
    ## n_seek_its=3,
    ## nGibbsSamples=3,
    ## Ksubset = 400,
    ## Knew = 400,
    ## block_gibbs_iterations=3,
    ## n_gibbs_burn_in_its=6,

profile_end <- Sys.time()

setwd(stitch_dir)
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


## use_sample_is_diploid = FALSE
## Time difference of 53.16302 secs
## mv profile.pdf profile.nodiploid.pdf

## use sample_is_diploid = TRUE
## Time difference of 33.62765 secs
## bash-4.2$ mv profile.pdf profile.usediploid.pdf


## Time difference of 41.5688 secsmv p
## bash-4.2$ mv profile.pdf profile.nodiploid.pdf

## Time difference of 29.70805 secs
## bash-4.2$ mv profile.pdf profile.nodiploid.nophasefile.pdf

##Time difference of 24.52199 secs
## bash-4.2$ mv profile.pdf profile.usediploid.nophasefile.pdf


## use diploid
## Time difference of 45.95948 secs
## no use diploid

## no diploid
## Time difference of 53.44265 secs


## with -O3 flag on
## Time difference of 1.377417 mins
## Time difference of 1.428216 min
## ./scripts/build-and-install.R && ./scripts/profile.R && mv profile.pdf profile.yesO3.pdf
## example ## g++ -std=gnu++11 -I"/well/davies/users/dcc832/bin/R-4.1.3/include" -DNDEBUG -I./ -I'/gpfs3/well/davies/users/dcc832/bin/R-4.1.3/library/Rcpp/include' -I'/gpfs3/well/davies/users/dcc832/bin/R-4.1.3/library/RcppArmadillo/include' -I/usr/local/include  -O3 -fpic  -g -O2  -c test.cpp -o test.o



## with no -O3 flag
## Time difference of 1.28888 mins
## ./scripts/profile.R && mv profile.pdf profile.O2.pdf
