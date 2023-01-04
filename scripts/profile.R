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




library("QUILT")
setwd("/well/davies/users/dcc832/hrc_test_2022_12_15/")
CHR="chr20"
REGIONSTART=10000001
REGIONEND=12000000
BUFFER=500000
reference_sample_file = "hg38_liftover_chr20.samples"
genetic_map_file = "CEU-chr20-final.b38.txt.gz"

use_hapMatcherR <- Sys.getenv("USE_HAPMATCHER_R")
OUTPUTDIR <- Sys.getenv("OUTPUTDIR")
profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)
profile_start <- Sys.time()

################## PROFILE HERE
QUILT(
    outputdir = OUTPUTDIR,
    chr=CHR,
    regionStart=REGIONSTART,
    regionEnd=REGIONEND,
    buffer=BUFFER,
    bamlist="bamlist.txt",
    posfile="posfile.txt",
    sampleNames_file="sampleNames.txt",
    nGibbsSamples=3,
    n_seek_its=2,
    use_mspbwt = FALSE,
    use_hapMatcherR = use_hapMatcherR
)
################## END OF PROFILE

profile_end <- Sys.time()
setwd(stitch_dir)
Rprof(NULL)
pd <- readProfileData(profout)

output_plot <- "profile.pdf"
pdf(output_plot, height = 24, width = 24)
par(mfrow = c(3, 1), oma=c(0, 0, 3, 0))
flameGraph(pd, order = "time", main = "Time")
flameGraph(pd, order = "alpha", main = "Alphabetically")
flameGraph(pd, order = "hot", main = "Hot path")
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
