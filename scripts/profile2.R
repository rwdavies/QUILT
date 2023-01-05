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




REBUILD <- as.logical(Sys.getenv("REBUILD"))
VERSION <- Sys.getenv("VERSION")
METHOD <- Sys.getenv("METHOD")
OBJECT <- Sys.getenv("OBJECT")
WHAT <- Sys.getenv("WHAT")

OUTPUTDIR <- file.path(
    "/well/davies/users/dcc832/hrc_test_2022_12_15/",
    paste0("output_", VERSION, "_", METHOD, "_", OBJECT)
)
OUTPUTFILE <- file.path(OUTPUTDIR, "RData/QUILT_prepared_reference.chr20.10000001.12000000.RData")
OUTPUTPDF <- paste0("~/proj/QUILT/scratch/profile_", VERSION, "_", METHOD, "_", OBJECT, ".pdf")


if (OBJECT == "regular") {
    use_hapMatcherR <- FALSE
} else {
    use_hapMatcherR <- TRUE
}


if (METHOD == "full") {
    
    use_mspbwt <- FALSE
    zilong <- FALSE
    use_pbwt_index <- FALSE

} else if (METHOD == "mspbwt") {

    use_mspbwt <- TRUE
    zilong <- FALSE
    use_pbwt_index <- FALSE
    
} else if (METHOD == "zilong") {

    use_mspbwt <- FALSE
    zilong <- TRUE
    use_pbwt_index <- TRUE
    
}

if (WHAT == "BUILD") {

if (VERSION == "QUILT1") {
    
    ##install.packages("~/proj/QUILT/releases/QUILT_1.0.4.tar.gz")
    stopifnot(sessionInfo()["otherPkgs"][[1]][["QUILT"]][["Version"]] == "1.0.4")
    
    if (REBUILD | !file.exists(OUTPUTFILE)) {
        QUILT_prepare_reference(
            outputdir = OUTPUTDIR,
            reference_exclude_samplelist_file = "excludelist.txt",
            chr=CHR,
            nGen=100,
            reference_haplotype_file="temp3.hap.gz",
            reference_legend_file="temp3.legend.gz",
            reference_sample_file=reference_sample_file,
            genetic_map_file=genetic_map_file,
            regionStart=REGIONSTART,
            regionEnd=REGIONEND,
            buffer=BUFFER
        )
    }
} else {

    ## install.packages("~/proj/QUILT/releases/QUILT_2.0.0.tar.gz")
    stopifnot(sessionInfo()["otherPkgs"][[1]][["QUILT"]][["Version"]] == "2.0.0")
              
    if (REBUILD | !file.exists(OUTPUTFILE)) {

        QUILT_prepare_reference(
            outputdir = OUTPUTDIR,
            reference_exclude_samplelist_file = "excludelist.txt",
            chr=CHR,
            nGen=100,
            reference_haplotype_file="temp3.hap.gz",
            reference_legend_file="temp3.legend.gz",
            reference_sample_file=reference_sample_file,
            genetic_map_file=genetic_map_file,
            regionStart=REGIONSTART,
            regionEnd=REGIONEND,
            buffer=BUFFER,
            use_mspbwt=use_mspbwt,
            use_pbwt_index=use_pbwt_index,
            reference_vcf_file="temp2.vcf.gz",
            use_hapMatcherR = use_hapMatcherR
        )
    }
    
}
    stop("Done building")

}

profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)
profile_start <- Sys.time()

################## PROFILE HERE
if (VERSION == "QUILT1") {

    QUILT(
        outputdir = OUTPUTDIR,
        chr=CHR,
        regionStart=REGIONSTART,
        regionEnd=REGIONEND,
        buffer=BUFFER,
        bamlist="bamlist.txt",
        posfile="posfile.txt",
        sampleNames_file="sampleNames.txt"
    )
    
} else {

    QUILT(
        outputdir = OUTPUTDIR,
        chr=CHR,
        regionStart=REGIONSTART,
        regionEnd=REGIONEND,
        buffer=BUFFER,
        bamlist="bamlist.txt",
        posfile="posfile.txt",
        sampleNames_file="sampleNames.txt",
        use_mspbwt = use_mspbwt,
        zilong = zilong,
        use_hapMatcherR = use_hapMatcherR,
        reference_vcf_file="temp2.vcf.gz"        
    )

}
################## END OF PROFILE

profile_end <- Sys.time()
setwd(stitch_dir)
Rprof(NULL)
pd <- readProfileData(profout)

pdf(OUTPUTPDF, height = 24, width = 24)
par(mfrow = c(3, 1), oma=c(0, 0, 3, 0))
flameGraph(pd, order = "time", main = "Time")
flameGraph(pd, order = "alpha", main = "Alphabetically")
flameGraph(pd, order = "hot", main = "Hot path")
dev.off()

print(profile_end - profile_start)

