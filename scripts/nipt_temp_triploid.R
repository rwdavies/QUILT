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
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))




## important global variables
ANALYSIS_DIR <- "/well/davies/users/dcc832/nipt_test_2023_09_04/"
REF_PREFIX <- "HRC" ## or use ONEKG for 1000 Genomes, should work I think
## REF_PREFIX <- "ONEKG"
CHR <- "chr20"
REGIONSTART <- 10000000
REGIONEND <- 15000000
BUFFER <- 500000
FETAL_FRACTION <- 0.1
SIM_COVERAGE <- 4.0
SEED <- 12
READ_LENGTH <- 150 ## NYGC 30X





## other stuff, set at the start
library("QUILT")
library("data.table")
options(scipen = 999)
REGIONNAME <- paste0(CHR, ".", REGIONSTART, ".", REGIONEND)
REGIONSUBSET <- paste0(CHR, ":", REGIONSTART - BUFFER, "-", REGIONEND + BUFFER)
outputdir <- "quilt_output"
outputdir2 <- "quilt_output2"
outputdir3 <- "quilt_output3"
## choose ref panel here
## REF_PREFIX <- "ONEKG"
##REF_PREFIX <- "hrc"
reference_haplotype_file <- paste0(REF_PREFIX, ".hap.gz")
reference_legend_file <- paste0(REF_PREFIX, ".legend.gz")
reference_sample_file <- paste0(REF_PREFIX, ".samples")
##hrc_reference_haplotype_file <- paste0(HRC_REF_PREFIX, ".hap.gz")
##hrc_reference_legend_file <- paste0(HRC_REF_PREFIX, ".legend.gz")
##hrc_reference_sample_file <- paste0(HRC_REF_PREFIX, ".samples")



## choose some temporary directory
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(ANALYSIS_DIR)




library("STITCHprivate")


x <- read.table(reference_sample_file, header = TRUE)

## unclear, might be slow
## not sure what else I can get out of this
## ARGH, re-write
ref_samples <- read.table(reference_sample_file, header = TRUE)
colnames(ref_samples) <- c("SAMPLE", "POP", "GROUP", "SEX")
ref_samples[1:2500, "POP"] <- "ONE" ## just do a subset for now
ref_samples[, "SEX"] <- rep("male", nrow(ref_samples))
temp_reference_sample_file <- tempfile()
write.table(ref_samples, file = temp_reference_sample_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

     
outputdir3 <- "quilt_output3"
STITCHprivate::STITCH(
    chr = CHR,
    bamlist = "bamlist.nipt.txt",
    posfile = "pos.txt",
    phasefile = "phasefile.txt",       
    outputdir = outputdir3,
    K = 10,
    nGen = 100,
    nCores = 1,
    method = "triploid-nipt",
    ffs = rep(FETAL_FRACTION, 1),
    niterations = 1, ## I think OK
    refillIterations = NA,
    shuffleHaplotypeIterations = NA,
    reference_haplotype_file = paste0(REF_PREFIX, ".hap.gz"),
    reference_legend_file = paste0(REF_PREFIX, ".legend.gz"),
    reference_sample_file = temp_reference_sample_file,
    reference_populations = "ONE",
    joint_mat_fetal_output = TRUE,
    regionStart = 12000000,
    regionEnd = 12500000,
    buffer = 100000,
    keepInterimFiles = TRUE
)

## OK seems to work, next time, would take ~10 times longer to generate some results
## but then could hopefully work without massive difficulty
##     reference_iterations = 10
##    genfile = data_package$genfile,
##    genFfile = data_package$genFfile,



## afterwards, load up sample reads
## load up the other stuff
## process through once maybe?
