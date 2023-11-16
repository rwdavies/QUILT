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


if ( 1 == 0 ) {

    curdir <- getwd()
    library("testthat")
    library("QUILT")
    dir <- "~/proj/QUILT/"
    setwd(paste0(dir, "/QUILT/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(curdir)

}


## method / checking todo
## - fix final phasing read label aggregation
## - average depth when run is too high, is there a bug in what I've done here
## - not clear if block gibbs is working optimally, need to check other fetal fraction / haps, and/or walk through examples carefully

## suboptimal things to be fixed in a real version
## - extraction of SNPs from haplotype file done lazilly at multi-allelic SNPs
## - not all reads are sub-sampled in building the new simulated NIPT bam, reads that don't intersect the SNPs are ignored
## - subsetting of reads in simulation doesn't account for single vs paired end properly

## important global variables
HOSTNAME <- "smew"
if (HOSTNAME == "smew") {
    ANALYSIS_DIR <- "/data/smew1/rdavies/nipt_test_2023_09_11/"
} else if (HOSTNAME == "bmrc") {
    ANALYSIS_DIR <- "/well/davies/users/dcc832/nipt_test_2023_09_04/"
}
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



## helper
run <- function(...) {
    system(paste0(...))
}
wget_or_rsync <- function(cur, link) {
    if (!file.exists(cur)) {
        run("wget ", link)
    } else {
        run("rsync -av --progress ", cur, " .")
    }
}
old_path <- function(file) {
    if (HOSTNAME == "bmrc") {
        file.path("/well/davies/users/dcc832/nipt_test_2022_12_06/", file)
    } else {
        file.path("/data/smew1/rdavies/nipt_test_2023_09_11/", file)
    }
}


filter_and_convert_ref_vcf <- function(input_vcf, ref_vcf_prefix) {
    ## subset and convert into .hap .legend file
    run(
        "bcftools view --output-file ", ref_vcf_prefix, ".temp.vcf.gz ",
        "--output-type z --min-alleles 2 --max-alleles 2 --types snps ", input_vcf,
        " ", REGIONSUBSET
    )
    ## ARGH it contains multiple entries per site, force sites to be different, ideally should choose based on MAF but meh
    run(
        "gunzip -c ", ref_vcf_prefix, ".temp.vcf.gz | ",
        "awk '{if(NR==1 || substr($0, 1, -1) == \"#\" || a != $2) {print $0}; a=$2}' | bgzip",
        " > ", ref_vcf_prefix, ".ready.vcf.gz"
    )
    system(paste0("tabix -f ", ref_vcf_prefix, ".ready.vcf.gz"))
    system(paste0("bcftools convert --haplegendsample ", ref_vcf_prefix, " ", ref_vcf_prefix, ".ready.vcf.gz"))
}
    
get_random_samples_and_haps <- function() {
    if (!file.exists("1000G_2504_high_coverage.sequence.index")) {
        run("wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index")
    }
    index <- read.table("1000G_2504_high_coverage.sequence.index")
    set.seed(SEED);
    samples <- sort(sample(which(index[, 24] == "CEU"), 3))
    samples_to_use <- index[samples, 30]
    haps_to_use <- sample(c(1, 2), 3, replace= TRUE)
    return(list(samples_to_use, haps_to_use))
}










##
## download and prepare reference information
##




if (REF_PREFIX == "HRC") {
    ## just use original files where they are
    if (HOSTNAME == "bmrc") {
        HRC_DIR <- file.path("/well/davies/shared/hrc_ega/sinan_2020_06_20")
    } else {
        HRC_DIR <- file.path("/data/smew1/rdavies/sinan_2020_06_20")
    }
    REF_VCF_FILE_TO_USE <- file.path(HRC_DIR, paste0("hg38_liftover_", CHR, ".vcf.gz"))
    REF_SAMPLE_FILE <- file.path(HRC_DIR, paste0("hg38_liftover_", CHR, ".samples"))
    filter_and_convert_ref_vcf(REF_VCF_FILE_TO_USE, REF_PREFIX)
} else {
    ## download 1000G for chr
    oneKG_vcf_name <- paste0("CCDG_14151_B01_GRM_WGS_2020-08-05_", CHR, ".filtered.shapeit2-duohmm-phased.vcf.gz")
    wget_or_rsync(old_path("CEU_omni_recombination_20130507.tar"), "ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/CEU_omni_recombination_20130507.tar")
    if (!file.exists("/well/davies/users/dcc832/nipt_test_2022_12_05/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz")) {
        system(paste0("wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/", oneKG_vcf_name))
        system(paste0("wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/", oneKG_vcf_name, ".tbi"))
    } else {
        run("rsync -av /well/davies/users/dcc832/nipt_test_2022_12_05/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz* .")
    }
    filter_and_convert_ref_vcf(oneKG_vcf_name, REF_PREFIX)
    run("sed -i 's/sample population group sex/SAMPLE POP GROUP SEX/g' ", REF_PREFIX, ".samples")    
}




## ugh, recombination rate
## download build37 and liftOver
wget_or_rsync(old_path("CEU_omni_recombination_20130507.tar"), "ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/CEU_omni_recombination_20130507.tar")
run("tar -xvf CEU_omni_recombination_20130507.tar")

wget_or_rsync(old_path("liftOver"), "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver")
run("chmod +x liftOver")

wget_or_rsync(old_path("hg19ToHg38.over.chain.gz"), "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz")
run("R -f ~/proj/QUILT/scripts/make_b38_recomb_map.R --args ", ANALYSIS_DIR, " CEU ", gsub("chr", "", CHR))


## prepare a region to analyse
QUILT_prepare_reference(
    outputdir = outputdir,
    chr = CHR, 
    nGen = 100,
    reference_vcf_file = paste0(REF_PREFIX, ".ready.vcf.gz"),
    reference_haplotype_file = reference_haplotype_file,
    reference_legend_file = reference_legend_file,
    reference_sample_file = reference_sample_file,
    genetic_map_file = file.path(ANALYSIS_DIR, "CEU", paste0("CEU-", CHR, "-final.b38.txt.gz")),
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER
)















##
## prepare simulated NIPT bam and truth information
##




## choose 3 haps at random from 3 different samples from the SAME population (otherwise it is too easy!)

out <- get_random_samples_and_haps()
samples_to_use <- out[[1]]
haps_to_use <- out[[2]]


## build truth data
load(file.path(outputdir, "RData", paste0("QUILT_prepared_reference.", CHR, ".", REGIONSTART, ".", REGIONEND, ".RData")))
i <- match(samples_to_use, reference_samples[, 1])
which_haps_to_get <- 2 * i - 1 + (haps_to_use - 1) - 1 ## 0-based 
## get haps
truth_haps <- t(STITCH::inflate_fhb_t(rhb_t, haps_to_get = which_haps_to_get, nSNPs = nrow(pos)))
##
phase <- matrix(paste0(truth_haps[, 1], "|", truth_haps[, 2], "|", truth_haps[, 3]), ncol = 1)
colnames(phase) <- "sim"
write.table(phase, "phasefile.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
## do diploid one with first two haps as well
phase <- matrix(paste0(truth_haps[, 1], "|", truth_haps[, 2]), ncol = 1)
colnames(phase) <- "sim"
write.table(phase, "phasefile.diploid.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

## posfile too
write.table(pos, "pos.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

## get bams for these samples
if (HOSTNAME == "bmrc") {
    wget_or_rsync("/well/davies/shared/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/GRCh38_full_analysis_set_plus_decoy_hla.fa", )
} else {
    wget_or_rsync("/data/smew1/rdavies/nipt_test_2023_09_11/GRCh38_full_analysis_set_plus_decoy_hla.fa", "ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/GRCh38_full_analysis_set_plus_decoy_hla.fa")
}
index <- read.table("1000G_2504_high_coverage.sequence.index")

for(sample in samples_to_use) {
    message(paste0("Downloading ", sample, ", ", date()))
    link <- index[match(sample, index[, 30]), 1]
    run("wget ", link, ".crai")
    ## samtools view -b -s ", 1 / (30 / SIM_COVERAGE), " | samtools sort 
    run("samtools view -b -T GRCh38_full_analysis_set_plus_decoy_hla.fa ", link, " ", REGIONSUBSET, " > ", sample, ".", REGIONNAME, ".bam")
    run("samtools index ", sample, ".", REGIONNAME, ".bam")
}

## get reads for these samples
write.table(
    matrix(paste0(samples_to_use, ".", REGIONNAME, ".bam"), ncol = 1),
    file = "bamlist.txt",
    row.names = FALSE,
    col.names = FALSE,
    sep = "",
    quote = FALSE
)

   
## get reads 
STITCH(
    outputdir = outputdir,
    chr = CHR,
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER,
    nGen = 100,
    K = 10,
    bamlist = "bamlist.txt",
    posfile = "pos.txt",
    generateInputOnly = TRUE,
    save_sampleReadsInfo = TRUE
)

## for each sample, get the reads!
source("~/proj/QUILT/QUILT/R/gibbs-nipt.R")
set.seed(SEED) ## do it again

for(i_sample in 1:length(samples_to_use)) {

    sample <- samples_to_use[i_sample]
    load(file.path(outputdir, "input", paste0("sample.",i_sample, ".input.", REGIONNAME, ".RData")))
    load(file.path(outputdir, "input", paste0("sampleReadsInfo.",i_sample, ".input.", REGIONNAME, ".RData")))

    ## get a sample of truth labels
    which_haps_to_get <- 2 *match(samples_to_use[i_sample], reference_samples[, 1]) + -1:0 - 1
    truth_haps_t <- STITCH::inflate_fhb_t(rhb_t, haps_to_get = which_haps_to_get, nSNPs = nrow(pos))
    ## now - determine which come from which haps - and get the appropriate hap
    hap_to_use <- haps_to_use[i_sample]
    readProbs <- assign_fetal_read_probabilities(sampleReads, truth_haps_t)
    read_labels <- 2 - as.integer(runif(length(sampleReads)) < (readProbs[1, ] / colSums(readProbs)))

    
    ## get the reads we want
    cov_wanted <- c(0.5, 0.5 - FETAL_FRACTION / 2, FETAL_FRACTION / 2)[i_sample] * SIM_COVERAGE
    reads_wanted <- ((REGIONEND + BUFFER) - (REGIONSTART - BUFFER)) / READ_LENGTH * cov_wanted
    ## but this is from ALL reads, so downsample as appropriate
    n_reads_in_STITCH_object <- sum(nchar(sampleReadsInfo[, "strand"]))
    total_reads_in_region <- system(paste0("samtools view ", sample, ".",REGIONNAME, ".bam | wc -l"), intern = TRUE)
    reads_wanted <- round(reads_wanted * (n_reads_in_STITCH_object / as.integer(total_reads_in_region)))
    message(paste0("For the ", i_sample, " sample, keep ", reads_wanted, " reads "))
    
    ## sample the reads!
    reads_to_keep <- sample(which(read_labels == hap_to_use), reads_wanted, replace = FALSE)
    file1 <- tempfile()
    file2 <- tempfile()
    file3 <- tempfile()
    cat(sampleReadsInfo[reads_to_keep, "qname"], file = file1, sep = "\n")

    run("samtools view ", sample, ".",REGIONNAME, ".bam > ", file2)
    run("samtools view -H ", sample, ".",REGIONNAME, ".bam > ", file3)
    ## add the reads
    run("awk 'FNR==NR {a[$1]; next} FNR > 1 && $1 in a' ", file1, " ", file2, " >> ", file3)
    ## put into temporary bam file
    run("cat ", file3, " | samtools view -b > ", sample, ".prepped.", REGIONNAME, ".bam")
    
}


## merge the samples together
to_merge <- sapply(1:length(samples_to_use), function(i_sample) {
    sample <- samples_to_use[i_sample]
    paste0(sample, ".prepped.", REGIONNAME, ".bam")
})
##run("samtools cat --output-fmt BAM ", paste0(to_merge, collapse = " "), " | samtools sort > nipt.", REGIONNAME, ".bam")
run("samtools cat ", paste0(to_merge, collapse = " "), " | samtools sort > nipt.", REGIONNAME, ".bam")

## fix header
file <- tempfile()
run("samtools view -H nipt.", REGIONNAME, ".bam | grep -v '^@RG' > ", file)
cat("@RG", "ID:someid", "PL:illumina", "PM:Unknown", "LB:NA12812", "DS:GRCh38", "SM:sim\n", file = file, append = TRUE, sep = "\t")
run("samtools reheader ", file, " nipt.", REGIONNAME, ".bam > nipt.", REGIONNAME, ".temp.bam")
run("mv nipt.", REGIONNAME, ".temp.bam nipt.", REGIONNAME, ".bam")
run("samtools index nipt.", REGIONNAME, ".bam")




















##
## try imputing using QUILT NIPT
## 
##


##
## make slightly different reference panel, remove all related people, as well as the three samples used here
##
## prepare a region to analyse

## exclude ref samples and anyone not in original 2504
##sample_to_keep <- setdiff(reference_samples[, 1], samples_to_use)
##samples_to_exclude <- reference_samples[!reference_samples[, 1] %in% sample_to_keep, 1]
write.table(matrix(samples_to_use, ncol = 1), file = "ref.samples.exclude.specific.txt", row.names =FALSE, col.names = FALSE, sep = "", quote = FALSE)
## 
QUILT_prepare_reference(
    outputdir = outputdir2,
    chr = CHR, 
    nGen = 100,
    reference_vcf_file = paste0(REF_PREFIX, ".ready.vcf.gz"),
    reference_legend_file = reference_legend_file,
    reference_sample_file = reference_sample_file,
    reference_exclude_samplelist_file = "ref.samples.exclude.specific.txt",
    genetic_map_file = file.path(ANALYSIS_DIR, "CEU", paste0("CEU-", CHR, "-final.b38.txt.gz")),
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER,
    impute_rare_common = FALSE,
    rare_af_threshold = 0.01,
    use_mspbwt = TRUE
)


## can either do full or downsampled versions
## not entirely sure 
run("samtools view -b -s 100.1666 nipt.", REGIONNAME, ".bam > nipt.downsampled.", REGIONNAME, ".bam")
run("samtools view -b -s 100.999 nipt.", REGIONNAME, ".bam > nipt.downsampled.", REGIONNAME, ".bam")
run("samtools index nipt.downsampled.", REGIONNAME, ".bam")
cat("nipt.downsampled.", REGIONNAME, ".bam\n", sep = "", file = "bamlist.nipt.txt")

## seems to work
fflist <- "fflist.nipt.txt"
cat(FETAL_FRACTION, "\n", sep = "", file = "fflist.nipt.txt")
outputdir2 <- "quilt_output2"
make_plots <- FALSE
make_plots_block_gibbs <- FALSE
## make_plots = TRUE
## make_plots_block_gibbs = TRUE
nGibbsSamples <- 5
n_seek_its <- 3
QUILT(
    outputdir = outputdir2,
    chr = CHR,
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER,
    nGen = 100,
    bamlist = "bamlist.nipt.txt",
    posfile = "pos.txt",
    prepared_reference_filename = file.path(outputdir2, "RData", paste0("QUILT_prepared_reference.", REGIONNAME, ".RData")),
    method = "nipt",
    fflist = fflist,
    nGibbsSamples = nGibbsSamples,
    phasefile = "phasefile.txt",   
    n_seek_its = n_seek_its,
    make_plots = TRUE,
    make_plots_block_gibbs = TRUE,
    impute_rare_common = FALSE,
    use_mspbwt = TRUE
)
## 

## OK TODO
## - are the plots working?
## - does it seem to get better?
## - 



## have done, merge code bases
## now, try it again, then look at block gibbs


## things to look at
## - fix final phasing read label aggregation (maybe?)
## - really doesn't look like block gibbs is working properly
## - not clear if block gibbs is working optimally, need to check other fetal fraction / haps, and/or walk through examples carefully

## - todo
## 1) merge code, make work
## 2) work on block gibbs

## definitely take a look at blocking, even right at the start
