ANALYSIS_DIR <- "/home/zilong/.tmp/quilt_nipt_test"
REF_PREFIX <- "HRC"
CHR <- "chr20"
REGIONSTART <- 10082590
REGIONEND <- 15082589
BUFFER <- 500000
FETAL_FRACTION <- 0.2
SIM_COVERAGE <- 2.0
SEED <- 12
READ_LENGTH <- 150 ## NYGC 30X


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
reference_vcf_file <- paste0(REF_PREFIX, ".vcf.gz")
reference_haplotype_file <- paste0(REF_PREFIX, ".hap.gz")
reference_legend_file <- paste0(REF_PREFIX, ".legend.gz")
reference_sample_file <- paste0(REF_PREFIX, ".samples")
genetic_map_file <- "/home/zilong/.tmp/datahub/genetic-map/CEU-chr20-final.b38.txt.gz"
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(ANALYSIS_DIR)

QUILT_prepare_reference(
    outputdir = outputdir,
    chr = CHR, 
    nGen = 100,
    reference_vcf_file = reference_vcf_file,
    genetic_map_file = genetic_map_file,
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER
)

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

## helper
run <- function(...) {
    system(paste0(...))
}

run2 <- function(...) {
    system(paste0(...), intern = T)
}
##
## prepare simulated NIPT bam and truth information
##

## choose 3 haps at random from 3 different samples from the SAME population (otherwise it is too easy!)
## out <- get_random_samples_and_haps()
## samples_to_use <- out[[1]]
## haps_to_use <- out[[2]]
samples_to_use <- c("NA12878", "NA12889", "NA12890")
haps_to_use <- sample(c(1, 2), 3, replace= TRUE)

## build truth data
load(file.path(outputdir, "RData", paste0("QUILT_prepared_reference.", CHR, ".", REGIONSTART, ".", REGIONEND, ".RData")))
reference_samples <- run2("bcftools query -l ", reference_vcf_file)

i <- match(samples_to_use, reference_samples)
samples_to_use <- reference_samples[i]
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
source("~/Projects/github/QUILT/QUILT/R/gibbs-nipt.R")
set.seed(SEED) ## do it again

for(i_sample in 1:length(samples_to_use)) {
    sample <- samples_to_use[i_sample]
    load(file.path(outputdir, "input", paste0("sample.",i_sample, ".input.", REGIONNAME, ".RData")))
    load(file.path(outputdir, "input", paste0("sampleReadsInfo.",i_sample, ".input.", REGIONNAME, ".RData")))
    ## get a sample of truth labels
    which_haps_to_get <- 2 *match(samples_to_use[i_sample], reference_samples) + -1:0 - 1
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

## can either do full or downsampled versions
## not entirely sure 
run("samtools view -b -s 100.1666 nipt.", REGIONNAME, ".bam > nipt.downsampled.", REGIONNAME, ".bam")
run("samtools view -b -s 100.999 nipt.", REGIONNAME, ".bam > nipt.downsampled.", REGIONNAME, ".bam")
run("samtools index nipt.downsampled.", REGIONNAME, ".bam")
cat("nipt.downsampled.", REGIONNAME, ".bam\n", sep = "", file = "bamlist.nipt.txt")
fflist <- "fflist.nipt.txt"
cat(FETAL_FRACTION, "\n", sep = "", file = "fflist.nipt.txt")

## exclude ref samples and anyone not in original 2504
##sample_to_keep <- setdiff(reference_samples[, 1], samples_to_use)
##samples_to_exclude <- reference_samples[!reference_samples[, 1] %in% sample_to_keep, 1]
write.table(matrix(samples_to_use, ncol = 1), file = "ref.samples.exclude.specific.txt", row.names =FALSE, col.names = FALSE, sep = "", quote = FALSE)

##
## reference_vcf_file_exclude_target <- paste0(REF_PREFIX, ".uniq.vcf.gz")
## run("bcftools view -S ^ref.samples.exclude.specific.txt ", reference_vcf_file, " -Oz -o ", reference_vcf_file_exclude_target)
## run("bcftools index -f ", reference_vcf_file_exclude_target)


QUILT_prepare_reference(
    outputdir = outputdir2,
    chr = CHR, 
    nGen = 100,
    reference_vcf_file = reference_vcf_file,
    reference_exclude_samplelist_file = "ref.samples.exclude.specific.txt",
    genetic_map_file = genetic_map_file,
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER,
    use_mspbwt = TRUE,
    impute_rare_common = TRUE,
    rare_af_threshold = 0.001
)

## seems to work
make_plots <- FALSE
make_plots_block_gibbs <- FALSE
nGibbsSamples <- 5
n_seek_its <- 3

QUILT(
    outputdir = outputdir2,
    prepared_reference_filename = file.path(outputdir2, "RData", paste0("QUILT_prepared_reference.", REGIONNAME, ".RData")),
    chr = CHR,
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER,
    nGen = 100,
    method = "nipt",
    bamlist = "bamlist.nipt.txt",
    posfile = "pos.txt",
    fflist = "fflist.nipt.txt",
    phasefile = "phasefile.txt",   
    nGibbsSamples = nGibbsSamples,
    n_seek_its = n_seek_its,
    use_mspbwt = TRUE,
    mspbwtL = 5,
    mspbwtM = 1,
    impute_rare_common = FALSE,
    rare_af_threshold = 0.001,
    make_plots = make_plots,
    make_plots_block_gibbs = make_plots_block_gibbs
)

# impute_rare_common = FALSE, use_mspbwt = FALSE
## [2024-02-28 13:27:12] Final imputation accuracy for sample sim mat r2:0.975, PSE:0.2%, disc:1.4%
## [2024-02-28 13:27:12] Final imputation accuracy for sample sim fet r2:0.842, PSE:0.2%, disc:15.7%
## [2024-02-28 13:27:12] Final imputation accuracy for sample sim hap r2 mt:0.935, mu:0.93, pt:0.616

# impute_rare_common = TRUE, rare_af_threshold = 0.001, use_mspbwt = FALSE
## [2024-02-28 13:33:05] Final imputation accuracy for sample sim mat r2:0.979, PSE:0.2%, disc:2% (all SNPs)
## [2024-02-28 13:33:05] Final imputation accuracy for sample sim fet r2:0.256, PSE:0.2%, disc:11.3% (all SNPs)
## [2024-02-28 13:33:05] Final imputation accuracy for sample sim hap r2 mt:0.762, mu:0.754, pt:0.556 (all SNPs)

## [2024-02-28 16:25:44] Final imputation accuracy for sample sim mat r2:0.894, PSE:2.4%, disc:17.3% (all SNPs)
## [2024-02-28 16:25:44] Final imputation accuracy for sample sim fet r2:0.286, PSE:1%, disc:26.1% (all SNPs)
## [2024-02-28 16:25:44] Final imputation accuracy for sample sim hap r2 mt:0.32, mu:0.224, pt:0.308 (all SNPs)

QUILT_prepare_reference(
    outputdir = outputdir3,
    chr = CHR, 
    nGen = 100,
    reference_vcf_file = reference_vcf_file,
    reference_exclude_samplelist_file = "ref.samples.exclude.specific.txt",
    genetic_map_file = genetic_map_file,
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER,
    use_mspbwt = TRUE,
    mspbwt_nindices = 4,
    impute_rare_common = TRUE,
    rare_af_threshold = 0.01
)


make_plots <- FALSE
make_plots_block_gibbs <- FALSE
nGibbsSamples <- 7
n_seek_its <- 3

QUILT(
    outputdir = outputdir3,
    prepared_reference_filename = file.path(outputdir3, "RData", paste0("QUILT_prepared_reference.", REGIONNAME, ".RData")),
    chr = CHR,
    regionStart= REGIONSTART,
    regionEnd = REGIONEND,
    buffer = BUFFER,
    nGen = 100,
    bamlist = "bamlist.nipt.txt",
    posfile = "pos.txt",
    method = "nipt",
    fflist = "fflist.nipt.txt",
    phasefile = "phasefile.txt",   
    nGibbsSamples = nGibbsSamples,
    n_seek_its = n_seek_its,
    use_mspbwt = TRUE,
    mspbwtM = 1,
    mspbwtL = 10,
    make_plots = make_plots,
    make_plots_block_gibbs = make_plots_block_gibbs
)


## investigate double_list_of_starting_read_labels 
ff <- 0.2
nReads <- 100
H <- sample(c(1, 2, 3), prob = c(0.5, 1 - ff / 2, ff / 2), nReads, replace = TRUE)
H
