#!/usr/bin/env Rscript

## simple script to trio phase at known sites

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}


curdir <- getwd()
library("data.table")
library("parallel")
source("scripts/phase_functions.R")
options(scipen = 999)
## source("~/proj/QUILT/scripts/phase_functions.R")



##
## IMPORTANT
## 
## the code below calls prepare_for_phase2.R
## which needs the following to work
## . ~/.conda_activate && conda activate shapeit4 && shapeit4
## so either modify prepare_for_phase2.R
## or put conda activate code (from miniconda) into ~/.conda_activate
## and run code like the following
##
## . ~/.conda_activate
## conda create --name shapeit4
## conda install shapeit4=4.0 -c bioconda


    ## ## set parameters here
    ## outputdir <- "/data/smew1/rdavies/ukbb_gel_2023_01_26/phasing_test"
    ## chr <- "chr20"
    ## regionStartMinusBuffer <- 20000000
    ## regionEndPlusBuffer <- 21000000
    ## REF_FA <- "/data/smew1/rdavies/ukbb_gel_2023_01_26/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    ## pop <- "CEU" ## what population to get
    ## n <- 10 ## how many samples to get
    ## invcf <- "/data/smew1/rdavies/ukbb_gel_2023_01_26/2023_06_19A/beagle.20.20000000.21000000.clean.cols.1.10.vcf.gz"
    ## nCores <- 16
    ## rebuild <- FALSE ## turn on to force re-building everything
    ## genetic_map_file <- "/data/smew1/rdavies/ukbb_gel_2023_01_26/CEU/CEU-chr20-final.b38.txt.gz"
    ## refpanel <- "/data/smew1/rdavies/hrc/hg38_liftover_clean_chr20.vcf.gz" ## for trio triple het phasing

args <- commandArgs(trailingOnly = TRUE)

print(args)

outputdir <- args[1]
regionStartMinusBuffer <- as.integer(args[2])
regionEndPlusBuffer <- as.integer(args[3])
pop <- args[4]
n <- as.integer(args[5])
invcf <- args[6]
rebuild <- as.logical(args[7])


chr <- "chr20"
nCores <- 16
genetic_map_file <- "/data/smew1/rdavies/ukbb_gel_2023_01_26/CEU/CEU-chr20-final.b38.txt.gz"
refpanel <- "/data/smew1/rdavies/hrc/hg38_liftover_clean_chr20.vcf.gz" ## for trio triple het phasing
REF_FA <- "/data/smew1/rdavies/ukbb_gel_2023_01_26/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"




## derived
dir.create(outputdir, showWarnings = FALSE)
setwd(outputdir)
region <- paste0("chr20", ":", regionStartMinusBuffer, "-", regionEndPlusBuffer)
dir.create("bams/", showWarnings = FALSE)


## get sample names
out <- get_1000_genomes_samples(pop = pop, n = n)
samples <- out[["samples"]]
data <- out[["data"]]
ped <-  out[["ped"]]
samples_needed_for_phasing <- out[["samples_needed_for_phasing"]]
links_needed_for_phasing <- out[["links_needed_for_phasing"]]    


## should already be made now
tsv <- gsub(".vcf.gz", ".tsv", invcf)
system(paste0("bcftools query -f'chr20\t%POS\t%REF,%ALT\n' ", invcf, " > ", tsv))


## change directory for indexfiles
clean_region <- paste0(gsub("chr", "", chr), ".", regionStartMinusBuffer, ".", regionEndPlusBuffer)


print(paste0("In directory:", getwd()))
##
## download the high coverage bams and call genotypes
## 
out <- mclapply(1:length(samples_needed_for_phasing), mc.cores = nCores, function(i_sample) {
    
    sample <- samples_needed_for_phasing[i_sample]
    ftp <- links_needed_for_phasing[i_sample]
    ftp <- gsub("ftp://", "https://", ftp)
    
    outbam <- file.path("bams", paste0(sample, ".HC.", clean_region, ".bam"))
    outbam_lc <- file.path("bams", paste0(sample, ".1.0.", clean_region, ".bam"))
    
    OUT <- file.path("bams", paste0(sample, ".HC.", clean_region, ".vcf.gz"))
    
    
    ## only rebuild if necessary
    if (rebuild | (!file.exists(OUT) | (!file.exists(outbam_lc)))) {
        system(paste0("samtools view -b -T ", REF_FA, " ", ftp, " ", region, " > ", outbam))
        system(paste0("samtools index ", outbam))
    }
    
    if (!file.exists(outbam_lc) | rebuild) {
        ## about 1X
        cmd <- paste0("samtools view -s 2023.033333 ", outbam, " -o ", outbam_lc)
        system(cmd)
        system(paste0("samtools index ", outbam_lc))
    }
    
    ## get truth genotypes at those sites
    OUT <- file.path("bams", paste0(sample, ".HC.", clean_region, ".vcf.gz"))
    if (!file.exists(OUT) | rebuild) {
        message(paste0("Call variants"))
        system(paste0("bcftools mpileup -f ", REF_FA, " -I -E -a 'FORMAT/DP' -T ", invcf, " -r ", chr, " ", outbam, " -Ou | bcftools call -Aim -C alleles -T ", tsv, " -Oz -o ", OUT))
        system(paste0("bcftools index -f ", OUT))
    }
    
    ## remove high coverage bam files
    if (file.exists(outbam)) {
        unlink(outbam)
    }
    
})


## check for errors
##te <- (sapply(out, class) == "try-error") | sapply(out, is.null)
##if (sum(te) > 0) {
##    print_message(out[[which(te)[1]]]) # print first error
##    stop("Something went wrong")
##}





##
## merge the genotypes, make genotype truth file
## 
specific_posfile <- file.path(outputdir, paste0("posfile.HC.", clean_region, ".txt"))
specific_genfile <- file.path(outputdir, paste0("genfile.HC.", clean_region, ".txt"))

if (rebuild | !file.exists(specific_genfile)) {
    
    tmp <- tempfile()
    cat(paste0("bams/", samples, ".HC.", clean_region, ".vcf.gz"), sep = "\n", file = tmp)
    
    ## make a single truth VCF
    system(paste0("bcftools merge -o bams/all.HC.", clean_region, ".vcf.gz --file-list ", tmp))
    
    ## make a phasefile and a genfile
    out <- load_and_convert_vcf(paste0("bams/all.HC.", clean_region, ".vcf.gz"))
    pos <- out[["pos"]]
    G <- out[["G"]]
    ## add rownames
    rownames(pos) <- paste(chr, pos[, 2], pos[, 3], pos[, 4], sep = "-")
    
    write.table(pos, file = specific_posfile, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    write.table(G, file = specific_genfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    gen2 <- G
    pos2 <- pos
    
}

        






##
## do trio phasing, and then scaffolded phasing of triple hets
## 

specific_phasefile <- file.path(outputdir, paste0("phasefile.HC.", clean_region, ".txt"))

if (rebuild | !file.exists(specific_phasefile)) {
    
    ## make a single VCF
    tmp <- tempfile()
    cat(paste0("bams/", samples_needed_for_phasing, ".HC.", clean_region, ".vcf.gz"), sep = "\n", file = tmp)
    vcffile <- paste0("bams/everyone.HC.", clean_region, ".vcf.gz")
    system(paste0("bcftools merge -o ", vcffile, " --file-list ", tmp))
    
    ## subset to sites of interest
    vcffile <- paste0("bams/everyone.HC.", clean_region, ".vcf.gz")
    vcf <- fread(cmd = paste0("gunzip -c ", vcffile, "| grep -v '##'"), data.table = FALSE)
    
    ## match
    vcf2 <- vcf
    fix <- (rowSums(vcf2[, -c(1:9)] == ".")) == (ncol(vcf2) - 9)
    print(table(fix)) ## 2%? a bit high?
    vcf2[fix, -c(1:9)] <- "./." ## not sure?
    
    ## write
    vcffile_clean <- paste0("bams/specific_clean_everyone.HC.", clean_region, ".vcf.gz")    
    system(paste0("bcftools view -h ", vcffile, " > ", gsub(".gz", "", vcffile_clean)))
    fwrite(vcf2, file = gsub(".gz", "", vcffile_clean), row.names = FALSE, col.names = FALSE, sep = "\t", append = TRUE)
    system(paste0("bgzip -f ", gsub(".gz", "", vcffile_clean)))
    system(paste0("bcftools index ", vcffile_clean))
    
    
    
    ## warning("Possibly clean up genotyping missing data somehow")
    ##
    ## AM HERE
    ## THIS SHOULD WORK
    ## PHASING A BIT SLOW BUT OTHERWISE OK
    ## 1) fix the fact that many of these SNPs I can't genotype properly
    ## 2) can't quite get phasing to work yet! fix that too
    ## SHOULD HAVE THE RIGHT INPUT, MORE OR LESS
    ## RE-CONFIGURE AND RUN USING CODE BELOW
    ## 
    
    
    ## do phasing fancy way
    
    ## ref panel
    chr <- "chr20"
    region <- paste0("chr20", ":", regionStartMinusBuffer, "-", regionEndPlusBuffer)
    outstem <- paste0("bams/everyone.HC.", clean_region)
    
    args <- c(
        vcffile = vcffile_clean,
        outstem = outstem,
        fid = paste0(ped[, 1], collapse = "_"),
        iid = paste0(ped[, 2], collapse = "_"),
        parent1 = paste0(ped[, 3], collapse = "_"),
        parent2 = paste0(ped[, 4], collapse = "_"),
        recomb_file = genetic_map_file,
        region = region,
        refpanel = refpanel,
        nthreads = 16
    )
    
    cmd <- paste0(
        "R -f ", file.path(curdir, "scripts/prepare_for_phase2.R"), " --args ",
        paste0(args, collapse = " "), sep = ""
    )
    system(cmd)
    
    
    ##
    ## now load and fix etc
    ##
    ## check out / chean up / make phasefile
    vcffile <- paste0(outstem, ".phased.vcf.gz")
    phased_vcf <- fread(cmd = paste0("gunzip -c ", vcffile, "| grep -v '##'"), data.table = FALSE)
    
    ## check out scaffold
    vcffile <- paste0(outstem, ".scaffold.vcf.gz")
    scaffold <- fread(cmd = paste0("gunzip -c ", vcffile, "| grep -v '##'"), data.table = FALSE)        
    
    
    ## fill in the blanks at first
    phase2 <- array(NA, c(nrow(gen2), ncol(gen2), 2))
    colnames(phase2) <- colnames(gen2)
    a1 <- paste("chr20", pos2[, 2], pos2[, 3], pos2[, 4], sep = "-")
    a2 <- paste("chr20", phased_vcf[, 2], phased_vcf[, 4], phased_vcf[, 5], sep = "-")
    t1 <- match(a2, a1)
    t2 <- setdiff(1:length(a1), t1)
    stopifnot(sum(is.na(t1)) == 0)
    stopifnot(length(intersect(t1, t2)) == 0)
    for(i in 1:ncol(gen2)) {
        samp <- colnames(gen2)[i]
        h1 <- as.integer(substr(phased_vcf[, samp], 1, 1))
        h2 <- as.integer(substr(phased_vcf[, samp], 3, 3))
        phase2[t1, samp, 1] <- h1
        phase2[t1, samp, 2] <- h2
        ## fill in missing using scaffold
        suppressWarnings(phase2[t2, samp, 1] <- as.integer(substr(scaffold[t2, samp], 1, 1)))
        suppressWarnings(phase2[t2, samp, 2] <- as.integer(substr(scaffold[t2, samp], 3, 3)))
    }
    
    print(paste0("Missing data proportion for phase:", sum(is.na(phase2)) / prod(dim(phase2)) * 100))
    phase2B <- array("", c(nrow(gen2), ncol(gen2)))
    colnames(phase2B) <- colnames(phase2)
    for(i in 1:ncol(phase2)) {
        phase2B[, i] <- paste0(phase2[, i, 1], "|", phase2[, i, 2])
    }
    
    clean_region <- paste0("20.", regionStartMinusBuffer, ".", regionEndPlusBuffer)    
    ## specific_phasefile <- file.path(output_date, paste0("phasefile.HC.", clean_region, ".txt"))
    write.table(phase2B, file = specific_phasefile, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
}
