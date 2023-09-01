simulate_hla_pseudo_acceptance_test_dataset <- function(
    hla_gene,
    ref_dir,
    ref_fasta,
    dict_entries,
    nSamples = 1,
    cov = 1,
    readLength = 150,
    insert_size = 450,
    nCores = 1
) {

    load(file_quilt_hla_specific_haplotypes(ref_dir, hla_gene))
    load(file_quilt_hla_phased_haplotypes(ref_dir, hla_gene))
    start <- pos[1, 2]
    end <- pos[nrow(pos), 2]
    has_a_snp <- rep(FALSE, max(pos[, 2] - start + 1))
    snp_matcher <- pos[, 2] - start + 1
    has_a_snp[snp_matcher] <- TRUE
    snp_count <- cumsum(has_a_snp)
    n <- max(snp_matcher)

    cigar <- paste0(readLength, "M")
    bq <- paste0(rep(":", readLength), collapse = "") ## 25 I think
    sam_head <- paste0(paste0(dict_entries[, 1], "\t", dict_entries[, 2], "\t", dict_entries[, 3]), collapse = "\n")
    
    results <- mclapply(1:nSamples, mc.cores = nCores, function(iSample) {
        
        which_haps <- sample(1:length(hlahaptypes), 2)
        which_alleles <- hlahaptypes[which_haps]
        
        ## 0 vs 1 based
        haps <- inflate_fhb_t(rhb_t, haps_to_get = which_haps - 1, nSNPs = nSNPs)
        hap1 <- rep("A", n)
        hap2 <- rep("A", n)
        hap1[snp_matcher] <- pos[cbind(1:nrow(pos), haps[1, ] + 3)]
        hap2[snp_matcher] <- pos[cbind(1:nrow(pos), haps[2, ] + 3)]        
        
        nPairedReads <- round(
            cov * (end - start) / 
            2 / readLength
        )

        H <- sample(c(1, 2), replace = TRUE, nPairedReads)
        reads <- lapply(1:nPairedReads, function(iRead) {
            ## choose start of first read
            r1s <- round(start + 2 * insert_size + runif(1) * (end - start - 4 * insert_size)) - start
            r1e <- r1s + readLength - 1
            ## 
            r2s <- r1e + insert_size
            r2e <- r2s + readLength - 1
            ## choose haplotype
            if (H[iRead] == 1) {
                seq1 <- paste0(hap1[r1s:r1e], collapse = "")
                seq2 <- paste0(hap1[r2s:r2e], collapse = "")
            } else {
                seq1 <- paste0(hap2[r1s:r1e], collapse = "")
                seq2 <- paste0(hap2[r2s:r2e], collapse = "")
            }
            ## now write out the bits
            r1 <- list(
                paste0("r00", iRead), "0", "chr6", r1s + start - 1, "60",
                cigar, "*", "0", "0",
                seq1, bq
            )
            r2 <- list(
                paste0("r00", iRead), "0", "chr6", r2s + start - 1, "60",
                cigar, "*", "0", "0",
                seq2, bq
            )
            list(r1, r2)
        })
        ## make one read per
        reads <- c(lapply(reads, function(x) x[[1]]), lapply(reads, function(x) x[[2]]))
        reads <- reads[order(sapply(reads, function(x) x[[4]]))]
        ##
        reads <- sapply(reads, function(x) paste0(x, collapse = "\t"))
        reads <- paste0(reads, collapse = "\n")
        
        
        sam <- paste0(
            "@HD\tVN:1.5\tSO:coordinate\n",
            sam_head, "\n",
            "@RG\tID:7369_8x15\tSM:sample", iSample, "\n",
            reads
        )

        ## normal sequencing reads OK
        bam_file <- make_simple_bam(
            file_stem = tempfile(),
            sam = sam
        )
        
        ## for those in the gene itself
        ## put into header
        
        to_return <- list(
            alleles = which_alleles,
            who_copied = which_haps,
            bam_file = bam_file
        )

    })

    ## normalize a bit
    alleles <- t(sapply(results, function(x) x$alleles))
    who_copied <- t(sapply(results, function(x) x$who_copied))
    bam_files <- sapply(results, function(x) x$bam_file)
    bamlist <- tempfile()
    write.table(
        matrix(bam_files, ncol = 1),
        file = bamlist,
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
    )

    list(
        alleles = alleles,
        who_copied = who_copied,
        bam_files = bam_files,
        bamlist = bamlist
    )

}


if ( 1 == 0 ) {
    
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


}


## can start using real data (run test if it exists)
## then can simulate what haplotypes to carry
## and reads (somehow)
## then can push through code to check

test_that("can run acceptance test on HLA", {

    skip("wip")

    
    ## this is quite a big file!
    ref_dir <- "/well/davies/users/dcc832/quilt_hla_packages/quilt_hla_reference_panel_files_2021_12_28_benchmarking_3.43/"

    
    hla_gene <- "A"
    ref_fasta <- "/well/davies/users/dcc832/single_imp/2020_10_25/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    
    if (!file.exists(file_quilt_hla_specific_haplotypes(ref_dir, hla_gene))) {
        skip("input data not available to test QUILT HLA")
    }

    if (!file.exists(ref_fasta)) {
        skip("input data not available to test QUILT HLA")
    }

    outputdir <- tempdir()
    dict_entries <- read.table("~/proj/QUILT/hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict", sep = "\t", comment.char = "", skip = 1)

    ## source("~/proj/QUILT/QUILT/tests/testthat/test-acceptance-hla.R")
    data_package <- simulate_hla_pseudo_acceptance_test_dataset(
        hla_gene = hla_gene,
        ref_dir = ref_dir,
        ref_fasta = ref_fasta,
        dict_entries = dict_entries,
        nSamples = 1,
        cov = 1,
        readLength = 150,
        insert_size = 450,
        nCores = 1
    )

    bamlist <- data_package$bamlist

    ## UGH am here
    ## some stupid error somehow!
    ## really would benefit from running faster!
    
    ## 
    QUILT_HLA(
        outputdir = outputdir,
        bamlist = bamlist,
        region = hla_gene,
        prepared_hla_reference_dir = ref_dir,
        quilt_hla_haplotype_panelfile = file.path(ref_dir, paste0("quilt.hrc.hla.", hla_gene, ".haplotypes.RData")),
        dict_file = file.path("~/proj/QUILT/hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict")
    )
    ## 
    
    ## check that the inferred type is correct

    
    
})
