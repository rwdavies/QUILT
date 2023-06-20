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

## n_snps <- 500
n_snps <- 50
chr <- 10
K <- 6
set.seed(919)
phasemaster <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
reads_span_n_snps <- 3
## want about 4X here
n_reads <- round(4 * n_snps / reads_span_n_snps)
data_package <- STITCH::make_acceptance_test_data_package(
    reads_span_n_snps = reads_span_n_snps,
    n_samples = 3,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 2,
    chr = chr,
    K = K,
    phasemaster = phasemaster
)
refpack <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 100,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster
)



test_that("QUILT can remove samples from the reference", {

    set.seed(0105)

    ## make file with samples to exclude
    reference_samples <- refpack$reference_samples[, 1]
    refpack$reference_exclude_samplelist_file <- tempfile()
    refpack$reference_samples_to_exclude <- reference_samples[sample(length(reference_samples), 10)]
    refpack$reference_samples_to_keep <- setdiff(reference_samples, refpack$reference_samples_to_exclude)
    write.table(
        matrix(refpack$reference_samples_to_exclude, ncol = 1),
        file = refpack$reference_exclude_samplelist_file,
        row.names = FALSE,
        col.names = FALSE,
        sep = "",
        quote = FALSE
    )
    rm(reference_samples)
    
    regionStart <- 11
    regionEnd <- 40
    buffer <- 5
    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
    i <- 1

    for(use_reference_vcf in c(FALSE, TRUE)) {

        if (use_reference_vcf) {
            reference_vcf_file <- refpack$reference_vcf_file
            reference_haplotype_file <- ""
            reference_legend_file <- ""
        } else {
            reference_vcf_file <- ""
            reference_haplotype_file <- refpack$reference_haplotype_file
            reference_legend_file <- refpack$reference_legend_file
        }
        
        for(selection_method in c("sample", "pop")) {
            
            if (selection_method == "sample") {
                reference_exclude_samplelist_file <- refpack$reference_exclude_samplelist_file
                reference_populations <- NA
            } else {
                reference_exclude_samplelist_file <- ""
                reference_populations <- "CEU"
            }
            
            outputdir <- STITCH::make_unique_tempdir()

            ## print(paste0("use_reference_vcf = ", use_reference_vcf, ", selection_method = ", selection_method))
            
            QUILT_prepare_reference(
                outputdir = outputdir,
                chr = data_package$chr,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                reference_vcf_file = reference_vcf_file,
                reference_haplotype_file = reference_haplotype_file,
                reference_legend_file = reference_legend_file,
                reference_sample_file = refpack$reference_sample_file,
                reference_exclude_samplelist_file = reference_exclude_samplelist_file,
                reference_populations = reference_populations,
                nGen = 100
            )
            
            ## check prepared stuff
            load(file_quilt_prepared_reference(outputdir, regionName))

            if (selection_method == "sample") {
                expect_equal(length(refpack$reference_samples_to_keep) * 2, nrow(hapMatcherR))
                expect_equal(reference_samples[, 1], refpack$reference_samples_to_keep)
            } else {
                w <- refpack$reference_samples[, 2] == "CEU"                
                expect_equal(reference_samples[, "POP"], rep("CEU", sum(w)))
                expect_equal(nrow(hapMatcherR), sum(w) * 2)
                expect_equal(reference_samples[, 1], refpack$reference_samples[w, 1])
            }

        }
    }

})

