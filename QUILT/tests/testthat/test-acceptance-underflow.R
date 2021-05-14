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

test_that("QUILT can adjust to underflow problems", {
    
    outputdir <- STITCH::make_unique_tempdir()

    reads_span_n_snps <- 20
    n_snps <- 200
    n_reads <- round(100 * n_snps / reads_span_n_snps)
    K <- 6 ## make em different
    set.seed(465)
    phasemaster <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
    chr <- "chr20"
    data_package_hc <- STITCH::make_acceptance_test_data_package(
        reads_span_n_snps = reads_span_n_snps,
        n_samples = 1,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 2,
        chr = chr,
        K = K,
        phasemaster = phasemaster,
        phred_bq = 25
    )

    ## deliberately use different haplotypes
    ## so results are unlikely to drive down probabilities
    phasemaster2 <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
    refpack <- STITCH::make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = 50,
        reference_populations = c("CEU", "GBR"),
        chr = chr,
        phasemaster = phasemaster2
    )
    
    regionName <- data_package_hc$chr
    QUILT(
        outputdir = outputdir,
        chr = data_package_hc$chr,
        bamlist = data_package_hc$bamlist,
        posfile = data_package_hc$posfile,
        nGibbsSamples = 3,
        n_seek_its = 1,
        nGen = 100,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        genetic_map_file = refpack$reference_genetic_map_file,
        downsampleToCov = 1000,        
        maxDifferenceBetweenReads = 1e10
    )
        
    expect_true(file.exists(file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz"))))
    
})
