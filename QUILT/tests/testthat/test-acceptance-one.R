if ( 1 == 0 ) {

    library("testthat")
    library("QUILT")
    dir <- "~/proj/QUILT/"
    setwd(paste0(dir, "/QUILT/R"))
    a <- dir(pattern = "*R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)


}


n_snps <- 50
chr <- 10
K <- 2 
set.seed(919)
phasemaster <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
data_package <- STITCH::make_acceptance_test_data_package(
    reads_span_n_snps = 3,
    n_samples = 3,
    n_snps = n_snps,
    n_reads = n_snps * 10,
    seed = 2,
    chr = chr,
    K = K,
    phasemaster = phasemaster
)
refpack <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 50,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster
)



test_that("SEW can phase a sample at high coverage with long reads", {

    outputdir <- STITCH::make_unique_tempdir()
    regionStart <- 11
    regionEnd <- 40
    buffer <- 5
    ##
    QUILT_prepare_reference(
        outputdir = outputdir,
        chr = data_package$chr,
        nGen = 100,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        genetic_map_file = refpack$reference_genetic_map_file,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )

    ## temporary check
    expect_true(file.exists(quilt_prepared_reference(outputdir, paste0(data_package$chr, ".", regionStart, ".", regionEnd))))

    ## run QUILT
    ## 
    ## library("testthat")
    ## library("QUILT")
    ## dir <- "~/proj/QUILT/"
    ## setwd(paste0(dir, "/QUILT/R"))
    ## a <- dir(pattern = "*R")
    ## b <- grep("~", a)
    ## if (length(b) > 0) {
    ##     a <- a[-b]
    ## }
    ## o <- sapply(a, source)
    ## ## 
    ## QUILT(
    ##     outputdir = outputdir,
    ##     chr = data_package$chr,
    ##     regionStart = regionStart,
    ##     regionEnd = regionEnd,
    ##     buffer = buffer
    ##     bamlist = data_package$bamlist
    ## )

})
