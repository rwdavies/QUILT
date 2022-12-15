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


n_snps <- 50
chr <- 10
K <- 20
set.seed(919)
phasemaster <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
reads_span_n_snps <- 3
n_samples <- 3
## want about 4X here
n_reads <- round(4 * n_snps / reads_span_n_snps)
data_package <- STITCH::make_acceptance_test_data_package(
    reads_span_n_snps = reads_span_n_snps,
    n_samples = n_samples,
    n_snps = n_snps,
    n_reads = n_reads,
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




test_that("QUILT can exclude individuals from the reference panel and then impute", {
    
    outputdir <- STITCH::make_unique_tempdir()
    regionStart <- 11
    regionEnd <- 40
    buffer <- 5
    
    reference_exclude_samplelist_file <- tempfile()
    cat(sort(sample(refpack$reference_samples[, 1], 95)), sep = "\n", file = reference_exclude_samplelist_file)
    
    QUILT_prepare_reference(
        outputdir = outputdir,
        chr = data_package$chr,
        nGen = 100,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        genetic_map_file = refpack$reference_genetic_map_file,
        reference_exclude_samplelist_file = reference_exclude_samplelist_file,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        expRate = 0.5
    )
    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
    expect_true(file.exists(file_quilt_prepared_reference(outputdir, regionName)))
        
    
})

