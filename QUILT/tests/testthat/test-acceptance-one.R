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
    n_samples_per_pop = 50,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster
)

test_that("QUILT can impute a few samples in a standard way", {
    
    outputdir <- STITCH::make_unique_tempdir()
    regionStart <- 11
    regionEnd <- 40
    buffer <- 5
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
    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
    expect_true(file.exists(file_quilt_prepared_reference(outputdir, regionName)))
    i <- 1
    
    for(i in 1:2) {

        if (i == 1) {phasefile <- "" }
        if (i == 2) {phasefile <- data_package$phasefile}
        
        QUILT(
            outputdir = outputdir,
            chr = data_package$chr,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,
            phasefile = phasefile,
            Ksubset = 100,
            Knew = 25,
            nGibbsSamples = 5,
            n_seek_its = 2,
            nCores = 1,
            RData_objects_to_save = "final_set_of_results",
            addOptimalHapsToVCF = TRUE
        )
        
    }

    which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)
    
    ## now evaluate versus truth!
    check_quilt_output(
        file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
        data_package = data_package,
        which_snps = which_snps,
        tol = 0.1,
        min_info = 0.9
    )

    ## check loaded stuff
    load(file_quilt_output_RData(outputdir, regionName))
    expect_equal(length(final_set_of_results), 3)
    
})



if (1 == 0) {


    dim(full_alphaHat_t)
    dim(full_betaHat_t)
    dim(full_gamma_t)
    dim(full_transMatRate_t_H)
    dim(rhb_t)
    dim(distinctHapsB)
    dim(distinctHapsIE)
    dim(hapMatcher)
    dim(full_gammaSmall_t)
    full_gammaSmall_cols_to_get
   

    regionStart = NA
    regionEnd = NA
    buffer = NA
    bamlist = ""
    cramlist = ""
    sampleNames_file = ""
    reference = ""
    nCores = 1
    nGibbsSamples = 7
    n_seek_its = 3
    Ksubset = 400
    Knew = 100
    K_top_matches = 5
    heuristic_match_thin = 0.1
    output_filename = NULL
    prepared_reference_filename = ""
    tempdir = NA
    bqFilter = as.integer(17)
    panel_size = NA
    posfile = ""
    genfile = ""
    phasefile = ""
    maxDifferenceBetweenReads = 1e10
    make_plots = FALSE
    verbose = TRUE
    shuffle_bin_radius = 5000
    iSizeUpperLimit = 1e6
    record_read_label_usage = TRUE
    record_interim_dosages = TRUE
    
}
