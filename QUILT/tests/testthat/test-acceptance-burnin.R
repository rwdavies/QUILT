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
K <- 6
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
n_samples_per_pop <- 50
refpack <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = n_samples_per_pop,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster
)
Ktotal <- 2 * 2 * n_samples_per_pop

test_that("QUILT can impute a few samples in a standard way, using a large panel", {

    set.seed(0105)
    
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

    QUILT(
        outputdir = outputdir,
        chr = data_package$chr,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        nGibbsSamples = 3,
        n_seek_its = 4,
        n_burn_in_seek_its = 2,
        nCores = 1,
        Ksubset = 50,
        Knew = 50,
        ## RData_objects_to_save = "super_out_read_labels",
        ## output_RData_filename = "~/temp.RData",
        override_default_params_for_small_ref_panel = FALSE
    )
    ## 

    which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)
    ## now evaluate versus truth!
    check_quilt_output(
        file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
        data_package = data_package,
        which_snps = which_snps,
        tol = 0.1,
        min_info = 0.9
    )
    
    ## now evaluate versus true phase!
    check_sew_phase(
        file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
        phase = data_package$phase,
        which_snps = which_snps
    )

})

