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
refpack <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 50,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster
)


test_that("QUILT can avoid using small_eHapsCurrent_tc", {

    files <- NULL
    
    for(use_small_eHapsCurrent_tc in c(FALSE, TRUE)) {

        set.seed(101)
        outputdir <- STITCH::make_unique_tempdir()
        regionName <- data_package$chr
        QUILT(
            outputdir = outputdir,
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            nGibbsSamples = 5,
            n_seek_its = 1,
            nGen = 100,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            genetic_map_file = refpack$reference_genetic_map_file,
            use_small_eHapsCurrent_tc = use_small_eHapsCurrent_tc
        )

        which_snps <- NULL
        
        ## now evaluate versus truth!
        check_quilt_output(
            file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
            data_package = data_package,
            which_snps = which_snps,
            tol = 0.1,
            min_info = 0.9
        )
        
        files <- c(files, file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")))
    }

    ## hmm, phasing swap? otherwise similar?
    ## this check infludes INFO so will test posteriors
    a <- read.table(files[1])
    b <- read.table(files[2])
    expect_equal(a[, 8], b[, 8])
    
    ## make sure they are the same
    ##a <- system(paste0("md5sum ", shQuote(files[1]), " | cut -f1 --delimiter=' '"), intern = TRUE)
    ##b <- system(paste0("md5sum ", shQuote(files[2]), " | cut -f1 --delimiter=' '"), intern = TRUE)
    ##expect_equal(a, b)
    
})
