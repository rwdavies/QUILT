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

n_snps <- 200
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
    n_samples_per_pop = 500,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster
)
set.seed(010)




test_that("QUILT can impute a few samples in a standard way using either normal, mspbwt, or zilong approaches", {

    outputdir <- STITCH::make_unique_tempdir()

    regionStart <- 11
    regionEnd <- 200 - 10
    buffer <- 5
    i_method <- 2

    for(i_method in 1:3) {

        for(use_hapMatcherR in c(TRUE, FALSE)) {

            if (i_method == 1) {
                zilong <- FALSE
                use_mspbwt <- TRUE
            } else if (i_method == 2) {
                zilong <- TRUE
                use_mspbwt <- FALSE
            } else {
                zilong <- FALSE
                use_mspbwt <- FALSE
            }

            print(paste0("i_method = ", i_method, ", use_hapMatcherR = ", use_hapMatcherR))

            QUILT(
                outputdir = outputdir,
                chr = data_package$chr,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                bamlist = data_package$bamlist,
                posfile = data_package$posfile,
                genfile = data_package$genfile,
                phasefile = data_package$phasefile,
                reference_vcf_file = refpack$reference_vcf_file,
                reference_haplotype_file = refpack$reference_haplotype_file,
                reference_legend_file = refpack$reference_legend_file,
                genetic_map_file = refpack$reference_genetic_map_file,
                nGibbsSamples = 5,
                n_seek_its = 3,
                nCores = 1,
                nGen = 100,
                use_mspbwt = use_mspbwt,
                zilong = zilong,
                use_hapMatcherR = use_hapMatcherR
            )

            regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
            which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)

            ## now evaluate versus truth!
            check_quilt_output(
                file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
                data_package = data_package,
                which_snps = which_snps,
                tol = 0.1,
                min_info = 0.9
            )

        }

    }

})

