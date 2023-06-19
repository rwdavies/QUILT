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


set.seed(1911919)


n_snps <- 200
K <- 6
n_big_haps <- 100 ## 1000
chr <- 10


## make a spectrum or rare and common SNPs
phasemaster1 <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
phasemaster2 <- phasemaster1[, sample(1:K, n_big_haps, replace = TRUE)]

## make about 10% common, and do not touch
n_common <- round(n_snps / 10)
common <- sample(1:n_snps, n_common)
## make about 90% rare, and make about 1% freq
phasemaster2[-common, ] <- 0
## make remaining frequency about 1%
phasemaster2[-common, ] <- sample(c(0, 1), (n_snps - n_common) * n_big_haps, prob = c(0.99, 0.01), replace = TRUE)


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
    phasemaster = phasemaster2
)
refpack <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 500,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster2
)
set.seed(010)
chr <- data_package$chr



test_that("QUILT can impute a few samples in a standard way using either normal, mspbwt, or zilong approaches, also using rare common idea", {

    options(warn=2)
    
    regionStart <- 11
    regionEnd <- 200 - 10
    buffer <- 5
    i_method <- 1
    impute_rare_common <- TRUE

    for(impute_rare_common in c(TRUE, FALSE)) {

        ## this is the different methods: zilong, mspbwt, none, etc
        for(i_method in 2:2) {

            if (i_method == 1) {
                zilong <- FALSE
                use_mspbwt <- TRUE
                use_hapMatcherR <- TRUE
            } else if (i_method == 2) {
                zilong <- TRUE
                use_mspbwt <- FALSE
                use_hapMatcherR <- TRUE                    
            } else if (i_method == 3) {
                zilong <- FALSE
                use_mspbwt <- FALSE
                use_hapMatcherR <- TRUE                    
            } else {
                zilong <- FALSE
                use_mspbwt <- FALSE
                use_hapMatcherR <- FALSE
            }

            ## this is whether to do in one go (i_approach = 1), or do prepare reference first (i_approach = 2)
            for(i_approach in 1:1) {

                print("FIX THIS - make second one work")

                print(paste0("impute_rare_common = ", impute_rare_common ,", i_method = ", i_method, ", i_approach = ", i_approach))

                outputdir <- STITCH::make_unique_tempdir()
                
                if (i_approach == 1) {

                    QUILT(
                        outputdir = outputdir,
                        chr = chr,
                        regionStart = regionStart,
                        regionEnd = regionEnd,
                        buffer = buffer,
                        bamlist = data_package$bamlist,
                        posfile = data_package$posfile,
                        genfile = data_package$genfile,
                        phasefile = data_package$phasefile,
                        reference_vcf_file= refpack$reference_vcf_file,
                        reference_haplotype_file = refpack$reference_haplotype_file,
                        reference_legend_file = refpack$reference_legend_file,
                        genetic_map_file = refpack$reference_genetic_map_file,
                        nGibbsSamples = 5,
                        n_seek_its = 3,
                        nCores = 1,
                        nGen = 100,
                        use_mspbwt = use_mspbwt,
                        zilong = zilong,
                        use_hapMatcherR = use_hapMatcherR,
                        impute_rare_common = impute_rare_common,
                        mspbwt_nindices = 1,
                        mspbwtB = 32L
                    )

                } else {

                    QUILT_prepare_reference(
                        outputdir = outputdir,
                        chr = chr,
                        regionStart = regionStart,
                        regionEnd = regionEnd,
                        buffer = buffer,
                        reference_vcf_file= refpack$reference_vcf_file,
                        reference_haplotype_file = refpack$reference_haplotype_file,
                        reference_legend_file = refpack$reference_legend_file,
                        genetic_map_file = refpack$reference_genetic_map_file,
                        nGen = 100,
                        use_mspbwt = use_mspbwt,
                        use_zilong = zilong,
                        use_hapMatcherR = use_hapMatcherR,
                        impute_rare_common = impute_rare_common,
                        mspbwt_nindices = 1,
                        mspbwtB = 32L
                    )

                    QUILT(
                        outputdir = outputdir,
                        chr = chr,
                        regionStart = regionStart,
                        regionEnd = regionEnd,
                        buffer = buffer,
                        bamlist = data_package$bamlist,
                        posfile = data_package$posfile,
                        genfile = data_package$genfile,
                        phasefile = data_package$phasefile,
                        nGibbsSamples = 5,
                        n_seek_its = 3,
                        nCores = 1,
                        use_mspbwt = use_mspbwt,
                        zilong = zilong,
                        impute_rare_common = impute_rare_common                        
                    )

                }
                
                regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
                which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)

                ## a <- file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz"))
                ## print(dim(a))
                ## print(length(which_snps))
                ## print(table(which_snps))
                
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

    }

})




