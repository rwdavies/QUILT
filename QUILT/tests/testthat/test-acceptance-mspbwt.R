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
    n_samples = 4,
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



test_that("QUILT can impute a few samples in a standard way using either normal, or mspbwt approaches, also using rare common idea", {

    ## options(warn=2)

    regionStart <- 11
    regionEnd <- 200 - 10
    buffer <- 5
    i_method <- 4
    impute_rare_common <- TRUE
    nCores <- 1
    rare_af_threshold <- 0.01

    for(impute_rare_common in c(TRUE, FALSE)) {

        ## this is the different methods: zilong, mspbwt A, none, mspbwt B
        for(i_method in 1:4) {

            if (i_method == 1) {
                zilong <- FALSE
                use_mspbwt <- TRUE
                use_hapMatcherR <- TRUE
                heuristic_approach <- "A"
                use_list_of_columns_of_A <- TRUE                
            } else if (i_method == 2) {
                zilong <- TRUE
                use_mspbwt <- FALSE
                use_hapMatcherR <- TRUE
                heuristic_approach <- "A"
                use_list_of_columns_of_A <- TRUE                
            } else if (i_method == 3) {
                zilong <- FALSE
                use_mspbwt <- FALSE
                use_hapMatcherR <- TRUE
                heuristic_approach <- "A"
                use_list_of_columns_of_A <- TRUE
            } else {
                zilong <- FALSE
                use_mspbwt <- TRUE
                use_hapMatcherR <- TRUE
                heuristic_approach <- "B"
                use_list_of_columns_of_A <- FALSE
            }

            ## this is whether to do in one go (i_approach = 1), or do prepare reference first (i_approach = 2)
            for(i_approach in 1:2) {

                ## print(paste0("impute_rare_common = ", impute_rare_common ,", i_method = ", i_method, ", i_approach = ", i_approach, ", ", date()))
                set.seed(19)

                outputdir <- STITCH::make_unique_tempdir()
                
                if (i_approach == 1) {

                    QUILT(
                        outputdir = outputdir,
                        chr = chr,
                        regionStart = regionStart,
                        regionEnd = regionEnd,
                        buffer = buffer,
                        bamlist = data_package$bamlist,
                        reference_vcf_file= refpack$reference_vcf_file,
                        reference_haplotype_file = refpack$reference_haplotype_file,
                        reference_legend_file = refpack$reference_legend_file,
                        genetic_map_file = refpack$reference_genetic_map_file,
                        nGibbsSamples = 7,
                        n_seek_its = 2,
                        nCores = nCores,
                        nGen = 100,
                        use_mspbwt = use_mspbwt,
                        use_hapMatcherR = use_hapMatcherR,
                        impute_rare_common = impute_rare_common,
                        mspbwt_nindices = 1,
                        rare_af_threshold = rare_af_threshold,
                        heuristic_approach = heuristic_approach,
                        use_list_of_columns_of_A = use_list_of_columns_of_A
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
                        impute_rare_common = impute_rare_common,
                        use_mspbwt = use_mspbwt,
                        use_hapMatcherR = use_hapMatcherR,
                        mspbwt_nindices = 1,
                        rare_af_threshold = rare_af_threshold,
                        use_list_of_columns_of_A = use_list_of_columns_of_A
                    )

                    QUILT(
                        outputdir = outputdir,
                        chr = chr,
                        regionStart = regionStart,
                        regionEnd = regionEnd,
                        buffer = buffer,
                        bamlist = data_package$bamlist,
                        nGibbsSamples = 7,
                        n_seek_its = 2,
                        nCores = nCores,
                        use_hapMatcherR = use_hapMatcherR,                        
                        use_mspbwt = use_mspbwt,
                        impute_rare_common = impute_rare_common,
                        heuristic_approach = heuristic_approach                        
                    )
                    ## posfile = data_package$posfile,
                    ## genfile = data_package$genfile,
                    ## phasefile = data_package$phasefile,
                    

                }
                
                regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
                which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)

                ## a <- read.table(file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")))
                ## print(dim(a))
                ## print(length(which_snps))
                ## print(table(which_snps))
                
                ## now evaluate versus truth!
                check_quilt_output(
                    file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
                    data_package = data_package,
                    which_snps = which_snps,
                    tol = 0.1,
                    min_info = 0.01,
                    check_info_only = TRUE
                )
                ## surprisingly low info for some of these, occurs naturally, for some less confident ones, given low sample size

            }

        }

    }

})




