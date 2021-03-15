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

test_that("QUILT can impute a few samples in a standard way, using a small panel, with or without a genetic map", {
    
    outputdir <- STITCH::make_unique_tempdir()
    regionStart <- 11
    regionEnd <- 40
    buffer <- 5

    for(genetic_map_file in c("", refpack$reference_genetic_map_file)) {
        
        QUILT_prepare_reference(
            outputdir = outputdir,
            chr = data_package$chr,
            nGen = 100,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            reference_sample_file = refpack$reference_sample_file,
            genetic_map_file = genetic_map_file,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            expRate = 0.5
        )
        regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
        expect_true(file.exists(file_quilt_prepared_reference(outputdir, regionName)))
        
        QUILT(
            outputdir = outputdir,
            chr = data_package$chr,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile
        )
        
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
    
})


test_that("QUILT can impute entire chromosomes", {
    
    outputdir <- STITCH::make_unique_tempdir()
    QUILT_prepare_reference(
        outputdir = outputdir,
        chr = data_package$chr,
        nGen = 100,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        genetic_map_file = refpack$reference_genetic_map_file
    )
    regionName <- data_package$chr
    expect_true(file.exists(file_quilt_prepared_reference(outputdir, regionName)))
    
    QUILT(
        outputdir = outputdir,
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        nGibbsSamples = 3,
        n_seek_its = 1
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
    
})




test_that("QUILT can use all combinations of posfile, genfile and phasfile in the right way, with or without buffers", {

    verbose <- FALSE
    i_buffer <- 1
    
    for(i_buffer in 1:2) {

        if (verbose) {
            print(paste0("i_buffer = ", i_buffer))
        }
        
        outputdir <- STITCH::make_unique_tempdir()
        if (i_buffer == 1) {
            regionStart <- 11
            regionEnd <- 40
            buffer <- 5
        } else {
            regionStart <- NA
            regionEnd <- NA
            buffer <- NA
        }
        
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
        if (is.na(regionStart)) {
            regionName <- data_package$chr
        } else {
            regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
        }

        expect_true(file.exists(file_quilt_prepared_reference(outputdir, regionName)))
        i <- 1
        
        for(i in 1:4) {

            if (verbose) {
                print(paste0("i = ", i))            
            }
            
            posfile <- ""
            genfile <- ""
            phasefile <- ""
            if (i == 1) {
                ## keep them blank
            } else if (i == 2) {
                posfile <- data_package$posfile
                genfile <- data_package$genfile
            } else if (i == 3) {
                posfile <- data_package$posfile
                phasefile <- data_package$phasefile
            } else if (i == 4) {
                posfile <- data_package$posfile            
                genfile <- data_package$genfile
                phasefile <- data_package$phasefile
            }
            ## capture messages from this to check it worked
            output <- testthat::capture_messages(
                QUILT(
                    outputdir = outputdir,
                    chr = data_package$chr,
                    regionStart = regionStart,
                    regionEnd = regionEnd,
                    buffer = buffer,
                    bamlist = data_package$bamlist,
                    posfile = posfile,
                    genfile = genfile,
                    phasefile = phasefile,
                    nGibbsSamples = 3
                )
            )
            n_imp_dos <- length(grep("Final imputation dosage", output))
            n_phase <- length(grep("Final phasing accuracy", output))            
            if (verbose) {
                print(paste0("n_imp_dos = ", n_imp_dos))
                print(paste0("n_phase = ", n_phase))
            }
            if (i == 1) {
                expect_equal(n_imp_dos, 0)
                expect_equal(n_phase, 0)
            } else if (i == 2) {
                expect_equal(n_imp_dos, n_samples)
                expect_equal(n_phase, 0)
            } else if (i == 3) {
                expect_equal(n_imp_dos, n_samples)
                expect_equal(n_phase, n_samples)
            } else if (i == 4) {
                expect_equal(n_imp_dos, n_samples)
                expect_equal(n_phase, n_samples)
            }
            
            if (!is.na(regionStart)) {
                which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)
            } else {
                which_snps <- NULL
            }
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
