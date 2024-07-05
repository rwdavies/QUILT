if ( 1 == 0 ) {
    
    library("testthat")
    library("QUILT")
    curdir <- getwd()
    dir <- "~/proj/QUILT/"
    setwd(paste0(dir, "/QUILT/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(curdir)


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

            a <- grep("Final imputation accuracy", output)
            if (length(a) > 0) {
                n_imp_dos <- length(grep("r2", output[a]))
                n_phase <- length(grep("PSE", output[a]))
            } else {
                n_imp_dos <- 0
                n_phase <- 0
            }
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




test_that("QUILT can impute in one step", {
    
    outputdir <- STITCH::make_unique_tempdir()

    regionName <- data_package$chr
    QUILT(
        outputdir = outputdir,
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        nGibbsSamples = 3,
        n_seek_its = 1,
        nGen = 100,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        genetic_map_file = refpack$reference_genetic_map_file
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


test_that("QUILT can use specified or inferred nMaxDH", {

    set.seed(3210)
    
    for(nMaxDH in c(NA, 2 ** 4 - 1, 2 ** 8 - 1, 34)) {
        
        outputdir <- STITCH::make_unique_tempdir()
        regionName <- data_package$chr
        QUILT(
            outputdir = outputdir,
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            nGibbsSamples = 3,
            n_seek_its = 1,
            nGen = 100,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            genetic_map_file = refpack$reference_genetic_map_file,
            nMaxDH = nMaxDH
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
    }
    
})


test_that("QUILT can save or not save prepared reference data", {

    for(i in 1:3) {
        if (i == 1) {
            prepared_reference_filename <- ""
            save_prepared_reference <- FALSE            
        } else if (i == 2) {
            prepared_reference_filename <- ""
            save_prepared_reference <- TRUE
        } else if (i == 3) {
            prepared_reference_filename <- tempfile()
            save_prepared_reference <- TRUE            
        }

        outputdir <- STITCH::make_unique_tempdir()
        regionName <- data_package$chr
        
        QUILT(
            outputdir = outputdir,
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            nGibbsSamples = 3,
            n_seek_its = 1,
            nGen = 100,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            genetic_map_file = refpack$reference_genetic_map_file,
            save_prepared_reference = save_prepared_reference,
            prepared_reference_filename = prepared_reference_filename
        )

        if (save_prepared_reference) {
            if (prepared_reference_filename == "") {
                expect_true(file.exists(file_quilt_prepared_reference(outputdir, regionName)))
            } else {
                expect_true(file.exists(prepared_reference_filename))
            }
        } else {
            expect_false(file.exists(file_quilt_prepared_reference(outputdir, regionName)))
        }


        which_snps <- NULL
        
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


test_that("QUILT can impute samples with very few reads", {

    for(n_reads in c(0, 1, 2, 4)) {

        n_snps <- 50
        n_samples <- 2
        data_package_local <- STITCH::make_acceptance_test_data_package(
            reads_span_n_snps = reads_span_n_snps,
            n_samples = n_samples,
            n_snps = n_snps,
            n_reads = n_reads,
            seed = 2,
            chr = chr,
            K = K,
            phasemaster = phasemaster
        )

        if (n_reads == 0) {
            ## change MQ to 0, easier than fixing fixing test drivers...
            j <- 1
            inbam <- data_package_local$bam_files[j]
            system(paste0("mv ", inbam, " ", inbam, ".temp"))
            cmd <- paste0(
                "samtools view -h ", inbam, ".temp | ",
                "sed ", shQuote("s/60/0/g"), "|",
                "samtools view -b -o ", inbam
            )
            system(cmd)
            system(paste0("samtools index ", inbam))
        }

        for(addOptimalHapsToVCF in c(FALSE, TRUE)) {

            outputdir <- STITCH::make_unique_tempdir()
            regionName <- data_package_local$chr
            QUILT(
                outputdir = outputdir,
                chr = data_package_local$chr,
                bamlist = data_package_local$bamlist,
                posfile = data_package_local$posfile,
                nGibbsSamples = 3,
                n_seek_its = 1,
                nGen = 100,
                reference_haplotype_file = refpack$reference_haplotype_file,
                reference_legend_file = refpack$reference_legend_file,
                genetic_map_file = refpack$reference_genetic_map_file,
                addOptimalHapsToVCF = addOptimalHapsToVCF,
                phasefile = data_package_local$phasefile
            )
            which_snps <- NULL
            
            ## now evaluate versus truth!
            check_quilt_output(
                file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
                data_package = data_package,
                which_snps = which_snps,
                tol = 1.1,
                min_info = 0,
                max_missingness = 1.1
            )
            
        }
    }
    
})


test_that("QUILT can use or not use eigen to impute", {

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
        reference_sample_file = refpack$reference_sample_file,
        genetic_map_file = refpack$reference_genetic_map_file,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        expRate = 0.5
    )
    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
    expect_true(file.exists(file_quilt_prepared_reference(outputdir, regionName)))

    for(use_eigen in c(FALSE, TRUE)) {
        
        QUILT(
            outputdir = outputdir,
            chr = data_package$chr,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            use_eigen = use_eigen
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


