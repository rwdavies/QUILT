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



test_that("can avoid using eHapsCurrent_tc in genProbs calculation", {

    ## OK am here
    ## want to test this function
    ## using either eHapsCurrent_tc, or other approach used for eMatRead_t
    ## should be more straightforward? depending on offset?

    K <- 100
    nMaxDH <- 80
    Ksmall <- 10 ## i.e. the subset
    nSNPs <- 200
    nReads <- 20
    Jmax <- 1000
    S <- 1
    ref_error <- 0.01
    maxDifferenceBetweenReads <- 1e100
    rescale_eMatRead_t <- TRUE
    snp_start_1_based <- 1
    snp_end_1_based <- nSNPs
    ## make some test data
    which_haps_to_use <- sort(sample(K, Ksmall)) ## 1-based
    out <- make_fb_test_package(
        K = K,
        nReads = nReads,
        nSNPs = nSNPs,
        gridWindowSize = 32
    )
    grid <- out$grid
    nGrids <- out$nGrids
    rhi <- t(round(out$eHapsCurrent_tc[, , 1]))
    rhi_t <- t(rhi)
    rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)
    
    ##
    f <- function() {
        gammaMT_t <- array(runif(K * nGrids), c(Ksmall, nGrids))
        for(iGrid in 1:nGrids) {
            gammaMT_t[, iGrid] <- gammaMT_t[, iGrid] / sum(gammaMT_t[, iGrid])
        }
        gammaMT_t
    }
    gammaMT_t <- f()
    gammaMU_t <- f()
    gammaP_t <- f()
    ## this is on subse
    small_eHapsCurrent_tc <- array(NA, c(Ksmall, nSNPs, 1))
    small_eHapsCurrent_tc[, , 1] <- rhi_t[which_haps_to_use, ]
    small_eHapsCurrent_tc[small_eHapsCurrent_tc == 1] <- 1 - ref_error
    small_eHapsCurrent_tc[small_eHapsCurrent_tc == 0] <- ref_error    


    ##
    ## original version, using eHapsCurrent_tc
    ##
    genProbsM_t <- array(0, c(3, nSNPs))
    genProbsF_t <- array(0, c(3, nSNPs))
    hapProbs_t <- array(0, c(3, nSNPs))
    
    ## ame here, write this first 
    ## 
    rcpp_calculate_gn_genProbs_and_hapProbs(
        genProbsM_t = genProbsM_t,
        genProbsF_t = genProbsF_t,
        hapProbs_t = hapProbs_t,
        s = 0,
        eHapsCurrent_tc = small_eHapsCurrent_tc,
        gammaMT_t = gammaMT_t,
        gammaMU_t = gammaMU_t,
        gammaP_t = gammaP_t,
        grid = grid,
        snp_start_1_based = snp_start_1_based,
        snp_end_1_based = snp_end_1_based,
        grid_offset = 0
    )



    
    ##
    ## new version
    ##
    out <- make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nSNPs,
        ref_error = ref_error
    )
    distinctHapsB <- out[["distinctHapsB"]]
    distinctHapsIE <- out[["distinctHapsIE"]]            
    hapMatcher <- out[["hapMatcher"]]
    eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
    eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]
    nrow_which_hapMatcher_0 <- out[["nrow_which_hapMatcher_0"]]
    
    genProbsM_t_new <- array(0, c(3, nSNPs))
    genProbsF_t_new <- array(0, c(3, nSNPs))
    hapProbs_t_new <- array(0, c(3, nSNPs))
    
    rcpp_calculate_gibbs_small_genProbs_and_hapProbs_using_binary_objects(
        genProbsM_t = genProbsM_t_new,
        genProbsF_t = genProbsF_t_new,
        hapProbs_t = hapProbs_t_new,
        gammaMT_t = gammaMT_t,
        gammaMU_t = gammaMU_t,
        gammaP_t = gammaP_t,
        hapMatcher = hapMatcher,
        distinctHapsIE = distinctHapsIE,
        which_haps_to_use = which_haps_to_use,
        ref_error = ref_error,
        rhb_t = rhb_t
    )

    expect_equal( genProbsM_t_new,  genProbsM_t)
    expect_equal( genProbsF_t_new,  genProbsF_t)
    expect_equal( hapProbs_t_new, hapProbs_t)
    
})

test_that("can avoid inflating fhb_t using eHapsCurrent_tc to make eMatRead_t", {

    ##
    ## make some dummy input
    ## want reads, ref panel, etc
    ##
    ## simulate some stuff
    K <- 100
    nMaxDH <- 80
    Ksmall <- 10 ## i.e. the subset
    nSNPs <- 200
    nReads <- 20
    Jmax <- 1000
    S <- 1
    ref_error <- 0.01
    maxDifferenceBetweenReads <- 1e100
    rescale_eMatRead_t <- TRUE
    
    ## make some test data
    which_haps_to_use <- sort(sample(K, Ksmall)) ## 1-based
    out <- make_fb_test_package(
        K = K,
        nReads = nReads,
        nSNPs = nSNPs,
        gridWindowSize = 32
    )
    rhi <- t(round(out$eHapsCurrent_tc[, , 1]))
    rhi_t <- t(rhi)
    rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)
    ##
    sampleReads <- out$sampleReads
    nSNPs <- out$nSNPs
    grid <- out$grid
    ## merge together some reads, make a bit wacky
    read1 <- sampleReads[[1]]
    read2 <- sampleReads[[round(nReads / 2)]]
    read3 <- sampleReads[[nReads]]
    sampleReads[[nReads + 1]] <- list(
        read1[[1]] + read2[[1]] + read3[[1]] + 3 - 1,
        3,
        matrix(c(read1[[3]], read2[[3]], read3[[3]]), ncol = 1),
        matrix(c(read1[[4]], read2[[4]], read3[[4]]), ncol = 1)
    )
    nReads <- length(sampleReads)
    
    ##
    ## here is original version, using inflation, resulting in eMatRead_t
    ##
    small_eHapsCurrent_tc <- array(0, c(Ksmall, nSNPs, S))    
    inflate_fhb_t_in_place(
        rhb_t = rhb_t,
        small_eHapsCurrent_tc,
        haps_to_get = which_haps_to_use - 1,
        nSNPs = nSNPs,
        ref_error = ref_error
    )
    eMatRead_t_old <- array(1, c(Ksmall, nReads))
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t_old,
        sampleReads = sampleReads,
        eHapsCurrent_tc = small_eHapsCurrent_tc,
        s = 0,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = Jmax,
        eMatHapOri_t = matrix(0, nrow = 0, ncol = 0),
        pRgivenH1 = vector(),
        pRgivenH2 = vector(),
        prev = 1,
        suppressOutput = 1,
        prev_section = "N/A",
        next_section = "N/A",
        run_pseudo_haploid = FALSE,
        rescale_eMatRead_t = rescale_eMatRead_t
    );
    
    ##
    ## here will be new version, based on subset, doing the same thing
    ##

    ## first, build useful things
    out <- make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nSNPs,
        ref_error = ref_error
    )
    distinctHapsB <- out[["distinctHapsB"]]
    distinctHapsIE <- out[["distinctHapsIE"]]            
    hapMatcher <- out[["hapMatcher"]]
    eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
    eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]
    nrow_which_hapMatcher_0 <- out[["nrow_which_hapMatcher_0"]]

    eMatRead_t <- R_make_eMatRead_t_for_gibbs_using_objects(
        sampleReads = sampleReads,
        hapMatcher = hapMatcher,
        grid = grid,
        rhb_t = rhb_t,
        distinctHapsIE = distinctHapsIE,
        ref_error = ref_error,
        which_haps_to_use = which_haps_to_use,
        normalize = TRUE,
        do_checks = TRUE,
        rhi_t = rhi_t
    )
    expect_equal(eMatRead_t, eMatRead_t_old)


    
    eMatRead_t_new <- array(1, c(Ksmall, nReads))
    Rcpp_make_eMatRead_t_for_gibbs_using_objects(
        eMatRead_t = eMatRead_t_new,
        sampleReads = sampleReads,
        hapMatcher = hapMatcher,
        grid = grid,
        rhb_t = rhb_t,
        distinctHapsIE = distinctHapsIE,
        ref_error = ref_error,
        which_haps_to_use = which_haps_to_use,
        rescale_eMatRead_t = rescale_eMatRead_t,
        Jmax = Jmax,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads
    )
    expect_equal(eMatRead_t_new, eMatRead_t_old)    


    
})






test_that("profile using riyan fish data", {

    skip("speed test")
    setwd("/data/smew1/rdavies/riyan_debug_2021_03_15")

    i_chr <- 2
    
    ## compare speed between old, new
    set.seed(i_chr)
    regionStart <- regionEnd <- buffer <- NA
    chr <- paste0("chr", i_chr)
    bams <- dir("bams")
    bams <- bams[-grep(".bai", bams)][1:2]
    write.table(
        matrix(paste0("bams/", bams), ncol = 1),
        file = "bamlist.txt",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = ""
    )
    
    ## 
    QUILT(
        outputdir = "/data/smew1/rdavies/riyan_debug_2021_03_15",
        chr = chr,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        bamlist = "bamlist.txt",
        override_default_params_for_small_ref_panel = TRUE,
        nCores = 1,
        nGibbsSamples = 2,
        block_gibbs_iterations = c(3, 6),
        n_gibbs_burn_in_its = 10
    )

    if (1 == 0) {

        load("~/temp/20191209_Plate1_10A.RData")
        
        temp_compare_two_versions(
                rhb_t,
                small_eHapsCurrent_tc,
                which_haps_to_use,
                nSNPs,
                ref_error,
                maxDifferenceBetweenReads,
                Jmax,
                hapMatcher,
                grid,
                distinctHapsIE,
                sampleReads
        )



    }

})


test_that("profile using HRC data", {

    skip("speed test")
    
    setwd("~/proj/QUILT/")
    dir <- "/data/smew1/rdavies/quilt_data/hrc_2021_04_20/2021_04_20_bams"
    f <- dir(dir)[grep(".1.0.bam", dir(dir))]
    f <- f[-grep("bai", f)]
    bamlist <- tempfile()
    write.table(matrix(paste0(dir, "/", f), ncol = 1), file = bamlist, row.names = FALSE, col.names = FALSE, quote = FALSE)
    QUILT(
        outputdir = "/data/smew1/rdavies/quilt_data/hrc_2021_04_20/",
        chr = "chr20",
        regionStart = 1,
        regionEnd = 2000000,
        buffer = 10000,
        bamlist= bamlist,
        nCores = 3
    )

    if (1 == 0) {

        temp_compare_two_versions <- function(
                                              rhb_t,
                                              small_eHapsCurrent_tc,
                                              which_haps_to_use,
                                              nSNPs,
                                              ref_error,
                                              maxDifferenceBetweenReads,
                                              Jmax,
                                              hapMatcher,
                                              grid,
                                              distinctHapsIE,
                                              sampleReads
                                              ) {

            K <- length(which_haps_to_use)
            rescale_eMatRead_t <- TRUE
            nReads <- length(sampleReads)
            
            f1 <- function() {
                ##
                S <- 1
                ## argh
                ## print(paste0("start = ", Sys.time()))
                inflate_fhb_t_in_place(
                    rhb_t,
                    small_eHapsCurrent_tc,
                    haps_to_get = which_haps_to_use - 1,
                    nSNPs = nSNPs,
                    ref_error = ref_error
                )
                eMatRead_t_old <- array(1, c(K, nReads))        
                rcpp_make_eMatRead_t(
                    eMatRead_t = eMatRead_t_old,
                    sampleReads = sampleReads,
                    eHapsCurrent_tc = small_eHapsCurrent_tc,
                    s = 0,
                    maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                    Jmax = Jmax,
                    eMatHapOri_t = matrix(0, nrow = 0, ncol = 0),
                    pRgivenH1 = vector(),
                    pRgivenH2 = vector(),
                    prev = 1,
                    suppressOutput = 1,
                    prev_section = "N/A",
                    next_section = "N/A",
                    run_pseudo_haploid = FALSE,
                    rescale_eMatRead_t = rescale_eMatRead_t
                )
                eMatRead_t_old
            }


            f2 <- function() {
                eMatRead_t_new <- array(1, c(K, nReads))
                Rcpp_make_eMatRead_t_for_gibbs_using_objects(
                    eMatRead_t = eMatRead_t_new,
                    sampleReads = sampleReads,
                    hapMatcher = hapMatcher,
                    grid = grid,
                    rhb_t = rhb_t,
                    distinctHapsIE = distinctHapsIE,
                    ref_error = ref_error,
                    which_haps_to_use = which_haps_to_use,
                    rescale_eMatRead_t = rescale_eMatRead_t,
                    Jmax = Jmax,
                    maxDifferenceBetweenReads = maxDifferenceBetweenReads
                )
                eMatRead_t_new
            }

            library("testthat")
            expect_equal(f1(), f2())

            library("microbenchmark")
            ## can be slower (on tall data, HRC sized)
            ## so can be faster (on wide data)
            print(microbenchmark(f1(), f2(), times = 2))
            ## so def slower. hmm

        }

        
        skip("woo")
        for(sample_name in c("NA12878", "NA12878HT", "NA12878ONT")) {
            print(sample_name)
            load(paste0("/dev/shm/rwdavies/", sample_name, ".RData")    )
            temp_compare_two_versions(
                rhb_t,
                small_eHapsCurrent_tc,
                which_haps_to_use,
                nSNPs,
                ref_error,
                maxDifferenceBetweenReads,
                Jmax,
                hapMatcher,
                grid,
                distinctHapsIE,
                sampleReads
            )
        }

        ## slower but no problems!

    }

})


test_that("profile genProbs making with or without eHapsCurrent_tc", {

    skip("test")
    
    for(sample_name in c("NA12878", "NA12878HT", "NA12878ONT")) {
        
            print(sample_name)
            load(paste0("/dev/shm/rwdavies/", sample_name, ".RData")    )
snp_start_1_based <- 1
    snp_end_1_based <- nSNPs    
class(hapMatcher) <- "integer"
    
            f1 <- function() {
        inflate_fhb_t_in_place(
            rhb_t,
            small_eHapsCurrent_tc,
            haps_to_get = which_haps_to_use - 1,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
        rcpp_calculate_gn_genProbs_and_hapProbs(
            genProbsM_t = genProbsM_t,
            genProbsF_t = genProbsF_t,
            hapProbs_t = hapProbs_t,
            s = 0,
            eHapsCurrent_tc = small_eHapsCurrent_tc,
            gammaMT_t = gammaMT_t,
            gammaMU_t = gammaMU_t,
            gammaP_t = gammaP_t,
            grid = grid,
            snp_start_1_based = snp_start_1_based,
            snp_end_1_based = snp_end_1_based,
            grid_offset = 0
        )
        list(genProbsM_t = genProbsM_t,
            genProbsF_t = genProbsF_t,
            hapProbs_t = hapProbs_t)
    }


        
        f2 <- function() {
            ## rhb_t_subset <- rhb_t[which_haps_to_use, ]
    rcpp_calculate_gibbs_small_genProbs_and_hapProbs_using_binary_objects(
        genProbsM_t = genProbsM_t_new,
        genProbsF_t = genProbsF_t_new,
        hapProbs_t = hapProbs_t_new,
        gammaMT_t = gammaMT_t,
        gammaMU_t = gammaMU_t,
        gammaP_t = gammaP_t,
        hapMatcher = hapMatcher,
        distinctHapsIE = distinctHapsIE,
        which_haps_to_use = which_haps_to_use,
        ref_error = ref_error,
        rhb_t = rhb_t
    )
        list(genProbsM_t = genProbsM_t_new,
            genProbsF_t = genProbsF_t_new,
            hapProbs_t = hapProbs_t_new)
    }

    genProbsM_t <- array(0, c(3, nSNPs))
    genProbsF_t <- array(0, c(3, nSNPs))
    hapProbs_t <- array(0, c(3, nSNPs))
    
    
    genProbsM_t_new <- array(0, c(3, nSNPs))
    genProbsF_t_new <- array(0, c(3, nSNPs))
    hapProbs_t_new <- array(0, c(3, nSNPs))

        inflate_fhb_t_in_place(
            rhb_t,
            small_eHapsCurrent_tc,
            haps_to_get = which_haps_to_use - 1,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
    
    f <- function() {
        K <- nrow(small_eHapsCurrent_tc)
        nGrids <- ncol(rhb_t)
        gammaMT_t <- array(runif(K * nGrids), c(K, nGrids))
        for(iGrid in 1:nGrids) {
            gammaMT_t[, iGrid] <- gammaMT_t[, iGrid] / sum(gammaMT_t[, iGrid])
        }
        gammaMT_t
    }
    gammaMT_t <- f()
    gammaMU_t <- f()
    gammaP_t <- f()
    library("testthat")
    expect_equal(f1(), f2(), tol = 1e-5)


    library("microbenchmark")
    ## can be slower (on tall data, HRC sized)
    ## so can be faster (on wide data)
    print(microbenchmark(f1(), f2(), times = 2))
    ## so DEF slower, much. hmm


        }

        ## other one
        
        ## slower but no problems!
        ## hmm, though does add up, when run repeatedly
        ## can I make faster, or is it just what it is?


})
