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



test_that("can avoid using eHapsCurrent_tc and gammas in genProbs calculation", {

    set.seed(12)

    ## OK am here
    ## want to test this function
    ## using either eHapsCurrent_tc, or other approach used for eMatRead_t
    ## should be more straightforward? depending on offset?

    K <- 100
    nMaxDH <- 80
    Ksmall <- 10 ## i.e. the subset
    nSNPs <- 500
    nCommonSNPs <- 100
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
    ## actuall re-do a bit to make some rarer and some more common
    snp_is_common <- rep(FALSE, nSNPs)
    snp_is_common[sample(nSNPs, nCommonSNPs)] <- TRUE   
    rhi_t[, !snp_is_common] <- matrix(sample(c(0, 1), (nSNPs - nCommonSNPs) * K, replace = TRUE, prob = c(0.9, 0.1)), nrow = K, ncol = nSNPs - nCommonSNPs)
    rhi <- t(rhi_t)
    rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)

    ##
    f <- function() {
        alphaHat_t <- array(runif(Ksmall * nGrids), c(Ksmall, nGrids))
        betaHat_t <- array(runif(Ksmall * nGrids), c(Ksmall, nGrids))
        gammaMT_t <- alphaHat_t * betaHat_t
        c <- colSums(gammaMT_t)
        gammaMT_t <- gammaMT_t / rep(c, each = nrow(gammaMT_t))
        stopifnot(abs(colSums(gammaMT_t)  - 1) < 1e-10)
        return(list(alphaHat_t, betaHat_t, c, gammaMT_t))
    }
    
    out1 <- f()
    alphaHat_t1 <- out1[[1]]
    betaHat_t1 <- out1[[2]]
    c1 <- out1[[3]]
    gamma_t1 <- out1[[4]]

    out2 <- f()
    alphaHat_t2 <- out2[[1]]
    betaHat_t2 <- out2[[2]]
    c2 <- out2[[3]]
    gamma_t2 <- out2[[4]]

    out3 <- f()
    alphaHat_t3 <- out3[[1]]
    betaHat_t3 <- out3[[2]]
    c3 <- out3[[3]]
    gamma_t3 <- out3[[4]]
    
    ## also
    gammaMT_t <- gamma_t1
    gammaMU_t <- gamma_t2
    gammaP_t <- gamma_t3
    
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
        grid_offset = 0,
        sample_is_diploid = FALSE
    )




    ##
    ## new version
    ##
    out <- STITCH::make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nSNPs,
        ref_error = ref_error,
        use_hapMatcherR = FALSE
    )
    distinctHapsB <- out[["distinctHapsB"]]
    distinctHapsIE <- out[["distinctHapsIE"]]
    hapMatcher <- out[["hapMatcher"]]
    eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
    eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]
    nrow_which_hapMatcher_0 <- out[["nrow_which_hapMatcher_0"]]
    eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]
    eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
    out <- STITCH::make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nSNPs,
        ref_error = ref_error,
        use_hapMatcherR = TRUE
    )
    hapMatcherR <- out[["hapMatcherR"]]

    for(use_eMatDH_special_symbols in c(TRUE, FALSE)) {

        for(use_hapMatcherR in c(FALSE, TRUE)) {

            for(calculate_gamma_on_the_fly in c(FALSE, TRUE)) {

                if (calculate_gamma_on_the_fly) {
                    alphaHat_t1 <- out1[[1]];     betaHat_t1 <- out1[[2]];     c1 <- out1[[3]];     gamma_t1 <- out1[[4]]
                    alphaHat_t2 <- out2[[1]];    betaHat_t2 <- out2[[2]];     c2 <- out2[[3]];     gamma_t2 <- out2[[4]]
                    alphaHat_t3 <- out3[[1]];     betaHat_t3 <- out3[[2]];     c3 <- out3[[3]];     gamma_t3 <- out3[[4]]
                    gammaMT <- array(0, c(1, 1)); gammaMU_t <- array(0, c(1, 1)); gammaP_t <- array(0, c(1, 1));                    
                } else {
                    alphaHat_t1 <- array(0, c(1, 1)); betaHat_t1 <- array(0, c(1, 1)); c1 <- array(0, 1)
                    alphaHat_t2 <- array(0, c(1, 1)); betaHat_t2 <- array(0, c(1, 1)); c2 <- array(0, 1)
                    alphaHat_t3 <- array(0, c(1, 1)); betaHat_t3 <- array(0, c(1, 1)); c3 <- array(0, 1)
                    gammaMT <- out1[[4]]; gammaMU_t <- out2[[4]]; gammaP_t <- out3[[4]]
                }

                genProbsM_t_new <- array(0, c(3, nSNPs))
                genProbsF_t_new <- array(0, c(3, nSNPs))
                hapProbs_t_new <- array(0, c(3, nSNPs))

                rcpp_calculate_gibbs_small_genProbs_and_hapProbs_using_binary_objects(
                    alphaHat_t1 = alphaHat_t1,
                    alphaHat_t2 = alphaHat_t2,
                    alphaHat_t3 = alphaHat_t3,     
                    betaHat_t1 = betaHat_t1,
                    betaHat_t2 = betaHat_t2, 
                    betaHat_t3 = betaHat_t3,
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    genProbsM_t = genProbsM_t_new,
                    genProbsF_t = genProbsF_t_new,
                    hapProbs_t = hapProbs_t_new,
                    gammaMT_t = gammaMT_t,
                    gammaMU_t = gammaMU_t,
                    gammaP_t = gammaP_t,
                    hapMatcher = hapMatcher,
                    hapMatcherR = hapMatcherR,
                    use_hapMatcherR = use_hapMatcherR,
                    distinctHapsB = distinctHapsB,
                    distinctHapsIE = distinctHapsIE,
                    which_haps_to_use = which_haps_to_use,
                    ref_error = ref_error,
                    rhb_t = rhb_t,
                    eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                    eMatDH_special_matrix = eMatDH_special_matrix,
                    use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                    calculate_gamma_on_the_fly = calculate_gamma_on_the_fly,
                    sample_is_diploid = FALSE
                )
                
                expect_equal( hapProbs_t_new, hapProbs_t)
                expect_equal( genProbsM_t_new, genProbsM_t)
                expect_equal( genProbsF_t_new, genProbsF_t)

            }

        }

    }


    
    ##
    ## now do the same thing but using rare common ideas
    ##
    nGrids <- ceiling(nSNPs / 32)
    grid <- floor((1:nSNPs) / 32)
    nCommonGrids <- ceiling(nCommonSNPs / 32)
    rhb <- array(0L, c(nCommonGrids, K))
    
    for(k in 1:K) {
        rhb[, k] <- rcpp_int_contract(rhi[snp_is_common, k])
    }
    rhb_t <- t(rhb)
    snp_is_rare_1_based <- which(!snp_is_common)        
    rare_per_hap_info <- sapply(1:K, function(k) {
        ## this is among rare SNPs
        snp_is_rare_1_based[which(rhi[!snp_is_common, k] == 1)]
    })
    common_snp_index <- integer(nSNPs)
    common_snp_index[which(snp_is_common)] <- 1:nCommonSNPs ## 1-based


    out <- STITCH::make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nCommonSNPs,
        ref_error = ref_error,
        use_hapMatcherR = TRUE
    )
    distinctHapsB <- out[["distinctHapsB"]]
    distinctHapsIE <- out[["distinctHapsIE"]]
    eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
    eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]
    nrow_which_hapMatcher_0 <- out[["nrow_which_hapMatcher_0"]]
    eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]
    eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
    hapMatcherR <- out[["hapMatcherR"]]
    

    genProbsM_t_new <- array(0, c(3, nSNPs))
    genProbsF_t_new <- array(0, c(3, nSNPs))
    hapProbs_t_new <- array(0, c(3, nSNPs))
    
    alphaHat_t1 <- out1[[1]];     betaHat_t1 <- out1[[2]];     c1 <- out1[[3]];     gamma_t1 <- out1[[4]]
    alphaHat_t2 <- out2[[1]];    betaHat_t2 <- out2[[2]];     c2 <- out2[[3]];     gamma_t2 <- out2[[4]]
    alphaHat_t3 <- out3[[1]];     betaHat_t3 <- out3[[2]];     c3 <- out3[[3]];     gamma_t3 <- out3[[4]]
    gammaMT <- array(0, c(1, 1)); gammaMU_t <- array(0, c(1, 1)); gammaP_t <- array(0, c(1, 1));                    


    ## this is the overall sites
    common_snp_overall0 <- which(snp_is_common) - 1
    common_snp_grids0 <- floor(common_snp_overall0 / 32)

    ## hopefully not a RAM monster
    ## create a new one that is per-site, and gives list of k
    rare_per_snp_info <- lapply(1:nSNPs, function(x) -1L)
    ## want this to be the index of the rare SNP in the overall
    snp_is_rare_1_based <- which(!snp_is_common)            
    for(k in 1:length(which_haps_to_use)) {
        ## expand back out to all SNPs
        snps <- rare_per_hap_info[[which_haps_to_use[[k]]]] ## this IS an index among all SNPs
        ## this is among all SNPs
        for(snp in snps) {
            rare_per_snp_info[[snp]] <- c(rare_per_snp_info[[snp]], k) ## keep 1-based
        }
    }
    
    rcpp_calculate_genProbs_and_hapProbs_final_rare_common(
        alphaHat_t1 = alphaHat_t1,
        alphaHat_t2 = alphaHat_t2,
        alphaHat_t3 = alphaHat_t3,     
        betaHat_t1 = betaHat_t1,
        betaHat_t2 = betaHat_t2, 
        betaHat_t3 = betaHat_t3,
        c1 = c1,
        c2 = c2,
        c3 = c3,
        genProbsM_t = genProbsM_t_new,
        genProbsF_t = genProbsF_t_new,
        hapProbs_t = hapProbs_t_new,
        gammaMT_t = gammaMT_t,
        gammaMU_t = gammaMU_t,
        gammaP_t = gammaP_t,
        hapMatcherR = hapMatcherR,
        distinctHapsB = distinctHapsB,
        distinctHapsIE = distinctHapsIE,
        which_haps_to_use = which_haps_to_use,
        ref_error = ref_error,
        rhb_t = rhb_t,
        eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
        eMatDH_special_matrix = eMatDH_special_matrix,
        use_eMatDH_special_symbols = use_eMatDH_special_symbols,
        rare_per_hap_info = rare_per_hap_info,
        common_snp_index = common_snp_index,
        snp_is_common = snp_is_common,
        rare_per_snp_info = rare_per_snp_info,
        sample_is_diploid = FALSE
    )

    expect_equal( hapProbs_t_new, hapProbs_t)
    expect_equal( genProbsM_t_new, genProbsM_t)
    expect_equal( genProbsF_t_new, genProbsF_t)
    
    ## bad_cols <- which(colSums(abs(hapProbs_t_new - hapProbs_t)) > 0.01)
    ## print("bad cols")
    ## print(head(bad_cols))
    ## print("first 10 cols")
    ## print(snp_is_common[1:10])

    ## print("lets check out those cols")
    ## print(hapProbs_t_new[1:2, 1:10])
    ## print(hapProbs_t[1:2, 1:10])
    ## print(hapProbs_t_new[1:2, 1:10] / hapProbs_t[1:2, 1:10])

    
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
    out <- STITCH::make_rhb_t_equality(
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
    eMatDH_special_symbols_list <- out[["eMatDH_special_symbols_list"]]
    nrow_which_hapMatcher_0 <- out[["nrow_which_hapMatcher_0"]]
    eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
    eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]

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


    for(use_eMatDH_special_symbols in c(FALSE, TRUE)) {

        for(use_hapMatcherR in c(FALSE, TRUE)) {

            out <- STITCH::make_rhb_t_equality(
                rhb_t = rhb_t,
                nMaxDH = nMaxDH,
                nSNPs = nSNPs,
                ref_error = ref_error,
                use_hapMatcherR = use_hapMatcherR
            )
            hapMatcher <- out[["hapMatcher"]]
            hapMatcherR <- out[["hapMatcherR"]]

            eMatRead_t_new <- array(1, c(Ksmall, nReads))
            Rcpp_make_eMatRead_t_for_gibbs_using_objects(
                eMatRead_t = eMatRead_t_new,
                sampleReads = sampleReads,
                hapMatcher = hapMatcher,
                hapMatcherR = hapMatcherR,
                use_hapMatcherR = use_hapMatcherR,
                grid = grid,
                rhb_t = rhb_t,
                distinctHapsIE = distinctHapsIE,
                eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                eMatDH_special_matrix = eMatDH_special_matrix,
                ref_error = ref_error,
                which_haps_to_use = which_haps_to_use,
                rescale_eMatRead_t = rescale_eMatRead_t,
                Jmax = Jmax,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                use_eMatDH_special_symbols = use_eMatDH_special_symbols
            )
            expect_equal(eMatRead_t_new, eMatRead_t_old)

        }

    }

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
        rhb_t = rhb_t,
        sample_is_diploid = FALSE
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
