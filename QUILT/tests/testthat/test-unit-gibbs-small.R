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
