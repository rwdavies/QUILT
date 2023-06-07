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



test_that("can make old eHapsCurrent_tc using rare common idea", {

    set.seed(1011010)
    
    ref_error <- 0.001
    maxDifferenceBetweenReads <- 1e6

    nReads <- 800
    nSNPs <- 2000
    nCommonSNPs <- 100
    K <- 1000
    Ksubset <- 50
    snp_is_common <- rep(FALSE, nSNPs)
    snp_is_common[sample(nSNPs, nCommonSNPs)] <- TRUE

    rhi <- array(0L, c(nSNPs, K))
    rhi[snp_is_common, ] <- sample(c(0L, 1L), K * nCommonSNPs, prob = c(0.7, 0.3), replace = TRUE)
    rhi[!snp_is_common, ] <- sample(c(0L, 1L), K * (nSNPs - nCommonSNPs), prob = c(0.99, 0.01), replace = TRUE)
    rhi_t <- t(rhi)

    ## now, make objects
    nGrids <- ceiling(nSNPs / 32)
    grid <- floor((1:nSNPs) / 32)
    nCommonGrids <- ceiling(nCommonSNPs / 32)
    rhb <- array(0L, c(nCommonGrids, K))
    for(k in 1:K) {
        rhb[, k] <- rcpp_int_contract(rhi[snp_is_common, k])
    }
    rhb_t <- t(rhb)


    ## make special objects
    out <- make_rhb_t_equality(
        rhb_t = rhb_t,
        nSNPs = nCommonSNPs,
        ref_error = ref_error,
        use_hapMatcherR = TRUE
    )
    hapMatcherR <- out[["hapMatcherR"]]
    eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]
    eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
    distinctHapsIE <- out[["distinctHapsIE"]]

    ## manually make special rare one too,
    snp_is_rare_1_based <- which(!snp_is_common)    
    rare_per_hap_info <- sapply(1:K, function(k) {
        ## for now store position among all SNPs
        snp_is_rare_1_based[which(rhi[!snp_is_common, k] == 1)]
    })


    which_haps_to_use <- sort(sample(1:K, Ksubset))
    
    old_eHapsCurrent_tc <- array(0, c(Ksubset, nSNPs, 1))
    old_eHapsCurrent_tc[, , 1][rhi_t[which_haps_to_use, ] == 0] <- ref_error
    old_eHapsCurrent_tc[, , 1][rhi_t[which_haps_to_use, ] == 1] <- (1 - ref_error)

    snp_is_common_1_based <- which(snp_is_common)

    new_eHapsCurrent_tc <- make_eHapsCurrent_tc_using_rare_and_common_stuff(
        hapMatcherR = hapMatcherR,
        distinctHapsIE = distinctHapsIE,
        eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
        eMatDH_special_matrix = eMatDH_special_matrix,
        rare_per_hap_info = rare_per_hap_info,
        snp_is_common = snp_is_common,
        which_haps_to_use = which_haps_to_use,
        snp_is_common_1_based = snp_is_common_1_based,
        Ksubset = Ksubset,
        ref_error = ref_error
    )

    expect_equal(new_eHapsCurrent_tc, old_eHapsCurrent_tc)

    ## ## get some fake-ish sample reads
    ## ## honestly they can be very fake here
    ## sampleReads <- lapply(1:nReads, function(ii) {
    ##     rl <- sample(1:5, 1) ## length
    ##     s <- sample(1:(nSNPs - rl), 1) ## start
    ##     e <- s + rl - 1 ## end
    ##     cg <- grid[round(mean(c(s, e)))] ## central grid
    ##     bq <- matrix(sample(c(-25:-20, 26:30), rl, replace = TRUE), ncol = 1)
    ##     u <- matrix(s:e - 1, ncol = 1)
    ##     list(
    ##         rl - 1,
    ##         cg,
    ##         bq,
    ##         u
    ##     )
    ## })
    ## sampleReads <- sampleReads[order(sapply(sampleReads, function(x) x[[2]]))]

    
    ## original_eMatRead_t <- array(1, c(Ksubset, nReads))
    ## rcpp_make_eMatRead_t(
    ##     eMatRead_t = original_eMatRead_t,
    ##     sampleReads = sampleReads,
    ##     eHapsCurrent_tc = eHapsCurrent_tc,
    ##     maxDifferenceBetweenReads = maxDifferenceBetweenReads,
    ##     s = 0,
    ##     Jmax = 100,
    ##     eMatHapOri_t = array(0, c(1, 1)),
    ##     pRgivenH1 = array(0, 1),
    ##     pRgivenH2 = array(0, 1),
    ##     prev = 0,
    ##     suppressOutput = 1,
    ##     prev_section = "wer",
    ##     next_section = "wer",
    ##     rescale_eMatRead_t = TRUE,
    ##     run_pseudo_haploid = FALSE
    ## )
    

    ## ## OK am here
    ## ## want to reproduce this using



    ## eMatRead_t <- array(1, c(Ksubset, nReads))    
    
    
})
