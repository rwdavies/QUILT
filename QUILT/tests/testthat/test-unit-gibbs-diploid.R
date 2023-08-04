if (1 == 0) {

    library("testthat"); library("STITCH"); library("rrbgen")
    dir <- "~/proj/QUILT/"
    setwd(paste0(dir, "/QUILT/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))


}



test_that("can identify all 1 eMatRead_t cols, and get 0-based indices of non-1 entries otherwise", {

    ## make sure a bunch of cases are here
    K <- 50
    nReads <- 200
    eMatRead_t <- sapply(1:nReads, function(iRead) {
        x <- rpois(1, 1) + 1
        y <- c(1, runif(x - 1))
        prob <- c(1, rep(0.1, x - 1))
        y[sample(1:x, K, prob = prob / sum(prob), replace =  TRUE)]
    })

    outR <- evaluate_read_variability(eMatRead_t)

    ## make sure everything OK here
    expect_equal(sort(unique(outR[["read_category"]])), 0:3)
    
    number_of_non_1_reads <- integer(nReads)
    indices_of_non_1_reads <- array(0L, c(K, nReads))
    read_category <- integer(nReads)
    rcpp_evaluate_read_variability(eMatRead_t, number_of_non_1_reads, indices_of_non_1_reads, read_category)

    expect_equal(outR[["number_of_non_1_reads"]], number_of_non_1_reads)
    expect_equal(outR[["indices_of_non_1_reads"]], indices_of_non_1_reads)
    expect_equal(outR[["read_category"]], read_category)

    ## print(outR[["indices_of_non_1_reads"]][1:5, 1:5])
    ## print(indices_of_non_1_reads[1:5, 1:5])

})



test_that("subtraction approach works", {

    set.seed(81)
    K <- 100
    alphaHat_m <- matrix(runif(K * 3), 3, K)
    betaHat_m <- matrix(runif(K * 3), 3, K)
    eMatRead_t <- matrix(c(1, runif(3))[sample(1:2, K, replace = TRUE)], K, 1)
    outR <- evaluate_read_variability(eMatRead_t)
    number_of_non_1_reads <- outR[[1]]
    indices_of_non_1_reads <- outR[[2]]
    read_category <- outR[[3]]
    
    ## 
    pC <- rowSums(alphaHat_m * betaHat_m)
    ab_m <- t(alphaHat_m * betaHat_m)

    for(iii in 1:6) {
        
        if (iii == 1) { h <- c(1, 2, 3) }
        if (iii == 2) { h <- c(1, 3, 2) }
        if (iii == 3) { h <- c(2, 1, 3) }
        if (iii == 4) { h <- c(2, 1, 3) }
        if (iii == 5) { h <- c(3, 1, 2) }
        if (iii == 6) { h <- c(3, 2, 1) }
        h_rC <- h[1]
        h_rA1 <- h[2]
        h_rA2 <- h[3]

        iRead <- 1
        
        ## check it out down here
        ## recall 0 = normal
        read_category <- 0L
        out0 <- evaluate_read_probabilities(alphaHat_m, betaHat_m, pC, read_category, iRead, h_rC, h_rA1, h_rA2, eMatRead_t, number_of_non_1_reads, indices_of_non_1_reads)

        ## rcpp 
        rcpp0_pA1 <- numeric(3)
        rcpp0_pA2 <- numeric(3)
        rcpp_evaluate_read_probabilities(alphaHat_m, betaHat_m, ab_m, pC, rcpp0_pA1, rcpp0_pA2, read_category, iRead - 1, h_rC - 1, h_rA1 - 1, h_rA2 - 1, eMatRead_t, number_of_non_1_reads, indices_of_non_1_reads, FALSE)
        
        ## 2 = subraction, same vals
        read_category <- 2L
        out2 <- evaluate_read_probabilities(alphaHat_m, betaHat_m, pC, read_category, iRead, h_rC, h_rA1, h_rA2, eMatRead_t, number_of_non_1_reads, indices_of_non_1_reads)

        ## rcpp
        rcpp2_pA1 <- numeric(3)
        rcpp2_pA2 <- numeric(3)
        rcpp_evaluate_read_probabilities(alphaHat_m, betaHat_m, ab_m, pC, rcpp2_pA1, rcpp2_pA2, read_category, iRead - 1, h_rC - 1, h_rA1 - 1, h_rA2 - 1, eMatRead_t, number_of_non_1_reads, indices_of_non_1_reads, FALSE)
        
        ## 3 = subraction
        read_category <- 3L
        out3 <- evaluate_read_probabilities(alphaHat_m, betaHat_m, pC, read_category, iRead, h_rC, h_rA1, h_rA2, eMatRead_t, number_of_non_1_reads, indices_of_non_1_reads)

        ## rcpp
        rcpp3_pA1 <- numeric(3)
        rcpp3_pA2 <- numeric(3)
        rcpp_evaluate_read_probabilities(alphaHat_m, betaHat_m, ab_m, pC, rcpp3_pA1, rcpp3_pA2, read_category, iRead - 1, h_rC - 1, h_rA1 - 1, h_rA2 - 1, eMatRead_t, number_of_non_1_reads, indices_of_non_1_reads, FALSE)

        expect_equal(out0, out2)
        expect_equal(out0, out3)
        
        expect_equivalent(out0[["pA1"]], rcpp0_pA1)
        expect_equivalent(out0[["pA2"]], rcpp0_pA2)
        
        expect_equivalent(out0[["pA1"]], rcpp2_pA1)
        expect_equivalent(out0[["pA2"]], rcpp2_pA2)

        expect_equivalent(out0[["pA1"]], rcpp3_pA1)
        expect_equivalent(out0[["pA2"]], rcpp3_pA2)

        
    }

})


test_that("can skip reads and more efficiently calculate probabilities for gibbs sampling, specifically diploid", {

    
    ## have test package
    ## have R version
    ## have Rcpp version
    ## try the two approaches
    set.seed(919)
    
    sample_is_diploid <- TRUE
    verbose <- FALSE
    S <- 1
    test_package <- make_quilt_fb_test_package(
        K = 100,
        nReads = 100,
        nSNPs = 32 * 30,
        gridWindowSize = 32,
        S = S,
        simple_ematread = TRUE
    )
    eMatRead_t <- test_package$list_of_eMatRead_t[[1]]
    outR <- evaluate_read_variability(eMatRead_t)
    
    ## require sufficiently diverse
    expect_equal(sort(unique(outR[["read_category"]])), 0:3)    
    
    ff <- 0
    sampleReads <- test_package[["sampleReads"]]
    nReads <- length(sampleReads)
    transMatRate_tc_H <- test_package[["transMatRate_tc_H"]]
    alphaMatCurrent_tc <- test_package[["alphaMatCurrent_tc"]]
    ## argh
    
    priorCurrent_m <- test_package[["priorCurrent_m"]]
    eHapsCurrent_tc <- test_package[["eHapsCurrent_tc"]]
    grid <- test_package$grid
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    
    ## it is probably better if fit is much poorer

    ## check that with no emissions, move one forward keeps sum the same
    true_H <- test_package$true_H[1]
    nReads <- length(sampleReads)
    wif <- sapply(sampleReads, function(x) x[[2]]) + 1
    
    double_list_of_starting_read_labels <- lapply(1:S, function(s) {
        lapply(1:1, function(i_sampling) {
            sample(c(1, 2), nReads, replace = TRUE)
        })
    })

    n_gibbs_starts <- 1
    n_gibbs_burn_in_its <- 0
    n_gibbs_sample_its <- 1


    ## here check in R that if I make everything the slow approach, it works the same
    ## note, also do this in C++ later
    super_out_R <- lapply(c(FALSE, TRUE), function(force_reset_read_category_zero) {
        set.seed(123)
        outR <- forwardBackwardGibbsNIPT(
            sampleReads = sampleReads,
            priorCurrent_m = priorCurrent_m,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            eHapsCurrent_tc = eHapsCurrent_tc,
            transMatRate_tc_H = transMatRate_tc_H,
            ff = ff,
            grid = grid,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            n_gibbs_sample_its = n_gibbs_sample_its,
            n_gibbs_starts = n_gibbs_starts,
            use_starting_read_labels = TRUE,
            verbose = verbose,
            double_list_of_starting_read_labels = double_list_of_starting_read_labels,
            force_reset_read_category_zero = force_reset_read_category_zero,
            return_alpha = TRUE,
            return_extra = TRUE
        )
        outR
    })

    expect_equal(super_out_R[[1]][["H"]], super_out_R[[2]][["H"]])
    outR <- super_out_R[[1]]

    param_list <- list(
        return_alpha = TRUE,
        return_extra = TRUE,
        return_genProbs = FALSE,
        return_gamma = FALSE,
        return_hapProbs = FALSE,
        return_p_store = FALSE,
        return_p1 = FALSE,
        return_gibbs_block_output = FALSE,
        return_advanced_gibbs_block_output = FALSE,
        use_starting_read_labels = TRUE,
        verbose = FALSE,
        run_fb_subset = FALSE,
        haploid_gibbs_equal_weighting = TRUE,
        gibbs_initialize_iteratively = FALSE,
        gibbs_initialize_at_first_read = TRUE,
        use_smooth_cm_in_block_gibbs = TRUE,
        use_small_eHapsCurrent_tc = TRUE,
        sample_is_diploid = TRUE,
        update_in_place = FALSE,
        do_shard_ff0_block_gibbs = TRUE,
        force_reset_read_category_zero = FALSE
    )


    wif0 <- as.integer(sapply(sampleReads, function(x) x[[2]]))
    grid_has_read <- rep(FALSE, nGrids)
    grid_has_read[wif0 + 1] <- TRUE

    double_list_of_starting_read_labelsX <-double_list_of_starting_read_labels
    a <-  double_list_of_starting_read_labelsX[[1]][[1]][1] 
    double_list_of_starting_read_labelsX[[1]][[1]][1] <- 10
    double_list_of_starting_read_labelsX[[1]][[1]][1] <- a

    set.seed(123)
    outRCPP <- rcpp_forwardBackwardGibbsNIPT(
        sampleReads = sampleReads,
        priorCurrent_m = priorCurrent_m,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        eHapsCurrent_tc = eHapsCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        ff = ff,
        grid = grid,
        n_gibbs_burn_in_its = n_gibbs_burn_in_its,
        n_gibbs_sample_its = n_gibbs_sample_its,
        n_gibbs_starts = n_gibbs_starts,
        blocks_for_output = array(0, c(1, 1)),
        double_list_of_starting_read_labels = double_list_of_starting_read_labels,
        param_list = param_list,
        suppressOutput = 1,
        wif0 = integer(1),
        L_grid = integer(1),
        alphaHat_t1 = array(0, c(K, nGrids)),
        alphaHat_t2 = array(0, c(K, nGrids)),
        alphaHat_t3 = array(0, c(K, nGrids)),
        betaHat_t1 = array(0, c(K, nGrids)),
        betaHat_t2 = array(0, c(K, nGrids)),
        betaHat_t3 = array(0, c(K, nGrids)),
        eMatGrid_t1 = array(0, c(K, nGrids)),
        eMatGrid_t2 = array(0, c(K, nGrids)),
        eMatGrid_t3 = array(0, c(K, nGrids)),
        gammaMT_t_local = array(0, c(K, nGrids)),
        gammaMU_t_local = array(0, c(K, nGrids)),
        gammaP_t_local = array(0, c(K, nGrids)),
        hapSum_tc = array(0, c(1, 1, 1)),
        eMatDH_special_matrix_helper = array(0, c(1, 1)),
        eMatDH_special_matrix = array(0,c(1, 1)),
        use_eMatDH_special_symbols = FALSE,
        distinctHapsB = array(0L, c(1, 1)),
        distinctHapsIE = array(0, c(1, 1)),
        hapMatcher = array(0L, c(1, 1)),
        hapMatcherR = array(as.raw(0), c(1, 1)),
        use_hapMatcherR = FALSE,
        rhb_t = array(0L, c(1, 1)),
        ref_error = 0,
        which_haps_to_use = 0:(K - 1),
        grid_has_read = grid_has_read,
        smooth_cm = numeric(1),
        skip_read_iteration = FALSE
    )

    ## getting there...
    outR2 <- evaluate_read_variability(eMatRead_t)
    outR3 <- evaluate_read_variability(outRCPP$eMatRead_t)
    
    ## print(head(cbind(
    ##     outR2$read_category,
    ##     double_list_of_starting_read_labels[[1]][[1]],
    ##     double_list_of_starting_read_labelsX[[1]][[1]],        
    ##     outR[["H"]],
    ##     outRCPP[["H"]]
    ## ), 10))

    init <- double_list_of_starting_read_labelsX[[1]][[1]]
    w1 <- outR2$read_category == 1    
    expect_equal(outR[["H"]][w1], outRCPP[["H"]][w1])
    expect_equal(outR[["H"]][w1], init[w1])
    expect_equal(outRCPP[["H"]][w1], init[w1])        

    w2 <- outR2$read_category == 2
    expect_equal(outR[["H"]][w2], outRCPP[["H"]][w2])

    ## print(outR[["H"]][w2])
    ## print(outRCPP[["H"]][w2])
    ## print(outR2[["read_category"]])

    
    ## for(i in 1:3) {
    ##     if (i == 1) {
    ##         print("start")
    ##         h <- double_list_of_starting_read_labelsX[[1]][[1]][w2]
    ##     } else if (i == 2) {
    ##         print("R")
    ##         h <- outR[["H"]][w2]
    ##     } else if (i == 3) {
    ##         print("Rcpp")
    ##         h <- outRCPP[["H"]][w2]
    ##     }
    ##     print(which(w2[which(h == 1)]))
    ## }
    ## ## first difference at read = 15
          
    ## ##expect_equal(outR[["H"]][w2], init[w2])
    ## ##expect_equal(outRCPP[["H"]][w2], init[w2])        

    ## print("rcpp")
    ## print(outRCPP$alphaHat_t1[1:10, 1])
    ## print(outRCPP$c1[1])

    ## print("R")
    ## print(outR[["alphaHat_t1"]][1:10, 1])
    ## print(outR[["c1"]][1])
    
    ## ## This is the grid of the read
    ## ## 
    ## wif0
    
    ## ## AM HERE
    ## ## STILL A PROBLEM FUCK
    ## ## SOMETHING ABOUT STARTING GRID2, THE INITIAL PROBABILITIES pC ARE DIFFERENT
    ## ## UNDERSTAND WHY

    ## ## just try to move forward
    ## iGrid <- 2
    ## alphaHat_t1 <- outR[["alphaHat_t1"]]
    ## eMatGrid_t1 <- outR[["eMatGrid_t1"]]    
    ## s <- 1
    ## c1 <- outR[["c1"]]
    ## alphaHat_t1 <- alpha_forward_one(iGrid, K, s, alphaHat_t1, transMatRate_tc_H, eMatGrid_t1, alphaMatCurrent_tc)
    ## alphaHat_t1[, iGrid] <- alphaHat_t1[, iGrid] * c1[iGrid]
    ## a <- 1 / sum(alphaHat_t1[, iGrid])
    ## c1[iGrid] <- c1[iGrid] * a
    ## alphaHat_t1[, iGrid] <- alphaHat_t1[, iGrid] * a
    ## alphaHat_t1[1:5, 2]
    ## sum(alphaHat_t1[, 2])
    
    
    ## ## can I check the two versions
    ## X_alphaHat_t1 <- outR[["alphaHat_t1"]]
    ## rcpp_alpha_forward_one_QUILT_faster(
    ##     s = 0,
    ##     iGrid = 1,
    ##     K = K,
    ##     alphaHat_t = alphaHat_t1,
    ##     transMatRate_tc_H = transMatRate_tc_H,
    ##     eMatGrid_t = eMatGrid_t1,
    ##     c = c1,
    ##     0,
    ##     grid_has_read = grid_has_read,
    ##     normalize = TRUE
    ## )
    ## sum(X_alphaHat_t1[, 2])

    ## print("this one")
    ## print(head(alphaHat_t1[, 2]))
    ## print("that one")
    ## print(head(X_alphaHat_t1[, 2]))
    
})
