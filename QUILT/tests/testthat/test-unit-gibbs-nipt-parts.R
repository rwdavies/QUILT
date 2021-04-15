test_that("small Gibbs forwards one step including profiling", {

    speed_test <- TRUE

    if (speed_test) {
        ## adjust to see differences
        K <- 1000
        nSNPs <- 5000
        nReads <- 100
    } else {
        ## ok, this seems to make sense
        ## simulate some stuff
        K <- 50
        nSNPs <- 200
        nReads <- 10
    }
    
    out <- make_quilt_fb_test_package(K = K, nReads = nReads, nSNPs = nSNPs, S = 1)
    nGrids <- out$nGrids
    priorCurrent_m <- out$priorCurrent_m
    sampleReads <- out[["sampleReads"]]
    alphaMatCurrent_tc <- out$alphaMatCurrent_tc    
    alphaMatCurrent_tc[] <- 1 / K

    eMatRead_t1 <- out[["list_of_eMatRead_t"]][[1]]
    eMatGrid_t1 <- out[["list_of_eMatGrid_t"]][[1]]

    transMatRate_tc_H <- out$transMatRate_tc_H
    betaHat_t1 <- array(1, c(K, nGrids))
    wif0 <- as.integer(sapply(out$sampleReads, function(x) x[[2]]))
    grid_has_read <- rep(FALSE, nGrids)
    grid_has_read[wif0 + 1] <- TRUE
    
    s <- 0
    minus_log_c_sum <- 0

    ##
    ## first, check things have not changed, including re-allocating reads
    ## 
    for(i_check in 1:2) {
        ##
        eMatRead_t1 <- out[["list_of_eMatRead_t"]][[1]]
        eMatGrid_t1 <- out[["list_of_eMatGrid_t"]][[1]]
        alphaHat_t1 <- array(1, c(K, nGrids))
        c1 <- array(1, nGrids)
        Rcpp_run_forward_haploid(
            alphaHat_t1,
            c1,
            eMatGrid_t1,
            alphaMatCurrent_tc,
            transMatRate_tc_H,
            priorCurrent_m,
            s
        )
        alphaHat_t1_new <- array(1, c(K, nGrids))
        alphaHat_t1_new[] <- alphaHat_t1[]
        c1_new <- array(1, nGrids)
        c1_new[] <- c1[]
        rcpp_reinitialize_in_iterations(
            s = 0,
            alphaHat_t = alphaHat_t1_new,
            c = c1_new,
            priorCurrent_m = priorCurrent_m,
            eMatGrid_t = eMatGrid_t1,
            K = K
        )
        for(iGrid0 in 1:(nGrids - 1)) {
            if (grid_has_read[iGrid0 + 1]) {
                iRead <- sample(which(wif0 == (iGrid0)), 1)
                if (runif(1) < 0.5) {
                    ## hmm, easiest is just to remove this read, that's the same
                    eMatGrid_t1[, iGrid0 + 1] <- eMatGrid_t1[, iGrid0 + 1] / eMatRead_t1[, iRead]
                    alphaHat_t1_new[, iGrid0 + 1] <- alphaHat_t1_new[, iGrid0 + 1] / eMatRead_t1[, iRead]
                    ## change
                    alphaConst <- 1 / sum(alphaHat_t1_new[, iGrid0 + 1])
                    c1_new[iGrid0 + 1] <- c1_new[iGrid0 + 1] * alphaConst
                    alphaHat_t1_new[, iGrid0 + 1] <- alphaHat_t1_new[, iGrid0 + 1] * alphaConst
                }
            }
            if (i_check == 1) {
                rcpp_alpha_forward_one(
                    s = 0,
                    iGrid = iGrid0,
                    K = K,
                    alphaHat_t = alphaHat_t1_new,
                    transMatRate_tc_H = transMatRate_tc_H,
                    eMatGrid_t = eMatGrid_t1,
                    alphaMatCurrent_tc = alphaMatCurrent_tc,
                    c = c1_new,
                    minus_log_c_sum = minus_log_c_sum,
                    normalize = TRUE
                )
            } else {
                rcpp_alpha_forward_one_QUILT_faster(
                    s = 0,
                    iGrid = iGrid0,
                    K = K,
                    alphaHat_t = alphaHat_t1_new,
                    transMatRate_tc_H = transMatRate_tc_H,
                    eMatGrid_t = eMatGrid_t1,
                    c = c1_new,
                    grid_has_read = grid_has_read,
                    minus_log_c_sum = minus_log_c_sum,
                    normalize = TRUE
                )
            }
        }
        ## now re-do
        alphaHat_t1_redone <- array(1, c(K, nGrids))
        c1_redone <- array(1, nGrids)
        Rcpp_run_forward_haploid(
            alphaHat_t1_redone,
            c1_redone,
            eMatGrid_t1,
            alphaMatCurrent_tc,
            transMatRate_tc_H,
            priorCurrent_m,
            s
        )
        expect_equal(alphaHat_t1_redone, alphaHat_t1_new)
        expect_equal(c1_redone, c1_new, tol = 1e-4)
    }


    ##
    ## finally, optional speed test
    ##

    f_original <- function() {
        for(iGrid in 1:(nGrids - 1)) {
            rcpp_alpha_forward_one(
                s = 0,
                iGrid = iGrid,
                K = K,
                alphaHat_t = alphaHat_t1,
                transMatRate_tc_H = transMatRate_tc_H,
                eMatGrid_t = eMatGrid_t1,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                c = c1,
                minus_log_c_sum = minus_log_c_sum,
                normalize = TRUE
            )
        }
    }
    f_modified <- function() {
        for(iGrid in 1:(nGrids - 1)) {
            rcpp_alpha_forward_one_QUILT_faster(
                s = 0,
                iGrid = iGrid,
                K = K,
                alphaHat_t = alphaHat_t1,
                transMatRate_tc_H = transMatRate_tc_H,
                eMatGrid_t = eMatGrid_t1,
                c = c1,
                minus_log_c_sum = minus_log_c_sum,
                grid_has_read = grid_has_read,                
                normalize = TRUE
            )
        }
    }

        
    if (speed_test) {
        ## compare against each other        
        print(microbenchmark::microbenchmark(
            f_original(),
            f_modified(),
            times = 20
        ))
    }

    

})


test_that("faster version of small Gibbs backwards works", {

    ## toggle on and off to check comparison of speed
    speed_test <- FALSE

    if (speed_test) {
        ## adjust to see differences
        K <- 1000
        nSNPs <- 5000
        nReads <- 100
    } else {
        ## ok, this seems to make sense
        ## simulate some stuff
        K <- 50
        nSNPs <- 200
        nReads <- 10
    }
    
    out <- make_quilt_fb_test_package(K = K, nReads = nReads, nSNPs = nSNPs)
    nGrids <- out$nGrids
    priorCurrent_m <- out$priorCurrent_m
    alphaMatCurrent_tc <- out$alphaMatCurrent_tc    
    alphaMatCurrent_tc[] <- 1 / K
    
    eMatGrid_t1 <- out[["list_of_eMatGrid_t"]][[1]]
    c1 <- array(1, nGrids)

    transMatRate_tc_H <- out$transMatRate_tc_H
    alphaHat_t1 <- array(1, c(K, nGrids))    
    betaHat_t1 <- array(1, c(K, nGrids))
    wif0 <- as.integer(sapply(out$sampleReads, function(x) x[[2]]))
    grid_has_read <- rep(FALSE, nGrids)
    grid_has_read[wif0 + 1] <- TRUE
    
    s <- 0
    Rcpp_run_forward_haploid(
        alphaHat_t1, c1,
        eMatGrid_t1,
        alphaMatCurrent_tc,
        transMatRate_tc_H,
        priorCurrent_m,
        s
    )
    
    f_original <- function() {
        Rcpp_run_backward_haploid(
            betaHat_t1,
            c1,
            eMatGrid_t1,
            alphaMatCurrent_tc,
            transMatRate_tc_H,
            s
        )
    }
    f_modified <- function() {
        Rcpp_run_backward_haploid_QUILT_faster(
            betaHat_t1,
            c1,
            eMatGrid_t1,
            transMatRate_tc_H,
            grid_has_read,
            s
        )
    }
    
    ## original
    betaHat_t1 <- array(1, c(K, nGrids))
    f_original()
    betaHat_t1_original <- betaHat_t1
    ## modified
    betaHat_t1 <- array(1, c(K, nGrids))
    f_modified()
    betaHat_t1_modified <- betaHat_t1
    ## 
    expect_equal(betaHat_t1_modified, betaHat_t1_original)

    if (speed_test) {
        ## compare against each other        
        print(microbenchmark::microbenchmark(
            f_original(),
            f_modified(),
            times = 20
        ))
    }

    

})
