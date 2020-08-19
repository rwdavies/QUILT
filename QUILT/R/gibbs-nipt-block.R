helper_block_gibbs_resampler <- function(
    H,
    blocked_snps,
    grid,
    wif0,
    ff,
    s,
    eHapsCurrent_tc,
    alphaMatCurrent_tc,
    priorCurrent_m,
    transMatRate_tc_H,
    sampleReads,
    L_grid,
    shuffle_bin_radius,
    smooth_cm,
    maxDifferenceBetweenReads = 1000,
    Jmax = 1000,
    language = "Rcpp",
    verbose = FALSE,
    do_checks = FALSE,
    block_approach = 1,
    block_gibbs_quantile_prob = 0.9,
    class_sum_cutoff = 0.06,
    use_smooth_cm_in_block_gibbs = FALSE
) {

    eMatRead_t <- array(1, c(dim(eHapsCurrent_tc)[1], length(sampleReads)))
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = eHapsCurrent_tc,
        s = s - 1,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = Jmax,
        eMatHapOri_t = array(0, c(1, 1)),
        pRgivenH1 = array(0),
        pRgivenH2 = array(0),
        prev = 0,
        suppressOutput = 1,
        prev_section = "",
        next_section = "",
        run_pseudo_haploid = FALSE,
        rescale_eMatRead_t = FALSE
    )
    ##
    fpp_stuff <- list(
        transMatRate_tc_H = transMatRate_tc_H,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        priorCurrent_m = priorCurrent_m ,
        eMatRead_t = eMatRead_t,
        s = s,
        sampleReads = sampleReads
    )
    initial_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)
    alphaHat_t1 <- initial_package[[1]][["alphaHat_t"]]
    alphaHat_t2 <- initial_package[[2]][["alphaHat_t"]]
    alphaHat_t3 <- initial_package[[3]][["alphaHat_t"]]
    betaHat_t1 <- initial_package[[1]][["betaHat_t"]]
    betaHat_t2 <- initial_package[[2]][["betaHat_t"]]
    betaHat_t3 <- initial_package[[3]][["betaHat_t"]]
    c1 <- initial_package[[1]][["c"]]
    c2 <- initial_package[[2]][["c"]]
    c3 <- initial_package[[3]][["c"]]
    eMatGrid_t1 <- initial_package[[1]][["eMatGrid_t"]]
    eMatGrid_t2 <- initial_package[[2]][["eMatGrid_t"]]
    eMatGrid_t3 <- initial_package[[3]][["eMatGrid_t"]]
    ##
    ##
    ##

    ##
    ## check attempt to identify breaks
    ##
    if (language == "R") {
        f <- R_define_blocked_snps_using_gamma_on_the_fly
    } else {
        f <- Rcpp_define_blocked_snps_using_gamma_on_the_fly
        s <- s - 1L
    }
    out <- f(
        alphaHat_t1 = alphaHat_t1,
        alphaHat_t2 = alphaHat_t2,
        alphaHat_t3 = alphaHat_t3,
        betaHat_t1 = betaHat_t1,
        betaHat_t2 = betaHat_t2,
        betaHat_t3 = betaHat_t3,
        c1 = c1,
        c2 = c2,
        c3 = c3,
        eMatGrid_t1 = eMatGrid_t1,
        eMatGrid_t2 = eMatGrid_t2,
        eMatGrid_t3 = eMatGrid_t3,
        transMatRate_tc_H = transMatRate_tc_H,
        shuffle_bin_radius = shuffle_bin_radius,
        L_grid = L_grid,
        grid = grid,
        s = s,
        block_gibbs_quantile_prob = block_gibbs_quantile_prob,
        verbose = verbose,
        use_smooth_cm_in_block_gibbs = use_smooth_cm_in_block_gibbs,
        smooth_cm = smooth_cm
    )
    attempted_blocked_snps <- out[["blocked_snps"]]
    if (language == "R") {
        ## 
    } else {
        f <- Rcpp_define_blocked_snps_using_gamma_on_the_fly
        s <- s + 1L
    }
    
    

    
    H_class <- calculate_H_class(
        eMatRead_t = eMatRead_t,
        alphaHat_t1 = alphaHat_t1,
        alphaHat_t2 = alphaHat_t2,
        alphaHat_t3 = alphaHat_t3,
        betaHat_t1 = betaHat_t1,
        betaHat_t2 = betaHat_t2,
        betaHat_t3 = betaHat_t3,
        ff = ff,
        wif0 = wif0,
        H = H,
        class_sum_cutoff = class_sum_cutoff
    )
    ##
    ## 
    ##
    n_blocks <- max(blocked_snps) + 1
    runif_block <- runif(n_blocks)    
    runif_total <- runif(n_blocks)
    nReads <- length(sampleReads)
    runif_proposed <- matrix(runif(nReads * 6), nrow = 6, ncol = nReads)
    if (language == "R") {
        f <- R_block_gibbs_resampler
        use_cpp_bits_in_R <- FALSE
    } else if (language == "R_with_Rcpp") {
        f <- R_block_gibbs_resampler
        use_cpp_bits_in_R <- TRUE
    } else if (language == "Rcpp") {
        f <- Rcpp_block_gibbs_resampler
        s <- s - 1L
        use_cpp_bits_in_R <- FALSE
    } else {
        stop("not a language!")
    }

    block_out <- f(
        alphaHat_t1 = alphaHat_t1,
        alphaHat_t2 = alphaHat_t2,
        alphaHat_t3 = alphaHat_t3,
        betaHat_t1 = betaHat_t1,
        betaHat_t2 = betaHat_t2,
        betaHat_t3 = betaHat_t3,
        c1 = c1,
        c2 = c2,
        c3 = c3,
        eMatGrid_t1 = eMatGrid_t1,
        eMatGrid_t2 = eMatGrid_t2,
        eMatGrid_t3 = eMatGrid_t3,
        H = H,
        H_class = H_class,
        eMatRead_t = eMatRead_t,
        blocked_snps = blocked_snps,
        runif_block = runif_block,
        runif_total = runif_total,
        runif_proposed = runif_proposed,
        grid = grid,
        wif0 = wif0,
        ff = ff,
        s = s,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        priorCurrent_m = priorCurrent_m,
        transMatRate_tc_H = transMatRate_tc_H,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = Jmax,
        do_checks = do_checks,
        initial_package = initial_package,
        verbose = verbose,
        fpp_stuff = fpp_stuff,
        use_cpp_bits_in_R = use_cpp_bits_in_R,
        block_approach = block_approach
    )

    if (language == "Rcpp") {
        s <- s + 1L
        ## need to package together
        block_out <- append(block_out, list(
            eMatGrid_t1 = eMatGrid_t1,
            eMatGrid_t2 = eMatGrid_t2,
            eMatGrid_t3 = eMatGrid_t3,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            H = H,
            alphaHat_t1 = alphaHat_t1,
            alphaHat_t2 = alphaHat_t2,
            alphaHat_t3 = alphaHat_t3,            
            betaHat_t1 = betaHat_t1,
            betaHat_t2 = betaHat_t2,
            betaHat_t3 = betaHat_t3
        ))
    } 
    
    ## do checks here
    final_package <- for_testing_get_full_package_probabilities(block_out[["H"]], fpp_stuff)
    expect_equal(block_out[["eMatGrid_t1"]], final_package[[1]][["eMatGrid_t"]])
    expect_equal(block_out[["eMatGrid_t2"]], final_package[[2]][["eMatGrid_t"]])
    expect_equal(block_out[["eMatGrid_t3"]], final_package[[3]][["eMatGrid_t"]])
    expect_equal(block_out[["c1"]], final_package[[1]][["c"]])
    expect_equal(block_out[["c2"]], final_package[[2]][["c"]])
    expect_equal(block_out[["c3"]], final_package[[3]][["c"]])
    expect_equal(block_out[["alphaHat_t1"]], final_package[[1]][["alphaHat_t"]])
    expect_equal(block_out[["alphaHat_t2"]], final_package[[2]][["alphaHat_t"]])
    expect_equal(block_out[["alphaHat_t3"]], final_package[[3]][["alphaHat_t"]])
    expect_equal(block_out[["betaHat_t1"]], final_package[[1]][["betaHat_t"]])
    expect_equal(block_out[["betaHat_t2"]], final_package[[2]][["betaHat_t"]])
    expect_equal(block_out[["betaHat_t3"]], final_package[[3]][["betaHat_t"]])
    ##
    block_out <- append(block_out, list(attempted_blocked_snps = attempted_blocked_snps))
    return(block_out)
}



R_gibbs_block_forward_one_master <- function(
    approach2_iRead,
    iGrid,
    s,
    alphaStore,
    log_cStore,
    rr,
    rr0,
    eMatGridLocal,
    eMatGridLocalc,
    transMatRate_tc_H,
    alphaMatCurrent_tc,
    priorCurrent_m,
    read_is_uninformative,
    block_approach,
    wif0,
    eMatRead_t,
    nReads,
    H,
    proposed_H,
    H_class,
    rlc,
    rlcM,
    runif_proposed,
    use_cpp_bits_in_R,
    do_checks,
    all_packages,
    cur_package,
    fpp_stuff,
    read_end_0_based
) {
    if (use_cpp_bits_in_R) {
        f <- Rcpp_gibbs_block_forward_one
        iGrid <- iGrid - 1
        s <- s - 1
        approach2_iRead[1] <- approach2_iRead[1] - 1L
    } else {
        f <- R_gibbs_block_forward_one
    }
    out <- f(
        approach2_iRead = approach2_iRead,
        iGrid = iGrid,
        s = s,
        alphaStore = alphaStore,
        log_cStore = log_cStore,
        rr = rr,
        rr0 = rr0,
        eMatGridLocal = eMatGridLocal,
        eMatGridLocalc = eMatGridLocalc,
        transMatRate_tc_H = transMatRate_tc_H,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        priorCurrent_m = priorCurrent_m,
        read_is_uninformative = read_is_uninformative,
        block_approach = block_approach,
        wif0 = wif0,
        eMatRead_t = eMatRead_t,
        nReads = nReads,
        H = H,
        proposed_H = proposed_H,
        H_class = H_class,
        rlc = rlc,
        rlcM = rlcM,
        runif_proposed = runif_proposed
    )
    if (use_cpp_bits_in_R) {
        iGrid <- iGrid + 1
        s <- s + 1
        approach2_iRead[1] <- approach2_iRead[1] + 1L
    } else {
        alphaStore <- out[["alphaStore"]]
        log_cStore <- out[["log_cStore"]]
        approach2_iRead <- out[["approach2_iRead"]]
        proposed_H <- out[["proposed_H"]]
        eMatGridLocalc <- out[["eMatGridLocalc"]]
    }
    ##
    ##
    ##
    if (do_checks && (block_approach == 1)) {
        for(ir in 1:6) {
            for(i in 1:3) {
                h <- rr[ir, i]
                expect_equal(all_packages[[ir]][[h]]$alphaHat_t[, iGrid], alphaStore[, h, ir])
                expect_equal(log(all_packages[[ir]][[h]]$c[iGrid]), log_cStore[iGrid, h, ir])
            }
        }
    }
    if (do_checks & (block_approach == 4)) {
        ## if (verbose) {
        ##     print_message("Check alpha for proposed H")
        ## }
        ## check alphaStores
        ## 2 and 5 breaking? after a 3
        rs <- (read_end_0_based + 2)
        re <- approach2_iRead[1] - 1
        w <- rs:re
        for(ir in 1:6) {
            Htemp <- H
            ## ONLY if change makes sense
            if (rs <= re) {
                Htemp[w] <- rr[ir, proposed_H[ir, w]]
            }
            cur_package <- for_testing_get_full_package_probabilities(Htemp, fpp_stuff)
            ## check alphas
            m <- cbind(cur_package[[rr[ir, 1]]]$eMatGrid_t[, iGrid], cur_package[[rr[ir, 2]]]$eMatGrid_t[, iGrid], cur_package[[rr[ir, 3]]]$eMatGrid_t[, iGrid])
            expect_equal(eMatGridLocalc[, , ir], m)
            ## 
            m <- cbind(cur_package[[1]]$alphaHat_t[, iGrid], cur_package[[2]]$alphaHat_t[, iGrid], cur_package[[3]]$alphaHat_t[, iGrid])
            expect_equal(alphaStore[, , ir], m)
        }
    }
    return(
        list(
            approach2_iRead = approach2_iRead,
            cur_package = cur_package,
            alphaStore = alphaStore,
            log_cStore = log_cStore,
            approach2_iRead = approach2_iRead,
            proposed_H = proposed_H,
            eMatGridLocalc = eMatGridLocalc 
        )
    )
}


## so normally, see below, just implement in the middle of other code
## here, for neater testing code, we just give it an initial (possibly messed up) H
## make it build the things it needs
R_block_gibbs_resampler <- function(
    alphaHat_t1,
    alphaHat_t2,
    alphaHat_t3,
    betaHat_t1,
    betaHat_t2,
    betaHat_t3,
    c1,
    c2,
    c3,
    eMatGrid_t1,
    eMatGrid_t2,
    eMatGrid_t3,
    H,
    H_class,
    eMatRead_t,    
    blocked_snps,
    runif_block,
    runif_total,
    runif_proposed,
    grid,
    wif0,
    ff,
    s,
    alphaMatCurrent_tc,
    priorCurrent_m,
    transMatRate_tc_H,
    log_cStore,
    maxDifferenceBetweenReads,
    Jmax,
    do_checks = FALSE,
    initial_package = NULL,
    verbose = FALSE,
    fpp_stuff = NULL,
    use_cpp_bits_in_R = FALSE,
    block_approach = 4
) {

    ##
    ##
    K <- dim(alphaMatCurrent_tc)[1]
    nSNPs <- length(grid)
    nGrids <- dim(alphaMatCurrent_tc)[2] + 1    
    S <- dim(alphaMatCurrent_tc)[3]
    nReads <- length(wif0)
    prior_probs <- c(0.5, (1 - ff) / 2, ff / 2)
    log_prior_probs <- log(prior_probs)
    read_is_uninformative <- logical(1)
    read_is_uninformative[] <- as.logical(NA)
    ## read_is_uninformative = NA,
    if (is.na(read_is_uninformative[1])) {
        read_is_uninformative <- array(FALSE, 2 * nReads)
        read_is_uninformative[seq(1, nReads, 2)] <- TRUE
    }
    ##
    if (use_cpp_bits_in_R) {
        f <- Rcpp_make_gibbs_considers
    } else {
        f <- R_make_gibbs_considers
    }
    out <- f(
        blocked_snps = blocked_snps,
        grid = grid,
        wif0 = wif0,
        nGrids = nGrids
    )
    ## yes, take out from both
    consider_reads_start_0_based <- out[["consider_reads_start_0_based"]]
    consider_reads_end_0_based  <- out[["consider_reads_end_0_based"]]
    consider_grid_start_0_based <- out[["consider_grid_start_0_based"]]
    consider_grid_end_0_based <- out[["consider_grid_end_0_based"]]
    consider_snp_start_0_based <- out[["consider_snp_start_0_based"]]
    consider_snp_end_0_based <- out[["consider_snp_end_0_based"]]
    consider_grid_where_0_based <- out[["consider_grid_where_0_based"]]
    n_blocks <- out[["n_blocks"]]
    ##
    ## add one check here
    ##
    stopifnot(consider_grid_start_0_based[1] == 0)
    stopifnot(consider_grid_end_0_based[n_blocks] == (nGrids - 1))
    for(i in 1:(length(consider_grid_start_0_based) - 1)) {
        d <- consider_grid_start_0_based[i + 1] - consider_grid_end_0_based[i]
        if (d != 1) {
            stop("bad making of consider grid")
        }
    }
    ##
    ## start, end,
    ## want 6 probabilities, change, sampled
    ## old probability, new probability
    block_results_columns <- c(
        "iBlock",
        "total",
        paste0("p", 1:6),
        "ir_chosen",
        "p_O1_given_H1_L",
        "p_O2_given_H2_L",
        "p_O3_given_H3_L",
        "p_O_given_H_L",
        "p_H_given_L",
        "p_H_given_O_L_up_to_C"
    )
    block_results <- matrix(0.0, nrow = n_blocks * 2, ncol = length(block_results_columns))
    colnames(block_results) <- block_results_columns
    ## 
    ## 
    rr <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(2, 3, 1),
        c(3, 1, 2),
        c(3, 2, 1)
    )
    rr <- matrix(as.integer(rr), ncol = 3)
    rr0 <- rr - 1L
    swap_list <- list(
        list(NA),
        list(c(2, 3)),
        list(c(2, 1)),
        list(c(1, 2), c(2, 3)),
        list(c(1, 2), c(1, 3)),
        list(c(1, 3))
    )
    ##
    p <- prior_probs
    rlc <- make_rlc(ff)
    rlcM <- array(NA, c(3, 6, 7))
    rlcM <- R_fill_rlcM(rlcM, rlc, rr)
    ## also get the reverse!
    ## rrR <- array(NA, c(6, 3))
    ## for(ir in 1:6) {
    ##     for(i in 1:3) {
    ##         rrR[ir, rr[ir, i]] <- i
    ##     }
    ## }
    ## rr means GOES TO that spot
    ## i.e. 4th row (c(2, 3, 1)) means
    ##   label 1 (1st entry) goes to the 2 spot
    ##   label 2 (2nd entry) goes to the 3 spot
    ##   label 3 (3rd entry) goes to the 1 spot
    ##
    if (do_checks) {
        ##
        all_packages <- lapply(1:6, function(ir) {
            iBlock <- 0        
            w <- wif0 %in% (consider_grid_start_0_based[iBlock + 1]:consider_grid_end_0_based[iBlock + 1])
            localH <- H
            localH[w] <- rr[ir, H[w]]
            return(for_testing_get_full_package_probabilities(localH, fpp_stuff))
        })
    } else {
        all_packages <- list(NA)
    }
    ## 
    ##
    logC_before <- numeric(3)
    ## remove first one from after, as we skip over that one
    logC_after <- c(sum(log(c1)), sum(log(c2)), sum(log(c3)))
    ##
    iBlock <- 0
    ever_changed <- as.numeric(0)
    ##
    sum_H <- numeric(3)
    nReads <- integer(1)
    nReads[1] <- length(H)
    ff_0_condition <- ff > 0
    ##
    ## begin of serious weirdness
    ##
    ## did not work sept 20 2019
    ## below, on second go after a pass through Rcpp, "0" was not working through seq or 0:X for loops including 0
    ## does work as an integer
    ##
    ## for(iRead in seq(0, nReads - 1)) {    
    for(iRead in seq(0L, as.integer(nReads - 1))) {
        sum_H[H[iRead + 1]] <- sum_H[H[iRead + 1]] + 1.0
    }
    ##
    ## end of serious weirdness
    ##
    ## H, and proposed_H, are ALWAYS 1-based
    proposed_H <- array(as.integer(-1), c(6, nReads)) 
    approach2_iRead <- integer(1)
    approach2_iRead[] <- 1L
    nReads <- length(H)
    ##
    ## initialize with nothing
    ##
    alphaStore <- array(0, c(K, 3, 6)) ## only do one of these. redo later
    alphaHatLocal <- betaHatLocal <- eMatGridLocal <- array(0, c(K, 3))
    eMatGridLocalc <- array(0, c(K, 3, 6))
    log_cStore <- array(0, c(nGrids, 3, 6))
    read_end_0_based <- -1 ## only matters for checks
    cur_package <- NULL
    ##
    ## loop!
    ##
    if (class(H[1]) != "integer") {
        stop("somehow H became non-integer before even starting")
    }
    
    for(iGrid in 1:nGrids) {    
        if (verbose) {
            print_message(paste0("iBlock = ", iBlock, ", iGrid = ", iGrid))
        }
        ##
        ## now, just go forward for each one!
        ##
        eMatGridLocal[, 1] <- eMatGrid_t1[, iGrid]
        eMatGridLocal[, 2] <- eMatGrid_t2[, iGrid]
        eMatGridLocal[, 3] <- eMatGrid_t3[, iGrid]
        ##
        out <- R_gibbs_block_forward_one_master(approach2_iRead = approach2_iRead, iGrid =iGrid, s = s, alphaStore = alphaStore, log_cStore = log_cStore, rr = rr, rr0 = rr0, eMatGridLocal = eMatGridLocal, eMatGridLocalc = eMatGridLocalc, transMatRate_tc_H = transMatRate_tc_H, alphaMatCurrent_tc = alphaMatCurrent_tc, priorCurrent_m = priorCurrent_m, read_is_uninformative = read_is_uninformative, block_approach = block_approach, wif0 = wif0, eMatRead_t = eMatRead_t, nReads = nReads, H = H, proposed_H = proposed_H, H_class = H_class, rlc = rlc, rlcM = rlcM, runif_proposed = runif_proposed, use_cpp_bits_in_R = use_cpp_bits_in_R, do_checks = do_checks, all_packages = all_packages, cur_package = cur_package, read_end_0_based = read_end_0_based, fpp_stuff = fpp_stuff)
        approach2_iRead <- out[["approach2_iRead"]]
        cur_package <- out[["cur_package"]]
        alphaStore <- out[["alphaStore"]]
        log_cStore <- out[["log_cStore"]]
        approach2_iRead <- out[["approach2_iRead"]]
        proposed_H <- out[["proposed_H"]]
        eMatGridLocalc <- out[["eMatGridLocalc"]]
        
        ##
        ## if this is when we consider the block change
        ##
        if (consider_grid_where_0_based[iGrid] > (- 1)) {
            iBlock <- consider_grid_where_0_based[iGrid]
            grid_start_0_based <- consider_grid_start_0_based[iBlock + 1]
            grid_end_0_based <- consider_grid_end_0_based[iBlock + 1]
            read_start_0_based <- consider_reads_start_0_based[iBlock + 1]
            read_end_0_based <- consider_reads_end_0_based[iBlock + 1]
            ## 
            ## now, build probabilities!
            ## need to use alphaBeta here!
            ##
            ## block relabelling
            ##
            if (do_checks & (block_approach == 1)) {
                if (0 < iBlock) {
                    if (verbose) {
                        print_message("Check everything before block relabelling")
                    }
                    cur_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)
                    ## here, only care about up to previous region
                    w <- 1:(consider_grid_end_0_based[iBlock + 1 - 1])
                    expect_equal(eMatGrid_t1[, w], cur_package[[1]][["eMatGrid_t"]][, w])
                    expect_equal(eMatGrid_t2[, w], cur_package[[2]][["eMatGrid_t"]][, w])
                    expect_equal(eMatGrid_t3[, w], cur_package[[3]][["eMatGrid_t"]][, w])
                    expect_equal(c1[w], cur_package[[1]][["c"]][w])
                    expect_equal(c2[w], cur_package[[2]][["c"]][w])
                    expect_equal(c3[w], cur_package[[3]][["c"]][w])
                    expect_equal(alphaHat_t1[, w], cur_package[[1]]$alphaHat_t[, w])
                    expect_equal(alphaHat_t2[, w], cur_package[[2]]$alphaHat_t[, w])
                    expect_equal(alphaHat_t3[, w], cur_package[[3]]$alphaHat_t[, w])
                }
            }
            if (use_cpp_bits_in_R) {
                f <- Rcpp_consider_block_relabelling
                s <- s - 1
                iGrid <- iGrid - 1
            } else {
                f <- R_consider_block_relabelling
            }
            if (class(H[1]) != "integer") {
                stop("somehow H became non-integer!")
            }
            out <- f(
                iBlock = iBlock,
                runif_block = runif_block,
                sum_H = sum_H,
                s = s,
                rr = rr,
                rr0 = rr0,
                ff = ff,
                log_prior_probs = log_prior_probs,
                logC_before = logC_before,
                logC_after = logC_after,
                verbose = verbose,
                swap_list = swap_list,
                eMatGridLocal = eMatGridLocal,
                betaHatLocal = betaHatLocal,
                iGrid = iGrid,
                grid_start_0_based = grid_start_0_based,
                grid_end_0_based = grid_end_0_based,
                read_start_0_based = read_start_0_based,
                read_end_0_based = read_end_0_based,
                wif0 = wif0,
                log_cStore = log_cStore,
                alphaStore = alphaStore,
                read_is_uninformative = read_is_uninformative,
                block_approach = block_approach,
                do_checks = do_checks,
                all_packages = all_packages,
                block_results = block_results,
                ever_changed = ever_changed,
                transMatRate_tc_H = transMatRate_tc_H,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                priorCurrent_m = priorCurrent_m,
                fpp_stuff = fpp_stuff,
                H = H,
                proposed_H = proposed_H,
                nReads = nReads,
                eMatRead_t = eMatRead_t,
                alphaHat_t1 = alphaHat_t1, betaHat_t1 = betaHat_t1, c1 = c1, eMatGrid_t1 = eMatGrid_t1,
                alphaHat_t2 = alphaHat_t2, betaHat_t2 = betaHat_t2, c2 = c2, eMatGrid_t2 = eMatGrid_t2,
                alphaHat_t3 = alphaHat_t3, betaHat_t3 = betaHat_t3, c3 = c3, eMatGrid_t3 = eMatGrid_t3
            )
            ##
            if (use_cpp_bits_in_R) {
                s <- s + 1
                iGrid <- iGrid + 1
            } else {
                eMatGrid_t1 <- out$eMatGrid_t1
                eMatGrid_t2 <- out$eMatGrid_t2
                eMatGrid_t3 <- out$eMatGrid_t3
                c1 <- out$c1
                c2 <- out$c2
                c3 <- out$c3
                alphaHat_t1 <- out$alphaHat_t1
                alphaHat_t2 <- out$alphaHat_t2
                alphaHat_t3 <- out$alphaHat_t3            
                betaHat_t1 <- out$betaHat_t1
                betaHat_t2 <- out$betaHat_t2
                betaHat_t3 <- out$betaHat_t3
                H <- out$H
                ever_changed <- out$ever_changed
                block_results <- out$block_results
                sum_H <- out$sum_H
            }
            if (do_checks) {
                if (verbose) {
                    print_message("Check everything after block relabelling, and before total relabelling")
                }
                cur_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)
                expect_equal(eMatGrid_t1, cur_package[[1]][["eMatGrid_t"]])
                expect_equal(eMatGrid_t2, cur_package[[2]][["eMatGrid_t"]])
                expect_equal(eMatGrid_t3, cur_package[[3]][["eMatGrid_t"]])
                expect_equal(c1[1:iGrid], cur_package[[1]][["c"]][1:iGrid])
                expect_equal(c2[1:iGrid], cur_package[[2]][["c"]][1:iGrid])
                expect_equal(c3[1:iGrid], cur_package[[3]][["c"]][1:iGrid])
                expect_equal(alphaHat_t1[, 1:iGrid], cur_package[[1]]$alphaHat_t[, 1:iGrid])
                expect_equal(alphaHat_t2[, 1:iGrid], cur_package[[2]]$alphaHat_t[, 1:iGrid])
                expect_equal(alphaHat_t3[, 1:iGrid], cur_package[[3]]$alphaHat_t[, 1:iGrid])
            }
            ##
            ## total relabelling
            ##
            ## weird results if use this on the fly, maybe due toe pass by reference? but hard to nail down
            if (ff_0_condition) {                 ## only do if ff > 0
                if (use_cpp_bits_in_R) {
                    f <- Rcpp_consider_total_relabelling
                } else {
                    f <- R_consider_total_relabelling
                }
                out <- f(
                    iBlock = iBlock,
                    rr = rr,
                    rr0 = rr0,
                    ff = ff,
                    log_prior_probs = log_prior_probs,
                    logC_before = logC_before,
                    logC_after = logC_after,
                    verbose = verbose,
                    swap_list = swap_list,
                    block_results = block_results,
                    runif_total = runif_total,
                    sum_H = sum_H,
                    H = H,
                    alphaHat_t1 = alphaHat_t1, betaHat_t1 = betaHat_t1, c1 = c1, eMatGrid_t1 = eMatGrid_t1,
                    alphaHat_t2 = alphaHat_t2, betaHat_t2 = betaHat_t2, c2 = c2, eMatGrid_t2 = eMatGrid_t2,
                    alphaHat_t3 = alphaHat_t3, betaHat_t3 = betaHat_t3, c3 = c3, eMatGrid_t3 = eMatGrid_t3
                )
                if (!use_cpp_bits_in_R) {
                    eMatGrid_t1 <- out$eMatGrid_t1
                    eMatGrid_t2 <- out$eMatGrid_t2
                    eMatGrid_t3 <- out$eMatGrid_t3
                    c1 <- out$c1
                    c2 <- out$c2
                    c3 <- out$c3
                    alphaHat_t1 <- out$alphaHat_t1
                    alphaHat_t2 <- out$alphaHat_t2
                    alphaHat_t3 <- out$alphaHat_t3            
                    betaHat_t1 <- out$betaHat_t1
                    betaHat_t2 <- out$betaHat_t2
                    betaHat_t3 <- out$betaHat_t3
                    H <- out$H
                    sum_H <- out$sum_H
                    logC_before <- out$logC_before
                    logC_after <- out$logC_after
                    block_results <- out$block_results
                }
                if (do_checks) {
                    if (verbose) {
                        print_message("Check everything after total relabelling for new H")
                    }
                    cur_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)
                    expect_equal(eMatGrid_t1, cur_package[[1]][["eMatGrid_t"]])
                    expect_equal(eMatGrid_t2, cur_package[[2]][["eMatGrid_t"]])
                    expect_equal(eMatGrid_t3, cur_package[[3]][["eMatGrid_t"]])
                    expect_equal(c1[1:iGrid], cur_package[[1]][["c"]][1:iGrid])
                    expect_equal(c2[1:iGrid], cur_package[[2]][["c"]][1:iGrid])
                    expect_equal(c3[1:iGrid], cur_package[[3]][["c"]][1:iGrid])
                    expect_equal(alphaHat_t1[, 1:iGrid], cur_package[[1]]$alphaHat_t[, 1:iGrid])
                    expect_equal(alphaHat_t2[, 1:iGrid], cur_package[[2]]$alphaHat_t[, 1:iGrid])
                    expect_equal(alphaHat_t3[, 1:iGrid], cur_package[[3]]$alphaHat_t[, 1:iGrid])
                    expect_equal(sum(H == 1), sum_H[1])
                    expect_equal(sum(H == 2), sum_H[2])
                    expect_equal(sum(H == 3), sum_H[3])
                }
            }
            ##
            ## if this is NOT the last block, reset, and checks
            ##
            if (((iBlock + 1) + 1) <= length(consider_grid_start_0_based)) {
                ## note above - iBlock + 1, first one, is 0 vs 1-based. second one is need to check future
                if (do_checks) {
                    ## reset for next time
                    ## iBlock <- 0 ## here, set!
                    if (verbose) {
                        print_message("Reset all packages")
                    }
                    all_packages <- lapply(1:6, function(ir) {
                        w <- wif0 %in% (consider_grid_start_0_based[(iBlock + 1) + 1]:consider_grid_end_0_based[(iBlock + 1) + 1])
                        localH <- H
                        localH[w] <- rr[ir, H[w]]
                        return(for_testing_get_full_package_probabilities(localH, fpp_stuff))
                    })
                } else {
                    all_packages <- NULL
                }
                if (use_cpp_bits_in_R) {
                    f <- Rcpp_reset_local_variables
                    iGrid <- iGrid - 1
                } else {
                    f <- R_reset_local_variables
                }
                out <- f(
                    iGrid = iGrid,
                    verbose = verbose,
                    alphaHatLocal = alphaHatLocal,
                    alphaStore = alphaStore,
                    alphaHat_t1 = alphaHat_t1,
                    alphaHat_t2 = alphaHat_t2,
                    alphaHat_t3 = alphaHat_t3,
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    log_cStore = log_cStore
                )
                if (use_cpp_bits_in_R) {
                    iGrid <- iGrid + 1
                } else {
                    alphaStore <- out[["alphaStore"]]
                    log_cStore <- out[["log_cStore"]]
                }
                if (do_checks) {
                    if (verbose) {
                        print_message("Check alphaStore, log_cStore after reset")
                    }
                    ## now once reset, should be the same as for current H
                    ## but here, no relabelling has yet taken place
                    ## so have "i" on the interior. they are all the same anyway
                    ## (the alpha with the all packages)
                    ## print(paste0("check iGrid=", iGrid, ", h = ", h, ", ir = ", ir, " i = ", i))
                    for(ir in 1:6) {
                        for(i in 1:3) {
                            h <- rr[ir, i]
                            expect_equal(all_packages[[ir]][[h]]$alphaHat_t[, iGrid], alphaStore[, h, ir])
                            expect_equal(log(all_packages[[ir]][[h]]$c[iGrid]), log_cStore[iGrid, h, ir])
                        }
                    }
                }
            }
            ##
            ## add to logC from previous run
            ##
            for(iGrid2 in grid_start_0_based:grid_end_0_based) {
                logC_before[1] <- logC_before[1] + log(c1[iGrid2 + 1])
                logC_before[2] <- logC_before[2] + log(c2[iGrid2 + 1])
                logC_before[3] <- logC_before[3] + log(c3[iGrid2 + 1])
            }
        }
        ##
        ## remove from this for now, this is the remaining C
        ##
        logC_after[1] <- logC_after[1] - log(c1[iGrid])
        logC_after[2] <- logC_after[2] - log(c2[iGrid])
        logC_after[3] <- logC_after[3] - log(c3[iGrid])
        ## 
    }
    
    ##
    ## re-run backward
    ## 
    betaHat_t1[, nGrids] <- c1[nGrids]
    Rcpp_run_backward_haploid(        
        betaHat_t1,
        c = c1,
        eMatGrid_t = eMatGrid_t1,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc = transMatRate_tc_H,
        s = s - 1
    )
    betaHat_t2[, nGrids] <- c2[nGrids]
    Rcpp_run_backward_haploid(        
        betaHat_t2,
        c = c2,
        eMatGrid_t = eMatGrid_t2,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc = transMatRate_tc_H,
        s = s - 1
    )
    betaHat_t3[, nGrids] <- c3[nGrids]
    Rcpp_run_backward_haploid(        
        betaHat_t3,
        c = c3,
        eMatGrid_t = eMatGrid_t3,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc = transMatRate_tc_H,
        s = s - 1
    )
    return(
        list(
            eMatGrid_t1 = eMatGrid_t1,
            eMatGrid_t2 = eMatGrid_t2,
            eMatGrid_t3 = eMatGrid_t3,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            alphaHat_t1 = alphaHat_t1,
            alphaHat_t2 = alphaHat_t2,
            alphaHat_t3 = alphaHat_t3,            
            betaHat_t1 = betaHat_t1,
            betaHat_t2 = betaHat_t2,
            betaHat_t3 = betaHat_t3,
            H = H,
            fpp_stuff = fpp_stuff,
            block_results = block_results,
            consider_snp_start_0_based = consider_snp_start_0_based,            
            consider_snp_end_0_based = consider_snp_end_0_based,
            consider_reads_start_0_based = consider_reads_start_0_based,
            consider_reads_end_0_based = consider_reads_end_0_based,
            alphaStore = alphaStore,
            log_cStore = log_cStore
        )
    )
}




R_define_blocked_snps_using_gamma_on_the_fly <- function(
    alphaHat_t1,
    alphaHat_t2,
    alphaHat_t3,
    betaHat_t1,
    betaHat_t2,
    betaHat_t3,
    c1,
    c2,
    c3,
    eMatGrid_t1,
    eMatGrid_t2,
    eMatGrid_t3,
    transMatRate_tc_H,
    shuffle_bin_radius,
    L_grid,
    grid,
    s,
    use_smooth_cm_in_block_gibbs,
    smooth_cm,
    block_gibbs_quantile_prob = 0.9,
    verbose = FALSE
) {
    nSNPs <- length(grid)
    nGrids <- length(c1)
    ## ## absolute difference between gamma
    ## diff <- array(0, c(3, nGrids - 1))
    ## g1a <- (alphaHat_t1[, 0 + 1] * betaHat_t1[, 0 + 1]) / c1[0 + 1]
    ## g2a <- (alphaHat_t2[, 0 + 1] * betaHat_t2[, 0 + 1]) / c2[0 + 1]
    ## g3a <- (alphaHat_t3[, 0 + 1] * betaHat_t3[, 0 + 1]) / c3[0 + 1]
    ## for(iGrid in 1:(nGrids - 1)) {
    ##     g1b <- (alphaHat_t1[, iGrid + 1] * betaHat_t1[, iGrid + 1]) / c1[iGrid + 1]
    ##     g2b <- (alphaHat_t2[, iGrid + 1] * betaHat_t2[, iGrid + 1]) / c2[iGrid + 1]
    ##     g3b <- (alphaHat_t3[, iGrid + 1] * betaHat_t3[, iGrid + 1]) / c3[iGrid + 1]
    ##     diff[1, iGrid - 1 + 1] <- sum( (g1b - g1a) ** 2)
    ##     diff[2, iGrid - 1 + 1] <- sum( (g2b - g2a) ** 2)
    ##     diff[3, iGrid - 1 + 1] <- sum( (g3b - g3a) ** 2)
    ##     g1a <- g1b
    ##     g2a <- g2b
    ##     g3a <- g3b        
    ## }
    ## rate <- colSums(diff)
    ## try jUpdate version
    diff2 <- array(0, c(3, nGrids - 1))
    for(iGrid in 0:(nGrids - 1 - 2)) {
        ## do not need P(O | lambda), the scaling takes care of that!
        diff2[1, iGrid + 1] <- 1 - sum(transMatRate_tc_H[1, iGrid + 1, s] * (
            alphaHat_t1[, iGrid + 1] *
            betaHat_t1[, iGrid + 1 + 1] *
            eMatGrid_t1[, iGrid + 1 + 1]
        ))
        diff2[2, iGrid + 1] <- 1 - sum(transMatRate_tc_H[1, iGrid + 1, s] * (
            alphaHat_t2[, iGrid + 1] *
            betaHat_t2[, iGrid + 1 + 1] *
            eMatGrid_t2[, iGrid + 1 + 1]
        ))
        diff2[3, iGrid + 1] <- 1 - sum(transMatRate_tc_H[1, iGrid + 1, s] * (
            alphaHat_t3[, iGrid + 1] *
            betaHat_t3[, iGrid + 1 + 1] *
            eMatGrid_t3[, iGrid + 1 + 1]
        ))
    }
    rate2 <- colSums(diff2)
    ## how does jUpate work? probably better?
    ## so want rough cut
    ## will see this more than once
    ## can probably be "too" blocky?
    smoothed_rate <- rcpp_make_smoothed_rate(
        sigma_rate = rate2,
        L_grid = L_grid,
        shuffle_bin_radius = shuffle_bin_radius,
        verbose = FALSE
    )
    ## break_threshold is at most 1, or max or the 50th break_threshold, or the
    break_thresh <- 1
    ## d <- quantile(smoothed_rate, probs = quantile_prob, na.rm = TRUE)
    d <- rcpp_simple_quantile(smoothed_rate, block_gibbs_quantile_prob);
    if (d < break_thresh) {
        break_thresh <- d
    }
    ##
    ## so kind of like with smoothing
    ## start with all available as those above break_threshold
    ## then keep going until done
    ##
    available <- smoothed_rate > break_thresh
    available[smoothed_rate < 0.01] <- FALSE
    available[is.na(smoothed_rate)] <- FALSE
    if (sum(available) == 0) {
        blocked_snps <- rep(0, length(grid))
        return(blocked_snps)
    }
    ## best <- which(available)[order(smoothed_rate[available], decreasing = TRUE)]
    nAvailable <- sum(available)
    best2 <- order(smoothed_rate, decreasing = TRUE)
    best <- best2[1:nAvailable]
    ##
    iBest <- 1
    to_keep <- c()
    for(iBest in 1:length(best)) {
        if (available[best[iBest]]) {
            ## now find start, end for this
            ## this is between grids "snp_best" and grid "snp_best" + 1
            snp_best <- best[iBest]
            ## now, consider where peak ends
            a <- max(snp_best - 1, 1)
            b <- min(snp_best + 1, nGrids)
            if (sum(available[a:b]) == 3) {
                ## no matter what, one further away
                ## return 1-based as well
                snp_left <- R_determine_where_to_stop(smoothed_rate, available, snp_best, break_thresh, nGrids, TRUE)
                snp_right <- R_determine_where_to_stop(smoothed_rate, available, snp_best, break_thresh, nGrids, FALSE)
                available[snp_left:snp_right] <- FALSE
            } else {
                available[a:b] <- FALSE
            }
            ## no matter what, keep
            to_keep <- c(to_keep, snp_best)
        }
    }
    if (min(to_keep) != 0) {
        to_keep <- c(to_keep, 0)
    }
    if (max(to_keep) != (nGrids - 1)) {
        to_keep <- c(to_keep, nGrids - 1)
    }
    blocks_to_consider <- sort(to_keep)
    ##
    if (!TRUE) {
        print(range(grid[which(true_blocked_snps == 0)]))
        print(range(grid[which(true_blocked_snps == 1)]))
        print(range(grid[which(true_blocked_snps == 2)]))
    }
    blocked_snps <- array(-1, nSNPs) ## argh
    blocked_grid <- array(-1, nGrids)
    ## 
    for(i in 0:(length(blocks_to_consider) - 1 - 1)) {
        a <- blocks_to_consider[i + 1] + 1 
        b <- blocks_to_consider[i + 1 + 1] + 1
        blocked_grid[a:b] <- i
    }
    ## 
    for(iSNP in 0:(nSNPs - 1)) {
        blocked_snps[iSNP + 1] <- blocked_grid[grid[iSNP + 1] + 1]    
    }
    ## now - if there are no reads in a block, get rid of
    ## iRead <- -1
    ## for(iSNP in 0:(nSNPs - 1)) {
    ##     b <- blocked_snps[iSNP + 1]
    ##     ## keep going until block is done, then known 
    ## }
    ##
    if (!TRUE) {
        print(range(grid[which(blocked_snps == 0)]))
        print(range(grid[which(blocked_snps == 1)]))
        print(range(grid[which(blocked_snps == 2)]))
    }
    return(
        list(
            blocked_snps = blocked_snps,
            break_thresh = break_thresh,
            smoothed_rate = smoothed_rate
        )
    )
}

for_testing_get_full_package_probabilities <- function(localH, fpp_stuff) {
    ## argh - setup
    transMatRate_tc_H <- fpp_stuff[["transMatRate_tc_H"]]
    alphaMatCurrent_tc <- fpp_stuff[["alphaMatCurrent_tc"]]
    priorCurrent_m <- fpp_stuff[["priorCurrent_m"]]
    eMatRead_t <- fpp_stuff[["eMatRead_t"]]
    s <- fpp_stuff[["s"]]
    sampleReads <- fpp_stuff[["sampleReads"]]
    ##
    K <- dim(alphaMatCurrent_tc)[1]
    nGrids <- dim(alphaMatCurrent_tc)[2] + 1    
    S <- dim(alphaMatCurrent_tc)[3]
    ## 
    package <- lapply(1:3, function(hap) {
        ## now, for first region specifically, introduce the options
        ##
        eMatGrid_t <- array(1, c(K, nGrids))
        rcpp_make_eMatGrid_t(eMatGrid_t = eMatGrid_t, eMatRead_t = eMatRead_t, H = localH, sampleReads = sampleReads, hap = hap, nGrids = nGrids, prev = 0, suppressOutput = 1, prev_section = "text", next_section = "", run_fb_grid_offset = 0, use_all_reads = FALSE, bound = FALSE, maxEmissionMatrixDifference =  1000, rescale = FALSE)
        ## initialize
        out <- initialize_gibbs_forward_backward(H = localH, hap = hap, s = s, sampleReads = sampleReads, priorCurrent_m = priorCurrent_m, alphaMatCurrent_tc = alphaMatCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, eMatRead_t = eMatRead_t, eMatGrid_t = eMatGrid_t)
        ## 
        return(append(out, list(eMatGrid_t = eMatGrid_t)))
    })
    log_p <- -sum(log(package[[1]][["c"]]) + log(package[[2]][["c"]]) + log(package[[3]][["c"]]))
    package <- append(package, list(log_p = log_p))
    return(package)
}





R_consider_block_relabelling <- function(
    iBlock,
    runif_block,
    sum_H,
    s,
    rr,
    rr0,
    ff,
    log_prior_probs,
    logC_before,
    logC_after,
    verbose,
    swap_list,
    eMatGridLocal,    
    betaHatLocal,
    iGrid,
    grid_start_0_based,
    grid_end_0_based,
    read_start_0_based,
    read_end_0_based,
    wif0,
    log_cStore,
    alphaStore,
    read_is_uninformative,
    block_approach,
    do_checks,
    all_packages,
    block_results,
    ever_changed,
    transMatRate_tc_H,
    alphaMatCurrent_tc,
    priorCurrent_m,
    fpp_stuff,
    H,
    proposed_H,
    nReads,
    eMatRead_t,
    alphaHat_t1,
    betaHat_t1,
    c1,
    eMatGrid_t1,
    alphaHat_t2,
    betaHat_t2,
    c2,
    eMatGrid_t2,
    alphaHat_t3,
    betaHat_t3,
    c3,
    eMatGrid_t3
) {
    ## 
    betaHatLocal[, 1] <- betaHat_t1[, iGrid]
    betaHatLocal[, 2] <- betaHat_t2[, iGrid]
    betaHatLocal[, 3] <- betaHat_t3[, iGrid]
    ##
    ## get this as well - inefficient
    ## 
    choice_log_probs_P <- array(0, 6)
    choice_log_probs_Pm <- array(0, c(6, 3))
    for(ir in 1:6) {
        for(i in 1:3) {
            ## get the c influence from within the grid
            logC_inside <- 0
            for(iGrid2 in grid_start_0_based:grid_end_0_based) {
                logC_inside <- logC_inside + sum(log_cStore[iGrid2 + 1, i, ir])
            }
            ## OK, so after includes original one
            ## inside one includes current c
            choice_log_probs_Pm[ir, i] <- 
                log(sum(alphaStore[, i, ir] * betaHatLocal[, i])) +
                - logC_before[i] +
                - logC_inside + 
                - logC_after[i]
            choice_log_probs_P[ir] <- choice_log_probs_P[ir] + choice_log_probs_Pm[ir, i]
        }
    }
    ## do check
    if (do_checks & (block_approach == 1)) {
        log_package_probs <- sapply(all_packages, function(x) x[["log_p"]])
        expect_equal(as.numeric(choice_log_probs_P), log_package_probs)
    }
    ## now add on real label probabilities
    if (block_approach == 1) {
        choice_log_probs_H <- calculate_block_read_label_probabilities(
            read_start_0_based = read_start_0_based,
            read_end_0_based = read_end_0_based,
            H = H,
            log_prior_probs = log_prior_probs,
            rr = rr
        )
    } else if (block_approach == 2) {
        choice_log_probs_H <- calculate_block_read_label_probabilities_using_read_informativeness(
            read_start_0_based = read_start_0_based,
            read_end_0_based = read_end_0_based,
            H = H,
            log_prior_probs = log_prior_probs,
            rr = rr,
            read_is_uninformative = read_is_uninformative
        )
    } else if (block_approach == 4) {
        choice_log_probs_H <- calculate_block_read_label_probabilities_using_proposed_H(
            read_start_0_based = read_start_0_based,
            read_end_0_based = read_end_0_based,
            proposed_H = proposed_H,
            log_prior_probs = log_prior_probs,
            rr = rr
        )            
    }
    ## now choose!
    choice_log_probs <- choice_log_probs_P + choice_log_probs_H
    ## now do sampling
    choice_log_probs <- choice_log_probs - max(choice_log_probs)
    ## cbind(choice_log_probs, choice_log_probs_P, choice_log_probs_H)
    choice_log_probs[choice_log_probs < (-100)] <- (-100)
    choice_probs <- exp(choice_log_probs)
    choice_probs[is.na(choice_probs)] <- 0
    if (ff == 0) {
        ##// only allow the two options
        choice_probs[2] <- 0
        choice_probs[4] <- 0
        choice_probs[5] <- 0
        choice_probs[6] <- 0
    }
    choice_probs <- choice_probs / sum(choice_probs)
    ##
    ## ir_chosen <- sample(1:6, size = 1, prob = choice_probs)
    chance <- runif_block[iBlock + 1]
    cumsum_choice_probs <- cumsum(choice_probs)
    ##
    ir_chosen <- 0
    for(i in 6:1) {
        if (chance < cumsum_choice_probs[i]) {
            ir_chosen <- i
        }
    }
    ## 
    if (verbose) {
        print_message(paste0("In block ", iBlock, ", see the following probabilities"))
        print_message(paste0(round(choice_probs, 4), collapse = ", "))
        print_message(paste0("Have selected block relabelling:", ir_chosen))
    }
    ## 
    ## store some probabilities about chosen change
    ##
    ibr <- 2 * iBlock + 1
    block_results[ibr, "iBlock"] <- iBlock
    block_results[ibr, "total"] <- 0
    block_results[ibr, paste0("p", 1:6)] <- choice_probs
    block_results[ibr, "ir_chosen"] <- ir_chosen
    block_results[ibr, c("p_O1_given_H1_L", "p_O2_given_H2_L", "p_O3_given_H3_L")] <- choice_log_probs_Pm[ir_chosen, ]
    block_results[ibr, "p_O_given_H_L"] <- sum(choice_log_probs_Pm[ir_chosen, ])
    ## this should be recorded for the ENTIRE set of reads NOT specific
    ## although the chosen set of 
    if (do_checks) {
        if (verbose) {
            print_message("Check H for block_results")
        }
        expect_equal(sum(H == 1), sum_H[1])
        expect_equal(sum(H == 2), sum_H[2])
        expect_equal(sum(H == 3), sum_H[3])
    }
    ## not entirely true for approach 2
    x <- 0
    for(i in 1:3) {
        if (sum_H[i] > 0) {
            x <- x + log_prior_probs[i] * sum_H[i]
        }
    }
    for(iRead in read_start_0_based:read_end_0_based) {
        h1 <- H[iRead + 1]        
        if (block_approach == 1 | block_approach == 2) {
            h2 <- h1
        } else if (block_approach == 4) {
            h2 <- proposed_H[ir_chosen, iRead + 1] 
        }
        x <- x - log_prior_probs[h1]
        x <- x + log_prior_probs[rr[ir_chosen, h2]]
    }
    block_results[ibr, "p_H_given_L"] <- x
    if (do_checks) {
        Htemp <- H
        for(iRead in read_start_0_based:read_end_0_based) {
            if (block_approach == 1 | block_approach == 2) {
                h <- H[iRead + 1]                
            } else if (block_approach == 4) {
                h <- proposed_H[ir_chosen, iRead + 1] 
            }
            Htemp[iRead + 1] <- rr[ir_chosen, h]
        }
        expect_equal(x, sum(log_prior_probs[Htemp]))
    }
    block_results[ibr, "p_H_given_O_L_up_to_C"] <- block_results[ibr, "p_O_given_H_L"] + block_results[ibr, "p_H_given_L"]
    ##
    if (!((ever_changed == 1)| (ir_chosen != 1))) {
        if (verbose) {
            print_message("No change warranted")
        }
    } else {
        if (ir_chosen != 1) {
            if (verbose) {
                print_message("Apply block relabelling")
                print_message(paste0("i.e. the following relabelling occurs:"))
                for(i in 1:3) {
                    print_message(paste0("previous label ", i, " -> ", rr[ir_chosen, i]))
                }
            }
        } else {
            if (verbose) {
                print_message("Reset alpha post-block checking")
            }
        }
        ever_changed[] <- 1
        ## apply the labelling!
        ## now, go back, and apply the switch
        iRead <- read_start_0_based + 1
        wif_read <- wif0[iRead]
        for(iGrid2 in (grid_start_0_based:grid_end_0_based + 1)) {
            if (verbose) {
                print(paste0("changing iGrid2 (1-based) = ", iGrid2))
            }
            ## rebuild eMatGridLocal using proposed H, if applicable!
            if (block_approach == 1 | block_approach == 2) {
                eMatGridLocal[, rr[ir_chosen, 1]] <- eMatGrid_t1[, iGrid2]
                eMatGridLocal[, rr[ir_chosen, 2]] <- eMatGrid_t2[, iGrid2]
                eMatGridLocal[, rr[ir_chosen, 3]] <- eMatGrid_t3[, iGrid2]
            } else if (block_approach == 4) {
                ##
                ## rebuild eMatGrid under
                ##
                eMatGridLocal[] <- 1
                while((iRead <= nReads) & (wif_read < (iGrid2 - 1))) {
                    iRead <- iRead + 1
                    if (iRead <= nReads) {
                        wif_read <- wif0[iRead]
                    }
                }
                while (iRead <= nReads & (wif_read == (iGrid2 - 1))) {
                    Hl <- proposed_H[ir_chosen, iRead] 
                    h <- rr[ir_chosen, Hl]
                    eMatGridLocal[, h] <- eMatGridLocal[, h] * eMatRead_t[, iRead]
                    ##
                    iRead <- iRead + 1
                    if (iRead <= nReads) {
                        wif_read <- wif0[iRead]
                    }
                }
            }
            ## 
            ## reset eMatGrid_t s
            ## 
            eMatGrid_t1[, iGrid2] <- eMatGridLocal[, 1]
            eMatGrid_t2[, iGrid2] <- eMatGridLocal[, 2]
            eMatGrid_t3[, iGrid2] <- eMatGridLocal[, 3]
            ##
            ## do the rest!
            ##
            if (iGrid2 == 1) {
                ## first grid!
                alphaHat_t1[, iGrid2] <- priorCurrent_m[, s] * eMatGrid_t1[, iGrid2]
                alphaHat_t2[, iGrid2] <- priorCurrent_m[, s] * eMatGrid_t2[, iGrid2]
                alphaHat_t3[, iGrid2] <- priorCurrent_m[, s] * eMatGrid_t3[, iGrid2]
            } else {
                ## normal!
                ## now go forward in traditional way
                alphaHat_t1[, iGrid2] <- eMatGrid_t1[, iGrid2] * (
                    transMatRate_tc_H[1, iGrid2 - 1, s] * alphaHat_t1[, iGrid2 - 1] + 
                    transMatRate_tc_H[2, iGrid2 - 1, s] * alphaMatCurrent_tc[, iGrid2 - 1, s]
                )
                alphaHat_t2[, iGrid2] <- eMatGrid_t2[, iGrid2] * (
                    transMatRate_tc_H[1, iGrid2 - 1, s] * alphaHat_t2[, iGrid2 - 1] + 
                    transMatRate_tc_H[2, iGrid2 - 1, s] * alphaMatCurrent_tc[, iGrid2 - 1, s]
                )
                alphaHat_t3[, iGrid2] <- eMatGrid_t3[, iGrid2] * (
                    transMatRate_tc_H[1, iGrid2 - 1, s] * alphaHat_t3[, iGrid2 - 1] + 
                    transMatRate_tc_H[2, iGrid2 - 1, s] * alphaMatCurrent_tc[, iGrid2 - 1, s]
                )
            }
            ## re-normalize
            c1[iGrid2] <- 1 / sum(alphaHat_t1[, iGrid2])
            alphaHat_t1[, iGrid2] <- c1[iGrid2] * alphaHat_t1[, iGrid2]
            ##
            c2[iGrid2] <- 1 / sum(alphaHat_t2[, iGrid2])
            alphaHat_t2[, iGrid2] <- c2[iGrid2] * alphaHat_t2[, iGrid2]
            ##
            c3[iGrid2] <- 1 / sum(alphaHat_t3[, iGrid2])
            alphaHat_t3[, iGrid2] <- c3[iGrid2] * alphaHat_t3[, iGrid2]
        }
        ## finally, re-do read labels
        for(iRead in read_start_0_based:read_end_0_based) {
            lost <- H[iRead + 1]
            if (block_approach == 1 | block_approach == 2) {
                gained <- rr[ir_chosen, lost]
            } else if (block_approach == 4) {
                gained <- rr[ir_chosen, proposed_H[ir_chosen, iRead + 1] ]
            }
            H[iRead + 1] <- gained
            ## one loses
            sum_H[gained] <- sum_H[gained] + 1
            ## one gains
            sum_H[lost] <- sum_H[lost] - 1
        }
        if (do_checks) {
            if (verbose) {
                print_message("Check H after reassignment")
            }
            expect_equal(sum(H == 1), sum_H[1])
            expect_equal(sum(H == 2), sum_H[2])
            expect_equal(sum(H == 3), sum_H[3])
        }
        ## 
    }
    return(
        list(
            eMatGrid_t1 = eMatGrid_t1,
            eMatGrid_t2 = eMatGrid_t2,
            eMatGrid_t3 = eMatGrid_t3,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            alphaHat_t1 = alphaHat_t1,
            alphaHat_t2 = alphaHat_t2,
            alphaHat_t3 = alphaHat_t3,            
            betaHat_t1 = betaHat_t1,
            betaHat_t2 = betaHat_t2,
            betaHat_t3 = betaHat_t3,
            H = H,
            ever_changed = ever_changed,
            block_results = block_results,
            sum_H = sum_H
        )
    )
}


R_consider_total_relabelling <- function(
    iBlock,
    rr,
    rr0,
    ff,
    log_prior_probs,
    logC_before,
    logC_after,
    verbose,
    swap_list,
    block_results,
    runif_total,
    sum_H,
    H,
    alphaHat_t1,
    betaHat_t1,
    c1,
    eMatGrid_t1,
    alphaHat_t2,
    betaHat_t2,
    c2,
    eMatGrid_t2,
    alphaHat_t3,
    betaHat_t3,
    c3,
    eMatGrid_t3,
    do_checks = FALSE
) {
    ##
    ## consider entire relabelling here
    ##
    choice_log_probs_Pm <- array(0, c(6, 3))
    choice_log_probs_P <- array(0, 6)
    for(ir in 1:6) {
        for(i in 1:3) {
            n <- sum_H[rr[ir, i]]
            if (n > 0) {
                choice_log_probs_Pm[ir, i] <- log_prior_probs[i] * n
            }
            choice_log_probs_P[ir] <- choice_log_probs_P[ir] + choice_log_probs_Pm[ir, i]
        }
    }
    choice_log_probs_P <- choice_log_probs_P - max(choice_log_probs_P)
    choice_log_probs_P[(choice_log_probs_P) < (-100)] <- -100
    choice_probs <- exp(choice_log_probs_P)
    if (ff == 0) {
        ## shouldn't happen, but anyway, force only allow 1 or 3
        choice_probs[2] <- 0
        choice_probs[4] <- 0
        choice_probs[5] <- 0
        choice_probs[6] <- 0        
    }
    choice_probs <- choice_probs / sum(choice_probs)
    ##
    chance <- runif_total[iBlock + 1]
    cumsum_choice_probs <- cumsum(choice_probs)
    ##
    ir_chosen <- 0
    for(i in 6:1) {
        if (chance < cumsum_choice_probs[i]) {
            ir_chosen <- i
        }
    }
    ##ir_chosen <- sample(1:6, 1, prob = choice_probs)
    ##
    if (verbose) {
        print_message(paste0("After block ", iBlock, ", on total relabelling, see the following probabilities"))
        print_message(paste0(round(choice_probs, 4), collapse = ", "))
        print_message(paste0("Have selected total relabelling:", ir_chosen, " using chance:", round(chance, 3)))
    }
    ## 
    ## store some probabilities about chosen change
    ##
    ibr <- 2 * iBlock + 1 + 1
    block_results[ibr, "iBlock"] <- iBlock
    block_results[ibr, "total"] <- 1
    block_results[ibr, paste0("p", 1:6)] <- choice_probs
    block_results[ibr, "ir_chosen"] <- ir_chosen
    ## these, but re-ordered
    b <- block_results[ibr - 1, c("p_O1_given_H1_L", "p_O2_given_H2_L", "p_O3_given_H3_L")]
    a <- array(0, 3)
    for(i in 1:3) {
        h <- rr[ir_chosen, i]
        a[h] <- b[i]
    }
    block_results[ibr, c("p_O1_given_H1_L", "p_O2_given_H2_L", "p_O3_given_H3_L")] <- b
    block_results[ibr, c("p_O_given_H_L")] <- sum(b)
    ## this is new set after change
    ## make this efficient!
    x <- 0
    for(i in 1:3) {
        h <- rr[ir_chosen, i]
        x <- x + log_prior_probs[h] * sum_H[i]
    }
    if (do_checks) {
        ## inefficient vs efficient way!
        Htemp <- rr[ir_chosen, H]
        expect_equal(table(H), sum_H)
        expect_equal(x, sum(log_prior_probs[Htemp]))
    }        
    block_results[ibr, "p_H_given_L"] <- x
    block_results[ibr, "p_H_given_O_L_up_to_C"] <- block_results[ibr, "p_O_given_H_L"] + block_results[ibr, "p_H_given_L"]
    ##
    if (ir_chosen == 1) {
        if (verbose) {
            print_message("No total relabelling to apply")
        }
    } else {
        if (verbose) {
            print_message("Apply total relabelling")
            print_message(paste0("i.e. the following relabelling occurs:"))
            for(i in 1:3) {
                print_message(paste0("previous label ", i, " -> ", rr[ir_chosen, i]))
            }
        }
        ## H, alpha, beta, eMatGrid, c
        H <- rr[ir_chosen, H]
        to_swap <- swap_list[[ir_chosen]]
        for(what in c("c", "alphaHat_t", "betaHat_t", "eMatGrid_t")) {
            for(i_swap in 1:length(to_swap)) {
                to_swap2 <- to_swap[[i_swap]]
                ## have c++ swap function for this too
                eval(parse(text = paste0("x=", what, to_swap2[1])))
                eval(parse(text = paste0(what, to_swap2[1], "=", what, to_swap2[2])))
                eval(parse(text = paste0(what, to_swap2[2], "=x")))
            }
        }
        a <- logC_before
        b <- logC_after
        c <- sum_H
        for(i_swap in 1:length(to_swap)) {
            ## can do more easily
            for(i in 1:3) {
                h <- rr[ir_chosen, i]
                logC_before[h] <- a[i]
                logC_after[h] <- b[i]
                sum_H[h] <- c[i]
            }
        }
    }
    return(
        list(
            eMatGrid_t1 = eMatGrid_t1,
            eMatGrid_t2 = eMatGrid_t2,
            eMatGrid_t3 = eMatGrid_t3,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            alphaHat_t1 = alphaHat_t1,
            alphaHat_t2 = alphaHat_t2,
            alphaHat_t3 = alphaHat_t3,            
            betaHat_t1 = betaHat_t1,
            betaHat_t2 = betaHat_t2,
            betaHat_t3 = betaHat_t3,
            H = H,
            logC_after = logC_after,
            logC_before = logC_before,
            ir_chosen = ir_chosen,
            block_results = block_results,
            sum_H = sum_H
        )
    )
}



calculate_block_read_label_probabilities <- function(
    read_start_0_based,
    read_end_0_based,
    H,
    log_prior_probs,
    rr
) {
    choice_log_probs_H <- array(0, 6)
    for(iRead in read_start_0_based:read_end_0_based) {
        h <- H[iRead + 1]
        for(ir in 1:6) {
            choice_log_probs_H[ir] <- choice_log_probs_H[ir] + log_prior_probs[rr[ir, h]]
        }
    }
    return(choice_log_probs_H)
}

calculate_block_read_label_probabilities_using_read_informativeness <- function(
    read_start_0_based,
    read_end_0_based,
    H,
    log_prior_probs,
    rr,
    read_is_uninformative
) {
    choice_log_probs_H <- array(0, 6)
    for(iRead in read_start_0_based:read_end_0_based) {
        if (!read_is_uninformative[iRead + 1]) {
            h <- H[iRead + 1]
            for(ir in 1:6) {
                choice_log_probs_H[ir] <- choice_log_probs_H[ir] + log_prior_probs[rr[ir, h]]
            }
        }
    }
    return(choice_log_probs_H)
}


calculate_block_read_label_probabilities_using_proposed_H <- function(
    read_start_0_based,
    read_end_0_based,
    proposed_H,
    log_prior_probs,
    rr
) {
    choice_log_probs_H <- array(0, 6)
    for(iRead in read_start_0_based:read_end_0_based) {
        for(ir in 1:6) {
            h <- rr[ir, proposed_H[ir, iRead + 1]]
            choice_log_probs_H[ir] <- choice_log_probs_H[ir] + log_prior_probs[h]
        }
    }
    return(choice_log_probs_H)
}



R_reset_local_variables <- function(
    iGrid,
    verbose,
    alphaHatLocal,
    alphaStore,
    alphaHat_t1,
    alphaHat_t2,
    alphaHat_t3,
    c1,
    c2,
    c3,
    log_cStore
) {
    if (verbose) {
        print_message("Reset alphaHatLocal, cStore")
    }
    ## reset alphaStore, cStore
    ## they are all using the current data
    ## check against new truth
    alphaHatLocal[, 1] <- alphaHat_t1[, iGrid]
    alphaHatLocal[, 2] <- alphaHat_t2[, iGrid]
    alphaHatLocal[, 3] <- alphaHat_t3[, iGrid]
    cLocal <- c(c1[iGrid], c2[iGrid], c3[iGrid])
    ## alphaStore[] <- 0
    for(ir in 1:6) {
        for(i in 1:3) {
            ## h <- rr[ir, i]
            h <- i
            alphaStore[, h, ir] <- alphaHatLocal[, i]
            log_cStore[iGrid, h, ir] <- log(cLocal[i])
        }
    }
    return(
        list(
            alphaStore = alphaStore,
            log_cStore = log_cStore
        )
    )
}


add_grey_background <- function(L_grid) {
    from <- floor(min(L_grid) / 1e6)
    to <- ceiling(max(L_grid) / 1e6)
    for(f in from:to) {
        abline(v = 1e6 * f, col = "grey")
    }
}




plot_attempt_to_reblock_snps <- function(
    outname,
    break_thresh,
    considers,
    grid_distances,
    L_grid,
    gibbs_block_output_list,
    smoothed_rate,
    L,
    block_results,
    uncertain_truth_labels,
    truth_labels,
    have_truth_haplotypes,
    sampleReads
) {
    ##
    xlim <- range(L_grid)
    consider_snp_start_0_based <- considers[["consider_snp_start_0_based"]]
    consider_snp_end_0_based <- considers[["consider_snp_end_0_based"]]
    consider_reads_start_0_based <- considers[["consider_reads_start_0_based"]]
    consider_reads_end_0_based <- considers[["consider_reads_end_0_based"]]
    ##
    before_gamma1_t <- gibbs_block_output_list[["before_gamma1_t"]]
    before_gamma2_t <- gibbs_block_output_list[["before_gamma2_t"]]
    after_gamma1_t <- gibbs_block_output_list[["after_gamma1_t"]]
    after_gamma2_t <- gibbs_block_output_list[["after_gamma2_t"]]
    before_read_labels <- gibbs_block_output_list[["before_read_labels"]]        
    after_read_labels <- gibbs_block_output_list[["after_read_labels"]]
    ## 
    ## width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 36), 200)
    png(outname, height = 10, width = 20, res = 100, units = "in")
    par(mfrow = c(7, 1))
    par(oma = c(0, 0, 5, 0))    
    grid_distances <- diff(L_grid)
    x <- L_grid[-1] - grid_distances    
    xleft <- L_grid[-length(L_grid)]
    xright <- L_grid[-1]
    ## for fbdstore
    midpoints <- L_grid[-1] - grid_distances / 2    
    xleft2 <- c(L_grid[1], midpoints)
    xright2 <- c(midpoints, L_grid[length(L_grid)])
    ##
    ## 1) find peaks to check
    ##
    ylim <- c(0, max(break_thresh, max(smoothed_rate, na.rm = TRUE)))
    ylim <- c(0, max(break_thresh, quantile(smoothed_rate, probs = c(0.99))))
    par(mar = c(0, 0, 3, 0))            
    plot(x = 0, y = 0, xlab = "Physical position", ylab = "Rate", main = "Location of shuffles to check", ylim = ylim, xlim = xlim)
    add_grey_background(L_grid)
    y <- smoothed_rate
    y[y > ylim[2]] <- ylim[2]
    lines(x = x, y = y, lwd = 2)
    for(iBlock in 0:(length(consider_snp_start_0_based) - 1)) {
        l <- L[consider_snp_start_0_based[iBlock + 1] + 1]
        r <- L[consider_snp_end_0_based[iBlock + 1] + 1]
        abline(v = l, col = "red")
        abline(v = r, col = "red")
        text(x = (l + r) / 2, y = ylim[2] - diff(ylim) * 0.1, labels = iBlock)
        text(x = (l + r) / 2, y = ylim[2] - diff(ylim) * 0.25, labels = round(block_results[iBlock + 1, "p1"], 2))
        text(x = (l + r) / 2, y = ylim[2] - diff(ylim) * 0.4, labels = round(block_results[iBlock + 1, "p3"], 2))
        ## if have (manual) read labels, add them
        rs <- consider_reads_start_0_based[iBlock + 1] + 1
        re <- consider_reads_end_0_based[iBlock + 1] + 1
        ## get counts
        ## mrll <- manual_read_labels[rs:re]
        ## labels <- sapply(c('v0', 'v3', 'v2', 'v23', 'v1', 'v13', 'v12', 'v123', "other"), function(v) {
        ##     sum(mrll == v, na.rm = TRUE)
        ## })
        ##labels <- sapply(0:7, function(v) {
        ##    sum(mrll == v, na.rm = TRUE)
        ##})
        ##names(labels) <- c("other", "v1", "v2", "v3", "v12", "v13", "v23", "v123")
        ##text(x = (l + r) / 2, y = ylim[1] + diff(ylim) * 0.1 + diff(ylim) * 0.6 * seq(0, 1, length.out = length(labels)), labels = paste0("class", names(labels), "=", labels))
        ## get col
        changed <- sum(block_results[which(block_results[, "iBlock"] == iBlock), "ir_chosen"] != 1) > 0
        if (changed) {
            col <- "green"
        } else {
            col <- "red"
        }
        ## rectangle underneath as well
        rect(
            xleft = l,
            xright = r,
            ybottom = -1, ytop = 0,
            col = col
        ) ##whichIsBest 0 = no switch
    }
    abline(h = break_thresh, col = "black", lwd = 2)    
    for(i_before in 1:2) {
        if (i_before == 1) {
            read_labels <- before_read_labels
            gamma1_t <- before_gamma1_t
            gamma2_t <- before_gamma2_t            
        } else {
            read_labels <- after_read_labels
            gamma1_t <- after_gamma1_t
            gamma2_t <- after_gamma2_t            
        }
        ## 
        ## 2) add in reads here
        ##
        par(mar = c(0, 0, 3, 0))            
        ylim <- c(0, 1)
        plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, cex = 1.5)
        if (have_truth_haplotypes) {
            truth <- truth_labels
            truth[uncertain_truth_labels] <- 0
        }
        y <- 0.1 + (read_labels - 1) / 2 + runif(length(read_labels)) / 4
        for(iRead in 1:length(sampleReads)) {
            u <- range(sampleReads[[iRead]][[4]])
            if (have_truth_haplotypes) {
                col <- c("black", "blue", "red")[truth[iRead] + 1]
            } else {
                col <- "black"
            }
            lwd <- 1
            segments(x0 = L[u[1] + 1], x1 = L[u[2] + 1], y0 = y[iRead], y1 = y[iRead], col = col, lwd = lwd)
        }
        ## 
        ##
        ## 3) plot gammas
        ##
        ## argh, make fake
        ##
        ## for now, plot gammas before
        ## 
        for(i_which in 1:2) {
            par(mar = c(0, 0, 3, 0))        
            scale_dosage <- 0
            colStore <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
            nCols <- length(colStore)
            if (i_which == 1) { gammaK_t <- gamma1_t;   main <- "Hap 1"}
            if (i_which == 2) { gammaK_t <- gamma2_t;   main <- "Hap 2" }
            ##
            K <- nrow(gammaK_t)
            ylim <- c(0, 1 + scale_dosage + scale_dosage)
            nGrids <- ncol(gammaK_t)
            backwards <- nGrids:1
            ##
            plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, cex = 1.5, main = main)
            x <- L_grid ## c(L_grid[1], L_grid) ## , L_grid[length(L_grid):1])
            xleft <- c(x[1] - (x[2] - x[1]) / 2, x[-length(x)])
            xright <- c(x[-1], x[length(x)] + (x[length(x)] - x[(length(x) - 1)]) / 2)
            m <- array(0, c(nGrids, K + 1))
            ## is this slow...
            for(i in 1:K) {
                m[, i + 1] <- m[, i] + gammaK_t[i, ]
            }
            ##
            for(j in 1:2) {
                for(i in K:1) {
                    ## can I o
                    ybottom <- m[, i]
                    ytop <- m[, i + 1]
                    if (max(ytop - ybottom) > 0.01) {
                        if (j == 1) {
                            rect(
                                xleft = xleft,
                                xright = xright,
                                ybottom = ybottom,
                                ytop = ytop,
                                border = NA,
                                col = colStore[(i %% nCols) + 1]
                            )
                        } else {
                            add_numbers(ytop, ybottom, x, i)
                        }
                    }
                }
            }
            for(iBlock in 0:(length(consider_snp_start_0_based) - 1)) {
                l <- L_grid[consider_snp_start_0_based[iBlock + 1] + 1]
                r <- L_grid[consider_snp_end_0_based[iBlock + 1] + 1]
                abline(v = l, col = "red")
                abline(v = r, col = "red")
                text(x = (l + r) / 2, y = ylim[2] - diff(ylim) * 0.1, labels = iBlock)
            }
        }
    }
    dev.off()
}




calculate_H_class <- function(
    eMatRead_t,
    alphaHat_t1,
    alphaHat_t2,
    alphaHat_t3,
    betaHat_t1,
    betaHat_t2,
    betaHat_t3,
    ff,
    wif0,
    H,
    class_sum_cutoff = 0.06
) {
    ## 
    prior_probs <- c(0.5, (1 - ff) / 2, ff / 2)
    p <- prior_probs    
    rlc <- matrix(0, nrow = 7, ncol = 3)
    rlc[1, ] <- c(1, 0, 0)
    rlc[2, ] <- c(0, 1, 0)
    rlc[3, ] <- c(0, 0, 1)
    rlc[4, ] <- c(p[1] / (p[1] + p[2]), p[2] / (p[1] + p[2]), 0)
    rlc[5, ] <- c(p[1] / (p[1] + p[3]), 0, p[3] / (p[1] + p[3]))
    rlc[6, ] <- c(0, p[2] / (p[2] + p[3]), p[3] / (p[2] + p[3]))
    rlc[7, ] <- p
    H_class <- integer(ncol(eMatRead_t))
    ## probs <- array(NA, c(nReads, 3))
    ## have normal probs and eMatGrid_t1, etc
    ## remove from current
    ## have eMatGrid_t1
    for(iRead in 1:ncol(eMatRead_t)) {
        ##
        iGrid <- wif0[iRead] + 1
        ##
        ## in an inefficient manner, calculate!
        ## remove from current
        alphaHat_m <- rbind(alphaHat_t1[, iGrid], alphaHat_t2[, iGrid], alphaHat_t3[, iGrid])
        betaHat_m <- rbind(betaHat_t1[, iGrid], betaHat_t2[, iGrid], betaHat_t3[, iGrid])
        ## 
        ## remove from current
        ##
        pC <- rowSums(alphaHat_m * betaHat_m)
        h_rC <- H[iRead] ## h_r current
        h_rA1 <- setdiff(1:3, H[iRead])[1] ## h_r alternate 1
        h_rA2 <- setdiff(1:3, H[iRead])[2] ## h_r alternate 1
        ## 
        pA1 <- pC ## read label becomes h_rA1
        pA2 <- pC ## read label becomes h_rA2
        ## now - blank out the two that switch
        ## default behaviour
        pA1[c(h_rC, h_rA1)] <- 0
        pA2[c(h_rC, h_rA2)] <- 0
        ## A1 loser
        pA1[h_rC] <- pA1[h_rC] + sum(alphaHat_m[h_rC, ] * betaHat_m[h_rC, ] / eMatRead_t[, iRead])
        ## A2 lower
        pA2[h_rC] <- pA1[h_rC]        
        ## A1 gainer
        pA1[h_rA1] <- pA1[h_rA1] + sum(alphaHat_m[h_rA1, ] * betaHat_m[h_rA1, ] * eMatRead_t[, iRead])
        ## A2 gainer
        pA2[h_rA2] <- pA2[h_rA2] + sum(alphaHat_m[h_rA2, ] * betaHat_m[h_rA2, ] * eMatRead_t[, iRead])
        prod_pC <- prod(pC) * prior_probs[h_rC]
        prod_pA1 <- prod(pA1) * prior_probs[h_rA1]
        prod_pA2 <- prod(pA2) * prior_probs[h_rA2]
        denom <- (prod_pC + prod_pA1 + prod_pA2)        
        norm_pC <- prod_pC / denom
        norm_pA1 <- prod_pA1 / denom
        norm_pA2 <- prod_pA2 / denom
        choice_probs <- c(norm_pC, norm_pA1, norm_pA2)
        ##
        ##
        ## now compare as usual
        ##
        ## determine how close it is to the three probs
        x <- array(NA, 3)
        x[h_rC] <- norm_pC
        x[h_rA1] <- norm_pA1
        x[h_rA2] <- norm_pA2
        ## probs[iRead, ] <- x
        ## check if those three probs are the default option probs
        ## within 0.01 of any option
        ## (what if they overlap? - future problem?)
        ## get label as 1 through 6 (or 7), or not!
        y <- array(NA, 7)
        for(i in 1:7) {
            ## check closest, and within range
            y[i] <- sum(abs(rlc[i, ] - x))
        }
        if ((min(y) < (class_sum_cutoff))) {
            H_class[iRead] <- as.integer(which.min(y))
        } else {
            H_class[iRead] <- as.integer(0)
        }
    }
    ## x <- sort(table(apply(round(probs[H_class == 0, ], 3), 1, paste0, collapse = ", ")))
    ## cbind(x, x)
    return(H_class)
}


R_gibbs_block_forward_one <- function(
    approach2_iRead,
    iGrid,
    s,
    alphaStore,
    log_cStore,
    rr,
    rr0,
    eMatGridLocal,
    eMatGridLocalc,    
    transMatRate_tc_H,
    alphaMatCurrent_tc,
    priorCurrent_m,
    read_is_uninformative,
    block_approach,
    wif0,
    eMatRead_t,
    nReads,
    H,
    proposed_H,
    H_class,
    rlc,
    rlcM,
    runif_proposed
) {
    ##
    ##
    if (block_approach == 2 | block_approach == 4) {
        if ((approach2_iRead[1] <= nReads)) {
            wif_read <- wif0[approach2_iRead[1]]
        }
        while ((approach2_iRead[1] <= nReads) & wif_read < (iGrid - 1)) {
            approach2_iRead[1] <- approach2_iRead[1] + 1
            if (approach2_iRead[1] <= nReads) {
                wif_read <- wif0[approach2_iRead[1]]
            }
        }
        if (block_approach == 2) {
            ## approach2_iRead[1] is 1-based in R
            while (approach2_iRead[1] <= (nReads - 1) & wif_read == (iGrid - 1)) {
                el <- eMatRead_t[, approach2_iRead[1]]
                Hl <- H[approach2_iRead[1]]
                if (read_is_uninformative[approach2_iRead[1]]) {
                    eMatGridLocal[, Hl] <- eMatGridLocal[, Hl] / el
                }
                approach2_iRead[1] <- approach2_iRead[1] + 1
                if (approach2_iRead[1] <= (nReads)) {
                    wif_read <- wif0[approach2_iRead[1]]
                }
            }
        } else if (block_approach == 4) {
            ## per-read, sample read
            ## add to appropriate
            ## then add to 
            ## opposite, start with BLANK. add back in based on sampling for the H
            eMatGridLocalc[] <- 1
            if ((approach2_iRead[1] <= nReads)) {
                wif_read <- wif0[approach2_iRead[1]]
            }
            while ((approach2_iRead[1] <= nReads) & wif_read == (iGrid - 1)) {
                ## old stuff
                Hori <- H[approach2_iRead[1]] ## current label
                Hc <- H_class[approach2_iRead[1]] ## current class
                eMatRead_t_col <- eMatRead_t[, approach2_iRead[1]]
                for(ir in 1:6) {
                    ## if it's Hc = 0, i.e. no assignment, keep as normal
                    if (Hc == 0) {
                        Hl <- as.integer(Hori) ## as.integer(rr[ir, Hori])
                    } else {
                        ## so, have class Hc
                        ## Hl <- as.integer(sample(c(1L, 2L, 3L), size = 1, prob = rlcM[Hc, ir, ]))
                        chance_prob <- runif_proposed[ir, approach2_iRead[1]]
                        if (chance_prob < rlcM[1, ir, Hc]) {
                            Hl <- 1
                        } else if (rlcM[1, ir, Hc] < chance_prob & chance_prob < (rlcM[1, ir, Hc] + rlcM[2, ir, Hc])) {
                            Hl <- 2
                        } else {
                            Hl <- 3
                        }
                    }
                    ## -->HERE<-- record label in "current" notation
                    ## pass it through rr later to get in "real" notation
                    proposed_H[ir, approach2_iRead[1]] <- Hl ## store as 1-based
                    ## -->HERE<-- use "past" label, as rr gets applied below
                    eMatGridLocalc[, Hl, ir] <- eMatGridLocalc[, Hl, ir] * eMatRead_t_col
                }
                approach2_iRead[1] <- approach2_iRead[1] + 1
                if (approach2_iRead[1] <= (nReads)) {
                    wif_read <- wif0[approach2_iRead[1]]
                }
            }
            ##
            ## continue!
            ##
        }
    }
    ##
    if (iGrid == 1) {
        for(ir in 1:6) {
            if (block_approach == 4) {
                eMatGridLocal <- eMatGridLocalc[, , ir]
            }
            for(i in 1:3) {
                h <- rr[ir, i]
                alphaStore[, h, ir] <- priorCurrent_m[, s] * eMatGridLocal[, i]
                d <- 1 / sum(alphaStore[, h, ir])
                log_cStore[1, h, ir] <- log(d)
                alphaStore[, h, ir] <- d * alphaStore[, h, ir]
            }
        }
    } else {
        for(ir in 1:6) {
            if (block_approach == 4) {
                eMatGridLocal <- eMatGridLocalc[, , ir]
            }
            for(i in 1:3) {
                h <- rr[ir, i]
                ## here, choose which eMatGrid_t to use
                ## now, can proceed!
                alphaStore[, h, ir] <- eMatGridLocal[, i] * (
                    transMatRate_tc_H[1, iGrid - 1, s] * alphaStore[, h, ir] +
                    transMatRate_tc_H[2, iGrid - 1, s] * alphaMatCurrent_tc[, iGrid - 1, s]
                )
                d <- 1 / sum(alphaStore[, h, ir])
                log_cStore[iGrid, h, ir] <- log(d)
                alphaStore[, h, ir ]<- d * alphaStore[, h, ir]
            }
        }
    }
    return(
        list(
            approach2_iRead = approach2_iRead,
            log_cStore = log_cStore,
            alphaStore = alphaStore,
            proposed_H = proposed_H,
            eMatGridLocalc = eMatGridLocalc
        )
    )
}

## class are
## 0 = original
## 1 = v1
## 2 = v2
## 3 = v3
## 4 = v12
## 5 = v13
## 6 = v23
## 7 = v123
check_agreements_with_truth_class <- function(true_H, true_H_class, H, rlc) {
    ##
    v <- array(0, c(7, 3)) ## what it can be
    rlc2 <- rlc
    rlc2[] <- as.integer(rlc2 > 0)
    check <- sapply(1:length(H), function(iRead) {
        hc <- true_H_class[iRead]
        if (hc == 0) {
            ## has to be the srR_ame
            return(true_H[iRead] == H[iRead])
        } else {
            ## can be a set
            ok <- which(rlc2[hc, ] == 1) ## these are OK
            return(H[iRead] %in% ok)
        }
    })
    x <- c(
        agree = sum(check),
        total = length(check),
        agree_class.gt0 = sum(check & (true_H_class > 0)),
        total_class.gt0 = sum((true_H_class > 0))
    )
    return(x["total_class.gt0"] - x["agree_class.gt0"])
}


## first dimension is the three options
## second dimension is the ir
## third dimension is the class
R_fill_rlcM <- function(rlcM, rlc, rr, verbose = FALSE) {
    ## rlcM <- array(NA, c(3, 6, 7))
    ## rlcI <- rlc
    ## make rlcI manually, based on assuming ff > 0
    rlcI <- array(integer(1), c(7, 3))
    rlcI[1, ] <- c(1L, 0L, 0L)
    rlcI[2, ] <- c(0L, 1L, 0L)
    rlcI[3, ] <- c(0L, 0L, 1L)
    rlcI[4, ] <- c(1L, 1L, 0L)
    rlcI[5, ] <- c(1L, 0L, 1L)
    rlcI[6, ] <- c(0L, 1L, 1L)
    rlcI[7, ] <- c(1L, 1L, 1L)            
    for(hcC in 1:7) {
        for(ir in 1:6) {
            if (verbose) {
                print(paste0("Currently working on current class_label hcC=", hcC, " and rr ir=", ir))
            }
            current_labels <- which(rlc[hcC, ] > 0)
            if (verbose) {
                print(paste0("Which labels this influences:", paste0(current_labels, collapse = ", ")))
            }
            future_labels <- rr[ir, current_labels]
            if (verbose) {
                print(paste0("With the remapping, this corresponds to the following relabels:", paste0(future_labels, collapse = ", ")))
            }
            x <- array(0, 3)
            x[future_labels] <- 1
            ## 
            hcF <- -1
            for(i in 1:7) {
                d <- 0
                for(j in 1:3) {
                    if ((rlcI[i, j] == 1) == (x[j] == 1)) {
                        d <- d + 1
                    }
                }
                if (d == 3) {
                    if (hcF != -1) {
                        stop("Faulty assumption in filling rlcM")
                    }
                    hcF <- i;
                }
            }
            if (hcF == -1) {
                stop("failed to set hcF")
            }
            ## hcF <- which(sapply(1:7, function(i) {sum(rlcI[i, ] == x) == 3}))
            if (verbose) {
                print(paste0("With the future class label of:", hcF))
            }
            if (length(hcF) > 1) {
                stop(paste0("bad logic with hcC = ", hcC, " and ir = ", ir))
            }
            rlcM[1:3, ir, hcC] <- rlc[hcF, rr[ir, ]]
            if (verbose) {
                print(paste0("So the current probabilities are:", paste0(round(rlc[hcC, ], 2), collapse = ", ")))
                print(paste0("And the future probabilities are:", paste0(round(rlc[hcF, rr[ir, ]], 2), collapse = ", ")))                
            }
        }
    }
    return(rlcM)
}



R_make_gibbs_considers <- function(
    blocked_snps,
    grid,
    wif0,
    nGrids,
    do_removal = TRUE, ## turn off for easy, partial checks
    verbose = FALSE
) {
    ##
    ##
    n_blocks <- max(blocked_snps) + 1
    nSNPs <- length(blocked_snps)
    nReads <- length(wif0)
    ##
    ## run a check / assert
    ## 
    for(iSNP in 1:(nSNPs - 1)) {
        if ((blocked_snps[iSNP + 1] - blocked_snps[iSNP]) > 1) {
            stop("problem with blocked snps")
        }
    }
    ##
    ##
    consider_snp_start_0_based <- integer(n_blocks)
    consider_snp_end_0_based <- integer(n_blocks)
    cur_block <- -1
    prev_block <- -1
    iBlock <- 1
    start <- 0
    for(iSNP in 1:nSNPs) {
        cur_block <- blocked_snps[iSNP]
        if (iSNP == nSNPs) {
            record <- TRUE
        } else if (blocked_snps[iSNP] < blocked_snps[iSNP + 1]) {
            record <- TRUE
        } else {
            record <- FALSE
        }
        if (record) {
            consider_snp_start_0_based[iBlock] <- start
            consider_snp_end_0_based[iBlock] <- iSNP - 1                        
            start <- iSNP - 1 + 1
            iBlock <- iBlock + 1
        }
    }
    ##
    ## grid
    ##
    consider_grid_start_0_based <- integer(n_blocks)
    consider_grid_end_0_based <- integer(n_blocks)
    for(iBlock in 1:n_blocks) {    
        consider_grid_start_0_based[iBlock] <- grid[consider_snp_start_0_based[iBlock] + 1]
        consider_grid_end_0_based[iBlock] <- grid[consider_snp_end_0_based[iBlock] + 1]
    }
    ##
    ## reads
    ##
    blocked_grid <- integer(nGrids)
    for(iBlock in 1:n_blocks) {
        a <- consider_grid_start_0_based[iBlock]
        b <- consider_grid_end_0_based[iBlock]
        for(i in a:b) {
            blocked_grid[i + 1] <- iBlock - 1
        }
    }
    ##
    ##
    ##
    consider_reads_start_0_based <- integer(n_blocks)
    consider_reads_start_0_based[] <- -1L
    consider_reads_end_0_based <- integer(n_blocks)
    consider_reads_end_0_based[] <- -1L
    ## here, all are 1-based
    previous_block_first_iRead <- 1
    previous_grid <- wif0[previous_block_first_iRead] + 1
    previous_block <- blocked_grid[previous_grid] + 1 ## 1-based
    ## you are here
    ## you need to slow these reads in correctly
    ## then put in cpp
    ## then put whole thing in cpp
    ## then do more tests!
    for(this_iRead in 2:(nReads)) {
        this_grid <- wif0[this_iRead] + 1
        this_block <- blocked_grid[this_grid] + 1
        if (this_iRead == nReads) {
            ## done, record
            consider_reads_start_0_based[this_block] <- previous_block_first_iRead - 1
            consider_reads_end_0_based[this_block] <- this_iRead - 1
        } else if (previous_block < this_block) {
            consider_reads_start_0_based[previous_block] <- previous_block_first_iRead - 1
            consider_reads_end_0_based[previous_block] <- (this_iRead - 1) - 1 ## second minus 1 is because it is before
            ## now reset previous
            previous_block_first_iRead <- this_iRead ## this read, 1-based
            previous_block <- blocked_grid[wif0[this_iRead] + 1] + 1
            previous_grid <- this_grid
        }
    }
    ##
    ## now, hard bit!
    ##
    if (do_removal) {
        remove <- consider_reads_start_0_based == -1
        w <- which(remove)
        if (sum(remove) > 0) {
            jBefore <- 1
            for(jNow in 1:length(w)) {
                if (jNow == length(w)) {
                    todo <- TRUE
                } else {
                    if ((w[jNow + 1] - w[jNow]) == 1) {
                        todo <- FALSE
                        jBefore <- jBefore - 1
                    } else {
                        todo <- TRUE
                    }
                }
                if (todo) {
                    if (verbose) {
                        print(paste0("todo, jBefore = ", jBefore, ", jNow = ", jNow))
                    }
                    s1 <- w[jBefore]
                    e1 <- w[jNow]
                    if (verbose) {
                        print(paste0("s1 = ", s1, ", e1 = ", e1))
                    }
                    x <- ceiling(0.5 * (consider_grid_start_0_based[s1] + consider_grid_end_0_based[e1]))
                    y <- ceiling(0.5 * (consider_snp_start_0_based[s1] + consider_snp_end_0_based[e1]))
                    ## if s1 is 1, then need to start at 0
                    if (s1 == 1) {
                        s1 <- 2
                        x <- 0
                        y <- 0
                    }
                    if (e1 == n_blocks) {
                        e1 <- e1 - 1
                        x <- consider_grid_end_0_based[n_blocks]
                        y <- consider_snp_end_0_based[n_blocks]
                    }
                    ## so push grid and read start / end
                    consider_grid_start_0_based[e1 + 1] <- x
                    consider_grid_end_0_based[s1 - 1] <- x - 1
                    ## same for snps
                    consider_snp_start_0_based[e1 + 1] <- y
                    consider_snp_end_0_based[s1 - 1] <- y - 1
                    jBefore <- jNow
                }
                jBefore <- jBefore + 1
            }
            ## 
            n_new_blocks <- sum(!remove)
            ##
            new_consider_reads_start_0_based <- integer(n_new_blocks)
            new_consider_reads_end_0_based <- integer(n_new_blocks)
            new_consider_grid_start_0_based <- integer(n_new_blocks)
            new_consider_grid_end_0_based <- integer(n_new_blocks)
            new_consider_snp_start_0_based <- integer(n_new_blocks)
            new_consider_snp_end_0_based <- integer(n_new_blocks)
            ## 
            i_prev_block <- 0
            for(iBlock in 1:n_blocks) {
                if (!remove[iBlock]) {
                    i_prev_block <- i_prev_block + 1
                    ##  reads
                    new_consider_reads_start_0_based[i_prev_block] <- consider_reads_start_0_based[iBlock]
                    new_consider_reads_end_0_based[i_prev_block] <- consider_reads_end_0_based[iBlock]
                    ## grid
                    new_consider_grid_start_0_based[i_prev_block] <- consider_grid_start_0_based[iBlock]
                    new_consider_grid_end_0_based[i_prev_block] <- consider_grid_end_0_based[iBlock]
                    ## snp
                    new_consider_snp_start_0_based[i_prev_block] <- consider_snp_start_0_based[iBlock]
                    new_consider_snp_end_0_based[i_prev_block] <- consider_snp_end_0_based[iBlock]
                }
            }
            ## reset
            consider_reads_start_0_based <- new_consider_reads_start_0_based
            consider_reads_end_0_based <- new_consider_reads_end_0_based
            consider_grid_start_0_based <- new_consider_grid_start_0_based
            consider_grid_end_0_based <- new_consider_grid_end_0_based
            consider_snp_start_0_based <- new_consider_snp_start_0_based
            consider_snp_end_0_based <- new_consider_snp_end_0_based
        }
    }
    ##
    ## finish
    ## 
    n_blocks <- length(consider_snp_end_0_based)
    consider_grid_where_0_based <- integer(nGrids)
    consider_grid_where_0_based[] <- -1L
    for(i in 1:n_blocks) {
        b <- consider_grid_end_0_based[i]
        consider_grid_where_0_based[b + 1] <- i - 1
    }
    return(
        list(
            consider_reads_start_0_based = consider_reads_start_0_based ,
            consider_reads_end_0_based  = consider_reads_end_0_based,
            consider_grid_start_0_based = consider_grid_start_0_based,            
            consider_grid_end_0_based = consider_grid_end_0_based,
            consider_snp_start_0_based = consider_snp_start_0_based,
            consider_snp_end_0_based = consider_snp_end_0_based,
            consider_grid_where_0_based = consider_grid_where_0_based,
            n_blocks = n_blocks
        )
    )
}











R_ff0_shard_block_gibbs_resampler <- function(
    alphaHat_t1,
    alphaHat_t2,
    alphaHat_t3,
    betaHat_t1,
    betaHat_t2,
    betaHat_t3,
    c1,
    c2,
    c3,
    eMatGrid_t1,
    eMatGrid_t2,
    eMatGrid_t3,
    H,
    eMatRead_t,    
    blocked_snps,
    grid,
    wif0,
    s,
    alphaMatCurrent_tc,
    priorCurrent_m,
    transMatRate_tc_H,
    do_checks = FALSE,
    initial_package = NULL,
    verbose = FALSE,
    fpp_stuff = NULL
) {

    ##
    ##
    K <- dim(alphaMatCurrent_tc)[1]
    nSNPs <- length(grid)
    nGrids <- dim(alphaMatCurrent_tc)[2] + 1    
    S <- dim(alphaMatCurrent_tc)[3]
    nReads <- length(wif0)


    ##
    ## compare arbitrary splits
    ##
    out <- Rcpp_make_gibbs_considers(
        blocked_snps = blocked_snps,
        grid = grid,
        wif0 = wif0,
        nGrids = nGrids
    )
    ## yes, take out from both
    consider_reads_start_0_based <- out[["consider_reads_start_0_based"]]
    consider_reads_end_0_based  <- out[["consider_reads_end_0_based"]]
    consider_grid_start_0_based <- out[["consider_grid_start_0_based"]]
    consider_grid_end_0_based <- out[["consider_grid_end_0_based"]]
    consider_snp_start_0_based <- out[["consider_snp_start_0_based"]]
    consider_snp_end_0_based <- out[["consider_snp_end_0_based"]]
    consider_grid_where_0_based <- out[["consider_grid_where_0_based"]]
    iGridConsider <- 0

    f_get_alt_prob_from_truth <- function(H, split_grid, wif0, fpp_stuff) {
        H_temp <- H
        H_temp[wif0 <= split_grid] <- 3 - H_temp[wif0 <= split_grid]
        alt_package <- for_testing_get_full_package_probabilities(H_temp, fpp_stuff)
        return(-sum(log(alt_package[[1]][["c"]]) + log(alt_package[[2]][["c"]])))
    }

    shard_block_columns <- c(
        "iBlock",
        "p_stay", "p_flip",
        "flip_mode",
        "p_O_stay", "p_O_flip"
    )
    shard_block_results <- matrix(0.0, nrow = length(consider_grid_end_0_based) - 1, ncol = length(shard_block_columns))
    colnames(shard_block_results) <- shard_block_columns
    ## 

    
    ##
    ## now do progression
    ##
    ## 
    minus_log_c1_sum <- 0
    minus_log_c2_sum <- 0
    in_flip_mode <- FALSE
    ## initialize normally
    iRead <- 0
    for(iGrid in 0:(nGrids - 1)) {
        ##
        ## normal forward one (includes initialization)
        ##
        if (iGrid == 0) {
            ## 
            alphaHat_t1[, iGrid + 1] <- priorCurrent_m[, s] * eMatGrid_t1[, iGrid + 1]
            c1[iGrid + 1] <- 1 / sum(alphaHat_t1[, iGrid + 1])
            alphaHat_t1[, iGrid + 1] <- alphaHat_t1[, iGrid + 1] * c1[iGrid + 1]
            ##
            alphaHat_t2[, iGrid + 1] <- priorCurrent_m[, s] * eMatGrid_t2[, iGrid + 1]
            c2[iGrid + 1] <- 1 / sum(alphaHat_t2[, iGrid + 1])
            alphaHat_t2[, iGrid + 1] <- alphaHat_t2[, iGrid + 1] * c2[iGrid + 1]
        } else {
            if (in_flip_mode) {
                x <- eMatGrid_t1[, iGrid + 1]
                eMatGrid_t1[, iGrid + 1] <- eMatGrid_t2[, iGrid + 1]
                eMatGrid_t2[, iGrid + 1] <- x
            }
            ## only thing that changes is eMatGrid_t
            rcpp_alpha_forward_one(s = 0, iGrid = iGrid, K = K, alphaHat_t = alphaHat_t1, transMatRate_tc_H = transMatRate_tc_H, eMatGrid_t = eMatGrid_t1, alphaMatCurrent_tc = alphaMatCurrent_tc, c = c1, minus_log_c_sum = minus_log_c1_sum, normalize = TRUE)
            rcpp_alpha_forward_one(s = 0, iGrid = iGrid, K = K, alphaHat_t = alphaHat_t2, transMatRate_tc_H = transMatRate_tc_H, eMatGrid_t = eMatGrid_t2, alphaMatCurrent_tc = alphaMatCurrent_tc, c = c2, minus_log_c_sum = minus_log_c2_sum, normalize = TRUE)
        }
        ##
        ##  go over read labels too, maybe flip them
        ##
        done <- FALSE
        while(!done) {
            if (iRead > (nReads - 1)) {
                done <- TRUE
            } else {
                if (wif0[iRead + 1] == iGrid) {
                    if (in_flip_mode) {
                        H[iRead + 1] <- 3 - H[iRead + 1]
                    }
                    iRead <- iRead + 1
                }
                if (iRead > (nReads - 1)) {
                    done <- TRUE
                } else {
                    if (wif0[iRead + 1] > iGrid) {
                        done <- TRUE
                    }
                }
            }
        }
        ##
        ## now do this bit
        ##
        if (consider_grid_where_0_based[iGrid + 1] > (- 1) && (iGrid < (nGrids - 1))) {
            iGridConsider <- consider_grid_where_0_based[iGrid + 1] ## still 0-based
            split_grid <- consider_grid_end_0_based[iGridConsider + 1]
            if (verbose) {
                print(paste0("Considering split_grid = ", split_grid))
            }
            if (do_checks) {
                ## first, get true current and alternate probs
                true_cur_prob <- f_get_alt_prob_from_truth(H, split_grid = -1, wif0, fpp_stuff)        
                true_alt_prob <- f_get_alt_prob_from_truth(H, split_grid, wif0, fpp_stuff)
                true_difference <- true_alt_prob - true_cur_prob
                current_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)
                w <- 1:(iGrid + 1)
                expect_equal(c1[w], current_package[[1]][["c"]][w])
                expect_equal(c2[w], current_package[[2]][["c"]][w])
                expect_equal(alphaHat_t1[, w], current_package[[1]][["alphaHat_t"]][, w])
                expect_equal(alphaHat_t2[, w], current_package[[2]][["alphaHat_t"]][, w])                
            }
            ##
            ## now do on fly version
            ## 
            ## alt
            ## now try to do "on fly"
            p_original <- -sum(log(c1) + log(c2))
            if (in_flip_mode) {
                x1 <- sum(alphaHat_t1[, split_grid + 1] * betaHat_t1[, split_grid + 1])
                x2 <- sum(alphaHat_t2[, split_grid + 1] * betaHat_t2[, split_grid + 1])
            } else {
                x1 <- sum(alphaHat_t1[, split_grid + 1] * betaHat_t2[, split_grid + 1])
                x2 <- sum(alphaHat_t2[, split_grid + 1] * betaHat_t1[, split_grid + 1])
            }
            calculated_difference <- log(x1) + log(x2) - log(c1[split_grid + 1]) - log(c2[split_grid + 1])
            if (in_flip_mode) {
                calculated_difference <- calculated_difference * -1
                p_original <- p_original - calculated_difference
            }
            p_alt <- p_original + calculated_difference
            ## yup, this looks right
            if (do_checks) {
                expect_equal(true_cur_prob, p_original)
                expect_equal(true_alt_prob, p_alt)
                expect_equal(true_difference, calculated_difference)
            }
            probs <- c(1, exp(calculated_difference))
            probs <- probs / sum(probs)
            if (verbose) {
                print(paste0("Remain the same prob is :", probs[1]))
            }
            ## now sample with respect to these
            in_flip_mode <- runif(1) > probs[1]
            ## record stuff now yyyyyyyyyeeeeeeeeeeeeeeeeessssssssssssss
            shard_block_results[iGridConsider + 1, "iBlock"] <- iGridConsider ## 0-based
            shard_block_results[iGridConsider + 1, c("p_stay", "p_flip")] <- probs
            shard_block_results[iGridConsider + 1, "flip_mode"] <- as.integer(in_flip_mode)
            shard_block_results[iGridConsider + 1, c("p_O_stay", "p_O_flip")] <- c(p_original, p_alt)
            ##block_results[iGridConsider + 1, "p_O1_given_H1_L"]
            ##
            if (verbose) {
                if (in_flip_mode) {
                    print(paste0("FLIP ME UP BRO"))
                } else {
                    print("no flipping for me thanks")
                }
            }
        }
    }
    
    ##
    ## re-run backward
    ## 
    betaHat_t1[, nGrids] <- c1[nGrids]
    Rcpp_run_backward_haploid(        
        betaHat_t1,
        c = c1,
        eMatGrid_t = eMatGrid_t1,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc = transMatRate_tc_H,
        s = s - 1
    )
    betaHat_t2[, nGrids] <- c2[nGrids]
    Rcpp_run_backward_haploid(        
        betaHat_t2,
        c = c2,
        eMatGrid_t = eMatGrid_t2,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc = transMatRate_tc_H,
        s = s - 1
    )
    betaHat_t3[, nGrids] <- c3[nGrids]
    Rcpp_run_backward_haploid(        
        betaHat_t3,
        c = c3,
        eMatGrid_t = eMatGrid_t3,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc = transMatRate_tc_H,
        s = s - 1
    )
    return(
        list(
            eMatGrid_t1 = eMatGrid_t1,
            eMatGrid_t2 = eMatGrid_t2,
            eMatGrid_t3 = eMatGrid_t3,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            alphaHat_t1 = alphaHat_t1,
            alphaHat_t2 = alphaHat_t2,
            alphaHat_t3 = alphaHat_t3,            
            betaHat_t1 = betaHat_t1,
            betaHat_t2 = betaHat_t2,
            betaHat_t3 = betaHat_t3,
            H = H,
            fpp_stuff = fpp_stuff,
            shard_block_results = shard_block_results,
            consider_snp_start_0_based = consider_snp_start_0_based,            
            consider_snp_end_0_based = consider_snp_end_0_based,
            consider_reads_start_0_based = consider_reads_start_0_based,
            consider_reads_end_0_based = consider_reads_end_0_based
        )
    )
}


