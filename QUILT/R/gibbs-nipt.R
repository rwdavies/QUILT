## currently will only work for S=1
## true gibbs sampling for NIPT (three haploid) setting
## borrow heavily from existing code!
forwardBackwardGibbsNIPT <- function(
    sampleReads,
    priorCurrent_m,
    alphaMatCurrent_tc,
    eHapsCurrent_tc,
    transMatRate_tc_H,
    ff,
    blocks_for_output,   ##  = array(0, c(1, 1))
    Jmax_local = 100,
    maxDifferenceBetweenReads = 1000,
    maxEmissionMatrixDifference = 1e10,
    run_fb_subset = FALSE,
    run_fb_grid_offset = 0,
    return_genProbs = TRUE,
    return_hapProbs = TRUE,
    return_gamma = TRUE,
    return_alpha = FALSE,
    grid = -1,
    snp_start_1_based = NA,
    snp_end_1_based = NA,
    generate_fb_snp_offsets = FALSE,
    alphaStart = array(0, 1),
    betaEnd = array(0, 1),
    suppressOutput = 0, ## for consistency
    n_gibbs_starts = 1,    
    n_gibbs_sample_its = 1,
    n_gibbs_burn_in_its = 1,    
    use_starting_read_labels = FALSE,
    return_p_store = FALSE,
    seed = NA,
    verbose = FALSE,
    true_H = NULL,
    double_list_of_starting_read_labels = NULL,
    seed_vector = -1,
    rescale_eMatRead_t = TRUE,
    bound_eMatGrid_t = FALSE,
    rescale_eMatGrid_t = FALSE,
    haploid_gibbs_equal_weighting = TRUE,
    return_extra = FALSE,
    gibbs_initialize_iteratively = FALSE,
    gibbs_initialize_at_first_read = TRUE,
    do_block_resampling = FALSE,
    artificial_relabel = -1,
    alphaHat_t1 = array(0, c(1, 1)),
    alphaHat_t2 = array(0, c(1, 1)),
    alphaHat_t3 = array(0, c(1, 1)),
    betaHat_t1 = array(0, c(1, 1)),
    betaHat_t2 = array(0, c(1, 1)),
    betaHat_t3 = array(0, c(1, 1)),
    hapSum_tc = array(0, c(1, 1, 1)),
    update_hapSum_tc = FALSE,
    record_read_set = FALSE,
    class_sum_cutoff = 0.06,
    wif0 = integer(1),
    L_grid = integer(1),
    prev_list_of_alphaBetaBlocks = NULL,
    i_snp_block_for_alpha_beta = 1,
    perform_block_gibbs = FALSE,
    shuffle_bin_radius = NULL,
    block_gibbs_iterations = NULL,
    return_gibbs_block_output = NULL,
    force_reset_read_category_zero = TRUE
) {
    if (!is.na(seed)) {
        set.seed(seed)
    }
    ##
    n_gibbs_full_its <- n_gibbs_burn_in_its + n_gibbs_sample_its
    ## initialize reads
    nReads <- length(sampleReads)
    nGrids <- ncol(transMatRate_tc_H) + 1
    K <- nrow(eHapsCurrent_tc)
    S <- dim(eHapsCurrent_tc)[3]
    s <- 1
    list_of_ending_read_labels <- list()
    ##
    if (return_p_store) {
        p_store_cols <- c(
            "p_1", "p_2", "p_3",
            "chance", "h_rC", "h_rN",
            "agreePer",
            "itType",
            "p_O1_given_H1_L",
            "p_O2_given_H2_L",
            "p_O3_given_H3_L",
            "p_O_given_H_L",
            "p_H_given_L",
            "p_H_given_O_L_up_to_C"
        )
        p_store <- array(0, c(n_gibbs_starts * n_gibbs_full_its * nReads, length(p_store_cols)))
        colnames(p_store) <- p_store_cols
    } else {
        p_store <- array(0, c(2, 2))
    }
    ##
    ## for re-weighting genProbs
    prior_probs <- c(0.5, (1 - ff) / 2, ff / 2)
    p <- prior_probs
    if (record_read_set) {
        ## read label classifier
        rlc <- make_rlc(ff)
        H_class <- integer(nReads)        
        ##
    } else {
        rlc <- matrix(0, nrow = 1, ncol = 1)
        H_class <- integer(1)
    }
    ## 
    if (n_gibbs_full_its > 0) {
        hg_ll_rescaled <- array(-1, n_gibbs_starts * n_gibbs_sample_its)
    } else {
        hg_ll_rescaled <- array(-1, n_gibbs_starts)
    }
    hg_log_mult <- array(-1, 2) ## only first used. two spots as weird Rcpp bug
    if (is.na(snp_start_1_based)) {
        snp_start_1_based <- 1
        snp_end_1_based <- length(grid)
    }
    ##
    ## this is done ocne
    ##
    eMatRead_t <- array(1, c(K, nReads))
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = eHapsCurrent_tc,
        s = s - 1,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = Jmax_local,
        eMatHapOri_t = array(0, c(1, 1)),
        pRgivenH1 = array(0),
        pRgivenH2 = array(0),
        prev = 0,
        suppressOutput = 1,
        prev_section = "",
        next_section = "",
        run_pseudo_haploid = FALSE,
        rescale_eMatRead_t = rescale_eMatRead_t
    )
    ## 
    ## evaluate how necessary they are
    ##
    outR <- evaluate_read_variability(eMatRead_t)
    number_of_non_1_reads <- outR[["number_of_non_1_reads"]]
    indices_of_non_1_reads <- outR[["indices_of_non_1_reads"]]
    read_category <- outR[["read_category"]]
    ##
    if (force_reset_read_category_zero) {
        read_category[read_category != 1L] <- 0L
    }
    ##
    to_return <- list()
    if (return_extra) {
        to_return[["per_iteration_store"]] <- lapply(1:(n_gibbs_starts * n_gibbs_sample_its), function(x) as.list(1:2))
    }
    ##to_return[["temp"]] <- lapply(1:n_gibbs_starts, function(x) as.list(1:3))
    ##
    ## initialize and iterate in its own function
    ##
    ## column names are P(O^1 | \lambda) for (1, 2, 3), then P(H | \lambda), then P(O, H | \lambda) aka the likelihood (for H) ll = log10(P(H | O, \lambda))
    cols <- c(
        "s", "i_samp", "i_it", "i_result_it",
        "p_O1_given_H1_L", "p_O2_given_H2_L", "p_O3_given_H3_L",
        "p_O_given_H_L", "p_H_given_L", "p_H_given_O_L_up_to_C",
        "p_set_H_given_L", "relabel"
    )
    per_it_likelihoods <- array(0, c(S * n_gibbs_starts * n_gibbs_full_its, length(cols)))
    colnames(per_it_likelihoods) <- cols
    ##
    if (!run_fb_subset) {
        n_outer <- n_gibbs_starts;
    } else {
        n_outer <- n_gibbs_starts * n_gibbs_sample_its;
    }
    for(i_outer in 1:n_outer) {
        if (!run_fb_subset) {
            i_gibbs_samplings <- i_outer; ## 1-based
        } else {
            i_gibbs_samplings <- floor((i_outer - 1) / n_gibbs_sample_its) + 1 ## 1-based
        }
        if (seed_vector[1] >= 0) {
            set.seed(seed_vector[i_gibbs_samplings])
        }
        ## why did this not screw things up?
        runif_reads <- runif(nReads * n_gibbs_full_its)
        if (gibbs_initialize_at_first_read) {
            first_read_for_gibbs_initialization <- 1
        } else {
            first_read_for_gibbs_initialization <- sample(nReads, 1)
        }
        ##
        if (use_starting_read_labels) {
            H <- double_list_of_starting_read_labels[[s]][[i_outer]]
        } else {
            H <- random_gibbs_nipt_read_labels(nReads, ff)
        }
        ## 
        ## initialize - yuck - should be its own function in R, but meh
        ##
        if (gibbs_initialize_iteratively) {
            ##
            ## here, the initialization is done on the fly, using the first iteration
            ## much simpler!
            ##
            nGrids <- ncol(alphaMatCurrent_tc) + 1
            eMatGrid_t1 <- array(1, c(K, nGrids))
            eMatGrid_t2 <- array(1, c(K, nGrids))
            eMatGrid_t3 <- array(1, c(K, nGrids))
            ## initialize alpha, beta suuuuper simple
            alphaHat_t1 <- array(1, c(K, nGrids)); betaHat_t1 <- array(1, c(K, nGrids)); c1 <- array(1, nGrids)
            alphaHat_t2 <- array(1, c(K, nGrids)); betaHat_t2 <- array(1, c(K, nGrids)); c2 <- array(1, nGrids)
            alphaHat_t3 <- array(1, c(K, nGrids)); betaHat_t3 <- array(1, c(K, nGrids)); c3 <- array(1, nGrids)
            ##
            ##
            ##
            Rcpp_run_forward_haploid(alphaHat_t = alphaHat_t1, c = c1, eMatGrid_t = eMatGrid_t1, alphaMatCurrent_tc = alphaMatCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, priorCurrent_m  = priorCurrent_m, s = s - 1, alphaStart = 0, run_fb_subset = FALSE , initialize_only = TRUE);
            Rcpp_run_forward_haploid(alphaHat_t = alphaHat_t2, c = c2, eMatGrid_t = eMatGrid_t2, alphaMatCurrent_tc = alphaMatCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, priorCurrent_m  = priorCurrent_m, s = s - 1, alphaStart = 0, run_fb_subset = FALSE , initialize_only = TRUE);
            Rcpp_run_forward_haploid(alphaHat_t = alphaHat_t3, c = c3, eMatGrid_t = eMatGrid_t3, alphaMatCurrent_tc = alphaMatCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, priorCurrent_m  = priorCurrent_m, s = s - 1, alphaStart = 0, run_fb_subset = FALSE , initialize_only = TRUE);
            ##
            ##
            ##
        } else {
            ##
            eMatGrid_t1 <- array(1, c(K, nGrids))
            eMatGrid_t2 <- array(1, c(K, nGrids))
            eMatGrid_t3 <- array(1, c(K, nGrids))
            ## 
            rcpp_make_eMatGrid_t(eMatGrid_t = eMatGrid_t1, eMatRead_t = eMatRead_t, H = H, sampleReads = sampleReads, hap = 1, nGrids = nGrids, prev = 0, suppressOutput = 1, prev_section = "text", next_section = "", run_fb_grid_offset = run_fb_grid_offset, use_all_reads = FALSE, bound = bound_eMatGrid_t, maxEmissionMatrixDifference =  maxEmissionMatrixDifference, rescale = rescale_eMatGrid_t)
            ##
            rcpp_make_eMatGrid_t(eMatGrid_t = eMatGrid_t2, eMatRead_t = eMatRead_t, H = H, sampleReads = sampleReads, hap = 2, nGrids = nGrids, prev = 0, suppressOutput = 1, prev_section = "text", next_section = "", run_fb_grid_offset = run_fb_grid_offset, use_all_reads = FALSE, bound = bound_eMatGrid_t, maxEmissionMatrixDifference =  maxEmissionMatrixDifference, rescale = rescale_eMatGrid_t)
            ## 
            rcpp_make_eMatGrid_t(eMatGrid_t = eMatGrid_t3, eMatRead_t = eMatRead_t, H = H, sampleReads = sampleReads, hap = 3, nGrids = nGrids, prev = 0, suppressOutput = 1, prev_section = "text", next_section = "", run_fb_grid_offset = run_fb_grid_offset, use_all_reads = FALSE, bound = bound_eMatGrid_t, maxEmissionMatrixDifference =  maxEmissionMatrixDifference, rescale = rescale_eMatGrid_t)
            ##
            ##
            ## 
            ##
            out <- initialize_gibbs_forward_backward(H = H, hap = 1, s = s, sampleReads = sampleReads, priorCurrent_m = priorCurrent_m, alphaMatCurrent_tc = alphaMatCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, eMatRead_t = eMatRead_t, eMatGrid_t = eMatGrid_t1)
            alphaHat_t1 <- out$alphaHat_t; betaHat_t1 <- out$betaHat_t; c1 <- out$c
            ##
            out <- initialize_gibbs_forward_backward(H = H, hap = 2, s = s, sampleReads = sampleReads, priorCurrent_m = priorCurrent_m, alphaMatCurrent_tc = alphaMatCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, eMatRead_t = eMatRead_t, eMatGrid_t = eMatGrid_t2)
            alphaHat_t2 <- out$alphaHat_t; betaHat_t2 <- out$betaHat_t; c2 <- out$c
            ##
            out <- initialize_gibbs_forward_backward(H = H, hap = 3, s = s, sampleReads = sampleReads, priorCurrent_m = priorCurrent_m, alphaMatCurrent_tc = alphaMatCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, eMatRead_t = eMatRead_t, eMatGrid_t = eMatGrid_t3)
            alphaHat_t3 <- out$alphaHat_t; betaHat_t3 <- out$betaHat_t; c3 <- out$c
        }
        ## 
        ## iterate, and possibly save!
        ## 
        prior_probs <- c(0.5, (1 - ff) / 2, ff / 2)
        ##
        if (n_gibbs_full_its > 0) {
            iteration_run_list <- 1:n_gibbs_full_its
        } else {
            iteration_run_list <- NULL
            i_result_it <- i_outer
            i_ever_it <- i_outer
            out <- save_various_gammas(return_gamma = return_gamma, return_genProbs = return_genProbs, return_hapProbs = return_hapProbs, return_extra = return_extra, H = H, alphaHat_t1 = alphaHat_t1, betaHat_t1 = betaHat_t1, c1 = c1, alphaHat_t2 = alphaHat_t2, betaHat_t2 = betaHat_t2, c2 = c2, alphaHat_t3 = alphaHat_t3, betaHat_t3 = betaHat_t3, c3 = c3, eHapsCurrent_tc = eHapsCurrent_tc, s = s, grid = grid, snp_start_1_based = snp_start_1_based, snp_end_1_based = snp_end_1_based, run_fb_grid_offset = run_fb_grid_offset, i_gibbs_samplings = i_gibbs_samplings, i_result_it = i_result_it, to_return = to_return, hg_log_mult = hg_log_mult, hg_ll_rescaled = hg_ll_rescaled, haploid_gibbs_equal_weighting = haploid_gibbs_equal_weighting, prior_probs = prior_probs, ff = ff)
            to_return <- out$to_return
            hg_log_mult <- out$hg_log_mult
            hg_ll_rescaled <- out$hg_ll_rescaled
        }
        ##
        for(iteration in iteration_run_list) {
            if (run_fb_subset) {
                i_result_it <- i_outer
                i_ever_it <- i_outer
            } else {
                if ((iteration) > n_gibbs_burn_in_its) {
                    i_result_it <- n_gibbs_sample_its * (i_gibbs_samplings - 1) + iteration - n_gibbs_burn_in_its;
                } else {
                    i_result_it <- -1;
                }
                ## this is e.g. 
                i_ever_it <- n_gibbs_full_its * (i_gibbs_samplings - 1) + iteration
            }
            if (verbose) {
                print(paste0("i_gibbs_samplings = ", i_gibbs_samplings))
                print(paste0("iteration = ", iteration))
                print(paste0("i_ever_it = ", i_ever_it))
                print(paste0("i_result_it = ", i_result_it))
            }
            ##
            ## do an iteration here
            ##
            out <- gibbs_nipt_one_iteration(s = s, iteration = iteration, sampleReads = sampleReads, priorCurrent_m = priorCurrent_m, transMatRate_tc_H = transMatRate_tc_H, alphaMatCurrent_tc = alphaMatCurrent_tc, n_gibbs_full_its = n_gibbs_full_its, nGrids = nGrids,K = K, H = H, eMatRead_t = eMatRead_t, runif_reads = runif_reads, alphaHat_t1 = alphaHat_t1, betaHat_t1 = betaHat_t1, c1 = c1, eMatGrid_t1 = eMatGrid_t1, alphaHat_t2 = alphaHat_t2, betaHat_t2 = betaHat_t2, c2 = c2, eMatGrid_t2 = eMatGrid_t2, alphaHat_t3 = alphaHat_t3, betaHat_t3 = betaHat_t3, c3 = c3, eMatGrid_t3 = eMatGrid_t3, p_store = p_store, verbose = verbose, return_p_store = return_p_store, true_H = true_H, i_gibbs_samplings = i_gibbs_samplings, ff = ff, bound_eMatGrid_t = bound_eMatGrid_t, per_it_likelihoods = per_it_likelihoods, gibbs_initialize_iteratively = gibbs_initialize_iteratively, first_read_for_gibbs_initialization = first_read_for_gibbs_initialization, do_block_resampling = do_block_resampling, artificial_relabel = artificial_relabel, prior_probs = prior_probs, i_result_it = i_result_it, record_read_set = record_read_set, H_class = H_class, rlc = rlc, class_sum_cutoff = class_sum_cutoff, number_of_non_1_reads = number_of_non_1_reads, indices_of_non_1_reads  = indices_of_non_1_reads, read_category = read_category)
            ## unpack
            alphaHat_t1 <- out$alphaHat_t1; betaHat_t1 <- out$betaHat_t1; c1 <- out$c1; eMatGrid_t1 <- out$eMatGrid_t1
            alphaHat_t2 <- out$alphaHat_t2; betaHat_t2 <- out$betaHat_t2; c2 <- out$c2; eMatGrid_t2 <- out$eMatGrid_t2
            alphaHat_t3 <- out$alphaHat_t3; betaHat_t3 <- out$betaHat_t3; c3 <- out$c3; eMatGrid_t3 <- out$eMatGrid_t3
            p_store <- out$p_store
            H <- out[["H"]]
            per_it_likelihoods <- out$per_it_likelihoods
            H_class <- out[["H_class"]]
            ##
            ## if this is a sampling iteration, save! add to gammas, etc
            ## 
            if (run_fb_subset | ((iteration) > n_gibbs_burn_in_its)) {
                ## if this is fb subset, we are doing this every time
                ## otherwise, it depends!
                ## 
                out <- save_various_gammas(return_gamma = return_gamma, return_genProbs = return_genProbs, return_hapProbs = return_hapProbs, return_extra = return_extra, H = H, alphaHat_t1 = alphaHat_t1, betaHat_t1 = betaHat_t1, c1 = c1, alphaHat_t2 = alphaHat_t2, betaHat_t2 = betaHat_t2, c2 = c2, alphaHat_t3 = alphaHat_t3, betaHat_t3 = betaHat_t3, c3 = c3, eHapsCurrent_tc = eHapsCurrent_tc, s = s, grid = grid, snp_start_1_based = snp_start_1_based, snp_end_1_based = snp_end_1_based, run_fb_grid_offset = run_fb_grid_offset, i_gibbs_samplings = i_gibbs_samplings, i_result_it = i_result_it, to_return = to_return, hg_log_mult = hg_log_mult, hg_ll_rescaled =hg_ll_rescaled, haploid_gibbs_equal_weighting = haploid_gibbs_equal_weighting, prior_probs = prior_probs, ff = ff)
                to_return <- out$to_return
                hg_log_mult <- out$hg_log_mult
                hg_ll_rescaled <- out$hg_ll_rescaled
                if (!run_fb_subset) {
                    list_of_ending_read_labels <- append(list_of_ending_read_labels, list(H = H))
                }
            }
        }
        ##
        ## end of the iterations
        ##
        ##
        to_return[["H"]] <- H ## this means certain things will be overridden
        to_return[["H_class"]] <- H_class
        if (return_p_store) {
            to_return[["p_store"]] <- p_store
        }
        if (return_alpha) {
            to_return[["alphaHat_t1"]] <- alphaHat_t1
            to_return[["alphaHat_t2"]] <- alphaHat_t2
            to_return[["alphaHat_t3"]] <- alphaHat_t3
            to_return[["betaHat_t1"]] <- betaHat_t1
            to_return[["betaHat_t2"]] <- betaHat_t2
            to_return[["betaHat_t3"]] <- betaHat_t3
            to_return[["c1"]] <- c1
            to_return[["c2"]] <- c2
            to_return[["c3"]] <- c3
        }
    }
    ## finish!
    for(object in
        c(
            "genProbsM_t", "genProbsF_t", "hapProbs_t",
            "gammaMT_t", "gammaMU_t", "gammaP_t",
            "gammaM_t", "gammaF_t"
        )) {
        to_return[[object]] <- to_return[[object]] / sum(exp(hg_ll_rescaled))
    }
    ## add likelihoods - always useful!
    to_return[["hg_ll_rescaled"]] <- hg_ll_rescaled
    to_return[["double_list_of_ending_read_labels"]] <- list(list_of_ending_read_labels)
    to_return[["per_it_likelihoods"]] <- per_it_likelihoods
    ##
    if (return_extra) {
        to_return[["eMatRead_t"]] <- eMatRead_t
        to_return[["eMatGrid_t1"]] <- eMatGrid_t1
        to_return[["eMatGrid_t2"]] <- eMatGrid_t2
        to_return[["eMatGrid_t3"]] <- eMatGrid_t3        
    }
    return(to_return)
}


save_various_gammas <- function(
    return_gamma,
    return_genProbs,
    return_hapProbs,    
    return_extra,
    H,
    alphaHat_t1,
    betaHat_t1,
    c1,
    alphaHat_t2,
    betaHat_t2,
    c2,
    alphaHat_t3,
    betaHat_t3,
    c3,
    eHapsCurrent_tc,
    s,
    grid,
    snp_start_1_based,
    snp_end_1_based,
    run_fb_grid_offset,
    i_gibbs_samplings,
    i_result_it,
    to_return,
    hg_log_mult,
    hg_ll_rescaled,
    haploid_gibbs_equal_weighting,
    prior_probs,
    ff,
    sample_is_diploid = FALSE
) {
    if (return_gamma | return_genProbs | return_hapProbs) {
        read_counts <- c(sum(H == 1), sum(H == 2), sum(H == 3))
        prob_matrix <- determine_label_probabilities(read_counts, ff)
        ## do separate?
        gamma1 <- (alphaHat_t1 %*% diag(1 / c1)) * betaHat_t1
        gamma2 <- (alphaHat_t2 %*% diag(1 / c2)) * betaHat_t2
        gamma3 <- (alphaHat_t3 %*% diag(1 / c3)) * betaHat_t3
        ## now go!
        gammaMT_t <-
            prob_matrix[1, 1] * gamma1 +
            prob_matrix[1, 2] * gamma2 +
            prob_matrix[1, 3] * gamma3
        gammaMU_t <-
            prob_matrix[2, 1] * gamma1 +
            prob_matrix[2, 2] * gamma2 +
            prob_matrix[2, 3] * gamma3
        gammaP_t <-
            prob_matrix[3, 1] * gamma1 +
            prob_matrix[3, 2] * gamma2 +
            prob_matrix[3, 3] * gamma3
        ##
        ## genProbs here
        ##
        nSNPs <- ncol(eHapsCurrent_tc)
        genProbsM_t <- array(0, c(3, nSNPs))
        genProbsF_t <- array(0, c(3, nSNPs))
        hapProbs_t <- array(0, c(3, nSNPs))
        rcpp_calculate_gn_genProbs_and_hapProbs(
            genProbsM_t = genProbsM_t,
            genProbsF_t = genProbsF_t,
            hapProbs_t = hapProbs_t,
            eHapsCurrent_tc = eHapsCurrent_tc,
            s = s - 1,
            gammaMT_t = gammaMT_t,
            gammaMU_t = gammaMU_t,
            gammaP_t = gammaP_t,            
            grid = grid,
            snp_start_1_based = snp_start_1_based,
            snp_end_1_based = snp_end_1_based,
            grid_offset = run_fb_grid_offset,
            sample_is_diploid = sample_is_diploid
        )
        ## 
        if (return_genProbs | return_hapProbs) {
            ##
            ll_current <- calculate_likelihoods_values(c1, c2, c3, H, nGrids = length(c1), prior_probs, ff)[6]
            ## 
            out <- fly_weighter(
                i = i_result_it,
                gCurrent = genProbsM_t,
                gAverage = to_return[["genProbsM_t"]],
                ll_current = ll_current,
                log_mult = hg_log_mult,
                ll_rescaled = hg_ll_rescaled,
                gCurrent2 = genProbsF_t,
                gAverage2 = to_return[["genProbsF_t"]],
                gCurrent3 = hapProbs_t,
                gAverage3 = to_return[["hapProbs_t"]],
                equal_weighting = haploid_gibbs_equal_weighting
            )
            to_return[["genProbsM_t"]] <- out$gAverage
            to_return[["genProbsF_t"]] <- out$gAverage2
            to_return[["hapProbs_t"]] <- out$gAverage3
            hg_log_mult <- out$log_mult
            hg_ll_rescaled <- out$ll_rescaled
            log_rescale <- out$log_rescale
            relative_difference <- out$relative_difference
            ##
            if (return_extra) {
                ## store every genProbs here
                to_return[["per_iteration_store"]][[i_result_it]]$genProbsM_t <- genProbsM_t
                to_return[["per_iteration_store"]][[i_result_it]]$genProbsF_t <- genProbsF_t
                to_return[["per_iteration_store"]][[i_result_it]]$hapProbs_t <- hapProbs_t                
            }
        }
        if (return_gamma | return_hapProbs) {
            ## ugh
            gammaM_t <- gammaMT_t + gammaMU_t
            gammaF_t <- gammaMT_t + gammaP_t
            ##
            ## need to decide how to normalize!
            to_return[["gammaMT_t"]] <- single_fly_weight(i_result_it, relative_difference, gammaMT_t, to_return[["gammaMT_t"]], log_rescale)
            to_return[["gammaMU_t"]] <- single_fly_weight(i_result_it, relative_difference, gammaMU_t, to_return[["gammaMU_t"]], log_rescale)
            to_return[["gammaP_t"]] <- single_fly_weight(i_result_it, relative_difference, gammaP_t, to_return[["gammaP_t"]], log_rescale)
            to_return[["gammaM_t"]] <- single_fly_weight(i_result_it, relative_difference, gammaM_t, to_return[["gammaM_t"]], log_rescale)
            to_return[["gammaF_t"]] <- single_fly_weight(i_result_it, relative_difference, gammaF_t, to_return[["gammaF_t"]], log_rescale)
        }
    }
    return(
        list(
            to_return = to_return,
            hg_log_mult = hg_log_mult,
            hg_ll_rescaled = hg_ll_rescaled
        )
    )
}




random_gibbs_nipt_read_labels <- function(nReads, ff) {
    x <- runif(nReads)
    H <- array(0L, nReads)
    ## mt = 0            -> 0.5
    ## mu = 0.5          -> 0.5 + ff / 2
    ## p  = 0.5 + ff / 2 -> 1
    for(iRead in 1:nReads) {
        if (x[iRead] < 0.5) {
            H[iRead] = 1
        } else if ((0.5 <= x[iRead]) & (x[iRead] < (0.5 + ff / 2))) {
            H[iRead] = 2
        } else {
            H[iRead] = 3
        }
    }
    return(H)
}


gibbs_nipt_one_iteration <- function(
    s,
    iteration,
    sampleReads,
    priorCurrent_m,
    transMatRate_tc_H,
    alphaMatCurrent_tc,
    n_gibbs_full_its,
    nGrids,
    K,
    H,
    eMatRead_t,
    runif_reads,
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
    p_store,
    verbose,
    return_p_store,
    true_H,
    i_gibbs_samplings,
    ff,
    bound_eMatGrid_t,
    per_it_likelihoods,
    gibbs_initialize_iteratively = FALSE,
    first_read_for_gibbs_initialization = 1,
    do_block_resampling = FALSE,
    artificial_relabel = -1,
    prior_probs,
    i_result_it,
    record_read_set,
    H_class,
    rlc,
    class_sum_cutoff,
    number_of_non_1_reads,
    indices_of_non_1_reads,
    read_category
) {
    ##
    nReads <- length(sampleReads)
    ##
    ##
    ##
    ##
    done_reads <- FALSE
    iRead <- 0
    iGrid <- 1
    for(iGrid in 1:nGrids) {
        if (verbose) {
            print(paste0("iGrid = ", iGrid))
        }
        ## build new alphaHat_m
        if (iGrid > 1) {
            ## does not include c here
            alphaHat_t1 <- alpha_forward_one(iGrid, K, s, alphaHat_t1, transMatRate_tc_H, eMatGrid_t1, alphaMatCurrent_tc)
            alphaHat_t2 <- alpha_forward_one(iGrid, K, s, alphaHat_t2, transMatRate_tc_H, eMatGrid_t2, alphaMatCurrent_tc)
            alphaHat_t3 <- alpha_forward_one(iGrid, K, s, alphaHat_t3, transMatRate_tc_H, eMatGrid_t3, alphaMatCurrent_tc)
            ## now, do previous normalization
            alphaHat_t1[, iGrid] <- alphaHat_t1[, iGrid] * c1[iGrid]
            a <- 1 / sum(alphaHat_t1[, iGrid])
            c1[iGrid] <- c1[iGrid] * a
            alphaHat_t1[, iGrid] <- alphaHat_t1[, iGrid] * a
            ##
            alphaHat_t2[, iGrid] <- alphaHat_t2[, iGrid] * c2[iGrid]
            a <- 1 / sum(alphaHat_t2[, iGrid])
            c2[iGrid] <- c2[iGrid] * a
            alphaHat_t2[, iGrid] <- alphaHat_t2[, iGrid] * a
            ##
            alphaHat_t3[, iGrid] <- alphaHat_t3[, iGrid] * c3[iGrid]
            a <- 1 / sum(alphaHat_t3[, iGrid])
            c3[iGrid] <- c3[iGrid] * a
            alphaHat_t3[, iGrid] <- alphaHat_t3[, iGrid] * a
        } else {
            rcpp_reinitialize_in_iterations(s - 1, alphaHat_t1, c1, priorCurrent_m, eMatGrid_t1, K)
            rcpp_reinitialize_in_iterations(s - 1, alphaHat_t2, c2, priorCurrent_m, eMatGrid_t2, K)
            rcpp_reinitialize_in_iterations(s - 1, alphaHat_t3, c3, priorCurrent_m, eMatGrid_t3, K)
        }
        ## 
        alphaHat_m <- rbind(alphaHat_t1[, iGrid], alphaHat_t2[, iGrid], alphaHat_t3[, iGrid])
        betaHat_m <- rbind(betaHat_t1[, iGrid], betaHat_t2[, iGrid], betaHat_t3[, iGrid])
        pC <- rowSums(alphaHat_m * betaHat_m) ## 
        if (verbose) {
            print(paste0("set pC = ", paste0(pC, collapse = ", ")))
        }
        ##
        iRead <- iRead + 1
        iReadStart <- iRead
        if (!done_reads) {
            read_wif_iRead <- sampleReads[[iRead]][[2]];
        } else {
            read_wif_iRead = -1;
        }
        ## set three booleans about progress type
        ##
        ##
        ##
        while((done_reads == FALSE) & (read_wif_iRead + 1) == iGrid) {
            ##
            ## don't bother if read is uninformative
            ##
            if (read_category[iRead] != 1L) {
                ##
                ## decide on what type of work is being done
                ##
                currently_doing_normal_progress <- FALSE
                currently_doing_gibbs_initialization <- FALSE
                currently_doing_pass_through <- FALSE
                if (!gibbs_initialize_iteratively) {
                    currently_doing_normal_progress <- TRUE
                } else {
                    if ((iRead < first_read_for_gibbs_initialization) & (iteration == 1)) {
                        currently_doing_pass_through <- TRUE
                    } else if ((first_read_for_gibbs_initialization <= iRead) & (iteration == 1)) {
                        currently_doing_gibbs_initialization <- TRUE
                    } else if ((iRead < first_read_for_gibbs_initialization) & (iteration == 2)) {
                        currently_doing_gibbs_initialization <- TRUE
                    } else {
                        currently_doing_normal_progress <- TRUE
                    }
                }
                ## in R, make a check, always, exactly 1 should be chosen
                x <- c(currently_doing_normal_progress, currently_doing_gibbs_initialization, currently_doing_pass_through)
                if (sum(x) != 1) {
                    print(x)
                    stop("Bad initialization choice")
                }
                ##
                if (verbose) {
                    print(paste0("------------------------ it = ", iteration, ", iRead = ", iRead))
                }
                ##
                if (currently_doing_normal_progress) {
                    ## i_remove = 0 is remove it from the haplotype that normally has it
                    ## i_remove = 1 is add it to the haplotype that is (temporarily) gaining it
                    ## these are 1-3 based as normal here
                    h_rC <- H[iRead] ## h_r current
                    h_rA1 <- setdiff(1:3, H[iRead])[1] ## h_r alternate 1
                    h_rA2 <- setdiff(1:3, H[iRead])[2] ## h_r alternate 1
                    if (verbose) {
                        print(paste0("h_rC = ", h_rC, ", h_rA1 = ", h_rA1, ", h_rA2 = ", h_rA2))
                    }
                    out <- evaluate_read_probabilities(
                        alphaHat_m = alphaHat_m,
                        betaHat_m = betaHat_m, 
                        pC = pC,
                        read_category = read_category,
                        iRead = iRead,
                        h_rC = h_rC,
                        h_rA1 = h_rA1,
                        h_rA2 = h_rA2,
                        eMatRead_t = eMatRead_t,
                        number_of_non_1_reads = number_of_non_1_reads,
                        indices_of_non_1_reads = indices_of_non_1_reads
                    )
                    pA1 <- out[["pA1"]]
                    pA2 <- out[["pA2"]]                    
                } else if (currently_doing_gibbs_initialization) {
                    h_rC <- 1
                    h_rA1 <- 2
                    h_rA2 <- 3
                    pA1 <- pC
                    pA2 <- pC
                    pC[h_rC] <- 0
                    pA1[h_rA1] <- 0
                    pA2[h_rA2] <- 0
                    for(k in 1:K) {
                        ## everything is multiplied!
                        pC[h_rC] <- pC[h_rC] + alphaHat_m[h_rC, k]  * betaHat_m[h_rC, k] * eMatRead_t[k, iRead]
                        pA1[h_rA1] <- pA1[h_rA1] + alphaHat_m[h_rA1, k]  * betaHat_m[h_rA1, k] * eMatRead_t[k, iRead]
                        pA2[h_rA2] <- pA2[h_rA2] + alphaHat_m[h_rA2, k]  * betaHat_m[h_rA2, k] * eMatRead_t[k, iRead]
                    }
                } else if (currently_doing_pass_through) {
                    h_rC <- 1
                    h_rA1 <- 2
                    h_rA2 <- 3
                    ## nothing!
                    pA1 <- pC
                    pA2 <- pC
                }
                ## NOTE - THIS works out to be the same thing as e.g. prod(pA1) *  (prior_probs[h_rA1] / prior_probs[h_rC])
                ## which make a bit more sense
                prod_pC <- prod(pC) * prior_probs[h_rC]
                prod_pA1 <- prod(pA1) * prior_probs[h_rA1]
                prod_pA2 <- prod(pA2) * prior_probs[h_rA2]
                ##
                denom <- (prod_pC + prod_pA1 + prod_pA2)
                norm_pC <- prod_pC / denom
                norm_pA1 <- prod_pA1 / denom
                norm_pA2 <- prod_pA2 / denom
                if (verbose) {
                    print(paste0("normalized_probabilities = ", paste0(c(norm_pC, norm_pA1, norm_pA2), collapse = ",")))
                    print(paste0("pC = ", paste0(pC, collapse = ",")))
                    print(paste0("denom = ", denom))
                }
                ##
                ##
                ##
                ## mt = 0            -> 0.5
                ## mu = 0.5          -> 0.5 + ff / 2
                ## p  = 0.5 + ff / 2 -> 1
                ##
                ## which one is the winner
                ##
                cumsum_flip_probs <- c(0, 0, 0) ## H 1,2,3 ordering, i.e. mt, mu, p ordering
                cumsum_flip_probs[h_rC] <- norm_pC
                cumsum_flip_probs[h_rA1] <- norm_pA1
                cumsum_flip_probs[h_rA2] <- norm_pA2
                cumsum_flip_probs <- cumsum(cumsum_flip_probs)
                ##
                chance <- runif_reads[nReads * (iteration - 1) + iRead]
                h_rN <- 0
                for(i in 3:1) {
                    if (chance < cumsum_flip_probs[i]) {
                        h_rN <- i
                    }
                }
                if (!(h_rN %in% c(1, 2, 3))) {
                    print(c("h_rC", "h_rA1", "h_rA2"))
                    print(c(h_rC, h_rA1, h_rA2))
                    print(cumsum_flip_probs)
                    print(chance)
                    stop("bad h_rN")
                }
                ## if (iRead == 4) {
                ##     print(alphaHat_m[, 1:5])
                ##     print(betaHat_m[, 1:5])
                ##     print(paste0("iRead = ", iRead, ", read_category(iRead) = ", read_category[iRead]))
                ##     print(paste0("h_rC = ", h_rC, ", h_rN = ", h_rN))
                ##     print(paste0("chance = ", chance))
                ##     ##     print(paste0("cumsum_flip_probs = ", cumsum_flip_probs, collapse = ", "))
                ##     print(paste0("norm_pC = ", paste0(norm_pC, collapse = ", ")))
                ##     print(paste0("norm_pA1 = ", paste0(norm_pA1, collapse = ", ")))
                ##     print(paste0("pC = ", paste0(pC, collapse = ", ")))
                ##     print(paste0("pA1 = ", paste0(pA1, collapse = ", ")))
                ##     ##     print(paste0("prior_probs = ", paste0(prior_probs, collapse = ", ")))
                ##     ##
                ## }
                ##
                if (
                ((h_rN != h_rC) | currently_doing_gibbs_initialization) & (!currently_doing_pass_through)
                ) {
                    ##
                    H[iRead] <- h_rN
                    ## update alphas
                    if (currently_doing_normal_progress) {
                        for(k in 0:(K - 1)) {
                            alphaHat_m[h_rC, k + 1] <- alphaHat_m[h_rC, k + 1] / eMatRead_t[k + 1, iRead]
                        }
                    }
                    for(k in 0:(K - 1)) {
                        alphaHat_m[h_rN, k + 1] <- alphaHat_m[h_rN, k + 1] * eMatRead_t[k + 1, iRead]
                    }
                    ## now - first - who loses
                    if (currently_doing_normal_progress) {
                        if (h_rC == 1) { for(k in 1:K) { eMatGrid_t1[k, iGrid] = eMatGrid_t1[k, iGrid] / eMatRead_t[k, iRead] } }
                        if (h_rC == 2) { for(k in 1:K) { eMatGrid_t2[k, iGrid] = eMatGrid_t2[k, iGrid] / eMatRead_t[k, iRead] } }
                        if (h_rC == 3) { for(k in 1:K) { eMatGrid_t3[k, iGrid] = eMatGrid_t3[k, iGrid] / eMatRead_t[k, iRead] } }
                    }
                    if (h_rN == 1) { for(k in 1:K) { eMatGrid_t1[k, iGrid] = eMatGrid_t1[k, iGrid] * eMatRead_t[k, iRead] } }
                    if (h_rN == 2) { for(k in 1:K) { eMatGrid_t2[k, iGrid] = eMatGrid_t2[k, iGrid] * eMatRead_t[k, iRead] } }
                    if (h_rN == 3) { for(k in 1:K) { eMatGrid_t3[k, iGrid] = eMatGrid_t3[k, iGrid] * eMatRead_t[k, iRead] } }
                    ## reset pC
                    if (currently_doing_normal_progress) {
                        for(i in 1:3) {
                            if (h_rC == 1) {
                                if (h_rN == 2) {pC[i] = pA1[i];}
                                if (h_rN == 3) {pC[i] = pA2[i];}
                            } else if (h_rC == 2) {
                                if (h_rN == 1) {pC[i] = pA1[i];}
                                if (h_rN == 3) {pC[i] = pA2[i];}
                            } else if (h_rC == 3) {
                                if (h_rN == 1) {pC[i] = pA1[i];}
                                if (h_rN == 2) {pC[i] = pA2[i];}
                            }
                        }
                    } else if (currently_doing_gibbs_initialization) {
                        ## otherwise, if gibbs_initialize
                        ## reset no matter what, based on new sampled read
                        for(i in 1:3) {
                            ## if h_rN is 1, we're good already
                            if (h_rN == 2) {pC[i] = pA1[i];}
                            if (h_rN == 3) {pC[i] = pA2[i];}
                        }
                    }
                    ##if (iRead <= 4) {
                    ##    print(paste0("the new PC is ", paste0(pC, collapse = ", ")))
                    ## }
                }
                ##
                ##
                ##
                if ((sum(pC > 1*10**50) > 0) | (sum(pC < 1*10**(-50)) > 0)) {
                    if (verbose) {
                        print(paste0("renormalize"))
                    }
                    ## inject back
                    alphaHat_t1[, iGrid] <- alphaHat_m[1, ]
                    alphaHat_t2[, iGrid] <- alphaHat_m[2, ]
                    alphaHat_t3[, iGrid] <- alphaHat_m[3, ]
                    ## update-c
                    a <- 1 / sum(alphaHat_t1[, iGrid])
                    c1[iGrid] <- c1[iGrid] * a
                    alphaHat_t1[, iGrid] <- alphaHat_t1[, iGrid] * a
                    print(paste0("a = ", a))
                    ##
                    a <- 1 / sum(alphaHat_t2[, iGrid])
                    c2[iGrid] <- c2[iGrid] * a
                    alphaHat_t2[, iGrid] <- alphaHat_t2[, iGrid] * a
                    ##
                    a <- 1 / sum(alphaHat_t3[, iGrid])
                    c3[iGrid] <- c3[iGrid] * a
                    alphaHat_t3[, iGrid] <- alphaHat_t3[, iGrid] * a
                    ## now re-get stuff
                    alphaHat_m <- rbind(alphaHat_t1[, iGrid], alphaHat_t2[, iGrid], alphaHat_t3[, iGrid])
                    betaHat_m <- rbind(betaHat_t1[, iGrid], betaHat_t2[, iGrid], betaHat_t3[, iGrid])
                    print(paste0("before probabilities = ", paste0(pC, collapse = ", ")))                
                    pC <- rowSums(alphaHat_m * betaHat_m) ##
                    print(paste0("renormalized_probabilities = ", paste0(pC, collapse = ", ")))
                    print(paste0("pC = ", paste0(pC, collapse = ",")))
                }
                ##
                ## record the read set it came from, if so inclined
                ##
                if (record_read_set) {
                    ## determine how close it is to the three probs
                    x <- array(NA, 3)
                    x[h_rC] <- norm_pC
                    x[h_rA1] <- norm_pA1
                    x[h_rA2] <- norm_pA2
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
                        H_class[iRead] <- which.min(y)
                    } else {
                        H_class[iRead] <- 0
                    }
                }
                ##
                if (return_p_store) {
                    w <-
                        (i_gibbs_samplings - 1) * n_gibbs_full_its * nReads +
                        (iteration - 1) * nReads + iRead
                    ##p_store[w, c("p_1", "p_2", "p_3")[h_rC]] <- norm_pC
                    ##p_store[w, c("p_1", "p_2", "p_3")[h_rA1]] <- norm_pA1
                    ##p_store[w, c("p_1", "p_2", "p_3")[h_rA2]] <- norm_pA2
                    p_store[w, h_rC] <- norm_pC
                    p_store[w, h_rA1] <- norm_pA1
                    p_store[w, h_rA2] <- norm_pA2
                    p_store[w, "chance"] <- chance
                    p_store[w, "h_rC"] <- h_rC
                    p_store[w, "h_rN"] <- h_rN
                    p_store[w, "agreePer"] <- sum(abs(H ==  true_H)) / length(H) * 100
                    p_store[w, "itType"] <- sum(as.integer((0:2) * c(currently_doing_pass_through, currently_doing_gibbs_initialization, currently_doing_normal_progress)))
                    ##
                    p_store[w, "p_O1_given_H1_L"] <- (-sum(log(c1)) + log(sum(alphaHat_m[1, ])))
                    p_store[w, "p_O2_given_H2_L"] <- (-sum(log(c2)) + log(sum(alphaHat_m[2, ])))
                    p_store[w, "p_O3_given_H3_L"] <- (-sum(log(c3)) + log(sum(alphaHat_m[3, ])))
                    p_store[w, "p_O_given_H_L"] <- p_store[w, "p_O1_given_H1_L"] + p_store[w, "p_O2_given_H2_L"] + p_store[w, "p_O3_given_H3_L"]
                    p_store[w, "p_H_given_L"] <- sum(log(prior_probs[H]))
                    p_store[w, "p_H_given_O_L_up_to_C"] <- p_store[w, "p_O_given_H_L"] + p_store[w, "p_H_given_L"]
                }
            }
            ##
            ## do this for all reads anyway
            ##
            iRead <- iRead + 1
            if (nReads < iRead) {
                done_reads <- TRUE
                read_wif_iRead <- -1;
            } else {
                read_wif_iRead <- sampleReads[[iRead]][[2]];
            }
        }
        iRead <- iRead - 1
        iReadEnd <- iRead
        if (verbose) {
            print(paste0("INJECT BACK"))
        }
        ## inject back into alphaHat_t1
        alphaHat_t1[, iGrid] <- alphaHat_m[1, ]
        alphaHat_t2[, iGrid] <- alphaHat_m[2, ]
        alphaHat_t3[, iGrid] <- alphaHat_m[3, ]
        ## update-c
        a <- 1 / sum(alphaHat_t1[, iGrid])
        c1[iGrid] <- c1[iGrid] * a
        alphaHat_t1[, iGrid] <- alphaHat_t1[, iGrid] * a
        ##
        a <- 1 / sum(alphaHat_t2[, iGrid])
        c2[iGrid] <- c2[iGrid] * a
        alphaHat_t2[, iGrid] <- alphaHat_t2[, iGrid] * a
        ##
        a <- 1 / sum(alphaHat_t3[, iGrid])
        c3[iGrid] <- c3[iGrid] * a
        alphaHat_t3[, iGrid] <- alphaHat_t3[, iGrid] * a
    }
    if (verbose) {
        print("FINALIZE ITERATION AND RESET")
    }
    betaHat_t1[, iGrid] <- c1[nGrids]
    betaHat_t2[, iGrid] <- c2[nGrids]
    betaHat_t3[, iGrid] <- c3[nGrids]
    Rcpp_run_backward_haploid(
        betaHat_t1,
        c = c1,
        eMatGrid_t = eMatGrid_t1,
        alphaMatCurrent_t = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        s = s - 1
    )
    Rcpp_run_backward_haploid(
        betaHat_t2,
        c = c2,
        eMatGrid_t = eMatGrid_t2,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        s = s - 1
    )
    Rcpp_run_backward_haploid(
        betaHat_t3,
        c = c3,
        eMatGrid_t = eMatGrid_t3,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        s = s - 1
    )
    ## try complete re-labelling here
    relabel <- 1
    if (do_block_resampling) {
        out <- consider_and_try_entire_relabelling(H = H, ff = ff, relabel = artificial_relabel)
        relabel <- out$relabel
        H <- out$H
        if (relabel > 1) {
            ## alpha
            out <- apply_relabel(alphaHat_t1, alphaHat_t2, alphaHat_t3, relabel)
            alphaHat_t1 <- out$m1; alphaHat_t2 <- out$m2; alphaHat_t3 <- out$m3
            ## beta
            out <- apply_relabel(betaHat_t1, betaHat_t2, betaHat_t3, relabel)
            betaHat_t1 <- out$m1; betaHat_t2 <- out$m2; betaHat_t3 <- out$m3
            ## eHatMapSNP
            out <- apply_relabel(eMatGrid_t1, eMatGrid_t2, eMatGrid_t3, relabel)
            eMatGrid_t1 <- out$m1; eMatGrid_t2 <- out$m2; eMatGrid_t3 <- out$m3
            ## c
            out <- apply_relabel(c1, c2, c3, relabel)
            c1 <- out$m1; c2 <- out$m2; c3 <- out$m3
        }
    }
    ##
    per_it_likelihoods <- add_to_per_it_likelihoods(s = s, per_it_likelihoods = per_it_likelihoods, i_gibbs_samplings = i_gibbs_samplings, iteration = iteration, i_result_it = i_result_it, n_gibbs_full_its =  n_gibbs_full_its, H = H, c1 = c1, c2 = c2, c3 = c3, ff = ff, prior_probs = prior_probs, relabel = relabel)
    return(
        list(
            alphaHat_t1 = alphaHat_t1,
            betaHat_t1 = betaHat_t1,
            c1 = c1,
            eMatGrid_t1 = eMatGrid_t1,
            alphaHat_t2 = alphaHat_t2,
            betaHat_t2 = betaHat_t2,
            c2 = c2,
            eMatGrid_t2 = eMatGrid_t2,
            alphaHat_t3 = alphaHat_t3,
            betaHat_t3 = betaHat_t3,
            c3 = c3,
            eMatGrid_t3 = eMatGrid_t3,
            p_store = p_store,
            H = H,
            H_class = H_class,
            per_it_likelihoods = per_it_likelihoods,
            currently_doing_normal_progress = currently_doing_normal_progress,
            currently_doing_gibbs_initialization = currently_doing_gibbs_initialization,
            currently_doing_pass_through = currently_doing_pass_through
        )
    )
}


determine_label_probabilities <- function(rc, ff, return_neutral = TRUE) {
    if (return_neutral) {
        prob_matrix <- array(0, c(3, 3))
        diag(prob_matrix) <- 1
        return(prob_matrix)
    }
    p <- c(0.5, (1 - ff) / 2, (ff / 2))
    log_p <- log(p)
    options <- array(0, 6)
    rc_mat <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(2, 3, 1),
        c(3, 1, 2),
        c(3, 2, 1)
    )
    for(i in 1:6) {
        options[i] <-
            rc[rc_mat[i, 1]] * log_p[1] +
            rc[rc_mat[i, 2]] * log_p[2] +
            rc[rc_mat[i, 3]] * log_p[3]
        if (i == 1) {
            max <- options[1]
        } else {
            if (sum(is.na(options)) > 0) {
                print(rc)
            }
            if (options[i] > max) {
                max <- options[i]
            }
        }
    }
    probs <- array(0, 6)
    prob_sum <- 0
    for(i in 1:6) {
        options[i] <- options[i] - max
        if (options[i] < (-100)) {
            options[i] <- (-100)
        }
        probs[i] <- exp(options[i])
        prob_sum <- prob_sum + probs[i]
    }
    for(i in 1:6) {
        probs[i] <- probs[i] / prob_sum
    }
    ## now - make matrix from this
    ## matrix is probability of that original labelling (row) being that option (column)
    prob_matrix <- array(0, c(3, 3))
    for(original_label in 1:3) {
        for(true_label in 1:3) {
            for(i in 1:6) {
                if (rc_mat[i, true_label] == original_label) {
                    prob_matrix[original_label, true_label] <-
                        prob_matrix[original_label, true_label] + probs[i]
                }
            }
        }
    }
    ## so e.g. prob_matrix[label is 1/2/3, is really mt/mu/p]
    return(prob_matrix)
}








## equal weighting - if TRUE, then ignore ll, everything just averaged
weighted_average_on_the_fly <- function(g, ll, log_mult_max = 40, equal_weighting = FALSE, language = "R") {
    if (language == "R") {
        local_fly_weighter <- fly_weighter
    } else if (language == "Rcpp") {
        local_fly_weighter <- rcpp_fly_weighter
    } else {
        stop("bad language")
    }
    ##
    ll_on_fly <- array(0, length(ll))
    ll_rescaled <- array(100, length(ll))
    ##
    for(i in 1:length(g)) {
        ll_current <- ll[i]
        ll_on_fly[i] <- ll_current
        gCurrent <- g[[i]]
        if (i == 1) {
            gAverage <- array(0, dim(gCurrent))
            log_mult <- array(-1, 2) ## for fucks sake R
        }
        ##
        out <- local_fly_weighter(
            i = i,
            gCurrent = gCurrent,
            gAverage = gAverage, ## for i=1, reset
            ll_current = ll_current,
            log_mult = log_mult,
            ll_rescaled = ll_rescaled,
            gCurrent2 = array(0, c(1, 1)), ## do not care
            gAverage2 = array(0, c(1, 1)),
            gCurrent3 = array(0, c(1, 1)), ## do not care
            gAverage3 = array(0, c(1, 1)),
            relative_difference = -1,
            log_rescale = -1,
            log_mult_max = log_mult_max,
            equal_weighting = equal_weighting
        )
        if (language == "R") {
            gAverage <- out$gAverage
            log_mult <- out$log_mult
            ll_rescaled <- out$ll_rescaled
        }
        ##
    }
    gAverage <- gAverage / sum(exp(ll_rescaled))
    return(gAverage)
}

## here, basically, the above is a skeleton
fly_weighter <- function(
    i,
    gCurrent,
    gAverage,
    ll_current,
    log_mult,
    ll_rescaled,
    gCurrent2,
    gAverage2,
    gCurrent3,
    gAverage3,
    relative_difference = -1,
    log_rescale = -1,
    log_mult_max = 40,
    equal_weighting = FALSE
) {
    if (equal_weighting) {
        ll_current <- 1
    }
    log_rescale <- NA
    relative_difference <- 1
    if (i == 1) {
        log_mult <- ll_current
        gAverage <- 1 * gCurrent ## running average
        gAverage2 <- 1 * gCurrent2 ## running average
        gAverage3 <- 1 * gCurrent3 ## running average
    } else {
        relative_log_difference <- ll_current - log_mult
        if (relative_log_difference > log_mult_max) {
            ##
            log_rescale <- log_mult - ll_current
            ll_rescaled[1:(i - 1)] <- ll_rescaled[1:(i - 1)] + log_rescale
            ## rescale
            gAverage <- gAverage * exp(log_rescale)
            gAverage2 <- gAverage2 * exp(log_rescale)
            gAverage3 <- gAverage3 * exp(log_rescale)            
            ## this is the new value
            log_mult <- ll_current
        }
        relative_difference <- exp(ll_current - log_mult)
        gAverage <- gAverage + relative_difference * gCurrent
        gAverage2 <- gAverage2 + relative_difference * gCurrent2
        gAverage3 <- gAverage3 + relative_difference * gCurrent3
    }
    ll_rescaled[i] <- ll_current - log_mult
    ## if (verbose) {
    ##     ## not
    ##     a <- head(colSums(weights[1:i] * t(sapply(g[1:i], I))) / sum(weights[1:i]), 1)
    ##     b <- gAverage[1] / sum(exp(ll_rescaled[1:i]))
    ##     print_message(paste0("truth vs estimate = ", paste0(c(a, b, a - b), collapse = ", ")))
    ## }
    return(
        list(
            gAverage = gAverage,
            gAverage2 = gAverage2,
            gAverage3 = gAverage3,
            log_mult = log_mult,
            ll_rescaled = ll_rescaled,
            log_rescale = log_rescale,
            relative_difference = relative_difference
        )
    )
}


single_fly_weight <- function(
    i,
    relative_difference,
    gCurrent,
    gAverage,
    log_rescale = NA    
) {
    if (i == 1) {
        gAverage <- gCurrent
    } else {
        if (!is.na(log_rescale)) {
            gAverage <- gAverage * exp(log_rescale)        
        }
        gAverage <- gAverage + relative_difference * gCurrent
    }
    return(gAverage)
}

## not fully polished
plot_label_usage_using_rect <- function(
    fbsoL,
    log_l,
    ff,
    order = NULL,
    true_H = NULL,
    true_H_ll = NULL,
    sampleReads = NULL,
    gibbs_first_read = 0,
    relabels = NULL
) {
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    z <- t(sapply(fbsoL[[1]]$list_of_ending_read_labels, I))
    labels <- log_l
    if (is.null(order)) {
        order <- 1:nrow(z)
    }
    z <- z[order, ]
    labels <- labels[order]
    ##
    x <- 1:(ncol(z) + 1)
    y <- 1:(nrow(z) + 1)
    xlim <- c(-5, max(x))
    ylim <- c(1, max(y) + 2)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, col = "white", axes = FALSE)
    names(cbPalette)[1:3] <- c("mt", "mu", "p")
    if (!is.null(sampleReads)) {
        wif <- sapply(sampleReads, function(x) x[[2]])
        a <- unique(wif)
        for(aa in a) {
            w <- which(wif == aa)
            abline(v = head(w, 1))
            abline(v = tail(w, 1) + 1)
        }
    }
    for(i_row in 1:(nrow(z))) {
        ## lol... anyway
        H <- z[i_row, ]
        ## simplest thing - choose "mt" as label that's best for mt
        hlpm <- determine_label_probabilities(table(H), ff)
        ## argh -
        best <- c(which.max(hlpm[, 1]), mu = which.max(hlpm[, 2]), p = which.max(hlpm[, 3]))
        best2 <- c(mt = which.max(hlpm[1, ]), mu = which.max(hlpm[2, ]), p = which.max(hlpm[3, ]))
        H2 <- c("mt", "mu", "p")[match(H, best)]
        H2 <- c("mt", "mu", "p")[H]
        col <- cbPalette[H2]
        rect(
            xleft = x[-length(x)], xright = x[-1],
            ybottom = i_row,
            ytop = i_row + 1,
            col = col,
            border = NA
        )
        if ((i_row == 1) & (gibbs_first_read > 0)) {
            rect(
                xleft = gibbs_first_read, xright = gibbs_first_read + 1,
                ybottom = i_row,
                ytop = i_row + 1,
                col = NA,
                border = "red"
            )
        }
        ##text(x = 20 + -1:1, y = 0.5 + i_row, label = best)
        ##text(x = 25 + -1:1, y = 0.5 + i_row, label = best2)
        ## if (i_row == 50) {
        ##     rect(xleft = x[1], xright = tail(x, 1),
        ##          ybottom = i_row,
        ##          ytop = i_row + 1,
        ##          border = "red")
        ## }
    }
    if (!is.null(true_H)) {
        i_row <- i_row + 1.5
        H <- true_H
        col <- cbPalette[H]
        rect(
            xleft = x[-length(x)], xright = x[-1],
            ybottom = i_row,
            ytop = i_row + 1,
            col = col,
            border = NA
        )
        rect(
            xleft = x[1], xright = tail(x, 1),
            ybottom = i_row,
            ytop = i_row + 1,
            border = "red"
        )
    }
    q <- diff(xlim) * 0.02
    y <- 0.5 + 1:length(labels)
    if (!is.null(true_H_ll)) {
        labels <- c(labels, true_H_ll) ## default is NULL
        y <- c(y, tail(y, 1) + 1.5)
    }
    ## add some labelling!
    text(x = xlim[1] + q, y = y, labels = round(labels, 2), cex = 1)
    y <- 0.5 + 1:nrow(z)
    y <- c(y, tail(y, 1) + diff(tail(y, 2)))
    if (!is.null(relabels)) {
        text(x = xlim[2] - 2 * q, y = y, labels = c(relabels, "relabels"), cex = 1)
    }
    text(x = xlim[2] - q, y = y, labels = c(order, "order"), cex = 1)
    ##
}

calc_prob_of_set_of_reads <- function(ff, rc, reorder = NULL) {
    if (!is.null(reorder)) {
        rc <- rc[reorder]
    }
    return(dmultinom(rc, prob = c(0.5, 0.5 - ff / 2, ff / 2), log = TRUE))
    ## n <- sum(rc)
    ## log_p <-
    ##     n * log(n) - n +
    ##     - (rc[1] * log(rc[1]) - rc[1]) +
    ##     - (rc[2] * log(rc[2]) - rc[2]) +
    ##     - (rc[3] * log(rc[3]) - rc[3]) +
    ##     rc[1] * log(0.5) +
    ##     rc[2] * log(0.5 - ff / 2) +
    ##     rc[3] * log(ff / 2)
    ## return(log_p)
}


## consider collapsed options
consider_collapsed_options <- function(a1, a2, b, rc, ff) {
    a = c(a1, a2);
    x <- c(sum(rc[a]), rc[b])
    prob <- c(0.5, 0.5 - ff / 2, ff / 2)
    p <- c(sum(prob[a]), prob[b])
    dmultinom(x = x, prob = p, log = TRUE)
}


get_weights_for_entire_relabelling <- function(rc, ff) {
    ## 1, 2, 3, / 1, 3, 2, etc
    reorder_log_probs <- c(
        calc_prob_of_set_of_reads(ff, rc = rc, reorder = c(1, 2, 3)),
        calc_prob_of_set_of_reads(ff, rc = rc, reorder = c(1, 3, 2)),
        calc_prob_of_set_of_reads(ff, rc = rc, reorder = c(2, 1, 3)),
        calc_prob_of_set_of_reads(ff, rc = rc, reorder = c(2, 3, 1)),
        calc_prob_of_set_of_reads(ff, rc = rc, reorder = c(3, 1, 2)),
        calc_prob_of_set_of_reads(ff, rc = rc, reorder = c(3, 2, 1))
    )
    max <- max(reorder_log_probs)
    weights <- exp(reorder_log_probs - max)
    weights <- weights / sum(weights)
    return(weights)
}


## so here we do a gibbs sample on a set of read labels
## in the first case, an entire relabelling
## i.e. see if a permutation of the current labels would be beneficial
consider_and_try_entire_relabelling <- function(H, ff, take_max = FALSE, relabel = -1) {
    ##
    rc <- c(sum(H == 1), sum(H == 2), sum(H == 3))
    ## so I can do this easily at the end of every iteration, will save some I think
    weights <- get_weights_for_entire_relabelling(rc, ff)
    ## "uncollapsing" much harder!
    if (relabel == -1) {
        if (!take_max) {
            relabel <- sample(1:6, prob = weights, size = 1)
        } else {
            relabel <- which.max(weights)
        }
    }
    if (relabel != 1) {
        ## make a note of what's been done
        if (relabel == 1) {reorderX <- c(1, 2, 3)}
        if (relabel == 2) {reorderX <- c(1, 3, 2)}
        if (relabel == 3) {reorderX <- c(2, 1, 3)}
        if (relabel == 4) {reorderX <- c(3, 1, 2)}
        if (relabel == 5) {reorderX <- c(2, 3, 1)}
        if (relabel == 6) {reorderX <- c(3, 2, 1)}
        ## apply weighting
        new_H <- H
        for(iRead in 1:length(H)) {
            new_H[iRead] <- reorderX[H[iRead]]
        }
        H <- new_H
    }
    return(list(H = H, relabel = relabel))
    ## would an "uncollapsing" help things
    ## uncollapsed_log_probs <- c(
    ##     consider_collapsed_options(1, 2, 3, rc, ff), ## all that really matters is third one
    ##     consider_collapsed_options(1, 3, 2, rc, ff), ## what we have here
    ##     consider_collapsed_options(2, 3, 1, rc, ff)
    ## )
}



apply_relabel <- function(m1, m2, m3, relabel) {
    if (relabel == 2) {
        ## reorder <- c(1, 3, 2)
        ## swap 2 and 3
        mT <- m2
        m2 <- m3
        m3 <- mT
    }
    if (relabel == 3) {
        ## reorder <- c(2, 1, 3)
        ## swap 2 and 1
        mT <- m2
        m2 <- m1
        m1 <- mT
    }
    if (relabel == 4) {
        ## reorderX <- c(3, 1, 2)
        mT <- m1
        m1 <- m2
        m2 <- m3
        m3 <- mT
    }
    if (relabel == 5) {
        ## reorderX <- c(2, 3, 1)
        mT <- m1
        m1 <- m3
        m3 <- m2
        m2 <- mT
    }
    if (relabel == 6) {
        ## reorder <- c(3, 2, 1)
        ## swap 3 and 1
        mT <- m3
        m3 <- m1
        m1 <- mT
    }
    return(
        list(
            m1 = m1,
            m2 = m2,
            m3 = m3
        )
    )
}


add_to_per_it_likelihoods <- function(
    s,
    per_it_likelihoods,
    i_gibbs_samplings,
    iteration,
    i_result_it,
    n_gibbs_full_its,
    H,
    c1,
    c2,
    c3,
    ff,
    prior_probs,
    relabel
) {
    i <- (i_gibbs_samplings - 1) * n_gibbs_full_its + iteration ## yup! (i_outer?)
    per_it_likelihoods[i, "s"] <- s
    per_it_likelihoods[i, "i_samp"] <- i_gibbs_samplings
    per_it_likelihoods[i, "i_it"] <- iteration
    per_it_likelihoods[i, "i_result_it"] <- i_result_it
    per_it_likelihoods[i, "p_O1_given_H1_L"] <- -sum(log(c1))
    per_it_likelihoods[i, "p_O2_given_H2_L"] <- -sum(log(c2))
    per_it_likelihoods[i, "p_O3_given_H3_L"] <- -sum(log(c3))
    per_it_likelihoods[i, "p_O_given_H_L"] <- sum(per_it_likelihoods[i, c("p_O1_given_H1_L", "p_O2_given_H2_L", "p_O3_given_H3_L")])
    per_it_likelihoods[i, "p_H_given_L"] <- sum(log(prior_probs[H]))    
    per_it_likelihoods[i, "p_H_given_O_L_up_to_C"] <- per_it_likelihoods[i, "p_O_given_H_L"] + per_it_likelihoods[i, "p_H_given_L"]
    rc <- c(sum(H == 1), sum(H == 2), sum(H == 3))
    per_it_likelihoods[i, "p_set_H_given_L"] <- calc_prob_of_set_of_reads(ff = ff, rc = rc)
    per_it_likelihoods[i, "relabel"] <- relabel
    return(per_it_likelihoods)
}


## what do I want to see out of this afterwards - log likehoods? accuracy per iteration / process?
get_temporary_experimental_variables <- function(temporary_experimental_var, i_loop) {
    if (temporary_experimental_var == 1) {
        ## always - fixed, coming in, n_gibbs_starts = 10 (or 1?)
        n_loop <- 3
        method_list <- c("triploid-nipt", "gibbs-nipt", "gibbs-nipt")
        K_list <- c(10, 100, 1700)
        n_gibbs_full_its_list <- c(NA, 10, 0)
        gridWindowSize_list <- c(10000, 10000, NA)
        nGen_list <- round(4 * 20000 / K_list)
    }
    return(
        list(
            n_loop = n_loop,
            K = K_list[i_loop],
            method = method_list[i_loop],
            n_gibbs_full_its = n_gibbs_full_its_list[i_loop],
            gridWindowSize = gridWindowSize_list[i_loop],
            nGen = nGen_list[i_loop]
        )
    )
}



calculate_likelihoods_consistent_with_the_truth <- function(
    phase,
    samples_with_phase,
    tempdir,
    regionName,
    priorCurrent_m,
    alphaMatCurrent_tc,
    eHapsCurrent_tc,
    transMatRate_tc_H,
    maxDifferenceBetweenReads,
    maxEmissionMatrixDifference,
    Jmax,
    ffs,
    bundling_info,
    grid,
    nCores,
    nLLs = 10,
    verbose = FALSE
) {
    ## this must exist and be 
    n_phase <- length(samples_with_phase)
    if (n_phase == 0) {
        return(NULL) ## should not be called though
    }
    print_message("Begin calculating truth likelihoods and labels given phase")    
    lls <- array(0, c(n_phase + 1, nLLs))
    ## iiSample is 1, 2, 3, ... among samples with phase
    ## iSample is there index among full set from sampleNames / canonical ordering
    iiSample <- 1
    S <- dim(eHapsCurrent_tc)[3]
    ##
    truth_likelihoods_and_labels <- mclapply(
        1:n_phase,
        mc.cores = nCores,
        FUN = function(iiSample) {
            nSNPs <- dim(phase)[[1]]
            truth_haps_t <- array(NA, c(3, nSNPs))
            iSample <- samples_with_phase[iiSample]
            if (verbose) {
                print_message(paste0("iiSample = ", iiSample))
                print_message(paste0("iSample = ", iSample))
            }
            ff <- ffs[iSample]
            ##
            if (verbose) {            
                print_message("load sample reads")
            }
            sampleReads <- get_sampleReads_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = iSample,
                bundling_info = bundling_info
            )$sampleReads
            ## assign those reads to haplotypes with probabilities
            for(i in 1:3) {
                truth_haps_t[i, ] <- phase[, iiSample, i, drop = FALSE]
            }
            readProbs_t <- assign_fetal_read_probabilities(
                sampleReads = sampleReads,
                truth_haps_t = truth_haps_t
            )
            ## convert those probabilities into probabilities they originate from sets of haplotypes
            ## e.g. given ff, and a read could come from matt and p, what is the probability it came from each
            if (verbose) {            
                print_message("get read groupings")
            }
            out <- get_read_groupings_given_fetal_fraction_and_cov(
                readProbs_t = readProbs_t,
                phase = phase,
                iiSample = iiSample,
                ff = ff
            )
            ## sample read labels consistent with the truth            
            if (verbose) {            
                print_message("sample read labels")
            }
            double_list_of_sampled_truth_read_labels <- lapply(1:S, function(s) {
                lapply(1:nLLs, function(iLL) {
                    H <- sample_H_for_NIPT_given_groupings(
                        groupings = out$groupings,
                        counts = out$counts,
                        ff = ff
                    )
                    return(H)
                })
            })
            ##
            if (verbose) {            
                print_message("calculate likelihoods")
            }
            ##
            out <- run_forward_backwards(
                sampleReads = sampleReads,
                method = "gibbs-nipt",
                priorCurrent_m = priorCurrent_m,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                eHapsCurrent_tc = eHapsCurrent_tc,
                transMatRate_tc_H = transMatRate_tc_H,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                maxEmissionMatrixDifference = maxEmissionMatrixDifference,
                Jmax = Jmax,
                n_gibbs_starts = nLLs, ## might be longer!
                n_gibbs_burn_in_its = 0,
                n_gibbs_sample_its = 0,
                double_list_of_starting_read_labels = double_list_of_sampled_truth_read_labels,
                ffs = ffs,
                iSample = iSample,
                suppressOutput = 1,
                return_gamma = FALSE,
                return_genProbs = FALSE,
                return_hapProbs = FALSE,
                grid = grid
            )
            ## extract likelihoods!
            return(
                list(
                    likelihoods = out[[1]]$per_it_likelihoods,
                    double_list_of_sampled_truth_read_labels = double_list_of_sampled_truth_read_labels
                )
            )
        }
    )
    check_mclapply_OK(truth_likelihoods_and_labels)
    print_message("Done calculating truth likelihoods and labels given phase")
    return(truth_likelihoods_and_labels)
}

assign_fetal_read_probabilities <- function(
    sampleReads,
    truth_haps_t
) {
    ## get phased haplotypes
    eMatRead_t <- array(1, c(nrow(truth_haps_t), length(sampleReads)))
    eHapsCurrent_tc <- array(NA, c(dim(truth_haps_t), 1))
    eHapsCurrent_tc[, , 1] <- truth_haps_t
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = eHapsCurrent_tc,
        s = 1 - 1,
        maxDifferenceBetweenReads = 1000000,
        Jmax = 100,
        eMatHapOri_t = array(0, c(1, 1)),
        pRgivenH1 = array(1),
        pRgivenH2 = array(1),
        prev = 1,
        suppressOutput = 1,
        prev_section = "text",
        next_section = "text",
        run_pseudo_haploid = FALSE,
        rescale_eMatRead_t = FALSE
    )
    return(eMatRead_t) ## is readProbs_t
}


get_read_groupings_given_fetal_fraction_and_cov <- function(
    readProbs_t,
    phase,
    iiSample,
    ff
) {
    ## should not be "real" but OK
    real_frp <- c(0.5, 0.5 - ff / 2, ff / 2)
    names(real_frp) <- c("matt", "matu", "pat")
    n_sample_reads <- ncol(readProbs_t)
    ## round - but preserve sum
    ## n_reads_to_sample <- preserve_round(n_sample_reads * real_frp)
    ## I think this is easier than I made it out to be
    matt <- readProbs_t[1, ] > 0.5
    matu <- readProbs_t[2, ] > 0.5
    pat <- readProbs_t[3, ] > 0.5
    matt[is.na(matt)] <- FALSE
    matu[is.na(matu)] <- FALSE
    pat[is.na(pat)] <- FALSE        
    ## 8 options
    ## order matt, matu, pat - first by first entry
    groupings <- cbind(
        all = (matt & matu & pat) | (!matt & !matu & !pat),
        matt_matu = matt & matu & !pat,
        matt_pat = matt & !matu & pat,
        matu_pat = !matt & matu & pat,
        matt_only = matt & !matu & !pat,
        matu_only = !matt & matu & !pat,           
        pat_only = !matt & !matu & pat
    )
    ## now, want a single, large, desired fraction of total
    ##
    ## figure out proportions, for mat1, mat2 and pat, then draw!
    ##
    groups <- c("matt", "matu", "pat")
    ## 
    ## OK - can be things split between two groups that have integer number of counts
    ## need to sample if too close?
    counts <- lapply(groups, function(who) {
        ##
        other1 <- setdiff(groups, who)[1]
        other2 <- setdiff(groups, who)[2]
        counts <- array(0, 4)
        names(counts) <- c(who, other1, other2, "all")
        return(counts)
    })
    names(counts) <- groups
    cg <- colSums(groupings)
    ## do the joint one - again, preserve sum
    all_split <- preserve_round(real_frp * cg["all"])
    ## do the single ones - easy
    for(who in groups) {
        counts[[who]][who] <- cg[paste0(who, "_only")]
        counts[[who]]["all"] <- all_split[who]
        ## do vs other here!
        ## need AGAIN to use preserve round
        for(j in 1:2) {
            other <- setdiff(groups, who)[j]
            f1 <- real_frp[who] / (real_frp[who] + real_frp[other])
            a <- paste0(sort(c(who, other)), collapse = "_")
            to_use <- preserve_round(cg[a] * c(f1, 1 - f1))
            ## do not really care about over-writing - will over-write both!
            counts[[who]][other] <- to_use[1]
            counts[[other]][who] <- to_use[2]
        }
    }
    ## check - this is fiddly - throw errow NOW if wrong
    check_counts_groupings(counts, groupings)
    ## also - check each entry
    n_reads_to_sample <- sapply(counts, sum)
    ##
    ## now, sample reads from these!
    ##
    reads <- sort(unlist(sapply(groups, function(who) {
        ##
        other1 <- setdiff(groups, who)[1]
        other2 <- setdiff(groups, who)[2]
        index_who_only <- groupings[, paste0(who, "_only")]
        index_who_other1 <- groupings[,  paste0(sort(c(who, other1)), collapse = "_")]
        index_who_other2 <- groupings[,  paste0(sort(c(who, other2)), collapse = "_")]
        index_all <- groupings[, "all"]        
        n_iwo <- sum(index_who_only)
        n_iwo1 <- sum(index_who_other1)
        n_iwo2 <- sum(index_who_other2)        
        n_all <- sum(index_all)
        ## 
        local_counts <- counts[[who]]
        if (sum(local_counts) == 0) {
            reads <- NULL
        } else {
            n_to_sample <- local_counts
            ## sample reads here
            reads <- sort(c(
                sample(which(index_who_only), min(n_to_sample[1], n_iwo)),
                sample(which(index_who_other1), min(n_to_sample[2], n_iwo1)),
                sample(which(index_who_other2), min(n_to_sample[3], n_iwo2)),
                sample(which(index_all), n_to_sample[4])
            ))
        }
        return(reads)
    })))
    ##
    if (length(reads) != sum(n_reads_to_sample)) {
        print(length(reads))
        print(sum(n_reads_to_sample))
        stop("Did not select reads properly")
    }
    ##
    return(
        list(
            reads = reads,
            groupings = groupings,
            counts = counts
        )
    )
}

check_counts_groupings <- function(counts, groupings) {
    a <- sum(unlist(counts))
    b <- nrow(groupings)
    if (a != b) {
        print(counts)
        print(paste0("counts has a sum of ", a))
        print(paste0("groupings is ", b, " long"))
        stop("Count association and groupings does not agree")
    }
    return(NULL)
}
    

## 
preserve_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    if (abs(sum(x) - sum(y)) > 0.1) {
        print(x)
        print(y)
        stop("preserve round has not worked")
    }
    return(y)
}



## two different things
## downsample and get new set of reads
## choose read labels given readProbs_t for all reads
sample_H_for_NIPT_given_groupings <- function(
    groupings,
    counts,
    ff
)  {
    ## OK - here we want labels
    nReads <- nrow(groupings)
    H <- array(0, nReads)
    groups <- names(counts)
    ## all - sample according to probabilities
    ## but H is in units of 1 = mt, 2 = mu, 3 = p
    prob <- c(0.5, 0.5 - ff / 2, ff / 2)
    names(prob) <- c("matt", "matu", "pt")
    ## 
    H[which(groupings[, "all"])] <- sample(
        1:3,
        size = sum(groupings[, "all"]),
        prob = prob,
        replace = TRUE
    )
    ## note - above is all = yes or all = no, so both are fine
    ## onlies
    H[which(groupings[, "matt_only"])] <- 1
    H[which(groupings[, "matu_only"])] <- 2
    H[which(groupings[, "pat_only"])] <- 3
    ## shared
    for(i in 1:2) {
        for(j in (i + 1):3) {
            who1 <- groups[i]
            who2 <- groups[j]
            a <- counts[[who1]][[who2]]
            b <- counts[[who2]][[who1]]
            if ((a + b) > 0) {
                prob <- c(a / (a + b), b / (a + b))
                if (sum(is.na(prob)) > 0) {
                    print("bad prob 2")                    
                    print(a)
                    print(b)
                    print(prob)
                }
                ## 
                H[which(groupings[, paste0(who1, "_", who2)])] <- sample(
                    c(i, j), size = a + b, prob = prob, replace = TRUE
                )
            }
        }
    }
    ##
    return(H)
}






initialize_gibbs_forward_backward <- function(
    H,
    hap,
    s,
    sampleReads,
    priorCurrent_m,
    alphaMatCurrent_tc,
    transMatRate_tc_H,
    eMatGrid_t,
    eMatRead_t = NULL,
    bound_eMatGrid_t = FALSE,
    rescale_eMatGrid_t = FALSE,
    maxEmissionMatrixDifference = 1000
) {
    nReads <- length(sampleReads)
    K <- dim(alphaMatCurrent_tc)[1]
    nGrids <- dim(alphaMatCurrent_tc)[2] + 1
    ## initialize
    c <- array(0, nGrids)
    alphaHat_t <- array(0, c(K, nGrids))
    betaHat_t <- array(0, c(K, nGrids))        
    ##
    if (is.null(eMatGrid_t)) {
        eMatGrid_t <- array(1, c(K, nReads))
        rcpp_make_eMatGrid_t(
            eMatGrid_t = eMatGrid_t,
            eMatRead_t = eMatRead_t,
            H = H,
            sampleReads = sampleReads,
            hap = hap,
            nGrids = nGrids,
            prev = 0,
            suppressOutput = 1,
            prev_section = "",
            next_section = "",
            use_all_reads = FALSE,
            bound_eMatGrid_t = bound_eMatGrid_t,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            rescale_eMatGrid_t = rescale_eMatGrid_t
        )
    }
    ## 
    ## alphaHat
    ##
    ## for(k in 1:K) {
    ##     alphaHat_t[k, 1] <- priorCurrent[k] * eMatHapSNP_t[k, 1]
    ## }
    ## c[1] <- 1 / sum(alphaHat_t[, 1])
    ## alphaHat_t[, 1] <- alphaHat_t[, 1] * c[1]
    Rcpp_run_forward_haploid(
        alphaHat_t = alphaHat_t,
        c = c,
        eMatGrid_t = eMatGrid_t,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        priorCurrent_m = priorCurrent_m,
        s = s - 1,
        alphaStart = 0,
        run_fb_subset = FALSE,
        initialize_only = FALSE
    )
    ##
    ## do betahat
    ##
    betaHat_t[, nGrids] <- c[nGrids]
    Rcpp_run_backward_haploid(        
        betaHat_t,
        c = c,
        eMatGrid_t = eMatGrid_t,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc = transMatRate_tc_H,
        s = s - 1
    )
    return(
        list(
            alphaHat_t = alphaHat_t,
            betaHat_t = betaHat_t,
            eMatGrid_t = eMatGrid_t,
            c = c
        )
    )
}


alpha_forward_one <- function(t, K, s, alphaHat_t, transMatRate_tc_H, eMatGrid_t, alphaMatCurrent_tc) {
    alphaConst <- 0;
    for(k in 1:K) {
        alphaConst <- alphaConst + alphaHat_t[k, t - 1]
    }
    alphaConst <- alphaConst * transMatRate_tc_H[2, t - 1, s]
    for(k in 1:K) {
        alphaHat_t[k, t] <-
            eMatGrid_t[k,t] *
            ( transMatRate_tc_H[1, t - 1, s] * alphaHat_t[k, t - 1] +  
              alphaConst * alphaMatCurrent_tc[k, t - 1, s])
    }
    return(alphaHat_t)
}



make_rlc <- function(ff) {
    if (length(ff) != 1) {
        stop("ff too long")
    }
    p <- c(0.5, (1 - ff) / 2, ff / 2)
    rlc <- matrix(0, nrow = 7, ncol = 3)
    rlc[1, ] <- c(1, 0, 0)
    rlc[2, ] <- c(0, 1, 0)
    rlc[3, ] <- c(0, 0, 1)
    rlc[4, ] <- c(p[1] / (p[1] + p[2]), p[2] / (p[1] + p[2]), 0)
    rlc[5, ] <- c(p[1] / (p[1] + p[3]), 0, p[3] / (p[1] + p[3]))
    rlc[6, ] <- c(0, p[2] / (p[2] + p[3]), p[3] / (p[2] + p[3]))
    rlc[7, ] <- p
    return(rlc)
}



determine_a_set_of_truth_labels_for_nipt <- function(
    sampleReads,
    truth_haps,
    phase,
    ff
) {
    readProbs_t <- assign_fetal_read_probabilities(
        sampleReads = sampleReads,
        truth_haps_t = t(truth_haps)
    )
    ## convert those probabilities into probabilities they originate from sets of haplotypes
    ## e.g. given ff, and a read could come from matt and p, what is the probability it came from each
    out <- get_read_groupings_given_fetal_fraction_and_cov(
        readProbs_t = readProbs_t,
        phase = phase,
        iiSample = 1,
        ff = ff
    )
    ## 
    truth_labels <- sample_H_for_NIPT_given_groupings(
        groupings = out$groupings,
        counts = out$counts,
        ff = ff
    )
    ## hmm, for now, make these exact
    groupings <- out$groupings
    uncertain_truth_labels <- rep(TRUE, nrow(out$groupings))
    uncertain_truth_labels[groupings[, "matt_only"]] <- FALSE
    uncertain_truth_labels[groupings[, "matu_only"]] <- FALSE
    uncertain_truth_labels[groupings[, "pat_only"]] <- FALSE   
    ## 
    list(
        truth_labels = truth_labels,
        uncertain_truth_labels = uncertain_truth_labels
    )
}


evaluate_read_variability <- function(eMatRead_t) {
    K <- nrow(eMatRead_t)
    nReads <- ncol(eMatRead_t)
    number_of_non_1_reads <- integer(nReads)
    indices_of_non_1_reads <- array(0, c(K, nReads))
    thresh <- 1 - 1e-12
    thresh2 <- as.integer(floor(K * 0.20)) ## threshold for "few" subtractions to do
    ## category 0 = normal (full mult)
    ## category 1 = can skip completely (all 1s)
    ## category 2 = subraction of specific entries, all with the same value (one non-1 value seen)
    ## category 3 = subraction of specific entries
    read_category <- integer(nReads)
    read_category[] <- 0L ## assume 
    for(iRead in 1:nReads) {
        c <- 1
        val <- -1
        more_than_two <- FALSE
        for(k in 1:K) {
            if (eMatRead_t[k, iRead] < (thresh)) {
                indices_of_non_1_reads[c, iRead] <- k - 1 ## make 0-based
                c <- c + 1
                if (val == -1) {
                    val <- eMatRead_t[k, iRead]
                } else {
                    if (eMatRead_t[k, iRead] != val) {
                        more_than_two <- TRUE
                    }
                }
            }
        }
        number_of_non_1_reads[iRead] <- c - 1
        if (number_of_non_1_reads[iRead] == 0) {
            read_category[iRead] <- 1L
        } else if (!more_than_two) {
            read_category[iRead] <- 2L
        } else if (number_of_non_1_reads[iRead] < thresh2) {
            read_category[iRead] <- 3L
        } else {
            read_category[iRead] <- 0L
        }
    }
    return(
        list(
            number_of_non_1_reads = number_of_non_1_reads,
            indices_of_non_1_reads = indices_of_non_1_reads,
            read_category = read_category
        )
    )
}



evaluate_read_probabilities <- function(
    alphaHat_m,
    betaHat_m,
    pC,
    read_category,
    iRead,
    h_rC,
    h_rA1,
    h_rA2,
    eMatRead_t,
    number_of_non_1_reads,
    indices_of_non_1_reads
) {
    ## need same three probabilities for flip 1, flip 2
    pA1 <- pC ## read label becomes h_rA1
    pA2 <- pC ## read label becomes h_rA2
    if (read_category[iRead] == 0L) {
        ## default behaviour
        pA1[c(h_rC, h_rA1)] <- 0
        pA2[c(h_rC, h_rA2)] <- 0
        K <- nrow(eMatRead_t)
        for(k in 1:K) {
            ## A1 loser
            pA1[h_rC]  <- pA1[h_rC]  + alphaHat_m[h_rC, k]  * betaHat_m[h_rC, k] / eMatRead_t[k, iRead]
            ## A1 gainer
            pA1[h_rA1] <- pA1[h_rA1] + alphaHat_m[h_rA1, k]  * betaHat_m[h_rA1, k] * eMatRead_t[k, iRead]
            ## A2 gainer
            pA2[h_rA2] <- pA2[h_rA2] + alphaHat_m[h_rA2, k]  * betaHat_m[h_rA2, k] * eMatRead_t[k, iRead]
        }
        ## A2 loser
        pA2[h_rC] <- pA1[h_rC]
        
    } else if (read_category[iRead] == 2L) {
        ##
        ## subtraction based approach, all the same
        ##
        ## setup already done with no change
        ## now remove and consider
        ## 
        ## now go back over, remove, and add back in
        ##
        val1 <- 0
        val2 <- 0
        val3 <- 0
        for(ik in 1:number_of_non_1_reads[iRead]) {
            k <- indices_of_non_1_reads[ik, iRead] + 1 ## make 1-based here            
            ## A1 loser
            val1 <- val1 + alphaHat_m[h_rC, k]  * betaHat_m[h_rC, k]
            ## A1 gainer
            val2 <- val2 + alphaHat_m[h_rA1, k]  * betaHat_m[h_rA1, k]
            ## A2 gainer
            val3 <- val3 + alphaHat_m[h_rA2, k]  * betaHat_m[h_rA2, k]
        }
        ## k at the end is fine
        pA1[h_rC] <-  pA1[h_rC]  + val1 * (1 / eMatRead_t[k, iRead] - 1)
        pA1[h_rA1] <- pA1[h_rA1] + val2 * (eMatRead_t[k, iRead] - 1)
        pA2[h_rA2] <- pA2[h_rA2]  + val3 * (eMatRead_t[k, iRead] - 1)
        ## A2 loser
        pA2[h_rC] <- pA1[h_rC]
    } else if (read_category[iRead] == 3L) {        
        ##
        ## subtraction based approach, all distinct
        ##
        ## setup already done with no change
        ## now remove and consider
        ## 
        ## now go back over, remove, and add back in
        ##
        for(ik in 1:number_of_non_1_reads[iRead]) {
            k <- indices_of_non_1_reads[ik, iRead] + 1 ## make 1-based here
            ## A1 loser
            pA1[h_rC]  <- pA1[h_rC]  + alphaHat_m[h_rC, k]  * betaHat_m[h_rC, k] * (1 / eMatRead_t[k, iRead] - 1)
            ## A1 gainer
            pA1[h_rA1] <- pA1[h_rA1] + alphaHat_m[h_rA1, k]  * betaHat_m[h_rA1, k] * (eMatRead_t[k, iRead] - 1)
            ## A2 gainer
            pA2[h_rA2] <- pA2[h_rA2] + alphaHat_m[h_rA2, k]  * betaHat_m[h_rA2, k] * (eMatRead_t[k, iRead] - 1)
        }
        ## A2 loser
        pA2[h_rC] <- pA1[h_rC]
    }
    return(list(pA1 = pA1, pA2 = pA2))
}

