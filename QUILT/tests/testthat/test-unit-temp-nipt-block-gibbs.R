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
    QUILT::rcpp_nth_partial_sort(runif(10), 3) ## I dunno
    
}

test_that("blarh", {

    skip("not for routine use")
    
    i_it <- 3
    nSNPs <- 18808
    
    ## load("/well/davies/users/dcc832/werAwerBwerC.3.18808.RData")
    load("/data/smew1/rdavies/nipt_test_2023_09_11/werAwerBwerC.3.18808.RData")

    
    suppressOutput <- 0
    shard_check_every_pair <- TRUE
   
    ## problem in add gammas and genProbs_t to output

    ##
    K <- length(which_haps_to_use)
    S <- 1
    ## print(paste0("start = ", Sys.time()))
    if (use_small_eHapsCurrent_tc & !use_provided_small_eHapsCurrent_tc) {
        inflate_fhb_t_in_place(
            rhb_t,
            small_eHapsCurrent_tc,
            haps_to_get = which_haps_to_use - 1,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
    }
    if (make_plots_block_gibbs) {
        return_gibbs_block_output <- TRUE
        return_advanced_gibbs_block_output <- TRUE
    }
    ##
    if (use_sample_is_diploid && ff == 0) {
        sample_is_diploid <- TRUE
    }else {
        sample_is_diploid <- FALSE
    }
    ## I dunno, could clearly be fixed
    if (is.null(eMatRead_t)) {
        pass_in_eMatRead_t <- FALSE
        eMatRead_t <- array(0, c(1, 1))
    } else {
        pass_in_eMatRead_t <- TRUE
    }
    param_list <- list(
        return_alpha = TRUE,
        return_extra = TRUE,
        return_genProbs = return_genProbs,
        return_gamma = as.logical(return_gamma | make_plots),
        return_hapProbs = as.logical(return_hapProbs | make_plots),
        return_p_store = return_p_store,
        return_p1 = return_p1,
        return_gibbs_block_output = return_gibbs_block_output,
        return_advanced_gibbs_block_output = return_advanced_gibbs_block_output,
        use_starting_read_labels = TRUE,
        verbose = verbose,
        run_fb_subset = FALSE,
        haploid_gibbs_equal_weighting = TRUE,
        gibbs_initialize_iteratively = gibbs_initialize_iteratively,
        gibbs_initialize_at_first_read = gibbs_initialize_at_first_read,
        use_smooth_cm_in_block_gibbs = use_smooth_cm_in_block_gibbs,
        use_small_eHapsCurrent_tc = use_small_eHapsCurrent_tc,
        sample_is_diploid = sample_is_diploid,
        update_in_place = FALSE,
        do_shard_block_gibbs = FALSE, ## for now
        shard_check_every_pair = TRUE,
        force_reset_read_category_zero = FALSE,
        disable_read_category_usage = disable_read_category_usage,
        calculate_gamma_on_the_fly = calculate_gamma_on_the_fly,
        rescale_eMatRead_t = rescale_eMatRead_t,
        pass_in_eMatRead_t = pass_in_eMatRead_t,
        make_eMatRead_t_rare_common = make_eMatRead_t_rare_common,
        pass_in_alphaBeta = TRUE,
        update_hapSum = FALSE,
        record_read_set = TRUE,
        perform_block_gibbs = perform_block_gibbs,
        use_eMatDH_special_symbols = use_eMatDH_special_symbols
    )
    ## use_provided_small_eHapsCurrent_tc = use_provided_small_eHapsCurrent_tc
    ## this should catch hopefully rare underflow problems and re-run the samples
    done_imputing <- FALSE
    n_imputing <- 0
    n_gibbs_burn_in_its <- small_ref_panel_gibbs_iterations
    n_gibbs_full_its <- n_gibbs_burn_in_its + n_gibbs_sample_its
    skip_read_iteration <- rep(FALSE, n_gibbs_full_its)
    if (small_ref_panel_skip_equally_likely_reads) {
        skip_read_iteration <- rep(TRUE, n_gibbs_full_its)
        skip_read_iteration[small_ref_panel_equally_likely_reads_update_iterations] <- FALSE
    }

    Ksubset <- 400
    nGrids <- 588
    nGrids 
                gammaMT_t_local <- array(0, c(Ksubset, nGrids))
                gammaMU_t_local <- array(0, c(Ksubset, nGrids))
                gammaP_t_local <- array(0, c(Ksubset, nGrids))
    

        out <- rcpp_forwardBackwardGibbsNIPT(
            sampleReads = sampleReads,
            eMatRead_t = eMatRead_t,
            priorCurrent_m = small_priorCurrent_m,
            alphaMatCurrent_tc = small_alphaMatCurrent_tc,
            eHapsCurrent_tc = small_eHapsCurrent_tc,
            transMatRate_tc_H = small_transMatRate_tc_H,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            use_hapMatcherR = use_hapMatcherR,
            distinctHapsB =distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            eMatDH_special_matrix = eMatDH_special_matrix,
            rhb_t = rhb_t,
            ref_error = ref_error,
            which_haps_to_use = which_haps_to_use,
            ff = ff,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax_local = Jmax,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            run_fb_grid_offset = FALSE,
            blocks_for_output = array(0, c(1, 1)),
            grid = grid,
            skip_read_iteration = skip_read_iteration,
            alphaHat_t1 = alphaHat_t1,
            alphaHat_t2 = alphaHat_t2,
            alphaHat_t3 = alphaHat_t3,
            betaHat_t1 = betaHat_t1,
            betaHat_t2 = betaHat_t2,
            betaHat_t3 = betaHat_t3,
            eMatGrid_t1 = eMatGrid_t1,
            eMatGrid_t2 = eMatGrid_t2,
            eMatGrid_t3 = eMatGrid_t3,
            gammaMT_t_local = gammaMT_t_local,
            gammaMU_t_local = gammaMU_t_local,
            gammaP_t_local = gammaP_t_local,
            hapSum_tc = array(0, c(1, 1, 1)),
            snp_start_1_based = -1,
            snp_end_1_based = -1,
            generate_fb_snp_offsets = FALSE,
            suppressOutput = 0,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            n_gibbs_sample_its = n_gibbs_sample_its,
            n_gibbs_starts = n_gibbs_starts,
            double_list_of_starting_read_labels = double_list_of_starting_read_labels,
            prev_list_of_alphaBetaBlocks = as.list(c(1, 2)),
            i_snp_block_for_alpha_beta = -1,
            do_block_resampling = FALSE, ## turn off for now
            seed_vector = 0,
            class_sum_cutoff = 0.06, ## what is this
            wif0 = wif0,
            grid_has_read = grid_has_read,
            shuffle_bin_radius = shuffle_bin_radius,
            L_grid = L_grid,
            block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
            smooth_cm = smooth_cm,
            param_list = param_list,
            block_gibbs_quantile_prob = block_gibbs_quantile_prob,
            artificial_relabel = -1,
            common_snp_index = common_snp_index,
            snp_is_common = snp_is_common,
            rare_per_hap_info = rare_per_hap_info,
            rare_per_snp_info = rare_per_snp_info
        )
        eMatRead_t <- out[["eMatRead_t"]]
        
        ## print(out$per_it_likelihoods[, -c(1, 2, 4, 5, 6, 7)])
        print(out$per_it_likelihoods[, c("i_it", "p_O_given_H_L", "p_H_given_L", "p_H_class_given_L")])
        
        ## print(names(out[["gibbs_block_output_list"]][[1]]))

        ## they get better but then worse?
        ## that seems really weird?
        ## is this because I haven't done read sampling afterwards based on that
        x1 <- out$gibbs_block_output_list[[1]]$gibbs_block_output$block_results        
        x2 <- out$gibbs_block_output_list[[2]]$gibbs_block_output$block_results
        x3 <- out$gibbs_block_output_list[[3]]$gibbs_block_output$block_results
        print(rbind(
            x1[c(1, nrow(x1) - 1), ],
            x2[c(1, nrow(x2) - 1), ],
            x3[c(1, nrow(x3) - 1), ]
        ))

        ## stop("WER")
        
        sapply(1:3, function(i) {
            m <- out[["gibbs_block_output_list"]][[i]][["shard_block_output"]][[1]]
            sum(diff(m[, "ir_chosen"]) != 0)
        })
        ## a bit better, but kind of crazy?


    truth_label_set <- determine_a_set_of_truth_labels_for_nipt(
                    sampleReads = sampleReads,
                    truth_haps = truth_haps,
                    phase = truth_haps,
                    ff = ff
                )
   truth_labels <-  truth_label_set$truth_labels
   uncertain_truth_labels <-  truth_label_set$uncertain_truth_labels

        fpp_stuff <- list(
            transMatRate_tc_H = small_transMatRate_tc_H,
            alphaMatCurrent_tc = small_alphaMatCurrent_tc,
            priorCurrent_m = small_priorCurrent_m ,
            eMatRead_t = eMatRead_t,
            s = 1,
            sampleReads = sampleReads
        )
    ## calculate true H class?
    initial_package <- for_testing_get_full_package_probabilities(truth_labels, fpp_stuff)
    truth_class <- calculate_H_class(
        eMatRead_t,
        alphaHat_t1 = initial_package[[1]]$alphaHat_t,
        alphaHat_t2 = initial_package[[2]]$alphaHat_t,
        alphaHat_t3 = initial_package[[3]]$alphaHat_t,
        betaHat_t1 = initial_package[[1]]$betaHat_t,
        betaHat_t2 = initial_package[[2]]$betaHat_t,
        betaHat_t3 = initial_package[[3]]$betaHat_t,
        ff = ff,
        wif0 = wif0,
        H = truth_labels,
        class_sum_cutoff = 1
    )

    
    ##     print(names(out))
        
    ##     print(table(out[["double_list_of_ending_read_labels"]][[1]]))
    ##     print(table(out[["double_list_of_ending_read_labels"]][[1]][[1]]))

    ##     x <- out[["gibbs_block_output_list"]][[1]]

        ## 
        ## how do I plot and/or know if it is good
        ##
        ## 

        ##
        ##try to plot block gibbs output
        ##
        f <- function(o, i) {
            x <- o[[paste0("alphaHat_t", i)]] * o[[paste0("betaHat_t", i)]]
            apply(x, 2, function(y) y / sum(y))
        }
        ## cheat, start with truth labels, add an artificial break or two, can I find it
        ## 


        for(i in 1:3) {
            outname <- paste0("~/robbie.test.", i, ".png")
            x <- out$gibbs_block_output_list[[1]]$gibbs_block_output[[1]]
            colnames(x)[colnames(x) == "ir_chosen"] <- "ir_left"
            x <- cbind(x, ir_right = x[, "ir_left"])
            shard_block_results <- x[seq(1, nrow(x), 2), ]
            before_read_labels <- out$gibbs_block_output_list[[i]]$before_read_labels
            after_read_labels <- out$gibbs_block_output_list[[i]]$after_read_labels                
            before_gamma1_t <- out$gibbs_block_output_list[[i]]$before_gamma1_t
            before_gamma2_t <- out$gibbs_block_output_list[[i]]$before_gamma2_t
            before_gamma3_t <- out$gibbs_block_output_list[[i]]$before_gamma3_t
            after_gamma1_t <- out$gibbs_block_output_list[[i]]$after_gamma1_t
            after_gamma2_t <- out$gibbs_block_output_list[[i]]$after_gamma2_t
            after_gamma3_t <- out$gibbs_block_output_list[[i]]$after_gamma3_t
            before_H_class <- out$gibbs_block_output_list[[i]]$before_H_class
            after_H_class <- out$gibbs_block_output_list[[i]]$after_H_class
            ## 
            plot_shard_block_output(
                shard_block_results = shard_block_results,
                before_gamma1_t = before_gamma1_t,
                before_gamma2_t = before_gamma2_t,
                before_gamma3_t = before_gamma3_t,
                after_gamma1_t = after_gamma1_t,
                after_gamma2_t = after_gamma2_t,
                after_gamma3_t = after_gamma3_t,
                before_read_labels = before_read_labels,
                after_read_labels = after_read_labels,
                nGrids = nGrids,
                outname = outname,
                L_grid = L_grid,
                L = L,
                uncertain_truth_labels = uncertain_truth_labels,
                truth_labels = truth_labels,
                have_truth_haplotypes = have_truth_haplotypes,
                sampleReads = sampleReads,
                wif0 = wif0,
                grid = grid,
                method = method,
                ff = ff,
                only_plot_confident_reads = TRUE,
                before_H_class = before_H_class,
                after_H_class = after_H_class
            )
        }

        

        
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


    ## make some packages here
 
   for(iii in 1:2) {

        if (iii == 1) {
            ## choose before labels, the ACTUAL ones we're looking at
            H <- out[["double_list_of_ending_read_labels"]][[1]][[1]]
            file <- "/well/davies/users/dcc832/werAwerBwerC.package.real.RData"            
        } else {
            ## choose truth labels, don't really expect a change!
            H <- truth_labels
            file <- "/well/davies/users/dcc832/werAwerBwerC.package.truth.RData"
        }
        package <- for_testing_get_full_package_probabilities(H, fpp_stuff, ff = ff)
        print("saving")
        save(
            H, ff, nGrids,
            package, eMatRead_t, wif0,
            sampleReads,
            truth_labels, uncertain_truth_labels,
            L, L_grid, grid, 
            small_alphaMatCurrent_tc,
            small_priorCurrent_m,
            small_transMatRate_tc_H,
            fpp_stuff,
            file = file
        )
    
        out_new <- R_shard_block_gibbs_resampler(
            alphaHat_t1 = package[[1]]$alphaHat_t,
            alphaHat_t2 = package[[2]]$alphaHat_t,
            alphaHat_t3 = package[[3]]$alphaHat_t,
            betaHat_t1 = package[[1]]$betaHat_t,
            betaHat_t2 = package[[2]]$betaHat_t,
            betaHat_t3 = package[[3]]$betaHat_t,
            c1 = package[[1]]$c,
            c2 = package[[2]]$c,
            c3 = package[[3]]$c,
            eMatGrid_t1 = package[[1]]$eMatGrid_t,
            eMatGrid_t2 = package[[2]]$eMatGrid_t,
            eMatGrid_t3 = package[[3]]$eMatGrid_t,
            H = H,
            ff = ff,
            eMatRead_t = eMatRead_t,
            grid = grid,
            wif0 = wif0,
            s = 1,
            alphaMatCurrent_tc = small_alphaMatCurrent_tc,
            priorCurrent_m = small_priorCurrent_m,
            transMatRate_tc_H = small_transMatRate_tc_H,
            do_checks = FALSE,
            verbose = FALSE,
            fpp_stuff = NULL,
            shard_check_every_pair = TRUE,
            H_class = package[["H_class"]]
        )

        print(package$log_p)
        print(-(sum(log(out_new$c1)) + sum(log(out_new$c2)) + sum(log(out_new$c3))))
    }
    

        ## VERY similar

        ## ACTUALLY a bit better!!! 


        ## what about the read probabilities, re-sampling, that kind of thing
        
        ## maybe just make a new simple plot here
        ## before shard gammas
        ## after ghard gammas
        ## read labels beforehand coloured
        ## real labels after coloured

        f <- function(o, i) {
            x <- o[[paste0("alphaHat_t", i)]] * o[[paste0("betaHat_t", i)]]
            apply(x, 2, function(y) y / sum(y))
        }

        ## cheat, start with truth labels, add an artificial break or two, can I find it

        ## before_read_labels <- out[["double_list_of_ending_read_labels"]][[1]][[1]]                
    before_gamma1_t <- f(out, 1)
    before_gamma2_t <- f(out, 2)
    before_gamma3_t <- f(out, 3)
    after_read_labels <- out_new[["H"]]
    shard_block_results <- out_new[["shard_block_results"]]
    after_gamma1_t <- f(out_new, 1)
    after_gamma2_t <- f(out_new, 2)
    after_gamma3_t <- f(out_new, 3)

    ## before_H_class <- out[["H_class"]]
    before_H_class <- rep(0, length(sampleReads))


    ##
    ## check out the triploid ones
    ##
    nReads <- length(sampleReads)
    ori_sampleReads <- sampleReads
    ## load("/well/davies/users/dcc832/temp.triploid1.RData")    
    load("/data/smew1/rdavies/nipt_test_2023_09_11/temp.triploid1.RData")
    
    u <- runif(nReads)
    triploid_read_labels <- rep(3, nReads)
    triploid_read_labels[u < (readLabelProbs[1, ] + readLabelProbs[2, ])] <- 2
    triploid_read_labels[u < (readLabelProbs[1, ])] <- 1    
                after_gamma1_t = gamma1
                after_gamma2_t = gamma2
                after_gamma3_t = gamma3

                shard_block_results <- array(1, c(nGrids - 1, 7))
                colnames(shard_block_results) <- c("ir_chosen", "p1", "p2", "p3", "p4", "p5", "p6")
                shard_block_results[1:3] <- c(1, 2, 1)
                shard_block_results[, 2:7] <- 1/6
                ## jesus christ?
                ## 
                initial_package <- for_testing_get_full_package_probabilities(truth_labels, fpp_stuff)
                initial_package$log_p
                before_gamma1_t <- initial_package[[1]]$gamma_t
                before_gamma2_t <- initial_package[[2]]$gamma_t
                before_gamma3_t <- initial_package[[3]]$gamma_t
                before_read_labels <- truth_labels
                ## 
                initial_package <- for_testing_get_full_package_probabilities(triploid_read_labels, fpp_stuff)
                initial_package$log_p
                after_read_labels <- triploid_read_labels
                after_gamma1_t <- initial_package[[1]]$gamma_t
                after_gamma2_t <- initial_package[[2]]$gamma_t
                after_gamma3_t <- initial_package[[3]]$gamma_t
                truth_class <- rep(1, length(sampleReads))
                
    
    
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
    outname <- "~/temp.png"

            plot_shard_block_output(
                shard_block_results = shard_block_results,
                before_gamma1_t = before_gamma1_t,
                before_gamma2_t = before_gamma2_t,
                before_gamma3_t = before_gamma3_t,
                after_gamma1_t = after_gamma1_t,
                after_gamma2_t = after_gamma2_t,
                after_gamma3_t = after_gamma3_t,
                before_read_labels = before_read_labels,
                after_read_labels = after_read_labels,
                nGrids = nGrids,
                outname = outname,
                L_grid = L_grid,
                L = L,
                uncertain_truth_labels = uncertain_truth_labels,
                truth_labels = truth_labels,
                have_truth_haplotypes = have_truth_haplotypes,
                sampleReads = sampleReads,
                wif0 = wif0,
                grid = grid,
                method = method,
                ff = ff,
                only_plot_confident_reads = TRUE,
                before_H_class = truth_class
            )       

            1

            ## do I need to colour them
            table(truth_labels, uncertain_truth_labels)
            
    ##     ## yup, has made some suggestions!
    ##     out_new$shard_block_results[, "ir_chosen"]

    ##     ##
    ##     sum(log(out[["c1"]])) + sum(log(out[["c2"]])) + sum(log(out[["c3"]]))
    ##     sum(log(out_new[["c1"]])) + sum(log(out_new[["c2"]])) + sum(log(out_new[["c3"]])) ## lower? good?

    ##     ## all of them are lower, which is promising!
    ##     sum(diff(apply(out[["alphaHat_t2"]] * out[["betaHat_t2"]], 2, which.max)) != 0) ## higher
    ##     sum(diff(apply(out_new[["alphaHat_t2"]] * out_new[["betaHat_t2"]], 2, which.max)) != 0) ## lower, so promising

    ##     sum(diff(apply(out[["alphaHat_t3"]] * out[["betaHat_t3"]], 2, which.max)) != 0) ## higher
    ##     sum(diff(apply(out_new[["alphaHat_t3"]] * out_new[["betaHat_t3"]], 2, which.max)) != 0) ## lower, so promising


            ## what about within a block
            table(before_H_class, truth_labels)

            ## try manually making some breaks, checking out likelihoods

        fpp_stuff <- list(
            transMatRate_tc_H = small_transMatRate_tc_H,
            alphaMatCurrent_tc = small_alphaMatCurrent_tc,
            priorCurrent_m = small_priorCurrent_m ,
            eMatRead_t = out$eMatRead_t,
            s = 1,
            sampleReads = sampleReads
        )
        H <- after_read_labels
        initial_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)        
        initial_package$log_p
    rx <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(3, 1, 2),
        c(2, 3, 1),
        c(3, 2, 1)
    )

        ## what if I flip this 
    ma <- mclapply(487 + -3:3, mc.cores = 4, function(iSplit) {
            ## keep up to that point, then do the splits after words, get likelihoods
            w <- (wif0 > iSplit)
        x <- sapply(1:6, function(ir) {
                H <- after_read_labels
                H[w] <- rx[ir, H[w]]
                initial_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)        
                initial_package$log_p
        })
        ## x + -x[1]
        x
        })
    t(sapply(ma, I)) ## yeah, one better, rest similar

    ## ARGH FUCK WAY OFF THERE

    ## CAN I FIND EASY INDUCED BREAKS?
    initial_package$log_p
    
    ## calculate true H class?
    initial_package <- for_testing_get_full_package_probabilities(truth_labels, fpp_stuff)
    truth_class <- calculate_H_class(
        eMatRead_t,
        alphaHat_t1 = initial_package[[1]]$alphaHat_t,
        alphaHat_t2 = initial_package[[2]]$alphaHat_t,
        alphaHat_t3 = initial_package[[3]]$alphaHat_t,
        betaHat_t1 = initial_package[[1]]$betaHat_t,
        betaHat_t2 = initial_package[[2]]$betaHat_t,
        betaHat_t3 = initial_package[[3]]$betaHat_t,
        ff = ff,
        wif0 = wif0,
        H = truth_labels,
        class_sum_cutoff = 0.06
    )

    table(truth_class, before_H_class)
    

    ## so ~83% of truth labels are uncertain
    table(uncertain_truth_labels)    
    table(truth_class, truth_labels)
    table(truth_class[!uncertain_truth_labels], truth_labels[!uncertain_truth_labels])

    
})
