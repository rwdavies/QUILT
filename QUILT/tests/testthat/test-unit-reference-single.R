if ( 1 == 0 ) {
    
    library("testthat")
    library("QUILT")
    dir <- "~/proj/QUILT/"
    setwd(paste0(dir, "/QUILT/R"))
    a <- dir(pattern = "*R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    
}


test_that("can build necessary components from make_rhb_t_equality", {

    ## make them mostly one of three options, with a few small changes
    K <- 500
    nSNPs <- 100
    reference_haps <- array(as.integer(runif(nSNPs * K) > 0.5), c(nSNPs, K))
    rhi <- reference_haps
    rhi_t <- t(rhi)
    rhb_t <- STITCH::make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)
    
    ref_error <- 0.01
    ref_one_minus_error <- 1 - ref_error

    for(nMaxDH in c(3, 255)) {
    
        ## make haplotype matching objects
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
        ##
        ## perform
        ## 
        expect_equal(sum(eMatDH_special_grid_which != 0), length(eMatDH_special_values_list))
        if (length(eMatDH_special_values_list) == 0) {
            expect_equal(nrow_which_hapMatcher_0, 0)
        } else {
            expect_equal(nrow_which_hapMatcher_0, sum(sapply(eMatDH_special_values_list, length)))
        }

    }

})


test_that("can run a single gl sample through reference haplotypes quickly with grid of 32", {

    ## switch on and off. off does thes tests, on checks cpp times
    speed_test <- FALSE
    if (!speed_test) {
        nSNPs <- 100
        L <- sort(sample(1:1000, 100))
        suppressOutput <- 1        
    } else {
        suppressOutput <- 0
        nSNPs <- 1000
        L <- sort(sample(1:10000, 1000))
        expect_equal(1, 1)
    }

    ## looks good
    out <- assign_positions_to_grid(
        L = L,
        grid32 = TRUE,
        gridWindowSize = NA
    )
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids
    snps_in_grid_1_based <- out$snps_in_grid_1_based
    cM_grid <- out$cM_grid

    ##
    dl <- diff(L_grid)
    expRate <- 1
    nGen <- 10
    sigmaCurrent <- exp(-nGen * expRate / 100 / 1000000 * dl)
    transMatRate_t <- rbind(sigmaCurrent, 1 - sigmaCurrent)


    use_eMatDH <- TRUE; i_setup <- 2

    gammaSmall_cols_to_get <- c(-1, 0, -1, 1) ## i.e. 0-based, what to put it in
    nSmallGammaGrids <- 2
    K_top_matches <- 5
    get_best_haps_from_thinned_sites <- TRUE
    
    ## build some haplotypes and encode them
    for(use_eMatDH in c(FALSE, TRUE)) {        
        for(i_setup in 1:3) {

            ## for cpp
            return_dosage <- TRUE
            return_betaHat_t <- TRUE
            return_gamma_t <- TRUE
            return_gammaSmall_t <- TRUE
            get_best_haps_from_thinned_sites <- TRUE

            ##
            ## i_setup = 1 --- 
            ## i_setup = 2 --- 
            ## i_setup = 3 --- same as above, but only get best_haps from thinned sites
            if (!suppressOutput) {            
                print(paste0("use_eMatDH = ", use_eMatDH, ", i_setup = ", i_setup))
            }
            if (i_setup == 1) {
                ## test gamma makes sense
                set.seed(9910)
                K <- 100
                reference_haps <- array(as.integer(runif(nSNPs * K) > 0.5), c(nSNPs, K))
            } else if (i_setup == 2  | i_setup == 3){
                ## test can easily work with large, mostly similar emissions
                set.seed(242)
                K <- 10000
                reference_haps <- array(NA, c(nSNPs, K))
                haps <- array(sample(c(0, 1), nSNPs * 3, replace = TRUE), c(3, nSNPs))
                for(k in 1:K) {
                    reference_haps[, k] <- haps[sample(1:3, 1), ]
                }
                ## tiny noise
                w <- runif(K * nSNPs) < 0.01
                reference_haps[w] <- 1 - reference_haps[w]
                if (i_setup == 3) {
                    ## particularly for cpp
                    return_dosage <- FALSE
                    return_betaHat_t <- FALSE
                    return_gamma_t <- FALSE
                    return_gammaSmall_t <- FALSE
                    get_best_haps_from_thinned_sites <- TRUE
                }
            }

            ## 
            best_haps_stuff_list <- list()
            best_haps_stuff_list <- as.list(1:sum(gammaSmall_cols_to_get >= 0))
            
            ## make them mostly one of three options, with a few small changes
            rhi <- reference_haps
            rhi_t <- t(rhi)
            rhb_t <- STITCH::make_rhb_t_from_rhi_t(rhi_t)
            rhb <- t(rhb_t)

            nMaxDH <- 2 ** 10 - 1
            ref_error <- 0.01
            ref_one_minus_error <- 1 - ref_error

            
            ## make haplotype matching objects
            out <- make_rhb_t_equality(
                rhb_t = rhb_t,
                nMaxDH = nMaxDH,
                nSNPs = nSNPs,
                ref_error = ref_error
            )
            distinctHapsB <- out[["distinctHapsB"]]
            distinctHapsIE <- out[["distinctHapsIE"]]            
            hapMatcher <- out[["hapMatcher"]]
            
            ## make some gls from this
            my_hap <- c(
                rhi_t[1, 1:36],
                rhi_t[2, 37:60],
                rhi_t[3, 61:100]
            )
            ## my_hap <- c(
            ##     rhi_t[1, 1:32],
            ##     rhi_t[2, 33:64],
            ##     rhi_t[3, 65:100]
            ## )
            ## draw some reads from this
            ## draw some reads from this
            nReads <- round(nSNPs * 1.5)
            u <- sort(sample(1:nSNPs, nReads, replace = TRUE))
            bq <- rep(-10, length(u))
            bq[my_hap[u] == 1] <- 10
            gl <- make_gl_from_u_bq(u, bq, nSNPs) 
            ## add a bit of noise
            ## w <- sample(1:length(bq), round(length(bq) / 10))
            ## bq[w] <- - bq[w]
            ## sample(c(-10, 10), nReads, replace = TRUE)
            ## round(t(rbind(gl[, 1:32], my_hap[1:32])), 3)
            ## from this, make prob of ref and prob of alt emissio for each
            ## my_hap <- 0.90 * my_hap + 0.10 * runif(nSNPs)
            rm(reference_haps, rhi_t, rhi); gc(reset = TRUE); gc(reset = TRUE)

            if (!speed_test) {
                ## check the same
                R_eMatDH <- build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error)
                Rcpp_eMatDH <- Rcpp_build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error)
                expect_equal(R_eMatDH, Rcpp_eMatDH)
            }

            K <- nrow(rhb_t)
            nGrids <- ncol(rhb_t)
            ##
            alphaHat_t <- array(0, c(K, nGrids))
            betaHat_t <- array(0, c(K, nGrids))
            gamma_t <- array(0, c(K, nGrids))
            gammaSmall_t <- array(0, c(K, nSmallGammaGrids))
            dosage <- numeric(nSNPs)
            
            ## so now, want to 
            if (!speed_test) {
                
                ## 
                a <- system.time(
                    outR <- R_haploid_dosage_versus_refs(
                        gl = gl,
                        alphaHat_t = alphaHat_t,
                        betaHat_t = betaHat_t,
                        gamma_t = gamma_t,
                        gammaSmall_t = gammaSmall_t,
                        dosage = dosage,
                        transMatRate_t = transMatRate_t,
                        rhb_t = rhb_t,
                        ref_error = ref_error,
                        use_eMatDH = use_eMatDH,
                        distinctHapsB = distinctHapsB,
                        distinctHapsIE = distinctHapsIE,
                        hapMatcher = hapMatcher,
                        return_extra = TRUE,
                        gammaSmall_cols_to_get = gammaSmall_cols_to_get,
                        get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
                        return_dosage = TRUE,
                        return_gammaSmall_t = TRUE,
                        return_gamma_t = TRUE,
                        K_top_matches = K_top_matches
                    )
                )

                if (!suppressOutput) {                
                    print("R speed")
                }
                R_gamma_t <- outR[["gamma_t"]]
                R_gammaSmall_t <- outR[["gammaSmall_t"]]
                R_dosage <- outR[["dosage"]]
                R_alphaHat_t <- outR[["alphaHat_t"]]
                R_betaHat_t <- outR[["betaHat_t"]]
                R_best_haps_stuff_list <- outR[["best_haps_stuff_list"]]
                expect_equal(colSums(R_gamma_t), rep(1, 4))
                expect_true((cor(R_dosage, my_hap) ** 2) > 0.7) ## is better for truer match against reference
                ## alternate
                
            }
            
            a <- system.time(
                Rcpp_haploid_dosage_versus_refs(
                    gl = gl,
                    alphaHat_t = alphaHat_t,
                    betaHat_t = betaHat_t,
                    gamma_t = gamma_t,
                    dosage = dosage,
                    transMatRate_t = transMatRate_t,
                    rhb_t = rhb_t,
                    ref_error = ref_error,
                    use_eMatDH = use_eMatDH,
                    distinctHapsB = distinctHapsB,
                    distinctHapsIE = distinctHapsIE,
                    hapMatcher = hapMatcher,
                    gammaSmall_t = gammaSmall_t,
                    gammaSmall_cols_to_get = gammaSmall_cols_to_get,
                    suppressOutput = suppressOutput,
                    K_top_matches = K_top_matches,
                    best_haps_stuff_list = best_haps_stuff_list,
                    get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
                    return_betaHat_t = return_betaHat_t,
                    return_dosage = return_dosage,
                    return_gamma_t = return_gamma_t,
                    return_gammaSmall_t = return_gammaSmall_t
                )
            )
            if (!suppressOutput) {
                print("cpp speed")
                print(a)
            }

            if (!speed_test) {
                if (!suppressOutput) {
                    print(use_eMatDH)
                }
                if (i_setup <= 2) {
                    expect_equal(R_alphaHat_t, alphaHat_t)
                    expect_equal(R_betaHat_t, betaHat_t)
                    expect_equal(R_gamma_t, gamma_t)
                    expect_equal(R_dosage, dosage)
                }
                expect_equal(R_best_haps_stuff_list, best_haps_stuff_list)
            }

        }
    }

    
})


test_that("can fast eMat", {


    ## simulate some stuff
    K <- 100
    nSNPs <- 200
    nReads <- 20
    out <- make_fb_test_package(K = K, nReads = nReads, nSNPs = nSNPs)

    rhi <- t(round(out$eHapsCurrent_tc[, , 1]))
    rhi_t <- t(rhi)
    rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)

    ref_error <- 0.01
    sampleReads <- out$sampleReads
    nSNPs <- out$nSNPs
    eHapsCurrent_tc <- array(0, c(K, nSNPs, 1))
    eHapsCurrent_tc[] <- rhi_t
    eHapsCurrent_tc[eHapsCurrent_tc == 0] <- ref_error
    eHapsCurrent_tc[eHapsCurrent_tc == 1] <- 1 - ref_error
    
    maxDifferenceBetweenReads <- 1e10
    normal_eMatRead_t <- array(1, c(K, nReads))
    rcpp_make_eMatRead_t(
        eMatRead_t = normal_eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = eHapsCurrent_tc,
        s = 0,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads ,
        Jmax = 1000,
        eMatHapOri_t = array(0, c(1, 1)),
        pRgivenH1 = array(0, 1),
        pRgivenH2 = array(0, 1),
        prev = 0,
        suppressOutput = 1,
        prev_section = "wer",
        next_section = "wer",
        run_pseudo_haploid = FALSE,
        rescale_eMatRead_t = FALSE
    )

    eMatRead_t <- make_eMatRead_t_using_binary(
        sampleReads = sampleReads,
        rhb_t = rhb_t,
        nSNPs = nSNPs,
        ref_error = ref_error,
        language = "R",
        n = 10
    )
    expect_equal(eMatRead_t, normal_eMatRead_t)    

    ## OK!


    Rcpp_eMatRead_t <- make_eMatRead_t_using_binary(
        sampleReads = sampleReads,
        rhb_t = rhb_t,
        nSNPs = nSNPs,
        ref_error = ref_error,
        language = "Rcpp"
    )

    expect_equal(eMatRead_t, Rcpp_eMatRead_t)    
    
})


test_that("profile", {

    skip("skip profiling")
    
    if (1 == 0) {

        load("~/Downloads/impute_develop.RData")
        nMaxDH <- 2 ** 8 - 1
        ref_error <- 1e-2
        out <- make_rhb_t_equality(
            rhb_t = rhb_t,
            nMaxDH = nMaxDH,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
        distinctHapsB <- out[["distinctHapsB"]]
        distinctHapsIE <- out[["distinctHapsIE"]]            
        hapMatcher <- out[["hapMatcher"]]
        save(nMaxDH, ref_error, distinctHapsB, distinctHapsIE, hapMatcher, nSNPs, transMatRate_t, gl, rhb_t, file = "~/Downloads/temp_quilt_speed_test.RData", compress = FALSE)

    }

    ## 
    ## load("~/Downloads/temp_quilt_speed_test.RData")
    load("~/Downloads/impute_develop.RData")
    transMatRate_t <- full_transMatRate_t_H
    K <- nrow(rhb_t)    
    nGrids <- ncol(rhb_t)
    i_hap <- 1
    u <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
    bq <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
    w <- bq != 0
    bq <- bq[w]
    u <- u[w]
    if (length(u) == 0) {
        print_message(paste0("Read label assignment includes no reads for haplotype ", i_hap))
    }
    ## 
    gl <- make_gl_from_u_bq(u, bq, nSNPs)
    ## so u is 1-based, so e.g. having 1-32 means grid 1, 33-64 means grid 2, etc
    ## so unique(floor(u / 32)) having 2 means yes there is a var in grid 2, etc
    grid_has_variant <- array(FALSE, nGrids)
    grid_has_variant[unique(floor(u / 32))] <- TRUE ## at least for 32-base grids

    which_hapMatcher_0 <- which(hapMatcher == 0, arr.ind = TRUE) 
    which_hapMatcher_0 <- which_hapMatcher_0 - 1 ## make both 0-based


                ## // now do additional check on stragglers
		## // if (which_hapMatcher_0(hapMatcher_0_position, 1) == iGrid) {
		## //     while(continue_hapMatcher_0 && which_hapMatcher_0(hapMatcher_0_position, 1) == iGrid) {
		## //         k = which_hapMatcher_0(hapMatcher_0_position, 0);
		## //         std::uint32_t tmp(rhb_t(k, iGrid));                
                ## //         //
                ## //         prob = 1;
                ## //         for(int b = 0; b < nSNPsLocal; b++) {
                ## //             dR = gl_local(0, b);
                ## //             dA = gl_local(1, b);
                ## //             if (tmp & (1<<b)) {
                ## //                 // alternate
                ## //                 prob *= (dR * ref_error + dA * ref_one_minus_error);
                ## //             } else {
                ## //                 prob *= (dR * ref_one_minus_error + dA * ref_error);
                ## //             }
                ## //         }
		## // 	alphaHat_t_col(k) = (jump_prob + not_jump_prob * alphaHat_t_col(k)) * prob;
		## // 	run_total += alphaHat_t_col(k);
		## // 	hapMatcher_0_position++;
		## // 	if (hapMatcher_0_position == n_which_hapMatcher_0) {
  		## // 	    continue_hapMatcher_0 = false;
		## // 	    hapMatcher_0_position = 0;			    
		## // 	}
		## //     }
		## // }
		## // if (iGrid > which_hapMatcher_0(hapMatcher_0_position, 1)) {
		## //   while(continue_hapMatcher_0 && (iGrid > which_hapMatcher_0(hapMatcher_0_position, 1))) {
		## // 	hapMatcher_0_position++;
		## // 	if (hapMatcher_0_position == n_which_hapMatcher_0) {
  		## // 	    continue_hapMatcher_0 = false;
		## // 	    hapMatcher_0_position = 0;
		## // 	}
		## //     }
		## // }
		## //
    
    
    ## nMaxDH <- 2 ** 10 - 1
    ## ref_error <- 1e-2

    ## so about a third or more of sites are invariant, so prob is constant
    ## what is the best way to get this, through u?



    ##
    ##
    ## profile forward
    ##
    ##
    if (1 == 0) {
        
        print("------------------")
        use_eMatDH <- TRUE
        ref_one_minus_error <- 1 - ref_error
        eMatDH <- Rcpp_build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error)
        alphaHat_t_new<- array(0, c(K, nGrids))
        alphaHat_t_original <- array(0, c(K, nGrids))    
        ##betaHat_t <- array(0, c(K, nGrids))
        ##gamma_t <- array(0, c(K, nGrids))
        betaHat_t <- array(0, c(1, 1))
        gamma_t <- array(0, c(1, 1))
        gammaSmall_t <- array(0, c(1, 1))
        print(paste0("----------new version----------------"))
        ## kind of want to be one bigger, 
        ##hapMatcher0 <- hapMatcher - 1
        ##hapMatcher0[hapMatcher0 == -1] <- nMaxDH - 1
        eMatDH_bigger <- rbind(0, eMatDH) ## make one bigger here for testing
        c1 <- array(0, c(1, nGrids))
        hapMatcher0 <- hapMatcher == 0
        ww <- seq(1, nGrids, length.out = max(1, round(heuristic_match_thin * nGrids)))
        full_gammaSmall_cols_to_get <- array(-1, nGrids)
        full_gammaSmall_cols_to_get[ww] <- 0:(length(ww) - 1)
        return_gammaSmall_t <- TRUE
        ##
        only_store_alpha_at_gamma_small <- FALSE
        ## try to store where the hapMatcher0 are
        
        Rcpp_haploid_reference_single_forward_version2(
            full_gammaSmall_cols_to_get,
            gl,
            alphaHat_t_new,
            c1,
            transMatRate_t,
            rhb_t,
            hapMatcher,
            eMatDH_bigger,
            nGrids,
            nSNPs,
            K,
            use_eMatDH,
            ref_error,
            only_store_alpha_at_gamma_small
        )
        alphaHat_t <- alphaHat_t_new
        c <- c1
        ## save(alphaHat_t, c, file = "~/Downloads/impute_develop_alpha.RData")
        print(range(c1))
        print(paste0("----------original version----------------"))
        c2 <- array(0, c(1, nGrids))
        Rcpp_haploid_reference_single_forward(
            full_gammaSmall_cols_to_get,
            gl,
            alphaHat_t_original,
            c2,
            transMatRate_t,
            rhb_t,
            hapMatcher,
            eMatDH,
            nGrids,
            nSNPs,
            K,
            use_eMatDH,
            ref_error,
            only_store_alpha_at_gamma_small
        )
        print(c1[1:10])
        print(c2[1:10])
        print(grid_has_variant[1:10])
        print(sum(log10(c1[-1])))
        print(sum(log10(c2[-1])))
        print(prod((c1[-1])))
        print(prod((c2[-1])))

        print(median(alphaHat_t_new[, nGrids]))
        print(median(alphaHat_t_original[, nGrids]))
        
        expect_equal(alphaHat_t_new[, nGrids], alphaHat_t_original[, nGrids])
        expect_equal(c1, c2)

        stop("for now that is it")

    }


    if (1 == 0) {

        print("----------------")
        load(file = "~/Downloads/impute_develop_alpha.RData")
        nMaxDH <- 2 ** 8 - 1        
        return_betaHat_t <- FALSE
        return_dosage <- FALSE
        return_gamma_t <- FALSE
        return_gammaSmall_t <- TRUE
        make_plots <- FALSE
        full_betaHat_t <- array(0, c(1, 1))
        if (make_plots) {
            full_gamma_t <- array(0, c(K, nGrids))
        } else {
            full_gamma_t <- array(0, c(1, 1))        
        }
        iSample <- 1
        ww <- seq(1, nGrids, length.out = max(1, round(heuristic_match_thin * nGrids)))
        full_gammaSmall_cols_to_get <- array(-1, nGrids)
        full_gammaSmall_cols_to_get[ww] <- 0:(length(ww) - 1)
        use_eMatDH <- TRUE
        ref_one_minus_error <- 1 - ref_error
        add_zero_row <- TRUE
        eMatDH_bigger <- Rcpp_build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error, add_zero_row)
        add_zero_row <- FALSE
        eMatDH_smaller <- Rcpp_build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error, add_zero_row)            

        print("--- backwards new version ---------------")
        dosage1 <- numeric(nSNPs)
        full_gammaSmall_t1 <-  array(0, c(K, length(ww)))
        Rcpp_haploid_reference_single_backward_version2(
            alphaHat_t,
            full_betaHat_t,
            full_gamma_t,
            full_gammaSmall_t1,
            full_gammaSmall_cols_to_get,
            dosage1,
            nGrids,
            transMatRate_t,
            eMatDH_bigger,
            hapMatcher,
            nSNPs,
            K,
            use_eMatDH,
            rhb_t,
            ref_error,
            gl,
            c,
            distinctHapsIE,
            return_betaHat_t,
            return_dosage,
            return_gamma_t,
            return_gammaSmall_t,
            nMaxDH
        )


        print("--- backwards old version ---------------")
        dosage2 <- numeric(nSNPs)
        full_gammaSmall_t2 <-  array(0, c(K, length(ww)))        
        Rcpp_haploid_reference_single_backward(
            alphaHat_t,
            full_betaHat_t,
            full_gamma_t,
            full_gammaSmall_t2,
            full_gammaSmall_cols_to_get,
            dosage2,    
            nGrids,
            transMatRate_t,
            eMatDH_smaller,
            hapMatcher,
            nSNPs,
            K,
            use_eMatDH,
            rhb_t,
            ref_error,
            gl,
            c,
            distinctHapsIE,
            return_betaHat_t,
            return_dosage,
            return_gamma_t,
            return_gammaSmall_t,
            nMaxDH
        )

        print(mean(dosage1))
        print(mean(dosage2))
        expect_equal(dosage1, dosage2)
        print(mean(full_gammaSmall_t1))
        print(mean(full_gammaSmall_t2))
        expect_equal(full_gammaSmall_t1, full_gammaSmall_t2)
        stop("END OF BETA BIT")        
    }





    
    
    for(i in 1:2) {

        use_eMatDH <- TRUE
        print(paste0("----------i = ", i))
        alphaHat_t <- array(0, c(K, nGrids))
        betaHat_t <- array(0, c(1, 1))
        gamma_t <- array(0, c(1, 1))
        ww <- seq(1, nGrids, length.out = max(1, round(heuristic_match_thin * nGrids)))
        full_gammaSmall_cols_to_get <- array(-1, nGrids)
        full_gammaSmall_cols_to_get[ww] <- 0:(length(ww) - 1)
        full_gammaSmall_t1 <-  array(0, c(K, length(ww)))
        
        ## can I run some arbitrary number of gls, need many alphaHats for instance
        ## can I run BOTH gls as the same time for instance? does this matter
        for(i_hap in 1:1) {
            ## 
            u <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
            bq <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
            w <- bq != 0
            bq <- bq[w]
            u <- u[w]
            if (length(u) == 0) {
                print_message(paste0("Read label assignment includes no reads for haplotype ", i_hap))
            }
            ## 
            gl <- make_gl_from_u_bq(u, bq, nSNPs)
            ## 
            dosage <- numeric(nSNPs)
            if (i == 1) {
                f <- Rcpp_haploid_dosage_versus_refs_version2
            } else {
                f <- Rcpp_haploid_dosage_versus_refs
            }
            f(
                gl = gl,
                alphaHat_t = alphaHat_t,
                betaHat_t = betaHat_t,
                gamma_t = gamma_t,
                dosage = dosage,
                transMatRate_t = transMatRate_t,
                rhb_t = rhb_t,
                ref_error = ref_error,
                use_eMatDH = use_eMatDH,
                distinctHapsB = distinctHapsB,
                distinctHapsIE = distinctHapsIE,                
                hapMatcher = hapMatcher,
                suppressOutput = 0,
                gammaSmall_t = full_gammaSmall_t1,
                return_dosage = FALSE,
                return_betaHat_t = FALSE,
                return_gamma_t = FALSE,
                return_gammaSmall_t = TRUE,
                gammaSmall_cols_to_get = full_gammaSmall_cols_to_get
            )
            
            if (i == 1) {dosage1 <- dosage}
            if (i == 2) {dosage2 <- dosage}
        }
    }


    expect_equal(dosage1, dosage2)
    print(mean(dosage1))
    print(mean(dosage2))

    if (1 == 0) {

        m <- apply(alphaHat_t, 2, function(x) length(unique(x)))
        m2 <- mclapply(1:(ncol(alphaHat_t) - 1), mc.cores = 4, function(icol) {
            a <- unique(alphaHat_t[, icol])
            b <- unique(alphaHat_t[, icol + 1])
            c <- paste0(match(alphaHat_t[, icol], a), "-", match(alphaHat_t[, icol + 1], b))
            length(unique(c))
            ## sum(table(alphaHat_t[, icol], alphaHat_t[, icol + 1]) != 0)
        })
        m2 <- unlist(m2)

        m3 <- cbind(
            diff = diff(m),
            max = m2,
            m[-1] - m2,
            m2 - m[-length(m)]
        )
        cbind(
            m,
            rbind(m3, c(NA, NA, NA, NA))
        )[1:40, ]
        
        ## AM HERE
        ## HOW MANY ARE CREATED AND DESTROYED
        ## CAN I CONSIDER A REDUCED REPRESENTATION
        ## WHERE I UPDATE THESE ON EVERY GO
        ## AND DO FEWER CALCULATIONS
            ## AM 
        
        ## OK, so how may created / lost per go
        
        ## OK so sometimes several thousands, and sometimes many fewer
        ## seriously can I use lower precision on the values I store...

        ## so often 50+ 
        x <- apply(hapMatcher, 2, function(x) length(unique(x)))
        sort(x)
        mean(x) ## often 80 ish haplotypes
        y <- apply(eMatDH, 2, function(x) length(unique(x)))
        sort(y)
        mean(y)
        mean(y[y != 1]) ## so more like 10 options on average

        ##
        ## use some sweet sweet real data
        ##
        iGrid <- 100
        ## sort(unique(exp(round(log10(alphaHat_t[, iGrid]), 3))))
        prev_alphaHat_col <- alphaHat_t[, iGrid]
        prev_alphaHat_master <- unique(alphaHat_t[, iGrid])
        prev_alphaHat_class <- match(prev_alphaHat_col, prev_alphaHat_master)
        ##
        ## next one
        ## 
        dh_col <- hapMatcher[, iGrid + 1]
        eMatDH_col <- eMatDH[, iGrid + 1]
        prob_col <- eMatDH_col[dh_col]
        ##
        ## original way to do it
        ##
        jump_prob <- full_transMatRate_t_H[2, iGrid] / K
        not_jump_prob <- full_transMatRate_t_H[1, iGrid]
        new_alphaHat_col_old <- (jump_prob + not_jump_prob * prev_alphaHat_col) * prob_col
        ##
        ## alternate way to do it
        ##
        class_compare <- array(0, c(K, nMaxDH + 1))
        to_blank <- array(0, c(K, 2))
        to_reclass <- array(0, c(K, nMaxDH + 1))
        to_reclass_count <- array(0, c(K))
        count <- 1
        for(k in 1:K) {
            if (class_compare[prev_alphaHat_class[k], dh_col[k]] == 0) {
                class_compare[prev_alphaHat_class[k], dh_col[k]] <- 1
                to_blank[count, 1] <- prev_alphaHat_class[k]
                to_blank[count, 2] <- dh_col[k]
                count <- count + 1
            }
        }
        ## now calculate these values only, in these locations
        for(i in 1:(count - 1)) {
            alpha_val <- prev_alphaHat_master[to_blank[i, 1]]
            prob_val <- eMatDH_col[to_blank[i, 2]]
            class_compare[to_blank[i, 1], to_blank[i, 2]] <-
                (jump_prob + not_jump_prob * alpha_val) * prob_val
        }
        ## now calculate new alphaHat column from this
        new_alphaHat_col_new <- numeric(K)
        for(k in 1:K) {
            new_alphaHat_col_new[k] <- class_compare[prev_alphaHat_class[k], dh_col[k]]
        }
        ## now figure out which ones have been used
        for(i in 1:(count - 1)) {
            ac <- to_blank[i, 1]
            dh <- to_blank[i, 2]
            to_reclass_count[ac] <- to_reclass_count[ac] + 1
            to_reclass[ac, to_reclass_count[ac]] <- dh
        }
        ## now re-class those that have been split?

        ## OK, from this, for those split, get split
        
        expect_equal(new_alphaHat_col_old, new_alphaHat_col_new)
        ## now need to re-set class
        for(i in 1:(count - 1)) {
            class_compare[to_blank[i, 1], to_blank[i, 2]] <- 0
        }
        ## OK this is now blank, that is good. count re-sets itself
        expect_equal(sum(class_compare != 0), 0)
        ## now need to update class for next round
        ## 

        
        


        
        
        
    }
    
})

    












test_that("prototype idea of classes for faster alphaHat", {

    skip("not sure this is the right idea")

    
    ## so have initial class
    ## and then have prob, some change
    ## re-make class in efficient way after calc
    set.seed(1)
    K <- 1000
    jump_prob <- 0.01 / K
    not_jump_prob <- 1 - 0.01
    ## previous class
    prev_class <- sample(1:10, K, replace = TRUE)
    prev_alphaHat_master <- runif(10)
    prev_alphaHat_master <- prev_alphaHat_master / 10 ## normalized
    prev_alphaHat_col <- prev_alphaHat_master[prev_class]
    ## now have 4 options for reference haplotype, ignore 0 value for now
    dh_col <- sample(1:4, K, replace = TRUE)
    eMatDH_col <- runif(4) ## and their outputs
    prob_col <- eMatDH_col[dh_col]

    ##
    ## original way to do it
    ##
    new_alphaHat_col <- (jump_prob + not_jump_prob * prev_alphaHat_col) * prob_col

    ##
    ## potential new way to do it
    ##
    ## so almost all the mass in low values, remember, these are ordered





})



test_that("can quickly in cpp get positiosn greater than certain value", {

    ## either an equivalency or not
    K_top_matches <- 10
    
    for(i in 1:1) {
        
        if (i == 1) {
            set.seed(100)
            x <- runif(100) ## so will be distinct
        } else {
            x <- seq(0, 1, length.out = 100)
            x[85:95] <- 0.90 ## so this is the one we need more of
        }

        K <- 100
        alphaHat_t <- array(1, c(K, 2))
        iGrid <- 1
        betaHat_t_col <- array(x, c(K, 1))
        gamma_t_col <- array(0, c(K, 1))

        for(j in 1:2) {
            if (j == 1) { f<- Rcpp_get_top_K_or_more_matches_while_building_gamma }
            if (j == 2) { f<- R_get_top_K_or_more_matches_while_building_gamma }            
            out <- f(
                alphaHat_t = alphaHat_t,
                betaHat_t_col = betaHat_t_col,
                gamma_t_col = gamma_t_col,
                iGrid = iGrid,
                K = K,
                K_top_matches = K_top_matches
            )
            if (j == 1) {out1 <- out}
            if (j == 2) {out2 <- out}
        }

        expect_equal(out1, out2)

    }
    
    ## so first pass, can get value, while going over and building
    ## second pass, can get those that meet value

})



test_that(" can select good haps properly", {

    set.seed(91)
    K <- 10000
    Knew <- 30
    K_top_matches <- 5
    Ksubset <- 400
    previously_selected_haplotypes <- sample(1:K, Ksubset - Knew)

    ## now on each hap, do some, add some overlap though, then study...
    new_haps <- list(
        list(90:95, 123:30, 250:260),
        list(310:315, 350:355, 360:370)
    )

    new_ones <- sort(everything_select_good_haps(
        Knew,
        K_top_matches,
        new_haps,
        previously_selected_haplotypes,
        K
    ))

    ## hm, no real test here, just dummy test, should not exceed 500
    expect_equal(max(new_ones) <= 500, TRUE)

})
