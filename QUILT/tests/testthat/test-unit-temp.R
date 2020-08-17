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

test_that("profile", {

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

    


