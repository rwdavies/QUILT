select_new_haps_internal <- function(hap)  {
    Z <- round(hap)
    Z[mono0] <- 0
    Z[mono1] <- 1
    ##
    out <- MatchZ_Algorithm5_Rcpp(
        X = fhb_t,
        a = indices[["a"]],
        u = indices[["u"]],
        v = indices[["v"]],
        c = indices[["c"]],
        d = indices[["d"]],
        Z = Z
    )
    rownames(out) <- 1:nrow(out)
    out2 <- out[rownames(unique(out[, 3:4])), ]
    out2 <- data.frame(out2)
    out2$diff <- out2[, 4] - out2[, 3]
    out <- data.frame(out)
    out$diff <- out[, 4] - out[, 3]
    return(list(out, out2))
}


get_more_new_good_haps <- function(a, n, great_candidates) {
    a <- a[is.na(match(a[, 2] + 1, great_candidates)), , drop = FALSE]
    if (nrow(a) == 0) {
        return(NULL)
    } else if (nrow(a) == 1) {
        return(a[, 2] + 1)
    }
    sample(x = a[, 2] + 1, size = min(n, nrow(a)), prob = a[, 5] / sum(a[, 5]))
}



select_new_haps_mspbwt <- function(
    hapProbs_t,
    hapMatcher,
    ms_indices,
    Knew,
    Kfull,
    all_symbols,
    nGrids
) {
    a <- lapply(1:2, function(ihap) {
        hap <- round(hapProbs_t[ihap, ])
        Zs <- rcpp_int_contract(hap)
        Zg <- mspbwt::map_Z_to_all_symbols(Zs, all_symbols)
        mtm <- mspbwt::Rcpp_ms_MatchZ_Algorithm5(
            X = hapMatcher,
            ms_indices = ms_indices,
            Z = Zg
        )
        mtm[, 2] <- mtm[, 2] + 1 ## make 1-based here
        mtm
    })
    ## choose based on length and frequency
    b <- rbind(a[[1]], a[[2]])
    key <- nGrids * b[, 3] + b[, 4]
    c <- match(key, unique(nGrids * b[, 3] + b[, 4]))
    weight <- 1 / table(c)[c]
    length <- b[, 4] - b[, 3]
    desire <- weight * length
    b <- b[order(-desire), ]
    vals <- unique(b[, 2])
    if (length(vals) >= Knew) {
        new_haps <- vals[1:Knew]
    } else {
        new_haps <- array(NA, Knew)
        new_haps[1:length(vals)] <- vals
        new_haps[-c(1:length(vals))] <- sample(setdiff(1:Kfull, vals), Knew - length(vals), replace = FALSE)
    }
    new_haps
}



select_new_haps_mspbwt_v2 <- function(
    hapProbs_t,
    hapMatcher,
    ms_indices,
    Knew,
    Kfull,
    all_symbols,
    nGrids
) {
    nIndices <- 3
    mtms <- lapply(1:nIndices, function(iIndex) {
        a <- lapply(1:2, function(ihap) {
            hap <- round(hapProbs_t[ihap, ])
            Zs <- rcpp_int_contract(hap) 
            w <- seq(iIndex, ncol(hapMatcher), nIndices)  
            Zg <- mspbwt::map_Z_to_all_symbols(Zs[w], ms_indices[[iIndex]][["all_symbols"]])
            mtm <- mspbwt::Rcpp_ms_MatchZ_Algorithm5(
              X = hapMatcher[, w],
              ms_indices = ms_indices[[iIndex]],
              Z = Zg
            )
            mtm[, 2] <- mtm[, 2] + 1 ## make 1-based here
            key <- nGrids * mtm[, 3] + mtm[, 4]
            length <- mtm[, 4] - mtm[, 3]
            mtm <- cbind(mtm, key, length)
            mtm
        })
        mtm <- rbind(a[[1]], a[[2]])
        mtm <- mtm[order(-mtm[, 6], mtm[, 5]), ]
    })
    ## all that now matters is haplotype being considered, and length
    mtm <- rbind(mtms[[1]], mtms[[2]], mtms[[3]])
    unique_haps <- unique(mtm[, 2])
    if (length(unique_haps) > Knew) {
        return(mtm[1:Knew, 2])
    } else if(length(unique_haps) <= Knew)  {
        new_haps <- array(NA, Knew)
        unique_haps <- mtm[, 2]
        new_haps[1:length(unique_haps)] <- unique_haps
        new_haps[-c(1:length(unique_haps))] <- sample(setdiff(1:Kfull, unique_haps), Knew - length(unique_haps), replace = FALSE)
        return(new_haps)
    } else {
        ## sort by length
        ## though darn this does nothing about region specificity!
        ## though not sure that is a problem
        mtm <- mtm[order(mtm[, "length"], mtm[, 2]), ]
        return(unique(mtm[, 2])[1:Knew])
    }
}


##             print("--------------wer----------------saving stuff--------------wer--------------")
##                 hapProbs_t = gibbs_iterate[["hapProbs_t"]]                
##                 hapMatcher = hapMatcher
##                 ms_indices = ms_indices
##                 Knew = Knew
##                 Kfull = nrow(rhb_t)
##             all_symbols = ms_indices[["all_symbols"]]
##             dir <- "/well/davies/users/dcc832/QUILT_nicola_testing_2022_08_26/"                
##             if (i_gibbs_sample == 1 & i_it == 1) {
##             save(
## hapProbs_t,
## hapMatcher,
## ms_indices,
## Knew,
## Kfull,
## all_symbols,
## nGrids  ,
## file = paste0(dir, "stuff.RData"))
##             }
##                 nGrids = nGrids
##             save(
## hapProbs_t,
## file = paste0(dir, "stuff.", i_gibbs_sample, ".", i_it, ".RData"))
##             print("--------------wer----------------done saving stuff--------------wer--------------")    







impute_using_split_reads_and_small_ref_panel <- function(
    H,
    which_haps_to_use,
    sampleReads,
    rhb_t,
    nSNPs,
    full_alphaHat_t,
    full_betaHat_t,
    full_gamma_t,
    full_gammaSmall_t,
    full_gammaSmall_cols_to_get,
    full_transMatRate_t_H,
    distinctHapsB,
    distinctHapsIE,
    hapMatcher,
    eMatDH_special_grid_which,
    eMatDH_special_values_list,
    ref_error,
    make_plots,
    outplotprefix,
    have_truth_haplotypes,
    truth_haps,
    truth_labels,
    uncertain_truth_labels,
    L_grid,
    L,
    inRegion2,
    cM_grid,
    ancAlleleFreqAll,
    return_good_haps,
    plot_description,
    Knew,
    previously_selected_haplotypes,
    sample_name,
    smooth_cm,
    regionStart,
    regionEnd,
    buffer,
    minGLValue,
    return_dosage = FALSE,
    return_betaHat_t = FALSE,
    return_gamma_t = FALSE,
    K_top_matches = 5,
    heuristic_match_thin = 0.01,
    suppressOutput = 1
    ) {






    
    K <- nrow(rhb_t)
    dosage <- numeric(nSNPs)
    nGrids <- ncol(rhb_t)
    ##
    ## for now - can be ruthless - don't worry about speed
    ##
    full_alphaHat_t <- array(0, c(length(which_haps_to_use), nGrids))
    
    hapMatcherL <- hapMatcher[which_haps_to_use, ]
    rhb_tL <- rhb_t[which_haps_to_use, ]


## save(
##     H,
##     which_haps_to_use,
##     sampleReads,
##     rhb_tL,
##     nSNPs,
##     full_alphaHat_t,
##     full_betaHat_t,
##     full_gamma_t,
##     full_gammaSmall_t,
##     full_gammaSmall_cols_to_get,
##     full_transMatRate_t_H,
##     distinctHapsB,
##     distinctHapsIE,
##     hapMatcherL,
##     eMatDH_special_grid_which,
##     eMatDH_special_values_list,
##     ref_error,
##     make_plots,
##     outplotprefix,
##     have_truth_haplotypes,
##     truth_haps,
##     truth_labels,
##     uncertain_truth_labels,
##     L_grid,
##     L,
##     inRegion2,
##     cM_grid,
##     ancAlleleFreqAll,
##     return_good_haps,
##     plot_description,
##     Knew,
##     previously_selected_haplotypes,
##     sample_name,
##     smooth_cm,
##     regionStart,
##     regionEnd,
##     buffer,
##     minGLValue,
##     return_dosage,
##     return_betaHat_t,
##     return_gamma_t,
##     K_top_matches ,
##     heuristic_match_thin ,
##     suppressOutput,
##     file = "~/temp.RData")

##     load("~/temp.RData")
    ##
    ##
    ##
    w <- which(eMatDH_special_grid_which > 0)
    if (length(w) > 0) {
        for(iw in 1:length(w)) {
            ww <- w[iw]
            a <- sum((which_haps_to_use - 1) %in% eMatDH_special_values_list[[iw]])
            if (a == 0) {
                eMatDH_special_grid_which[ww] <- 0
            } else {
                eMatDH_special_values_list[[iw]] <- which(hapMatcherL[, ww] == 0) - 1L
            }
        }
    }
    
    ##
    ##
    ## 
    return_gammaSmall_t <- FALSE
    get_best_haps_from_thinned_sites <- FALSE
    best_haps_stuff_list <- list()
    return_dosage <- TRUE
    return_betaHat_t <- FALSE
    return_gamma_t = FALSE
    K_top_matches = 5
    heuristic_match_thin = 0.01
    suppressOutput = 1
    ##
    ##
    ##
    for(i_hap in 1:2) {
        ##
        u <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
        bq <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
        w <- bq != 0
        bq <- bq[w]
        u <- u[w]
        if (length(u) == 0) {
            print_message(paste0("Read label assignment includes no reads for haplotype ", i_hap))
        }
        gl <- make_gl_from_u_bq(u, bq, nSNPs, minGLValue = minGLValue)
        use_eMatDH <- TRUE
        c <-  array(1, c(nGrids))  ## more useful for debugging
        Rcpp_haploid_dosage_versus_refs(
            gl = gl,
            alphaHat_t = full_alphaHat_t,
            betaHat_t = full_betaHat_t,
            c = c,
            gamma_t = full_gamma_t,
            dosage = dosage,
            transMatRate_t = full_transMatRate_t_H,
            rhb_t = rhb_tL,
            ref_error = ref_error,
            use_eMatDH = use_eMatDH,
            distinctHapsB = distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            hapMatcher = hapMatcherL,
            eMatDH_special_grid_which = eMatDH_special_grid_which,
            eMatDH_special_values_list = eMatDH_special_values_list,
            suppressOutput = suppressOutput,
            return_dosage = return_dosage,
            return_betaHat_t = return_betaHat_t,
            return_gamma_t = return_gamma_t,
            return_gammaSmall_t = return_gammaSmall_t,
            gammaSmall_t = full_gammaSmall_t,
            gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
            get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
            best_haps_stuff_list = best_haps_stuff_list,
            K_top_matches = K_top_matches,
            always_normalize = FALSE,
            is_version_2 = TRUE,
            normalize_emissions = TRUE
        )
        ## do some checks here
        if ((min(dosage) < -1e-5) | (1 + 1e-5) < max(dosage)) {
            print(range(dosage))
            stop("Dosage observed outside of range of 0 to 1 on forward-backward full iteration. Something has gone wrong. Please report this")
        }
        dosageNew <- numeric(nSNPs)
        dosageNew <- dosage[]
        if (i_hap == 1) { dosage1 <- dosageNew}
        if (i_hap == 2) { dosage2 <- dosageNew}
    }
    
    to_return <- list(
        dosage1 = dosage1,
        dosage2 = dosage2
    )
    return(to_return)
}

