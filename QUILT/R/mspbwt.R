build_mspbwt_indices <- function(
   hapMatcher,
   hapMatcherR,
   mspbwt_nindices,
   use_hapMatcherR,
   all_symbols,
   use_list_of_columns_of_A
) {
    if (use_hapMatcherR) {
        ncol <- ncol(hapMatcherR)
        nrow <- nrow(hapMatcherR)        
    } else {
        ncol <- ncol(hapMatcher)
        nrow <- nrow(hapMatcher)
    }
    ## auto-choose?
    if (nrow > 1e5) {
        egs <- 1000
    } else {
        egs <- 100
    }
    ms_indices <- lapply(1:mspbwt_nindices, function(iIndex) {
        w <- seq(iIndex, ncol, mspbwt_nindices)
        if (use_hapMatcherR) {
            X1C <- hapMatcherR[, w, drop = FALSE]
        } else {
            X1C <- hapMatcher[, w, drop = FALSE]
        }
        out <- mspbwt::Rcpp_ms_BuildIndices_Algorithm5(
            X1C = X1C,
            all_symbols = all_symbols[w],
            indices = list(),
            verbose = FALSE,
            egs = egs
        )
        ## turn off d
        out[["d"]] <- matrix(1L, 1, 1)
        if (use_list_of_columns_of_A) {
            ## nuke most of a
            a <- out[["a"]]
            list_of_columns_of_A <- as.list(1:ncol(a))
            if (ncol(a) < 10) {
                cols <- 1:ncol(a) 
            } else {
                cols <- sort(unique(c(1, 3, 5, seq(1, ncol(a), 10))))
            }
            for(i in cols) {
                list_of_columns_of_A[[i]] <- a[, i]
            }
            out[["list_of_columns_of_A"]] <- list_of_columns_of_A
            out[["a"]] <- NULL
        }
        return(out)
    })
    ms_indices
}

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
    hapMatcherR,
    use_hapMatcherR,
    ms_indices,
    Knew,
    Kfull
) {
    iIndex <- 1
    ihap <- 1
    nIndices <- length(ms_indices)
    if (use_hapMatcherR) {
        nGrids <- ncol(hapMatcherR)
    } else {
        nGrids <- ncol(hapMatcher)
    }
    a <- lapply(1:2, function(ihap) {
        hap <- round(hapProbs_t[ihap, ])
        Zs <- rcpp_int_contract(hap)
        mtms <- lapply(1:nIndices, function(iIndex) {
            which_grids <- seq(iIndex, nGrids, nIndices)
            Z_local <- mspbwt::map_Z_to_all_symbols(Zs[which_grids], ms_indices[[iIndex]][["all_symbols"]])
            mtm <- mspbwt::Rcpp_ms_MatchZ_Algorithm5(
                X = hapMatcher,
                XR = hapMatcherR,
                use_XR = use_hapMatcherR,
                ms_indices = ms_indices[[iIndex]],
                Z = Z_local,
                cols_to_use0 = as.integer(which_grids - 1L),
                use_cols_to_use0 = TRUE,
                verbose = FALSE,
                min_length = 3
            )
            mtm[, 2] <- mtm[, 2] + 1 ## make 1-based here
            key <- nGrids * mtm[, 3] + mtm[, 4]
            length <- mtm[, 4] - mtm[, 3]
            mtm <- cbind(mtm, key, length)
            mtm
        })
        mtm <- mtms[[1]]
        if (length(mtms) > 1) {
            for(j in 2:length(mtms)) {
                mtm <- rbind(mtm, mtms[[j]])
            }
        }
        mtm
    })
    mtm <- rbind(a[[1]], a[[2]])
    ## order everything
    mtm <- mtm[order(-mtm[, "length"], mtm[, "key"]), , drop = FALSE]
    unique_haps <- unique(mtm[, "indexB0"])
    if (length(unique_haps) == 0) {
        ## special fluke case
        ## likely driven by very small regions we are trying to impute, with few / no matches above the min length above
        new_haps <- sample(1:Kfull, Knew)
        return(new_haps)
    } else if(length(unique_haps) <= Knew)  {
        ##new_haps <- array(NA, Knew)
        ## new_haps[1:length(unique_haps)] <- unique_haps
        ## new_haps[-c(1:length(unique_haps))] <- sample(setdiff(1:Kfull, unique_haps), Knew - length(unique_haps), replace = FALSE)
        ##
        ## so in this faster version
        ## oversample all we could possibly want. then take the new ones, plus some needed new ones
        ##
        new_haps <- unique(c(
            unique_haps,
            sample(Kfull, length(unique_haps) + Knew, replace = FALSE)
        ))[1:Knew]
        return(new_haps)
    } else {
        ## so this doesn't do anything about region specificity
        ## as long as there are buffers it should be pretty OK
        ## it will favour the longest matches, and take one per key
        ## if it exhausts that, it will take other unique long hones
        unique_keys <- unique(mtm[, "key"])
        unique_haps_at_unique_keys <- unique(mtm[match(unique_keys, mtm[, "key"]), "indexB0"])
        if (length(unique_haps_at_unique_keys) > Knew) {
            return(unique_haps_at_unique_keys[1:Knew])
        } else {
            ## otherwise, take unique ones, then next best ones, from length down
            return(c(unique_haps_at_unique_keys, setdiff(unique_haps, unique_haps_at_unique_keys))[1:Knew])
        }
    }
}









select_new_haps_mspbwt_v3 <- function(
    hapProbs_t,
    hapMatcher,
    hapMatcherR,
    use_hapMatcherR,
    ms_indices,
    Knew,
    Kfull,
    mspbwtL,
    mspbwtM,
    heuristic_approach,
    method
) {    

    ## save(
    ## hapProbs_t,
    ## hapMatcher,
    ## hapMatcherR,
    ## use_hapMatcherR,
    ## ms_indices,
    ## Knew,
    ## Kfull,
    ## mspbwtL,
    ## mspbwtM,
    ## heuristic_approach,
    ## method,
    ## file = "~/temp.RData")
    ## stop("WER - werwer")
    
    iIndex <- 1
    ihap <- 1
    nIndices <- length(ms_indices)
    if (use_hapMatcherR) {
        nGrids <- ncol(hapMatcherR)
        K <- nrow(hapMatcherR)
    } else {
        nGrids <- ncol(hapMatcher)
        K <- nrow(hapMatcher)
    }
    if (!is.null(ms_indices[[1]]$list_of_columns_of_A)) {
        use_list_of_columns_of_A <- TRUE
    } else {
        use_list_of_columns_of_A <- FALSE
    }
    
    out <- lapply(1:nrow(hapProbs_t), function(ihap) {
        
        hap <- round(hapProbs_t[ihap, ])
        Zs <- rcpp_int_contract(hap)
        
        mtms <- lapply(1:nIndices, function(iIndex) {
            
            ## print_message(paste0("ihap = ", ihap, ", iIndex = ", iIndex))
            which_grids <- seq(iIndex, nGrids, nIndices)
            Z_local <- mspbwt::map_Z_to_all_symbols(Zs[which_grids], ms_indices[[iIndex]][["all_symbols"]])
            if (use_hapMatcherR) {
                X <- matrix(0, 1, 1)
                XR <- hapMatcherR
            } else {
                X <- hapMatcher
                XR <- matrix(0, 1, 1)
            }
            
            if (heuristic_approach == "A") {

                ## call this one approach A
                ## message(paste0("ihap = ", ihap, ", iIndex = ", iIndex, ", ", date()))
                mtm <- mspbwt::Rcpp_find_good_matches_without_a(
                    Z = Z_local,
                    all_symbols = ms_indices[[iIndex]][["all_symbols"]],
                    usge_all = ms_indices[[iIndex]][["usge_all"]],
                    egs = ms_indices[[iIndex]][["egs"]],
                    pbwtL = mspbwtL,
                    pbwtM = mspbwtM,
                    hapMatcherR = hapMatcherR,
                    which_snps_in_hapMatcherR = which_grids,
                    verbose = FALSE,
                    list_of_columns_of_A = ms_indices[[iIndex]][["list_of_columns_of_A"]],
                    use_list_of_columns_of_A = use_list_of_columns_of_A,
                    K = K
                )
                
                ## I think this is right
                colnames(mtm) <- c("start0", "index0", "len1")
                mtm <- cbind(mtm[, "index0"], mtm[, "start0"] + 1, mtm[, "start0"] + mtm[, "len1"], mtm[, "len1"])
                colnames(mtm) <- c("index0", "start1", "end1", "len1")
                
            } else {

                if (is.null(ms_indices[[iIndex]][["a"]])) {
                    stop("You need to set use_list_of_columns_of_A to FALSE to use this option")
                }

                mtm <- mspbwt::Rcpp_ms_MatchZ_Algorithm5(
                    X = X,
                    XR = XR,
                    use_XR = use_hapMatcherR,
                    ms_indices = ms_indices[[iIndex]],
                    Z = Z_local,
                    mspbwtM = mspbwtM,
                    mspbwtL = mspbwtL,
                    do_up_and_down_scan = TRUE,
                    cols_to_use0 = as.integer(which_grids - 1L),
                    use_cols_to_use0 = TRUE,
                    verbose = FALSE,
                    have_d = FALSE,
                    cap_scan_count = max(100L, mspbwtL) ## don't bother doing a crazy number
                )
                
                
            }
            ## )[["uppy_downy_reporter"]]
            ## change to 1-based
            mtm[, "index0"] <- mtm[, "index0"] + 1
            colnames(mtm)[colnames(mtm) == "index0"] <- "index1"
            ## order so the same index and start comes first, with longest first
            mtm <- mtm[order(mtm[, 1], -mtm[, "end1"], -mtm[, "start1"]), ]            
            x <- c(FALSE, diff(mtm[, "index1"]) == 0 & diff(mtm[, "start1"]) == 0)
            ## y <- !duplicated(paste0(mtm[, "index1"], "-", mtm[, "start1"]))

            ##print("in-A")
            ##save(mtm, x, file = "~/temp.RData")
            mtm <- mtm[!x, ]
            ##print("out-A")
            
            ## 
            key <- nGrids * mtm[, "start1"] + mtm[, "end1"]
            length <- mtm[, "len1"]
            mtm <- cbind(mtm, key, n = iIndex) ## , length)
            mtm
        })
        mtm <- mtms[[1]]
        if (length(mtms) > 1) {
            for(j in 2:length(mtms)) {
                mtm <- rbind(mtm, mtms[[j]])
            }
        }
        ## 
        ## mtm <- mtm[order(mtm[, 1], -mtm[, "end1"], -mtm[, "start1"]), ]
        ## order by length
        mtm <- mtm[order(-mtm[, "len1"]), ]
        mtm
    })
    ## check max number
    if (method == "diploid") {
        unique_haps <- unique(c(out[[1]][, 1], out[[2]][, 1]))
    } else {
        unique_haps <- unique(c(out[[1]][, 1], out[[2]][, 1], out[[3]][, 1]))
    }
    if (length(unique_haps) == 0) {
        ## special fluke case
        ## likely driven by very small regions we are trying to impute, with few / no matches above the min length above
        new_haps <- sample(1:Kfull, Knew)
        return(new_haps)
    } else if (length(unique_haps) <= Knew)  {
        ##
        ## so in this faster version
        ## oversample all we could possibly want. then take the new ones, plus some needed new ones
        ##
        new_haps <- unique(c(
            unique_haps,
            sample(Kfull, length(unique_haps) + Knew, replace = FALSE)
        ))[1:Knew]
        return(new_haps)
    } else {
        ##
        ## heuristically, prioritize based on length and new-ness
        ## do this for each of the two haps
        ##
        results <- lapply(out, function(mtm) {
            ## m <- max(mtm[, "end1"])
            weight <- numeric(nrow(mtm))
            cur_sum <- numeric(max(mtm[, "end1"]))
            cur_sum[] <- 1
            for(i in 1:nrow(mtm)) {
                s <- mtm[i, "start1"]
                e <- mtm[i, "end1"]
                weight[i] <- (e - s + 1) * 1 / sum(cur_sum[s:e])
                cur_sum[s:e] <- cur_sum[s:e] + 1
            }
            ##
            o <- order(-weight)
            if (is.null(dim(mtm))) {
                save(o, mtm, weight, out, file = "~/debug.RData")
                stop("there is an error because mtm has no dimension")
            }
            if (length(weight) != nrow(mtm)) {
                save(o, mtm, weight, out, file = "~/debug.RData")
                stop("there is an error because weight and mtm have different sizes, please see ~/debug.RData")
            }
            if (min(o) < 1 | max(o) > length(weight)) {
                save(o, mtm, weight, out, file = "~/debug.RData")
                stop("there is an error because order of weight went wrong, please see ~/debug.RData")
            }
            mtm <- mtm[o, ]
            mtm[, "index1"]
        })
        ## pad out one of them
        x <- results[[1]]
        y <- results[[2]]
        if (method == "nipt") {
            z <- results[[3]]
        } else {
            z <- NULL
        }
        a <- max(c(length(x), length(y), length(z)))
        x <- c(x, rep(NA, a - length(x)))
        y <- c(y, rep(NA, a - length(y)))
        if (method == "nipt") {
            z <- c(z, rep(NA, a - length(z)))
        }
        if (method == "diploid") {
            unique_ordered_haps <- unique(c(t(cbind(x, y))))
        } else {
            unique_ordered_haps <- unique(c(t(cbind(x, y, z))))
        }
        unique_ordered_haps <- unique_ordered_haps[!is.na(unique_ordered_haps)]
        ## want to interleave, then make unique, then take first Knew
        ## print(paste("select", length(unique_ordered_haps), " unique haps after post-selection 2"))
        ## 
        if (length(unique_ordered_haps) >= Knew) {
            return(unique_ordered_haps[1:Knew])
        } else {
            ## add in some other (potentially) duplicated haps
            new_haps <- c(setdiff(unique_ordered_haps, unique_haps), unique_haps)[1:Knew]
            print(paste("select", length(unique(new_haps)), " unique haps after post-selection 3"))
            return(new_haps)
        }
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
        if (use_eigen) {
            arma_alphaHat_t <- array(0, c(1, 1))
            eigen_alphaHat_t <- full_alphaHat_t
        } else {
            arma_alphaHat_t <- full_alphaHat_t
            eigen_alphaHat_t <- array(0, c(1, 1))
        }
        Rcpp_haploid_dosage_versus_refs(
            gl = gl,
            arma_alphaHat_t = arma_alphaHat_t,
            eigen_alphaHat_t = eigen_alphaHat_t,
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

