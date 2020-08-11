if (1 == 0) {
    
    dl <- diff(L_grid)
    expRate <- 1
    nGen <- 10
    sigmaCurrent <- exp(-nGen * expRate / 100 / 1000000 * dl)
    transMatRate_t <- rbind(sigmaCurrent, 1 - sigmaCurrent)
    K <- nrow(rhb_t)
    nGrids <- ncol(rhb_t)
    ## 
    alphaHat_t <- array(0, c(K, nGrids))
    betaHat_t <- array(0, c(K, nGrids))
    gamma_t <- array(0, c(K, nGrids))
}


make_gl_from_u_bq <- function(u, bq, nSNPs) {
    gl <- array(1, c(2, nSNPs))
    if (length(u) == 0) {
        return(gl)
    }
    probs <- convertScaledBQtoProbs(matrix(bq, ncol = 1))
    ## 
    for(i in 1:length(u)) {
        gl[, u[i]] <- gl[, u[i]] * probs[i, ]
    }
    return(gl)
}


build_eMatDH <- function(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error) {
    nMaxDH <- nrow(distinctHapsB)
    eMatDH <- array(0, c(nMaxDH, nGrids))
    for(iGrid in 0:(nGrids - 1)) {
        ## 
        s <- 32 * iGrid + 1 ## 1-based start
        e <- min(32 * (iGrid + 1), nSNPs) ## 1-based end
        nSNPsLocal <- e - s + 1
        gl_local <- gl[, s:e, drop = FALSE]        
        ##
        for(k in 0:(nMaxDH - 1)) {
            ref_hapLocal <- STITCH::int_expand(distinctHapsB[k + 1, iGrid + 1], nSNPs = nSNPsLocal)
            ## 
            prob <- 1
            for(b in 0:(nSNPsLocal - 1)) {
                x <- ref_hapLocal[b + 1]
                dR <- gl_local[1, b + 1] ## ref
                dA <- gl_local[2, b + 1] ## alt
                if (x == 0) { ## if ref is 0
                    prob <- prob * (dR * ref_one_minus_error + dA * ref_error)
                } else {
                    prob <- prob * (dR * ref_error + dA * ref_one_minus_error)
                }
            }
            eMatDH[k + 1, iGrid + 1] <- prob
        }
    }
    return(eMatDH)
}


get_prob_for_k <- function(rhb_t, k, iGrid, nSNPsLocal, ref_error, ref_one_minus_error, gl_local) {
    ref_hapLocal <- STITCH::int_expand(rhb_t[k, iGrid], nSNPs = nSNPsLocal)
    ## 
    prob <- 1
    for(b in 0:(nSNPsLocal - 1)) {
        x <- ref_hapLocal[b + 1]
        dR <- gl_local[1, b + 1] ## ref
        dA <- gl_local[2, b + 1] ## alt
        if (x == 0) { ## if ref is 0
            prob <- prob * (dR * ref_one_minus_error + dA * ref_error)
        } else {
            prob <- prob * (dR * ref_error + dA * ref_one_minus_error)
        }
    }
    return(prob)
}


R_haploid_dosage_versus_refs <- function(
    gl,
    alphaHat_t,
    betaHat_t,
    gamma_t,
    dosage,
    transMatRate_t,
    rhb_t,
    ref_error,
    use_eMatDH,
    distinctHapsB,
    distinctHapsIE,
    hapMatcher,
    return_extra = FALSE
) {
    ## run one sample haplotype against potentially very many other haplotypes
    K <- nrow(alphaHat_t)
    nGrids <- ncol(alphaHat_t)
    nSNPs <- ncol(gl)
    one_over_K <- 1 / K
    ref_one_minus_error <- 1 - ref_error
    c <- array(0, nGrids)
    if (use_eMatDH) {
        nMaxDH <- nrow(distinctHapsB)
    } else {
        nMaxDH <- 1
    }
    ##
    ## build emissionGrid container version
    ##
    if (use_eMatDH) {
        eMatDH <- build_eMatDH(
            distinctHapsB = distinctHapsB,
            gl = gl,
            nGrids = nGrids,
            nSNPs = nSNPs,
            ref_error = ref_error,
            ref_one_minus_error = ref_one_minus_error
        )
    } else {
        eMatDH <- array(0, c(1, 1))
    }
    ##
    ## initialize alphaHat_t
    ## 
    iGrid <- 0
    s <- 32 * iGrid + 1 ## 1-based start
    e <- min(32 * (iGrid + 1), nSNPs) ## 1-based end
    nSNPsLocal <- e - s + 1
    gl_local <- gl[, s : e, drop = FALSE]
    for(k in 0:(K - 1)) {
        if (use_eMatDH) {
            dh <- hapMatcher[k + 1, iGrid + 1]
        } else {
            dh <- 0
        }
        if (dh > 0) {
            prob <- eMatDH[dh, iGrid + 1]
        } else {
            prob <- get_prob_for_k(rhb_t, k + 1, iGrid + 1, nSNPsLocal, ref_error, ref_one_minus_error, gl_local)
        }
        alphaHat_t[k + 1, iGrid + 1] <- prob * one_over_K
    }
    c[iGrid + 1] <- 1 / sum(alphaHat_t[, iGrid + 1])
    alphaHat_t[, iGrid + 1] <- alphaHat_t[, iGrid + 1] * c[iGrid + 1]
    ##
    ## run forward algorithm
    ##
    for(iGrid in 1:(nGrids - 1)) {
        ## 
        jump_prob <- transMatRate_t[2, iGrid] / K
        s <- 32 * iGrid + 1 ## 1-based start
        e <- min(32 * (iGrid + 1), nSNPs) ## 1-based end
        nSNPsLocal <- e - s + 1
        gl_local <- gl[, s : e, drop = FALSE]        
        ##
        for(k in 0:(K - 1)) {
            if (use_eMatDH) {
                dh <- hapMatcher[k + 1, iGrid + 1] ## this is 1-based in R
            } else {
                dh <- 0
            }
            if (dh > 0) {
                prob <- eMatDH[dh, iGrid + 1]
            } else {
                prob <- get_prob_for_k(rhb_t, k + 1, iGrid + 1, nSNPsLocal, ref_error, ref_one_minus_error, gl_local)
            }
            alphaHat_t[k + 1, iGrid + 1] <- (jump_prob + (1 - jump_prob) * alphaHat_t[k + 1, iGrid + 1 - 1]) * prob
        }
        ## normalize
        c[iGrid + 1] <- 1 / sum(alphaHat_t[, iGrid + 1])
        alphaHat_t[, iGrid + 1] <- alphaHat_t[, iGrid + 1] * c[iGrid + 1]
    }
    ##
    ## run backward algorithm
    ##
    ematcol <- array(1, K)
    iGrid <- nGrids - 1
    for(iGrid in (nGrids - 1):0) {
        if (iGrid == (nGrids - 1)) {
            betaHat_t_col <- rep(1, K) ## c[nGrids - 1 + 1]
        } else {
            ## this is from 
            jump_prob <- transMatRate_t[2, iGrid + 1] / K
            ## I think we want previous emissions?
            s <- 32 * (iGrid + 1) + 1 ## 1-based start
            e <- min(32 * ((iGrid + 1) + 1), nSNPs) ## 1-based end
            nSNPsLocal <- e - s + 1
            gl_local <- gl[, s : e, drop = FALSE]
            for(k in 0:(K - 1)) {
                if (use_eMatDH) {
                    dh <- hapMatcher[k + 1, iGrid + 1 + 1]
                } else {
                    dh <- 0
                }
                if (dh > 0) {
                    prob <- eMatDH[dh, iGrid + 1 + 1]
                } else {
                    prob <- get_prob_for_k(rhb_t, k + 1, iGrid + 1 + 1, nSNPsLocal, ref_error, ref_one_minus_error, gl_local)
                }
                ematcol[k + 1] <- prob
            }
            ##
            e_times_b <- betaHat_t_col * ematcol
            val <- jump_prob * sum(e_times_b)
            betaHat_t_col <- ((1 - jump_prob) * e_times_b + val)
        }
        ##
        ## now second bit, store of finish off
        ##
        gamma_t_col <- (alphaHat_t[, iGrid + 1] * betaHat_t_col)
        ## do this bit here
        s <- 32 * iGrid + 1 ## 1-based start
        e <- min(32 * (iGrid + 1), nSNPs) ## 1-based end
        nSNPsLocal <- e - s + 1
        dosageL <- array(0, nSNPsLocal)        
        if (use_eMatDH) {
            matched_gammas <- array(0, nMaxDH)
            for(k in 0:(K - 1)) {
                gk <- gamma_t_col[k + 1]
                dh <- hapMatcher[k + 1, iGrid + 1] ## this is 1-based in R
                if (dh > 0) {
                    matched_gammas[dh] <- matched_gammas[dh] + gk
                } else {
                    ref_hapLocal <- STITCH::int_expand(rhb_t[k + 1, iGrid + 1], nSNPs = nSNPsLocal)
                    ref_hapLocal[ref_hapLocal == 0] <- ref_error
                    ref_hapLocal[ref_hapLocal == 1] <- (1 - ref_error)
                    dosageL <- dosageL + gk * ref_hapLocal
                }
            }
            for(b in 0:(nSNPsLocal - 1)) {
                for(dh in 0:(nMaxDH - 1)) {
                    dosageL[b + 1] <- dosageL[b + 1] + matched_gammas[dh + 1] * distinctHapsIE[dh + 1, s + b] ## s is 1-based here
                }
            }
            dosage[s:e] <- dosageL
        } else {
            for(k in 0:(K - 1)) {
                ref_hapLocal <- STITCH::int_expand(rhb_t[k + 1, iGrid + 1], nSNPs = nSNPsLocal)
                ref_hapLocal[ref_hapLocal == 0] <- ref_error
                ref_hapLocal[ref_hapLocal == 1] <- (1 - ref_error)
                dosageL <- dosageL + gamma_t_col[k + 1] * ref_hapLocal
            }
            dosage[s:e] <- dosageL            
        }
        ## finish up
        betaHat_t_col <- betaHat_t_col * c[iGrid + 1]
        if (return_extra) {
            betaHat_t[, iGrid + 1] <- betaHat_t_col
            gamma_t[, iGrid + 1] <- (alphaHat_t[, iGrid + 1] * betaHat_t[, iGrid + 1]) / c[iGrid + 1]
        }
    }
    ##
    ## make gamma
    ##
    ## compare to each other and dosageX
    return(
        list(
            gamma_t = gamma_t,
            dosage = dosage,
            alphaHat_t = alphaHat_t,
            betaHat_t = betaHat_t
        )
    )
}





make_rhb_t_equality <- function(rhb_t, nSNPs, ref_error, nMaxDH) {
    K <- nrow(rhb_t)
    nGrids <- ncol(rhb_t)
    ## --- hapMatcher
    ## matrix K x nGrids
    ## 0 = no match
    ## i is match to ith haplotype in distinctHaps i.e. i
    hapMatcher <- array(0L, c(K, nGrids))
    ## --- distinctHapsB
    ## matrix with nMaxDH x nGrids
    ## matrix with the distinct haplotypes
    distinctHapsB <- array(0, c(nMaxDH, nGrids)) ## store encoded binary
    for(iGrid in 1:nGrids) {
        ## can safely ignore the end, it will be zeros what is not captured
        a <- table(rhb_t[, iGrid], useNA = "always")
        a <- a[order(-a)]
        names_a <- as.integer(names(a))
        w <- names_a[1:min(length(names_a), nMaxDH)]
        distinctHapsB[1:length(w), iGrid] <- w
        ## match against
        hapMatcher[, iGrid] <- as.integer(match(rhb_t[, iGrid], distinctHapsB[, iGrid]))
        hapMatcher[which(is.na(hapMatcher[, iGrid])), iGrid] <- 0L
    }
    ## inflate them too, they're pretty small
    distinctHapsIE <- array(0L, c(nMaxDH, nSNPs)) ## inflated, with ref_error
    for(iGrid in 0:(nGrids - 1)) {
        s <- 32 * iGrid + 1 ## 1-based start
        e <- min(32 * (iGrid + 1), nSNPs) ## 1-based end
        nSNPsLocal <- e - s + 1
        for(k in 1:nMaxDH) {
            distinctHapsIE[k, s:e] <- rcpp_int_expand(distinctHapsB[k, iGrid + 1], nSNPsLocal)
        }
    }
    distinctHapsIE[distinctHapsIE == 0] <- ref_error
    distinctHapsIE[distinctHapsIE == 1] <- 1 - ref_error
    return(
        list(
            distinctHapsB = distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            hapMatcher = hapMatcher
        )
    )
}



make_eMatRead_t_using_binary <- function(
    sampleReads,
    rhb_t,
    nSNPs,
    ref_error,
    language = "Rcpp",
    n = 5000
) {
    nReads <- length(sampleReads)
    K <- nrow(rhb_t)
    eMatRead_t <- array(1, c(K, nReads))
    ## get read start and end too
    bq <- as.numeric(unlist(sapply(sampleReads, function(x) x[[3]])))
    u <- as.integer(unlist(sapply(sampleReads, function(x) x[[4]])))
    ##
    ps <- array(NA, c(2, length(bq)))
    ##
    w <- bq < 0
    eps <- 10 ** (bq[w] / 10)
    ps[1, w] <- 1 - eps
    ps[2, w] <- eps / 3
    w <- bq > 0
    eps <- 10 ** (-bq[w] / 10)
    ps[1, w] <- eps / 3
    ps[2, w] <- (1 - eps)
    ## 
    nr <- sapply(sampleReads, function(x) x[[1]]) + 1
    start <- cumsum(c(1, nr))[-(1 + length(nr))]
    end <- cumsum(nr)
    ## 
    rhb <- t(rhb_t)
    ceil_K_n <- ceiling(K / n)
    if (language == "R") {
        eMatRead_t <- internal_make_eMatRead_t_using_binary(eMatRead_t, rhb, K, nSNPs, u, ps, nReads, start, end, nr, ref_error, ceil_K_n, n)
    } else if (language == "Rcpp") {
        rcpp_internal_make_eMatRead_t_using_binary(eMatRead_t, rhb, K, nSNPs, u, ps, nReads, start, end, nr, ref_error, ceil_K_n, n)
    }
    return(eMatRead_t)
}

internal_make_eMatRead_t_using_binary <- function(eMatRead_t, rhb, K, nSNPs, u, ps, nReads, start, end, nr, ref_error, ceil_K_n, n) {
    ## could get where each read starts and ends?
    ref_one_minus_error <- 1 - ref_error
    for(i in 0:(ceil_K_n - 1)) {
        sh <- n * (i) ## 0-based
        eh <- n * (i + 1) - 1 ## 0-based
        if (eh > (K - 1)) {
            eh <- K - 1
        }
        KL <- eh - sh + 1 ## 1-based, length
        ## so run from s to e, 1-based
        haps  <- inflate_fhb(
            rhb = rhb,
            haps_to_get = sh:eh,
            nSNPs = nSNPs
        )
        for(ik in 0:(KL - 1)) {
            hap <- haps[, ik + 1]
            for(iRead in 0:(nReads - 1)) {
                s <- start[iRead + 1] - 1
                e <- end[iRead + 1] - 1
                prob <- 1
                for(j in 0:(nr[iRead + 1] - 1)) {
                    pR <- ps[1, s + j + 1]
                    pA <- ps[2, s + j + 1]
                    if (hap[u[s + j + 1] + 1] == 0) {
                        prob <- prob * (pR * ref_one_minus_error + pA * ref_error)
                    } else {
                        prob <- prob * (pA * ref_one_minus_error + pR * ref_error)
                    }
                }
                eMatRead_t[(sh + ik) + 1, iRead + 1] <- prob
            }
        }
    }
    return(eMatRead_t)
}


calculate_eMatRead_t_usig_rhb_t <- function(
    sampleReads,
    rhb_t,
    rescale_eMatRead_t = TRUE
) {
    nReads <- length(sampleReads)
    K <- nrow(rhb_t)
    eMatRead_t_full  <- array(1, c(K, nReads))
    s <- 0
    ## slow but whatevs, will fix
    n <- 5000
    eMatRead_t_full  <- array(1, c(K, nReads))    
    for(i in 1:ceiling(K / n)) {
        print(paste0(i, " / ", ceiling(K / n)))
        s <- n * (i - 1) + 1
        e <- min(n * i, K)
        KL <- e - s + 1
        ehc <- array(0, c(KL, nSNPs, 1))
        ## do it in parts
        ehc[, , 1] <- inflate_fhb_t(
            rhb_t,
            haps_to_get = s:e - 1,
            nSNPs
        )
        eMatRead_t_partial  <- array(1, c(KL, nReads))            
        gc(reset = TRUE); gc(reset = TRUE);
        ## 
        rcpp_make_eMatRead_t(
            eMatRead_t = eMatRead_t_partial,
            sampleReads = sampleReads,
            eHapsCurrent_tc = ehc,
            s = 0,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax = 10000,
            eMatHapOri_t = array(0, c(1, 1)), ## ugh                                                      
            pRgivenH1 = array(0),
            pRgivenH2 = array(0),
            run_pseudo_haploid = FALSE,
            prev = 0,
            suppressOutput = 1,
            prev_section = "text",
            next_section = "text",
            rescale_eMatRead_t = rescale_eMatRead_t
        )
        eMatRead_t_full[s:e, ] <- eMatRead_t_partial
        if ((i %% 10) == 0) {
            gc(reset = TRUE);    gc(reset = TRUE);
        }
    }
    gc(reset = TRUE);    gc(reset = TRUE);
    return(eMatRead_t_full)
}