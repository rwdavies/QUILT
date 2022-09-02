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
        key <- nGrids * mtm[, 3] + mtm[, 4]
        length <- mtm[, 4] - mtm[, 3]
        mtm <- cbind(mtm, key, length)
        mtm
    })
    mtm <- rbind(a[[1]], a[[2]])
    mtm <- mtm[order(-mtm[, 6], mtm[, 5]), ]
    unique_keys <- unique(mtm[, "key"])
    if (length(unique_keys) > Knew) {
        return(mtm[1:Knew, 2])
    } else if(nrow(mtm) <= Knew)  {
        new_haps <- array(NA, Knew)
        vals <- b[, 2]
        new_haps[1:length(vals)] <- vals
        new_haps[-c(1:length(vals))] <- sample(setdiff(1:Kfull, vals), Knew - length(vals), replace = FALSE)
        return(new_haps)
    } else {
        ## now need to do something smart
        ## there are too many to take all
        ## but too few unique to just only choose those
        ## so we order to be first unique, then second unique, etc
        mtm <- cbind(mtm, ikey = 0)
        for(i in 2:nrow(mtm)) {
            if (mtm[i, "key"] == mtm[i - 1, "key"]) {
                mtm[i, "ikey"] <- mtm[i - 1, "ikey"] + 1
            }
        }
        mtm <- mtm[order(mtm[, "ikey"], mtm[, "length"], mtm[, "key"]), ]
        return(mtm[1:Knew, 2])
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
