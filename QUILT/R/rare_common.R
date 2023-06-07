make_eHapsCurrent_tc_using_rare_and_common_stuff <- function(
    hapMatcherR,
    distinctHapsIE,
    eMatDH_special_matrix_helper,
    eMatDH_special_matrix,
    rare_per_hap_info,
    snp_is_common,
    which_haps_to_use,
    snp_is_common_1_based,
    Ksubset,
    ref_error
) {
    one_minus_ref_error <- 1 - ref_error
    Ksubset <- length(which_haps_to_use)
    nSNPs <- length(snp_is_common)
    nCommonSNPs <- length(snp_is_common_1_based)
    nCommonGrids <- ceiling(nCommonSNPs / 32)
    ## should be easy in the first instance
    eHapsCurrent_tc <- array(ref_error, c(Ksubset, nSNPs, 1))    
    ## build eHapsCurrent_tc here
    for(iGrid in 1:nCommonGrids) {
        s <- 1 + 32 * (iGrid - 1)
        e <- min(32 * iGrid, nCommonSNPs)
        w <- snp_is_common_1_based[s:e]
        for(i_k in 1:Ksubset) {
            k <- which_haps_to_use[i_k]
            i <- as.integer(hapMatcherR[k, iGrid])
            if (i > 0) {
                ## note, could make distinctHapsB pretty easily if needed for RAM
                eHapsCurrent_tc[i_k, w, 1] <- distinctHapsIE[i, s:e]
            } else {
                b <- rcpp_simple_binary_matrix_search(
                    val = k - 1,
                    mat = eMatDH_special_matrix,
                    s1 = eMatDH_special_matrix_helper[iGrid, 1],
                    e1 = eMatDH_special_matrix_helper[iGrid, 2]
                )
                x <- rcpp_int_expand(b, nSNPs = e - s + 1)
                eHapsCurrent_tc[i_k, w[x == 1], 1] <- one_minus_ref_error
            }
        }
    }
    ## now to rare ones
    for(i_k in 1:length(which_haps_to_use)) {
        k <- which_haps_to_use[i_k]
        eHapsCurrent_tc[i_k, rare_per_hap_info[[k]], 1] <- one_minus_ref_error
    }
    eHapsCurrent_tc
}
