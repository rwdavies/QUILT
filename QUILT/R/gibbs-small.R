R_make_eMatRead_t_for_gibbs_using_objects <- function(
    sampleReads,
    hapMatcher,
    grid,
    rhb_t,
    distinctHapsIE,
    ref_error,
    which_haps_to_use,
    normalize = TRUE,
    do_checks = FALSE,
    rhi_t = NULL
) {
    ##
    Ksmall <- length(which_haps_to_use)
    nReads <- length(sampleReads)
    ## row = haplotypes, column = grid
    eMatRead_t <- array(1, c(Ksmall, nReads))
    for(iRead in 1:nReads) {
        sampleRead <- sampleReads[[iRead]]
        J <- sampleRead[[1]]
        bq <- sampleRead[[3]]
        u <- sampleRead[[4]]
        iGrid0 <- grid[u[1] + 1]
        haps_at_grid <- hapMatcher[which_haps_to_use, iGrid0 + 1]
        iGrid0_prev <- iGrid0
        for(j in 1:(J + 1)) {
            if(bq[j] < 0) {
                eps <- 10 ** (bq[j] / 10)
                pR <- 1 - eps
                pA <- eps / 3
            } else {
                eps <- 10 ** (-bq[j] / 10)
                pR <- eps / 3
                pA <- 1 - eps
            }
            jj <- u[j]
            iGrid0 <- grid[u[j] + 1]
            if (iGrid0 != iGrid0_prev) {
                haps_at_grid <- hapMatcher[which_haps_to_use, iGrid0 + 1]
            }
            iGrid0_prev <- iGrid0
            for(k in 1:Ksmall) {
                if (haps_at_grid[k] > 0) {
                    e <- distinctHapsIE[haps_at_grid[k], u[j] + 1]
                } else {
                    ## OK, I need to get this one - do fast way first? rare?
                    val <- rcpp_int_expand(rhb_t[which_haps_to_use[k], iGrid0 + 1], 32)[u[j] - iGrid0 * 32 + 1]
                    ## check it?
                    if (do_checks) {
                        if (rhi_t[which_haps_to_use[k], u[j] + 1] != val) {
                            stop("bad!")
                        }
                    }
                    if (val == 1) {
                        e <- 1 - ref_error
                    } else {
                        e <- ref_error
                    }
                }
                eMatRead_t[k, iRead] <- eMatRead_t[k, iRead] * (e * pA + (1 - e) * pR)
            }
        }
        if (normalize) {
            eMatRead_t[, iRead] <- eMatRead_t[, iRead] / max(eMatRead_t[, iRead])
        }
    }
    eMatRead_t
}




