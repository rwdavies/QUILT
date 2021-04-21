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




temp_compare_two_versions <- function(
            rhb_t,
            small_eHapsCurrent_tc,
            which_haps_to_use,
            nSNPs,
            ref_error,
            maxDifferenceBetweenReads,
            Jmax,
            hapMatcher,
            grid,
            distinctHapsIE,
            sampleReads
) {

    K <- length(which_haps_to_use)
    rescale_eMatRead_t <- TRUE
    nReads <- length(sampleReads)
    
    f1 <- function() {
        ##
        S <- 1
        ## argh
        ## print(paste0("start = ", Sys.time()))
        inflate_fhb_t_in_place(
            rhb_t,
            small_eHapsCurrent_tc,
            haps_to_get = which_haps_to_use - 1,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
        eMatRead_t_old <- array(1, c(K, nReads))        
        rcpp_make_eMatRead_t(
            eMatRead_t = eMatRead_t_old,
            sampleReads = sampleReads,
            eHapsCurrent_tc = small_eHapsCurrent_tc,
            s = 0,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax = Jmax,
            eMatHapOri_t = matrix(0, nrow = 0, ncol = 0),
            pRgivenH1 = vector(),
            pRgivenH2 = vector(),
            prev = 1,
            suppressOutput = 1,
            prev_section = "N/A",
            next_section = "N/A",
            run_pseudo_haploid = FALSE,
            rescale_eMatRead_t = rescale_eMatRead_t
        )
        eMatRead_t_old
    }


    f2 <- function() {
        eMatRead_t_new <- array(1, c(K, nReads))
        Rcpp_make_eMatRead_t_for_gibbs_using_objects(
            eMatRead_t = eMatRead_t_new,
            sampleReads = sampleReads,
            hapMatcher = hapMatcher,
            grid = grid,
            rhb_t = rhb_t,
            distinctHapsIE = distinctHapsIE,
            ref_error = ref_error,
            which_haps_to_use = which_haps_to_use,
            rescale_eMatRead_t = rescale_eMatRead_t,
            Jmax = Jmax,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads
        )
        eMatRead_t_new
    }

    library("testthat")
    expect_equal(f1(), f2())

    library("microbenchmark")
    ## can be slower (on tall data, HRC sized)
    ## so can be faster (on wide data)
    print(microbenchmark(f1(), f2(), times = 2))
    ## so def slower. hmm

}
