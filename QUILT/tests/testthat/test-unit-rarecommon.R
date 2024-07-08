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

}



test_that("can make old eHapsCurrent_tc using rare common idea", {

    set.seed(1011010)
    
    ref_error <- 0.001
    maxDifferenceBetweenReads <- 1e6

    nReads <- 800
    nSNPs <- 2000
    nCommonSNPs <- 100
    K <- 1000
    Ksubset <- 50
    snp_is_common <- rep(FALSE, nSNPs)
    snp_is_common[sample(nSNPs, nCommonSNPs)] <- TRUE
    common_snp_index <- integer(nSNPs)
    common_snp_index[which(snp_is_common)] <- 1:nCommonSNPs ## 1-based
    

    rhi <- array(0L, c(nSNPs, K))
    rhi[snp_is_common, ] <- sample(c(0L, 1L), K * nCommonSNPs, prob = c(0.7, 0.3), replace = TRUE)
    rhi[!snp_is_common, ] <- sample(c(0L, 1L), K * (nSNPs - nCommonSNPs), prob = c(0.99, 0.01), replace = TRUE)
    rhi_t <- t(rhi)

    ## now, make objects
    nGrids <- ceiling(nSNPs / 32)
    grid <- floor((1:nSNPs) / 32)
    nCommonGrids <- ceiling(nCommonSNPs / 32)
    rhb <- array(0L, c(nCommonGrids, K))
    for(k in 1:K) {
        rhb[, k] <- rcpp_int_contract(rhi[snp_is_common, k])
    }
    rhb_t <- t(rhb)


    ## make special objects
    out <- STITCH::make_rhb_t_equality(
        rhb_t = rhb_t,
        nSNPs = nCommonSNPs,
        ref_error = ref_error,
        use_hapMatcherR = TRUE,
        verbose = FALSE
    )
    
    hapMatcherR <- out[["hapMatcherR"]]
    eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]
    eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
    distinctHapsIE <- out[["distinctHapsIE"]]

    ## manually make special rare one too,
    snp_is_rare_1_based <- which(!snp_is_common)    
    rare_per_hap_info <- sapply(1:K, function(k) {
        ## for now store position among all SNPs
        snp_is_rare_1_based[which(rhi[!snp_is_common, k] == 1)]
    })


    which_haps_to_use <- sort(sample(1:K, Ksubset))
    
    old_eHapsCurrent_tc <- array(0, c(Ksubset, nSNPs, 1))
    old_eHapsCurrent_tc[, , 1][rhi_t[which_haps_to_use, ] == 0] <- ref_error
    old_eHapsCurrent_tc[, , 1][rhi_t[which_haps_to_use, ] == 1] <- (1 - ref_error)

    snp_is_common_1_based <- which(snp_is_common)

    new_eHapsCurrent_tc <- make_eHapsCurrent_tc_using_rare_and_common_stuff(
        hapMatcher = array(0, c(1, 1)),
        hapMatcherR = hapMatcherR,
        distinctHapsIE = distinctHapsIE,
        use_hapMatcherR = TRUE,
        eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
        eMatDH_special_matrix = eMatDH_special_matrix,
        rare_per_hap_info = rare_per_hap_info,
        snp_is_common = snp_is_common,
        which_haps_to_use = which_haps_to_use,
        snp_is_common_1_based = snp_is_common_1_based,
        Ksubset = Ksubset,
        ref_error = ref_error
    )

    expect_equal(new_eHapsCurrent_tc, old_eHapsCurrent_tc)

    ## ## get some fake-ish sample reads
    ## ## honestly they can be very fake here
    sampleReads <- lapply(1:nReads, function(ii) {
        rl <- sample(1:5, 1) ## length
        s <- sample(1:(nSNPs - rl), 1) ## start
        e <- s + rl - 1 ## end
        cg <- grid[round(mean(c(s, e)))] ## central grid
        bq <- matrix(sample(c(-25:-20, 26:30), rl, replace = TRUE), ncol = 1)
        u <- matrix(s:e - 1, ncol = 1)
        list(
            rl - 1,
            cg,
            bq,
            u
        )
    })
    sampleReads <- sampleReads[order(sapply(sampleReads, function(x) x[[2]]))]

    ## generate using eHapsCurrent_tc
    original_eMatRead_t <- array(1, c(Ksubset, nReads))
    rcpp_make_eMatRead_t(
        eMatRead_t = original_eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = new_eHapsCurrent_tc,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        s = 0,
        Jmax = 100,
        eMatHapOri_t = array(0, c(1, 1)),
        pRgivenH1 = array(0, 1),
        pRgivenH2 = array(0, 1),
        prev = 0,
        suppressOutput = 1,
        prev_section = "wer",
        next_section = "wer",
        rescale_eMatRead_t = TRUE,
        run_pseudo_haploid = FALSE
    )

    original_eMatRead_t[c(12, 15, 42, 50), 1] ## OK so maybe

    round(new_eHapsCurrent_tc[, 16:19, 1] )
    
    sampleReads[[1]][[3]]
    
    
    ## now generate using special objects (can't!)
    ## new_eMatRead_t <- array(1, c(Ksubset, nReads))
    ## print("wer1")
    ## Rcpp_make_eMatRead_t_for_gibbs_using_objects(
    ##     sampleReads = sampleReads,
    ##     eMatRead_t = new_eMatRead_t,
    ##     hapMatcher = array(0, c(1, 1)),
    ##     hapMatcherR = hapMatcherR,
    ##     use_hapMatcherR = TRUE,
    ##     grid = grid,
    ##     rhb_t = array(as.raw(0), c(1, 1)),
    ##     distinctHapsIE = distinctHapsIE,
    ##     eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
    ##     eMatDH_special_matrix = eMatDH_special_matrix,
    ##     ref_error = ref_error,
    ##     which_haps_to_use = which_haps_to_use,
    ##     Jmax = 100,
    ##     maxDifferenceBetweenReads = maxDifferenceBetweenReads,
    ##     use_eMatDH_special_symbols = TRUE,
    ##     rescale_eMatRead_t = TRUE
    ## );
    ## print("wer2")    

    ## expect_equal(original_eMatRead_t, new_eMatRead_t)

    ## OK, now from here, verify I can quickly figure out which reads I need to interrogate

    check_eMatRead_t <- determine_which_reads_require_k_haps(
        sampleReads = sampleReads,
        nSNPs = nSNPs,
        rare_per_hap_info = rare_per_hap_info,
        which_haps_to_use = which_haps_to_use
    )

    ## do slower check
    slow_check_eMatRead_t <- array(1, c(Ksubset, nReads))
    rphiL <- rare_per_hap_info[which_haps_to_use]    
    x <- unlist(sapply(1:length(rphiL), function(k) {
        rep(k, length(rphiL[[k]]))
    }))
    y <- unlist(rphiL)
    for(iRead in 1:nReads) {
        u <- sampleReads[[iRead]][[4]] + 1 ## 1-based
        a <- !is.na(match(y, u))
        expect_one <- rep(TRUE, Ksubset)
        if (sum(a) > 0) {
            expect_one[x[a]] <- FALSE
        }
        slow_check_eMatRead_t[, iRead] <- expect_one
    }

    expect_equal(slow_check_eMatRead_t, check_eMatRead_t)

    ## hopefully not a RAM monster
    ## create a new one that is per-site, and gives list of k
    rare_per_snp_info <- lapply(1:nSNPs, function(x) -1L)
    ## want this to be the index of the rare SNP in the overall
    snp_is_rare_1_based <- which(!snp_is_common)            
    for(k in 1:length(which_haps_to_use)) {
        ## expand back out to all SNPs
        snps <- rare_per_hap_info[[which_haps_to_use[[k]]]] ## this IS an index among all SNPs
        ## this is among all SNPs
        for(snp in snps) {
            rare_per_snp_info[[snp]] <- c(rare_per_snp_info[[snp]], k) ## keep 1-based
        }
    }
    
    ##
    ## now use this one and expand it 
    ##
    aabbaa <- check_eMatRead_t[, 13]
    final_eMatRead_t <- check_eMatRead_t
    final_eMatRead_t[] <- 1
    Rcpp_make_eMatRead_t_for_final_rare_common_gibbs_using_objects(
        eMatRead_t = final_eMatRead_t,
        rare_per_hap_info = rare_per_hap_info,
        common_snp_index = common_snp_index,
        sampleReads = sampleReads,
        snp_is_common = snp_is_common,
        hapMatcherR = hapMatcherR,
        grid = grid,
        distinctHapsIE = distinctHapsIE,
        eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
        eMatDH_special_matrix = eMatDH_special_matrix,
        ref_error = ref_error,
        which_haps_to_use = which_haps_to_use,
        rescale_eMatRead_t = TRUE,        
        Jmax = 100,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        rare_per_snp_info = rare_per_snp_info
    )

    
    expect_equal(final_eMatRead_t, original_eMatRead_t)

    ## print("ori values")
    ## 
    ## print(w)
    ## ## wk1 <- c(7, 23, 26, 50)    
    ## wk1 <- 5:8
    ## wk1 <- 21:26
    ## ## print(rare_per_hap_info[which_haps_to_use[wk1]])
    ## ## Maybe check out 6, 7, 8?
    ## print(m[wk1, ])
    ## print(m[c(7, 23, 26, 50), ])
    ## ## print(m)
    ## u <- sampleReads[[13]][[4]] + 1
    ## ##print(sampleReads[[13]][[3]]) ## base qualities, it is ref ref
    ## ##print(snp_is_common[u])     ## the second one is common
    ## print("new eHaps values")
    ## print(round(new_eHapsCurrent_tc[wk1, u, 1]))
    
    ## ## print(distinctHapsIE[as.integer(hapMatcherR[which_haps_to_use[wk1], 1]), 3]    )
    ## print("what hapMatcher values should be")
    ## print(sapply(as.integer(hapMatcherR[which_haps_to_use[wk1], 1]),
    ##        function(a) {
    ##     if (a == 0) {
    ##         return(NA)
    ##     }
    ##     distinctHapsIE[a, 3]
    ## }))

    ## ## print(as.integer(hapMatcherR[which_haps_to_use[wk1], 1]))
    
    
})



if (1 == 0) {

    sampleReads <- allSNP_sampleReads
    nReads <- length(sampleReads)
    nSNPs <- nrow(pos_all)
    
    ## 1 is normal
    ## 0 is special
    eMatRead_t <- array(1, c(Ksubset, nReads))
    
    ## rare_per_hap_info is 1-based
    ## want knowledge of for each hap, 
    rphiL <- rare_per_hap_info[which_haps_to_use]

    ## build an index of what sites have a k
    has_a_k <- logical(nSNPs)
    for(k in 1:Ksubset) {
        has_a_k[rphiL[[k]]] <- TRUE
    }
    n_affected_SNPs <- sum(has_a_k)
    ## now get their positions
    snp_n_for_has_a_k <- integer(nSNPs)
    snp_n_for_has_a_k[has_a_k] <- 1:n_affected_SNPs
    
    ## among those sites, figure out which k are involved
    which_k <- vector("list", n_affected_SNPs)
    for(k in 1:Ksubset) {
        for(n in snp_n_for_has_a_k[rphiL[[k]]]) {
            which_k[[n]] <- c(which_k[[n]], k)
        }
    }
    
    ## now find reads
    read_has_a_k <- logical(nReads)
    for(iRead in 1:nReads) {
        u <- sampleReads[[iRead]][[4]] + 1 ## 1-based
        if (sum(has_a_k[u]) > 0) {
            ## get these (could be more than 1)
            x <- snp_n_for_has_a_k[u]
            y <- x[x > 0]
            for(yy in y) {
                ks <- which_k[[yy]]
                eMatRead_t[ks, iRead] <- 0
            }
        }
    }
    
    ## verify this here by doing everything manually
    x <- unlist(sapply(1:length(rphiL), function(k) {
        rep(k, length(rphiL[[k]]))
    }))
    y <- unlist(rphiL)
    
    for(iRead in 1:nReads) {
        u <- sampleReads[[iRead]][[4]] + 1 ## 1-based
        a <- !is.na(match(y, u))
        expect_one <- rep(TRUE, Ksubset)
        if (sum(a) > 0) {
            expect_one[x[a]] <- FALSE
        }
        stopifnot(length(expect_one) == Ksubset)
        stopifnot(eMatRead_t[, iRead] == expect_one)
        which(!expect_one)
        which(eMatRead_t[, iRead] == 0)
    }

    ## AWW YEAH BABY
    ## ROCK FROM THIS
    
}



