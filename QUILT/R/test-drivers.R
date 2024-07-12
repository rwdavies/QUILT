check_quilt_output <- function(
    file,
    data_package,
    which_snps = NULL,
    tol = 0.1,
    min_info = 0.9,
    max_missingness = 0.1,
    check_info_only = FALSE
) {
    vcf <- read.table(
        file,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    ## check things in turn
    STITCH::check_vcf_info_scores(vcf, min_info)
    if (check_info_only) {
        return(NULL)
    }
    ## check columns
    N <- ncol(vcf) - 9
    nSNPs <- nrow(vcf)
    expect_equal(length(unique(vcf[, 9])), 1)
    truth_phase <- data_package$phase
    per_sample_entries <- strsplit(vcf[1, 9], ":")[[1]]
    for(iSample in 1:N) {
        ## check things in turn
        per_col_vcf <- vcf[, iSample + 9]
        if (!is.null(which_snps)) {
            sample_truth_haps <- truth_phase[which_snps, iSample, ]
        } else {
            sample_truth_haps <- truth_phase[, iSample, ]
        }
        sample_truth_gen <- sample_truth_haps[, 1] + sample_truth_haps[, 2]
        sample_results <- t(sapply(strsplit(per_col_vcf, ":"), I))
        colnames(sample_results) <- per_sample_entries
        ## check not too much missingness i.e.b non-confident data
        expect_true((sum(sample_results[ , "GP"] == "./.") / nSNPs) <= max_missingness)
        ## check genotypes that exist
        gt <- c(NA, NA, 0, 1, 1, 2)[match(sample_results[, "GT"], c("./.", ".|.", "0|0", "0|1", "1|0", "1|1"))]
        ##print(length(gt))
        ##print(length(sample_truth_gen))
        x <- gt != sample_truth_gen
        ## can't all be 0
        if (sum(!is.na(x)) > 0) {
            expect_true(sum(x, na.rm = TRUE) / sum(!is.na(x)) <= tol)
        }
        ## don't continue if all the dosages are also missing
        if (sum(sample_results[, "DS"] != ".") > 0) {
            ##
            ## check genotype posteriors sum to 1 and are between 0 and 1
            ##
            gps <- t(sapply(strsplit(sample_results[, "GP"], ","), as.numeric))
            expect_true(0.998 <= min(rowSums(gps)))
            expect_true(max(rowSums(gps)) <= 1.002)
            expect_true(0 <= min(gps))
            expect_true(max(gps) <= 1)
            ##
            ## check dosages are close to truth
            ##
            expect_true(max(abs(as.numeric(sample_results[, "DS"]) - sample_truth_gen)) < tol)
            ## now check haplotype dosages - hmm, not sure if safe, could be recombs
            ## should ideally be PSE based, but oh well, these are small tests
            observed_haps <- t(sapply(strsplit(unlist(strsplit(sample_results[, "HD"], ":")), ","), I))
            val1 <- max(c(
                max(abs(as.numeric(observed_haps[, 1]) - sample_truth_haps[, 1])),
                max(abs(as.numeric(observed_haps[, 2]) - sample_truth_haps[, 2]))
            ))
            val2 <- max(c(
                max(abs(as.numeric(observed_haps[, 1]) - sample_truth_haps[, 2])),
                max(abs(as.numeric(observed_haps[, 2]) - sample_truth_haps[, 1]))
            ))
            check_value <- (val1 < tol) | (val2 < tol)
            if (!check_value) {
                print("observed then truth")
                print(paste0("val1 = ", val1, ", val2 = ", val2, ", tol = ", tol))
                print(table(observed_haps[, 1], sample_truth_haps[, 1]))
                print(table(observed_haps[, 1], sample_truth_haps[, 2]))
                print(table(observed_haps[, 2], sample_truth_haps[, 1]))
                print(table(observed_haps[, 2], sample_truth_haps[, 2]))
                print(data.frame(observed_haps, sample_truth_haps))
                print(val1)
                print(val2)
            }
            expect_true(check_value)
        }
    }
    return(NULL)
}


## here, we check phasing is exact
check_sew_phase <- function(file, phase, which_snps = NULL) {
    vcf <- read.table(
        file,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    ## GT is first column
    gt <- sapply(strsplit(vcf[, 10], ":"), function(x) x[1])
    gt1 <- substr(gt, 1, 1)
    gt2 <- substr(gt, 3, 3)
    ##
    if (is.null(which_snps)) {
        t1 <- phase[, 1, 1]
        t2 <- phase[, 1, 2]
    } else {
        t1 <- phase[which_snps, 1, 1]
        t2 <- phase[which_snps, 1, 2]
    }
    ## basically, need
    hap1_check <- (sum(gt1 != t1) == 0) | (sum(gt2 != t1) == 0)
    hap2_check <- (sum(gt1 != t2) == 0) | (sum(gt2 != t2) == 0)
    if (!hap1_check | !hap2_check) {
        print(paste0("phasing results are:", gt))
        print(paste0("truth is:", apply(phase[, 1, ], 1, paste, collapse = "|")))
    }
    expect_equal(hap1_check, TRUE)
    expect_equal(hap2_check, TRUE)
}





#' @export
make_quilt_fb_test_package <- function(
    K = 4,
    nReads = 8,
    nSNPs = 10,
    S = 2,
    gridWindowSize = 3,
    method = "diploid",
    ff = NA,
    seed = 4916,
    eHapsMin = 0.1,
    return_eMatGridTri_t = TRUE,
    L = NULL,
    bq_mult = 30,
    randomize_sample_read_length = FALSE,
    simple_ematread = FALSE
) {
    if (method == "diploid") {
        n_haps <- 2
        hap_probs <- c(0.5, 0.5)
    } else if (method == "triploid-nipt") {
        n_haps <- 3
        if (is.na(ff)) {
            stop("make ff not NA")
        }
        hap_probs <- c(0.5, (1 - ff) / 2, ff / 2)
    } else {
        stop("this only supports diploid (general?) and triploid-nipt")
    }
    set.seed(seed)
    if (is.null(L)) {
        L <- 1:nSNPs
    }
    maxDifferenceBetweenReads <- 1e10
    maxEmissionMatrixDifference <- 1e50
    Jmax_local <- 1e3
    ##
    out <- assign_positions_to_grid(
        L = L,
        gridWindowSize = gridWindowSize
    )
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids
    snps_in_grid_1_based <- out$snps_in_grid_1_based
    ## make almost the same
    eHapsCurrent_tc <- array(NA, c(K, nSNPs, S))
    alphaMatCurrent_tc <- array(NA, c(K, nGrids - 1, S))
    ##
    for(s in 1:S) {
        eHapsCurrent_tc[, , s] <-
            0.1 * array(runif(K * nSNPs), c(K, nSNPs)) +
            0.9 * array(sample(c(0, 1), K * nSNPs, replace = TRUE), c(K, nSNPs))
        if ((method == "triploid-nipt") & (K >= 3)) {
            eHapsCurrent_tc[1, ,s] <- rep(c(eHapsMin, 1 - eHapsMin), each = 1, len = nSNPs)
            eHapsCurrent_tc[2, ,s] <- rep(c(eHapsMin, 1 - eHapsMin), each = 2, len = nSNPs)
            eHapsCurrent_tc[3, ,s] <- (1 - eHapsMin)
        }
        if (simple_ematread && gridWindowSize == 32) {
            for(iGrid in 1:nGrids) {
                w <- 1:32 + 32 * (iGrid - 1)
                x <- rpois(1, 3) + 1 ## how many distinct haps
                y <- matrix(sample(c(0, 1), x * 32, prob = c(0.90, 0.10), replace = TRUE), x, 32)
                prob <- c(1, rep(0.1, x - 1))                
                w2 <- sample(1:x, K, replace =  TRUE, prob = prob / sum(prob))
                eHapsCurrent_tc[, w, s] <- y[w2, ]
            }
        }
        ## m <- array(runif(K * (nGrids - 1)), c(K, (nGrids - 1)))
        alphaMatCurrent_tc[, , s] <- 1 / K
    }
    sigmaCurrent_m <- array(0.9 + 0.1 * runif((nGrids - 1) * S), c(nGrids - 1, S))
    priorCurrent_m <- array(1 / K, c(K, S))
    ##
    transMatRate_tc_D <- get_transMatRate_m("diploid", sigmaCurrent_m)
    transMatRate_tc_H <- get_transMatRate_m("diploid-inbred", sigmaCurrent_m)
    ## sample reads
    cr <- sort(sample(1:nGrids, nReads, replace = TRUE)) - 1
    true_H <- sample(1:n_haps, nReads, replace = TRUE, prob = hap_probs)
    ## choose them from haps 1 and 2
    sampleReads <- lapply(
        1:length(cr),
        function(ii) {
        i <- cr[ii] ## 0-based
        w <- which(grid == i)
        g <- w[sample(1:length(w), size = 1)] - 1
        x <- c(g - 2, g - 1, g, g + 1, g + 2)
        x <- x[x %in% 0:(nSNPs - 1)] ## 0-based
        if (randomize_sample_read_length) {
            ## choose two points at random, make them start and end
            a <- sort(sample(1:length(x), 2, replace = TRUE))
            x <- x[a[1]:a[2]]
        }
        i_hap <- true_H[ii]
        bq <- matrix(round(bq_mult * (eHapsCurrent_tc[i_hap, x + 1, 1] - 0.5)), ncol = 1)
        list(
            length(x) - 1,
            i,
            bq,
            matrix(x,ncol=1,nrow=length(x))
        )
    })
    N <- 1
    ## almost certainly want to pre-declare
    list_of_eMatRead_t <- lapply(0:(S - 1), function(s) {
        eMatRead_t <- array(1, c(K, nReads))
        rcpp_make_eMatRead_t(
            eMatRead_t = eMatRead_t,
            sampleReads = sampleReads,
            eHapsCurrent_tc = eHapsCurrent_tc,
            s = s,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax = Jmax_local,
            eMatHapOri_t = array(0, c(1, 1)), ## ugh
            pRgivenH1 = array(0),
            pRgivenH2 = array(0),
            run_pseudo_haploid = FALSE,
            prev = 0,
            suppressOutput = 1,
            prev_section = "text",
            next_section = "text"
        )
        return(eMatRead_t)
    })
    ##
    list_of_eMatGridTri_t <- NULL
    if (method == "triploid-nipt" && return_eMatGridTri_t) {
        list_of_eMatGridTri_t <- lapply(1:S, function(s) {
            eMatGridTri_t <- rcpp_make_and_bound_eMatGridTri_t(
                eMatRead_t = list_of_eMatRead_t[[s]],
                ff = ff,
                sampleReads = sampleReads,
                nReads = nReads,
                K = K,
                nGrids = nGrids,
                maxEmissionMatrixDifference = maxEmissionMatrixDifference,
                run_fb_grid_offset = 0
            )
        })
    }
    ##
    list_of_eMatGrid_t <- lapply(0:(S - 1), function(s) {
        eMatGrid_t <- array(1, c(K, nGrids))
        rcpp_make_eMatGrid_t(
            eMatGrid_t = eMatGrid_t,
            eMatRead_t = list_of_eMatRead_t[[s + 1]],
            H = 1,
            sampleReads = sampleReads,
            hap = 1,
            nGrids = nGrids,
            run_fb_grid_offset = 0,
            use_all_reads = TRUE,
            bound = TRUE,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            prev = 0,
            suppressOutput = 1,
            prev_section = "text",
            next_section = "text"
        )
        return(eMatGrid_t)
    })
    ##
    return(
        list(
            eHapsCurrent_tc = eHapsCurrent_tc,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            sigmaCurrent_m = sigmaCurrent_m,
            priorCurrent_m = priorCurrent_m,
            transMatRate_tc_D = transMatRate_tc_D,
            transMatRate_tc_H = transMatRate_tc_H,
            S = S,
            K = K,
            gridWindowSize = gridWindowSize,
            nReads = nReads,
            nSNPs = nSNPs,
            L = L,
            N = N,
            grid = grid,
            grid_distances = grid_distances,
            L_grid = L_grid,
            nGrids = nGrids,
            snps_in_grid_1_based = snps_in_grid_1_based,
            sampleReads = sampleReads,
            true_H = true_H,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax_local = Jmax_local,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            list_of_eMatRead_t = list_of_eMatRead_t,
            list_of_eMatGrid_t = list_of_eMatGrid_t,
            list_of_eMatGridTri_t = list_of_eMatGridTri_t
        )
    )
}



#' @export
make_reference_single_test_package <- function(
    nSNPs = 500,
    K = 1000,
    seed = 4916,
    L = NULL,
    expRate = 100,
    nGen = 10,
    nMaxDH = 2 ** 8 - 1,
    ref_error = 0.01,
    gammaSmall_cols_to_get = NULL
) {
    set.seed(seed)
    if (is.null(L)) {
        L <- sort(sample(1:(nSNPs * 10), nSNPs))
    }
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
    dl <- diff(L_grid)
    sigmaCurrent <- exp(-nGen * expRate / 100 / 1000000 * dl)
    transMatRate_t <- rbind(sigmaCurrent, 1 - sigmaCurrent)
    ## in each grid, make a few, then a tail
    reference_haps <- array(0, c(nSNPs, K))
    for(iGrid in 1:(nGrids - 1)) {
        ##
        ## this is designed to be semi-realistic
        ## choose some number of haplotypes to sample from (y)
        ## then choose some smaller number of haplotypes to make random (z)
        ## then y are made from z with some small alterations
        ## like deep branches vs short branches
        ##
        y <- rpois(n = K, lambda = 10)
        nLocal <- 4
        if (iGrid == 3 | iGrid == 10) {
            ##
            y[sample(1:K, nMaxDH + 20, replace = FALSE)] <- 1:(nMaxDH + 20)
            ##
            nLocal <- 6
        }
        n <- max(y) + 1
        y2 <- array(NA, c(32, n))
        ## here, not that deep
        z <- max(rpois(n = 1, lambda = 6), nLocal)
        ## choose each one with a different frequency
        for(iSNP in 1:32) {
            if (iGrid == 3 | iGrid == 10) {
                af <- rbeta(n = 1, 1, 1)
            } else {
                af <- rbeta(n = 1, 0.1, 0.1)
            }
            y2[iSNP, 1:z] <- sample(c(0, 1), z, replace = TRUE, prob = c(1 - af, af))
        }
        ##
        for(i in (z + 1):n) {
            y2[, i] <- y2[, sample(1:z, 1)]
            y2[sample(1:32, nLocal), i] <- sample(c(0, 1), nLocal, replace = TRUE)
        }
        x <- y2[, y + 1]
        reference_haps[32 * (iGrid - 1) + 1:32, ] <- x ## rows are SNPs
    }
    rhi <- reference_haps
    rhi_t <- t(rhi)
    rhb_t <- STITCH::make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)
    K <- nrow(rhb_t)
    nGrids <- ncol(rhb_t)
    out <- STITCH::make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nSNPs,
        ref_error = ref_error,
        use_hapMatcherR = FALSE
    )
    distinctHapsB <- out[["distinctHapsB"]]
    distinctHapsIE <- out[["distinctHapsIE"]]
    hapMatcher <- out[["hapMatcher"]]
    eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
    eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]
    eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
    eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]
    ## get hapMatcherR as well
    out <- STITCH::make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nSNPs,
        ref_error = ref_error,
        use_hapMatcherR = TRUE
    )
    hapMatcherR <- out[["hapMatcherR"]]
    my_hap <- c(
        rhi_t[1, 1:36],
        rhi_t[2, 37:60],
        rhi_t[3, 61:250],
        rhi_t[4, 257:400],
        rhi_t[5, 401:nSNPs]
    )
    ## make about 0.25X coverage
    nReads <- round(nSNPs * 0.25)
    u <- sort(sample(1:nSNPs, nReads, replace = TRUE))
    ## though specifically remove some regions, make have no variants
    u <- u[!(u %in% 200:300)]
    ##
    bq <- rep(-10, length(u))
    bq[my_hap[u] == 1] <- 10
    gl <- make_gl_from_u_bq(u, bq, nSNPs)
    ## choose some of these
    if (is.null(gammaSmall_cols_to_get)) {
        gammaSmall_cols_to_get <- integer(nGrids)
        gammaSmall_cols_to_get[] <- -1L
        w <- seq(2, nGrids, 4)
        gammaSmall_cols_to_get[w] <- as.integer(0:(length(w) - 1))
    }
    return(
        list(
            distinctHapsB = distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            rhb_t = rhb_t,
            gl = gl,
            transMatRate_t = transMatRate_t,
            ref_error = ref_error,
            eMatDH_special_values_list = eMatDH_special_values_list,
            eMatDH_special_grid_which = eMatDH_special_grid_which,
            eMatDH_special_matrix = eMatDH_special_matrix,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            gammaSmall_cols_to_get = gammaSmall_cols_to_get,
            truth_hap = my_hap
        )
    )
}

