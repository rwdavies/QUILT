check_quilt_output <- function(
    file,
    data_package,
    which_snps = NULL,
    tol = 0.1,
    min_info = 0.9,
    max_missingness = 0.1
) {
    vcf <- read.table(
        file,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    ## check things in turn
    STITCH::check_vcf_info_scores(vcf, min_info)
    ## check columns
    N <- ncol(vcf) - 9
    nSNPs <- nrow(vcf)
    truth_phase <- data_package$phase
    for(iSample in 1:N) {
        ## check things in turn
        per_col_vcf <- vcf[, iSample + 9]
        sample_truth_haps <- truth_phase[which_snps, iSample, ]
        sample_truth_gen <- sample_truth_haps[, 1] + sample_truth_haps[, 2]
        sample_results <- t(sapply(strsplit(per_col_vcf, ":"), I))
        ## check not too much missingness i.e. non-confident data
        expect_true((sum(sample_results[ ,1] == "./.") / nSNPs) < max_missingness)
        ## check genotypes that exist
        gt <- c(NA, 0, 1, 2)[match(sample_results[, 1], c("./.", "0/0", "0/1", "1/1"))]
        x <- gt != sample_truth_gen
        expect_true(sum(x, na.rm = TRUE) / sum(!is.na(x)) < tol)
        ## now check dosages
        expect_true(max(abs(as.numeric(sample_results[, 3]) - sample_truth_gen)) < tol)
        ## now check haplotype dosages - hmm, not sure if safe, could be recombs
        ## should ideally be PSE based, but oh well, these are small tests
        observed_haps <- t(sapply(strsplit(unlist(strsplit(sample_results[, 4], ":")), ","), I))
        val1 <- max(c(
            max(abs(as.numeric(observed_haps[, 1]) - sample_truth_haps[, 1])),
            max(abs(as.numeric(observed_haps[, 2]) - sample_truth_haps[, 2]))
        ))
        val2 <- max(c(
            max(abs(as.numeric(observed_haps[, 1]) - sample_truth_haps[, 2])),
            max(abs(as.numeric(observed_haps[, 2]) - sample_truth_haps[, 1]))
        ))
        expect_true((val1 < tol) | (val2 < tol))
    }
    return(NULL)
}


## here, we check phasing is exact
check_sew_phase <- function(vcf, phase, which_snps = NULL) {
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
    randomize_sample_read_length = FALSE
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
        m <- array(runif(K * (nGrids - 1)), c(K, (nGrids - 1)))
        alphaMatCurrent_tc[, , s] <- t(t(m) / colSums(m))
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

