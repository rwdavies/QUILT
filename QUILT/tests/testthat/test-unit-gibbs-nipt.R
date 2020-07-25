if (1 == 0) {

    library("testthat"); library("STITCH"); library("rrbgen")
    ##    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    ##    dir <- "~/Google Drive/STITCH/"
    dir <- "~/proj/STITCH-private/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))


}




test_that("read set probabilities agree", {

    H <- sample(c(1, 2, 3), 100, replace = TRUE, prob = c(0.5, 0.25, 0.25))
    ff <- 0.25
    rc <- c(sum(H == 1), sum(H == 2), sum(H == 3))
    expect_equal(
        calc_prob_of_set_of_reads(ff, rc),
        rcpp_calc_prob_of_set_of_reads(ff, rc)
    )

    ## expect to be similar!
    exp(calc_prob_of_set_of_reads(ff, rc))
    n <- 100
    rc <- c(sum(H == 1), sum(H == 2), sum(H == 3))
    factorial(n) / factorial(rc[1]) /factorial(rc[2]) / factorial(rc[3]) *
        0.5 ** rc[1] * 0.25 ** rc[2] *  0.25 ** rc[3]
    dmultinom(rc, prob = c(0.5, 0.25, 0.25))

    ## off by 2 orders or magnitude?
    exp(calc_prob_of_set_of_reads(ff, rc = table(H)))

    ff <- 0.20
    H1 <- c(rep(1, 65), rep(2, 27), rep(3, 8))
    H2 <- c(rep(1, 28), rep(2, 37), rep(3, 35))
    calc_prob_of_set_of_reads(ff, rc = table(H1))
    calc_prob_of_set_of_reads(ff, rc = table(H2))

})



test_that("can re-label read labels where entirely changing is much better", {

    ## so e.g. we have our labelling of 1, 2, 3
    ## but compared to truth, mt = 3, mu = 2, p = 1
    ## given f, we should be able to switch that easily
    ff <- 0.2
    H <- c(rep(1, 34), rep(2, 61), rep(3, 5))
    rc <- c(34, 61, 5)

    expect_equal(
        get_weights_for_entire_relabelling(rc, ff),
        rcpp_get_weights_for_entire_relabelling(rc, ff)
    )

    set.seed(72)
    outR <- consider_and_try_entire_relabelling(H, ff)

    H_for_cpp <- vector("integer", length(H))
    H_for_cpp[] <- as.integer(H[])

    set.seed(72)
    outRcpp <- rcpp_consider_and_try_entire_relabelling(H_for_cpp, ff)

    ## print(rbind(H, outR[["H"]], H_for_cpp))
    expect_equal(outR[["relabel"]], outRcpp)
    expect_equal(as.integer(outR[["H"]]), as.integer(H_for_cpp))

    ## test that can apply relabel appropriately
    rcORI <- round(100 * c(0.5, 0.5 - ff / 2, ff / 2))

    for(i in 1:6) {

        ## start with an H
        HORI <- c(rep(1, rcORI[1]), rep(2, rcORI[2]), rep(3, rcORI[3]))
        ## choose one of the 6 permutations
        if (i == 1) {reorder <- c(1, 2, 3)}
        if (i == 2) {reorder <- c(1, 3, 2)}
        if (i == 3) {reorder <- c(2, 1, 3)}
        if (i == 4) {reorder <- c(2, 3, 1)}
        if (i == 5) {reorder <- c(3, 1, 2)}
        if (i == 6) {reorder <- c(3, 2, 1)}
        ## i.e. this particular re-ordering has been applied
        H <- reorder[HORI]
        ## force this to happen as well
        out <- consider_and_try_entire_relabelling(H, ff, take_max = TRUE)
        relabel <- out$relabel
        expect_equal(out$H, HORI)
        H2 <- H
        H2[1] <- H2[1] + 1; H2[1] <- H2[1] - 1
        H2 <- as.integer(H2)
        rcpp_consider_and_try_entire_relabelling(H2, ff, relabel = relabel)
        expect_equal(H2, HORI)
        ##
        mt <- as.numeric(names(table(H))[table(H) == rcORI[1]])
        mu <- as.numeric(names(table(H))[table(H) == rcORI[2]])
        p <- as.numeric(names(table(H))[table(H) == rcORI[3]])
        ##m1o <- matrix(1, 1, 1); m2o <- matrix(2, 1, 1); m3o <- matrix(3, 1, 1)
        ##m1 <- matrix(mt, 1, 1); m2 <- matrix(mu, 1, 1); m3 <- matrix(p, 1, 1)
        m1o <- array(mt, c(2, 2)); m2o <- array(mu, c(2, 2)); m3o <- array(p, c(2, 2))
        ##
        m1 <- array(1, c(2, 2)); m2 <- array(2, c(2, 2)); m3 <- array(3, c(2, 2))
        m1v <- array(1, 2); m2v <- array(2, 2); m3v <- array(3, 2)
        ##
        out <- apply_relabel(m1, m2, m3, relabel)
        expect_equal(out$m1, m1o)
        expect_equal(out$m2, m2o)
        expect_equal(out$m3, m3o)
        ## c++
        rcpp_apply_mat_relabel(m1, m2, m3, relabel)
        expect_equal(m1, m1o)
        expect_equal(m2, m2o)
        expect_equal(m3, m3o)
        ## vec
        rcpp_apply_vec_relabel(m1v, m2v, m3v, relabel)
        expect_equal(m1v, array(mt, 2))
        expect_equal(m2v, array(mu, 2))
        expect_equal(m3v, array(p, 2))

    }

    H <- c(rep(1, 45), rep(2, 25), rep(3, 30))
    rc <- c(45, 25, 30)
    ff <- 0.10538
    expect_equal(
        get_weights_for_entire_relabelling(rc, ff),
        rcpp_get_weights_for_entire_relabelling(rc, ff)
    )

    set.seed(72)
    outR <- consider_and_try_entire_relabelling(H, ff)
    H_for_cpp <- vector("integer", length(H))
    H_for_cpp[] <- as.integer(H[])

    set.seed(72)
    outRcpp <- rcpp_consider_and_try_entire_relabelling(H_for_cpp, ff)
    expect_equal(outR[["relabel"]], outRcpp)
    expect_equal(as.integer(outR[["H"]]), as.integer(H_for_cpp))

})


test_that("gibbs-nipt normalization scheme makes sense", {

    ## so idea is that say the likelihhods are, in e space,
    for(i_l in 1:13) {
        set.seed(i_l * 100)
        log_mult_max <- 10
        equal_weighting <- FALSE
        if (i_l == 1) { ll <- c(-47, -43)}
        if (i_l == 2) { ll <- c(-47, -43, -42, -46, -50)}
        if (i_l == 3) { ll <- c(-100, -40, -42)}
        if (i_l == 4) { ll <- c(-100, -100, -100 + log_mult_max - 1, -100 + log_mult_max + 1, -100 + log_mult_max -2)}
        if (i_l == 5) { ll <- c(-100, -100, -100 + log_mult_max + 1, -100 + log_mult_max - 1, -100 + log_mult_max +2)}
        if (i_l == 6) { ll <- c(-100, -100, -100 + log_mult_max - 1, -100 + log_mult_max + 1, -100 + 2 * log_mult_max - 2, -100 + 2 * log_mult_max - 2, -100 + 2 * log_mult_max + 2)}
        if (i_l == 7) { ll <- c(-100, -10, -100, 10) }
        if (i_l == 8) { ll <- c(-10, -100, -10, 100) }
        if (i_l == 9) { ll <- -(ceiling(log_mult_max * 1.1) * (10:1))}
        if (i_l == 10) { ll <- c(-100, -100, -100 + log_mult_max - 1, -100 + log_mult_max + 1, -100 + 2 * log_mult_max - 2, -100 + 2 * log_mult_max - 2, -100 + 2 * log_mult_max - 2)}
        if (i_l == 11) { log_mult_max <- 2; ll <- rep(-100, 100)}
        if (i_l == 12) { log_mult_max <- 40; equal_weighting = TRUE; ll <- c(-100, -100, -100 + log_mult_max - 1, -100 + log_mult_max + 1, -100 + 2 * log_mult_max - 2, -100 + 2 * log_mult_max - 2, -100 + 2 * log_mult_max + 2)}
        if (i_l == 13) { log_mult_max <- 40; ll <- c(-327.4957, -282.7341, -315.1111, -367.5477, -433.4326, -321.4748, -401.8183, -361.8846, -385.7613, 247.5394) }
        ll <- ll - max(ll) ## for now!
        x <- exp(ll)
        weights <- x / sum(x)
        if (equal_weighting) {
            weights <- rep(1 / length(x), length(x))
        }
        expect_equal(sum(weights), 1)
        for(j in 1:4) {
            ## and for each one, row-wise, we have values
            if (j == 1) { g <- lapply(1:length(ll), function(x) rep(1, 10))}
            if (j == 2) { g <- lapply(1:length(ll), function(x) runif(10)) }
            if (j == 3) { g <- lapply(1:length(ll), function(x) rep(1e4, 10)) }
            if (j == 4) { g <- lapply(1:length(ll), function(x) rep(1 / weights[x], 10))}
            ## re-shape to matrix
            g <- lapply(g, function(x) array(x, c(2, 5)))
            ##
            gA <- weights[1] * g[[1]]
            for(i in 2:length(g)) { gA <- gA + weights[i] * g[[i]]}
            gA2 <- weighted_average_on_the_fly(g, ll, log_mult_max, equal_weighting = equal_weighting, language = "R")
            ## gA3
            gA3 <- weighted_average_on_the_fly(g, ll, log_mult_max, equal_weighting = equal_weighting, language = "Rcpp")
            ## yay
            expect_equal(gA, gA2)
            expect_equal(gA, gA3)
        }
    }
    ## also check that equal weighter works

})


test_that("gibbs-nipt works in R and CPP and gives the same thing", {

    ##skip("wip")
    verbose <- FALSE
    ff <- 0.2
    s <- 1
    S <- 1
    test_package <- make_fb_test_package(
        K = 4,
        nReads = 10,
        nSNPs = 3,
        gridWindowSize = NA,
        method = "triploid-nipt",
        ff = ff,
        seed = 100,
        eHapsMin = 0.01,
        S = S
    )
    grid <- test_package$grid
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    ## it is probably better if fit is much poorer
    eHapsCurrent_tc[] <- array(runif(prod(dim(eHapsCurrent_tc))), dim(eHapsCurrent_tc))
    KKK <- K * K * K
    sampleReads <- test_package$sampleReads[1]
    ## check that with no emissions, move one forward keeps sum the same
    true_H <- test_package$true_H[1]
    nReads <- length(sampleReads)
    wif <- sapply(sampleReads, function(x) x[[2]]) + 1

    set.seed(99)

    to_run <- 1:8
    
    for(i_run in to_run) {

        ##
        if (verbose) {
            print(paste0("-------i_run = ", i_run, " ---------------"))
        }

        haploid_gibbs_equal_weighting <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)[i_run]
        n_gibbs_burn_in_its <- c(1, 1, 2, 2, 1, 1, 2, 2)[i_run]
        n_gibbs_sample_its <- c(1, 1, 2, 2, 1, 1, 2, 2)[i_run]
        n_gibbs_starts <- c(1, 1, 1, 1, 2, 2, 2, 2)[i_run]

        double_list_of_starting_read_labels <- lapply(1:S, function(s) {
            lapply(1:n_gibbs_starts, function(i_sampling) {
                x <- runif(nReads);
                H <- array(0L, nReads)
                for(iRead in 1:nReads) {
                    if (x[iRead] < 0.5) {
                        H[iRead] = 1
                    } else if ((0.5 <= x[iRead]) & (x[iRead] < (0.5 + ff / 2))) {
                        H[iRead] = 2
                    } else {
                        H[iRead] = 3
                    }
                }
                return(H)
            })
        })

        set.seed(123)        
        outR <- forwardBackwardGibbsNIPT(
            sampleReads = sampleReads,
            priorCurrent_m = priorCurrent_m,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            eHapsCurrent_tc = eHapsCurrent_tc,
            transMatRate_tc_H = transMatRate_tc_H,
            ff = ff,
            blocks_for_output = array(0, c(1, 1)),
            grid = grid,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            n_gibbs_sample_its = n_gibbs_sample_its,
            n_gibbs_starts = n_gibbs_starts,
            use_starting_read_labels = TRUE,
            return_p_store = TRUE,
            true_H = true_H,
            return_alpha = TRUE,
            verbose = verbose,
            double_list_of_starting_read_labels = double_list_of_starting_read_labels,
            return_extra = TRUE,
            haploid_gibbs_equal_weighting = haploid_gibbs_equal_weighting,
            record_read_set = TRUE
        )

        set.seed(123)
        outRCPP <- rcpp_forwardBackwardGibbsNIPT(
            sampleReads = sampleReads,
            priorCurrent_m = priorCurrent_m,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            eHapsCurrent_tc = eHapsCurrent_tc,
            transMatRate_tc_H = transMatRate_tc_H,
            ff = ff,
            blocks_for_output = array(0, c(1, 1)),
            grid = grid,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            n_gibbs_sample_its = n_gibbs_sample_its,
            n_gibbs_starts = n_gibbs_starts,
            use_starting_read_labels = TRUE,
            suppressOutput = 1,
            verbose = verbose,
            return_p_store = TRUE,
            return_alpha = TRUE,
            return_gamma = TRUE,
            return_extra = TRUE,
            double_list_of_starting_read_labels = double_list_of_starting_read_labels,
            haploid_gibbs_equal_weighting = haploid_gibbs_equal_weighting,
            alphaHat_t1 = array(0, c(1, 1)),
            alphaHat_t2 = array(0, c(1, 1)),
            alphaHat_t3 = array(0, c(1, 1)),
            betaHat_t1 = array(0, c(1, 1)),
            betaHat_t2 = array(0, c(1, 1)),
            betaHat_t3 = array(0, c(1, 1)),
            hapSum_tc = array(0, c(1, 1, 1)),
            record_read_set = TRUE,
            wif0 = integer(1),
            L_grid = integer(1)
        )

        expect_equal(outR[["H_class"]], outRCPP[["H_class"]])
        expect_equal(outR$eMatRead_t, outRCPP$eMatRead_t)
        expect_equal(outR$eMatGrid_t1, outRCPP$eMatGrid_t1)
        expect_equal(outR$eMatGrid_t2, outRCPP$eMatGrid_t2)
        expect_equal(outR$eMatGrid_t3, outRCPP$eMatGrid_t3)
        expect_equal(outR$alphaHat_t1, outRCPP$alphaHat_t1)
        expect_equal(outR$alphaHat_t2, outRCPP$alphaHat_t2)
        expect_equal(outR$alphaHat_t3, outRCPP$alphaHat_t3)
        expect_equal(outR$betaHat_t1, outRCPP$betaHat_t1)
        expect_equal(outR$betaHat_t2, outRCPP$betaHat_t2)
        expect_equal(outR$betaHat_t3, outRCPP$betaHat_t3)
        
        expect_equal(as.numeric(outR$c1), as.numeric(outRCPP$c1))        
        expect_equal(as.numeric(outR$c2), as.numeric(outRCPP$c2))
        expect_equal(as.numeric(outR$c3), as.numeric(outRCPP$c3))

        if (verbose) {
            print("wif")
            print(wif)
            print("true_H")
            print(true_H)
            print("----------colsum difference")
            print(colSums(abs(outR$alphaHat_t1 - outRCPP$alphaHat_t1)))
            print("----------head(outR_pstore)")
            print(head(outR$p_store, 10))
            print("----------head(outRCPP_pstore)")
            print(head(outRCPP$p_store, 10))
        }

        expect_equal(c(outR$H), c(outRCPP$H)   )
        ## check ending read labels vs p_store
        ## expect_equal(c(outR$hg_ll_rescaled), c(outRCPP$hg_ll_rescaled))
        expect_equal(outR$likelihoods, outRCPP$likelihoods)
        expect_equal(outR$double_list_of_ending_read_labels[[s]], outRCPP$double_list_of_ending_read_labels[[s]])
        ## all columns agree, except
        for(col in setdiff(colnames(outR$p_store), "agreePer")) {
            expect_equal(max(abs(outR$p_store[, col] - outRCPP$p_store[, col])) < 1e-4, TRUE)
        }
        
        ##
        for(i_start in 1:n_gibbs_starts) {
            for(i_sample_it in 1:n_gibbs_sample_its) {
                i_result_it <- n_gibbs_sample_its * (i_start - 1) + i_sample_it
                ##
                offset1 <- nReads * (n_gibbs_burn_in_its + n_gibbs_sample_its) * (i_start - 1)
                offset2 <- nReads * n_gibbs_burn_in_its
                w <- offset1 + offset2 + nReads * (i_sample_it - 1) + 1:nReads
                expect_equal(c(outR$double_list_of_ending_read_labels[[s]][[i_result_it]]), as.numeric(outR$p_store[w, "h_rN"]))
                expect_equal(c(outRCPP$double_list_of_ending_read_labels[[s]][[i_result_it]]), as.numeric(outRCPP$p_store[w, "h_rN"]))
            }
        }
        
        ## end of working normally
        ## check weighting
        if (haploid_gibbs_equal_weighting) {
            n <- (n_gibbs_starts * n_gibbs_sample_its)
            weights <- rep(1 / n, n)
        } else {
            ll <- outR$per_it_likelihoods[outR$per_it_likelihoods[, "i_result_it"] >= 0, "p_H_given_O_L_up_to_C"]
            ll <- ll - max(ll)
            weights <- exp(ll)
            weights <- weights / sum(weights)
        }
        ##
        
        for(what in c("genProbsM_t", "genProbsF_t")) {
            object_average <- array(0, dim(outR$per_iteration_store[[1]][[what]]))
            for(i_outer in 1:(n_gibbs_starts * n_gibbs_sample_its)) {
                object <- outR$per_iteration_store[[i_outer]][[what]]
                object_average <- object_average + weights[i_outer] * object
            }
            expect_equal(outR[[what]], object_average)
        }
        expect_equal(outR$per_it_likelihoods, outRCPP$per_it_likelihoods)
        ##
        expect_equal(outR$genProbsM_t, outRCPP$genProbsM_t)
        expect_equal(outR$genProbsF_t, outRCPP$genProbsF_t)
        print(outR$gammaMT_t)
        print(outRCPP$gammaMT_t)
        expect_equal(outR$gammaMT_t, outRCPP$gammaMT_t)
        expect_equal(outR$gammaMU_t, outRCPP$gammaMU_t)
        expect_equal(outR$gammaP_t, outRCPP$gammaP_t)
        expect_equal(outR$hapProbs_t, outRCPP$hapProbs_t)

        ##
        if (1 == 0) {
            ## check hapProbs is OK with genProbs
            a <- c(outR$hapProbsMT)
            b <- c(outR$hapProbsMU)
            c <- c(outR$hapProbsP)
            print(head(a))
            genProbsM_t <- rbind(
            (1 - a) * (1 - b),
            (1 - a) * b + (1 - b) * a,
            a * b
            )
            genProbsF_t <- rbind(
            (1 - a) * (1 - c),
            (1 - a) * c + (1 - c) * a,
            a * c
            )
            print(paste0("-------i_run = ", i_run, " ---------------"))
            print("-----------")
            print(outR$genProbsM_t)
            print(genProbsM_t)
            print("==========")
            print(outR$genProbsF_t)
            print(genProbsF_t)
            print("...........")
            ## hapProbs average out
            (outR[["temp"]][[1]]$hapProbsMT + outR[["temp"]][[2]]$hapProbsMT) / 2 - a
            (outR[["temp"]][[1]]$hapProbsMU + outR[["temp"]][[2]]$hapProbsMU) / 2 - b
            (outR[["temp"]][[1]]$hapProbsP + outR[["temp"]][[2]]$hapProbsP) / 2 - c
            ## OK - gammas are the same!
            (outR[["temp"]][[1]][["gammaMT_t"]] + outR[["temp"]][[2]][["gammaMT_t"]]) / 2 - outR$gammaMT_t
            ## OK - this makes sense
            (colSums(outR[["temp"]][[1]][["gammaMT_t"]] * eHapsCurrent_t) +
             colSums(outR[["temp"]][[2]][["gammaMT_t"]] * eHapsCurrent_t)) / 2 - a
            ## genProbs averaging makes sense
            outR$genProbsM_t - (outR[["temp"]][[1]]$genProbsM_t + outR[["temp"]][[2]]$genProbsM_t) / 2
            ##
            a <- c(outR[["temp"]][[2]]$hapProbsMT)
            b <- c(outR[["temp"]][[2]]$hapProbsMU)
            genProbsM_t <- rbind(
            (1 - a) * (1 - b),
            (1 - a) * b + (1 - b) * a,
            a * b
            )
            outR[["temp"]][[2]]$genProbsM_t - genProbsM_t
            ##
            outR$genProbsM_t -  genProbsM_t ## still differences! substantial! on one at least
            expect_equal(genProbsM_t, outR$genProbsM_t)
            expect_equal(genProbsF_t, outR$genProbsF_t)
        }

    }

})



test_that("gibbs-nipt initializing iteratively produces same alpha and beta, and same in R as Rcpp", {

    ## need this as Rcpp uses older version. eventually this will break I imagine when Rcpp is updated
    suppressWarnings(RNGversion("3.5.0"))

    ##skip("wip")
    ff <- 0.2
    test_packages <- lapply(c(NA, 3), function(gridWindowSize) {
        make_fb_test_package(
            K = 20,
            nReads = 50,
            nSNPs = 50,
            gridWindowSize = gridWindowSize,
            method = "triploid-nipt",
            ff = ff,
            seed = 100,
            eHapsMin = 0.01,
            S = 1
        )
    })
    names(test_packages) <- c("NA", "3")

    test_language <- function(
        language,
        gibbs_initialize_iteratively = FALSE,
        gibbs_initialize_at_first_read = TRUE,
        n_gibbs_burn_in_its = 0,
        n_gibbs_sample_its = 1,
        artificial_relabel = -1,
        do_block_resampling = FALSE
    ) {
        ##
        ## gridWindowSize is in global of this environment
        ##
        if (is.na(gridWindowSize)) {
            gridWindowSize <- "NA"
        }
        test_package <- test_packages[[as.character(gridWindowSize)]]
        grid <- test_package$grid
        nSNPs <- test_package$nSNPs
        nGrids <- test_package$nGrids
        K <- test_package$K
        transMatRate_tc_H <- test_package$transMatRate_tc_H
        alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
        priorCurrent_m <- test_package$priorCurrent_m
        eHapsCurrent_tc <- test_package$eHapsCurrent_tc
        KKK <- K * K * K
        sampleReads <- test_package$sampleReads
        ## check that with no emissions, move one forward keeps sum the same
        true_H <- test_package$true_H
        nReads <- length(sampleReads)
        wif <- sapply(sampleReads, function(x) x[[2]]) + 1
        n_gibbs_starts <- 1
        ##
        if (language == "R") {
            fbgn <- forwardBackwardGibbsNIPT
        } else if ((language == "Rcpp")) {
            fbgn <- rcpp_forwardBackwardGibbsNIPT
        }
        ##
        ##  here we do two passes
        ## for the first one, run, get labels at the end plus alpha
        ## for the second one, given the labels, get the alpha
        ## the goal is to make sure the forward/backwards alpha use is OK
        ##
        ## first one - initialize, somehow, get read labels
        verbose <- FALSE
        set.seed(72)
        running_seeds <- c(24)
        starting_read_labels <- sample(c(1, 2, 3), replace = TRUE, size = length(sampleReads))
        initial_starting_read_labels <- starting_read_labels
        double_list_of_starting_read_labels <- list(list(starting_read_labels))
        rm(starting_read_labels)
        out1 <- fbgn(
            sampleReads = sampleReads, priorCurrent_m = priorCurrent_m, alphaMatCurrent_tc = alphaMatCurrent_tc, eHapsCurrent_tc = eHapsCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, ff = ff, blocks_for_output = array(0, c(1, 1)), grid = grid, n_gibbs_starts = n_gibbs_starts, return_p_store = TRUE, return_alpha = TRUE, verbose = verbose, double_list_of_starting_read_labels = double_list_of_starting_read_labels, seed_vector = running_seeds, return_extra = TRUE, alphaHat_t1 = array(0, c(1, 1)), alphaHat_t2 = array(0, c(1, 1)), alphaHat_t3 = array(0, c(1, 1)), betaHat_t1 = array(0, c(1, 1)), betaHat_t2 = array(0, c(1, 1)), betaHat_t3 = array(0, c(1, 1)), hapSum_tc = array(0, c(1, 1, 1)), wif0 = integer(1), L_grid = integer(1),
            use_starting_read_labels = TRUE,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its ,
            n_gibbs_sample_its = n_gibbs_sample_its,
            gibbs_initialize_iteratively = gibbs_initialize_iteratively,
            gibbs_initialize_at_first_read = gibbs_initialize_at_first_read,
            artificial_relabel = artificial_relabel,
            do_block_resampling = do_block_resampling,
            record_read_set = TRUE
        )
        ##
        ## second - pass those labels through
        ## note - n_gibbs_starts is taken as number of read labels passed
        ##
        double_list_of_starting_read_labels <- out1$double_list_of_ending_read_labels
        ## will determine
        out2 <- fbgn(
            sampleReads = sampleReads, priorCurrent_m = priorCurrent_m, alphaMatCurrent_tc = alphaMatCurrent_tc, eHapsCurrent_tc = eHapsCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, ff = ff, blocks_for_output = array(0, c(1, 1)), grid = grid, return_p_store = TRUE, return_alpha = TRUE, verbose = verbose, double_list_of_starting_read_labels = double_list_of_starting_read_labels, seed_vector = running_seeds, return_extra = TRUE, alphaHat_t1 = array(0, c(1, 1)), alphaHat_t2 = array(0, c(1, 1)), alphaHat_t3 = array(0, c(1, 1)), betaHat_t1 = array(0, c(1, 1)), betaHat_t2 = array(0, c(1, 1)), betaHat_t3 = array(0, c(1, 1)), hapSum_tc = array(0, c(1, 1, 1)),  wif0 = integer(1), L_grid = integer(1),
            use_starting_read_labels = TRUE,
            n_gibbs_starts = length(double_list_of_starting_read_labels[[1]]), ## might be longer!
            n_gibbs_burn_in_its = 0,
            n_gibbs_sample_its = 0, ## what
            gibbs_initialize_iteratively = FALSE,
            do_block_resampling = FALSE,
            suppressOutput = 1
        )
        ##
        for(w in c("alphaHat_t1", "alphaHat_t2", "alphaHat_t3", "c1", "c2", "c3")) {
            expect_equal(out1[[w]], out2[[w]])
        }
        return(list(out1 = out1, out2 = out2, initial_starting_read_labels = initial_starting_read_labels))
    }

    ## so here, test initialization
    ## also, test if just the same read labels are given, get same gamma, etc

    for(gridWindowSize in c(NA, 3)) {

        outRA <- test_language(language = "R", gibbs_initialize_iteratively = FALSE)
        outRcppA <- test_language(language = "Rcpp", gibbs_initialize_iteratively = FALSE)
        
        outRB <- test_language(language = "R", gibbs_initialize_iteratively = TRUE, gibbs_initialize_at_first_read = TRUE, n_gibbs_burn_in_its = 0, n_gibbs_sample_its = 1)
        outRcppB <- test_language(language = "Rcpp", gibbs_initialize_iteratively = TRUE, gibbs_initialize_at_first_read = TRUE, n_gibbs_burn_in_its = 0, n_gibbs_sample_its = 1)
        outRC <- test_language(language = "R", gibbs_initialize_iteratively = TRUE, gibbs_initialize_at_first_read = FALSE, n_gibbs_burn_in_its = 1, n_gibbs_sample_its = 1)
        outRcppC <- test_language(language = "Rcpp", gibbs_initialize_iteratively = TRUE, gibbs_initialize_at_first_read = FALSE, n_gibbs_burn_in_its = 1, n_gibbs_sample_its = 1)

        expect_equal(as.numeric(outRA[["out1"]]$c1), as.numeric(outRcppA[["out1"]]$c1))
        expect_equal(as.numeric(outRA[["out1"]]$c1), as.numeric(outRcppA[["out1"]]$c1))
        expect_equal(as.numeric(outRB[["out1"]]$c1), as.numeric(outRcppB[["out1"]]$c1))

        ## this may fail on R >= 3.6.0
        ## can see with e.g. suppressWarnings(RNGversion("3.5.0"))
        ## hopefully Rcpp figures it out
        expect_equal(as.numeric(outRC[["out1"]]$c1), as.numeric(outRcppC[["out1"]]$c1))

        ## also test can do artificial re-labeling for tricky cases
        for(i in 1:6) {
            outRD <- test_language(language = "R", gibbs_initialize_iteratively = FALSE, artificial_relabel = i, do_block_resampling = TRUE)
            outRcppD <- test_language(language = "Rcpp", gibbs_initialize_iteratively = FALSE, artificial_relabel = i, do_block_resampling = TRUE)
        }
    }

    ## reset!
    suppressWarnings(RNGversion("3.6.0"))

})








test_that("speed X", {

    skip("optional speedtest")
    library(microbenchmark)

     iGrid <- 100
     K <- 100
     nGrids <- 500
     alphaHat_t <- array(runif(nGrids * K), c(K, nGrids))
     transMatRate_t_H <- array(runif(nGrids * 2), c(2, nGrids))
     eMatHapSNP_t <- array(runif(K * nGrids), c(K, nGrids))
     alphaMat_t <- array(runif(K * nGrids), c(K, nGrids))
     c <- array(runif(nGrids), nGrids)

     print(microbenchmark(
         rcpp_alpha_forward_oneORIGINAL(
             iGrid = iGrid,
             K = K,
             alphaHat_t = alphaHat_t,
             transMatRate_t_H = transMatRate_t_H,
             eMatHapSNP_t = eMatHapSNP_t,
             alphaMat_t = alphaMat_t,
             c = c,
             normalize = TRUE
         ),
         rcpp_alpha_forward_one(
             iGrid = iGrid,
             K = K,
             alphaHat_t = alphaHat_t,
             transMatRate_t_H = transMatRate_t_H,
             eMatHapSNP_t = eMatHapSNP_t,
             alphaMat_t = alphaMat_t,
             c = c,
             normalize = TRUE
         ),
         times = 100
     ))


})


test_that("can do multiple gibbs samplings", {

    run_speed_test <- FALSE

    if (run_speed_test) {
        n_gibbs_starts <- 2
        suppressOutput <- 0
        K <- 2000
        nReads <- 100
        nSNPs <- 1000
        temp_file <- "temp.gibbs.nipt.RData"
        n_gibbs_burn_in_its <- 4
        n_gibbs_sample_its <- 1
    } else {
        n_gibbs_starts <- 2
        suppressOutput <- 1
        K <- 20
        nSNPs <- 10
        nReads <- 50
        temp_file <- "XX.RData"
        n_gibbs_burn_in_its <- 0
        n_gibbs_sample_its <- 1
    }


    S <- 1
    ff <- 0.2
    rebuild <- TRUE
    if (rebuild | ((!run_speed_test) | (run_speed_test & !file.exists(temp_file)))) {
        test_package <- make_fb_test_package(
            K = K,
            nReads = nReads,
            nSNPs = nSNPs,
            gridWindowSize = NA,
            method = "triploid-nipt",
            ff = ff,
            seed = 100,
            eHapsMin = 0.01,
            return_eMatGridTri_t = FALSE,
            S = S
        )
        if (run_speed_test) {
            save(test_package, file = temp_file, compress = FALSE)
        }
    } else {
        load(temp_file)
    }

    grid <- test_package$grid
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    KKK <- K * K * K
    sampleReads <- test_package$sampleReads
    eMatHap_t <- test_package$eMatHap_t
    ## check that with no emissions, move one forward keeps sum the same
    true_H <- test_package$true_H
    nReads <- length(sampleReads)
    wif <- sapply(sampleReads, function(x) x[[2]]) + 1


    ##
    double_list_of_starting_read_labels <- lapply(1:S, function(s) {
        lapply(1:n_gibbs_starts, function(i) {
            set.seed(100 * i)
            x <- runif(nReads);
            H <- array(0L, nReads)
            for(iRead in 1:nReads) {
                if (x[iRead] < 0.5) {
                    H[iRead] = 1
                } else if ((0.5 <= x[iRead]) & (x[iRead] < (0.5 + ff / 2))) {
                    H[iRead] = 2
                } else {
                    H[iRead] = 3
                }
            }
            return(H)
        })
    })

    running_seeds <- c(24, 991, 5, 10, 20, 52, 2, 90, 4, 2292)
    if (n_gibbs_starts > length(running_seeds)) {
        stop("need more seeds")
    }

    test_by_language <- function(language) {
        set.seed(900)
        ## print(paste0("===============", language, "==============="))
        if (language == "R") {
            fbgn <- forwardBackwardGibbsNIPT
        } else if ((language == "Rcpp")) {
            fbgn <- rcpp_forwardBackwardGibbsNIPT
        } else if ((language == "RcppOriginal")) {
            fbgn <- rcpp_forwardBackwardGibbsNIPT_ORIGINAL
        }
        ## yuck - doesn't deal with sampling at all
        out <- lapply(1:n_gibbs_starts, function(i_sampling) {
            set.seed(running_seeds[i_sampling])
            out <- fbgn(
                sampleReads = sampleReads,
                priorCurrent_m = priorCurrent_m,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                eHapsCurrent_tc = eHapsCurrent_tc,
                transMatRate_tc_H = transMatRate_tc_H,
                ff = ff,
                blocks_for_output = array(0, c(1, 1)),
                grid = grid,
                n_gibbs_burn_in_its = n_gibbs_burn_in_its,
                n_gibbs_sample_its = n_gibbs_sample_its,
                use_starting_read_labels = TRUE,
                return_gamma = TRUE,
                return_genProbs = TRUE,
                return_hapProbs = TRUE,
                return_alpha = TRUE,
                double_list_of_starting_read_labels = list(list(double_list_of_starting_read_labels[[1]][[i_sampling]])),
                alphaHat_t1 = array(0, c(1, 1)),
                alphaHat_t2 = array(0, c(1, 1)),
                alphaHat_t3 = array(0, c(1, 1)),
                betaHat_t1 = array(0, c(1, 1)),
                betaHat_t2 = array(0, c(1, 1)),
                betaHat_t3 = array(0, c(1, 1)),
                hapSum_tc = array(0, c(1, 1, 1)),
                wif0 = integer(1),
                L_grid = integer(1)
            )
            return(out)
        })
        out_manual <- list(
            genProbsM_t = out[[1]]$genProbsM_t,
            genProbsF_t = out[[1]]$genProbsF_t,
            gammaMT_t = out[[1]]$gammaMT_t,
            gammaMU_t = out[[1]]$gammaMU_t,
            gammaP_t = out[[1]]$gammaP_t,
            everything = out
        )
        for(j in 2:n_gibbs_starts) {
            out_manual[["genProbsM_t"]] <- out_manual[["genProbsM_t"]] + out[[j]]$genProbsM_t
            out_manual[["genProbsF_t"]] <- out_manual[["genProbsF_t"]] + out[[j]]$genProbsF_t
            out_manual[["gammaMT_t"]] <- out_manual[["gammaMT_t"]] + out[[j]]$gammaMT_t
            out_manual[["gammaMU_t"]] <- out_manual[["gammaMU_t"]] + out[[j]]$gammaMU_t
            out_manual[["gammaP_t"]] <- out_manual[["gammaP_t"]] + out[[j]]$gammaP_t
        }
        out_manual[["genProbsM_t"]] <- out_manual[["genProbsM_t"]] / n_gibbs_starts
        out_manual[["genProbsF_t"]] <- out_manual[["genProbsF_t"]] / n_gibbs_starts
        out_manual[["gammaMT_t"]] <- out_manual[["gammaMT_t"]] / n_gibbs_starts
        out_manual[["gammaMU_t"]] <- out_manual[["gammaMU_t"]] / n_gibbs_starts
        out_manual[["gammaP_t"]] <- out_manual[["gammaP_t"]] / n_gibbs_starts
        out_auto <- fbgn(
            sampleReads = sampleReads,
            priorCurrent_m = priorCurrent_m,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            eHapsCurrent_tc = eHapsCurrent_tc,
            transMatRate_tc_H = transMatRate_tc_H,
            ff = ff,
            blocks_for_output = array(0, c(1, 1)),
            grid = grid,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            n_gibbs_sample_its = n_gibbs_sample_its,
            use_starting_read_labels = TRUE,
            double_list_of_starting_read_labels = double_list_of_starting_read_labels,
            n_gibbs_starts = n_gibbs_starts,
            seed_vector = running_seeds,
            return_gamma = TRUE,
            return_alpha = TRUE,
            suppressOutput = suppressOutput,
            verbose = FALSE,
            haploid_gibbs_equal_weighting = TRUE,
            alphaHat_t1 = array(0, c(1, 1)),
            alphaHat_t2 = array(0, c(1, 1)),
            alphaHat_t3 = array(0, c(1, 1)),
            betaHat_t1 = array(0, c(1, 1)),
            betaHat_t2 = array(0, c(1, 1)),
            betaHat_t3 = array(0, c(1, 1)),
            hapSum_tc = array(0, c(1, 1, 1)),
            wif0 = integer(1),
            L_grid = integer(1)
        )
        ## this is basically all I care about
        return(
            list(
                auto = out_auto,
                manual = out_manual
            )
        )
    }

    if (run_speed_test) {
        print("-----------OLD------------")
        results_RcppOriginal <- test_by_language("RcppOriginal")
        print("-----------NEW------------")
    }

    results_Rcpp <- test_by_language("Rcpp")
    expect_equal(1, 1) ## just so not empty

    if (!run_speed_test) {

        results_R <- test_by_language("R")
        ## also check they agree across
        expect_equal(results_R$auto$gammaMT_t, results_R$manual[["gammaMT_t"]])
        expect_equal(results_R$auto$gammaMU_t, results_R$manual[["gammaMU_t"]])
        expect_equal(results_R$auto$gammaP_t, results_R$manual[["gammaP_t"]])
        expect_equal(results_R$auto$genProbsM_t, results_R$manual[["genProbsM_t"]])
        expect_equal(results_R$auto$genProbsF_t, results_R$manual[["genProbsF_t"]])

        ## tests within Rcpp
        expect_equal(results_Rcpp$auto$gammaMT_t, results_Rcpp$manual[["gammaMT_t"]])
        expect_equal(results_Rcpp$auto$gammaMU_t, results_Rcpp$manual[["gammaMU_t"]])
        expect_equal(results_Rcpp$auto$gammaP_t, results_Rcpp$manual[["gammaP_t"]])
        expect_equal(results_Rcpp$auto$genProbsM_t, results_Rcpp$manual[["genProbsM_t"]])
        expect_equal(results_Rcpp$auto$genProbsF_t, results_Rcpp$manual[["genProbsF_t"]])

        ## manual gamma Rcpp vs R not working
        expect_equal(
            results_R$manual$everything[[1]]$alphaHat_t1,
            results_Rcpp$manual$everything[[1]]$alphaHat_t1
        )
        ## test
        expect_equal(results_R$auto$gammaMT_t, results_Rcpp$auto$gammaMT_t)
        expect_equal(results_R$auto$gammaMU_t, results_Rcpp$auto$gammaMU_t)
        expect_equal(results_R$auto$gammaP_t, results_Rcpp$auto$gammaP_t)
        expect_equal(results_R$auto$genProbsM_t, results_Rcpp$auto$genProbsM_t)
        expect_equal(results_R$auto$genProbsF_t, results_Rcpp$auto$genProbsF_t)

    }

})


test_that("can determine labels", {

    return_neutral <- FALSE
    rc <- c(14, 17, 19)
    ff <- 0.2
    expect_equal(dim(determine_label_probabilities(rc, ff, return_neutral)     ), c(3, 3))
    expect_equal(dim(rcpp_determine_label_probabilities(rc, ff, return_neutral)     ), c(3, 3))
    expect_equal(determine_label_probabilities(rc, ff, return_neutral), rcpp_determine_label_probabilities(rc, ff, return_neutral) )

})


test_that("label values do not suck", {

    return_neutral <- FALSE
    rc <- c(13733, 12944, 10095)
    ff <- 0.2
    determine_label_probabilities(rc, ff, return_neutral)
    expect_equal(determine_label_probabilities(rc, ff, return_neutral), rcpp_determine_label_probabilities(rc, ff, return_neutral) )

})


test_that("averaging logic makes sense", {


    ## so idea is Gibbs sampling
    ## like P(G | O, \lambda)
    ## = \sum_{h in \scrict{H}} P(G, H | O, \lambda)
    ## = \sum_{h in \scrict{H}} P(G | H, O, \lambda) P(H | O, \lambda)
    ## \prop \sum_{h*} P(G | H, O, \lambda)
    ## where h* ~ P(H | O, \lambda)
    ## big question - real Gibbs is fine
    ## BUT we might get stuck, due to initialization
    ## and not be able to move between them properly
    ## understand that here, and how to move, or re-estimate importance of individual spots

    haploid_fit <- function(whichReads) {
        out <- forwardBackwardHaploid(
            sampleReads = sampleReads[whichReads],
            eHapsCurrent_tc = eHapsCurrent_tc,
            grid = grid,
            priorCurrent_m = priorCurrent_m,
            transMatRate_tc_H = transMatRate_tc_H,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            alphaHat_t = array(0, c(1, 1)),
            betaHat_t = array(0, c(1, 1)),
            gamma_t = array(0, c(1, 1)),
            eMatGrid_t = array(0, c(1, 1)),            
            maxDifferenceBetweenReads = 1000,
            maxEmissionMatrixDifference = 1e10,
            Jmax = 10,
            suppressOutput = 1,
            model = as.integer(9),
            gammaSum0_tc = array(0, c(1, 1, 1)),
            gammaSum1_tc = array(0, c(1, 1, 1)),
            alphaMatSum_tc = array(0, c(1, 1, 1)),
            hapSum_tc = array(0, c(1, 1, 1)),
            priorSum_m = array(0, c(1, 1)),
            pRgivenH1_m = array(0, c(1, 1)),
            pRgivenH2_m = array(0, c(1, 1)),
            run_pseudo_haploid = FALSE,
            blocks_for_output = array(0, c(1, 1)),
            prev_list_of_alphaBetaBlocks = list(),
            return_gamma = TRUE
        )
        ll <- -sum(log(out$c))
        return(list(out = out, ll =ll))
    }
    diploid_fit <- function() {
        out <- forwardBackwardDiploid(
            sampleReads = sampleReads,
            grid = grid,
            priorCurrent_m = priorCurrent_m,
            transMatRate_tc_D = transMatRate_tc_D,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            eHapsCurrent_tc = eHapsCurrent_tc,
            maxDifferenceBetweenReads = 1000,
            maxEmissionMatrixDifference = 1e10,
            Jmax = 10,
            suppressOutput = 1,
            run_fb_subset = FALSE,
            gammaSum0_tc = array(0, c(1, 1, 1)),
            gammaSum1_tc = array(0, c(1, 1, 1)),
            alphaMatSum_tc = array(0, c(1, 1, 1)),
            hapSum_tc = array(0, c(1, 1, 1)),
            priorSum_m = array(0, c(1, 1)),
            alphaHat_t = array(0, c(1, 1)),
            betaHat_t = array(0, c(1, 1)),
            gamma_t = array(0, c(1, 1)),
            eMatGrid_t = array(0, c(1, 1)),                        
            blocks_for_output = array(0, c(1, 1))    ,
            return_extra = TRUE,
            return_genProbs = TRUE,
            return_gamma = TRUE,
            snp_start_1_based = 1,
            snp_end_1_based = ncol(eHapsCurrent_tc),
            prev_list_of_alphaBetaBlocks = list(),
            rescale_eMatGrid_t = FALSE,
            rescale_eMatRead_t = FALSE
        )
        return(list(out = out, ll = -sum(log(out$c))))
    }



    ## set.seed(100)
    ## n_snps <- 10
    ## hap_probs <- c(0.5, 0.5) ## 3 haps
    ## grid <- 0:(n_snps - 1)
    ## nGrids <- length(grid)
    ## K <- 2
    ## S <- 1
    ## sampleRead <- function(x, grid) {
    ##     y <- x + c(0, 1, 2, 3)
    ##     y <- y[y < (n_snps -1)]
    ##     return(list(
    ##         length(y) -1, round(median(grid[y + 1])),
    ##         matrix(sample(c(10, -10), length(y), replace = TRUE), ncol = 1),
    ##         matrix(y, ncol = 1)
    ##     ))
    ## }
    ## sampleReads <- lapply(1:(n_snps - 3), sampleRead, grid = grid)
    ## eHapsCurrent_tc <- array(runif(K * n_snps * S), c(K, n_snps, S))
    ## ## ignore transitions?
    ## nReads <- length(sampleReads)
    ## sigmaCurrent <- rep(0.99, nGrids - 1)
    ## alphaMatCurrent_t <- array(1 / K, c(K, nGrids - 1))
    ## priorCurrent <- runif(K) / K
    ## priorCurrent <- priorCurrent / sum(priorCurrent)
    ## transMatRate_t_H <- get_transMatRate(method = "diploid-inbred", sigmaCurrent)
    ## transMatRate_t_D <- get_transMatRate(method = "diploid", sigmaCurrent)
    ## ##
    set.seed(100)
    S <- 1
    s <- 1
    nSNPs <- 10
    ff <- 0
    K <- 2
    nReads <- 8
    test_package <- make_fb_test_package(
        K = K,
        nReads = nReads,
        nSNPs = nSNPs,
        gridWindowSize = NA,
        method = "triploid-nipt",
        ff = ff,
        seed = 100,
        eHapsMin = 0.01,
        return_eMatGridTri_t = FALSE,
        S = S
    )
    grid <- test_package$grid
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    transMatRate_tc_D <- test_package$transMatRate_tc_D
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    KKK <- K * K * K
    sampleReads <- test_package$sampleReads
    ## check that with no emissions, move one forward keeps sum the same
    true_H <- test_package$true_H
    nReads <- length(sampleReads)
    wif <- sapply(sampleReads, function(x) x[[2]]) + 1
    hap_probs <- c(0.5, 0.5) ## 3 haps

    i_snp_pg <- 5 ## for P(G), what snp
    i_g_pg <- 1 ## for P(G), which genotype

    p_H_given_L <- function(H) {
        return(sum(log(hap_probs[H])))
    }
    p_G_given_O_H_L <- function(whichReads) {
        haploid_fit(H == 1)$ll + haploid_fit(H == 2)$ll
    }

    H <- sample(c(1, 2), size = nReads, prob = hap_probs, replace = TRUE)

    ## calculate for all possibly H
    out <- array(NA, c(2 ** nReads, nReads + 7))
    colnames(out) <- c(
        paste0("H_r", 1:nReads),
        "p_O1_given_H1_L",
        "p_O2_given_H2_L",
        "p_H_given_L",
        "p_O_given_H_L",
        "p_O_H_given_L",
        "p_H_given_O_L",
        "p_G_given_O_H_lambda"
    )
    for(i in 0:(2 ** nReads - 1)) {
        x <- intToBits(i)
        H <- as.integer(x)[1:nReads] + 1
        p_gh_given_o_h1_lambda <- colSums(haploid_fit(H == 1)$out$list_of_gamma_t[[s]] * eHapsCurrent_tc[, , s])
        p_gh_given_o_h2_lambda <- colSums(haploid_fit(H == 2)$out$list_of_gamma_t[[s]] * eHapsCurrent_tc[, , s])
        a <- p_gh_given_o_h1_lambda[i_snp_pg]
        b <- p_gh_given_o_h2_lambda[i_snp_pg]
        p_g_given_o_h_lambda <- c((1 - a) * (1 - b), a * (1 - b) + b * (1 - a), a * b)[i_g_pg]
        out[i + 1, 1:nReads] <- H
        out[i + 1, "p_O1_given_H1_L"] <- haploid_fit(H == 1)$ll
        out[i + 1, "p_O2_given_H2_L"] <- haploid_fit(H == 2)$ll
        out[i + 1, "p_H_given_L"] <- p_H_given_L(H)
        out[i + 1, "p_O_given_H_L"] <- out[i + 1, "p_O1_given_H1_L"] + out[i + 1, "p_O2_given_H2_L"]
        out[i + 1, "p_O_H_given_L"] <- out[i + 1, "p_H_given_L"] + out[i + 1, "p_O1_given_H1_L"] + out[i + 1, "p_O2_given_H2_L"]
        ##
        out[i + 1, "p_G_given_O_H_lambda"] <- p_g_given_o_h_lambda
    }
    ## remove constant - who cares
    C1 <- out[, "p_O_H_given_L"]
    a <- max(C1)
    C2 <- C1 - a
    out[, "p_H_given_O_L"] <- exp(C2) / sum(exp(C2))

    p_O_given_L_haploid <- sum(exp(out[, "p_O_H_given_L"]))
    p_O_given_L_haploid_alternate <- sum(exp(out[, "p_H_given_L"]) * exp(out[, "p_O_given_H_L"]))
    expect_equal(p_O_given_L_haploid, p_O_given_L_haploid_alternate)

    p_O_given_L_diploid <- exp(diploid_fit()$ll)
    expect_equal(p_O_given_L_haploid, p_O_given_L_diploid)

    ##
    ## P(O | \lambda)
    ##
    ## OK - NOW - with truly random sampling
    ## this SUCKs - they are all EQUALLY LIKELY
    nRep <- 1e6
    z <- sapply(sample(x = 1:nrow(out), size = nRep, prob = exp(out[, "p_H_given_L"]), replace = TRUE), function(i_H) {
        return(as.numeric(out[i_H, "p_O_given_H_L"]))
    })
    ##
    x <- seq(3, 5, length.out = 20)
    y <- sapply(x, function(a) {
        mean(exp(z[1:round(10 ** a)]))
    })
    plot(x, y)
    abline(h = p_O_given_L_diploid) ## yup visual inspection, looks fine
    a <- max(z); z <- z - a; exp(a) * mean(exp(z))

    ##
    ## P(G_{5} = 0 | O, \lambda)
    ##
    genProbsD <- diploid_fit()$out$genProbs[i_g_pg, i_snp_pg]
    genProbsH <- sum(out[, "p_G_given_O_H_lambda"] * out[, "p_H_given_O_L"])
    expect_equal(genProbsD, genProbsH)

    z <- sapply(sample(x = 1:nrow(out), size = nRep, prob = out[, "p_H_given_O_L"], replace = TRUE), function(i_H) {
        return(as.numeric(out[i_H, "p_G_given_O_H_lambda"]))
    })

    ##
    x <- seq(3, 6, length.out = 20)
    y <- sapply(x, function(a) {
        mean(z[1:round(10 ** a)])
    })
    plot(x, y)
    abline(h = genProbsD) ## yup, works


})





if(1 == 0) {
test_that("properly simplified gibbs-nipt forward and backward", {

    ## OK - promising enough!
    skip("wip")
    ff <- 0.2
    test_package <- make_fb_test_package(
        K = 20,
        nReads = 50,
        nSNPs = 100,
        gridWindowSize = 20,
        method = "triploid-nipt",
        ff = ff,
        seed = 100,
        eHapsMin = 0.01
    )
    grid <- test_package$grid
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_t_H <- test_package$transMatRate_t_H
    alphaMatCurrent_t <- test_package$alphaMatCurrent_t
    priorCurrent <- test_package$priorCurrent
    eHapsCurrent_t <- test_package$eHapsCurrent_t
    KKK <- K * K * K
    sampleReads <- test_package$sampleReads
    eMatHap_t <- test_package$eMatHap_t
    ## check that with no emissions, move one forward keeps sum the same
    true_H <- test_package$true_H
    nReads <- length(sampleReads)
    wif <- sapply(sampleReads, function(x) x[[2]]) + 1 ## 1-based???
    ## so reads should be very easy to tell apart
    ##

    outT <- forwardBackwardTriploid(
        sampleReads = sampleReads,
        priorCurrent = priorCurrent,
        alphaMatCurrent_t = alphaMatCurrent_t,
        eHapsCurrent_t = eHapsCurrent_t,
        transMatRate_t_H = transMatRate_t_H,
        ff = ff,
        grid = grid,
        blocks_for_output = array(0, c(1, 1))
    )
    TgenProbsM_t <- outT$genProbsM_t
    TgenProbsF_t <- outT$genProbsF_t
    ##
    ## OK - if I get gammas out, can I do this
    ## assign reads in a preliminary fashion to these?
    ## can I try to place reads based on these?
    gammaKMT_t <- outT$gammaKMT_t
    gammaKMU_t <- outT$gammaKMU_t
    gammaKP_t <- outT$gammaKP_t

    get_H_estimate <- function(noise = 0) {
        ## OK - so have the marginal gammas
        readProbs <- array(0, c(3, nReads)) ## mt, mu, p
        estimated_H <- array(0, nReads)
        for(iRead in 1:nReads) {
            curGrid <- wif[iRead]
            pMT <- sum(gammaKMT_t[, curGrid] * eMatHap_t[, iRead])
            pMU <- sum(gammaKMU_t[, curGrid] * eMatHap_t[, iRead])
            pP <-  sum(gammaKP_t[, curGrid] * eMatHap_t[, iRead])
            readProbs[, iRead] <- c(pMT, pMU, pP) + noise
            readProbs[, iRead] <- readProbs[, iRead] / sum(readProbs[, iRead])
            ## choose one
            estimated_H[iRead] <- sample(c(1, 2, 3), size = 1, prob = readProbs[, iRead])
        }
        ## ridiculous - table(true_H, estimated_H)
        ## better table(true_H, estimated_H)
        table(estimated_H, true_H) ## not bad
        return(estimated_H)
    }

    N <- 10
        result <- array(0, c(N, 5))
        AgenProbsM_t <- array(0, c(3, ncol(TgenProbsM_t)))
        AgenProbsF_t <- array(0, c(3, ncol(TgenProbsM_t)))
        LL <- array(0, N)
    for(i in 1:N) {
        use_starting_read_labels <- TRUE
        ## estimate
        estimated_H <- get_H_estimate(noise = 0.2) ## ugh, too accurate?
        ## table(estimated_H, true_H)
        ## estimated_H <- true_H
        ## random
        ## how do I nudge?
        out <- rcpp_forwardBackwardGibbsNIPT(
            sampleReads = sampleReads,
            priorCurrent = priorCurrent,
            alphaMatCurrent_t = alphaMatCurrent_t,
            eHapsCurrent_t = eHapsCurrent_t,
            transMatRate_t_H = transMatRate_t_H,
            ff = ff,
            blocks_for_output = array(0, c(1, 1)),
            alphaHat_t1 = array(0, c(1, 1)),
            alphaHat_t2 = array(0, c(1, 1)),
            alphaHat_t3 = array(0, c(1, 1)),
            betaHat_t1 = array(0, c(1, 1)),
            betaHat_t2 = array(0, c(1, 1)),
            betaHat_t3 = array(0, c(1, 1)),
            hapSum_tc = array(0, c(1, 1, 1)),
            grid = grid,
            wif0 = integer(1),
            n_gibbs_burn_in_its = 1,
            n_gibbs_sample_its = 1,
            use_starting_read_labels = use_starting_read_labels,
            return_p_store = TRUE,
            starting_read_labels = estimated_H
        )
    ##
    ## not working very well!
    ## better initialization?
    table(out$H, true_H)
    LL[i] <- log(out$p1) + log(out$p2) + log(out$p3)
    gammaM_t <- out$gammaM_t
    gammaF_t <- out$gammaF_t
    genProbsM_t <- out$genProbsM_t
    genProbsF_t <- out$genProbsF_t
    AgenProbsM_t <- AgenProbsM_t + genProbsM_t
    AgenProbsF_t <- AgenProbsF_t + genProbsF_t
    ##
    ## is the gibbs getting "stuck"? or blocking...
    ##
    dM <- (AgenProbsM_t[2, ] + 2 * AgenProbsM_t[3, ]) / i
    dF <- (AgenProbsF_t[2, ] + 2 * AgenProbsF_t[3, ]) / i
    dM_triploid <- (TgenProbsM_t[2, ] + 2 * TgenProbsM_t[3, ])
    dF_triploid <- (TgenProbsF_t[2, ] + 2 * TgenProbsF_t[3, ])
    dM_truth <- colSums(eHapsCurrent_t[c(1, 2), ])
    dF_truth <- colSums(eHapsCurrent_t[c(1, 3), ])
    ## compare against truth here
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    par(mfrow = c(2, 1))
    ## maternal
    plot(x = 0, y = 0, xlim = c(1, length(grid)), ylim = c(0, 2), main = i)
    lines(dM_triploid, col = "green", lwd = 2)
    lines(dM_truth, col = "black")
    lines(dM, col = "red", lwd = 2)
    plot(x = 0, y = 0, xlim = c(1, length(grid)), ylim = c(0, 2))
    lines(dF_triploid, col = "green", lwd = 2)
    lines(dF_truth, col = "black", lwd = 1)
        lines(dF, col = "red", lwd = 2)
        legend("bottomleft", c("triploid", "truth", "gibbs"), col = c("green", "black", "red"), lwd = 2)
    ## r2
    result[i, ] <- c(
        LL[i],
        mean(abs(dM - dM_truth)),
        mean(abs(dM_triploid - dM_truth)),
        mean(abs(dF - dF_truth)),
        mean(abs(dF_triploid - dF_truth))
    )
        ## OK this is working yay
    Sys.sleep(0.1)
    }




})
}


test_that("gibbs-nipt works with S > 1", {

    S <- 3
    verbose <- FALSE
    ff <- 0.2
    n_gibbs_starts <- 3
    test_package <- make_fb_test_package(
        K = 4,
        nReads = 10,
        nSNPs = 3,
        gridWindowSize = NA,
        method = "triploid-nipt",
        ff = ff,
        seed = 100,
        eHapsMin = 0.01,
        S = S
    )
    grid <- test_package$grid
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    priorCurrent_m <- test_package$priorCurrent_m
    priorCurrent_m[] <- array(runif(prod(dim(priorCurrent_m))), dim(priorCurrent_m))
    priorCurrent_m <- t(t(priorCurrent_m) / rowSums(t(priorCurrent_m)))
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    ## it is probably better if fit is much poorer
    eHapsCurrent_tc[] <- array(runif(prod(dim(eHapsCurrent_tc))), dim(eHapsCurrent_tc))
    KKK <- K * K * K
    sampleReads <- test_package$sampleReads
    ## check that with no emissions, move one forward keeps sum the same
    true_H <- test_package$true_H[1]
    nReads <- length(sampleReads)
    wif <- sapply(sampleReads, function(x) x[[2]]) + 1

    double_list_of_starting_read_labels <- initialize_read_labels_from_random(
        method = "gibbs-nipt",
        ffs = ff,
        iSample = 1,
        n_gibbs_starts = n_gibbs_starts,
        sampleReads = sampleReads,
        S = S
    )
    ## needs 6 entries
    ## 3 gibbs samplings and S is 2
    seed_vector <- 1:(S * n_gibbs_starts) * 10
    
    set.seed(123)
    out_auto <- run_forward_backwards(
        method = "gibbs-nipt",
        sampleReads = sampleReads,
        priorCurrent_m = priorCurrent_m,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        eHapsCurrent_tc = eHapsCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        ff = ff,
        return_genProbs = TRUE,
        double_list_of_starting_read_labels = double_list_of_starting_read_labels,
        n_gibbs_starts = n_gibbs_starts,
        grid = grid,
        seed_vector = seed_vector,
        return_p_store = TRUE
    )
    out_manual <- lapply(1:S, function(s) {
        set.seed(123)
        run_forward_backwards(
            method = "gibbs-nipt",
            sampleReads = sampleReads,
            priorCurrent_m = priorCurrent_m[, s, drop = FALSE],
            alphaMatCurrent_tc = alphaMatCurrent_tc[, , s, drop = FALSE],
            eHapsCurrent_tc = eHapsCurrent_tc[, , s, drop = FALSE],
            transMatRate_tc_H = transMatRate_tc_H[, , s, drop = FALSE],
            ff = ff,
            return_genProbs = TRUE,
            seed_vector = seed_vector[S * (s - 1) + 1:n_gibbs_starts],
            double_list_of_starting_read_labels = double_list_of_starting_read_labels[s],
            n_gibbs_starts = n_gibbs_starts,
            grid = grid,
            suppressOutput = 1,
            return_p_store = TRUE
        )
    })

    ## check ending read labels
    for(s in 1:S) {
        expect_equal(
            out_auto[[1]]$double_list_of_ending_read_labels[[s]],
            out_manual[[s]][[1]]$double_list_of_ending_read_labels[[1]]
        )
    }
    
    manual_genProbsM_t <- out_manual[[1]][[1]]$genProbsM_t
    manual_genProbsF_t <- out_manual[[1]][[1]]$genProbsF_t    
    for(s in 2:S) {
        manual_genProbsM_t <- manual_genProbsM_t + out_manual[[s]][[1]]$genProbsM_t
        manual_genProbsF_t <- manual_genProbsF_t + out_manual[[s]][[1]]$genProbsF_t        
    }
    manual_genProbsM_t <- manual_genProbsM_t / S
    manual_genProbsF_t <- manual_genProbsF_t / S    
    ##
    expect_equal(
        out_auto[[1]]$genProbsM_t,
        manual_genProbsM_t
    )
    expect_equal(
        out_auto[[1]]$genProbsF_t,
        manual_genProbsF_t
    )


})


