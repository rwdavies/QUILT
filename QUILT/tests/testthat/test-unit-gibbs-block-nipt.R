if (1 == 0) {

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
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))



}



test_that("can do get_log_p_H_class in R and Rcpp", {

    H_class <- sample(0:7, 1000, replace = TRUE)
    ff <- 0.2
    
    expect_equal(
        rcpp_get_log_p_H_class(H_class, ff),
        get_log_p_H_class(H_class, ff)
    )

})


test_that("can do get_log_p_H_class2 in R and Rcpp", {

    expect_equal(
        get_log_p_H_class2(5, 10, 15, 20, 25, 30, 0.2),
        rcpp_get_log_p_H_class2(5, 10, 15, 20, 25, 30, 0.2)
    )

})



test_that("can confirm that can calculate read label probabilities using H class in R and Rcpp", {

    H_class <- sample(0:7, 1000, replace = TRUE)
    read_start_0_based <- 100
    read_end_0_based <- 200
    
    ff <- 0.2
    prior_probs <- c(0.5, (1 - ff) / 2, ff / 2)
    log_prior_probs <- log(prior_probs)
    rr <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(2, 3, 1),
        c(3, 1, 2),
        c(3, 2, 1)
    )
    rr <- matrix(as.integer(rr), ncol = 3)

    choice_log_probs_H <- calculate_block_read_label_probabilities_using_H_class(
        read_start_0_based = read_start_0_based,
        read_end_0_based = read_end_0_based,
        H_class = H_class,
        log_prior_probs = log_prior_probs,
        rr = rr,
        ff = ff
    )

    rcpp_choice_log_probs_H <- rcpp_calculate_block_read_label_probabilities_using_H_class(
        read_start_0_based = read_start_0_based,
        read_end_0_based = read_end_0_based,
        H_class = H_class,
        log_prior_probs = log_prior_probs,
        rr = rr,
        ff = ff
    )

    expect_equal(
        choice_log_probs_H,
        rcpp_choice_log_probs_H
    )

})



test_that("idea of one_based_swap can carry over between R and Rcpp", {

    rx <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(3, 1, 2),
        c(2, 3, 1),
        c(3, 2, 1)
    )
    
    for(ir_chosen in 1:6) {

        set.seed(119911)
        H <- sample(1:3, 100, replace = TRUE)
        H_class <- sample(0:7, 100, replace = TRUE)

        one_based_swap <- c(1, 1 + rx[ir_chosen, 1:3], 8 - rx[ir_chosen, 3:1], 8)

        H_new <- integer(100)
        H_class_new <- integer(100)

        ## in R
        for(iRead in 0:99) {
            H_class_new[iRead + 1] <- one_based_swap[H_class[iRead + 1] + 1] - 1
            H_new[iRead + 1] <- one_based_swap[H[iRead + 1] + 1] - 1L
        }

        ## in Rcpp (test function to confirm logic)
        ## these will OVERWRITE the original versions
        ## hence the tests below
        rcpp_test_one_based_swap(rx, ir_chosen, H, H_class)

        expect_equal(H, H_new)
        expect_equal(H_class, H_class_new)

    }

})

test_that("can do single quantile hack in cpp using arma", {
    x <- runif(1000)
    q <- 0.90
    r1 <- rcpp_simple_quantile(x, q)
    r2 <- x[order(x)[as.integer(q * length(x)) + 1]]
    expect_equal(r1, r2)
})


test_that("can make rlcM", {

    for(ff in c(0.2, 0, 1, 0.5)) {
        
        rlc <- make_rlc(ff)
        rr <- rbind(
            c(1, 2, 3),
            c(1, 3, 2),
            c(2, 1, 3),
            c(2, 3, 1),
            c(3, 1, 2),
            c(3, 2, 1)
        )
        rr0 <- rr - 1
        rlcM <- array(-1, c(3, 6, 7))
        rlcM <- R_fill_rlcM(rlcM, rlc, rr)
        
        rlcM2 <- array(-1, c(3, 6, 7))    
        Rcpp_fill_rlcM(rlcM2, rlc, rr0)
        
         expect_equal(rlcM, rlcM2)

    }
        
})


test_that("can reconfigure blocked_snps", {

    set.seed(9919)
    nSNPs <- 100
    L <- sort(sample(nSNPs * 10, nSNPs))
    gridWindowSize <- 10
    out <- assign_positions_to_grid(
        L = L,
        gridWindowSize = gridWindowSize
    )
    grid <- out$grid
    nGrids <- out[["nGrids"]]
    nReads <- 200

    ## make up blocked_snps
    blocked_snps <- array(6, nSNPs)
    blocked_snps[0 <= grid & grid <= 10] <- 0 ## no reads
    blocked_snps[11 <= grid & grid <= 15] <- 1
    blocked_snps[16 <= grid & grid <= 20] <- 2
    blocked_snps[21 <= grid & grid <= 22] <- 3 ## no reads
    blocked_snps[23 <= grid & grid <= 24] <- 4 ## no reads
    blocked_snps[25 <= grid & grid <= 40] <- 5
    blocked_snps[41 <= grid & grid <= 50] <- 6 ## no reads
    blocked_snps[51 <= grid & grid <= 60] <- 7
    blocked_snps[61 <= grid] <- 8 ## no reads
    
    language <- "R"

    ## sample some reads, include two consecutive that are NA, and a single
    wif0 <- sort(sample(setdiff(grid, c(0:10, 11:15, 21:24)), 200, replace = TRUE))
    languages <- c("R", "Rcpp", "R")
    do_checks <- c(FALSE, FALSE, TRUE)
    i_option <- 1
    verbose <- FALSE
    
    out <- lapply(1:3, function(i_option) {

        language <- languages[i_option]
        if (verbose) {
            print(paste0("---------- ", language, " ---------"))
        }
        
        do_check <- do_checks[i_option]
        do_removal <- !do_check
        if (language == "R") {
            f <- R_make_gibbs_considers
        } else if (language == "Rcpp") {
            f <- Rcpp_make_gibbs_considers
        }
        
        out <- f(
            blocked_snps = blocked_snps,
            grid = grid,
            wif0 = wif0,
            nGrids = nGrids,
            do_removal = do_removal
        )
        consider_reads_start_0_based <- out[["consider_reads_start_0_based"]]
        consider_reads_end_0_based  <- out[["consider_reads_end_0_based"]]
        consider_grid_start_0_based <- out[["consider_grid_start_0_based"]]
        consider_grid_end_0_based <- out[["consider_grid_end_0_based"]]
        consider_snp_start_0_based <- out[["consider_snp_start_0_based"]]
        consider_snp_end_0_based <- out[["consider_snp_end_0_based"]]
        consider_grid_where_0_based <- out[["consider_grid_where_0_based"]]
        n_blocks <- out[["n_blocks"]]

        ##
        ## build everything easy (slow, R) way
        ##
        if (do_check) {
            check_consider_snp_start_0_based <- as.integer(sapply(0:(n_blocks - 1), function(i) {    min(which(blocked_snps == i))}) - 1)
            check_consider_snp_end_0_based <- as.integer(sapply(0:(n_blocks - 1), function(i) {    max(which(blocked_snps == i))}) - 1)        
            check_consider_grid_start_0_based <- grid[check_consider_snp_start_0_based + 1]
            check_consider_grid_end_0_based <- grid[check_consider_snp_end_0_based + 1]        
            n_blocks <- max(blocked_snps + 1)
            check_consider_reads_start_0_based <- integer(n_blocks)
            check_consider_reads_end_0_based <- integer(n_blocks)
            for(i in 1:n_blocks) {
                a <- check_consider_grid_start_0_based[i]
                b <- check_consider_grid_end_0_based[i]
                w <- wif0 %in% (a:b)
                if (sum(w) == 0) {
                    ## later, re-ject this block!
                    check_consider_reads_start_0_based[i] <- -1
                    check_consider_reads_end_0_based[i] <- -1          
                } else {
                    check_consider_reads_start_0_based[i] <- min(which(w)) - 1
                    check_consider_reads_end_0_based[i] <- max(which(w)) - 1
                }
            }
            ##
            expect_equal(check_consider_reads_start_0_based, consider_reads_start_0_based)
            expect_equal(check_consider_reads_end_0_based, consider_reads_end_0_based)
            ##
            expect_equal(check_consider_grid_start_0_based, consider_grid_start_0_based)
            expect_equal(check_consider_grid_end_0_based, consider_grid_end_0_based)
            ##
            expect_equal(check_consider_snp_start_0_based, consider_snp_start_0_based)
            expect_equal(check_consider_snp_end_0_based, consider_snp_end_0_based)
        }

        ##
        ## check grids for no spaces
        ##
        d <- (consider_grid_start_0_based[-1] - 1) - (consider_grid_end_0_based[-n_blocks])
        expect_equal(sum(abs(d)), 0)
        d <- (consider_snp_start_0_based[-1] - 1) - (consider_snp_end_0_based[-n_blocks])
        expect_equal(sum(abs(d)), 0)

        return(out)

    })

    ## check same between languages
    for(what in c("consider_reads_start_0_based", "consider_reads_end_0_based", "consider_grid_start_0_based", "consider_grid_end_0_based", "consider_snp_start_0_based", "consider_snp_end_0_based", "consider_grid_where_0_based", "n_blocks")) {
        ## print(what)
        ## print(out[[1]][[what]]) ## R
        ## print(out[[2]][[what]]) ## Rcpp
        if (verbose) {
            print(what)
        }
        expect_equal(out[[1]][[what]], out[[2]][[what]])
    }

    check_non_overlapping_and_span <- function(x, y, n) {
        check <- array(NA, n)
        for(i in 1:length(x)) {
            w <- (x[i] + 1):(y[i] + 1)
            expect_equal(sum(!is.na(check[w])), 0)
            check[w] <- i
        }
        expect_equal(sum(is.na(check)), 0)
        return(NULL)
    }
    check_non_overlapping_and_span(out[[1]][["consider_reads_start_0_based"]], out[[1]][["consider_reads_end_0_based"]], nReads)
    check_non_overlapping_and_span(out[[1]][["consider_grid_start_0_based"]], out[[1]][["consider_grid_end_0_based"]], nGrids)
    check_non_overlapping_and_span(out[[1]][["consider_snp_start_0_based"]], out[[1]][["consider_snp_end_0_based"]], nSNPs)

})



test_that("can perform block gibbs when ff is 0", {

    set.seed(199)
    do_checks <- TRUE
    verbose <- FALSE
    s <- 1
    S <- 1
    ff <- 0
    n_gibbs_starts <- 3
    nSNPs <- 200
    gridWindowSize <- 2000
    L <- sort(sample(100000, nSNPs))
    shuffle_bin_radius <- 2000
    class_sum_cutoff <- 0.06
    ##
    test_package <- make_quilt_fb_test_package(
        K = 4,
        nReads = nSNPs * 4,
        nSNPs = nSNPs,
        gridWindowSize = gridWindowSize,
        method = "triploid-nipt",
        ff = ff,
        seed = 101,
        S = 1,
        L = L,
        eHapsMin = 0.01, ## make more confident
        bq_mult = 40,
        randomize_sample_read_length = TRUE,
        return_eMatGridTri_t = FALSE
    )
    grid <- test_package$grid
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    ## for here, these are the same (no prior)
    alphaMatCurrent_tc[] <- 1 / K
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    true_H <- as.integer(test_package$true_H) ## always keep as 1-based integer, EVEN in cpp
    maxDifferenceBetweenReads <- 1000
    Jmax <- 1000
    grid_distances <- test_package$grid_distances
    L_grid <- test_package$L_grid
    eMatRead_t <- test_package$list_of_eMatRead_t[[1]]    
    ##
    sampleReads <- test_package$sampleReads
    ## re-size the reads to feature more random length
    wif0 <- as.integer(sapply(sampleReads, function(x) x[[2]]))
    grid_has_read <- rep(FALSE, nGrids)
    grid_has_read[wif0 + 1] <- TRUE
    
    nReads <- length(sampleReads)
    ##
    rr <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(2, 3, 1),
        c(3, 1, 2),
        c(3, 2, 1)
    )
    rr <- matrix(as.integer(rr), ncol = 3)
    rr0 <- rr - 1L    
    ir_test <- 3
    rlc <- make_rlc(ff)
    ## 
    initial_package <- for_testing_get_full_package_probabilities(true_H, list( eHapsCurrent_tc = eHapsCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, alphaMatCurrent_tc = alphaMatCurrent_tc,priorCurrent_m = priorCurrent_m ,eMatRead_t = eMatRead_t,s = s,sampleReads = sampleReads))
    ##
    ##
    ##
    true_H_class <- calculate_H_class(
        eMatRead_t = eMatRead_t,
        alphaHat_t1 = initial_package[[1]][["alphaHat_t"]],
        alphaHat_t2 = initial_package[[2]][["alphaHat_t"]],
        alphaHat_t3 = initial_package[[3]][["alphaHat_t"]],
        betaHat_t1 = initial_package[[1]][["betaHat_t"]],
        betaHat_t2 = initial_package[[2]][["betaHat_t"]],
        betaHat_t3 = initial_package[[3]][["betaHat_t"]],
        ff = ff,
        wif0 = wif0,
        H = true_H,
        class_sum_cutoff = class_sum_cutoff
    )
    H <- integer(length(true_H))
    H[] <- true_H[]
    blocked_snps <- integer(nSNPs)
    blocked_snps[grid <= 20] <- 0L
    blocked_snps[20 < grid] <- 1L
    languages_to_test <- c("R", "R_with_Rcpp", "Rcpp")
    block_approach <- 1
    block_gibbs_quantile_prob <- 0.9
    ## language <- "R_with_Rcpp"
    ## language <- "Rcpp"    
    verbose <- FALSE
    use_smooth_cm_in_block_gibbs <- FALSE
    smooth_cm <- numeric(1)
    
    for(block_approach in c(1, 4)) {

        if (verbose) {
            print(paste0("---------- block_approach = ", block_approach, " ----------"))
        }
        ## out <- mclapply(languages_to_test, mc.cores = 1, function(language) {
        for(language in languages_to_test) {
            if (verbose) {
                print(paste0("------------", language, "------------"))
            }
            block_out <- helper_block_gibbs_resampler(
                H = H,
                blocked_snps = blocked_snps,
                grid = grid,
                wif0 = wif0,
                grid_has_read = grid_has_read,
                ff = ff,
                s = s,
                eHapsCurrent_tc = eHapsCurrent_tc,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                priorCurrent_m = priorCurrent_m,
                transMatRate_tc_H = transMatRate_tc_H,
                sampleReads = sampleReads,
                L_grid = L_grid,
                shuffle_bin_radius = shuffle_bin_radius,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                Jmax = Jmax,
                do_checks = do_checks,
                language = language,
                verbose = verbose,
                block_approach = block_approach,
                class_sum_cutoff = class_sum_cutoff,
                block_gibbs_quantile_prob = block_gibbs_quantile_prob,
                smooth_cm = smooth_cm,
                use_smooth_cm_in_block_gibbs = use_smooth_cm_in_block_gibbs
            )
            ## expect all jumps to be in 1 or 3
            expect_equal(0, sum(!(block_out$block_results[, "ir_chosen"] %in% c(0, 1, 3)))) ## 0 for skipped
        }
       
    }


})















test_that("can perform block gibbs", {

    ## run a few times to converge
    ## artificially break them
    ## then run!
    set.seed(11001)
    do_checks <- TRUE
    verbose <- FALSE
    s <- 1
    S <- 1
    ff <- 0.25
    n_gibbs_starts <- 3
    nSNPs <- 200
    gridWindowSize <- 2000
    L <- sort(sample(100000, nSNPs))
    shuffle_bin_radius <- 2000
    class_sum_cutoff <- 0.06
    ##
    test_package <- make_quilt_fb_test_package(
        K = 4,
        nReads = nSNPs * 4,
        nSNPs = nSNPs,
        gridWindowSize = gridWindowSize,
        method = "triploid-nipt",
        ff = ff,
        seed = 101,
        S = 1,
        L = L,
        eHapsMin = 0.01, ## make more confident
        bq_mult = 40,
        randomize_sample_read_length = TRUE,
        return_eMatGridTri_t = FALSE
    )
    
    grid <- test_package$grid
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    alphaMatCurrent_tc[] <- 1 / K
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    true_H <- as.integer(test_package$true_H) ## always keep as 1-based integer, EVEN in cpp
    maxDifferenceBetweenReads <- 1000
    Jmax <- 1000
    grid_distances <- test_package$grid_distances
    L_grid <- test_package$L_grid
    eMatRead_t <- test_package$list_of_eMatRead_t[[1]]    
    ##
    block_gibbs_quantile_prob <- 0.9
    sampleReads <- test_package$sampleReads
    ## re-size the reads to feature more random length
    wif0 <- as.integer(sapply(sampleReads, function(x) x[[2]]))
    grid_has_read <- rep(FALSE, nGrids)
    grid_has_read[wif0 + 1] <- TRUE
    nReads <- length(sampleReads)
    ##
    rr <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(2, 3, 1),
        c(3, 1, 2),
        c(3, 2, 1)
    )
    rr <- matrix(as.integer(rr), ncol = 3)
    rr0 <- rr - 1L    
    rlc <- make_rlc(ff)
    rx <- rbind(
        c(1L, 2L, 3L),
        c(1L, 3L, 2L),
        c(2L, 1L, 3L),
        c(3L, 1L, 2L),
        c(2L, 3L, 1L),
        c(3L, 2L, 1L)
    )
    ## 
    initial_package <- for_testing_get_full_package_probabilities(true_H, list( eHapsCurrent_tc = eHapsCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, alphaMatCurrent_tc = alphaMatCurrent_tc,priorCurrent_m = priorCurrent_m ,eMatRead_t = eMatRead_t,s = s,sampleReads = sampleReads))
    ##
    ##
    ##
    true_H_class <- calculate_H_class(
        eMatRead_t = eMatRead_t,
        alphaHat_t1 = initial_package[[1]][["alphaHat_t"]],
        alphaHat_t2 = initial_package[[2]][["alphaHat_t"]],
        alphaHat_t3 = initial_package[[3]][["alphaHat_t"]],
        betaHat_t1 = initial_package[[1]][["betaHat_t"]],
        betaHat_t2 = initial_package[[2]][["betaHat_t"]],
        betaHat_t3 = initial_package[[3]][["betaHat_t"]],
        ff = ff,
        wif0 = wif0,
        H = true_H,
        class_sum_cutoff = 1
    )
    use_smooth_cm_in_block_gibbs <- FALSE
    smooth_cm <- numeric(1)

    

    ## this is supposed to be like real data
    ## imagine some small blocks off in the middle of lots of blocksn
    ## want some small block in the middle to test on
    ir_test_options <- c(1:6, 1:6, 1:6) ## this is the ir_chosen to flip
    ## 1 = original
    ## 2 = do not count some H
    ## 4 = sample H for all options with probabilities
    ## 6 = sample H_class (?then possibly re-sample H? or do organically?)
    block_approach_options <- rep(c(1, 4, 6), each = 6)
    to_run <- 1:length(ir_test_options)
    ir_test <- 5

    ## for now - no chance it works first time
    to_run <- c(1:6, 13:18) ## block_approach 1 and -6
    ## to_run <- 17
    
    ## i_option <- 14
    
    ## weirdly runs into problem when not multi-cored
    ## some very weird stuff going on!
    ## , mc.cores = length(to_run),
    ## out <- lapply(to_run, function(i_option) {
    
    for(i_option in to_run) {
        
        ## 
        verbose <- FALSE
        do_checks <- FALSE
        ir_test <- ir_test_options[i_option]
        block_approach <- block_approach_options[i_option]
        if (verbose) {
            print(paste0("---------- ir_test = ", ir_test, " ----------"))
            print(paste0("---------- block_approach = ", block_approach, " ----------"))
        }
        ##
        ## recast blocked_snps
        ## 
        languages_to_test <- c("R", "R_with_Rcpp", "Rcpp")
        language <- "R"
        
        out <- mclapply(languages_to_test, mc.cores = 1, function(language) {

            if (verbose) {
                print(paste0("------------", language, "------------"))
            }
            set.seed(4)
            do_checks <- FALSE
            ## if (language == "R_with_Rcpp") {
            ##     do_checks <- FALSE
            ## }
            ## OK, am here
            ## ir_test = 3 is failing once Rcpp runs
            ## fuck it, for simplicity, make it first block that's messed up
            H <- integer(length(true_H))
            H[] <- true_H[]
            n <- 2
            blocked_snps <- floor(n * grid / (max(grid) + 0.1))
            ## block_to_consider <- ((n / 2) - 1)
            block_to_consider <- 0
            a <- (((grid)[blocked_snps == block_to_consider]))
            range_to_flip <- min(a):max(a)
            ##

            one_based_swap <- as.integer(c(1, 1 + rr[ir_test, 1:3], 8 - rr[ir_test, 3:1], 8))
            ## re-sample with respect to rr and true_H
            for(iRead in 1:nReads) {
                ## argh, I'm not sure I've done this right
                hci <- true_H_class[iRead]
                if (hci == 0) {
                    hci <- 7
                }
                hco <- one_based_swap[hci + 1] - 1
                ## now sample according to that new class
                hi <- as.integer(sample(c(1L, 2L, 3L), size = 1, prob = rlc[hci, ]))
                ho <- as.integer(sample(c(1L, 2L, 3L), size = 1, prob = rlc[hco, ]))
                ## 
                if (wif0[iRead] %in% range_to_flip) {
                    H[iRead] <- ho
                } else {
                    H[iRead] <- hi
                }
                ## ## so based on the class, determine the new probabilities
                ## if (hlc == 0) {
                ##     ##h <- sample(c(1L, 2L, 3L), prob = rlc[7, ])
                ##     hlc <- 7 ## same thing
                ## }
                ## ## sample a label based on those probabilities
                ## hl <- as.integer(sample(c(1L, 2L, 3L), size = 1, prob = rlc[hlc, ]))
                ## 
                ## if (hlc == 0) {
                ##     hl <- rr[ir_test, h]
                ## } else {
                ##     hl <- as.integer(sample(c(1L, 2L, 3L), size = 1, prob = rlc[h, ]))
                ## }
            }

            ## ## check that this makes sense
            ## initial_package <- for_testing_get_full_package_probabilities(H, list( eHapsCurrent_tc = eHapsCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, alphaMatCurrent_tc = alphaMatCurrent_tc,priorCurrent_m = priorCurrent_m ,eMatRead_t = eMatRead_t,s = s,sampleReads = sampleReads))
            ## H_class <- calculate_H_class(
            ##     eMatRead_t = eMatRead_t,
            ##     alphaHat_t1 = initial_package[[1]][["alphaHat_t"]],
            ##     alphaHat_t2 = initial_package[[2]][["alphaHat_t"]],
            ##     alphaHat_t3 = initial_package[[3]][["alphaHat_t"]],
            ##     betaHat_t1 = initial_package[[1]][["betaHat_t"]],
            ##     betaHat_t2 = initial_package[[2]][["betaHat_t"]],
            ##     betaHat_t3 = initial_package[[3]][["betaHat_t"]],
            ##     ff = ff,
            ##     wif0 = wif0,
            ##     H = H,
            ##     class_sum_cutoff = 0.06
            ## )
            ## ## 
            ## w <- wif0 %in% range_to_flip
            ## table(H_class[w], H[w])
            ## table(true_H_class[w], H_class[w])            
            ## prior_probs <- c(0.5, (1 - ff) / 2, ff / 2)
            ## log_prior_probs <- log(prior_probs)
            ## calculate_block_read_label_probabilities_using_H_class(0, 420, H_class,log_prior_probs, rr, ff)
            ## ## happiest with 4? what?
            
            ## now, check class OK
            check_agreements_with_truth_class(true_H, true_H_class, H, rlc)
            check_agreements_with_truth_class(true_H, true_H_class, true_H, rlc)
            
            ## determine H_class
            block_out <- helper_block_gibbs_resampler(
                H = H,
                blocked_snps = blocked_snps,
                grid = grid,
                wif0 = wif0,
                grid_has_read = grid_has_read,
                ff = ff,
                s = s,
                eHapsCurrent_tc = eHapsCurrent_tc,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                priorCurrent_m = priorCurrent_m,
                transMatRate_tc_H = transMatRate_tc_H,
                sampleReads = sampleReads,
                L_grid = L_grid,
                shuffle_bin_radius = shuffle_bin_radius,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                Jmax = Jmax,
                do_checks = do_checks,
                language = language,
                verbose = FALSE,
                block_approach = block_approach,
                class_sum_cutoff = 1,
                smooth_cm = smooth_cm,
                use_smooth_cm_in_block_gibbs = use_smooth_cm_in_block_gibbs,
                block_gibbs_quantile_prob = block_gibbs_quantile_prob                
            )
            block_out$block_results[c(1, 3), ]
            
            ## AM HERE
            ## FOR ir_test = 5, the first block is not suggesting what I know it should do
            ## can I fix this after lunch?
            
            
            ## print(block_out$block_results[, c("p1", "p2", "ir_chosen")])
            
            ## check H (class) here
            x <- check_agreements_with_truth_class(true_H, true_H_class, block_out[["H"]], rlc)
            expect_equal(as.logical(x < 7), TRUE) ## weird?

            ## print(block_out$block_results[2, 1:4])
            ## expect new results pretty good
            if (verbose) {
                print(paste0("language=", language))
                print(table(H, true_H))
                print(table(block_out[["H"]], true_H)) ## won't be "exact", will be close!
                print(block_out[["block_results"]])
            }
            ## expect_equal(sum(abs((block_out[["H"]] - true_H) > 0)) < 5, TRUE)
            ## expect_equal(as.integer(block_out[["block_results"]][2, "ir_chosen"]), ir_test)
            return(block_out)
        })
        check_mclapply_OK(out)
        names(out) <- languages_to_test

        if (length(out) > 1) {
            for(j in 2:length(out)) {
                expect_equal(out[[1]]$H, out[[j]]$H)
                expect_equal(out[[1]]$block_results, out[[j]]$block_results)
                ## for(what in c("alphaStore", "log_cStore", "alphaStore2", "log_cStore2")) {
                ## for(what in c("alphaStore")) {
                ##     print(paste0("------------", what, "----------"))
                ##     expect_equal(out[[1]][[what]], out[[j]][[what]])
                ## }
                for(what in c("alphaHat_t", "betaHat_t", "c", "eMatGrid_t")) {
                    for(i in 1:3) {
                        expect_equal(out[[1]][[paste0(what, i)]], out[[j]][[paste0(what, i)]])
                    }
                }
            }
        }
        1
       
    }

    ## so in this one, not updating H_class properly?
    ## })
    ## check_mclapply_OK(out)

})




test_that("can perform ff shard block gibbs including at all sites", {

    for(ff in c(0.2)) {    

        set.seed(199)
        do_checks <- TRUE
        verbose <- FALSE
        s <- 1
        S <- 1
        n_gibbs_starts <- 1
        nSNPs <- 1000
        gridWindowSize <- 2000
        L <- sort(sample(200000, nSNPs))
        shuffle_bin_radius <- 2000
        class_sum_cutoff <- 0.06
        ##
        test_package <- make_quilt_fb_test_package(
            K = 4,
            nReads = nSNPs * 4,
            nSNPs = nSNPs,
            gridWindowSize = gridWindowSize,
            method = "triploid-nipt",
            ff = ff,
            seed = 101,
            S = 1,
            L = L,
            eHapsMin = 0.01, ## make more confident
            bq_mult = 40,
            randomize_sample_read_length = TRUE,
            return_eMatGridTri_t = FALSE
        )
        grid <- test_package$grid
        nSNPs <- test_package$nSNPs
        nGrids <- test_package$nGrids
        K <- test_package$K
        transMatRate_tc_H <- test_package$transMatRate_tc_H
        alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
        priorCurrent_m <- test_package$priorCurrent_m
        eHapsCurrent_tc <- test_package$eHapsCurrent_tc
        true_H <- as.integer(test_package$true_H) ## always keep as 1-based integer, EVEN in cpp
        maxDifferenceBetweenReads <- 1000
        maxEmissionMatrixDifference <- 1e6
        Jmax <- 1000
        grid_distances <- test_package$grid_distances
        L_grid <- test_package$L_grid
        eMatRead_t <- test_package$list_of_eMatRead_t[[1]]
        
        ##
        sampleReads <- test_package$sampleReads
        ## re-size the reads to feature more random length
        wif0 <- as.integer(sapply(sampleReads, function(x) x[[2]]))
        nReads <- length(sampleReads)
        ##
        initial_package <- for_testing_get_full_package_probabilities(true_H, list( eHapsCurrent_tc = eHapsCurrent_tc, transMatRate_tc_H = transMatRate_tc_H, alphaMatCurrent_tc = alphaMatCurrent_tc,priorCurrent_m = priorCurrent_m ,eMatRead_t = eMatRead_t,s = s,sampleReads = sampleReads))

        ## also get H_class here
        
        ##
        ## so here make a 10 blocks
        ## real defined below
        ## 
        H <- integer(length(true_H))
        H[] <- true_H[] ## make clone not copy by reference
        ## am here, can I get H_class somehow, don't I have code for this
        blocked_snps <- integer(nSNPs)
        for(i in 0:9) {
            blocked_snps[(10 * i <= grid) & (grid < (10 * (i + 1)))] <- as.integer(i)
        }
        ## 
        rr <- rbind(
            c(1, 2, 3),
            c(1, 3, 2),
            c(2, 1, 3),
            c(2, 3, 1),
            c(3, 1, 2),
            c(3, 2, 1)
        )
        rr <- matrix(as.integer(rr), ncol = 3)
        if (ff == 0) {
            irs <- c(1, 3, 1, 3)
        } else {
            irs <- c(1, 2, 4, 6)
        }
        ## introduce artificial splits
        ir <- irs[1]
        w <- (wif0 >= 0) & (wif0 <= 19)
        H[w] <- rr[ir, H[w]]
        ##
        ir <- irs[2]
        w <- (wif0 >= 20) & (wif0 <= 29)
        H[w] <- rr[ir, H[w]]
        ## 
        ir <- irs[3]
        w <- (wif0 >= 30) & (wif0 <= 59)
        H[w] <- rr[ir, H[w]]
        ##
        ir <- irs[4]
        w <- (wif0 >= 60) & (wif0 <= 99)
        H[w] <- rr[ir, H[w]]
        languages_to_test <- c("R", "R_with_Rcpp", "Rcpp")
        block_approach <- 1
        block_gibbs_quantile_prob <- 0.9
        ## language <- "R_with_Rcpp"
        ## language <- "Rcpp"    
        use_smooth_cm_in_block_gibbs <- FALSE
        smooth_cm <- numeric(1)
        ## 
        ## do not need block approaches
        ## 
        fpp_stuff <- list(
            transMatRate_tc_H = transMatRate_tc_H,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            priorCurrent_m = priorCurrent_m ,
            eMatRead_t = eMatRead_t,
            s = s,
            sampleReads = sampleReads
        )
        initial_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)
        alphaHat_t1 <- initial_package[[1]][["alphaHat_t"]]
        alphaHat_t2 <- initial_package[[2]][["alphaHat_t"]]
        alphaHat_t3 <- initial_package[[3]][["alphaHat_t"]]
        betaHat_t1 <- initial_package[[1]][["betaHat_t"]]
        betaHat_t2 <- initial_package[[2]][["betaHat_t"]]
        betaHat_t3 <- initial_package[[3]][["betaHat_t"]]
        c1 <- initial_package[[1]][["c"]]
        c2 <- initial_package[[2]][["c"]]
        c3 <- initial_package[[3]][["c"]]
        eMatGrid_t1 <- initial_package[[1]][["eMatGrid_t"]]
        eMatGrid_t2 <- initial_package[[2]][["eMatGrid_t"]]
        eMatGrid_t3 <- initial_package[[3]][["eMatGrid_t"]]

        H_class <- calculate_H_class(
            eMatRead_t = eMatRead_t,
            alphaHat_t1 = initial_package[[1]][["alphaHat_t"]],
            alphaHat_t2 = initial_package[[2]][["alphaHat_t"]],
            alphaHat_t3 = initial_package[[3]][["alphaHat_t"]],
            betaHat_t1 = initial_package[[1]][["betaHat_t"]],
            betaHat_t2 = initial_package[[2]][["betaHat_t"]],
            betaHat_t3 = initial_package[[3]][["betaHat_t"]],
            ff = ff,
            wif0 = wif0,
            H = H,
            class_sum_cutoff = class_sum_cutoff
        )
        
        ##
        ##
        ##
        ori_H <- H
        verbose <- FALSE
        i_run <- 3
        
        for(i_run in 1:4) {

            ##print("")
            ##print(paste0("------------"))        
            ##print(paste0("  ff = ", ff, ", i_run = ", i_run))
            ##print(paste0("------------"))
            ##print("")
            
            if (i_run == 1) {
                language  <- "R"
                f <- R_shard_block_gibbs_resampler
                s <- 1
                shard_check_every_pair <- FALSE
                do_checks <- TRUE
            } else if (i_run == 2) {
                language <- "Rcpp"
                f <- Rcpp_shard_block_gibbs_resampler
                s <- 0
                shard_check_every_pair <- FALSE
                do_checks <- FALSE
            } else if (i_run == 3) {
                language  <- "R"
                f <- R_shard_block_gibbs_resampler
                s <- 1
                shard_check_every_pair <- TRUE
                do_checks <- FALSE            
            } else if (i_run == 4) {
                language <- "Rcpp"
                f <- Rcpp_shard_block_gibbs_resampler
                s <- 0
                shard_check_every_pair <- TRUE
                do_checks <- FALSE
            }
            if (verbose) {
                print(paste0("=======language = ", language, "========="))
            }

            set.seed(10)
            ## re-init
            H <- ori_H
            H[1] <- 2L * H[1]
            H[1] <- H[1] / 2L
            H <- as.integer(H) ## ARGH
            initial_package <- for_testing_get_full_package_probabilities(H, fpp_stuff)
            alphaHat_t1 <- initial_package[[1]][["alphaHat_t"]]
            alphaHat_t2 <- initial_package[[2]][["alphaHat_t"]]
            alphaHat_t3 <- initial_package[[3]][["alphaHat_t"]]
            betaHat_t1 <- initial_package[[1]][["betaHat_t"]]
            betaHat_t2 <- initial_package[[2]][["betaHat_t"]]
            betaHat_t3 <- initial_package[[3]][["betaHat_t"]]
            c1 <- initial_package[[1]][["c"]]
            c2 <- initial_package[[2]][["c"]]
            c3 <- initial_package[[3]][["c"]]
            eMatGrid_t1 <- initial_package[[1]][["eMatGrid_t"]]
            eMatGrid_t2 <- initial_package[[2]][["eMatGrid_t"]]
            eMatGrid_t3 <- initial_package[[3]][["eMatGrid_t"]]
            
            
            block_out <- f(
                alphaHat_t1 = alphaHat_t1,
                alphaHat_t2 = alphaHat_t2,
                alphaHat_t3 = alphaHat_t3,
                betaHat_t1 = betaHat_t1,
                betaHat_t2 = betaHat_t2,
                betaHat_t3 = betaHat_t3,
                c1 = c1,
                c2 = c2,
                c3 = c3,
                eMatGrid_t1 = eMatGrid_t1,
                eMatGrid_t2 = eMatGrid_t2,
                eMatGrid_t3 = eMatGrid_t3,
                H = H,
                H_class = H_class,
                ff = ff,
                eMatRead_t = eMatRead_t,    
                blocked_snps = blocked_snps,
                grid = grid,
                wif0 = wif0,
                s = s,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                priorCurrent_m = priorCurrent_m,
                transMatRate_tc_H = transMatRate_tc_H,
                do_checks = do_checks,
                verbose = verbose,
                fpp_stuff = fpp_stuff,
                shard_check_every_pair = shard_check_every_pair
            )

            ## initial_package = initial_package,
            
            ##
            ## OK, that DID work, neat
            ##
            if (language == "Rcpp") {
                block_out <- append(block_out, list(eMatGrid_t1 = eMatGrid_t1, eMatGrid_t2 = eMatGrid_t2, eMatGrid_t3 = eMatGrid_t3, c1 = c1, c2 = c2, c3 = c3,H = H, alphaHat_t1 = alphaHat_t1, alphaHat_t2 = alphaHat_t2, alphaHat_t3 = alphaHat_t3, betaHat_t1 = betaHat_t1, betaHat_t2 = betaHat_t2, betaHat_t3 = betaHat_t3))
            }
            H <- block_out$H

            ## print(table(ori_H, H))
            ## print(table(H, true_H))
            ## print(round(block_out$shard_block_results, 3))            
            ## print(block_out[["shard_block_results"]][, "ir_chosen"])


            ## want it to match up 1 to 1
            ##

            for(i in 1:3) {
                if (sum(true_H == i) > 0) {
                    expect_equal(length(unique(table(H[true_H == i]))), 1)
                }
            }
            
            final_package <- for_testing_get_full_package_probabilities(block_out[["H"]], fpp_stuff)
            expect_equal(block_out[["eMatGrid_t1"]], final_package[[1]][["eMatGrid_t"]])
            expect_equal(block_out[["eMatGrid_t2"]], final_package[[2]][["eMatGrid_t"]])
            expect_equal(block_out[["c1"]], final_package[[1]][["c"]])
            expect_equal(block_out[["c2"]], final_package[[2]][["c"]])
            expect_equal(block_out[["alphaHat_t1"]], final_package[[1]][["alphaHat_t"]])
            expect_equal(block_out[["alphaHat_t2"]], final_package[[2]][["alphaHat_t"]])
            expect_equal(block_out[["betaHat_t1"]], final_package[[1]][["betaHat_t"]])
            expect_equal(block_out[["betaHat_t2"]], final_package[[2]][["betaHat_t"]])
            if (i_run %in% 1:2) {
                if (language == "R") {    block_out_R <- block_out}
                if (language == "Rcpp") { block_out_Rcpp <- block_out}
            }
            
        }

        ## only tests the second two, the blocks
        expect_equal(block_out_R[["shard_block_results"]], block_out_Rcpp[["shard_block_results"]])
        
        ## print("R")
        ## print(block_out_R[["shard_block_results"]][1:3, ])
        ## print("rcpp")
        ## print(block_out_Rcpp[["shard_block_results"]][1:3, ])

    }
    
})


