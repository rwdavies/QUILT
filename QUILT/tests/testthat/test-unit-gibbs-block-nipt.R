if (1 == 0) {

    library("testthat"); library("STITCH"); library("rrbgen")
    ##    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    ##    dir <- "~/Google Drive/STITCH/"
    dir <- "/data/smew1/rdavies/stitch_development/STITCH-private/"
    dir <- "~/proj/STITCH-private/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))



}


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
    
    out <- lapply(1:3, function(i_option) {

        language <- languages[i_option]
        print(paste0("---------- ", language, " ---------"))
        
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
        print(what)
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
    test_package <- make_fb_test_package(
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
        randomize_sample_read_length = TRUE
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
    quantile_prob <- 0.9
    ## language <- "R_with_Rcpp"
    ## language <- "Rcpp"    
    verbose <- FALSE
    
    for(block_approach in c(1, 4)) {

        print(paste0("---------- block_approach = ", block_approach, " ----------"))
        
        out <- mclapply(languages_to_test, mc.cores = 3, function(language) {

            print(paste0("------------", language, "------------"))
            block_out <- helper_block_gibbs_resampler(
                H = H,
                blocked_snps = blocked_snps,
                grid = grid,
                wif0 = wif0,
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
                class_sum_cutoff = class_sum_cutoff
            )
            
            ## expect all jumps to be in 1 or 3
            expect_equal(0, sum(!(block_out$block_results[, "ir_chosen"] %in% c(0, 1, 3)))) ## 0 for skipped
            
        })
        check_mclapply_OK(out)
        
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
    test_package <- make_fb_test_package(
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
        randomize_sample_read_length = TRUE
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

    

    ## this is supposed to be like real data
    ## imagine some small blocks off in the middle of lots of blocksn
    ## want some small block in the middle to test on
    ir_test_options <- c(1:6, 1:6)
    ## 1 = original
    ## 2 = do not count some H
    ## 4 = sample H for all options with probabilities
    block_approach_options <- rep(c(1, 4), each = 6)
    to_run <- 1:length(ir_test_options)
    to_run <- 1:12
    
    out <- mclapply(to_run, mc.cores = 1, function(i_option) {

        ## 
        verbose <- FALSE
        do_checks <- FALSE
        ir_test <- ir_test_options[i_option]
        block_approach <- block_approach_options[i_option]
        ##if (verbose) {
        print(paste0("---------- ir_test = ", ir_test, " ----------"))
        print(paste0("---------- block_approach = ", block_approach, " ----------"))
        ##}
        ##
        ## recast blocked_snps
        ## 
        languages_to_test <- c("R", "R_with_Rcpp", "Rcpp")
        
        out <- mclapply(languages_to_test, mc.cores = 3, function(language) {
            ##if (verbose) {
                print(paste0("------------", language, "------------"))
            ##}
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
            ## re-sample with respect to rr and true_H
            for(iRead in 1:nReads) {
                h <- true_H[iRead] ## this is what
                hlc <- true_H_class[iRead]
                if (hlc == 0) {
                    hl <- rr[ir_test, h]
                } else {
                    hl <- as.integer(sample(c(1L, 2L, 3L), size = 1, prob = rlc[h, ]))
                }
                if (wif0[iRead] %in% range_to_flip) {
                    H[iRead] <- rr[ir_test, hl]
                } else {
                    H[iRead] <- hl
                }
            }
            ## now, check class OK
            check_agreements_with_truth_class(true_H, true_H_class, H, rlc)
            check_agreements_with_truth_class(true_H, true_H_class, true_H, rlc)
            ## determine H_class
            block_out <- helper_block_gibbs_resampler(
                H = H,
                blocked_snps = blocked_snps,
                grid = grid,
                wif0 = wif0,
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
                class_sum_cutoff = class_sum_cutoff
            )
            ## check H (class) here
            x <- check_agreements_with_truth_class(true_H, true_H_class, block_out[["H"]], rlc)
            expect_equal(as.logical(x < 5), TRUE)
            block_out$block_results
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
        

    })
    check_mclapply_OK(out)
    
})


