if ( 1 == 0 ) {

    library("testthat")
    library("QUILT")
    dir <- "~/proj/QUILT/"
    setwd(paste0(dir, "/QUILT/R"))
    a <- dir(pattern = "*R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)

}


test_that("can bound make_gl_from_u_bq", {

    ## caused by overflow issues
    ## want to bound at least a little bit!

    set.seed(10)
    nSNPs <- 100
    bq <- sample(c(-50, -10, 10, 50), nSNPs, replace = TRUE)
    u <- 1:(nSNPs)

    ## no change
    m1 <- make_gl_from_u_bq(u, bq, nSNPs, minGLValue = 0) ## original

    ## check Rcpp makes sense
    minGLValue <- 1e-2
    gl <- array(NA, dim(m1))
    gl[] <- m1[]
    to_fix <- as.integer(which(colSums(gl < minGLValue) > 0))
    for(i_col in to_fix) {
        a <- gl[1, i_col]
        b <- gl[2, i_col]
        if (a > b) {
            b <- b / a
            a <- 1
            if (b < minGLValue) {
                b <- minGLValue
            }
        } else {
            a <- a / b
            b <- 1
            if (a < minGLValue) {
                a <- minGLValue
            }
        }
        gl[1, i_col] <- a
        gl[2, i_col] <- b
    }

    gl2 <- array(NA, dim(m1))
    gl2[] <- m1[]
    Rcpp_make_gl_bound(gl2, minGLValue, to_fix - 1);

    expect_equal(gl, gl2)

    ## now do bounding
    m2 <- make_gl_from_u_bq(u, bq, nSNPs, minGLValue = 1e-2)

    expect_true(sum(m1 == 1e-2) == 0)    ## confirm that m1 has not changed
    expect_false(sum(m2 == 1e-2) == 0)   ## confirm that m2 has changed

    ## check that rcpp works too

})



test_that(" can select good haps properly", {

    set.seed(91)
    K <- 10000
    Knew <- 30
    K_top_matches <- 5
    Ksubset <- 400
    previously_selected_haplotypes <- sample(1:K, Ksubset - Knew)

    ## now on each hap, do some, add some overlap though, then study...
    new_haps <- list(
        list(90:95, 123:30, 250:260),
        list(310:315, 350:355, 360:370)
    )

    new_ones <- sort(everything_select_good_haps(
        Knew,
        K_top_matches,
        new_haps,
        previously_selected_haplotypes,
        K
    ))

    ## hm, no real test here, just dummy test, should not exceed 500
    expect_equal(max(new_ones) <= 500, TRUE)

})


test_that("can quickly in cpp get positions greater than certain value", {

    ## either an equivalency or not
    K_top_matches <- 10

    for(i in 1:1) {

        if (i == 1) {
            set.seed(100)
            x <- runif(100) ## so will be distinct
        } else {
            x <- seq(0, 1, length.out = 100)
            x[85:95] <- 0.90 ## so this is the one we need more of
        }

        K <- 100
        alphaHat_t <- array(1, c(K, 2))
        iGrid <- 1
        betaHat_t_col <- array(x, c(K, 1))
        gamma_t_col <- array(0, c(K, 1))

        for(j in 1:2) {
            if (j == 1) { f<- Rcpp_get_top_K_or_more_matches_while_building_gamma }
            if (j == 2) { f<- R_get_top_K_or_more_matches_while_building_gamma }
            out <- f(
                alphaHat_t = alphaHat_t,
                betaHat_t_col = betaHat_t_col,
                gamma_t_col = gamma_t_col,
                iGrid = iGrid,
                K = K,
                K_top_matches = K_top_matches
            )
            if (j == 1) {out1 <- out}
            if (j == 2) {out2 <- out}
        }

        expect_equal(out1, out2)

    }

    ## so first pass, can get value, while going over and building
    ## second pass, can get those that meet value

})




test_that("can do simple binary search", {

    vec <- sort(sample(100000, 1000))
    vec2 <- sort(sample(100000, 500))
    vec3 <- sort(sample(100000, 300))
    locations <- c(1, 37, 50, 900, 1000)

    ##
    out <- t(sapply(1:length(locations), function(i) {
        index <- locations[i]
        val <- vec[index]
        c(
            simple_binary_search(val, vec),
            rcpp_simple_binary_search(val, vec) + 1
        )
    }))

    expect_equal(out[, 1], out[, 2])
    expect_equal(out[, 1], match(vec[locations], vec))

    ## length 1
    expect_equal(1, simple_binary_search(7, 7))
    expect_equal(0, rcpp_simple_binary_search(7, 7))

    ## same thing but pad
    for(i_pad in 1:3) {

        if (i_pad == 1) {
            ## top and bottom
            mat <- rbind(cbind(vec2, vec2 + 10), cbind(vec, vec + 10), cbind(vec3, vec3 + 10))
            s1 <- 501; e1 <- 1500
        } else if (i_pad == 2) {
            ## top only
            mat <- rbind(cbind(vec2, vec2 + 10), cbind(vec, vec + 10))
            s1 <- 501; e1 <- 1500
        } else {
            ## bottom only
            mat <- rbind(cbind(vec, vec + 10), cbind(vec3, vec3 + 10))
            s1 <- 1; e1 <- 1000
        }

        ##
        out <- t(sapply(1:length(locations), function(i) {
            index <- locations[i]
            val <- vec[index]
            c(
                simple_binary_matrix_search(val, mat, s1, e1),
                rcpp_simple_binary_matrix_search(val, mat, s1, e1)
            )
        }))

        expect_equal(out[, 1], out[, 2])
        expect_equal(out[, 1], vec[locations] + 10)

    }

})



test_that("can build necessary components from make_rhb_t_equality", {


    set.seed(9920)

    ## make them mostly one of three options, with a few small changes
    K <- 500
    nSNPs <- 100

    reference_haps <- array(0L, c(nSNPs, K))
    ## introduce some shapiness to this
    small <- array(as.integer(runif(nSNPs * 10) > 0.5), c(nSNPs, 10))
    for(k in 1:K) {
        reference_haps[, k] <- small[, sample(10, 1)]
        if ((k %% 10) == 0) {
            small[sample(nSNPs, 1), sample(10, 1)] <- sample(c(0L, 1L), 1)
        }
    }
    ##
    rhi <- reference_haps
    rhi_t <- t(rhi)
    rhb_t <- STITCH::make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)

    ref_error <- 0.01
    ref_one_minus_error <- 1 - ref_error
    nMaxDH <- 6

    for(nMaxDH in c(3, 6, 255, NA)) {

        for(use_hapMatcherR in c(FALSE, TRUE)) {

            ## make haplotype matching objects
            out <- STITCH::make_rhb_t_equality(
                rhb_t = rhb_t,
                nMaxDH = nMaxDH,
                nSNPs = nSNPs,
                ref_error = ref_error,
                use_hapMatcherR = use_hapMatcherR
            )

            distinctHapsB <- out[["distinctHapsB"]]
            distinctHapsIE <- out[["distinctHapsIE"]]
            hapMatcher <- out[["hapMatcher"]]
            hapMatcherR <- out[["hapMatcherR"]]
            eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
            eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]
            eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
            eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]
            nrow_which_hapMatcher_0 <- out[["nrow_which_hapMatcher_0"]]
            ##
            ## perform
            ##
            expect_equal(sum(eMatDH_special_grid_which != 0), length(eMatDH_special_values_list))
            if (length(eMatDH_special_values_list) == 0) {
                expect_equal(nrow_which_hapMatcher_0, 0)
            } else {
                expect_equal(nrow_which_hapMatcher_0, sum(sapply(eMatDH_special_values_list, length)))
            }

            ##
            ## check can re-build
            ##
            for(i in 1:2) {
                if (i == 1) {
                    f <- simple_binary_matrix_search
                } else {
                    f <- rcpp_simple_binary_matrix_search
                }

                rebuilt_rhb_t <- array(as.integer(1), c(nrow(rhb_t), ncol(rhb_t)))
                for(k in 1:nrow(rhb_t)) {
                    for(iGrid in 1:ncol(rhb_t)) {
                        if (!use_hapMatcherR) {
                            i <- hapMatcher[k, iGrid]
                        } else {
                            i <- as.integer(hapMatcherR[k, iGrid])
                        }
                        if (i > 0) {
                            b <- distinctHapsB[i, iGrid]
                        } else {
                            b <- f( ## test both R and Rcpp versions
                                val = k - 1,
                                mat = eMatDH_special_matrix,
                                s1 = eMatDH_special_matrix_helper[iGrid, 1],
                                e1 = eMatDH_special_matrix_helper[iGrid, 2]
                            )
                        }
                        rebuilt_rhb_t[k, iGrid] <- b
                    }
                }
            }

            expect_equal(rhb_t, rebuilt_rhb_t)

        }

    }

})


test_that("can fast eMat", {


    ## simulate some stuff
    K <- 100
    nSNPs <- 200
    nReads <- 20
    out <- make_fb_test_package(K = K, nReads = nReads, nSNPs = nSNPs)

    rhi <- t(round(out$eHapsCurrent_tc[, , 1]))
    rhi_t <- t(rhi)
    rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)

    ref_error <- 0.01
    sampleReads <- out$sampleReads
    nSNPs <- out$nSNPs
    eHapsCurrent_tc <- array(0, c(K, nSNPs, 1))
    eHapsCurrent_tc[] <- rhi_t
    eHapsCurrent_tc[eHapsCurrent_tc == 0] <- ref_error
    eHapsCurrent_tc[eHapsCurrent_tc == 1] <- 1 - ref_error

    maxDifferenceBetweenReads <- 1e10
    normal_eMatRead_t <- array(1, c(K, nReads))
    rcpp_make_eMatRead_t(
        eMatRead_t = normal_eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = eHapsCurrent_tc,
        s = 0,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads ,
        Jmax = 1000,
        eMatHapOri_t = array(0, c(1, 1)),
        pRgivenH1 = array(0, 1),
        pRgivenH2 = array(0, 1),
        prev = 0,
        suppressOutput = 1,
        prev_section = "wer",
        next_section = "wer",
        run_pseudo_haploid = FALSE,
        rescale_eMatRead_t = FALSE
    )

    eMatRead_t <- make_eMatRead_t_using_binary(
        sampleReads = sampleReads,
        rhb_t = rhb_t,
        nSNPs = nSNPs,
        ref_error = ref_error,
        language = "R",
        n = 10
    )
    expect_equal(eMatRead_t, normal_eMatRead_t)

    ## OK!


    Rcpp_eMatRead_t <- make_eMatRead_t_using_binary(
        sampleReads = sampleReads,
        rhb_t = rhb_t,
        nSNPs = nSNPs,
        ref_error = ref_error,
        language = "Rcpp"
    )

    expect_equal(eMatRead_t, Rcpp_eMatRead_t)

})




test_that("can avoid normalizing alphaHat and betaHat throughout forward algorithm, but recover correct gamma when asked, as well as dosages", {

    ##
    ## setup
    ##
    ##
    speed_test <- FALSE
    rebuild <- FALSE
    use_hapMatcherR <- FALSE
    if (!speed_test) {
        nSNPs <- 500
        K <- 100
        nMaxDH <- 20
        languages <- c("R", "Rcpp", "Rcpp2")
        suppressOutput <- 1
    } else {
        nSNPs <- 20000
        K <- 30000
        nMaxDH <- 100
        languages <- c("Rcpp", "Rcpp2")
        suppressOutput <- 0
    }
    file <- paste0("~/Downloads/test_package.", speed_test, ".RData")
    if (!speed_test | !file.exists(file) | rebuild) {
        message("building")
        test_package <- make_reference_single_test_package(
            nSNPs = nSNPs,
            K = K,
            nMaxDH = nMaxDH
        )
        if (speed_test) {
            save(test_package, file = file, compress = FALSE)
        }
    } else {
        print("loading")
        load(file)
    }
    distinctHapsB <- test_package[["distinctHapsB"]]
    distinctHapsIE <- test_package[["distinctHapsIE"]]
    hapMatcher <- test_package[["hapMatcher"]]
    hapMatcherR <- test_package[["hapMatcherR"]]
    rhb_t <- test_package[["rhb_t"]]
    gl <- test_package[["gl"]]
    transMatRate_t <- test_package[["transMatRate_t"]]
    ref_error <- test_package[["ref_error"]]
    eMatDH_special_grid_which <- test_package[["eMatDH_special_grid_which"]]
    eMatDH_special_values_list <- test_package[["eMatDH_special_values_list"]]
    eMatDH_special_matrix <- test_package[["eMatDH_special_matrix"]]
    eMatDH_special_matrix_helper <- test_package[["eMatDH_special_matrix_helper"]]
    gammaSmall_cols_to_get <- test_package[["gammaSmall_cols_to_get"]]
    K <- nrow(rhb_t)
    nSNPs <- ncol(gl)
    nGrids <- ncol(rhb_t)
    nSmallGammaGrids <- sum(gammaSmall_cols_to_get != (-1))
    truth_hap <- test_package[["truth_hap"]]
    get_best_haps_from_thinned_sites <- TRUE
    K_top_matches <- 5
    use_eMatDH <- TRUE
    always_normalize <- TRUE
    language <- "Rcpp"
    min_emission_prob_normalization_threshold <- 1e-100 ## make much smaller than usual
    use_eMatDH_special_symbols <- TRUE
    is_version_2 <- TRUE
    normalize_emissions <- FALSE

    ##
    ## run normal version in R, or Rcpp, and make sure the same
    ##
    master_f <- function(
                         always_normalize,
                         language,
                         is_version_2,
                         normalize_emissions  = FALSE,
                         use_eMatDH_special_symbols = FALSE,
                         use_hapMatcherR = FALSE,
                         use_eigen = FALSE
                         ) {

        if (language == "R") {
            f <- R_haploid_dosage_versus_refs
        } else if (language == "Rcpp") {
            f <- Rcpp_haploid_dosage_versus_refs
        } else {
            f <- Rcpp_haploid_dosage_versus_refs
        }
        if (speed_test) {
            print(
                paste0(
                    "-------",
                    "language = ", language, ", ",
                    "always_normalize = ", always_normalize, " ",
                    "is_version_2 = ", is_version_2, " ",
                    "normalize_emissions = ", normalize_emissions, " ",
                    "use_eigen = ", use_eigen
                )
            )
        }
        ##

        stopifnot(!(use_eigen && is_version_2))
        if (use_eigen) {
            arma_alphaHat_t <- array(0, c(1, 1))
            eigen_alphaHat_t <- array(0, c(K, nGrids))
        } else  {
            arma_alphaHat_t <- array(0, c(K, nGrids))
            eigen_alphaHat_t <- array(0, c(1, 1))
        } 


        betaHat_t <- array(0, c(K, nGrids))
        c <-  array(1, c(nGrids))
        gamma_t <- array(0, c(K, nGrids))
        gammaSmall_t <- array(0, c(K, nSmallGammaGrids))
        dosage <- numeric(nSNPs)
        best_haps_stuff_list <- as.list(1:sum(gammaSmall_cols_to_get >= 0))

        if (use_eMatDH_special_symbols) {
            rhb_t <- matrix(as.integer(1), 1, 1) ## nuke!
        }

        ##
        out <- f(
            gl = gl,
            arma_alphaHat_t = arma_alphaHat_t,
            eigen_alphaHat_t = eigen_alphaHat_t,
            betaHat_t = betaHat_t,
            c = c,
            gamma_t = gamma_t,
            gammaSmall_t = gammaSmall_t,
            dosage = dosage,
            transMatRate_t = transMatRate_t,
            rhb_t = rhb_t,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            eMatDH_special_matrix = eMatDH_special_matrix,
            use_eMatDH_special_symbols = use_eMatDH_special_symbols,
            ref_error = ref_error,
            use_eMatDH = use_eMatDH,
            distinctHapsB = distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            use_hapMatcherR = use_hapMatcherR,
            return_extra = TRUE,
            gammaSmall_cols_to_get = gammaSmall_cols_to_get,
            get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
            return_dosage = TRUE,
            return_gammaSmall_t = TRUE,
            return_gamma_t = TRUE,
            K_top_matches = K_top_matches,
            always_normalize = always_normalize,
            min_emission_prob_normalization_threshold = min_emission_prob_normalization_threshold,
            best_haps_stuff_list = best_haps_stuff_list,
            eMatDH_special_values_list = eMatDH_special_values_list,
            eMatDH_special_grid_which = eMatDH_special_grid_which,
            suppressOutput = suppressOutput,
            normalize_emissions = normalize_emissions,
            use_eigen = use_eigen,
            is_version_2 = is_version_2
        )
        ## 
        if (language == "R") {
            alphaHat_t <- out[["alphaHat_t"]]
            c <- out[["c"]]
            gammaSmall_t <- out[["gammaSmall_t"]]
            dosage <- out[["dosage"]]
            betaHat_t <- out[["betaHat_t"]]
        }

        to_return <- list(
            betaHat_t = betaHat_t,
            c = c,
            gammaSmall_t = gammaSmall_t,
            dosage = dosage
        )

        if (use_eigen) {
            to_return <- append(to_return, list(alphaHat_t = eigen_alphaHat_t))
        } else {
            to_return <- append(to_return, list(alphaHat_t = arma_alphaHat_t))
        }
        to_return
    }

    ## note - is_version_2 rendered obsolete here, so will always use the single available version!
    out_R_always_normalize <- master_f(TRUE, "R", TRUE)
    out_Rcpp_always_normalize <- master_f(TRUE, "Rcpp", FALSE)
    out_Rcpp2_always_normalize <- master_f(TRUE, "Rcpp", TRUE)
    out_Rcpp3_always_normalize <- master_f(TRUE, "Rcpp", TRUE, use_eMatDH_special_symbols = TRUE)
    out_Rcpp4_always_normalize <- master_f(TRUE, "Rcpp", TRUE, use_eMatDH_special_symbols = TRUE, use_hapMatcherR = TRUE)
    out_Rcpp5_always_normalize <- master_f(TRUE, "Rcpp", is_version_2 = FALSE, use_eMatDH_special_symbols = TRUE, use_hapMatcherR = TRUE, use_eigen = TRUE)
    
    out_R_seldom_normalize <- master_f(FALSE, "R", TRUE)
    out_Rcpp_seldom_normalize <- master_f(FALSE, "Rcpp", FALSE)
    out_Rcpp2_seldom_normalize <- master_f(FALSE, "Rcpp", TRUE)
    out_Rcpp3_seldom_normalize <- master_f(FALSE, "Rcpp", TRUE, use_eMatDH_special_symbols = TRUE)
    out_Rcpp4_seldom_normalize <- master_f(FALSE, "Rcpp", TRUE, use_eMatDH_special_symbols = TRUE, use_hapMatcherR = TRUE)
    out_Rcpp5_seldom_normalize <- master_f(FALSE, "Rcpp", is_version_2 = FALSE, use_eMatDH_special_symbols = TRUE, use_hapMatcherR = TRUE, use_eigen = TRUE)    


    ##
    ## do checks if not just for speed
    ##
    if (!speed_test) {
        ##
        ## check results for R versions
        ##
        expect_equal(out_R_always_normalize[["dosage"]], out_R_seldom_normalize[["dosage"]])
        expect_equal(out_R_always_normalize[["gammaSmall_t"]], out_R_seldom_normalize[["gammaSmall_t"]])
        expect_equal(sum(log(out_R_always_normalize[["c"]])), sum(log(out_R_seldom_normalize[["c"]])))
        expect_true(sum(abs(out_R_always_normalize[["c"]] - out_R_seldom_normalize[["c"]]) > 0) > 0)
    }
    if (!speed_test) {
        ##
        ## check cpp version here too
        ##
        ##
        expect_equal(sum(log(out_Rcpp_always_normalize[["c"]])), sum(log(out_Rcpp_seldom_normalize[["c"]])))
        expect_true(sum(abs(out_Rcpp_always_normalize[["c"]] - out_Rcpp_seldom_normalize[["c"]]) > 0) > 0)
        expect_equal(sum(log(out_R_always_normalize[["c"]])), sum(log(out_Rcpp_always_normalize[["c"]])))
        expect_equal(sum(log(out_R_always_normalize[["c"]])), sum(log(out_Rcpp_seldom_normalize[["c"]])))
        ## check dosages vs each other and R
        expect_equal(out_R_always_normalize[["dosage"]], out_Rcpp_always_normalize[["dosage"]], tolerance = 1e-5)
        expect_equal(out_R_always_normalize[["dosage"]], out_Rcpp_seldom_normalize[["dosage"]], tolerance = 1e-5)
        ## check gammaSmall vs each other and R
        expect_equal(out_R_always_normalize[["gammaSmall_t"]], out_Rcpp_always_normalize[["gammaSmall_t"]], tolerance = 1e-5)
        expect_equal(out_R_always_normalize[["gammaSmall_t"]], out_Rcpp_seldom_normalize[["gammaSmall_t"]], tolerance = 1e-5)
    }


    ##
    ## check alternate cpp version
    ##
    if (!speed_test) {
        ##
        ## c
        ##
        expect_equal(sum(log(out_Rcpp2_always_normalize[["c"]])), sum(log(out_Rcpp2_seldom_normalize[["c"]])))
        expect_true(sum(abs(out_Rcpp2_always_normalize[["c"]] - out_Rcpp2_seldom_normalize[["c"]]) > 0) > 0)
        expect_equal(sum(log(out_R_always_normalize[["c"]])), sum(log(out_Rcpp2_always_normalize[["c"]])))
        expect_equal(sum(log(out_R_always_normalize[["c"]])), sum(log(out_Rcpp2_seldom_normalize[["c"]])))
        expect_equal(sum(log(out_R_always_normalize[["c"]])), sum(log(out_Rcpp3_seldom_normalize[["c"]])))
        expect_equal(sum(log(out_R_always_normalize[["c"]])), sum(log(out_Rcpp4_seldom_normalize[["c"]])))
        expect_equal(sum(log(out_R_always_normalize[["c"]])), sum(log(out_Rcpp5_seldom_normalize[["c"]])))        
        ##
        ## dosages
        ##
        expect_equal(out_R_always_normalize[["dosage"]], out_Rcpp2_always_normalize[["dosage"]], tolerance = 1e-5)
        expect_equal(out_R_always_normalize[["dosage"]], out_Rcpp2_seldom_normalize[["dosage"]], tolerance = 1e-5)
        expect_equal(out_Rcpp_always_normalize[["dosage"]], out_Rcpp2_always_normalize[["dosage"]], tolerance = 1e-5)
        expect_equal(out_Rcpp_always_normalize[["dosage"]], out_Rcpp3_always_normalize[["dosage"]], tolerance = 1e-5)
        expect_equal(out_Rcpp_always_normalize[["dosage"]], out_Rcpp4_always_normalize[["dosage"]], tolerance = 1e-5)
        expect_equal(out_Rcpp_always_normalize[["dosage"]], out_Rcpp5_always_normalize[["dosage"]], tolerance = 1e-5)        
        ##
        ## gammas
        ##
        expect_equal(out_R_always_normalize[["gammaSmall_t"]], out_Rcpp2_always_normalize[["gammaSmall_t"]])
        expect_equal(out_R_always_normalize[["gammaSmall_t"]], out_Rcpp2_seldom_normalize[["gammaSmall_t"]])
        expect_equal(out_Rcpp_always_normalize[["gammaSmall_t"]], out_Rcpp2_always_normalize[["gammaSmall_t"]])
        expect_equal(out_Rcpp_always_normalize[["gammaSmall_t"]], out_Rcpp3_always_normalize[["gammaSmall_t"]])
        expect_equal(out_Rcpp_always_normalize[["gammaSmall_t"]], out_Rcpp4_always_normalize[["gammaSmall_t"]])
        expect_equal(out_Rcpp_always_normalize[["gammaSmall_t"]], out_Rcpp5_always_normalize[["gammaSmall_t"]])        
    }


    ##
    ## check dosages and gamma for 4 Rcpp options, without normalization
    ##
    if (!speed_test) {

        out_Rcpp_baseline <- master_f(always_normalize = TRUE, language = "Rcpp", is_version_2 = FALSE, normalize_emissions = FALSE)
        for(always_normalize in c(FALSE, TRUE)) {
            for(normalize_emissions in c(FALSE, TRUE)) {
                for(i_version_eigen in 1:3) {
                    is_version_2 <- c(FALSE, TRUE, FALSE)[i_version_eigen]
                    use_eigen <- c(FALSE, FALSE, TRUE)[i_version_eigen]
                    ## print(
                    ##     paste0(
                    ##         "-------",
                    ##         "language = ", language, ", ",
                    ##         "always_normalize = ", always_normalize, " ",
                    ##         "is_version_2 = ", is_version_2, " ",
                    ##         "normalize_emissions = ", normalize_emissions
                    ##     )
                    ## )
                    out_Rcpp_comparison <- master_f(always_normalize = always_normalize, language = "Rcpp", is_version_2 = is_version_2, normalize_emissions = normalize_emissions, use_eigen = use_eigen)
                    expect_equal(out_Rcpp_baseline[["gammaSmall_t"]], out_Rcpp_comparison[["gammaSmall_t"]], tolerance = 1e-3) ## yuck, scary? but checked on larger data?
                    expect_equal(out_Rcpp_baseline[["dosage"]],       out_Rcpp_comparison[["dosage"]], tolerance = 1e-3)
                }
            }
        }
    }


    ##
    ## extra stuff for later (possibly)
    ##
    if (1 == 0) {

        x <- colSums(out_Rcpp2_always_normalize_c_inexact[["gammaSmall_t"]])
        expect_equal(x, rep(1, length(x)))

        expect_equal(colSums(out_Rcpp2_always_normalize_c_inexact[["gammaSmall_t"]]), 1)
        alphaHat_t1 <- out_Rcpp2_always_normalize_c_exact[["alphaHat_t"]]
        c1 <- out_Rcpp2_always_normalize_c_exact[["c"]]
        alphaHat_t2 <- out_Rcpp2_always_normalize_c_inexact[["alphaHat_t"]]
        c2 <- out_Rcpp2_always_normalize_c_inexact[["c"]]
        for(iGrid in which(gammaSmall_cols_to_get > (-1))) {
            x <- unique(round(log(alphaHat_t1[, iGrid] / alphaHat_t2[, iGrid]), 5))
            expect_equal(length(x), 1)
            ## expect_equal(x, exp(sum(log(c1[1:iGrid])) - sum(log(c2[1:iGrid]))))
            expect_equal(
                exp(x),
                exp(sum(log(c1[1:iGrid])) - sum(log(c2[1:iGrid])))
            )
        }

        expect_equal(out_Rcpp2_always_normalize_c_exact[["gammaSmall_t"]], out_Rcpp2_always_normalize_c_inexact[["gammaSmall_t"]])

        expect_equal(out_R_always_normalize[["dosage"]],    out_Rcpp2_always_normalize_c_exact[["dosage"]])
        expect_equal(out_R_always_normalize[["dosage"]],    out_Rcpp2_always_normalize_c_inexact[["dosage"]])
        expect_equal(out_R_always_normalize[["dosage"]],    out_Rcpp2_seldom_normalize_c_exact[["dosage"]])
        expect_equal(out_R_always_normalize[["dosage"]],    out_Rcpp2_seldom_normalize_c_inexact[["dosage"]])

    }

})




test_that("can run a single gl sample through reference haplotypes quickly with grid of 32", {

    ## switch on and off. off does thes tests, on checks cpp times
    speed_test <- FALSE
    if (!speed_test) {
        nSNPs <- 100
        L <- sort(sample(1:1000, 100))
        suppressOutput <- 1
    } else {
        suppressOutput <- 0
        nSNPs <- 1000
        L <- sort(sample(1:10000, 1000))
        expect_equal(1, 1)
    }

    ## looks good
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

    ##
    dl <- diff(L_grid)
    expRate <- 10 ## increase rate - easier to see jump issues
    nGen <- 10
    sigmaCurrent <- exp(-nGen * expRate / 100 / 1000000 * dl)
    transMatRate_t <- rbind(sigmaCurrent, 1 - sigmaCurrent)


    use_eMatDH <- TRUE; i_setup <- 2

    gammaSmall_cols_to_get <- c(-1, 0, -1, 1) ## i.e. 0-based, what to put it in
    nSmallGammaGrids <- 2
    K_top_matches <- 5
    get_best_haps_from_thinned_sites <- TRUE

    use_eMatDH <- TRUE
    i_setup <- 1
    always_normalize <- FALSE

    ## build some haplotypes and encode them
    for(i_setup in 1:3) {

        for(always_normalize in c(FALSE, TRUE)) {

            setup_result <- list()

            for(use_eMatDH in c(FALSE, TRUE)) {

                ## for cpp
                return_dosage <- TRUE
                return_betaHat_t <- TRUE
                return_gamma_t <- TRUE
                return_gammaSmall_t <- TRUE
                get_best_haps_from_thinned_sites <- TRUE

                ##
                ## i_setup = 1 ---
                ## i_setup = 2 ---
                ## i_setup = 3 --- same as above, but only get best_haps from thinned sites
                if (!suppressOutput) {
                    print(paste0("use_eMatDH = ", use_eMatDH, ", i_setup = ", i_setup))
                }
                if (i_setup == 1) {
                    ## test gamma makes sense
                    set.seed(9910)
                    K <- 100
                    reference_haps <- array(as.integer(runif(nSNPs * K) > 0.5), c(nSNPs, K))
                } else if (i_setup == 2  | i_setup == 3){
                    ## test can easily work with large, mostly similar emissions
                    set.seed(242)
                    K <- 10000
                    reference_haps <- array(NA, c(nSNPs, K))
                    haps <- array(sample(c(0, 1), nSNPs * 3, replace = TRUE), c(3, nSNPs))
                    for(k in 1:K) {
                        reference_haps[, k] <- haps[sample(1:3, 1), ]
                    }
                    ## tiny noise
                    w <- runif(K * nSNPs) < 0.01
                    reference_haps[w] <- 1 - reference_haps[w]
                    if (i_setup == 3) {
                        ## particularly important for cpp
                        return_dosage <- FALSE
                        return_betaHat_t <- FALSE
                        return_gamma_t <- FALSE
                        return_gammaSmall_t <- FALSE
                        get_best_haps_from_thinned_sites <- TRUE
                    }
                }

                ##
                best_haps_stuff_list <- list()
                best_haps_stuff_list <- as.list(1:sum(gammaSmall_cols_to_get >= 0))

                ## make them mostly one of three options, with a few small changes
                rhi <- reference_haps
                rhi_t <- t(rhi)
                rhb_t <- STITCH::make_rhb_t_from_rhi_t(rhi_t)
                rhb <- t(rhb_t)

                ## force use of special
                nMaxDH  <- 100
                ref_error <- 0.01
                ref_one_minus_error <- 1 - ref_error


                ## make haplotype matching objects
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
                if (!use_eMatDH) {
                    use_eMatDH_special_symbols <- FALSE
                } else {
                    use_eMatDH_special_symbols <- TRUE
                }
                use_hapMatcherR <- FALSE ## possibly change this!
                out <- STITCH::make_rhb_t_equality(
                    rhb_t = rhb_t,
                    nMaxDH = nMaxDH,
                    nSNPs = nSNPs,
                    ref_error = ref_error,
                    use_hapMatcherR = TRUE
                )
                hapMatcherR <- out[["hapMatcherR"]]


                ## make some gls from this
                my_hap <- c(
                    rhi_t[1, 1:36],
                    rhi_t[2, 37:60],
                    rhi_t[3, 61:100]
                )
                ## my_hap <- c(
                ##     rhi_t[1, 1:32],
                ##     rhi_t[2, 33:64],
                ##     rhi_t[3, 65:100]
                ## )
                ## draw some reads from this
                ## draw some reads from this
                nReads <- round(nSNPs * 1.5)
                u <- c(1, sort(sample(1:nSNPs, nReads, replace = TRUE)), nSNPs) ## include start and end
                bq <- rep(-10, length(u))
                bq[my_hap[u] == 1] <- 10
                gl <- make_gl_from_u_bq(u, bq, nSNPs)
                ## add a bit of noise
                ## w <- sample(1:length(bq), round(length(bq) / 10))
                ## bq[w] <- - bq[w]
                ## sample(c(-10, 10), nReads, replace = TRUE)
                ## round(t(rbind(gl[, 1:32], my_hap[1:32])), 3)
                ## from this, make prob of ref and prob of alt emissio for each
                ## my_hap <- 0.90 * my_hap + 0.10 * runif(nSNPs)
                rm(reference_haps, rhi_t, rhi); gc(reset = TRUE); gc(reset = TRUE)

                if (!speed_test) {
                    ## check the same
                    R_eMatDH <- build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error)
                    Rcpp_eMatDH <- Rcpp_build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error)
                    expect_equal(R_eMatDH, Rcpp_eMatDH)
                }

                K <- nrow(rhb_t)
                nGrids <- ncol(rhb_t)
                ##
                alphaHat_t <- array(0, c(K, nGrids))
                betaHat_t <- array(0, c(K, nGrids))
                c <-  array(1, c(nGrids))
                gamma_t <- array(0, c(K, nGrids))
                gammaSmall_t <- array(0, c(K, nSmallGammaGrids))
                dosage <- numeric(nSNPs)

                ## so now, want to
                if (!speed_test) {

                    ##
                    a <- system.time(
                        outR <- R_haploid_dosage_versus_refs(
                            gl = gl,
                            arma_alphaHat_t = alphaHat_t,
                            betaHat_t = betaHat_t,
                            c = c,
                            gamma_t = gamma_t,
                            gammaSmall_t = gammaSmall_t,
                            dosage = dosage,
                            transMatRate_t = transMatRate_t,
                            rhb_t = rhb_t,
                            ref_error = ref_error,
                            use_eMatDH = use_eMatDH,
                            distinctHapsB = distinctHapsB,
                            distinctHapsIE = distinctHapsIE,
                            hapMatcher = hapMatcher,
                            return_extra = TRUE,
                            gammaSmall_cols_to_get = gammaSmall_cols_to_get,
                            get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
                            return_dosage = TRUE,
                            return_gammaSmall_t = TRUE,
                            return_gamma_t = TRUE,
                            K_top_matches = K_top_matches
                        )
                    )

                    if (!suppressOutput) {
                        print("R speed")
                    }
                    R_gamma_t <- outR[["gamma_t"]]
                    R_gammaSmall_t <- outR[["gammaSmall_t"]]
                    R_dosage <- outR[["dosage"]]
                    R_alphaHat_t <- outR[["alphaHat_t"]]
                    R_betaHat_t <- outR[["betaHat_t"]]
                    R_best_haps_stuff_list <- outR[["best_haps_stuff_list"]]
                    R_c <- outR[["c"]]
                    expect_equal(colSums(R_gamma_t), rep(1, 4))
                    expect_true((cor(R_dosage, my_hap) ** 2) > 0.7) ## is better for truer match against reference
                    ## alternate

                }

                a <- system.time(
                    Rcpp_haploid_dosage_versus_refs(
                        gl = gl,
                        arma_alphaHat_t = alphaHat_t,
                        eigen_alphaHat_t = array(0, c(1, 1)),
                        betaHat_t = betaHat_t,
                        c = c,
                        gamma_t = gamma_t,
                        dosage = dosage,
                        transMatRate_t = transMatRate_t,
                        rhb_t = rhb_t,
                        ref_error = ref_error,
                        use_eMatDH = use_eMatDH,
                        distinctHapsB = distinctHapsB,
                        distinctHapsIE = distinctHapsIE,
                        eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                        eMatDH_special_matrix = eMatDH_special_matrix,
                        use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                        eMatDH_special_grid_which = eMatDH_special_grid_which,
                        eMatDH_special_values_list = eMatDH_special_values_list,
                        hapMatcher = hapMatcher,
                        hapMatcherR = hapMatcherR,
                        use_hapMatcherR = use_hapMatcherR,
                        gammaSmall_t = gammaSmall_t,
                        gammaSmall_cols_to_get = gammaSmall_cols_to_get,
                        suppressOutput = suppressOutput,
                        K_top_matches = K_top_matches,
                        best_haps_stuff_list = best_haps_stuff_list,
                        get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
                        return_betaHat_t = return_betaHat_t,
                        return_dosage = return_dosage,
                        return_gamma_t = return_gamma_t,
                        return_gammaSmall_t = return_gammaSmall_t,
                        normalize_emissions = FALSE,
                        always_normalize = always_normalize,
                        use_eigen = FALSE
                    )
                )
                if (!suppressOutput) {
                    print("cpp speed")
                    print(a)
                }

                if (!speed_test) {
                    if (!suppressOutput) {
                        print(use_eMatDH)
                    }
                    if (i_setup <= 2) {
                        if (always_normalize) {
                            expect_equal(R_alphaHat_t, alphaHat_t)
                            expect_equal(R_c, c)
                            expect_equal(R_betaHat_t, betaHat_t)
                        }
                        expect_equal(R_gamma_t, gamma_t)
                        expect_equal(R_dosage, dosage)
                    }
                    expect_equal(R_best_haps_stuff_list, best_haps_stuff_list)
                }

                setup_result[[as.character(use_eMatDH)]] <- list(
                    alphaHat_t = alphaHat_t,
                    betaHat_t = betaHat_t,
                    gamma_t = gamma_t,
                    dosage = dosage,
                    c = c
                )
            }

            ## compare results for a setup per eMatDH
            for(what in c("alphaHat_t", "betaHat_t", "gamma_t", "dosage", "c")) {
                m1 <- setup_result[["FALSE"]][[what]]
                m2 <- setup_result[["TRUE"]][[what]]
                ## don't check alpha and beta if always_normize is on
                if (
                    !(!always_normalize && (what == "alphaHat_t" | what == "betaHat_t" | what == "c")) &&
                    !(what == "alphaHat_t" & !return_gamma_t)
                ) {
                    if (suppressWarnings((max(m1 / m2 - 1, na.rm = TRUE) > 1e-8))) {
                        print(i_setup)
                        print(what)
                    }
                    expect_equal(setup_result[["FALSE"]][[what]], setup_result[["TRUE"]][[what]])
                }
            }

        }

    }

})




test_that("profile", {

    skip("for profiling")

    if (1 == 0) {

        load("~/Downloads/impute_develop.RData")
        nMaxDH <- 2 ** 8 - 1
        ref_error <- 1e-2
        out <- STITCH::make_rhb_t_equality(
            rhb_t = rhb_t,
            nMaxDH = nMaxDH,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
        distinctHapsB <- out[["distinctHapsB"]]
        distinctHapsIE <- out[["distinctHapsIE"]]
        hapMatcher <- out[["hapMatcher"]]
        eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
        eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]
        save(
            ancAlleleFreqAll, cM_grid, distinctHapsB, distinctHapsIE, full_gammaSmall_cols_to_get, full_transMatRate_t_H, H, hapMatcher, heuristic_match_thin, inRegion2, K_top_matches, Knew, L, L_grid, nSNPs, previously_selected_haplotypes, ref_error, rhb_t, sampleReads,
            eMatDH_special_values_list, eMatDH_special_grid_which,
            nMaxDH,
            file = "~/Downloads/impute_develop.RData", compress = FALSE
        )

    }

    ##
    ## load("~/Downloads/temp_quilt_speed_test.RData")
    ##
    load("~/Downloads/impute_develop.RData")
    transMatRate_t <- full_transMatRate_t_H
    K <- nrow(rhb_t)
    nGrids <- ncol(rhb_t)
    i_hap <- 1
    u <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
    bq <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
    w <- bq != 0
    bq <- bq[w]
    u <- u[w]
    if (length(u) == 0) {
        print_message(paste0("Read label assignment includes no reads for haplotype ", i_hap))
    }
    ##
    gl <- make_gl_from_u_bq(u, bq, nSNPs)
    ## so u is 1-based, so e.g. having 1-32 means grid 1, 33-64 means grid 2, etc
    ## so unique(floor(u / 32)) having 2 means yes there is a var in grid 2, etc
    grid_has_variant <- array(FALSE, nGrids)
    grid_has_variant[unique(floor(u / 32))] <- TRUE ## at least for 32-base grids


    ##
    ## new version
    ##
    use_eMatDH <- TRUE
    ref_one_minus_error <- 1 - ref_error
    eMatDH <- Rcpp_build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error)
    betaHat_t <- array(0, c(1, 1))
    gamma_t <- array(0, c(1, 1))
    gammaSmall_t <- array(0, c(1, 1))
    ww <- seq(1, nGrids, length.out = max(1, round(heuristic_match_thin * nGrids)))
    full_gammaSmall_cols_to_get <- array(-1, nGrids)
    full_gammaSmall_cols_to_get[ww] <- 0:(length(ww) - 1)
    return_gammaSmall_t <- FALSE
    ##
    only_store_alpha_at_gamma_small <- TRUE
    ## try to store where the hapMatcher0 are
    add_zero_row <- TRUE
    eMatDH_bigger <- Rcpp_build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error, add_zero_row)

    ##
    ## get first alpha and c
    ##
    iGrid <- 0
    s <- 32 * iGrid + 1 ## 1-based start
    e <- min(32 * (iGrid + 1), nSNPs) ## 1-based end
    nSNPsLocal <- e - s + 1
    gl_local <- gl[, s : e, drop = FALSE]
    first_alpha <- array(0, K)
    one_over_K <- 1 / K
    for(k in 0:(K - 1)) {
        if (use_eMatDH) {
            dh <- hapMatcher[k + 1, iGrid + 1]
        } else {
            dh <- 0
        }
        if (dh > 0) {
            prob <- eMatDH[dh, iGrid + 1]
        } else {
            prob <- get_prob_for_k(rhb_t, k + 1, iGrid + 1, nSNPsLocal, ref_error, ref_one_minus_error, gl_local)
        }
        first_alpha[k + 1] <- prob * one_over_K
    }
    first_c <- 1 / sum(first_alpha)
    first_alpha <- first_c * first_alpha

    ##
    ## profile forward
    ##

    if (1 == 1) {

        print("-----------------------")
        print("--------forward--------")
        print("-----------------------")

        alphaHat_t_new <- array(0, c(K, nGrids))
        alphaHat_t_new[, 1] <- first_alpha
        c1 <- array(1, c(1, nGrids))
        alphaHat_t_original <- array(0, c(K, nGrids))
        alphaHat_t_original[, 1] <- first_alpha
        c2 <- array(1, c(1, nGrids))
        c1[1] <- first_c
        c2[1] <- first_c

        ## OK, is there some way I can avoid more multiplications!

        ## run them normally too for checking!
        print(paste0("----------new version----------------"))
        nReps <- 10
        print(system.time({
            for(i in 1:nReps) {
                Rcpp_haploid_reference_single_forward_version2(
                    full_gammaSmall_cols_to_get,
                    gl,
                    alphaHat_t_new,
                    c1,
                    transMatRate_t,
                    rhb_t,
                    hapMatcher,
                    eMatDH_bigger,
                    nGrids,
                    nSNPs,
                    K,
                    use_eMatDH,
                    ref_error,
                    only_store_alpha_at_gamma_small,
                    always_normalize = FALSE,
                    min_emission_prob_normalization_threshold = 1e-100,
                    eMatDH_special_grid_which,
                    eMatDH_special_values_list
                )
            }}))

        print(paste0("----------original version----------------"))
        print(system.time({
            for(i in 1:nReps) {
                Rcpp_haploid_reference_single_forward(
                    full_gammaSmall_cols_to_get,
                    gl,
                    alphaHat_t_original,
                    c2,
                    transMatRate_t,
                    rhb_t,
                    hapMatcher,
                    eMatDH_bigger,
                    nGrids,
                    nSNPs,
                    K,
                    use_eMatDH,
                    ref_error,
                    only_store_alpha_at_gamma_small,
                    always_normalize = FALSE,
                    min_emission_prob_normalization_threshold = 1e-100
                )
            }}))


        print(paste0("--------done forward testing version----------"))
        if (1 == 0) {
            print(c1[1:10])
            print(c2[1:10])
            print(c1[1:10] - c2[1:10])
            print("aaaaaaaaaaa")
            print(alphaHat_t_new[1:10, 2])
            print(alphaHat_t_original[1:10, 2])
            print(alphaHat_t_new[1:10, 2] - alphaHat_t_original[1:10, 2])
            print("bbbbbbbbbbbb")
            print(grid_has_variant[1:10])
            print(sum(log10(c1[-1])))
            print(sum(log10(c2[-1])))
            print(prod((c1[-1])))
            print(prod((c2[-1])))
            print(median(alphaHat_t_new[, nGrids]))
            print(median(alphaHat_t_original[, nGrids]))
        }

        ## is this true
        expect_equal(
            sum(log(c1)),
            sum(log(c2))
        )

        ## check relationship between alphas is the same, at the end, when normalized
        expect_equal(
            alphaHat_t_new[, nGrids],
            alphaHat_t_original[, nGrids]
        )

        ## also, want alphas to be the same, throughout!
        for(iGrid in which(full_gammaSmall_cols_to_get > (-1))) {
            x <- unique(round(log(alphaHat_t_new[, iGrid] / alphaHat_t_original[, iGrid]), 5))
            expect_equal(length(x), 1)
            ## expect_equal(x, exp(sum(log(c1[1:iGrid])) - sum(log(c2[1:iGrid]))))
            expect_equal(
                exp(x),
                exp(sum(log(c1[1:iGrid])) - sum(log(c2[1:iGrid]))),
                tol = 1e-5
            )
        }

        ## stop("WER")

    }


    if (1 == 1) {

        print("-----------------------")
        print("--------backward-------")
        print("-----------------------")
        ## load(file = "~/Downloads/impute_develop_alpha.RData")
        ##nMaxDH <- 2 ** 8 - 1

        for(i_option in 1:2) {

            return_betaHat_t <- FALSE
            return_dosage <- FALSE
            return_gamma_t <- FALSE
            return_gammaSmall_t <- FALSE
            make_plots <- FALSE
            return_gammaSmall_t <- FALSE

            if (i_option == 1) {
                only_store_alpha_at_gamma_small <- TRUE
                get_best_haps_from_thinned_sites <- TRUE
                return_dosage <- FALSE
                best_haps_stuff_list1 <- as.list(1:sum(full_gammaSmall_cols_to_get >= 0))
                best_haps_stuff_list2 <- as.list(1:sum(full_gammaSmall_cols_to_get >= 0))
                print("--------- getting best K ---------")
            } else if (i_option == 2) {
                only_store_alpha_at_gamma_small <- FALSE
                get_best_haps_from_thinned_sites <- FALSE
                return_dosage <- TRUE
                best_haps_stuff_list1 <- list()
                best_haps_stuff_list2 <- list()
                print("--------- getting dosage ---------")
                return_gamma_t <- FALSE
            }
            if (return_gamma_t == TRUE) {
                gamma_t1 <- array(0, c(K, nGrids))
                gamma_t2 <- array(0, c(K, nGrids))
                only_store_alpha_at_gamma_small <- FALSE
            } else {
                gamma_t1 <- array(0, c(1, 1))
                gamma_t2 <- array(0, c(1, 1))
                full_gamma_t <- array(0, c(1, 1))
            }


            full_betaHat_t <- array(0, c(1, 1))
            dosage1 <- numeric(nSNPs)
            dosage2 <- numeric(nSNPs)


            ##
            ## run forward to initialize
            ##
            alphaHat_t <- array(0, c(K, nGrids))
            alphaHat_t[, 1] <- first_alpha
            c <- array(1, c(1, nGrids))
            c[1] <- first_c

            Rcpp_haploid_reference_single_forward_version2(
                full_gammaSmall_cols_to_get,
                gl,
                alphaHat_t,
                c,
                transMatRate_t,
                rhb_t,
                hapMatcher,
                eMatDH_bigger,
                nGrids,
                nSNPs,
                K,
                use_eMatDH,
                ref_error,
                only_store_alpha_at_gamma_small = only_store_alpha_at_gamma_small,
                always_normalize = FALSE,
                min_emission_prob_normalization_threshold = 1e-100,
                eMatDH_special_grid_which = eMatDH_special_grid_which,
                eMatDH_special_values_list = eMatDH_special_values_list
            )


            print("--- backwards new version ---------------")
            nReps <- 10
            print(system.time({
                for(i in 1:nReps) {
            Rcpp_haploid_reference_single_backward_version2(
                alphaHat_t,
                betaHat_t,
                gamma_t1,
                gammaSmall_t,
                best_haps_stuff_list1,
                full_gammaSmall_cols_to_get,
                dosage1,
                nGrids,
                transMatRate_t,
                eMatDH_bigger,
                hapMatcher,
                nSNPs,
                K,
                use_eMatDH,
                rhb_t,
                ref_error,
                gl,
                c,
                distinctHapsIE,
                return_betaHat_t,
                return_dosage,
                return_gamma_t,
                return_gammaSmall_t,
                get_best_haps_from_thinned_sites,
                nMaxDH,
                K_top_matches = 5,
                eMatDH_special_grid_which = eMatDH_special_grid_which,
                eMatDH_special_values_list = eMatDH_special_values_list
            )
                }
           }))


            alphaHat_t <- array(0, c(K, nGrids))
            alphaHat_t[, 1] <- first_alpha
            c <- array(1, c(1, nGrids))
            c[1] <- first_c

            Rcpp_haploid_reference_single_forward(
                full_gammaSmall_cols_to_get,
                gl,
                alphaHat_t,
                c,
                transMatRate_t,
                rhb_t,
                hapMatcher,
                eMatDH_bigger,
                nGrids,
                nSNPs,
                K,
                use_eMatDH,
                ref_error,
                only_store_alpha_at_gamma_small,
                always_normalize = FALSE,
                min_emission_prob_normalization_threshold = 1e-100
            )


            print("--- backwards old version ---------------")
            print(system.time({
                for(i in 1:nReps) {
            Rcpp_haploid_reference_single_backward(
                alphaHat_t,
                betaHat_t,
                gamma_t2,
                gammaSmall_t,
                best_haps_stuff_list2,
                full_gammaSmall_cols_to_get,
                dosage2,
                nGrids,
                transMatRate_t,
                eMatDH_bigger,
                hapMatcher,
                nSNPs,
                K,
                use_eMatDH,
                rhb_t,
                ref_error,
                gl,
                c,
                distinctHapsIE,
                return_betaHat_t,
                return_dosage,
                return_gamma_t,
                return_gammaSmall_t,
                get_best_haps_from_thinned_sites,
                nMaxDH,
                K_top_matches = 5
            )
                }
            }))

            if (i_option == 1) {
                expect_equal(
                    sapply(best_haps_stuff_list1, function(x) x[[1]]),
                    sapply(best_haps_stuff_list2, function(x) x[[1]])
                )
            } else {
                ##expect_equal(dosage1, dosage2)
            }
            if (return_gamma_t) {
                ##max(abs(colSums(gamma_t1) - 1)) ## yuck
                ##max(abs(colSums(gamma_t2) - 1)) ## yup
                ## expect_equal(max(abs(colSums(gamma_t2) - 1)))
                ## expect_equal(sum(colSums(gamma_t2) != 1), 0)
                ## expect_equal(sum(colSums(gamma_t2) != 1), 0)
                ## expect_equal(gamma_t1[1:100, ], gamma_t2[1:100, ])
            }

        }

    }


    stop("WER")



    if (1 == 1) {

        alphaHat_t <- array(0, c(K, nGrids))
        betaHat_t <- array(0, c(1, 1))
        gamma_t <- array(0, c(1, 1))
        ww <- seq(1, nGrids, length.out = max(1, round(heuristic_match_thin * nGrids)))
        gammaSmall_cols_to_get <- array(-1, nGrids)
        gammaSmall_cols_to_get[ww] <- 0:(length(ww) - 1)
        gammaSmall_t <-  array(0, c(K, length(ww)))

        ## what type of iteration
        for(i_what in 1:2) {

            if (i_what == 1) {
                print("----------store gammas at subset only -----------")
                return_gammaSmall_t <- TRUE
                return_dosage <- FALSE
                return_good_haps <- TRUE
                get_best_haps_from_thinned_sites <- TRUE
                best_haps_stuff_list <- as.list(1:sum(gammaSmall_cols_to_get >= 0))
            } else {
                print("----------calculate dosages -----------")
                return_gammaSmall_t <- FALSE
                return_dosage <- TRUE
                return_good_haps <- FALSE
                best_haps_stuff_list <- list()
                get_best_haps_from_thinned_sites <- FALSE
            }

            for(i_version in 2:1) {

                i_hap <- 1
                ## for(i_hap in 1:1) {
                alphaHat_t[] <- 0
                gammaSmall_t[] <- 0
                use_eMatDH <- TRUE
                ##
                u <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
                bq <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
                w <- bq != 0
                bq <- bq[w]
                u <- u[w]
                if (length(u) == 0) {
                    print_message(paste0("Read label assignment includes no reads for haplotype ", i_hap))
                }
                ##
                gl <- make_gl_from_u_bq(u, bq, nSNPs)
                ##
                dosage <- numeric(nSNPs)
                if (i_version == 2) {
                    is_version_2 <- TRUE
                    print(paste0("----------- new version --------------"))
                } else {
                    is_version_2 <- FALSE
                    print(paste0("----------- old version --------------"))
                }
                Rcpp_haploid_dosage_versus_refs(
                    gl = gl,
                    alphaHat_t = alphaHat_t,
                    betaHat_t = betaHat_t,
                    gamma_t = gamma_t,
                    dosage = dosage,
                    transMatRate_t = transMatRate_t,
                    rhb_t = rhb_t,
                    ref_error = ref_error,
                    use_eMatDH = use_eMatDH,
                    distinctHapsB = distinctHapsB,
                    distinctHapsIE = distinctHapsIE,
                    hapMatcher = hapMatcher,
                    suppressOutput = 0,
                    gammaSmall_t = gammaSmall_t,
                    return_dosage = return_dosage,
                    return_betaHat_t = FALSE,
                    return_gamma_t = FALSE,
                    return_gammaSmall_t = return_gammaSmall_t,
                    gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                    get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
                    best_haps_stuff_list = best_haps_stuff_list,
                    eMatDH_special_grid_which = eMatDH_special_grid_which,
                    eMatDH_special_values_list = eMatDH_special_values_list,
                    K_top_matches = 5
                )
                if (i_version == 1) {dosage1 <- dosage; gammaSmall_t1 <- gammaSmall_t}
                if (i_version == 2) {dosage2 <- dosage; gammaSmall_t2 <- gammaSmall_t}

            }

            if (return_dosage) {
                expect_equal(dosage1, dosage2)
                expect_true(mean(dosage1) != 0)
            } else {
                expect_equal(gammaSmall_t1, gammaSmall_t2)
                expect_true(mean(gammaSmall_t2) != 0)
            }

        }



    }

    if (1 == 0) {

        m <- apply(alphaHat_t, 2, function(x) length(unique(x)))
        m2 <- mclapply(1:(ncol(alphaHat_t) - 1), mc.cores = 4, function(icol) {
            a <- unique(alphaHat_t[, icol])
            b <- unique(alphaHat_t[, icol + 1])
            c <- paste0(match(alphaHat_t[, icol], a), "-", match(alphaHat_t[, icol + 1], b))
            length(unique(c))
            ## sum(table(alphaHat_t[, icol], alphaHat_t[, icol + 1]) != 0)
        })
        m2 <- unlist(m2)

        m3 <- cbind(
            diff = diff(m),
            max = m2,
            m[-1] - m2,
            m2 - m[-length(m)]
        )
        cbind(
            m,
            rbind(m3, c(NA, NA, NA, NA))
        )[1:40, ]

        ## AM HERE
        ## HOW MANY ARE CREATED AND DESTROYED
        ## CAN I CONSIDER A REDUCED REPRESENTATION
        ## WHERE I UPDATE THESE ON EVERY GO
        ## AND DO FEWER CALCULATIONS
            ## AM

        ## OK, so how may created / lost per go

        ## OK so sometimes several thousands, and sometimes many fewer
        ## seriously can I use lower precision on the values I store...

        ## so often 50+
        x <- apply(hapMatcher, 2, function(x) length(unique(x)))
        sort(x)
        mean(x) ## often 80 ish haplotypes
        y <- apply(eMatDH, 2, function(x) length(unique(x)))
        sort(y)
        mean(y)
        mean(y[y != 1]) ## so more like 10 options on average

        ##
        ## use some sweet sweet real data
        ##
        iGrid <- 100
        ## sort(unique(exp(round(log10(alphaHat_t[, iGrid]), 3))))
        prev_alphaHat_col <- alphaHat_t[, iGrid]
        prev_alphaHat_master <- unique(alphaHat_t[, iGrid])
        prev_alphaHat_class <- match(prev_alphaHat_col, prev_alphaHat_master)
        ##
        ## next one
        ##
        dh_col <- hapMatcher[, iGrid + 1]
        eMatDH_col <- eMatDH[, iGrid + 1]
        prob_col <- eMatDH_col[dh_col]
        ##
        ## original way to do it
        ##
        jump_prob <- full_transMatRate_t_H[2, iGrid] / K
        not_jump_prob <- full_transMatRate_t_H[1, iGrid]
        new_alphaHat_col_old <- (jump_prob + not_jump_prob * prev_alphaHat_col) * prob_col
        ##
        ## alternate way to do it
        ##
        class_compare <- array(0, c(K, nMaxDH + 1))
        to_blank <- array(0, c(K, 2))
        to_reclass <- array(0, c(K, nMaxDH + 1))
        to_reclass_count <- array(0, c(K))
        count <- 1
        for(k in 1:K) {
            if (class_compare[prev_alphaHat_class[k], dh_col[k]] == 0) {
                class_compare[prev_alphaHat_class[k], dh_col[k]] <- 1
                to_blank[count, 1] <- prev_alphaHat_class[k]
                to_blank[count, 2] <- dh_col[k]
                count <- count + 1
            }
        }
        ## now calculate these values only, in these locations
        for(i in 1:(count - 1)) {
            alpha_val <- prev_alphaHat_master[to_blank[i, 1]]
            prob_val <- eMatDH_col[to_blank[i, 2]]
            class_compare[to_blank[i, 1], to_blank[i, 2]] <-
                (jump_prob + not_jump_prob * alpha_val) * prob_val
        }
        ## now calculate new alphaHat column from this
        new_alphaHat_col_new <- numeric(K)
        for(k in 1:K) {
            new_alphaHat_col_new[k] <- class_compare[prev_alphaHat_class[k], dh_col[k]]
        }
        ## now figure out which ones have been used
        for(i in 1:(count - 1)) {
            ac <- to_blank[i, 1]
            dh <- to_blank[i, 2]
            to_reclass_count[ac] <- to_reclass_count[ac] + 1
            to_reclass[ac, to_reclass_count[ac]] <- dh
        }
        ## now re-class those that have been split?

        ## OK, from this, for those split, get split

        expect_equal(new_alphaHat_col_old, new_alphaHat_col_new)
        ## now need to re-set class
        for(i in 1:(count - 1)) {
            class_compare[to_blank[i, 1], to_blank[i, 2]] <- 0
        }
        ## OK this is now blank, that is good. count re-sets itself
        expect_equal(sum(class_compare != 0), 0)
        ## now need to update class for next round
        ##








    }

})












test_that("prototype idea of classes for faster alphaHat", {

    skip("not sure this is the right idea")


    ## so have initial class
    ## and then have prob, some change
    ## re-make class in efficient way after calc
    set.seed(1)
    K <- 1000
    jump_prob <- 0.01 / K
    not_jump_prob <- 1 - 0.01
    ## previous class
    prev_class <- sample(1:10, K, replace = TRUE)
    prev_alphaHat_master <- runif(10)
    prev_alphaHat_master <- prev_alphaHat_master / 10 ## normalized
    prev_alphaHat_col <- prev_alphaHat_master[prev_class]
    ## now have 4 options for reference haplotype, ignore 0 value for now
    dh_col <- sample(1:4, K, replace = TRUE)
    eMatDH_col <- runif(4) ## and their outputs
    prob_col <- eMatDH_col[dh_col]

    ##
    ## original way to do it
    ##
    new_alphaHat_col <- (jump_prob + not_jump_prob * prev_alphaHat_col) * prob_col

    ##
    ## potential new way to do it
    ##
    ## so almost all the mass in low values, remember, these are ordered





})





