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




test_that("can avoid normalizing alphaHat and betaHat throughout forward algorithm, but recover correct gamma when asked, as well as dosages", {


    ##
    ## setup
    ##
    ##
    speed_test <- FALSE
    rebuild <- FALSE
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
    file <- paste0("test_package.", speed_test, ".RData")
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
        message("loading")
        load(file)
    }
    distinctHapsB <- test_package[["distinctHapsB"]]
    distinctHapsIE <- test_package[["distinctHapsIE"]]
    hapMatcher <- test_package[["hapMatcher"]]
    rhb_t <- test_package[["rhb_t"]]
    gl <- test_package[["gl"]]
    transMatRate_t <- test_package[["transMatRate_t"]]
    ref_error <- test_package[["ref_error"]]
    eMatDH_special_grid_which <- test_package[["eMatDH_special_grid_which"]]
    eMatDH_special_values_list <- test_package[["eMatDH_special_values_list"]]
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
    language <- "R"
    min_emission_prob_normalization_threshold <- 1e-100 ## make much smaller than usual
     
    ##
    ## run normal version in R, or Rcpp, and make sure the same
    ##
    master_f <- function(always_normalize, language, is_version_2) {
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
                    "is_version_2 = ", is_version_2
                )
            )
        }
        ##
        alphaHat_t <- array(0, c(K, nGrids))
        betaHat_t <- array(0, c(K, nGrids))
        c <-  array(1, c(nGrids))
        gamma_t <- array(0, c(K, nGrids))
        gammaSmall_t <- array(0, c(K, nSmallGammaGrids))
        dosage <- numeric(nSNPs)
        best_haps_stuff_list <- as.list(1:sum(gammaSmall_cols_to_get >= 0))
        ##
        out <- f(
            gl = gl,
            alphaHat_t = alphaHat_t,
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
            K_top_matches = K_top_matches,
            always_normalize = always_normalize,
            min_emission_prob_normalization_threshold = min_emission_prob_normalization_threshold,
            best_haps_stuff_list = best_haps_stuff_list,
            eMatDH_special_values_list = eMatDH_special_values_list,
            eMatDH_special_grid_which = eMatDH_special_grid_which,
            is_version_2 = is_version_2,
            suppressOutput = suppressOutput
        )
        if (language == "R") {
            alphaHat_t <- out[["alphaHat_t"]]
            c <- out[["c"]]
            gammaSmall_t <- out[["gammaSmall_t"]]
            dosage <- out[["dosage"]]
            betaHat_t <- out[["betaHat_t"]]
        }
        return(
            list(
                alphaHat_t = alphaHat_t,
                betaHat_t = betaHat_t,
                c = c,
                gammaSmall_t = gammaSmall_t,
                dosage = dosage
            )
        )
    }
    out_R_always_normalize <- master_f(TRUE, "R", NA)
    out_Rcpp_always_normalize <- master_f(TRUE, "Rcpp", FALSE)
    out_Rcpp2_always_normalize <- master_f(TRUE, "Rcpp", TRUE)
    
    out_R_seldom_normalize <- master_f(FALSE, "R", NA)
    out_Rcpp_seldom_normalize <- master_f(FALSE, "Rcpp", FALSE)
    out_Rcpp2_seldom_normalize <- master_f(FALSE, "Rcpp", TRUE)

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
        ##
        ## dosages
        ##
        expect_equal(out_R_always_normalize[["dosage"]], out_Rcpp2_always_normalize[["dosage"]])
        expect_equal(out_R_always_normalize[["dosage"]], out_Rcpp2_seldom_normalize[["dosage"]])
        expect_equal(out_Rcpp_always_normalize[["dosage"]], out_Rcpp2_always_normalize[["dosage"]])
        ##
        ## gammas
        ##
        expect_equal(out_R_always_normalize[["gammaSmall_t"]], out_Rcpp2_always_normalize[["gammaSmall_t"]])
        expect_equal(out_R_always_normalize[["gammaSmall_t"]], out_Rcpp2_seldom_normalize[["gammaSmall_t"]])
        expect_equal(out_Rcpp_always_normalize[["gammaSmall_t"]], out_Rcpp2_always_normalize[["gammaSmall_t"]])
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
