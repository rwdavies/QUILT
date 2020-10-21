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




test_that("can avoid normalizing alphaHat throughout forward backwards, but recover correct gamma when asked, as well as dosages", {


    ##
    ## setup
    ##
    ##
    ##
    test_package <- make_reference_single_test_package(
        nSNPs = 500,
        K = 100
    )
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
    
    for(always_normalize in c(TRUE, FALSE)) {

        for(language in c("R", "Rcpp")) {            

            if (language == "R") {
                f <- R_haploid_dosage_versus_refs
            } else {
                f <- Rcpp_haploid_dosage_versus_refs
            }

            ##
            alphaHat_t <- array(0, c(K, nGrids))
            betaHat_t <- array(0, c(K, nGrids))
            c <-  array(1, c(nGrids))
            gamma_t <- array(0, c(K, nGrids))
            gammaSmall_t <- array(0, c(K, nSmallGammaGrids))
            dosage <- numeric(nSNPs)
            best_haps_stuff_list <- as.list(1:sum(gammaSmall_cols_to_get >= 0))

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
                return_gamma_t = FALSE,
                K_top_matches = K_top_matches,
                always_normalize = always_normalize,
                min_emission_prob_normalization_threshold = min_emission_prob_normalization_threshold,
                best_haps_stuff_list = best_haps_stuff_list,
                eMatDH_special_values_list = eMatDH_special_values_list,
                eMatDH_special_grid_which = eMatDH_special_grid_which
            )

            if (always_normalize) {
                if (language == "R") {
                    R_alphaHat_t_always_norm <- out[["alphaHat_t"]]                
                    R_c_always_norm <- out[["c"]]
                    R_gammaSmall_t_always_norm <- out[["gammaSmall_t"]]
                    R_dosage_always_norm <- out[["dosage"]]
                } else {
                    Rcpp_dosage_always_norm <- dosage
                    Rcpp_gammaSmall_t_always_norm <- gammaSmall_t
                    Rcpp_c_always_norm <- c
                }
            } else {
                if (language == "R") {
                    R_alphaHat_t_seldom_norm <- out[["alphaHat_t"]]
                    R_c_seldom_norm <- out[["c"]]
                    R_gammaSmall_t_seldom_norm <- out[["gammaSmall_t"]]
                    R_dosage_seldom_norm <- out[["dosage"]]
                } else {
                    Rcpp_dosage_seldom_norm <- dosage
                    Rcpp_gammaSmall_t_seldom_norm <- gammaSmall_t
                    Rcpp_c_seldom_norm <- c                    
                }
            }

        }

    }

    ##
    ## check results for R versions
    ##
    expect_equal(R_dosage_always_norm, R_dosage_seldom_norm)
    expect_equal(R_gammaSmall_t_always_norm, R_gammaSmall_t_seldom_norm)
    expect_equal(sum(log(R_c_always_norm)), sum(log(R_c_seldom_norm)))
    expect_true(sum(abs(R_c_always_norm - R_c_seldom_norm)) > 0)

    ##
    ## check cpp version here too
    ##
    expect_equal(sum(log(Rcpp_c_always_norm)), sum(log(Rcpp_c_seldom_norm)))
    ## expect_true(sum(abs(Rcpp_c_always_norm - Rcpp_c_seldom_norm)) > 0)
    

    ## check cpp versions here!

    if (1 == 0) {

        cbind(log(c_always_norm), log(c_seldom_norm))

        colSums(gammaSmall_t_R_always_norm)
        colSums(gammaSmall_t_R_seldom_norm)
        
        iGrid <- 3
        x <- alphaHat_t_always_norm[, iGrid] / prod(c_always_norm[1:iGrid])
        y <- alphaHat_t_seldom_norm[, iGrid] / prod(c_seldom_norm[1:iGrid])
        table(round(log10(x / y), 2))
        head(cbind(x, y))

    }
    
})

