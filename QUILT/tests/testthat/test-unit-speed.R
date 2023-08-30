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


test_that("speed test for cpp stuff", {

    skip("wip")
    
    out <- make_quilt_fb_test_package(
        K = 400,
        nReads = 32000,
        nSNPs = 32000,
        S = 1,
        gridWindowSize = 32
    )

    sampleReads <- out$sampleReads
    eMatRead_t <- out$list_of_eMatRead_t[[1]]    
    K <- out$K
    nGrids <- out$nGrids
    alphaHat_t1 <- array(0, c(K, nGrids))
    betaHat_t1 <- array(0, c(K, nGrids))
    alphaHat_t2 <- array(0, c(K, nGrids))
    betaHat_t2 <- array(0, c(K, nGrids))

    print(microbenchmark::microbenchmark(
                    rcpp_test1(
                        alphaHat_t1,
                        alphaHat_t2,
                        betaHat_t1,
                        betaHat_t2,
                        sampleReads,
                        eMatRead_t    
                    ),
                    rcpp_test2(
                        alphaHat_t1,
                        alphaHat_t2,
                        betaHat_t1,
                        betaHat_t2,
                        sampleReads,
                        eMatRead_t    
                    )           ,
                    rcpp_test3(
                        alphaHat_t1,
                        alphaHat_t2,
                        betaHat_t1,
                        betaHat_t2,
                        sampleReads,
                        eMatRead_t    
                    )           ,
                    rcpp_test4(
                        alphaHat_t1,
                        alphaHat_t2,
                        betaHat_t1,
                        betaHat_t2,
                        sampleReads,
                        eMatRead_t
                    )           ,
                    rcpp_test5(
                        alphaHat_t1,
                        alphaHat_t2,
                        betaHat_t1,
                        betaHat_t2,
                        sampleReads,
                        eMatRead_t
                    )           ,
                    times = 5L))

    

})
