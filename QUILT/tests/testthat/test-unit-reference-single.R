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

test_that("can run a single gl sample through reference haplotypes quickly with grid of 32", {

    ## switch on and off. off does thes tests, on checks cpp times
    speed_test <- FALSE
    if (!speed_test) {
        nSNPs <- 100
        L <- sort(sample(1:1000, 100))
    } else {
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
    expRate <- 1
    nGen <- 10
    sigmaCurrent <- exp(-nGen * expRate / 100 / 1000000 * dl)
    transMatRate_t <- rbind(sigmaCurrent, 1 - sigmaCurrent)


    use_eMatDH <- TRUE; i_setup <- 2

    gammaSmall_cols_to_get <- c(-1, 0, -1, 1) ## i.e. 0-based, what to put it in
    nSmallGammaGrids <- 2    
    
    ## build some haplotypes and encode them
    for(use_eMatDH in c(FALSE, TRUE)) {
        for(i_setup in 1:2) {

            print(paste0("use_eMatDH = ", use_eMatDH, ", i_setup = ", i_setup))
            if (i_setup == 1) {
                ## test gamma makes sense
                set.seed(9910)
                K <- 100
                reference_haps <- array(as.integer(runif(nSNPs * K) > 0.5), c(nSNPs, K))
            } else {
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
            }
            
            ## make them mostly one of three options, with a few small changes
            rhi <- reference_haps
            rhi_t <- t(rhi)
            rhb_t <- STITCH::make_rhb_t_from_rhi_t(rhi_t)
            rhb <- t(rhb_t)

            nMaxDH <- 2 ** 10 - 1
            ref_error <- 0.01
            ref_one_minus_error <- 1 - ref_error

            
            ## make haplotype matching objects
            out <- make_rhb_t_equality(
                rhb_t = rhb_t,
                nMaxDH = nMaxDH,
                nSNPs = nSNPs,
                ref_error = ref_error
            )
            distinctHapsB <- out[["distinctHapsB"]]
            distinctHapsIE <- out[["distinctHapsIE"]]            
            hapMatcher <- out[["hapMatcher"]]
            
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
            u <- sort(sample(1:nSNPs, nReads, replace = TRUE))
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
            gamma_t <- array(0, c(K, nGrids))
            gammaSmall_t <- array(0, c(K, nSmallGammaGrids))
            dosage <- numeric(nSNPs)
            
            ## so now, want to 
            if (!speed_test) {
                
                ## 
                a <- system.time(
                    outR <- R_haploid_dosage_versus_refs(
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
                        return_extra = TRUE
                    )
                )
                
                print("R speed")
                ## print(a)
                R_gamma_t <- outR[["gamma_t"]]
                R_dosage <- outR[["dosage"]]
                R_alphaHat_t <- outR[["alphaHat_t"]]
                R_betaHat_t <- outR[["betaHat_t"]]    
                expect_equal(colSums(R_gamma_t), rep(1, 4))
                expect_true((cor(R_dosage, my_hap) ** 2) > 0.7) ## is better for truer match against reference
                ## alternate
                
            }
            
            a <- system.time(
            Rcpp_haploid_dosage_versus_refs(gl = gl,alphaHat_t = alphaHat_t,betaHat_t = betaHat_t,gamma_t = gamma_t,dosage = dosage,transMatRate_t = transMatRate_t,rhb_t = rhb_t,ref_error = ref_error,use_eMatDH = use_eMatDH,distinctHapsB = distinctHapsB,distinctHapsIE = distinctHapsIE, hapMatcher = hapMatcher,gammaSmall_t = gammaSmall_t,gammaSmall_cols_to_get = gammaSmall_cols_to_get, suppressOutput = 0)
            )
            print("cpp speed")
            print(a)

            ## also check components are the same
            if (use_eMatDH) {

                ## check can just minimally get dosage
                dosageX <- numeric(nSNPs)
                Rcpp_haploid_dosage_versus_refs(
                    gl = gl, transMatRate_t = transMatRate_t, rhb_t = rhb_t,ref_error = ref_error,use_eMatDH = use_eMatDH,distinctHapsB = distinctHapsB,distinctHapsIE = distinctHapsIE, hapMatcher = hapMatcher, alphaHat_t = alphaHat_t,
                    betaHat_t = array(0, c(1, 1)),
                    gamma_t = array(0, c(1, 1)),
                    dosage = dosageX,
                    gammaSmall_t = gammaSmall_t,
                    return_dosage = TRUE,
                    return_betaHat_t = FALSE,
                    return_gamma_t = FALSE,
                    return_gammaSmall_t = FALSE,
                    gammaSmall_cols_to_get = array(0, c(1, 1))
                )
                expect_equal(dosage, dosageX)
                ## check can just minimally get small gamma
                gammaSmall_tX <- array(0, c(K, nSmallGammaGrids))
                Rcpp_haploid_dosage_versus_refs(
                    gl = gl, transMatRate_t = transMatRate_t, rhb_t = rhb_t,ref_error = ref_error,use_eMatDH = use_eMatDH,distinctHapsB = distinctHapsB,distinctHapsIE = distinctHapsIE, hapMatcher = hapMatcher, alphaHat_t = alphaHat_t,
                    betaHat_t = array(0, c(1, 1)),
                    gamma_t = array(0, c(1, 1)),
                    dosage = numeric(1),
                    gammaSmall_t = gammaSmall_tX,
                    return_dosage = FALSE,
                    return_betaHat_t = FALSE,
                    return_gamma_t = FALSE,
                    return_gammaSmall_t = TRUE,
                    gammaSmall_cols_to_get = gammaSmall_cols_to_get
                )
                expect_equal(gammaSmall_tX, R_gamma_t[, gammaSmall_cols_to_get >= 0])
                
                

            }
            
            if (!speed_test) {
                print(use_eMatDH)
                expect_equal(R_alphaHat_t, alphaHat_t)
                expect_equal(R_betaHat_t, betaHat_t)
                expect_equal(R_gamma_t, gamma_t)
                expect_equal(R_dosage, dosage)
            }

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


test_that("profile", {

    if (1 == 0) {

        load("/well/davies/users/dcc832/temp.RData")

        nMaxDH <- 2 ** 10 - 1
        ref_error <- 1e-2
        out <- make_rhb_t_equality(
            rhb_t = rhb_t,
            nMaxDH = nMaxDH,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
        distinctHapsB <- out[["distinctHapsB"]]
        distinctHapsIE <- out[["distinctHapsIE"]]            
        hapMatcher <- out[["hapMatcher"]]
        save(distinctHapsB, distinctHapsIE, hapMatcher, nSNPs, transMatRate_t, gl, rhb_t, file = "/well/davies/users/dcc832/temp.RData", compress = FALSE)    

    }

    skip("wer")

    load("/well/davies/users/dcc832/temp.RData")
    K <- nrow(rhb_t)    
    nGrids <- ncol(rhb_t)
    nMaxDH <- 2 ** 10 - 1
    ref_error <- 1e-2
    i <- 1
    
    for(i in 1:1) {

        if (i == 1) {use_eMatDH <- TRUE}
        print(paste0("----------i = ", i))
        alphaHat_t <- array(0, c(K, nGrids))
        ##betaHat_t <- array(0, c(K, nGrids))
        ##gamma_t <- array(0, c(K, nGrids))
        betaHat_t <- array(0, c(1, 1))
        gamma_t <- array(0, c(1, 1))
        gammaSmall_t <- array(0, c(1, 1))
        a <- system.time(
            for(i in 1:2) {
            dosage <- numeric(nSNPs)
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
                return_dosage = FALSE,
                return_betaHat_t = FALSE,
                return_gamma_t = FALSE,
                return_gammaSmall_t = FALSE,
                gammaSmall_cols_to_get = array(0, c(1, 1))
            )
            }
        )
        
        if (i == 1) {dosage1 <- dosage}
        if (i == 2) {dosage2 <- dosage}
        print(a)
    }
    expect_equal(dosage1, dosage2)
    ## yup, faster. 5.2 seconds to 4.5
    
    which(is.na(dosage))

    ## yup, faster, thought not massively
    ## 1 - 
    ## 2 - 
    ## 3 - 
    
    ## OK, good start!
    ## can I use vector extraction, get this down more

    
})

    









test_that("prototype idea of classes for faster alphaHat", {

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
