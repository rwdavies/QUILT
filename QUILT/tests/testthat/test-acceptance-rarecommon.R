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


set.seed(919)

n_snps <- 1200
K <- 6
n_big_haps <- 1000
chr <- 10


## make a spectrum or rare and common SNPs
phasemaster1 <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
phasemaster2 <- phasemaster1[, sample(1:K, n_big_haps, replace = TRUE)]

## make about 10% common, and do not touch
n_common <- round(n_snps / 10)
common <- sample(1:n_snps, n_common)
## make about 90% rare, and make about 1% freq
phasemaster2[-common, ] <- 0
## make remaining frequency about 1%
phasemaster2[-common, ] <- sample(c(0, 1), (n_snps - n_common) * n_big_haps, prob = c(0.99, 0.01), replace = TRUE)


reads_span_n_snps <- 3
## want about 4X here
n_reads <- round(4 * n_snps / reads_span_n_snps)
data_package <- STITCH::make_acceptance_test_data_package(
    reads_span_n_snps = reads_span_n_snps,
    n_samples = 3,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 2,
    chr = chr,
    K = K,
    phasemaster = phasemaster2
)
refpack <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = n_big_haps,
    reference_populations = c("GBR"),
    chr = chr,
    phasemaster = phasemaster2
)






test_that("can parse and use input hap VCF only, and use rare vs common idea", {

    
    af_cutoff <- 0.01 ## here, 1%. note maf < 1% corresponding to af > 99% is very rare so could ignore for simplicity?
    ## thinds robbie needs for gibbs
    use_hapMatcherR <- TRUE
    nMaxDH <- NA ## automatically inferred
    ref_error <- 0.001

    
    ##
    ## write data parser
    ## input: phased vcf, maf threshold
    ## output:
    ##  - at common sites: zilong indices, robbie things (e.g. hapMatcher)
    ##  - at rare sites: an object that somehow can be used to rebuild

    ##
    ## this stuff probably in C++, streaming, efficient, etc
    ## question: can we incorporate skipping of variants we don't like e.g. if multi-allelic, if position seen before, if not SNPs, etc
    ## this would be very desirable

    
    ## todo: convert this to nice streamable vcf parsing somehow?
    haps <- read.table(refpack$reference_vcf_file, stringsAsFactors = FALSE)
    LAll <- haps[, 2]
    h <- haps[, -(1:9)]
    h <- as.matrix(h)
    K <- ncol(h) * 2
    nSNPsAll <- nrow(h)
    rhi <- array(as.integer(NA), c(nSNPsAll, K)) ## reference haplotype integers
    rhi[, seq(1, K, 2)] <- matrix(as.integer(substr(h, 1, 1)), nrow = nSNPsAll)
    rhi[, seq(2, K, 2)] <- matrix(as.integer(substr(h, 3, 3)), nrow = nSNPsAll)
    ## filtering 
    ref_alleleCount <- cbind(rowSums(rhi), K, rowSums(rhi) / K)
    snp_is_common <- !(ref_alleleCount[, 3] < af_cutoff)
    L <- LAll[snp_is_common]
    nSNPs <- sum(snp_is_common)
    nGrids <- ceiling(nSNPs / 32)
    grid <- as.integer(floor((1:nSNPs) / 32))

    ## not sure of the right way to do this on the fly? with rhb or rhb_t?
    rhb <- matrix(0L, nrow = nGrids, ncol = K)
    for(iGrid in 1:nGrids) {
        s <- min(32 * (iGrid - 1) + 1, nSNPs)
        e <- min(32 * iGrid, nSNPs)
        for(k in 1:K) {
            rhb[iGrid, k] <- rcpp_int_contract(rhi[s:e, k])
        }
    }
    ## transpose (do after? during? probably can just work with rhb_t during)
    rhb_t <- t(rhb)





    
    ##
    ## ????? how to store rare stuff ?????
    ## ? how to make properly RAM efficient?
    rare_per_hap_info <- sapply(1:K, function(k) {
        which(rhi[!snp_is_common, k] == 1)
    })

    


    ## define mapping of all SNPs into grid
    ## grid is defined for common SNPs
    ## 
    gridAll <- as.integer(rep(NA, nSNPsAll))
    ## gridAll[snp_is_common] <- grid
    which_snp_is_common <- which(snp_is_common)
    for(iGrid in 1:nGrids) {
        ## start and end (for common)
        s <- min(32 * (iGrid - 1) + 1, nSNPs) ## first common SNP in grid
        e <- min(32 * iGrid, nSNPs) ## last common SNP in grid
        ## all SNPs in between
        sa <- which_snp_is_common[s] ## start in all set
        ea <- which_snp_is_common[e] ## end in all set
        ## all these can fill in
        gridAll[sa:ea] <- as.integer(iGrid - 1)
        if (iGrid == 1) {
            gridAll[1:sa] <- 0L
        } else if (iGrid == nGrids) {
            gridAll[ea:nSNPsAll] <- as.integer(nGrids - 1)
        }
        if (iGrid > 1) {
            ## check and fill backwards, if needed
            if ((sa - prev_ea) != 1) {
                ## need to fill in, just do to closer end, based on physical distance
                for(iSNP in (prev_ea + 1):(sa - 1)) {
                    stopifnot(is.na(gridAll[iSNP]))
                    if ((LAll[e] - LAll[iSNP]) > (LAll[iSNP] > LAll[s])) {
                        gridAll[iSNP] <- as.integer(iGrid - 1)
                    } else {
                        gridAll[iSNP] <- as.integer(iGrid - 2)
                    }
                }
            }
        }
        prev_ea <- ea
    }
    ## 
    gridRare <- gridAll[!snp_is_common]
    ## also specifically store start and end
    grid_rare_start_and_end <- sapply(tapply(X = 1:length(gridRare), INDEX = gridRare, range), I)

    
    ## 
    ## what I need from rhb_t for the common stuff
    ## 
    out <- make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nSNPs,
        ref_error = ref_error,
        use_hapMatcherR = use_hapMatcherR
    )
    distinctHapsB <- out[["distinctHapsB"]]
    distinctHapsIE <- out[["distinctHapsIE"]]
    hapMatcherR <- out[["hapMatcherR"]]
    eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]
    eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
    rm(rhb_t, out, hapMatcher)
    gc(reset = TRUE); gc(reset = TRUE);













    


    ##
    ## example of what I would do with imputation
    ##
    

    ## use this example Ksubset
    Ksubset <- 40
    which_haps_to_use <- sample(1:K, Ksubset) ## 1-based
    
    ## now suppose we have the following posteriors
    gamma1 <- matrix(runif(nGrids * Ksubset), nrow = Ksubset, ncol = nGrids)
    gamma1 <- gamma1 / rep(colSums(gamma1), each = nrow(gamma1))
    gamma2 <- matrix(runif(nGrids * Ksubset), nrow = Ksubset, ncol = nGrids)
    gamma2 <- gamma2 / rep(colSums(gamma2), each = nrow(gamma2))
    

    ## we can build dosages as follows
    hapDosageRare <- numeric(nSNPsAll - nSNPs)
    hapDosageCommon <- numeric(nSNPs)

    ## ??? inflate haps of interest
    ## ??? probably too slow, could do on the fly somehow down below maybe?
    rhi_t_common <- matrix(integer(1), nrow = Ksubset, ncol = nSNPsAll - nSNPs)
    for(k in 1:Ksubset) {
        rhi_t_common[k, rare_per_hap_info[[which_haps_to_use[k]]]] <- 1L
    }
    
    for(iGrid in 1:nGrids) {
        g1 <- gamma1[, iGrid]
        g2 <- gamma2[, iGrid]
        s <- min(32 * (iGrid - 1) + 1, nSNPs) ## first common SNP in grid
        e <- min(32 * iGrid, nSNPs) ## last common SNP in grid
        nLocalSNPs <- as.integer(e - s + 1)
        ## for rare SNPs
        sr <- grid_rare_start_and_end[1, iGrid] ## start rare
        er <- grid_rare_start_and_end[2, iGrid] ## start rare        
        for(ik in 1:Ksubset) {
            ##
            ## common
            ##
            k <- which_haps_to_use[ik]
            i <- as.integer(hapMatcherR[k, iGrid])
            if (i > 0) {
                ie <- distinctHapsIE[i, s:e]
            } else {
                b <- rcpp_simple_binary_matrix_search( ## test both R and Rcpp versions
                    val = k - 1,
                    mat = eMatDH_special_matrix,
                    s1 = eMatDH_special_matrix_helper[iGrid, 1],
                    e1 = eMatDH_special_matrix_helper[iGrid, 2]
                )
                ie <- rcpp_int_expand(b, nSNPs = nLocalSNPs)
            }
            hapDosageCommon[s:e] <- hapDosageCommon[s:e] + ie
            ##
            ## rare ???
            ##
            hapDosageRare[sr:er] <- hapDosageRare[sr:er] + rhi_t_common[ik, sr:er]
            ##
        }
    }
    hapDosageCommon <- hapDosageCommon / Ksubset
    hapDosageRare <- hapDosageRare / Ksubset

    
    ## combine
    hapDosage <- numeric(nSNPs)
    hapDosage[snp_is_common] <- hapDosageCommon[snp_is_common]
    hapDosageRare[!snp_is_common] <- hapDosageRare[!snp_is_common]    
    
    
})
