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

