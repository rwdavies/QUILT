## just use this one function
## then figure out what I need to do below that, if / as necessary
calculate_pse_and_r2_master <- function(
    method,
    have_truth_haplotypes,
    have_truth_genotypes,
    truth_haps,
    truth_gen,
    hap1 = NULL,
    hap2 = NULL,
    hap3 = NULL,
    impute_rare_common = FALSE,
    checking_all_snps = FALSE,
    verbose = FALSE,
    inRegion2 = NULL,
    af = NULL,
    dosage = NULL,
    mat_dosage = NULL,
    fet_dosage = NULL,
    prefix = ""
) {

    ## if you have neither, leave! nothing to check
    if (!have_truth_haplotypes && !have_truth_genotypes) {
        return(NULL)
    }
    
    out <- NULL
    
    if (method == "diploid") {
        out <- calculate_pse_and_r2_during_gibbs(
            have_truth_haplotypes = have_truth_haplotypes,
            have_truth_genotypes = have_truth_genotypes,
            inRegion2 = inRegion2,
            truth_haps = truth_haps,
            truth_gen = truth_gen,
            hap1 = hap1,
            hap2 = hap2,
            af = af,
            verbose = verbose,
            impute_rare_common = impute_rare_common,
            checking_all_snps = checking_all_snps,
            dosage = dosage,
            prefix = prefix
        )

    } else {

        ## meh don't bother writing just truth gen? kind of weird anyway
        
        if (have_truth_haplotypes) {

            calculate_pse_and_r2_during_gibbs_nipt(
                inRegion2 = inRegion2,
                hap1 = hap1,
                hap2 = hap2,
                hap3 = hap3,
                truth_haps = truth_haps,
                af = af,
                verbose = verbose,
                mat_dosage = mat_dosage,
                fet_dosage = fet_dosage,
                prefix = prefix,
                impute_rare_common = impute_rare_common,
                checking_all_snps = checking_all_snps
            )
        }

    }

    out
    
}





calculate_pse_and_r2_during_gibbs <- function(
    have_truth_haplotypes,
    have_truth_genotypes,
    inRegion2,
    truth_haps,
    truth_gen,
    hap1,
    hap2,
    af,
    verbose = FALSE,
    impute_rare_common = FALSE,
    checking_all_snps = FALSE,
    dosage = NULL,
    prefix = ""
) {
    ## wow, confident
    w <- (inRegion2)
    if (have_truth_haplotypes) {
        g <- truth_haps[inRegion2, 1] + truth_haps[inRegion2, 2]
    } else if (have_truth_genotypes) {
        g <- truth_gen[inRegion2, 1]
    }
    ## scaled version
    if (is.null(dosage)) {
        d <- (hap1 + hap2)[inRegion2]
    } else {
        d <- dosage[inRegion2]
    }
    r2 <-  round(cor(d - 2 * af[inRegion2], g - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
    fake_LL <- 1:sum(w)
    ##
    if (impute_rare_common) {
        if (checking_all_snps) {
            suffix <- " (all SNPs)"
        } else {
            suffix <- " (common SNPs only)"
        }
    } else {
        suffix <- ""
    }
    ## 
    if (have_truth_haplotypes) {
        values <- modified_calculate_pse(test = round(cbind(hap1, hap2)[w, ]), truth = truth_haps[w, ], LL = fake_LL)$values
        pse <- values["phase_errors_def1"] / values["phase_sites_def1"]
        pse <- round(100 * pse, 1)
        disc <- round(100 * values["disc_errors"] / values["dist_n"], 1)
    } else {
        pse <- NA
        disc <- NA
    }
    if (verbose) {
        if (have_truth_haplotypes) {
            print_message(paste0(prefix, "r2:", r2, ", PSE:", pse, "%, disc:", disc, "%", suffix))
        } else {
            print_message(paste0(prefix, "r2:", r2, suffix))            
        }
    }
    return(c(r2 = r2, pse = as.numeric(pse), disc = as.numeric(disc)))
}



## calculate_pse_and_r2_rare_common <- function(
##     hap1_all,
##     hap2_all,
##     have_truth_haplotypes,
##     truth_haps_all,
##     have_truth_genotypes,
##     truth_gen_all,
##     special_rare_common_objects,
##     verbose
## ) {
    
##     ref_alleleCount_all <- special_rare_common_objects[["ref_alleleCount_all"]]
##     af <- ref_alleleCount_all[, 3]
##     inRegion2 <- special_rare_common_objects[["inRegion2"]]
    
##     if (have_truth_haplotypes) {
##         x <- calculate_pse_and_r2_during_gibbs(
##             inRegion2 = inRegion2,
##             hap1 = hap1_all,
##             hap2 = hap2_all,
##             truth_haps = truth_haps_all,
##             af = af,
##             verbose = verbose,
##             impute_rare_common = TRUE,
##             all_snps = TRUE
##         )
##     } else if (have_truth_genotypes) {
##         r2 <-  round(cor((dosage_all)[inRegion2] - 2 * af[inRegion2], truth_gen_all[inRegion2, ] - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
##         print_message(paste0("Current accuracy for this gibbs sample for ", sample_name, " for all SNPs, r2:", r2))
##     }
    
## }





calculate_pse_and_r2_during_gibbs_nipt <- function(
    inRegion2,
    hap1,
    hap2,
    hap3,
    truth_haps,
    af,
    verbose = FALSE,
    mat_dosage = NULL,
    fet_dosage = NULL,
    prefix = "",
    impute_rare_common = FALSE,
    checking_all_snps = FALSE
) {
    if (impute_rare_common) {
        if (checking_all_snps) {
            suffix <- " (all SNPs)"
        } else {
            suffix <- " (common SNPs only)"
        }
    } else {
        suffix <- ""
    }
    to_return <- NULL
    ## wow, confident
    w <- (inRegion2)
    ## scaled version
    fake_LL <- 1:sum(w)
    for(i_hap in 1:2) {
        if (i_hap == 1) {
            who <- "mat "
            hapA <- hap1
            hapB <- hap2
            th1 <- 1
            th2 <- 2
            if (is.null(mat_dosage)) {
                dosage <- hapA + hapB
            } else {
                dosage <- mat_dosage
            }
        } else {
            who <- "fet "
            hapA <- hap1
            hapB <- hap3
            th1 <- 1
            th2 <- 3
            if (is.null(mat_dosage)) {
                dosage <- hapA + hapB
            } else {
                dosage <- fet_dosage
            }
        }
        values <- modified_calculate_pse(test = round(cbind(hapA, hapB)[w, ]), truth = truth_haps[w, c(th1, th2)], LL = fake_LL)$values
        g <- truth_haps[inRegion2, th1] + truth_haps[inRegion2, th2]
        r2 <-  round(cor((dosage)[inRegion2] - 2 * af[inRegion2], g - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
        pse <- round(100 * values["phase_errors_def1"] / values["phase_sites_def1"], 1)
        disc <- round(100 * values["disc_errors"] / values["dist_n"], 1)
        if (verbose) {
            print_message(paste0(prefix, who, "r2:", r2, ", PSE:", pse, "%, disc:", disc, "%", suffix))
        }
        to_return[paste0("r2_", who)] <- r2
        to_return[paste0("pse_", who)] <- pse
        to_return[paste0("disc_", who)] <- disc
    }
    if (verbose) {
        ##
        r2_1 <- round(cor(hap1[inRegion2] - af[inRegion2], truth_haps[inRegion2, 1] - af[inRegion2]) ** 2, 3)
        r2_2 <- round(cor(hap2[inRegion2] - af[inRegion2], truth_haps[inRegion2, 2] - af[inRegion2]) ** 2, 3)
        r2_3 <- round(cor(hap3[inRegion2] - af[inRegion2], truth_haps[inRegion2, 3] - af[inRegion2]) ** 2, 3) 
        print_message(paste0(prefix, "hap r2 mt:", r2_1, ", mu:", r2_2, ", pt:", r2_3, suffix))
    }
    return(to_return)
}






## final_phasing_accuracy_calculation <- function(
##     have_truth_haplotypes,
##     have_truth_genotypes,
##     truth_haps,
##     inRegion2,
##     dosage,
##     af,
##     phasing_haps,
##     gen,
##     sample_name
## ) {
##     if (have_truth_haplotypes) {    
##         w <- (inRegion2)
##         g <- truth_haps[inRegion2, 1] + truth_haps[inRegion2, 2]
##         r2 <-  round(cor((dosage)[inRegion2] - 2 * af[inRegion2], g - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
##         ##
##         x <- calculate_pse_and_r2_during_gibbs(inRegion2 = inRegion2, hap1 = phasing_haps[, 1], hap2 = phasing_haps[, 2], truth_haps = truth_haps, af = af, verbose = FALSE)
##         print_message(paste0("Final imputation dosage accuracy for sample ", sample_name, ", r2:", r2))
##         print_message(paste0("Final phasing accuracy for sample ", sample_name, ", pse:", x["pse"], ", disc(%):", x["disc"], "%"))
##     } else if (have_truth_genotypes) {
##         r2 <-  round(cor((dosage)[inRegion2] - 2 * af[inRegion2], gen[inRegion2, sample_name] - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
##         print_message(paste0("Final imputation dosage accuracy for sample ", sample_name, ", r2:", r2))
##     }
## }


## final_phasing_accuracy_calculation_rare_common <- function(
##     have_truth_haplotypes,
##     have_truth_genotypes,
##     truth_haps_all,
##     dosage_all,
##     hap1_all,
##     hap2_all,
##     gen_all,
##     sampleNames,
##     iSample,
##     special_rare_common_objects,
##     sample_name
## ) {
##     ref_alleleCount_all <- special_rare_common_objects[["ref_alleleCount_all"]]
##     af <- ref_alleleCount_all[, 3]
##     inRegion2 <- special_rare_common_objects[["inRegion2"]]
##     if (have_truth_haplotypes) {    
##         w <- (inRegion2)
##         g <- truth_haps_all[inRegion2, 1] + truth_haps_all[inRegion2, 2]
##         r2 <-  round(cor((dosage_all)[inRegion2] - 2 * af[inRegion2], g - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
##         ##
##         x <- calculate_pse_and_r2_during_gibbs(inRegion2 = inRegion2, hap1 = hap1_all, hap2 = hap2_all, truth_haps = truth_haps_all, af = af, verbose = FALSE)
##         print_message(paste0("Final imputation dosage accuracy for sample ", sample_name, ", r2:", r2))
##         print_message(paste0("Final phasing accuracy for sample ", sample_name, ", pse:", x["pse"], ", disc(%):", x["disc"], "%"))
##     } else if (have_truth_genotypes) {
##         r2 <-  round(cor((dosage_all)[inRegion2] - 2 * af[inRegion2], gen_all[inRegion2, sampleNames[iSample]] - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
##         print_message(paste0("Final imputation dosage accuracy for sample ", sample_name, ", r2:", r2))
##     }
## }
