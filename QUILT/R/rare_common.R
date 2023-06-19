make_eHapsCurrent_tc_using_rare_and_common_stuff <- function(
    hapMatcherR,
    distinctHapsIE,
    eMatDH_special_matrix_helper,
    eMatDH_special_matrix,
    rare_per_hap_info,
    snp_is_common,
    which_haps_to_use,
    snp_is_common_1_based,
    Ksubset,
    ref_error
) {
    one_minus_ref_error <- 1 - ref_error
    Ksubset <- length(which_haps_to_use)
    nSNPs <- length(snp_is_common)
    nCommonSNPs <- length(snp_is_common_1_based)
    nCommonGrids <- ceiling(nCommonSNPs / 32)
    ## should be easy in the first instance
    eHapsCurrent_tc <- array(ref_error, c(Ksubset, nSNPs, 1)) ## big! re-write to pass it around
    ## build eHapsCurrent_tc here
    for(iGrid in 1:nCommonGrids) {
        s <- 1 + 32 * (iGrid - 1)
        e <- min(32 * iGrid, nCommonSNPs)
        w <- snp_is_common_1_based[s:e]
        for(i_k in 1:Ksubset) {
            k <- which_haps_to_use[i_k]
            i <- as.integer(hapMatcherR[k, iGrid])
            if (i > 0) {
                ## note, could make distinctHapsB pretty easily if needed for RAM
                eHapsCurrent_tc[i_k, w, 1] <- distinctHapsIE[i, s:e]
            } else {
                b <- rcpp_simple_binary_matrix_search(
                    val = k - 1,
                    mat = eMatDH_special_matrix,
                    s1 = eMatDH_special_matrix_helper[iGrid, 1],
                    e1 = eMatDH_special_matrix_helper[iGrid, 2]
                )
                x <- rcpp_int_expand(b, nSNPs = e - s + 1)
                eHapsCurrent_tc[i_k, w[x == 1], 1] <- one_minus_ref_error
            }
        }
    }
    ## now to rare ones
    for(i_k in 1:length(which_haps_to_use)) {
        k <- which_haps_to_use[i_k]
        eHapsCurrent_tc[i_k, rare_per_hap_info[[k]], 1] <- one_minus_ref_error
    }
    eHapsCurrent_tc
}




get_initial_read_labels <- function(
    pos_all,
    allSNP_sampleReads,
    snp_is_common,
    hap1,
    hap2,
    maxDifferenceBetweenReads
) {
    nReads <- length(allSNP_sampleReads)
    eHapsCurrent_tc <- array(0.5, c(2, nrow(pos_all), 1))
    eHapsCurrent_tc[1, snp_is_common, 1] <- hap1
    eHapsCurrent_tc[2, snp_is_common, 1] <- hap2
    eMatRead_t <- array(1, c(2, nReads))
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = allSNP_sampleReads,
        eHapsCurrent_tc = eHapsCurrent_tc,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        s = 0,
        Jmax = 100,
        eMatHapOri_t = array(0, c(1, 1)),
        pRgivenH1 = array(0, 1),
        pRgivenH2 = array(0, 1),
        prev = 0,
        suppressOutput = 1,
        prev_section = "wer",
        next_section = "wer",
        rescale_eMatRead_t = TRUE,
        run_pseudo_haploid = FALSE
    )
    H <- as.integer(runif(nReads) < (eMatRead_t[1, ] / colSums(eMatRead_t))) + 1
    H
}

impute_final_gibbs_with_rare_common <- function(
    special_rare_common_objects,
    special_rare_common_objects_per_core,                                                
    allSNP_sampleReads,
    hap1,
    hap2,
    pos_all,
    maxDifferenceBetweenReads,
    hapMatcherR,
    distinctHapsIE,
    eMatDH_special_matrix_helper,
    eMatDH_special_matrix,
    Ksubset,
    ref_error,
    which_haps_to_use,
    small_ref_panel_gibbs_iterations,
    small_ref_panel_block_gibbs_iterations,
    allSNP_wif0,
    allSNP_grid_has_read,
    make_plots,
    outplotprefix,
    have_truth_haplotypes,
    truth_haps_all,
    have_truth_genotypes,
    truth_gen_all,
    truth_labels,
    uncertain_truth_labels,
    shuffle_bin_radius,
    make_plots_block_gibbs,
    sample_name,
    regionStart,
    regionEnd,
    buffer,
    i_it,
    i_gibbs_sample,
    ff0_shard_check_every_pair,
    sampleNames,
    iSample,
    phase_all
) {

    snp_is_common <- special_rare_common_objects[["snp_is_common"]]
    rare_per_hap_info <- special_rare_common_objects[["rare_per_hap_info"]]
    nSNPs <- nrow(pos_all)
    
    H <- get_initial_read_labels(
        pos_all = pos_all,
        allSNP_sampleReads = allSNP_sampleReads,
        snp_is_common = snp_is_common,
        hap1 = hap1,
        hap2 = hap2,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads
    )   

    snp_is_common_1_based <- which(snp_is_common)
    small_eHapsCurrent_tc <- make_eHapsCurrent_tc_using_rare_and_common_stuff(
        hapMatcherR = hapMatcherR,
        distinctHapsIE = distinctHapsIE,
        eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
        eMatDH_special_matrix = eMatDH_special_matrix,
        rare_per_hap_info = rare_per_hap_info,
        snp_is_common = snp_is_common,
        which_haps_to_use = which_haps_to_use,
        snp_is_common_1_based = snp_is_common_1_based,
        Ksubset = Ksubset,
        ref_error = ref_error
    )
    

    if (have_truth_haplotypes) {
        ## match
        s <- sampleNames[iSample]
        truth_haps <- cbind(phase_all[, s, 1], phase_all[, s, 2])
    } else {
        truth_haps <- NULL
    }

    
    use_small_eHapsCurrent_tc <- TRUE
    use_provided_small_eHapsCurrent_tc <- TRUE

    ## get objects out for convenience
    alphaHat_t1 <- special_rare_common_objects_per_core[["alphaHat_t1"]]
    betaHat_t1 <- special_rare_common_objects_per_core[["betaHat_t1"]]
    eMatGrid_t1 <- special_rare_common_objects_per_core[["eMatGrid_t1"]]
    alphaHat_t2 <- special_rare_common_objects_per_core[["alphaHat_t2"]]
    betaHat_t2 <- special_rare_common_objects_per_core[["betaHat_t2"]]
    eMatGrid_t2 <- special_rare_common_objects_per_core[["eMatGrid_t2"]]
    alphaHat_t3 <- special_rare_common_objects_per_core[["alphaHat_t3"]]
    betaHat_t3 <- special_rare_common_objects_per_core[["betaHat_t3"]]
    eMatGrid_t3 <- special_rare_common_objects_per_core[["eMatGrid_t3"]]
    gammaMT_t_local <- special_rare_common_objects_per_core[["gammaMT_t_local"]]
    gammaMU_t_local <- special_rare_common_objects_per_core[["gammaMU_t_local"]]
    gammaP_t_local <- special_rare_common_objects_per_core[["gammaP_t_local"]]

    ## global ones
    smooth_cm <- special_rare_common_objects[["smooth_cm"]]
    small_transMatRate_tc_H <- special_rare_common_objects[["small_transMatRate_tc_H"]]
    small_alphaMatCurrent_tc <- special_rare_common_objects[["small_alphaMatCurrent_tc"]]
    small_priorCurrent_m <- special_rare_common_objects[["small_priorCurrent_m"]]
    ref_alleleCount_all <- special_rare_common_objects[["ref_alleleCount_all"]]
    ancAlleleFreqAll <- ref_alleleCount_all[, 3]
    grid <- special_rare_common_objects[["grid"]]
    L_grid <- special_rare_common_objects[["L_grid"]]
    L <- special_rare_common_objects[["L"]]
    inRegion2 <- special_rare_common_objects[["inRegion2"]]    

    ## 
    double_list_of_starting_read_labels <- list(list(H))

    gibbs_iterate <- impute_one_sample(
        nSNPs = nSNPs,
        sampleReads = allSNP_sampleReads,
        small_eHapsCurrent_tc = small_eHapsCurrent_tc,
        small_transMatRate_tc_H = small_transMatRate_tc_H,
        alphaHat_t1 = alphaHat_t1,
        betaHat_t1 = betaHat_t1,
        eMatGrid_t1 = eMatGrid_t1,
        alphaHat_t2 = alphaHat_t2,
        betaHat_t2 = betaHat_t2,
        eMatGrid_t2 = eMatGrid_t2,
        alphaHat_t3 = alphaHat_t3,
        betaHat_t3 = betaHat_t3,
        eMatGrid_t3 = eMatGrid_t3,
        gammaMT_t_local = gammaMT_t_local,
        gammaMU_t_local = gammaMU_t_local,
        gammaP_t_local = gammaP_t_local,
        small_alphaMatCurrent_tc = small_alphaMatCurrent_tc,
        small_priorCurrent_m = small_priorCurrent_m,
        smooth_cm = smooth_cm,
        which_haps_to_use = which_haps_to_use,
        small_ref_panel_gibbs_iterations = small_ref_panel_gibbs_iterations,
        double_list_of_starting_read_labels = double_list_of_starting_read_labels,
        small_ref_panel_block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
        wif0 = allSNP_wif0,
        grid_has_read = allSNP_grid_has_read,
        make_plots = make_plots,
        plot_description = paste0("it", i_it, ".full.gibbs"),
        ancAlleleFreqAll = ancAlleleFreqAll,
        grid = grid,
        L_grid = L_grid,
        L = L,
        inRegion2 = inRegion2,
        cM_grid = cM_grid,
        outplotprefix = outplotprefix,
        have_truth_haplotypes = have_truth_haplotypes,
        truth_haps = truth_haps_all,
        have_truth_genotypes = have_truth_genotypes,
        truth_gen = truth_gen_all,
        truth_labels = truth_labels,
        uncertain_truth_labels = uncertain_truth_labels,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        shuffle_bin_radius = shuffle_bin_radius,
        make_plots_block_gibbs = make_plots_block_gibbs,
        sample_name = sample_name,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        use_small_eHapsCurrent_tc = use_small_eHapsCurrent_tc,
        use_provided_small_eHapsCurrent_tc = use_provided_small_eHapsCurrent_tc,
        use_sample_is_diploid = TRUE,
        i_it = i_it,
        i_gibbs_sample = i_gibbs_sample,
        ff0_shard_check_every_pair = ff0_shard_check_every_pair,
        suppressOutput = 1
    )
    
    ## AM HERE
    hapProbs_t <- gibbs_iterate$hapProbs_t
    hap1 <- hapProbs_t[1, ]
    hap2 <- hapProbs_t[2, ]

    return(
        list(
            hap1 = hap1,
            hap2 = hap2
        )
    )

}



calculate_pse_and_r2_rare_common <- function(
    hap1_all,
    hap2_all,
    have_truth_haplotypes,
    truth_haps_all,
    have_truth_genotypes,
    truth_gen_all,
    special_rare_common_objects,
    verbose
) {
    
    ref_alleleCount_all <- special_rare_common_objects[["ref_alleleCount_all"]]
    af <- ref_alleleCount_all[, 3]
    inRegion2 <- special_rare_common_objects[["inRegion2"]]
    
    if (have_truth_haplotypes) {
        x <- calculate_pse_and_r2_during_gibbs(
            inRegion2 = inRegion2,
            hap1 = hap1_all,
            hap2 = hap2_all,
            truth_haps = truth_haps_all,
            af = af,
            verbose = verbose,
            impute_rare_common = TRUE,
            all_snps = TRUE
        )
    } else if (have_truth_genotypes) {
        r2 <-  round(cor((dosage_all)[inRegion2] - 2 * af[inRegion2], truth_gen_all[inRegion2, ] - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
        print_message(paste0("Current accuracy for this gibbs sample for ", sample_name, " for all SNPs, r2:", r2))
    }
    
}



final_phasing_accuracy_calculation_rare_common <- function(
    have_truth_haplotypes,
    have_truth_genotypes,
    truth_haps_all,
    dosage_all,
    hap1_all,
    hap2_all,
    gen_all,
    sampleNames,
    iSample,
    special_rare_common_objects,
    sample_name
) {
    ref_alleleCount_all <- special_rare_common_objects[["ref_alleleCount_all"]]
    af <- ref_alleleCount_all[, 3]
    inRegion2 <- special_rare_common_objects[["inRegion2"]]
    if (have_truth_haplotypes) {    
        w <- (inRegion2)
        g <- truth_haps_all[inRegion2, 1] + truth_haps_all[inRegion2, 2]
        r2 <-  round(cor((dosage_all)[inRegion2] - 2 * af[inRegion2], g - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
        ##
        x <- calculate_pse_and_r2_during_gibbs(inRegion2 = inRegion2, hap1 = hap1_all, hap2 = hap2_all, truth_haps = truth_haps_all, af = af, verbose = FALSE)
        print_message(paste0("Final imputation dosage accuracy for sample ", sample_name, ", r2:", r2))
        print_message(paste0("Final phasing accuracy for sample ", sample_name, ", pse:", x["pse"], ", disc(%):", x["disc"], "%"))
    } else if (have_truth_genotypes) {
        r2 <-  round(cor((dosage_all)[inRegion2] - 2 * af[inRegion2], gen_all[inRegion2, sampleNames[iSample]] - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
        print_message(paste0("Final imputation dosage accuracy for sample ", sample_name, ", r2:", r2))
    }
}
