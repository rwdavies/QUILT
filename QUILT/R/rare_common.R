make_eHapsCurrent_tc_using_rare_and_common_stuff <- function(
    hapMatcher,
    hapMatcherR,
    use_hapMatcherR,
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
            if (use_hapMatcherR) {
                i <- as.integer(hapMatcherR[k, iGrid])
            } else {
                i <- as.integer(hapMatcher[k, iGrid])
            }
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
    hap3 = NULL,
    maxDifferenceBetweenReads,
    ff = NULL
) {
    nReads <- length(allSNP_sampleReads)
    nhap <- 2
    if (!is.null(hap3)) {
        nhap <- 3
    }
    eHapsCurrent_tc <- array(0.5, c(nhap, nrow(pos_all), 1))
    eHapsCurrent_tc[1, snp_is_common, 1] <- hap1
    eHapsCurrent_tc[2, snp_is_common, 1] <- hap2
    if (nhap == 3) {
        eHapsCurrent_tc[3, snp_is_common, 1] <- hap3
    }
    eMatRead_t <- array(1, c(nhap, nReads))
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
    if (nhap == 2) {
        H <- as.integer(runif(nReads) < (eMatRead_t[1, ] / colSums(eMatRead_t))) + 1
    } else {
        out <- get_read_groupings_given_fetal_fraction_and_cov(readProbs_t = eMatRead_t, phase = NULL, iiSample = NA, ff = ff)
        H <- sample_H_for_NIPT_given_groupings(groupings = out$groupings, counts = out$counts, ff = ff)
    }
    H
}

impute_final_gibbs_with_rare_common <- function(
    special_rare_common_objects,
    special_rare_common_objects_per_core,                                                
    allSNP_sampleReads,
    hap1,
    hap2,
    hap3,
    pos_all,
    maxDifferenceBetweenReads,
    hapMatcher,
    hapMatcherR,
    use_hapMatcherR,
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
    truth_labels_all,
    uncertain_truth_labels_all,
    shuffle_bin_radius,
    make_plots_block_gibbs,
    sample_name,
    regionStart,
    regionEnd,
    buffer,
    i_it,
    i_gibbs_sample,
    shard_check_every_pair,
    sampleNames,
    iSample,
    phase_all,
    ff,
    method
) {


    ## save(
    ## special_rare_common_objects,
    ## special_rare_common_objects_per_core,                                                
    ## allSNP_sampleReads,
    ## hap1,
    ## hap2,
    ## pos_all,
    ## maxDifferenceBetweenReads,
    ## hapMatcher,
    ## hapMatcherR,
    ## use_hapMatcherR,
    ## distinctHapsIE,
    ## eMatDH_special_matrix_helper,
    ## eMatDH_special_matrix,
    ## Ksubset,
    ## ref_error,
    ## which_haps_to_use,
    ## small_ref_panel_gibbs_iterations,
    ## small_ref_panel_block_gibbs_iterations,
    ## allSNP_wif0,
    ## allSNP_grid_has_read,
    ## make_plots,
    ## outplotprefix,
    ## have_truth_haplotypes,
    ## truth_haps_all,
    ## have_truth_genotypes,
    ## truth_gen_all,
    ## truth_labels_all,
    ## uncertain_truth_labels_all,
    ## shuffle_bin_radius,
    ## make_plots_block_gibbs,
    ## sample_name,
    ## regionStart,
    ## regionEnd,
    ## buffer,
    ## i_it,
    ## i_gibbs_sample,
    ## shard_check_every_pair,
    ## sampleNames,
    ## iSample,
    ## phase_all,
    ## method,
    ## file = "/data/smew1/rdavies/temp/123.RData", compress = FALSE)
    ## stop("WER")


    snp_is_common <- special_rare_common_objects[["snp_is_common"]]
    rare_per_hap_info <- special_rare_common_objects[["rare_per_hap_info"]]
    nSNPs <- nrow(pos_all)

    if (method == "nipt") { 
        H <- get_initial_read_labels(
            pos_all = pos_all,
            allSNP_sampleReads = allSNP_sampleReads,
            snp_is_common = snp_is_common,
            hap1 = hap1,
            hap2 = hap2,
            hap3 = hap3,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            ff = ff
        )
    } else {
        H <- get_initial_read_labels(
            pos_all = pos_all,
            allSNP_sampleReads = allSNP_sampleReads,
            snp_is_common = snp_is_common,
            hap1 = hap1,
            hap2 = hap2,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads
        )
    }
    
    ## snp_is_common_1_based <- which(snp_is_common)
    common_snp_index <- integer(nSNPs)
    common_snp_index[which(snp_is_common)] <- 1L:as.integer(sum(snp_is_common))
    
    ## small_eHapsCurrent_tc <- make_eHapsCurrent_tc_using_rare_and_common_stuff(
    ##     hapMatcher = hapMatcher,        
    ##     hapMatcherR = hapMatcherR,
    ##     use_hapMatcherR = use_hapMatcherR,
    ##     distinctHapsIE = distinctHapsIE,
    ##     eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
    ##     eMatDH_special_matrix = eMatDH_special_matrix,
    ##     rare_per_hap_info = rare_per_hap_info,
    ##     snp_is_common = snp_is_common,
    ##     which_haps_to_use = which_haps_to_use,
    ##     snp_is_common_1_based = snp_is_common_1_based,
    ##     Ksubset = Ksubset,
    ##     ref_error = ref_error
    ## )
    
    ## do this to get quick indicator on whether we need to examine this k and this read for rare haps
    ## not perfect overall but should be reasonably fast, and use much less RAM than making small_eHapsCurrent_tc

    ## do not need this anymore
    ## eMatRead_t <- determine_which_reads_require_k_haps(
    ##     sampleReads = allSNP_sampleReads,
    ##     nSNPs = nSNPs,
    ##     rare_per_hap_info = rare_per_hap_info,
    ##     which_haps_to_use = which_haps_to_use
    ## )
    Ksubset <- length(which_haps_to_use)
    nReads <- length(allSNP_sampleReads)
    ## 1 is normal
    ## 0 is special
    eMatRead_t <- array(1, c(Ksubset, nReads))
    
    

    if (have_truth_haplotypes) {
        ## match
        s <- sampleNames[iSample]
        truth_haps <- cbind(phase_all[, s, 1], phase_all[, s, 2])
    } else {
        truth_haps <- NULL
    }


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




##print("test it out")
##    Rcpp_make_eMatRead_t_for_final_rare_common_gibbs_using_objects(eMatRead_t, rare_per_hap_info, common_snp_index, snp_is_common, allSNP_sampleReads,hapMatcherR,grid,distinctHapsIE,eMatDH_special_matrix_helper,eMatDH_special_matrix,ref_error,which_haps_to_use,rescale_eMatRead_t = TRUE,Jmax = 100,maxDifferenceBetweenReads);
##print("end of test it out")    


    ## hopefully not a RAM monster
    ## create a new one that is per-site, and gives list of k
    rare_per_snp_info <- lapply(1:nSNPs, function(x) -1L)
    ## want this to be the index of the rare SNP in the overall
    for(k in 1:length(which_haps_to_use)) {
        ## expand back out to all SNPs
        snps <- rare_per_hap_info[[which_haps_to_use[[k]]]] ## this IS an index among all SNPs
        ## this is among all SNPs
        for(snp in snps) {
            rare_per_snp_info[[snp]] <- c(rare_per_snp_info[[snp]], k) ## keep 1-based
        }
    }
    

    gibbs_iterate <- impute_one_sample(
        nSNPs = nSNPs,
        hapMatcherR = hapMatcherR,
        sampleReads = allSNP_sampleReads,
        distinctHapsIE = distinctHapsIE,
        eMatDH_special_matrix = eMatDH_special_matrix,
        eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
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
        outplotprefix = outplotprefix,
        have_truth_haplotypes = have_truth_haplotypes,
        truth_haps = truth_haps_all,
        have_truth_genotypes = have_truth_genotypes,
        truth_gen = truth_gen_all,
        truth_labels = truth_labels_all,
        uncertain_truth_labels = uncertain_truth_labels_all,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        shuffle_bin_radius = shuffle_bin_radius,
        make_plots_block_gibbs = make_plots_block_gibbs,
        sample_name = sample_name,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        use_small_eHapsCurrent_tc = FALSE, ## set this to FALSE now
        use_sample_is_diploid = TRUE,
        i_it = i_it,
        i_gibbs_sample = i_gibbs_sample,
        shard_check_every_pair = shard_check_every_pair,
        suppressOutput = 1,
        verbose = FALSE,
        return_extra = FALSE,
        eMatRead_t = eMatRead_t,
        make_eMatRead_t_rare_common = TRUE, ## need this to use eMatRead_t properly
        common_snp_index = common_snp_index,
        snp_is_common = snp_is_common,
        rare_per_hap_info = rare_per_hap_info,
        rare_per_snp_info = rare_per_snp_info,
        ff = ff,
        method = method
    )

    

    hapProbs_t <- gibbs_iterate$hapProbs_t
    hap1 <- hapProbs_t[1, ]
    hap2 <- hapProbs_t[2, ]
    hap3 <- NULL
    if (method == "nipt") {
        hap3 <- hapProbs_t[3, ]
    }

    return(
        list(
            hap1 = hap1,
            hap2 = hap2,
            hap3 = hap3
        )
    )

}









determine_which_reads_require_k_haps <- function(
    sampleReads,
    nSNPs,
    rare_per_hap_info,
    which_haps_to_use
) {
    ##
    Ksubset <- length(which_haps_to_use)
    nReads <- length(sampleReads)
    ## 1 is normal
    ## 0 is special
    eMatRead_t <- array(1, c(Ksubset, nReads))
    ## rare_per_hap_info is 1-based
    ## want knowledge of for each hap, 
    ## build an index of what sites have a k
    has_a_k <- logical(nSNPs)
    for(k in 1:Ksubset) {
        x <- rare_per_hap_info[[which_haps_to_use[[k]]]]
        for(i in 1:length(x)) {
            has_a_k[x[i]] <- TRUE
        }
    }
    n_affected_SNPs <- sum(has_a_k)
    ## now get their positions
    snp_n_for_has_a_k <- integer(nSNPs)
    snp_n_for_has_a_k[has_a_k] <- 1:n_affected_SNPs
    ## among those sites, figure out which k are involved
    which_k <- vector("list", n_affected_SNPs)
    for(k in 1:Ksubset) {
        x <- rare_per_hap_info[[which_haps_to_use[[k]]]]
        for(n in snp_n_for_has_a_k[x]) {        
            which_k[[n]] <- c(which_k[[n]], k)
        }
    }
    ## now find reads
    read_has_a_k <- logical(nReads)
    for(iRead in 1:nReads) {
        u <- sampleReads[[iRead]][[4]] + 1 ## 1-based
        if (sum(has_a_k[u]) > 0) {
            ## get these (could be more than 1)
            x <- snp_n_for_has_a_k[u]
            y <- x[x > 0]
            for(yy in y) {
                ks <- which_k[[yy]]
                eMatRead_t[ks, iRead] <- 0
            }
        }
    }
    eMatRead_t
}
