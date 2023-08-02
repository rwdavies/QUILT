## zilong mspbwt32.cpp
## @param msp is xptr to msPBWT object in C++
select_new_haps_zilong_msp <- function(
    hapProbs_t,
    igibbs,
    outputdir,
    Kfull,
    Knew,
    msp,
    mspbwtB,
    mspbwtL,
    mspbwtM,
    nGrids
) {
  ##
  out <- lapply(1:2, function(x) {
    hap <- round(hapProbs_t[x, ])
    res <- as.data.frame(mspbwt_report(msp, mspbwtB, hap, mspbwtL, mspbwtM))
    res <- res[order(-res$lens),]
    res <- res[!duplicated(res[,c("haps")]),]
    res$starts <- res$ends - res$lens + 1
    res$keys <- res$starts * nGrids + res$ends
    res$haps <- res$haps + 1
    res
  })
  ## res <- do.call(rbind.data.frame, res)
  unique_haps <- unique(c(out[[1]][, "haps"], out[[2]][, "haps"]))
  ## interhaps <- intersect(res[res$ihap == 1, "haps"], res[res$ihap == 2, "haps"])
  ## saveRDS(res, file = file.path(outputdir, paste0("which_haps_to_use.i", igibbs, ".zilong.rds")))
  print_message(paste("select", length(unique_haps), " unique haps by mpbwt query before post-selection"))
  ## return(unique_haps[1:Knew])
  if (length(unique_haps) == 0) {
    new_haps <- sample(1:Kfull, Knew)
    return(new_haps)
  } else if (length(unique_haps) <= Knew) {
    ## cannot take a sample larger than the population when 'replace = FALSE'
    new_haps <- unique(c(
      unique_haps,
      sample(Kfull, min(length(unique_haps) + Knew, Knew), replace = FALSE)
    ))[1:Knew]
    print_message(paste("select", length(unique(new_haps)), " unique haps after post-selection 1"))
    return(new_haps)
  } else {
    ##
    ## heuristically, prioritize based on length and new-ness
    ## do this for each of the two haps
    ##
    results <- lapply(out, function(mtm) {
      m <- max(mtm[, "ends"])
      weight <- numeric(m)
      cur_sum <- numeric(m)
      cur_sum[] <- 1
      for(i in 1:nrow(mtm)) {
        s <- mtm[i, "starts"]
        e <- mtm[i, "ends"]
        weight[i] <- (e - s + 1) * 1 / sum(cur_sum[s:e])
        cur_sum[s:e] <- cur_sum[s:e] + 1
      }
      ##
      o <- order(-weight)
      mtm <- mtm[o, ]
      mtm[, "haps"]
    })
    ## pad out one of them
    x <- results[[1]]
    y <- results[[2]]
    if (length(x) > length(y)) {
      y <- c(y, rep(NA, length(x) - length(y)))
    } else {
      x <- c(x, rep(NA, length(y) - length(x)))
    }
    unique_ordered_haps <- unique(c(t(cbind(x, y))))
    unique_ordered_haps <- unique_ordered_haps[!is.na(unique_ordered_haps)]
    print_message(paste("select", length(unique(unique_ordered_haps)), " unique haps after post-selection 2"))
    if (length(unique_ordered_haps) >= Knew) {
      return(unique_ordered_haps[1:Knew])
    } else {
      ## add in some other (potentially) duplicated haps
      new_haps <- c(setdiff(unique_ordered_haps, unique_haps), unique_haps)[1:Knew]
      print(paste("select", length(unique(new_haps)), " unique haps after post-selection 3"))
      return(new_haps)
    }
  }
}


## zilong mspbwt32.cpp
## @param msp is xptr to msPBWT object in C++
select_new_haps_zilong_msp_robbie_version <- function(
    hapProbs_t,
    igibbs,
    outputdir,
    Kfull,
    Knew,
    msp,
    mspbwtB,
    mspbwtL,
    mspbwtM,
    nGrids,
    aggregated = FALSE
) {
    res <- lapply(1:2, function(x) {
        hap <- round(hapProbs_t[x, ])
        res <- mspbwt_report(msp, hap, mspbwtL, mspbwtB, aggregated = aggregated)
        return(res)
    })
    ## merge, filter, order
    ends <- c(res[[1]][["ends"]], res[[2]][["ends"]])
    lens <- c(res[[1]][["lens"]], res[[2]][["lens"]])
    res <- cbind(
        haps = c(res[[1]][["haps"]], res[[2]][["haps"]]),
        starts = ends - lens + 1,
        ends = ends,
        lens = lens,
        n = c(res[[1]][["n"]], res[[2]][["n"]])
    )
    res <- res[res[, "lens"] >= max(mspbwtM, 0), , drop = FALSE]
    res <- cbind(res, keys = res[, "starts"] * nGrids + res[, "ends"])
    ## remove those that are the same thing (effectively)
    res <- res[!duplicated((res[, "n"] + 1) * (res[, "starts"] + 1) + nGrids * (res[, "haps"])), ]
    ## remove those that are the same key
    res <- res[!duplicated(res[, "keys"]), ]
    ## can make this linear rather than an order
    ## I'm too lazy to do it now, see if it comes back with a profile
    ## print("make me linear")
    res <- res[order(-res[, "lens"]), ]
    ## 
    unique_haps <- unique(res[, "haps"])
    ## 
    if (length(unique_haps) == 0) {
        new_haps <- sample(1:Kfull, Knew)
        return(new_haps)
    } else if (length(unique_haps) <= Knew) {
        ## cannot take a sample larger than the population when 'replace = FALSE'
        new_haps <- unique(c(
            unique_haps,
            sample(Kfull, min(length(unique_haps) + Knew, Knew), replace = FALSE)
        ))[1:Knew]
        print(paste("select", length(unique(new_haps)), " unique haps after post-selection 1"))
        return(new_haps)
    } else {
        ##
        ## heuristically, prioritize based on length and new-ness
        ## 
        m <- max(res[, "ends"])
        weight <- numeric(m)
        cur_sum <- numeric(m)
        cur_sum[] <- 1
        for(i in 1:nrow(res)) {
            s <- res[i, "starts"]
            e <- res[i, "ends"]
            weight[i] <- (e - s + 1) * 1 / sum(cur_sum[s:e])
            cur_sum[s:e] <- cur_sum[s:e] + 1
        }
        ##
        unique_ordered_haps <- unique(res[order(-weight), "haps"])
        print(paste("select", length(unique_ordered_haps), " unique haps after post-selection 2"))
        ## 
        if (length(unique_ordered_haps) >= Knew) {
            return(unique_ordered_haps[1:Knew])
        } else {
            ## add in some other (potentially) duplicated haps
            new_haps <- c(setdiff(unique_ordered_haps, unique_haps), unique_haps)[1:Knew]
            print(paste("select", length(unique(new_haps)), " unique haps after post-selection 3"))
            return(new_haps)
        }
  }
}





    

get_and_impute_one_sample <- function(
    rhb_t,
    outputdir,
    nGibbsSamples,
    n_seek_its,
    n_burn_in_seek_its,
    full_alphaHat_t,
    full_betaHat_t,
    full_gamma_t,
    full_gammaSmall_t,
    full_gammaSmall_cols_to_get,
    full_transMatRate_t_H,
    small_transMatRate_tc_H,
    alphaHat_t1,
    betaHat_t1,
    eMatGrid_t1,
    alphaHat_t2,
    betaHat_t2,
    eMatGrid_t2,
    alphaHat_t3,
    betaHat_t3,
    eMatGrid_t3,
    gammaMT_t_local,
    gammaMU_t_local,
    gammaP_t_local,
    small_alphaMatCurrent_tc,
    small_priorCurrent_m,
    small_eHapsCurrent_tc,
    bam_files,
    L,
    pos,
    chr,
    tempdir,
    regionName,
    regionStart,
    regionEnd,
    buffer,
    gen,
    phase,
    gen_all,
    phase_all,
    iSample,
    grid,
    ancAlleleFreqAll,
    L_grid,
    verbose,
    shuffle_bin_radius,
    Ksubset,
    Knew,
    K_top_matches,
    heuristic_match_thin,
    record_interim_dosages,
    have_truth_haplotypes,
    have_truth_genotypes, 
    bqFilter,
    record_read_label_usage,
    sampleNames,
    smooth_cm,
    iSizeUpperLimit,
    maxDifferenceBetweenReads,
    make_plots,
    ref_error,
    distinctHapsB,
    distinctHapsIE,
    eMatDH_special_matrix_helper,
    eMatDH_special_matrix,
    use_eMatDH_special_symbols,
    hapMatcher,
    hapMatcherR,
    use_hapMatcherR,
    eMatDH_special_grid_which,
    eMatDH_special_values_list,
    inRegion2,
    cM_grid,
    af,
    use_bx_tag,
    bxTagUpperLimit,
    addOptimalHapsToVCF,
    make_plots_block_gibbs,
    estimate_bq_using_truth_read_labels,
    chrStart,
    chrEnd,
    gamma_physically_closest_to,
    hla_run,
    downsampleToCov,
    minGLValue,
    minimum_number_of_sample_reads,
    print_extra_timing_information,
    small_ref_panel_gibbs_iterations,
    small_ref_panel_block_gibbs_iterations,
    plot_per_sample_likelihoods,
    use_small_eHapsCurrent_tc,
    output_gt_phased_genotypes,
    mspbwtB,
    mspbwtL,
    mspbwtM,
    zilong,
    msp,
    use_mspbwt,
    ms_indices,
    use_splitreadgl,
    use_sample_is_diploid,
    plot_p1,
    small_ref_panel_skip_equally_likely_reads,
    small_ref_panel_equally_likely_reads_update_iterations,
    ff0_shard_check_every_pair,
    use_eigen,
    pos_all,
    special_rare_common_objects,
    special_rare_common_objects_per_core,    
    impute_rare_common,
    make_heuristic_plot,
    heuristic_approach
) {

    sample_name <- sampleNames[iSample]
    nSNPs <- nrow(pos)
    nGrids <- ncol(distinctHapsB)
    K <- nrow(hapMatcher)
    suppressOutput <- !print_extra_timing_information

    if (impute_rare_common) {

        loadBamAndConvert(
            iBam = iSample,
            L = pos_all[, 2],
            pos = pos_all,
            nSNPs = nrow(pos_all),
            bam_files = bam_files,
            iSizeUpperLimit = iSizeUpperLimit,
            bqFilter = bqFilter,
            chr = chr,
            N = length(sampleNames),
            downsampleToCov = downsampleToCov,
            sampleNames = sampleNames,
            inputdir = tempdir,
            regionName = regionName,
            tempdir = tempdir,
            chrStart = chrStart,
            chrEnd = chrEnd,
            chrLength = NA,
            save_sampleReadsInfo = TRUE,
            use_bx_tag = use_bx_tag,
            bxTagUpperLimit = bxTagUpperLimit,
            default_sample_no_read_behaviour = "return_null"
        )

        load(file_sampleReads(tempdir, iSample, regionName))
        load(file_sampleReadsInfo(tempdir, iSample, regionName))
        allSNP_sampleReads <- sampleReads
        allSNP_sampleReadsInfo <- sampleReadsInfo
        rm(sampleReads, sampleReadsInfo)

        allSNP_sampleReads <- snap_sampleReads_to_grid(
            sampleReads = allSNP_sampleReads,
            grid = special_rare_common_objects[["grid"]]
        )

        nGrids_all <- length(special_rare_common_objects$L_grid)
        allSNP_wif0 <- as.integer(sapply(allSNP_sampleReads, function(x) x[[2]]))
        allSNP_grid_has_read <- rep(FALSE, nGrids_all)
        allSNP_grid_has_read[unique(allSNP_wif0) + 1] <- TRUE

        if (have_truth_haplotypes) {
            s <- sampleNames[iSample]
            if (!(s %in% dimnames(phase_all)[[2]])) {
                stop("Something went wrong with phase naming")
            }
            truth_haps_all <- cbind(phase_all[, s, 1], phase_all[, s, 2])
        } else {
            truth_haps_all <- NULL
        }

        if (have_truth_genotypes) {
            s <- sampleNames[iSample]
            if (!(s %in% colnames(gen))) {
                stop("Something went wrong with gen naming")
            }
            truth_gen_all <- gen_all[, s, drop = FALSE]
        } else {
            truth_gen_all <- NULL
        }

        if (!is.null(truth_haps_all)) {
            truth_label_set <- determine_a_set_of_truth_labels(
                sampleReads = allSNP_sampleReads,
                truth_hap1 = truth_haps_all[, 1],
                truth_hap2 = truth_haps_all[, 2],
                maxDifferenceBetweenReads = maxDifferenceBetweenReads
            )
            truth_labels_all <- truth_label_set[["truth_labels"]]
            uncertain_truth_labels_all <- truth_label_set[["uncertain_truth_labels"]]
            rm(truth_label_set)
        } else {
            truth_labels_all <- NULL
            uncertain_truth_labels_all <- NULL
        }
        
        nSNPs_all <- nrow(pos_all)        
        dosage_all <- numeric(nSNPs_all)
        gp_t_all <- array(0, c(3, nSNPs_all))


    }

    ##
    ## sample read stuff - work off bam file!
    ##
    loadBamAndConvert(
        iBam = iSample,
        L = L,
        pos = pos,
        nSNPs = nSNPs,
        bam_files = bam_files,
        iSizeUpperLimit = iSizeUpperLimit,
        bqFilter = bqFilter,
        chr = chr,
        N = length(sampleNames),
        downsampleToCov = downsampleToCov,
        sampleNames = sampleNames,
        inputdir = tempdir,
        regionName = regionName,
        tempdir = tempdir,
        chrStart = chrStart,
        chrEnd = chrEnd,
        chrLength = NA,
        save_sampleReadsInfo = TRUE,
        use_bx_tag = use_bx_tag,
        bxTagUpperLimit = bxTagUpperLimit,
        default_sample_no_read_behaviour = "return_null"
    )

    load(file_sampleReads(tempdir, iSample, regionName))
    load(file_sampleReadsInfo(tempdir, iSample, regionName))

    if (length(sampleReads) < minimum_number_of_sample_reads) {
        print_message(paste0("Sample number ", iSample, " with sample name ", sampleNames[iSample], " has ", length(sampleReads), " reads which is fewer than the minimum ", minimum_number_of_sample_reads, ". This sample will therefore not be imputed and all results will be set to missing"))
        ## dummy up output here?
        ## note - useful, see writers.R,
        ## FORMAT <- "GT:GP:DS:HD" ## genotypes, posteriors, dosages, haploid-dosages
        per_sample_vcf_col <- "./.:.,.,.:.:.,."
        if (addOptimalHapsToVCF & have_truth_haplotypes) {
            FORMAT <- "GT:GP:DS:HD:OHD" ## (phased) genotypes, genotype posteriors, dosages, haploid-dosages, optimal haploid-dosages
            per_sample_vcf_col <- paste0(per_sample_vcf_col, ":.,.")
        }
        to_return <- list(
            sample_was_imputed = FALSE,
            per_sample_vcf_col = per_sample_vcf_col
        )
        return(to_return)
    }

    sample_alleleCount <- get_alleleCount(sampleReads, nrow(pos))
    print_message(paste0("The average depth of this sample is:", mean(sample_alleleCount[, 2])))
    print_message(paste0("There are ", length(sampleReads), " reads under consideration"))

    sampleReads <- snap_sampleReads_to_grid(
        sampleReads = sampleReads,
        grid = grid
    )

    i_gibbs_sample <- 1
    dosage <- numeric(nSNPs)
    nDosage <- 0
    gp_t <- array(0, c(3, nSNPs))
    wif0 <- as.integer(sapply(sampleReads, function(x) x[[2]]))
    grid_has_read <- rep(FALSE, nGrids)
    grid_has_read[wif0 + 1] <- TRUE

    nReads <- length(sampleReads)
    super_out_hap_dosages <- as.list(1:nGibbsSamples)
    super_out_read_labels <- as.list(1:nGibbsSamples)
    super_out_dosage_matrix <- as.list(1:nGibbsSamples)
    read_label_matrix_all <- array(NA, c(nReads, nGibbsSamples)) ## need this for phasing - store final read labels in it
    read_label_matrix_conf <- array(FALSE, c(nReads, nGibbsSamples)) ## need this for phasing - store whether these are confident reads

    return_gamma_t <- FALSE ## can turn on later, if required

    phasing_read_labels <- NULL
    phasing_haps <- NULL
    phasing_dosage <- NULL
    pse_mat <- NULL
    if (have_truth_haplotypes) {
        pse_mat <- array(NA, c(n_seek_its * (nGibbsSamples + 1), 6))
        colnames(pse_mat) <- c("ig", "ii", "pi", "r2", "pse", "disc")
        ## match
        s <- sampleNames[iSample]
        if (!(s %in% dimnames(phase)[[2]])) {
            stop("Something went wrong with phase naming")
        }
        truth_haps <- cbind(phase[, s, 1], phase[, s, 2])
    } else {
        truth_haps <- NULL
    }

    if (have_truth_genotypes) {
        s <- sampleNames[iSample]
        if (!(s %in% colnames(gen))) {
            stop("Something went wrong with gen naming")
        }
        truth_gen <- gen[, s, drop = FALSE]
    } else {
        truth_gen <- NULL
    }
    

    if (hla_run) {
        ## print_message("SPECIAL HLA CODE SIMON")
        gamma_total <- array(0, nrow(hapMatcher))
        list_of_gammas <- as.list(1:nGibbsSamples)
    }

    ## can just set this here
    if (is.null(phase) & !is.null(gen)) {
        have_truth_genotypes <- TRUE ## make UNIQUE to genotype having
    } else {
        have_truth_genotypes <- FALSE
    }

    ## leave NULL, add to it later if needed
    if (plot_per_sample_likelihoods) {
        for_likelihood_plotting <- as.list(1:(nGibbsSamples + 1))
    }

    if (plot_p1) {
        record_read_label_usage <- TRUE ## check this too
        p1_store <- lapply(1:(nGibbsSamples + 1), function(x) {
            as.list(1:n_seek_its)
        })
        read_store <- p1_store
        return_p1 <- TRUE
    } else {
        return_p1 <- FALSE
    }

    ## don't need this for routine use - or do better matching!
    ## truth_g <- as.integer(truth_gen[, sampleNames[iSample]])
    for(i_gibbs_sample in 1:(nGibbsSamples + 1)) {

        if (i_gibbs_sample == (nGibbsSamples + 1)) {
            print_message("Phasing it")
            phasing_it <- TRUE
        } else {
            phasing_it <- FALSE
        }

        outplotprefix <- file.path(outputdir, "plots", paste0("haps.", sample_name, ".", regionName, "_igs.", i_gibbs_sample, "."))

        if (record_read_label_usage & !phasing_it) {
            ## init, truth, after each gibbs, after each seek
            cols <- c("minp", "maxp", "pos", "truth", "uncertain", "init", paste0("gibbs", 1:n_seek_its))
            read_label_matrix <- array(NA, c(length(sampleReads), length(cols)))
            colnames(read_label_matrix) <- cols
            rownames(read_label_matrix) <- sampleReadsInfo[, "qname"]
            read_label_matrix[, "pos"] <- L_grid[sapply(sampleReads, function(x) x[[2]]) + 1]
            x <- t(sapply(sampleReads, function(sampleRead) {L[range(sampleRead[[4]]) + 1]}))
            read_label_matrix[, "minp"] <- x[, 1]
            read_label_matrix[, "maxp"] <- x[, 2]
        }


        ##
        ## possibly do truth
        ## note: ignore this
        if (have_truth_haplotypes) {

            if (verbose) {
                print_message(paste0("i_gibbs=", i_gibbs_sample, ", truth"))
            }
            truth_label_set <- determine_a_set_of_truth_labels(
                sampleReads = sampleReads,
                truth_hap1 = truth_haps[, 1],
                truth_hap2 = truth_haps[, 2],
                maxDifferenceBetweenReads = maxDifferenceBetweenReads
            )
            truth_labels <- truth_label_set[["truth_labels"]]
            uncertain_truth_labels <- truth_label_set[["uncertain_truth_labels"]]

            if ((i_gibbs_sample == 1) && estimate_bq_using_truth_read_labels) {
                bq_result <- estimate_bq(truth_labels = truth_labels, sampleReads = sampleReads, truth_haps = truth_haps)
                print(bq_result)
            }

            truth_all <- impute_using_everything(
                eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                eMatDH_special_matrix = eMatDH_special_matrix,
                use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                H = truth_labels,
                sampleReads = sampleReads,
                rhb_t = rhb_t,
                nSNPs = nSNPs,
                full_alphaHat_t = full_alphaHat_t,
                full_betaHat_t = full_betaHat_t,
                full_gamma_t = full_gamma_t,
                full_transMatRate_t_H  = full_transMatRate_t_H,
                distinctHapsB = distinctHapsB,
                distinctHapsIE = distinctHapsIE,
                hapMatcher = hapMatcher,
                hapMatcherR = hapMatcherR,
                use_hapMatcherR = use_hapMatcherR,
                eMatDH_special_grid_which = eMatDH_special_grid_which,
                eMatDH_special_values_list = eMatDH_special_values_list,
                ref_error = ref_error,
                make_plots = make_plots,
                outplotprefix = outplotprefix,
                have_truth_haplotypes = have_truth_haplotypes,
                have_truth_genotypes = have_truth_genotypes,
                truth_haps = truth_haps,
                truth_gen = truth_gen,
                truth_labels = truth_labels,
                uncertain_truth_labels = uncertain_truth_labels,
                L_grid = L_grid,
                L = L,
                inRegion2 = inRegion2,
                cM_grid = cM_grid,
                ancAlleleFreqAll = ancAlleleFreqAll,
                full_gammaSmall_t = full_gammaSmall_t,
                full_gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                return_good_haps = FALSE,
                plot_description = "0.truth",
                return_dosage = TRUE,
                sample_name = sample_name,
                smooth_cm = smooth_cm,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                minGLValue = minGLValue,
                suppressOutput = suppressOutput,
                use_eigen = use_eigen
            )

        } else {
            truth_labels <- NULL
            uncertain_truth_labels <- NULL
        }


        ##
        ## prepare output matrix
        ##
        if (record_interim_dosages) {
            cols <- paste0(rep(c("gibbs", "all"), n_seek_its), rep(1:n_seek_its, each = 2))
            if (have_truth_haplotypes) {
                cols <- c("truth", cols)
            }
            dosage_matrix <- array(NA, c(nSNPs, length(cols)))
            colnames(dosage_matrix) <- cols
            if (have_truth_haplotypes) {
                dosage_matrix[, "truth"] <- truth_all[["dosage1"]] + truth_all[["dosage2"]]
            }
        }

        if (record_read_label_usage) {
            if (have_truth_haplotypes) {
                read_label_matrix[, "truth"] <- truth_labels
                read_label_matrix[, "uncertain"] <- as.numeric(uncertain_truth_labels)
            }
        }
        i_it <- 1

        ##
        ## now loop, first on subset, then on all
        ##
        for(i_it in 1:n_seek_its) {

            ## here it is 1
            n_gibbs_starts <- 1
            if (i_it == 1 & !phasing_it) {
                which_haps_to_use <- sort(sample(1:K, Ksubset))
                double_list_of_starting_read_labels <- list(
                    lapply(1:n_gibbs_starts, function(i) {
                        H <- sample(c(1, 2), length(sampleReads), replace = TRUE)
                        return(H)
                    })
                )
                gibbs_initialize_iteratively <- TRUE
            } else {
                ## if phasing it will suck these up! same for "which_haps_to_use" (though is that bad?)
                ## which_haps_to_use <- which_haps_to_use
                double_list_of_starting_read_labels <- list(list(read_labels))
                gibbs_initialize_iteratively <- FALSE
            }

            if (i_it == 1 & record_read_label_usage) {
                read_label_matrix[, "init"] <- double_list_of_starting_read_labels[[1]][[1]]
            }

            if (verbose) {
                print_message(paste0("i_gibbs=", i_gibbs_sample, ", i_it = ", i_it, " small gibbs"))
            }

            if (record_interim_dosages) {
                return_genProbs <- TRUE
                return_hapProbs <- TRUE
            } else {
                return_genProbs <- FALSE
                return_hapProbs <- FALSE
            }

            ## print_message(paste0("i_gibbs, which_haps_to_use = ", paste(which_haps_to_use, collapse = ",")))
            if (zilong | use_mspbwt) {
                return_hapProbs <- TRUE
            }

            gibbs_iterate <- impute_one_sample(
                eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                eMatDH_special_matrix = eMatDH_special_matrix,
                use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                distinctHapsB = distinctHapsB,
                distinctHapsIE = distinctHapsIE,
                hapMatcher = hapMatcher,
                hapMatcherR = hapMatcherR,
                use_hapMatcherR = use_hapMatcherR,
                rhb_t = rhb_t,
                ref_error = ref_error,
                nSNPs = nSNPs,
                sampleReads = sampleReads,
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
                n_gibbs_starts = n_gibbs_starts,
                small_ref_panel_gibbs_iterations = small_ref_panel_gibbs_iterations,
                n_gibbs_sample_its = 1,
                double_list_of_starting_read_labels = double_list_of_starting_read_labels,
                small_ref_panel_block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
                perform_block_gibbs = TRUE,
                make_plots = make_plots,
                wif0 = wif0,
                grid_has_read = grid_has_read,
                plot_description = paste0("it", i_it, ".gibbs"),
                ancAlleleFreqAll = ancAlleleFreqAll,
                grid = grid,
                L_grid = L_grid,
                L = L,
                inRegion2 = inRegion2,
                cM_grid = cM_grid,
                outplotprefix = outplotprefix,
                have_truth_haplotypes = have_truth_haplotypes,
                truth_haps = truth_haps,
                have_truth_genotypes = have_truth_genotypes,
                truth_gen = truth_gen,
                truth_labels = truth_labels,
                uncertain_truth_labels = uncertain_truth_labels,
                verbose = FALSE,
                maxEmissionMatrixDifference = 1e100,
                return_p_store = FALSE,
                return_p1 = return_p1,
                return_extra = FALSE,
                return_genProbs = return_genProbs,
                return_hapProbs = return_hapProbs,
                return_gibbs_block_output = FALSE,
                gibbs_initialize_iteratively = gibbs_initialize_iteratively,
                gibbs_initialize_at_first_read = FALSE,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                rescale_eMatRead_t = TRUE,
                rescale_eMatGrid_t = FALSE,
                Jmax = 10000,
                suppressOutput = suppressOutput,
                shuffle_bin_radius = shuffle_bin_radius,
                make_plots_block_gibbs = make_plots_block_gibbs,
                sample_name = sample_name,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                use_small_eHapsCurrent_tc = use_small_eHapsCurrent_tc,
                use_sample_is_diploid = use_sample_is_diploid,
                small_ref_panel_skip_equally_likely_reads = small_ref_panel_skip_equally_likely_reads,
                small_ref_panel_equally_likely_reads_update_iterations = small_ref_panel_equally_likely_reads_update_iterations,
                i_it = i_it,
                i_gibbs_sample = i_gibbs_sample,
                ff0_shard_check_every_pair = ff0_shard_check_every_pair,
                zilong = zilong
            )

            if (plot_p1) {
                p1_store[[i_gibbs_sample]][[i_it]] <- gibbs_iterate[["p1"]]
                read_store[[i_gibbs_sample]][[i_it]] <- gibbs_iterate[["pH"]]
            }

            if (hla_run) {
                ## final phasing it, save gamma
                if (
                (i_it == n_seek_its)
                ) {
                    ## print_message("HLA SPECIAL CODE SIMON")
                    return_gamma_t <- TRUE
                    full_gamma_t <- array(0, c(K, nGrids))
                } else {
                    return_gamma_t <- FALSE
                }
            }

            if (plot_per_sample_likelihoods) {
                if (i_it == 1) {
                    for_likelihood_plotting[[i_gibbs_sample]] <- as.list(1:n_seek_its)
                }
                for_likelihood_plotting[[i_gibbs_sample]][[i_it]] <- gibbs_iterate$per_it_likelihoods
            }

            ## am here
            ## continue to work on this. weird error below hmm...
            if (verbose) {
                print_message(paste0("i_gibbs=", i_gibbs_sample, ", i_it = ", i_it, " full"))
            }

            if (i_it < n_seek_its) {
                return_good_haps <- TRUE
            } else {
                ## keep on, keep selecting good haps, phasing sucks these up
                return_good_haps <- TRUE
            }

            read_labels <- gibbs_iterate$double_list_of_ending_read_labels[[1]][[1]]
            previously_selected_haplotypes <- sample(which_haps_to_use, Ksubset - Knew)
            # is equavaliant to hapProbs_t
            return_dosage <- (have_truth_haplotypes | record_interim_dosages | (i_it > n_burn_in_seek_its))

            igibbs <- (i_gibbs_sample - 1) * n_seek_its + i_it ## 1-based

            which_haps_to_use_zilong_A <- NULL
            which_haps_to_use_zilong_B <- NULL
            which_haps_to_use_zilong_mspbwt <- NULL
            which_haps_to_use_quilt1 <- NULL
            
            if (zilong | (!use_mspbwt && make_heuristic_plot)) {

                ##print(paste0("zilong = ", zilong))
                ##print(paste0("(!use_mspbwt && make_heuristic_plot) = ", (!use_mspbwt && make_heuristic_plot)))

                Kfull <- nrow(hapMatcher)
                hap1 <- gibbs_iterate$hapProbs_t[1, ]
                hap2 <- gibbs_iterate$hapProbs_t[2, ]
                
                if (heuristic_approach == "A" | make_heuristic_plot) {

                    which_haps_to_use <- select_new_haps_zilong_msp(
                        gibbs_iterate$hapProbs_t,
                        igibbs = igibbs,
                        outputdir = outputdir,
                        Kfull = Kfull,
                        Knew = Knew,
                        msp = msp,
                        mspbwtB = mspbwtB,
                        mspbwtL = mspbwtL,
                        mspbwtM = mspbwtM,
                        nGrids = nGrids
                    )
                    
                    which_haps_to_use_zilong_A <- which_haps_to_use
                    
                }

                if (heuristic_approach == "B" | make_heuristic_plot) {

                    which_haps_to_use <- select_new_haps_zilong_msp_robbie_version(
                        gibbs_iterate$hapProbs_t,
                        igibbs = igibbs,
                        outputdir = outputdir,
                        Kfull = Kfull,
                        Knew = Knew,
                        msp = msp,
                        mspbwtB = mspbwtB,
                        mspbwtL = mspbwtL,
                        mspbwtM = mspbwtM,
                        nGrids = nGrids
                    )

                    which_haps_to_use_zilong_B <- which_haps_to_use

                }
                
            }

            if (use_mspbwt | (!zilong && make_heuristic_plot)) {

                ## for testing purposes
                if (use_splitreadgl) {
                    impute_all <- impute_using_split_reads_and_small_ref_panel(
                        H = read_labels,
                        which_haps_to_use = which_haps_to_use,
                        sampleReads = sampleReads,
                        rhb_t = rhb_t,
                        nSNPs = nSNPs,
                        full_alphaHat_t = full_alphaHat_t,
                        full_betaHat_t = full_betaHat_t,
                        full_gamma_t = full_gamma_t,
                        full_gammaSmall_t = full_gammaSmall_t,
                        full_gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                        full_transMatRate_t_H = full_transMatRate_t_H,
                        distinctHapsB = distinctHapsB,
                        distinctHapsIE = distinctHapsIE,
                        hapMatcher = hapMatcher,
                        eMatDH_special_grid_which = eMatDH_special_grid_which,
                        eMatDH_special_values_list = eMatDH_special_values_list,
                        ref_error = ref_error,
                        make_plots = FALSE, ## these plots are pretty useless? as plots?
                        outplotprefix = outplotprefix,
                        have_truth_haplotypes = have_truth_haplotypes,
                        truth_haps = truth_haps,
                        truth_labels = truth_labels,
                        uncertain_truth_labels = uncertain_truth_labels,
                        L_grid = L_grid,
                        L = L,
                        inRegion2 = inRegion2,
                        cM_grid = cM_grid,
                        ancAlleleFreqAll = ancAlleleFreqAll,
                        plot_description = paste0("it", i_it, ".full"),
                        return_good_haps = return_good_haps,
                        Knew = Knew,
                        return_dosage  = return_dosage,
                        previously_selected_haplotypes = previously_selected_haplotypes,
                        K_top_matches = K_top_matches,
                        heuristic_match_thin = heuristic_match_thin,
                        return_gamma_t = return_gamma_t,
                        sample_name = sample_name,
                        smooth_cm = smooth_cm ,
                        regionStart = regionStart,
                        regionEnd = regionEnd,
                        buffer = buffer,
                        minGLValue = minGLValue,
                        suppressOutput = suppressOutput
                    )
                    hap1 <- impute_all[["dosage1"]]
                    hap2 <- impute_all[["dosage2"]]
                } else {
                    hap1 <- gibbs_iterate$hapProbs_t[1, ]
                    hap2 <- gibbs_iterate$hapProbs_t[2, ]
                }

                ## this version here
                hapProbs_t <- rbind(hap1, hap2)
                Kfull <- nrow(hapMatcher)

                ## file <- paste0("~/temp.", i_gibbs_sample, ".", i_it, ".RData")
                ## print(paste0("saving to file:", file))
                ## save(
                ##     hapProbs_t,
                ##     use_hapMatcherR,
                ##     Knew,
                ##     Kfull,
                ##     mspbwtL,
                ##     mspbwtM,
                ##     file = file
                ## )

                if (heuristic_approach == "A" | make_heuristic_plot) {

                    which_haps_to_use <- select_new_haps_mspbwt_v3(
                        hapProbs_t = hapProbs_t,
                        hapMatcher = hapMatcher,
                        hapMatcherR = hapMatcherR,
                        use_hapMatcherR = use_hapMatcherR,
                        ms_indices = ms_indices,
                        Knew = Knew,
                        Kfull = Kfull,
                        mspbwtL = mspbwtL,
                        mspbwtM = mspbwtM,
                        heuristic_approach = "A"
                    )
                    which_haps_to_use_mspbwt_A <- which_haps_to_use

                }

                if (heuristic_approach == "B" | make_heuristic_plot) {

                    which_haps_to_use <- select_new_haps_mspbwt_v3(
                        hapProbs_t = hapProbs_t,
                        hapMatcher = hapMatcher,
                        hapMatcherR = hapMatcherR,
                        use_hapMatcherR = use_hapMatcherR,
                        ms_indices = ms_indices,
                        Knew = Knew,
                        Kfull = Kfull,
                        mspbwtL = mspbwtL,
                        mspbwtM = mspbwtM,
                        heuristic_approach = "B"
                    )
                    which_haps_to_use_mspbwt_B <- which_haps_to_use

                }
                
                

            }

            if ((!zilong && !use_mspbwt) | make_heuristic_plot) {
                
                impute_all <- impute_using_everything(
                    eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                    eMatDH_special_matrix = eMatDH_special_matrix,
                    use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                    H = read_labels,
                    sampleReads = sampleReads,
                    rhb_t = rhb_t,
                    nSNPs = nSNPs,
                    full_alphaHat_t = full_alphaHat_t,
                    full_betaHat_t = full_betaHat_t,
                    full_gamma_t = full_gamma_t,
                    full_gammaSmall_t = full_gammaSmall_t,
                    full_gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                    full_transMatRate_t_H = full_transMatRate_t_H,
                    distinctHapsB = distinctHapsB,
                    distinctHapsIE = distinctHapsIE,
                    hapMatcher = hapMatcher,
                    hapMatcherR = hapMatcherR,
                    use_hapMatcherR = use_hapMatcherR,
                    eMatDH_special_grid_which = eMatDH_special_grid_which,
                    eMatDH_special_values_list = eMatDH_special_values_list,
                    ref_error = ref_error,
                    make_plots = FALSE, ## these plots are pretty useless? as plots?
                    outplotprefix = outplotprefix,
                    have_truth_haplotypes = have_truth_haplotypes,
                    have_truth_genotypes = have_truth_genotypes,
                    truth_haps = truth_haps,
                    truth_gen = truth_gen,
                    truth_labels = truth_labels,
                    uncertain_truth_labels = uncertain_truth_labels,
                    L_grid = L_grid,
                    L = L,
                    inRegion2 = inRegion2,
                    cM_grid = cM_grid,
                    ancAlleleFreqAll = ancAlleleFreqAll,
                    plot_description = paste0("it", i_it, ".full"),
                    return_good_haps = return_good_haps,
                    Knew = Knew,
                    return_dosage  = return_dosage,
                    previously_selected_haplotypes = previously_selected_haplotypes,
                    K_top_matches = K_top_matches,
                    heuristic_match_thin = heuristic_match_thin,
                    return_gamma_t = return_gamma_t,
                    sample_name = sample_name,
                    smooth_cm = smooth_cm ,
                    regionStart = regionStart,
                    regionEnd = regionEnd,
                    buffer = buffer,
                    minGLValue = minGLValue,
                    suppressOutput = suppressOutput,
                    use_eigen = use_eigen
                )
                which_haps_to_use <- c(previously_selected_haplotypes, impute_all$new_haps)
                hap1 <- impute_all[["dosage1"]]
                hap2 <- impute_all[["dosage2"]]

                which_haps_to_use_quilt1 <- which_haps_to_use
            }

            if (make_heuristic_plot) {

                hapProbs_t <- gibbs_iterate$hapProbs_t 
                
                compare_heuristic_approaches(
                    hapProbs_t = hapProbs_t,
                    which_haps_to_use_quilt1 = which_haps_to_use_quilt1,
                    which_haps_to_use_zilong_A = which_haps_to_use_zilong_A,
                    which_haps_to_use_zilong_B = which_haps_to_use_zilong_B,
                    which_haps_to_use_mspbwt_A = which_haps_to_use_mspbwt_A,
                    which_haps_to_use_mspbwt_B = which_haps_to_use_mspbwt_B,                    
                    hapMatcherR = hapMatcherR,
                    hapMatcher = hapMatcher,
                    use_hapMatcherR = use_hapMatcherR,
                    distinctHapsB = distinctHapsB,
                    eMatDH_special_matrix = eMatDH_special_matrix,
                    eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                    outputdir = outputdir,
                    sample_name = sample_name,
                    regionName = regionName,
                    i_gibbs_sample = i_gibbs_sample,
                    i_it = i_it,
                    nGrids = nGrids
                )

            }

            
            if (record_interim_dosages) {
                dosage_matrix[, paste0("gibbs", i_it)] <- get_dosages_from_fbsoL(gibbs_iterate)
                dosage_matrix[, paste0("all", i_it)] <- hap1 + hap2
            }

            if (record_read_label_usage) {
                read_label_matrix[, paste0("gibbs", i_it)] <- read_labels
            }

            ## if not a phasing iteration, and past burn in, store!
            if (!phasing_it && (i_it > n_burn_in_seek_its)) {
                ## print(paste0("temptemp - ", i_gibbs_sample, " ", i_it, " ", n_seek_its, " ", n_burn_in_seek_its))
                dosage <- dosage + hap1 + hap2
                gp_t <- gp_t +
                    rbind((1 - hap1) * (1 - hap2), (1 - hap1) * hap2 + hap1 * (1 - hap2), hap1 * hap2)
                nDosage <- nDosage + 1

                if (have_truth_haplotypes) {
                    w <- i_it + n_seek_its * (i_gibbs_sample - 1)
                    x <- calculate_pse_and_r2_during_gibbs(inRegion2 = inRegion2, hap1 = hap1, hap2 = hap2, truth_haps = truth_haps, af = af, verbose = verbose, impute_rare_common = impute_rare_common, all_snps = FALSE)
                    pse_mat[w, ] <- c(i_gibbs_sample, i_it, as.integer(phasing_it), x)
                } else if (have_truth_genotypes) {
                    r2 <-  round(cor((dosage / nDosage)[inRegion2] - 2 * af[inRegion2], gen[inRegion2, sample_name] - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
                    print_message(paste0("Current accuracy for this gibbs sample for ", sample_name, ", r2:", r2))
                }

            }


        }


        if (impute_rare_common) {

            if (verbose) {
                print_message(paste0("i_gibbs=", i_gibbs_sample, ", i_it = ", i_it, " small gibbs (all SNPs)"))
            }

            out_rare_common <- impute_final_gibbs_with_rare_common(
                special_rare_common_objects = special_rare_common_objects,
                special_rare_common_objects_per_core = special_rare_common_objects_per_core,
                allSNP_sampleReads = allSNP_sampleReads,
                hap1 = hap1,
                hap2 = hap2,
                pos_all = pos_all,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                hapMatcher = hapMatcher,                
                hapMatcherR = hapMatcherR,
                use_hapMatcherR = use_hapMatcherR,
                distinctHapsIE = distinctHapsIE,
                eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                eMatDH_special_matrix = eMatDH_special_matrix,
                Ksubset = Ksubset,
                ref_error = ref_error,
                which_haps_to_use = which_haps_to_use,
                small_ref_panel_gibbs_iterations = small_ref_panel_gibbs_iterations,
                small_ref_panel_block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
                allSNP_wif0 = allSNP_wif0,
                allSNP_grid_has_read = allSNP_grid_has_read,
                make_plots = make_plots,
                outplotprefix = outplotprefix,
                have_truth_haplotypes = have_truth_haplotypes,
                truth_haps_all = truth_haps_all,
                have_truth_genotypes = have_truth_genotypes,
                truth_gen_all = truth_gen_all,
                truth_labels_all = truth_labels_all,
                uncertain_truth_labels_all = uncertain_truth_labels_all,
                shuffle_bin_radius = shuffle_bin_radius,
                make_plots_block_gibbs = make_plots_block_gibbs,
                sample_name = sample_name,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                i_it = i_it,
                i_gibbs_sample = i_gibbs_sample,
                ff0_shard_check_every_pair = ff0_shard_check_every_pair,
                sampleNames = sampleNames,
                iSample = iSample,
                phase_all = phase_all
            )

            hap1_all <- out_rare_common[["hap1"]]
            hap2_all <- out_rare_common[["hap2"]]

            if (!phasing_it && (i_it > n_burn_in_seek_its)) {

                dosage_all <- dosage_all + hap1_all + hap2_all
                gp_t_all <- gp_t_all + rbind(
                    (1 - hap1_all) * (1 - hap2_all),
                    (1 - hap1_all) * hap2_all + hap1_all * (1 - hap2_all),
                    hap1_all * hap2_all
                )
                ## nDosage <- nDosage + 1

                calculate_pse_and_r2_rare_common(                
                    hap1_all = hap1_all,
                    hap2_all = hap2_all,
                    have_truth_haplotypes = have_truth_haplotypes,
                    truth_haps_all = truth_haps_all,
                    have_truth_genotypes = have_truth_genotypes,
                    truth_gen_all = truth_gen_all,
                    special_rare_common_objects = special_rare_common_objects,
                    verbose = verbose
                )
                


            }
            

        }

        if (!phasing_it) {
            ## for phasing bit
            read_label_matrix_all[, i_gibbs_sample] <- read_labels
            read_label_matrix_conf[, i_gibbs_sample] <- assess_ability_of_reads_to_be_confident(
                hap1 = hap1,
                hap2 = hap2,
                sampleReads = sampleReads,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                minrp = 0.95,
                minmp = 0.95
            )
            if (record_read_label_usage) {
                super_out_read_labels[[i_gibbs_sample]] <- read_label_matrix
            }
            if (record_interim_dosages) {
                super_out_dosage_matrix[[i_gibbs_sample]] <- dosage_matrix
            }
        }

        ## save last
        if (i_gibbs_sample == nGibbsSamples) {
            can_hap <- nGibbsSamples
            out_best_labels <- determine_best_read_label_so_far(
                read_label_matrix_all = read_label_matrix_all,
                read_label_matrix_conf = read_label_matrix_conf,
                nReads = length(sampleReads),
                nGibbsSamples = nGibbsSamples,
                verbose = verbose,
                can_hap = can_hap
            )
            if (verbose) {
                x <- out_best_labels$flip_matrix[, can_hap]
                if (length(x) > 0) {
                    print_message(paste0("There are ", sum(x), " out of ", length(x), " regions that have been flipped by consensus"))
                }
            }
            read_labels <- out_best_labels$read_labels
        }

        if (phasing_it) {
            if (!impute_rare_common) {
                ## just save relevant stuff here
                out <- recast_haps(hd1 = hap1, hd2 = hap2, gp = t(gp_t))
                ## over-write here
                hap1 <- out$hd1
                hap2 <- out$hd2
                ## 
                phasing_haps <- cbind(hap1, hap2)
                ## 
                phasing_dosage <- hap1 + hap2
                phasing_read_labels <- read_labels
            } else {
                ## just save relevant stuff here
                out <- recast_haps(hd1 = hap1_all, hd2 = hap2_all, gp = t(gp_t_all))
                ## over-write here
                hap1_all <- out$hd1
                hap2_all <- out$hd2
                phasing_haps_all <- cbind(hap1_all, hap2_all)
            }
        }

        if (hla_run) {
            gamma1 <- impute_all$fbsoL$gammaMT_t
            gamma2 <- impute_all$fbsoL$gammaMU_t
            if (is.na(gamma_physically_closest_to)) {
                iGrid <- round(nGrids / 2)
            } else {

            }
            gamma1 <- gamma1[, iGrid]
            gamma2 <- gamma2[, iGrid]
            if (!phasing_it) {
                gamma_total <- gamma_total + gamma1 + gamma2
            }
            if (i_gibbs_sample <= nGibbsSamples) {
                list_of_gammas[[i_gibbs_sample]] <- list(gamma1, gamma2)
            }
        }

    }

    if (plot_p1) {
        ##save(p1_store, read_store, super_out_read_labels, file = "~/temp.RData")
        p1 <- p1_store[[1]][[1]]
        print(paste0(sum(p1 == 0) / (prod(dim(p1)))))
        if (small_ref_panel_skip_equally_likely_reads) {
            save(p1_store, file = "~/temp.RData")
        }
    }

    if(have_truth_haplotypes) {
        imputed_truth_haplotypes <- cbind(truth_all[["dosage1"]], truth_all[["dosage2"]])
    } else {
        imputed_truth_haplotypes <- NULL
    }

    dosage <- dosage / nDosage
    gp_t <- gp_t / nDosage

    if (impute_rare_common) {
        dosage_all <- dosage_all / nDosage
        gp_t_all <- gp_t_all / nDosage
    }

    
    ##
    ## print final accuracies
    ##
    if (impute_rare_common) {
        final_phasing_accuracy_calculation_rare_common(
            have_truth_haplotypes = have_truth_haplotypes,
            have_truth_genotypes = have_truth_genotypes,
            truth_haps_all = truth_haps_all,
            dosage_all = dosage_all,
            hap1_all = hap1_all,
            hap2_all = hap2_all,
            gen_all = gen_all,
            sampleNames = sampleNames,
            iSample = iSample,
            special_rare_common_objects = special_rare_common_objects,
            sample_name = sample_name
        )     
    } else {    
        final_phasing_accuracy_calculation(have_truth_haplotypes = have_truth_haplotypes, have_truth_genotypes = have_truth_genotypes, truth_haps = truth_haps, inRegion2 = inRegion2, dosage = dosage, af = af, phasing_haps = phasing_haps, gen = gen, sample_name = sample_name)
    }



    ## optionally plot here
    if (plot_per_sample_likelihoods) {
        plot_of_likelihoods_across_samplings_and_seek_its(
            outputdir = outputdir,
            for_likelihood_plotting = for_likelihood_plotting,
            nGibbsSamples = nGibbsSamples,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            block_gibbs_iterations = block_gibbs_iterations,
            n_seek_its = n_seek_its,
            sample_name = sample_name,
            regionName = regionName,
            truth_likelihoods = NULL,
            ylab = "p_H_given_O_L_up_to_C"
        )
    }


    ## perform switch over here
    if (impute_rare_common) {

        ##save(sampleReads, pos, pos_all, allSNP_sampleReads, file = "~/temp.RData")
        ##stop("WER")
        gp_t <- gp_t_all
        phasing_haps <- phasing_haps_all
        sampleReads <- allSNP_sampleReads
        nSNPs <- ncol(gp_t)


    }

    
    eij <- round(gp_t[2, ] + 2 * gp_t[3, ], 3) ## prevent weird rounding issues
    fij <- round(gp_t[2, ] + 4 * gp_t[3, ], 3) ##

    max_gen <- get_max_gen_rapid(gp_t)

    ##
    ## for allele count, from buildAlleleCount_subfunction in STITCH
    ##
    a <- unlist(sapply(sampleReads,function(x) x[[3]]))
    b <- unlist(sapply(sampleReads,function(x) x[[4]]))
    bqProbs <- STITCH::convertScaledBQtoProbs(matrix(a,ncol=1))
    ## y is numeric, z = integer counts (can have 0),
    c1 <- increment2N(
        y = as.numeric(bqProbs[,1]),
        z = as.numeric(b),
        yT = as.integer(nrow(bqProbs)),
        xT = as.integer(nSNPs - 1)
    )
    c2 <- increment2N(
        y = as.numeric(bqProbs[,2]),
        z = as.numeric(b),
        yT = as.integer(nrow(bqProbs)),
        xT = as.integer(nSNPs - 1)
    )
    per_sample_alleleCount <- cbind(c2, c1 + c2)


    ## make vcf character thing - NOTE - this is a HACK of previous one - but works OK here!

    if (addOptimalHapsToVCF & have_truth_haplotypes) {
        add_x_2_cols <- TRUE
        x_t <- t(imputed_truth_haplotypes)
    } else {
        add_x_2_cols <- FALSE
        x_t <- matrix()
    }

    per_sample_vcf_col <- STITCH::rcpp_make_column_of_vcf(
        gp_t = gp_t,
        use_read_proportions = FALSE,
        use_state_probabilities = TRUE,
        read_proportions = matrix(),
        q_t = t(phasing_haps),
        add_x_2_cols = add_x_2_cols,
        x_t = x_t
    )
    ## annoying but shouldn't be that slow I hope
    if (output_gt_phased_genotypes) {
        per_sample_vcf_col <-
            paste0(
                round(phasing_haps[, 1]), "|", round(phasing_haps[, 2]),
                substring(per_sample_vcf_col, first = 4, last = 100L)
            )
    }

    ## what to return
    to_return <- list(
        sample_was_imputed = TRUE,
        eij = eij,
        fij = fij,
        max_gen = max_gen,
        per_sample_alleleCount = per_sample_alleleCount,
        per_sample_vcf_col = per_sample_vcf_col,
        super_out_hap_dosages = super_out_hap_dosages,
        super_out_read_labels = super_out_read_labels,
        super_out_dosage_matrix = super_out_dosage_matrix
    )

    ## phasing_read_labels = phasing_read_labels,
    ## imputed_truth_haplotypes = imputed_truth_haplotypes,
    ## pse_mat = pse_mat,
    ## dosage = dosage,
    ## super_out_dosage_matrix = super_out_dosage_matrix,
    ## phasing_haps = phasing_haps,
    ## phasing_dosage = phasing_dosage,
    ## gp_t = gp_t,
    ##

    if (hla_run) {
        ## print_message("More HLA SIMON code")
        to_return[["gamma1"]] <- gamma1
        to_return[["gamma2"]] <- gamma2
        to_return[["gamma_total"]] <- gamma_total
        to_return[["list_of_gammas"]] <- list_of_gammas
    }

    return(to_return)


}


#' @export
modified_calculate_pse <- function(
    test,
    truth,
    LL,
    seed = NULL
) {
    ## for testing
    if (is.null(seed) == FALSE)
        set.seed(seed)
    rownames(test) <- 1:nrow(test)
    which_sites <-
        rowSums(truth == 0 | truth == 1, na.rm = TRUE) == 2 &
        rowSums(truth, na.rm = TRUE) == 1 &
        rowSums(is.na(truth)) == 0
    truth <- truth[which_sites, ]
    test <- test[which_sites, ]
    if (nrow(test) == 0)
        return(NA)
    ## as these sites, discrepency as well
    disc <- sum(round(rowSums(test)) != 1)
    test <- round(test)
    ## round test data for now. choose hets at random
    ## for homs
    test[, 1] <- as.integer(round(test[, 1]))
    test[, 2] <- as.integer(round(test[, 2]))
    testO <- test
    truthO <- truth
    for(i_option in 1:2) {
        test <- testO
        truth <- truthO
        ## specifically remove from consideration double phase switch errors
        w <- rowSums(test) == 1
        w2 <- which(diff(abs(test[w,1] - truth[w,1])) != 0)
        to_remove <- NULL
        if (length(w2) > 0) {
            w3 <- which(diff(w2) == 1)
            if (length(w3) > 0) {
                for(a in w3) {
                    c <- w2[c(a, a + 1)]
                    to_remove <- c(to_remove, which(w)[c])
                }
            }
        }
        ## double pse are two consecutive
        if (length(to_remove) > 0) {
            test <- test[-to_remove, ]
            truth <- truth[-to_remove, ]
        }
        ##
        if (i_option == 1) {
            ## only consider non-discrepent sites
            ## chose best start
            if (test[1, 1] != truth[1, 1])
                test <- test[, c(2, 1)]
            ## calculate number of differences
            w <- rowSums(test) == 1
            if (sum(w) == 0) {
                print("Test has no hets! possibly an error or homo over region, possibly no record dosages turned on in impute_all")
                switches1 <- cbind(i1 = NA, i2 = NA, l1 = NA, l2 = NA)
                phase_errors_def1 <- 0
                phase_sites_def1 <- 0
            } else {
                y <- diff(abs(test[w,1] - truth[w,1])) != 0
                phase_errors_def1 <- sum(y)
                phase_sites_def1 <- sum(w) - 1
                s <- as.integer(rownames(test[w, , drop = FALSE][c(as.logical(y), FALSE), , drop = FALSE]))
                e <- as.integer(rownames(test[w, , drop = FALSE][c(FALSE, as.logical(y)), , drop = FALSE]))
                switches1 <- cbind(i1 = s, i2 = e, l1 = LL[s], l2 = LL[e])
            }
        }
        if (i_option == 2) {
            choose_at_random <- which(rowSums(test) != 1)
            if (length(choose_at_random) > 0) {
                test[choose_at_random, ] <- 0
                r <- sample(
                    c(1, 2),
                    length(choose_at_random),
                    replace = TRUE
                )
                test[cbind(choose_at_random, r)] <- 1
            }
            ## chose best start
            if (test[1, 1] != truth[1, 1])
                test <- test[, c(2, 1)]
            ## calculate number of differences
            phase_errors_def2 <- sum(diff(abs(test[,1] - truth[,1])) != 0)
            phase_sites_def2 <- nrow(test) - 1
        }
    }
    ##
    return(
        list(
            values = c(
                phase_errors_def1 = phase_errors_def1,
                phase_sites_def1 = phase_sites_def1,
                phase_errors_def2 = phase_errors_def2,
                phase_sites_def2 = phase_sites_def2,
                disc_errors = disc,
                dist_n = nrow(test)
            ),
            switches1 = switches1
        )
    )
}


calculate_pse_and_r2_during_gibbs <- function(inRegion2, hap1, hap2, truth_haps, af, verbose = FALSE, impute_rare_common = FALSE, all_snps = FALSE) {
    ## wow, confident
    w <- (inRegion2)
    g <- truth_haps[inRegion2, 1] + truth_haps[inRegion2, 2]
    ## scaled version
    r2 <-  round(cor((hap1 + hap2)[inRegion2] - 2 * af[inRegion2], g - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
    fake_LL <- 1:sum(w)
    values <- modified_calculate_pse(test = round(cbind(hap1, hap2)[w, ]), truth = truth_haps[w, ], LL = fake_LL)$values
    ## discordance
    pse <- values["phase_errors_def1"] / values["phase_sites_def1"]
    pse <- round(100 * pse, 1)
    disc <- round(100 * values["disc_errors"] / values["dist_n"], 1)
    if (impute_rare_common) {
        if (all_snps && verbose) {
            print_message(paste0("r2:", r2, ", PSE:", pse, "%, disc:", disc, "% (all SNPs)"))
        } else {
            print_message(paste0("r2:", r2, ", PSE:", pse, "%, disc:", disc, "% (common SNPs only)"))
        }
    } else {
        if (verbose) {
            print_message(paste0("r2:", r2, ", PSE:", pse, "%, disc:", disc, "%"))
        }
    }
    return(c(r2 = r2, pse = as.numeric(pse), disc = as.numeric(disc)))
}



assess_ability_of_reads_to_be_confident <- function(hap1, hap2, sampleReads, maxDifferenceBetweenReads, minrp = 0.95, minmp = 0.95) {
    ##
    p <- calculate_eMatRead_t_vs_two_haplotypes(
        sampleReads,
        hap1 = hap1,
        hap2 = hap2,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads ,
        rescale_eMatRead_t = FALSE
    )
    p1 <- p[1, ]
    p2 <- p[2, ]
    ## maxmimum
    mp <- p1
    mp[p2 > p1] <- p2[p2 > p1]
    ## ratio
    rp <- p1 / (p1 + p2)
    rp[is.na(rp)] <- 0.5 ## if BOTH 0 i.e. suuuper unlikely
    rp[rp < 0.5] <- 1 - rp[rp < 0.5]
    conf <- (rp > minrp)  ## & (mp >
    return(conf)
}


## for(i_gibbs_sample in 1:nGibbsSamples) {
##     hap_matrix <- super_out_hap_dosages[[i_gibbs_sample]]
##     read_label_matrix <- super_out_read_labels[[i_gibbs_sample]]
##     ##
##     gp_t <- rbind(
##     (1 - h1) * (1 - h2),
##     (1 - h1) * (h2) + (h1) * (1 - h2),
##     (h1) * (h2)
##     )
##     eij <- round(gp_t[2, ] + 2 * gp_t[3, ], 3) ## prevent weird rounding issues
##     fij <- round(gp_t[2, ] + 4 * gp_t[3, ], 3) ##
##     thetaHat <- eij / 2
##     info <- 1 - (fij - eij**2) / 2 * thetaHat * (1-thetaHat)
## }



determine_best_read_label_so_far <- function(
    read_label_matrix_all,
    read_label_matrix_conf,
    nReads,
    nGibbsSamples,
    verbose,
    can_hap = 1
 ) {
    ##
    ## for default
    read_labels <- as.integer(read_label_matrix_all[, can_hap])
    flip_matrix <- array(FALSE, c(0, nGibbsSamples)) ##
    default_out <- list(
        read_labels = read_labels,
        flip_matrix = flip_matrix
    )
    ## how to assess "switch error rate"?
    ## do the reads partition themselves nicely
    a <- read_label_matrix_all
    a[!read_label_matrix_conf] <- NA
    a <- a[rowSums(is.na(a)) == 0, , drop = FALSE]
    if (nrow(a) < 10) {
        ## does not really matter
        if (verbose) {
            print_message("Insufficiently many confident reads for aggregating across runs")
        }
        return(default_out)
    }
    ## consider differences vs can_hap
    ## can_hap <- 1
    can <- a[, can_hap]
    a <- a - can
    s <- which(diff(rowSums(abs(a))) != 0) ## where they change, e.g. if a 39 here, it switches between 39 and 40
    ## probably some minimum number here
    if (length(s) == 0) {
        ## doesn't really need a message, the runs are the same at confident reads
        ##if (verbose) {
        ##    print_message("No differences between runs, Insufficiently many confident reads for aggregating across runs")
        ## }
        return(default_out)
    }
    ## so, e.g, want majority vote in each?
    ## what if I do a majority vote option
    ##
    s1 <- c(1, s + 1) #[-(length(s) + 1)]
    e1 <- c(s1[-1] - 1, nrow(a))
    ##
    flip_matrix <- array(FALSE, c(length(s1), nGibbsSamples))
    a2 <- a
    a2[] <- NA
    ## fuck it, just flip em
    for(i in 2:length(s1)) {
        s2 <- s1[i]
        prev <- a[e1[i - 1], ]
        cur <- a[s1[i], ]
        changed <- which(cur != 0)
        ## changed <- which(prev != cur)
        w <- s1[i]:nrow(a)
        if (length(changed) > 0) {
            if (length(changed) <= (nGibbsSamples / 2)) {
                ## trust canonical
                for(c1 in changed) {
                    reverted <- a[w, c1] + can[w]
                    reverted <- 3 - reverted
                    a[s1[i]:nrow(a), c1] <- reverted - can[w]
                }
                stopifnot(!(can_hap %in% changed))
            } else {
                ## here we change canonical, then any others identified
                changed <- which(cur == 0)
                stopifnot((can_hap %in% changed))
                for(c1 in changed) {
                    reverted <- a[w, c1] + can[w]
                    reverted <- 3 - reverted
                    a[w, c1] <- reverted - can[w]
                }
                reverted <- a[w, ] + can[w]
                ## now change canonical
                can[w] <- 3 - can[w]
                ## now re-set all of a
                a[w, ] <- reverted - can[w]
            }
        }
        flip_matrix[i, changed] <- TRUE
    }
    ## now - try using one of these and re-run and check phase
    read_labels <- as.integer(read_label_matrix_all[, can_hap])
    to_flip <- which(flip_matrix[, can_hap])
    if (length(to_flip) > 0) {
        for(i in to_flip) {
            read_labels[s1[i]:nReads] <- 3 - read_labels[s1[i]:nReads]
        }
    }
    return(
        list(
            read_labels = read_labels,
            flip_matrix = flip_matrix
        )
    )
}



if (1 == 0) {

    head(out$dosage_matrix)

    breaks2 <- sort(unique(c(
        c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
        c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
        c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
        seq(0.1, 0.5, length.out = 5),
        1 - seq(0.1, 0.5, length.out = 5),
        1 - c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1)
    )))
    d <- hap1 + hap2
    d <- dosage_matrix[, "truth"]
    r2_by_freq(breaks2, af, truth_g, d, flip = FALSE)


    dosage_matrix[, "gibbs0"]
    m <- round(cbind(
        r2_by_freq(breaks, af, truth_g, dosage_matrix[, "gibbs0"])[, c("n", "nA")],
        apply(dosage_matrix, 2, function(dosage) {
            r2_by_freq(breaks, af, truth_g, dosage, flip = TRUE)[, "simple"]
        })
    ), 3)
    m



}

determine_a_set_of_truth_labels <- function(sampleReads, truth_hap1, truth_hap2, maxDifferenceBetweenReads ) {
    ## do not count where both sites missing
    w <- is.na(truth_hap1) & is.na(truth_hap2)
    truth_hap1[w] <- 0.5
    truth_hap2[w] <- 0.5
    eMatRead_truth_t <- calculate_eMatRead_t_vs_two_haplotypes(
        sampleReads,
        hap1 = truth_hap1,
        hap2 = truth_hap2,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        rescale_eMatRead_t = FALSE
    )
    p1 <- eMatRead_truth_t[1, ]
    p2 <- eMatRead_truth_t[2, ]
    weird <- (colSums(is.na(eMatRead_truth_t)) > 0) | (colSums(eMatRead_truth_t == Inf) > 0) | (colSums(eMatRead_truth_t == -Inf) > 0)
    ## if both 0, set to 1 for both
    if (sum(p1 == 0 | p2 == 0) > 0) {
        ##
        w <- (p1 == 0) & (p2 == 0) & (!weird)
        p1[w] <- 1
        p2[w] <- 1
        w <- (p1 == 0) & (p2 !=0) & (!weird)
        p1[w] <- 1e-100
        p2[w][p2[w] < p1[w]] <- 1e-90
        w <- (p1 != 0) & (p2 ==0) & (!weird)
        p2[w] <- 1e-100
        p1[w][p1[w] < p2[w]] <- 1e-90
    }
    ##
    nReads <- length(sampleReads)
    H <- as.integer(runif(nReads) < (p2 / (p1 + p2))) + 1
    H[is.na(H)] <- sample(1:2, sum(is.na(H)), replace = TRUE)
    truth_labels <- H ## at least, one of them!
    fold_diff <- p1 / p2
    fold_diff[weird] <- 1 ## hmm
    fold_diff[fold_diff > 1] <- (p2 / p1)[fold_diff > 1]
    uncertain_truth_labels <- (fold_diff > 0.5) ## |         (p1 < 1e-5 & p2 < 1e-5)
    return(
        list(
            truth_labels = truth_labels,
            uncertain_truth_labels = uncertain_truth_labels
        )
    )
}



impute_using_everything <- function(
    eMatDH_special_matrix_helper,
    eMatDH_special_matrix,
    use_eMatDH_special_symbols,
    H,
    sampleReads,
    rhb_t,
    nSNPs,
    full_alphaHat_t,
    full_betaHat_t,
    full_gamma_t,
    full_gammaSmall_t,
    full_gammaSmall_cols_to_get,
    full_transMatRate_t_H,
    distinctHapsB,
    distinctHapsIE,
    hapMatcher,
    hapMatcherR,
    use_hapMatcherR,
    eMatDH_special_grid_which,
    eMatDH_special_values_list,
    ref_error,
    make_plots,
    outplotprefix,
    have_truth_haplotypes,
    have_truth_genotypes,
    truth_haps,
    truth_gen,
    truth_labels,
    uncertain_truth_labels,
    L_grid,
    L,
    inRegion2,
    cM_grid,
    ancAlleleFreqAll,
    return_good_haps,
    plot_description,
    Knew,
    previously_selected_haplotypes,
    sample_name,
    smooth_cm,
    regionStart,
    regionEnd,
    buffer,
    minGLValue,
    return_dosage = FALSE,
    return_betaHat_t = FALSE,
    return_gamma_t = FALSE,
    K_top_matches = 5,
    heuristic_match_thin = 0.01,
    suppressOutput = 1,
    use_eigen = FALSE
) {

    if (use_hapMatcherR) {
        K <- nrow(hapMatcherR)
    } else {
        K <- nrow(hapMatcher)
    }
    
    dosage <- numeric(nSNPs)
    nGrids <- ncol(distinctHapsB)
    if (return_good_haps) {
        return_gammaSmall_t <- FALSE
        get_best_haps_from_thinned_sites <- TRUE
        best_haps_stuff_list <- as.list(1:sum(full_gammaSmall_cols_to_get >= 0))
    } else {
        ##gammaSmall_t <-  array(0, c(1, 1))
        ##gammaSmall_cols_to_get <- array(0, 1)
        return_gammaSmall_t <- FALSE
        get_best_haps_from_thinned_sites <- FALSE
        best_haps_stuff_list <- list()
    }
    ##
    new_haps <- as.list(1:2)
    if (make_plots | return_gamma_t) {
        fbsoL <- as.list(1:3)
        fbsoL$list_of_gammas <- as.list(1:2)
        fbsoL$hapProbs_t <- array(NA, c(2, nSNPs))
        return_gamma_t <- TRUE
    } else {
        fbsoL <- NULL
    }
    for(i_hap in 1:2) {
        ##
        u <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
        bq <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
        w <- bq != 0
        bq <- bq[w]
        u <- u[w]
        if (length(u) == 0) {
            print_message(paste0("Read label assignment includes no reads for haplotype ", i_hap))
        }
        gl <- make_gl_from_u_bq(u, bq, nSNPs, minGLValue = minGLValue)
        use_eMatDH <- TRUE
        c <-  array(1, c(nGrids))  ## more useful for debugging
        if (use_eigen) {
            arma_alphaHat_t <- array(0, c(1, 1))
            eigen_alphaHat_t <- full_alphaHat_t
        } else {
            arma_alphaHat_t <- full_alphaHat_t
            eigen_alphaHat_t <- array(0, c(1, 1))
        }
        Rcpp_haploid_dosage_versus_refs(
            gl = gl,
            arma_alphaHat_t = arma_alphaHat_t,
            eigen_alphaHat_t = eigen_alphaHat_t,
            betaHat_t = full_betaHat_t,
            c = c,
            gamma_t = full_gamma_t,
            dosage = dosage,
            transMatRate_t = full_transMatRate_t_H,
            rhb_t = rhb_t,
            ref_error = ref_error,
            use_eMatDH = use_eMatDH,
            distinctHapsB = distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            use_hapMatcherR = use_hapMatcherR,
            eMatDH_special_grid_which = eMatDH_special_grid_which,
            eMatDH_special_values_list = eMatDH_special_values_list,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            eMatDH_special_matrix = eMatDH_special_matrix,
            use_eMatDH_special_symbols = use_eMatDH_special_symbols,
            suppressOutput = suppressOutput,
            return_dosage = return_dosage,
            return_betaHat_t = return_betaHat_t,
            return_gamma_t = return_gamma_t,
            return_gammaSmall_t = return_gammaSmall_t,
            gammaSmall_t = full_gammaSmall_t,
            gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
            get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
            best_haps_stuff_list = best_haps_stuff_list,
            K_top_matches = K_top_matches,
            always_normalize = FALSE,
            is_version_2 = !use_eigen, 
            normalize_emissions = TRUE,
            use_eigen = use_eigen
        )
        ## do some checks here
        if ((min(dosage) < -1e-5) | (1 + 1e-5) < max(dosage)) {
            print(range(dosage))
            stop("Dosage observed outside of range of 0 to 1 on forward-backward full iteration. Something has gone wrong. Please report this")
        }
        dosageNew <- numeric(nSNPs)
        dosageNew <- dosage[]
        if (i_hap == 1) { dosage1 <- dosageNew}
        if (i_hap == 2) { dosage2 <- dosageNew}
        ##
        ## if (return_good_haps) {
        ##     new_haps <- everything_per_hap_prepare_haps(
        ##         full_gammaSmall_t = full_gammaSmall_t,
        ##         K_top_matches = K_top_matches,
        ##         new_haps = new_haps,
        ##         i_hap = i_hap
        ##     )
        ## }
        if (return_good_haps) {
            new_haps[[i_hap]] <- everything_per_hap_rejig_haps(best_haps_stuff_list)
        }
        if ( return_gamma_t | make_plots) {
            ## requires gamma
            fbsoL$hapProbs_t[i_hap, ] <- dosageNew
            gamma_tNew <- array(0, c(K, nGrids))
            gamma_tNew[] <- full_gamma_t
            if (i_hap == 1) { fbsoL$gammaMT_t <- gamma_tNew}
            if (i_hap == 2) { fbsoL$gammaMU_t <- gamma_tNew}
        }
    }

    if (return_good_haps) {
        ## select some new haplotypes if possible!
        ## try all the depth-up-to-5 ones first, then go on
        new_haps <- everything_select_good_haps(
            Knew = Knew,
            K_top_matches = K_top_matches,
            new_haps = new_haps,
            previously_selected_haplotypes = previously_selected_haplotypes,
            K = K
        )
    }
    to_return <- list(
        dosage1 = dosage1,
        dosage2 = dosage2,
        new_haps = new_haps,
        fbsoL = fbsoL
    )

    if (make_plots) {
        plot_single_gamma_dosage(
            sampleReads = sampleReads,
            fbsoL = fbsoL,
            L_grid = L_grid,
            L = L,
            cM_grid = cM_grid,
            inRegion2 = inRegion2,
            ancAlleleFreqAll = ancAlleleFreqAll,
            outname = paste0(outplotprefix, plot_description, ".png"),
            method = "gibbs-nipt",
            haps = truth_haps,
            truth_gen = truth_gen,
            truth_labels = truth_labels,
            have_truth_haplotypes = have_truth_haplotypes,
            have_truth_genotypes = have_truth_genotypes,
            uncertain_truth_labels = uncertain_truth_labels,
            sample_name = sample_name,
            smooth_cm = smooth_cm,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer
        )
    }

    if (1 == 0) {
        check_accuracy_from_all(cheat_all)
    }
    return(to_return)
}



## do much less this time
everything_per_hap_rejig_haps <- function(
    best_haps_stuff_list
) {
    new_haps <- lapply(best_haps_stuff_list, function(x) {
        top_matches <- x[["top_matches"]] + 1 ## turn from 0 to 1 based for R stuff
        top_matches_values <- x[["top_matches_values"]]
        return(top_matches[order(-top_matches_values)])
    })
    return(new_haps)
}


##
## todo, make this work within the forward backward bit, then can post-process in R more quickly
##
everything_per_hap_prepare_haps <- function(
    full_gammaSmall_t,
    K_top_matches,
    new_haps,
    i_hap
) {
    ## OK, so this is slow, the quantile bit argh
    ## how can I make this fast?
    ## can I write c++ to do two passes over this
    ## otherwise
    n <- ncol(full_gammaSmall_t)
    to_out <- as.list(1:n)
    for(i in 1:n) {
        ## for these samples, return them with this
        x <- full_gammaSmall_t[, i]
        y <- R_get_top_K_or_more_matches(x, K_top_matches)
        to_out[[i]] <- y
    }
    ## now choose
    new_haps[[i_hap]] <- to_out
    return(new_haps)
}

R_get_top_K_or_more_matches <- function(x, K_top_matches) {
    val <- -rcpp_nth_partial_sort(-x, as.integer(K_top_matches))[K_top_matches]
    y <- which(x >= val) ## 1-based, surprisingly slow
    y <- y[order(-x[y])]
    return(y)
}

R_get_top_K_or_more_matches_while_building_gamma <- function(
    alphaHat_t,
    betaHat_t_col,
    gamma_t_col,
    iGrid,
    K,
    K_top_matches
) {
    ##
    top_K_values <- numeric(K_top_matches)
    beats_value <- 0;
    num_top_K <- 1
    for(k in 0:(K - 1)) {
        gamma_t_col[k + 1] <- alphaHat_t[k + 1, iGrid + 1] * betaHat_t_col[k + 1];
        if (gamma_t_col[k + 1] == top_K_values[0 + 1]) {
            num_top_K <- num_top_K + 1;
        } else if (gamma_t_col[k + 1] > top_K_values[0 + 1]) {
            num_top_K <- num_top_K + 1;
            beats_value <- 0;
            for(j in 0:(K_top_matches - 1)) {
                if (gamma_t_col[k + 1] > top_K_values[j + 1]) {
                    beats_value <- j
                }
            }
            if (beats_value > 0) {
                for(i in 0:(beats_value)) {
                    top_K_values[i + 1] = top_K_values[i + 1 + 1]
                }
            }
            top_K_values[beats_value + 1] <- gamma_t_col[k + 1]
        }
    }
    ##
    top_matches <- integer(num_top_K)
    top_matches_values <- numeric(num_top_K)
    count <- 0
    for(k in 0:(K - 1)) {
        if (gamma_t_col[k + 1] >= top_K_values[0 + 1]) {
            top_matches[count + 1] <- k + 1
            top_matches_values[count + 1] <- gamma_t_col[k + 1]
            count <- count + 1;
        }
    }
    ## note, can just leave for R in this way
    return(
        list(
            top_matches = top_matches[1:count] - 1, ## make 0-based
            top_matches_values = top_matches_values[1:count]
        )
    )
    ##top_matches[order(top_matches_values , decreasing = TRUE)])
}



everything_select_good_haps <- function(
    Knew,
    K_top_matches,
    new_haps,
    previously_selected_haplotypes,
    K
) {
    i <- 1
    to_keep <- NULL
    left <- Knew
    done <- FALSE
    while(!done) {
        if (i <= K_top_matches) {
            new <- unique(unlist(sapply(new_haps, function(x) lapply(x, function(y) y[i]))))
        } else {
            ## here you just take all because exhausted
            new <- unique(unlist(new_haps))
            done <- TRUE
        }
        new <- setdiff(new, previously_selected_haplotypes)
        new <- setdiff(new, to_keep)
        ## if (length(new) == 0) {
        ##     ## should not happen, have run out! then add some at random
        ##     print("AA")
        ##     toadd <- Knew - length(to_keep)
        ##     done <- TRUE
        ##     to_keep <- c(to_keep, sample(setdiff(setdiff(1:K, previously_selected_haplotypes), to_keep), toadd))
        if (length(new) < (Knew - length(to_keep))) {
            to_keep <- c(to_keep, new)
            i <- i + 1
        } else {
            toadd <- Knew - length(to_keep)
            to_keep <- c(to_keep, new[sample(1:length(new), toadd)])
            done <- TRUE
        }
    }
    ##     new_haps <- to_keep
    if (length(to_keep) < Knew) {
        ## if not enough, add some at random
        to_keep <- c(to_keep, sample(setdiff(1:K, c(to_keep, previously_selected_haplotypes)), Knew - length(to_keep)))
    }
    if (length(to_keep) != Knew) {
        print(paste0("length(to_keep) = ", length(to_keep)))
        print(paste0("Knew = ", Knew))
        save(Knew, K_top_matches, new_haps, previously_selected_haplotypes, K, file = "~/temp.RData")
        stop("Have returned too many haps, see debug file")
    }
    return(to_keep)
}


impute_one_sample <- function(
    eMatDH_special_matrix_helper = array(0, c(1, 1)),
    eMatDH_special_matrix = array(0,c(1, 1)),
    use_eMatDH_special_symbols = FALSE,
    distinctHapsB = array(0L, c(1, 1)),
    distinctHapsIE = array(0, c(1, 1)),
    hapMatcher = array(0L, c(1, 1)),
    hapMatcherR = array(as.raw(0), c(1, 1)),
    use_hapMatcherR = FALSE,
    rhb_t = array(0L, c(1, 1)),
    ref_error = 0.001,
    nSNPs,
    sampleReads,
    small_eHapsCurrent_tc,
    small_transMatRate_tc_H,
    alphaHat_t1,
    betaHat_t1,
    eMatGrid_t1,
    alphaHat_t2,
    betaHat_t2,
    eMatGrid_t2,
    alphaHat_t3,
    betaHat_t3,
    eMatGrid_t3,
    gammaMT_t_local,
    gammaMU_t_local,
    gammaP_t_local,
    small_alphaMatCurrent_tc,
    small_priorCurrent_m,
    smooth_cm,
    which_haps_to_use,
    n_gibbs_starts = 1,
    small_ref_panel_gibbs_iterations,
    n_gibbs_sample_its = 1,
    double_list_of_starting_read_labels,
    small_ref_panel_block_gibbs_iterations,
    perform_block_gibbs = TRUE,
    make_plots,
    maxDifferenceBetweenReads,
    wif0,
    grid_has_read,
    verbose = FALSE,
    shuffle_bin_radius,
    outplotprefix,
    plot_description,
    ancAlleleFreqAll,
    grid,
    L_grid,
    L,
    inRegion2,
    cM_grid,
    have_truth_haplotypes,
    truth_haps,
    have_truth_genotypes,
    truth_gen,
    truth_labels,
    sample_name,
    regionStart,
    regionEnd,
    buffer,
    uncertain_truth_labels,
    small_ref_panel_skip_equally_likely_reads = FALSE,
    small_ref_panel_equally_likely_reads_update_iterations = c(1,2,3,6,9,15),
    return_p_store = FALSE,
    return_p1 = FALSE,
    return_extra = FALSE,
    return_genProbs = TRUE,
    return_hapProbs = TRUE,
    return_gamma = FALSE,
    return_gibbs_block_output = FALSE,
    return_advanced_gibbs_block_output = FALSE,
    gibbs_initialize_iteratively = FALSE,
    gibbs_initialize_at_first_read = FALSE,
    maxEmissionMatrixDifference = 1e100,
    rescale_eMatRead_t = TRUE,
    rescale_eMatGrid_t = FALSE,
    Jmax = 10000,
    suppressOutput = 1,
    use_smooth_cm_in_block_gibbs = TRUE,
    block_gibbs_quantile_prob = 0.95,
    make_plots_block_gibbs = FALSE,
    use_small_eHapsCurrent_tc = TRUE,
    use_provided_small_eHapsCurrent_tc = FALSE,
    use_sample_is_diploid = FALSE,
    i_it = NA,
    i_gibbs_sample = NA,
    ff0_shard_check_every_pair = FALSE,
    zilong = FALSE
) {
    ##
    K <- length(which_haps_to_use)
    S <- 1
    ## print(paste0("start = ", Sys.time()))
    if (use_small_eHapsCurrent_tc & !use_provided_small_eHapsCurrent_tc) {
        inflate_fhb_t_in_place(
            rhb_t,
            small_eHapsCurrent_tc,
            haps_to_get = which_haps_to_use - 1,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
    }
    if (make_plots_block_gibbs) {
        return_gibbs_block_output <- TRUE
        return_advanced_gibbs_block_output <- TRUE
    }
    ##
    ## ugh, alphaMatCurrent_tc is a CONSTANT
    ##
    ## param_list <- list(
    ##     return_alpha = FALSE,
    ##     return_p_store = FALSE,
    ##     return_extra = FALSE,
    ##     return_genProbs = TRUE,
    ##     return_hapProbs = TRUE,
    ##     return_gamma = TRUE,
    ##     return_gibbs_block_output = FALSE,
    ##     return_advanced_gibbs_block_output = FALSE,
    ##     use_starting_read_labels = FALSE,
    ##     verbose = FALSE,
    ##     run_fb_subset = FALSE
    ## )
    if (use_sample_is_diploid) {
        sample_is_diploid <- TRUE
    }else {
        sample_is_diploid <- FALSE
    }
    param_list <- list(
        return_alpha = FALSE,
        return_extra = FALSE,
        return_genProbs = return_genProbs,
        return_gamma = as.logical(return_gamma | make_plots),
        return_hapProbs = as.logical(return_hapProbs | make_plots),
        return_p_store = return_p_store,
        return_p1 = return_p1,
        return_gibbs_block_output = return_gibbs_block_output,
        return_advanced_gibbs_block_output = return_advanced_gibbs_block_output,
        use_starting_read_labels = TRUE,
        verbose = verbose,
        run_fb_subset = FALSE,
        haploid_gibbs_equal_weighting = TRUE,
        gibbs_initialize_iteratively = gibbs_initialize_iteratively,
        gibbs_initialize_at_first_read = gibbs_initialize_at_first_read,
        use_smooth_cm_in_block_gibbs = use_smooth_cm_in_block_gibbs,
        use_small_eHapsCurrent_tc = use_small_eHapsCurrent_tc,
        sample_is_diploid = sample_is_diploid,
        update_in_place = FALSE,
        do_shard_ff0_block_gibbs = TRUE
    )
    ## use_provided_small_eHapsCurrent_tc = use_provided_small_eHapsCurrent_tc
    ## this should catch hopefully rare underflow problems and re-run the samples
    done_imputing <- FALSE
    n_imputing <- 0
    n_gibbs_burn_in_its <- small_ref_panel_gibbs_iterations
    n_gibbs_full_its <- n_gibbs_burn_in_its + n_gibbs_sample_its
    skip_read_iteration <- rep(FALSE, n_gibbs_full_its)
    if (small_ref_panel_skip_equally_likely_reads) {
        skip_read_iteration <- rep(TRUE, n_gibbs_full_its)
        skip_read_iteration[small_ref_panel_equally_likely_reads_update_iterations] <- FALSE
    }
    while(!done_imputing) {

##         print("SAVE ME")
##         save(
##             sampleReads,
## small_priorCurrent_m,
## small_alphaMatCurrent_tc,
## small_eHapsCurrent_tc,
## small_transMatRate_tc_H,
##             hapMatcher ,
##             hapMatcherR ,
##             use_hapMatcherR,
##             distinctHapsB ,
##             distinctHapsIE,
##             eMatDH_special_matrix_helper,
##             eMatDH_special_matrix ,
##             use_eMatDH_special_symbols,
##             rhb_t ,
##             ref_error, 
##             which_haps_to_use, 
##             maxDifferenceBetweenReads,
##             Jmax,
##             maxEmissionMatrixDifference,
##             grid ,
##             skip_read_iteration,
##             alphaHat_t1,
##             alphaHat_t2,
##             alphaHat_t3,
##             betaHat_t1 ,
##             betaHat_t2 ,
##             betaHat_t3 ,
##             eMatGrid_t1,
##             eMatGrid_t2,
##             eMatGrid_t3,
##             gammaMT_t_local, 
##             gammaMU_t_local ,
##             gammaP_t_local ,
##             suppressOutput,
##             n_gibbs_burn_in_its,
##             n_gibbs_sample_its,
##             n_gibbs_starts,
##             double_list_of_starting_read_labels,
##             perform_block_gibbs,
##             wif0 ,
##             grid_has_read,
##             shuffle_bin_radius,
##             L_grid ,
## small_ref_panel_block_gibbs_iterations,
##             rescale_eMatRead_t,
##             smooth_cm ,
##             param_list,
##             block_gibbs_quantile_prob,
##         ff0_shard_check_every_pair,
## file = "/dev/shm/rwdavies/temp.RData", compress = FALSE)
##         stop("WER")

        ## AM HERE
        ## AFTER THIS TRY LOADING UP
        ## SEE IF I CAN TWEAK WITH A-B TESTING

        
        out <- rcpp_forwardBackwardGibbsNIPT(
            sampleReads = sampleReads,
            priorCurrent_m = small_priorCurrent_m,
            alphaMatCurrent_tc = small_alphaMatCurrent_tc,
            eHapsCurrent_tc = small_eHapsCurrent_tc,
            transMatRate_tc_H = small_transMatRate_tc_H,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            use_hapMatcherR = use_hapMatcherR,
            distinctHapsB =distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            eMatDH_special_matrix = eMatDH_special_matrix,
            use_eMatDH_special_symbols = use_eMatDH_special_symbols,
            rhb_t = rhb_t,
            ref_error = ref_error,
            which_haps_to_use = which_haps_to_use,
            ff = 0,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax_local = Jmax,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            run_fb_grid_offset = FALSE,
            blocks_for_output = array(0, c(1, 1)),
            grid = grid,
            pass_in_alphaBeta = TRUE,
            skip_read_iteration = skip_read_iteration,
            alphaHat_t1 = alphaHat_t1,
            alphaHat_t2 = alphaHat_t2,
            alphaHat_t3 = alphaHat_t3,
            betaHat_t1 = betaHat_t1,
            betaHat_t2 = betaHat_t2,
            betaHat_t3 = betaHat_t3,
            eMatGrid_t1 = eMatGrid_t1,
            eMatGrid_t2 = eMatGrid_t2,
            eMatGrid_t3 = eMatGrid_t3,
            gammaMT_t_local = gammaMT_t_local,
            gammaMU_t_local = gammaMU_t_local,
            gammaP_t_local = gammaP_t_local,
            hapSum_tc = array(0, c(1, 1, 1)),
            snp_start_1_based = -1,
            snp_end_1_based = -1,
            generate_fb_snp_offsets = FALSE,
            suppressOutput = suppressOutput,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            n_gibbs_sample_its = n_gibbs_sample_its,
            n_gibbs_starts = n_gibbs_starts,
            double_list_of_starting_read_labels = double_list_of_starting_read_labels,
            prev_list_of_alphaBetaBlocks = as.list(c(1, 2)),
            i_snp_block_for_alpha_beta = -1,
            do_block_resampling = FALSE, ## turn off for now
            perform_block_gibbs = perform_block_gibbs,
            seed_vector = 0,
            update_hapSum = FALSE, ## do not bother for now
            class_sum_cutoff = 0.06, ## what is this
            record_read_set = TRUE, ## needed for block gibbs
            wif0 = wif0,
            grid_has_read = grid_has_read,
            shuffle_bin_radius = shuffle_bin_radius,
            L_grid = L_grid,
            block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
            rescale_eMatRead_t = rescale_eMatRead_t,
            smooth_cm = smooth_cm,
            param_list = param_list,
            block_gibbs_quantile_prob = block_gibbs_quantile_prob,
            artificial_relabel = -1,
            ff0_shard_check_every_pair = ff0_shard_check_every_pair
        )
        if (out[["underflow_problem"]]) {
            new_maxDifferenceBetweenReads <- max(1, maxDifferenceBetweenReads / 10)
            print_message(paste0("Underflow problem observed for sample ", sample_name, ". Re-setting maxDifferenceBetweenReads from ", maxDifferenceBetweenReads, " (log10 of ", log10(maxDifferenceBetweenReads), ") to ", new_maxDifferenceBetweenReads, " (log10 of ", log10(new_maxDifferenceBetweenReads), "). If this problem persists please reset maxDifferenceBetweenReads or downsampleToCov to lower coverage and/or reduce the impact of individual reads"))
            maxDifferenceBetweenReads <- new_maxDifferenceBetweenReads
            n_imputing <- n_imputing + 1
            if (n_imputing > 10) {
                stop(paste0("There were 5 consecutive underflow problems for sample ", sample_name, ". This is likely due to high coverage (at least locally) and/or reads that are very long. Please try setting downsampleToCov to a lower value and/or reduce maxDifferenceBetweenReads to decrease the probability of this happening again"))
            }
        } else {
            done_imputing <- TRUE
        }
    }
    ##
    genProbs_t <- out$genProbsM_t
    hapProbs_t <- out$happrobs_t
    dosage <- genProbs_t[2, ] + 2 * genProbs_t[3, ]
    out$dosage <- dosage
    ##
    if (make_plots) {
        plot_single_gamma_dosage(
            sampleReads = sampleReads,
            fbsoL = out,
            L_grid = L_grid,
            L = L,
            cM_grid = cM_grid,
            inRegion2 = inRegion2,
            ancAlleleFreqAll = ancAlleleFreqAll,
            outname = paste0(outplotprefix, plot_description, ".png"),
            method = "gibbs-nipt",
            haps = truth_haps,
            truth_labels = truth_labels,
            have_truth_haplotypes = have_truth_haplotypes,
            have_truth_genotypes = have_truth_genotypes,
            truth_gen = truth_gen,
            uncertain_truth_labels = uncertain_truth_labels,
            sample_name = sample_name,
            smooth_cm = smooth_cm,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            which_haps_to_use = which_haps_to_use
        )
        ## would like a plot here
    }
    if (make_plots_block_gibbs) {
        nGrids <- ncol(alphaHat_t1)
        for(n_block_it_to_plot in 1:length(small_ref_panel_block_gibbs_iterations)) {
            outname <- paste0(outplotprefix, "block", plot_description, ".n", n_block_it_to_plot, ".png")
            n_block_it_to_plot = n_block_it_to_plot
            plot_attempt_to_reblock_snps(
                out = out,
                nGrids = nGrids,
                block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
                outname = outname,
                L_grid = L_grid,
                L = L,
                uncertain_truth_labels = uncertain_truth_labels,
                truth_labels = truth_labels,
                have_truth_haplotypes = have_truth_haplotypes,
                sampleReads = sampleReads,
                n_block_it_to_plot = n_block_it_to_plot,
                wif0 = wif0,
                grid = grid
            )
        }
    }
    return(out)
}




get_alleleCount <- function(sampleReads , nSNPs) {
    alleleCount <- array(0,c(nSNPs, 3))
    a=unlist(sapply(sampleReads,function(x) x[[3]]))
    b=unlist(sapply(sampleReads,function(x) x[[4]]))
    bqProbs=convertScaledBQtoProbs(matrix(a,ncol=1))
    ## y is numeric, z = integer counts (can have 0),
    c1 <- increment2N(
        y = as.numeric(bqProbs[,1]),
        z = as.numeric(b),
        yT = as.integer(nrow(bqProbs)),
        xT = as.integer(nSNPs - 1)
    )
    c2 <- increment2N(
        y = as.numeric(bqProbs[,2]),
        z = as.numeric(b),
        yT = as.integer(nrow(bqProbs)),
        xT = as.integer(nSNPs - 1)
    )
    alleleCount[,1]=alleleCount[,1]+c2 # fixed nov 6 2015 - was backward
    alleleCount[,2]=alleleCount[,2]+c1+c2
    return(alleleCount)
}


#' @export
r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE) {
    if (flip) {
        w <- af > 0.5
        af[w] <- 1 - af[w]
        truthG[w] <- 2 - truthG[w]
        testDS[w] <- 2 - testDS[w]
    }
    if (!is.null(which_snps)) {
        af <- af[which_snps]
        truthG <- truthG[which_snps]
        testDS <- testDS[which_snps]
    }
    x <- cut(af, breaks = breaks)
    cors_per_af <- tapply(1:length(x), x, function(w) {
        c(
            n = length(w),
            nA = sum(truthG[w], na.rm = TRUE),
            simple = cor(truthG[w], testDS[w], use = 'pairwise.complete') ** 2,
            norm = cor(truthG[w] - 2 * af[w], testDS[w] - 2 * af[w], use = 'pairwise.complete') ** 2
        )
    })
    cors_per_af <- t(sapply(cors_per_af[!sapply(cors_per_af, is.null)], I))
    return(cors_per_af)
}


rsync_to_rescomp <- function(ext = "png") {
    analysis_date <- basename(getwd()) ## probably
    system(paste0("rsync --progress -av plots/*", ext, " rescompNew:/gpfs3/well/davies/users/dcc832/single_imp_test/", analysis_date, "/plots/ "))
}



calc_metric3 <- function(
    eMatRead_t,
    rare_thresh = 0.99,
    hap_thresh = 0.10,
    read_length_thresh = 0.80
) {
    ## among longest 20% of reads
    ## count as rare those occuring with probability of at least that of 99th percentile
    ## but found in fewer than 5% of samples
    x2 <- apply(eMatRead_t, 2, function(x) {
        rare <- x >= quantile(x, probs = rare_thresh, na.rm = TRUE)
        return(sum(rare))
    })
    is_rare <- x2 <= hap_thresh * nrow(eMatRead_t)
    is_long <- nr >= quantile(nr, probs = read_length_thresh)
    which_reads_to_use <- which(is_rare & is_long)
    ## count the occurences of these things
    a <- unlist(lapply(which_reads_to_use, function(iRead) {
        x <- eMatRead_t[, iRead]
        rare <- x >= quantile(x, probs = rare_thresh)
        return(which(rare))
    }))
    x <- increment2N(
        as.integer(length(a)),
        as.integer(nrow(eMatRead_t)),
        as.integer(rep(1, length(a))),
        as.integer(a)
    )[-1] ## first entry is 0-match
    return(
        list(
            x = x,
            which_reads_to_use = which_reads_to_use
        )
    )
    ## alt, among these long, rare reads
    ## look for reads ful-filling a need
}


calculate_eMatRead_t_usig_rhb_t <- function(
    sampleReads,
    rhb_t,
    rescale_eMatRead_t = TRUE
) {
    nReads <- length(sampleReads)
    K <- nrow(rhb_t)
    eMatRead_t_full  <- array(1, c(K, nReads))
    s <- 0
    ## slow but whatevs, will fix
    n <- 5000
    eMatRead_t_full  <- array(1, c(K, nReads))
    for(i in 1:ceiling(K / n)) {
        print(paste0(i, " / ", ceiling(K / n)))
        s <- n * (i - 1) + 1
        e <- min(n * i, K)
        KL <- e - s + 1
        ehc <- array(0, c(KL, nSNPs, 1))
        ## do it in parts
        ehc[, , 1] <- inflate_fhb_t(
            rhb_t,
            haps_to_get = s:e - 1,
            nSNPs
        )
        eMatRead_t_partial  <- array(1, c(KL, nReads))
        gc(reset = TRUE); gc(reset = TRUE);
        ##
        rcpp_make_eMatRead_t(
            eMatRead_t = eMatRead_t_partial,
            sampleReads = sampleReads,
            eHapsCurrent_tc = ehc,
            s = 0,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax = 10000,
            eMatHapOri_t = array(0, c(1, 1)), ## ugh
            pRgivenH1 = array(0),
            pRgivenH2 = array(0),
            run_pseudo_haploid = FALSE,
            prev = 0,
            suppressOutput = 1,
            prev_section = "text",
            next_section = "text",
            rescale_eMatRead_t = rescale_eMatRead_t
        )
        eMatRead_t_full[s:e, ] <- eMatRead_t_partial
        if ((i %% 10) == 0) {
            gc(reset = TRUE);    gc(reset = TRUE);
        }
    }
    gc(reset = TRUE);    gc(reset = TRUE);
    return(eMatRead_t_full)
}



calculate_eMatRead_t_some_haplotypes <- function(
    sampleReads,
    rhb_t,
    which_haps_to_use,
    maxDifferenceBetweenReads,
    rescale_eMatRead_t = TRUE
) {
    nReads <- length(sampleReads)
    K <- length(which_haps_to_use) ## 1-based
    eMatRead_t  <- array(1, c(K, nReads))
    s <- 0
    ## slow but whatevs, will fix
    n <- 1000
    eMatRead_t  <- array(1, c(K, nReads))
    ehc <- array(0, c(K, nSNPs, 1))
    ## do it in parts
    ehc[, , 1] <- inflate_fhb_t(
        rhb_t,
        haps_to_get = which_haps_to_use - 1,
        nSNPs
    )
    gc(reset = TRUE); gc(reset = TRUE);
    ##
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = ehc,
        s = 0,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = 1000,
        eMatHapOri_t = array(0, c(1, 1)), ## ugh
        pRgivenH1 = array(0),
        pRgivenH2 = array(0),
        run_pseudo_haploid = FALSE,
        prev = 0,
        suppressOutput = 1,
        prev_section = "text",
        next_section = "text",
        rescale_eMatRead_t = rescale_eMatRead_t
    )
    gc(reset = TRUE);    gc(reset = TRUE);
    gc(reset = TRUE);    gc(reset = TRUE);
    return(eMatRead_t)
}



calculate_eMatRead_t_vs_two_haplotypes <- function(
    sampleReads,
    hap1,
    hap2,
    maxDifferenceBetweenReads,
    rescale_eMatRead_t = TRUE
) {
    nReads <- length(sampleReads)
    K <- 2
    eMatRead_t <- array(1, c(K, nReads))
    s <- 0
    ## expand whole thing - yuck - fix this later!
    nSNPs <- length(hap1)
    ehc <- array(0, c(K, nSNPs, 1))
    ehc[1, , 1] <- hap1
    ehc[2, , 1] <- hap2
    ## gc(reset = TRUE); gc(reset = TRUE);
    ##
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = ehc,
        s = s,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = 1000,
        eMatHapOri_t = array(0, c(1, 1)), ## ugh
        pRgivenH1 = array(0),
        pRgivenH2 = array(0),
        run_pseudo_haploid = FALSE,
        prev = 0,
        suppressOutput = 1,
        prev_section = "text",
        next_section = "text",
        rescale_eMatRead_t = rescale_eMatRead_t
    )
    rm(ehc)
    ## gc(reset = TRUE);    gc(reset = TRUE);
    return(eMatRead_t)
}






##
## check accuracy quickly
##
check_accuracy_from_all <- function(all) {
    dosage1 <- all$dosage1
    dosage2 <- all$dosage2
    r2s <- c(
        cor(na12878_hap1, dosage1) ** 2,
        cor(na12878_hap2, dosage1) ** 2,
        cor(na12878_hap1, dosage2) ** 2,
        cor(na12878_hap2, dosage2) ** 2
    )
    testDS <- dosage1 + dosage2
    breaks <- sort(unique(c(
        seq(0, 1, length.out = 11)
    )))
    r2sb <- r2_by_freq(breaks, af, truth_g, testDS)
    r2sc <- r2_by_freq(c(0, 1), af, truth_g, testDS)
    afX <- af
    afX[afX > 0.5] <- (1 - af)[afX > 0.5]
    breaks <- sort(unique(c(
        c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
        c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
        c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
        seq(0.1, 0.5, length.out = 5)
    )))
    r2sm <- r2_by_freq(breaks, afX, truth_g, testDS)
    r2sm2 <- r2_by_freq(breaks, af, truth_g, testDS, flip = TRUE)
    breaks2 <- sort(unique(c(
        c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
        c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
        c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
        seq(0.1, 0.5, length.out = 5),
        1 - seq(0.1, 0.5, length.out = 5),
        1 - c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
        1 - c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
        1 - c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1)
    )))
    afX <- af
    r2sn <- r2_by_freq(breaks2, afX, truth_g, testDS)
    to_return <- list(
        r2s = r2s,
        r2sb = r2sb,
        r2sc = r2sc,
        r2sm = r2sm,
        r2sm2 = r2sm2,
        r2sn = r2sn
    )
    to_return
}






if (1 == 0) {

        ## check them inidividually
    ## check them in aggregate
    f <- function(output) {
        d1 <-
            super_out[[1]][[output]][["out"]][[1]]$dosage +
            super_out[[1]][[output]][["out"]][[2]]$dosage
        d2 <-
            super_out[[2]][[output]][["out"]][[1]]$dosage +
            super_out[[2]][[output]][["out"]][[2]]$dosage
        d3 <-
            super_out[[3]][[output]][["out"]][[1]]$dosage +
            super_out[[3]][[output]][["out"]][[2]]$dosage
        d4 <-
            super_out[[4]][[output]][["out"]][[1]]$dosage +
            super_out[[4]][[output]][["out"]][[2]]$dosage
        breaks <- sort(unique(c(
            seq(0, 0.01, length.out = 6),
            seq(0, 0.05, length.out = 6),
            seq(0, 1, length.out = 11)
        )))
        ## just 1
        cbind(
            t(r2_by_freq(breaks, af, truth_g, d1))[, c("n","simple")],
            t(r2_by_freq(breaks, af, truth_g, (d1 + d2) / 2))[, "simple"],
            t(r2_by_freq(breaks, af, truth_g, (d1 + d2 + d3 + d4) / 4))[, "simple"]
        )
    }

    100 * f("cheat_all")
    100 * f("impute_all1")
    100 * f("impute_all2")
    100 * f("impute_all3")

    100 * f("cheat_all")
    100 * f("impute_all1")
    100 * f("impute_all2")
    100 * f("impute_all3")

    ## however - how well does this do for UNLINKED reads
    ## to some extent, how much do I care? (somewhat?)
    ## can try the above, but initialize differently somehow?
    ## 400 haplotypes at random?
    ## hopefully OK?



}


    get_dosages_from_fbsoL <- function(fbsoL) {
        return(1 * fbsoL$genProbsM_t[2, ] + 2 * fbsoL$genProbsM_t[3, ])
    }






estimate_bq <- function(truth_labels, sampleReads, truth_haps) {
    ## do my own bqsr versus truth
    H <- truth_labels
    i_hap <- 1
    u1 <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
    bq1 <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
    i_hap <- 2
    u2 <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
    bq2 <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
    ##
    bqs <- c(bq1, bq2)
    vals <- c(truth_haps[u1, 1], truth_haps[u2, 2])
    ##
    cut_bqs <- cut(bqs, breaks = c(-Inf, -30, -20, -1, 0, 19, 29, Inf))
    results <- t(sapply(c("(-Inf,-30]", "(-30,-20]", "(-20,-1]", "(0,19]", "(19,29]", "(29, Inf]"), function(bqlevel) {
        w <- cut_bqs == bqlevel
        w0 <- cut_bqs == bqlevel & vals == 0
        w1 <- cut_bqs == bqlevel & vals == 1
        phred_stated <- round(exp(1 / sum(w) * sum(sum(log(abs(bqs[w]))))), 1)
        return(c(sum(w0), sum(w1), phred_stated))
    }))
    colnames(results)[1:3] <- c("N_ref", "N_alt", "phred_stated")
    rownames(results) <- c("(-Inf,-30]", "[-29,-20]", "[-19, 1]", "[1, 19]", "[20, 29]", "[30, Inf)")
    results <- cbind(results, accuracy = c(results[1:3, 1] / rowSums(results[1:3, 1:2]), results[4:6, 2] / rowSums(results[4:6, 1:2])))
    results <- cbind(results, phred_obs = -round(log10(1 - results[, "accuracy"]) * 10, 1))
    results <- results[, colnames(results) != "accuracy"]
    return(results)
}




## force phased haplotypes to agree with genotype posterior
recast_haps <- function(hd1, hd2, gp) {
    gt1 <- round(hd1) + round(hd2)
    max_val <- gp[, 1]
    gt3 <- rep(0, nrow(gp))
    for(i in 2:3) {
        w <- gp[, i] > max_val
        gt3[w] <- i - 1
        max_val[w] <- gp[w, i]
    }
    
    ##
    to_change <- which(gt3 != gt1)
    ## easy one
    hd1[to_change][gt3[to_change] == 0] <- 0
    hd2[to_change][gt3[to_change] == 0] <- 0
    hd1[to_change][gt3[to_change] == 2] <- 1
    hd2[to_change][gt3[to_change] == 2] <- 1
    ## hard one
    a1 <- hd1[to_change][gt3[to_change] == 1]
    a2 <- hd2[to_change][gt3[to_change] == 1]
    hd1[to_change][gt3[to_change] == 1][a1 > a2] <- 1
    hd2[to_change][gt3[to_change] == 1][a1 > a2] <- 0
    hd1[to_change][gt3[to_change] == 1][a1 <= a2] <- 0
    hd2[to_change][gt3[to_change] == 1][a1 <= a2] <- 1
    return(list(hd1 = hd1, hd2 = hd2))
}



final_phasing_accuracy_calculation <- function(
    have_truth_haplotypes,
    have_truth_genotypes,
    truth_haps,
    inRegion2,
    dosage,
    af,
    phasing_haps,
    gen,
    sample_name
) {
    if (have_truth_haplotypes) {    
        w <- (inRegion2)
        g <- truth_haps[inRegion2, 1] + truth_haps[inRegion2, 2]
        r2 <-  round(cor((dosage)[inRegion2] - 2 * af[inRegion2], g - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
        ##
        x <- calculate_pse_and_r2_during_gibbs(inRegion2 = inRegion2, hap1 = phasing_haps[, 1], hap2 = phasing_haps[, 2], truth_haps = truth_haps, af = af, verbose = FALSE)
        print_message(paste0("Final imputation dosage accuracy for sample ", sample_name, ", r2:", r2))
        print_message(paste0("Final phasing accuracy for sample ", sample_name, ", pse:", x["pse"], ", disc(%):", x["disc"], "%"))
    } else if (have_truth_genotypes) {
        r2 <-  round(cor((dosage)[inRegion2] - 2 * af[inRegion2], gen[inRegion2, sample_name] - 2 * af[inRegion2], use = "pairwise.complete.obs") ** 2, 3)
        print_message(paste0("Final imputation dosage accuracy for sample ", sample_name, ", r2:", r2))
    }
}

