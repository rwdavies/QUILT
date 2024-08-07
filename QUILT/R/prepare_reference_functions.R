make_face_vcf_with_sites_list <- function(outputdir, regionName, output_sites_filename, pos_to_use) {
    print_message("Make VCF with sites list")
    if (is.na(output_sites_filename)) {
        output_sites_filename <- file.path(outputdir, paste0("quilt.sites.", regionName, ".vcf.gz"))
    } else {
        STITCH::validate_output_filename(output_sites_filename, output_format = "bgvcf")
    }
    cat(
        "##fileformat=VCFv4.2",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
        sep = "\n",
        file = gsub(".gz", "", output_sites_filename)
    )
    to_out <- cbind(
        pos_to_use[ , 1],
        pos_to_use[, 2],
        ".",
        pos_to_use[, 3],
        pos_to_use[, 4],
        1000,
        "PASS",
        'DP=1000',
        "GT",
        "0/0"
    )
    write.table(to_out, file = gsub(".gz", "", output_sites_filename), append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    check_program_dependency("bgzip")
    check_program_dependency("tabix")
    system(paste0("bgzip -f ", gsub(".gz", "", output_sites_filename)))
    system(paste0("tabix -f ", output_sites_filename))
    print_message("Done making VCF with sites list")
    NULL
}



remove_sites_from_pos_to_use <- function(region_exclude_file, pos_to_use, chr) {
    if (!file.exists(region_exclude_file)) {
        stop(paste0("Cannot find region_exclude_file:", region_exclude_file))
    }
    print_message("Loading list of regions to exclude")
    regionstoexclude <- read.table(region_exclude_file, header = TRUE, as.is = TRUE)
    regionstoexclude <- regionstoexclude[regionstoexclude[, 2] == chr,, drop = FALSE]
    if (nrow(regionstoexclude) == 0) {
        warning("No regions to exclude based on region_exclude_file. Perhaps the chr isn't the same?")
    } else {
        regionstoexclude <- as.matrix(regionstoexclude)
        excluderegions <- matrix(as.double(regionstoexclude[, 3:4]),ncol = 2)
        rownames(excluderegions) <- regionstoexclude[, 1]
        ## now do the removal
        keep <- array(TRUE, nrow(pos_to_use))
        oldsum <- 0
        for(i in 1:nrow(excluderegions)){
            to_remove <-
                pos_to_use[,"POS"] >= excluderegions[i, 1] &
                pos_to_use[,"POS"] <= excluderegions[i, 2]
            keep[to_remove] <- FALSE
            ##sentence=paste("Excluding ",sum(keep==0)-oldsum," SNPs from region ",rownames(excluderegions)[i])
            ##print(sentence)
            ##oldsum <- sum(keep == 0)
        }
        pos_to_use <- pos_to_use[keep, ]
        print_message(paste0("Excluded ", length(keep) - sum(keep), " out of ", length(keep), " SNPs from ", nrow(excluderegions), " regions"))
    }
    pos_to_use
}



get_pos_to_use_from_reference_legend_file <- function(reference_legend_file, chr, regionStart, regionEnd, buffer) {
    ref_legend <- data.table::fread(cmd = paste0("gunzip -c ", shQuote(reference_legend_file)), data.table = FALSE)
    colnames(ref_legend) <- c("CHR", "POS", "REF", "ALT")
    ref_legend[, 1] <- chr
    if (!is.na(regionStart)) {
        a <- unique(ref_legend[
        (ref_legend[, "POS"] >= (regionStart - buffer)) &
        (ref_legend[, "POS"] <= (regionEnd   + buffer)), "POS"])
        pos_to_use <- ref_legend[match(a, ref_legend[, "POS"]), ]
    } else {
        pos_to_use <- ref_legend
    }
    pos_to_use
}



get_sigmaCurrent_m <- function(nGen, cM_grid, L_grid, expRate, minRate, maxRate, nGrids, S) {
    rateBetweenGrids <- nGen * diff(cM_grid) / 100
    expRateBetweenGrids <- nGen * diff(L_grid) / 1e6 * (expRate / 100)
    minRateBetweenGrids <- nGen * diff(L_grid) / 1e6 * (minRate / 100)
    maxRateBetweenGrids <- nGen * diff(L_grid) / 1e6 * (maxRate / 100)    
    m <- cbind(expRateBetweenGrids, rateBetweenGrids, minRateBetweenGrids)
    m[, 1] <- m[, 1] * diff(L_grid)
    m[, 2] <- m[, 2] * diff(L_grid)
    m[, 3] <- m[, 3] * diff(L_grid)
    m <- cbind(diff(L_grid), m)
    ##
    w <- rateBetweenGrids < minRateBetweenGrids
    rateBetweenGrids[w] <- minRateBetweenGrids[w]
    print_message(paste0("There are ", sum(w), " regions out of ", length(w), " below minimum recombination rate, setting them to minimum rate"))    
    w <- rateBetweenGrids > maxRateBetweenGrids
    rateBetweenGrids[w] <- maxRateBetweenGrids[w]
    print_message(paste0("There are ", sum(w), " regions out of ", length(w), " above maximum recombination rate, setting them to maximum rate"))    
    sigmaCurrent_m <- array(exp(-rateBetweenGrids), c(nGrids - 1, S))
    sigmaCurrent_m
}






rare_common_validate_region_to_impute_when_using_regionStart <- function(
    L,
    regionStart,
    regionEnd,
    buffer,
    snp_is_common
) {
    for (i in 1:3) {
        if (i == 1) {
            w1 <- ((regionStart - buffer) <= L) & (L < regionStart)
            x <- paste0((regionStart - buffer), " <= position < ", regionStart)
        } else if (i == 2) {
            w1 <- (regionStart <= L) & (L <= regionEnd)
            x <- paste0(regionStart, " <= position <= ", regionEnd)
        } else if (i == 3) {
            w1 <- (regionEnd < L) & (L <= (regionEnd + buffer))
            x <- paste0(regionEnd, " < position <= ", regionEnd + buffer)
        }
        s1 <- sum(w1 & snp_is_common)
        s2 <- sum(w1 & !snp_is_common)
        print_message(paste0("There are ", s1, " common and ", s2, " rare (", s1 + s2, " total) variants in the ", 
            c("left buffer", "central", "right buffer")[i], " region ", 
            x))
        if (i == 2) {
            nCentralSNPs <- s1 + s2
        }
    }
    if (length(L) < 2) {
        stop("There are fewer than 2 SNPs to impute, i.e. there is 1 SNP to impute. In this case, imputation is really just genotyping. STITCH could support genotyping but does not, and note that this kind of defeats the point of imputation. Please use your favourite genotyper e.g. GATK to genotype these SNPs. If you strongly disagree please file a bug report and this can be re-examined")
    }
    if (nCentralSNPs < 1) {
        stop("There are insufficient SNPs to impute in the central region")
    }
    return(NULL)
}


get_transMatRate_tc_H_and_smooth_cm <- function(sigmaCurrent_m, L_grid, shuffle_bin_radius) {
    small_transMatRate_tc_H <- get_transMatRate_m("pseudoHaploid", sigmaCurrent_m)
    X <- get_transMatRate_m("pseudoHaploid", sigmaCurrent_m)
    full_transMatRate_t_H <- array(NA, c(dim(X)[1], dim(X)[2]))
    ## ARGH R dropping dimensions
    full_transMatRate_t_H[, ] <- X[, , 1]
    rate2 <- -log(small_transMatRate_tc_H[1, , 1]) * 100
    smooth_cm <- rcpp_make_smoothed_rate(rate2, L_grid, shuffle_bin_radius, verbose = FALSE);
    smooth_cm <- smooth_cm / max(smooth_cm)
    return(
        list(
            smooth_cm = smooth_cm,
            small_transMatRate_tc_H = small_transMatRate_tc_H,
            full_transMatRate_t_H = full_transMatRate_t_H
        )
    )
}



prepare_full_objects_for_rare_common <- function(
    pos_all,
    expRate,
    minRate,
    maxRate,
    nGen,
    shuffle_bin_radius,
    genetic_map,
    rare_per_hap_info,
    snp_is_common,
    Ksubset,
    regionStart,
    regionEnd,
    ref_alleleCount_all
) {

    ## here with respect to ALL SNPs    
    L <- pos_all[, 2]
 
    out2 <- assign_positions_to_grid(
        L = L,
        grid32 = TRUE,
        cM = NULL
    )
    nGrids <- out2$nGrids
    grid <- out2$grid
    L_grid <- as.integer(out2$L_grid)
    dl <- diff(L_grid)

    cM_grid <- match_genetic_map_to_L(
        genetic_map = genetic_map,
        L = L_grid,
        expRate = expRate
    )
    
    sigmaCurrent_m <- get_sigmaCurrent_m(
        nGen,
        cM_grid,
        L_grid,
        expRate,
        minRate,
        maxRate,
        nGrids,
        S = 1
    )
    out <- get_transMatRate_tc_H_and_smooth_cm(
        sigmaCurrent_m = sigmaCurrent_m,
        L_grid = L_grid,
        shuffle_bin_radius = shuffle_bin_radius
    )

    ## yuck I still use these
    S <- 1
    small_priorCurrent_m <- array(1 / Ksubset, c(Ksubset, S))
    small_alphaMatCurrent_tc <- array(1 / Ksubset, c(Ksubset, nGrids - 1, S))

    if (is.na(regionStart)) {
        inRegion2 <- array(TRUE, length(L))
    } else {
        inRegion2 <- (regionStart <= L) & (L <= regionEnd)
    }
    
    special_rare_common_objects <- list(
        rare_per_hap_info = rare_per_hap_info,
        snp_is_common = snp_is_common,
        small_priorCurrent_m = small_priorCurrent_m,
        small_alphaMatCurrent_tc = small_alphaMatCurrent_tc,
        grid = grid,
        L = L,
        L_grid = L_grid,
        inRegion2 = inRegion2,
        cM_grid = cM_grid,
        ref_alleleCount_all = ref_alleleCount_all
    )
    special_rare_common_objects <- append(special_rare_common_objects, out)

    special_rare_common_objects
    
}
