make_and_write_output_file <- function(
    output_filename,
    sampleNames,
    nSNPs,
    N,
    pos,
    nCores,
    complete_set_of_results,
    inRegion2,
    sampleRanges,
    addOptimalHapsToVCF,
    output_gt_phased_genotypes,
    method
) {

    print_message("Begin making and writing output file")

    output_unbgzipped <- gsub(".gz", "", output_filename)

    ## make header
    if (method == "diploid") {
        make_and_write_quilt_header(
            output_vcf_header_file = output_unbgzipped,
            sampleNames = sampleNames,
            addOptimalHapsToVCF = addOptimalHapsToVCF,
            output_gt_phased_genotypes = output_gt_phased_genotypes
        )
    } else {
        make_and_write_quilt_nipt_header(
            output_vcf_header_file = output_unbgzipped,
            sampleNames = sampleNames
        )
    }

    ##
    ## re-build final annotation stuff across all cores
    ##
    hweCount <- array(0, c(nSNPs, 3))
    infoCount <- array(0, c(nSNPs, 2))
    afCount <- array(0, nSNPs)
    alleleCount <- array(0, c(nSNPs, 2))    
    for(i in 1:length(complete_set_of_results)) {
        afCount <- afCount + complete_set_of_results[[i]]$afCount
        hweCount <- hweCount + complete_set_of_results[[i]]$hweCount
        infoCount <- infoCount + complete_set_of_results[[i]]$infoCount
        alleleCount <- alleleCount + complete_set_of_results[[i]]$alleleCount
    }
    ## finalize
    alleleCount <- cbind(alleleCount, alleleCount[, 1] / alleleCount[, 2])
    thetaHat <- infoCount[, 1] / 2 / N
    denom <- 2 * N * thetaHat * (1-thetaHat)
    info <- 1 - infoCount[, 2] / denom
    ## block out those where thetaHat is really close to 0 or 1
    ## when very rare
    info[(round(thetaHat, 2) == 0) | (round(thetaHat, 2) == 1)] <- 1
    info[info < 0] <- 0
    estimatedAlleleFrequency <- afCount / N
    hwe <- STITCH::generate_hwe_on_counts(hweCount, nSNPs, nCores)
    
    ##
    ## summarize info bit
    ##
    if (method == "nipt") {
        FORMAT <- "GT:MGP:MDS:FGP:FDS"        
    }  else {
        if (addOptimalHapsToVCF) {
            FORMAT <- "GT:GP:DS:HD:OHD" ## genotypes, posteriors, dosages, haploid-dosages, optimal haploid-dosages
        } else {
            FORMAT <- "GT:GP:DS:HD" ## genotypes, posteriors, dosages, haploid-dosages
        }
    }
    INFO <- paste0(
        "EAF=", round(estimatedAlleleFrequency, 5), ";",
        "INFO_SCORE=", round(info, 5), ";",
        "HWE=", formatC(hwe, format = "e", digits = 2), ";",
        "ERC=", round(alleleCount[, 1], 5), ";",
        "EAC=", round(alleleCount[, 2] - alleleCount[, 1], 5), ";",
        "PAF=", round(alleleCount[, 3], 5)
    )
    vcf_matrix_to_out <- data.frame(matrix(
        data = NA,
        nrow = nSNPs,
        ncol = N + 9
    ))
    vcf_matrix_to_out[, 1] <- pos[, 1]
    vcf_matrix_to_out[, 2] <- pos[, 2]
    vcf_matrix_to_out[, 3] <- "."
    vcf_matrix_to_out[, 4] <- pos[, 3]
    vcf_matrix_to_out[, 5] <- pos[, 4]
    vcf_matrix_to_out[, 6] <- "."
    vcf_matrix_to_out[, 7] <- "PASS"
    vcf_matrix_to_out[, 8] <- INFO
    vcf_matrix_to_out[, 9] <- FORMAT

    ## should work
    for(i_core in 1:length(sampleRanges)) {
        sampleRange <- sampleRanges[[i_core]]
        sampleRange <- sampleRange[1]:sampleRange[2]
        for(ii in 1:length(sampleRange)) {
            iSample <- sampleRange[ii]
            vcf_matrix_to_out[, iSample + 9] <- complete_set_of_results[[i_core]][["results_across_samples"]][[ii]][["per_sample_vcf_col"]]
        }
    }

    ## ugh - inefficient
    vcf_matrix_to_out <- vcf_matrix_to_out[inRegion2, , drop = FALSE]
    
    data.table::fwrite(
        vcf_matrix_to_out,
        file = output_unbgzipped,
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE,
        append = TRUE,
        nThread = nCores
    )

    print_message("bgzip output file")
    if (length(grep("~", output_unbgzipped)) > 0) {
        ## not entirely sure why this system call isn't working otherwise
        check_system_OK(system(paste0("bgzip --threads ", nCores, " -f ", output_unbgzipped), intern = TRUE))
        print(system(paste0("tabix -f ", output_filename), intern = TRUE))
        check_system_OK(system(paste0("tabix -f ", output_filename), intern = TRUE))
    } else {
        check_system_OK(system(paste0("bgzip --threads ", nCores, " -f ", shQuote(output_unbgzipped)), intern = TRUE))
        check_system_OK(system(paste0("tabix -f ", shQuote(output_filename)), intern = TRUE))
    }
    print_message("Done making and writing output file")

    return(NULL)

}





get_per_snp_annot <- function(
    output_format,
    reference_panel_SNPs,                              
    estimatedAlleleFrequency = NULL,
    info = NULL,
    hwe = NULL,
    alleleCount = NULL,
    return_annotation_only = FALSE    
) {
    INFO <- NULL
    if (output_format == "bgvcf") {
        if (!return_annotation_only) {
            INFO <- paste0(
                "EAF=", round(estimatedAlleleFrequency, 5), ";",
                "INFO_SCORE=", round(info, 5), ";",
                "HWE=", formatC(hwe, format = "e", digits = 2), ";",
                "ERC=", round(alleleCount[, 1], 5), ";",
                "EAC=", round(alleleCount[, 2] - alleleCount[, 1], 5), ";",
                "PAF=", round(alleleCount[, 3], 5)
            )
        }
        annot_header <- paste0(
            '##INFO=<ID=INFO_SCORE,Number=.,Type=Float,Description="Info score">\n',
            '##INFO=<ID=EAF,Number=.,Type=Float,Description="Estimated allele frequency">\n',
            '##INFO=<ID=HWE,Number=.,Type=Float,Description="Hardy-Weinberg p-value">\n',
            '##INFO=<ID=ERC,Number=.,Type=Float,Description="Estimated number of copies of the reference allele from the pileup">\n',
            '##INFO=<ID=EAC,Number=.,Type=Float,Description="Estimated number of copies of the alternate allele from the pileup">\n',
            '##INFO=<ID=PAF,Number=.,Type=Float,Description="Estimated allele frequency using the pileup of reference and alternate alleles">\n'
        )
        if (sum(reference_panel_SNPs) > 0) {
            if (!return_annotation_only) {
                INFO <- paste0(INFO, ";REF_PANEL=", as.integer(reference_panel_SNPs))
            }
            annot_header <- paste0(
                annot_header,
                '##INFO=<ID=REF_PANEL,Number=.,Type=Integer,Description="Whether a SNP was (1) or was not (0) found in the reference panel during imputation">\n'
            )
        }
    } else {
        if (!return_annotation_only) {
            INFO <- data.frame(
                "EAF" =  round(estimatedAlleleFrequency, 5),
                "INFO_SCORE" = round(info, 5),
                "HWE" = formatC(hwe, format = "e", digits = 2),
                "ERC" = round(alleleCount[, 1], 5),
                "EAC" = round(alleleCount[, 2] - alleleCount[, 1], 5),
                "PAF" = round(alleleCount[, 3], 5)
            )
        }
        annot_header <- c(
            '#INFO_SCORE Imputation info score, same as IMPUTE info measure I_A',
            '#EAF Estimated allele frequency after imputation',
            '#HWE Hardy-Weinberg p-value',
            '#ERC Estimated number of copies of the reference allele from the pileup',
            '#EAC Estimated number of copies of the alternate allele from the pileup',
            '#PAF Estimated allele frequency using the pileup of reference and alternate alleles'
        )
        if (sum(reference_panel_SNPs) > 0) {
            if (!return_annotation_only) {
                INFO$REF_PANEL <- as.integer(reference_panel_SNPs)
            }
            annot_header <- c(
                annot_header,
                '#REF_PANEL Whether a SNP was (1) or was not (0) found in the reference panel during imputation'
            )
        }
    }
    return(
        list(
            INFO = INFO,
            annot_header = annot_header
        )
    )
}


make_and_write_quilt_header <- function(
    output_vcf_header_file,
    sampleNames,
    addOptimalHapsToVCF,
    output_gt_phased_genotypes
) {
    ## annotations now constant
    annot_header <- paste0(
        '##INFO=<ID=INFO_SCORE,Number=.,Type=Float,Description="Info score">\n',
        '##INFO=<ID=EAF,Number=.,Type=Float,Description="Estimated allele frequency">\n',
        '##INFO=<ID=HWE,Number=.,Type=Float,Description="Hardy-Weinberg p-value">\n',
        '##INFO=<ID=ERC,Number=.,Type=Float,Description="Estimated number of copies of the reference allele from the pileup">\n',
        '##INFO=<ID=EAC,Number=.,Type=Float,Description="Estimated number of copies of the alternate allele from the pileup">\n',
        '##INFO=<ID=PAF,Number=.,Type=Float,Description="Estimated allele frequency using the pileup of reference and alternate alleles">\n'
    )
    ##
    if (output_gt_phased_genotypes) {
        gt_annot <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">\n'
    } else {
        gt_annot <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Most likely genotype, given posterior probability of at least 0.90">\n'
    }
    header <- paste0(
        '##fileformat=VCFv4.0\n',
        annot_header,        
        gt_annot,
        '##FORMAT=<ID=GP,Number=3,Type=Float,Description="Posterior genotype probability of 0/0, 0/1, and 1/1">\n',
        '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Diploid dosage">\n',
        '##FORMAT=<ID=HD,Number=2,Type=Float,Description="Haploid dosages">\n'
    )
    if (addOptimalHapsToVCF) {
        header <- paste0(header, '##FORMAT=<ID=OHD,Number=2,Type=Float,Description="Optimal haploid dosages">\n')
    }
    header2 <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(sampleNames, collapse = "\t", sep="\t"), sep="\t")
    cat(header, header2, "\n", sep="", file = output_vcf_header_file)
    return(NULL)
}

make_and_write_quilt_nipt_header <- function(
    output_vcf_header_file,
    sampleNames
) {
    ## annotations now constant
    annot_header <- paste0(
        '##INFO=<ID=INFO_SCORE,Number=.,Type=Float,Description="Info score from maternal genotype posteriors">\n',
        '##INFO=<ID=EAF,Number=.,Type=Float,Description="Estimated allele frequency">\n',
        '##INFO=<ID=HWE,Number=.,Type=Float,Description="Hardy-Weinberg p-value from maternal genotypes">\n',
        '##INFO=<ID=ERC,Number=.,Type=Float,Description="Estimated number of copies of the reference allele from the pileup">\n',
        '##INFO=<ID=EAC,Number=.,Type=Float,Description="Estimated number of copies of the alternate allele from the pileup">\n',
        '##INFO=<ID=PAF,Number=.,Type=Float,Description="Estimated allele frequency using the pileup of reference and alternate alleles">\n'
    )
    ##
    gt_annot <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes in order of maternal transmitted, maternal untransmitted, and fetal transmitted">\n'
    header <- paste0(
        '##fileformat=VCFv4.0\n',
        annot_header,        
        gt_annot,
        '##FORMAT=<ID=MGP,Number=3,Type=Float,Description="Maternal Posterior genotype probability of 0/0, 0/1, and 1/1">\n',
        '##FORMAT=<ID=MDS,Number=1,Type=Float,Description="Maternal Diploid dosage">\n',
        '##FORMAT=<ID=FGP,Number=3,Type=Float,Description="Maternal Posterior genotype probability of 0/0, 0/1, and 1/1">\n',
        '##FORMAT=<ID=FDS,Number=1,Type=Float,Description="Maternal Diploid dosage">\n'
    )
    header2 <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(sampleNames, collapse = "\t", sep="\t"), sep="\t")
    cat(header, header2, "\n", sep="", file = output_vcf_header_file)
    return(NULL)
}
