make_and_write_output_file <- function(
    output_filename,
    sampleNames,
    nSNPs,
    N,
    pos,
    nCores,
    complete_set_of_results,
    inRegion2,
    tempdir,
    regionName,
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
    ## only keep the annotation of SNPs in the output region
    ## the per-sample columns were written to disk (already subset to inRegion2)
    ## by each sample after imputation, so here they are pasted together in
    ## blocks of SNPs, rather than holding all N x nSNPs strings in memory
    annot <- list(
        pos[inRegion2, 1],
        format(pos[inRegion2, 2], scientific = FALSE, trim = TRUE),
        ".",
        pos[inRegion2, 3],
        pos[inRegion2, 4],
        ".",
        "PASS",
        INFO[inRegion2],
        FORMAT
    )
    per_sample_files <- file_quilt_per_sample_vcf_col(tempdir, 1:N, regionName)
    paste_vcf_col_files(
        files = per_sample_files,
        output_file = output_unbgzipped,
        nLines = sum(inRegion2),
        annot = annot,
        append = TRUE,
        tempdir = tempdir
    )
    unlink(per_sample_files)

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



## write one sample's column of the VCF (only SNPs in inRegion2) to disk
write_per_sample_vcf_col <- function(per_sample_vcf_col, file, inRegion2) {
    if (length(per_sample_vcf_col) == 1) {
        ## e.g. sample not imputed, all missing
        per_sample_vcf_col <- rep(per_sample_vcf_col, sum(inRegion2))
    } else {
        if (length(per_sample_vcf_col) != length(inRegion2)) {
            stop(paste0("Per-sample VCF column has ", length(per_sample_vcf_col), " entries but there are ", length(inRegion2), " SNPs"))
        }
        per_sample_vcf_col <- per_sample_vcf_col[inRegion2]
    }
    con <- gzfile(file, open = "w", compression = 1)
    on.exit(close(con))
    writeLines(per_sample_vcf_col, con)
    return(NULL)
}



## paste files of per-sample VCF columns together line by line (tab separated),
## optionally preceded by the annotation columns in annot (list of vectors of
## length nLines or 1). Lines are streamed through in C++ without creating R
## strings, to keep memory low and to be fast. If there are more files than
## can be open at once (by default up to 2000, or the system limit on open
## files if lower), they are first pasted together in groups into
## intermediate files
paste_vcf_col_files <- function(
    files,
    output_file,
    nLines,
    annot = NULL,
    append = FALSE,
    tempdir = dirname(output_file),
    max_open = NULL
) {
    if (nLines == 0) {
        if (!append) {
            file.create(output_file)
        }
        return(NULL)
    }
    if (!is.null(annot)) {
        annot_file <- tempfile(pattern = "quilt_vcf_annot.", tmpdir = tempdir, fileext = ".txt")
        on.exit(unlink(annot_file))
        writeLines(do.call(paste, c(annot, sep = "\t")), annot_file)
        files <- c(annot_file, files)
    }
    if (is.null(max_open)) {
        ## leave some room for other open files
        max_open <- max(
            Rcpp_get_max_open_files(min(length(files), 2000) + 100) - 100,
            10
        )
    }
    if (length(files) > max_open) {
        groups <- split(files, ceiling(seq_along(files) / max_open))
        group_files <- tempfile(pattern = paste0("quilt_vcf_cols.", seq_along(groups), "."), tmpdir = tempdir, fileext = ".txt.gz")
        on.exit(unlink(group_files), add = TRUE)
        for(i in seq_along(groups)) {
            Rcpp_paste_line_files(
                files = groups[[i]],
                output_file = group_files[i],
                nLines = nLines,
                append = FALSE,
                gz_output = TRUE
            )
        }
        files <- group_files
    }
    Rcpp_paste_line_files(
        files = files,
        output_file = output_file,
        nLines = nLines,
        append = append,
        gz_output = grepl(".gz$", output_file)
    )
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
        '##FORMAT=<ID=FGP,Number=3,Type=Float,Description="Fetal Posterior genotype probability of 0/0, 0/1, and 1/1">\n',
        '##FORMAT=<ID=FDS,Number=1,Type=Float,Description="Fetal Diploid dosage">\n'
    )
    header2 <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(sampleNames, collapse = "\t", sep="\t"), sep="\t")
    cat(header, header2, "\n", sep="", file = output_vcf_header_file)
    return(NULL)
}
