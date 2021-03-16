file_quilt_prepared_reference <- function(outputdir, regionName) {
    output_file <- file.path(
        outputdir, "RData", paste0("QUILT_prepared_reference.", regionName, ".RData")
    )
    return(output_file)
}



file_quilt_output_RData <- function(outputdir, regionName) {
    output_file <- file.path(
        outputdir, "RData", paste0("QUILT_output.", regionName, ".RData")
    )
    return(output_file)
}


file_quilt_hla_panelfile <- function(outputdir, hlaregion) {
    quilt_panelfile <- file.path(
        outputdir,
        paste0("quilt.hla.", hlaregion, ".haplotypes.RData")
    )
    return(quilt_panelfile)
}

file_quilt_hla_all_alleles_kmers <- function(outputdir) {
    quilt_hla_all_alleles_kmers <- file.path(
        outputdir,
        "HLAallalleleskmers.out"
    )
    return(quilt_hla_all_alleles_kmers)
}

file_quilt_hla_full_alleles_filled_in <- function(outputdir, hla_region) {
    quilt_hla_full_alleles_filled_in <- file.path(
        outputdir,
        paste0("HLA", hla_region, "fullallelesfilledin.out")
    )
    return(quilt_hla_full_alleles_filled_in)
}

file_quilt_hla_full <- function(outputdir, hla_region) {
    quilt_hla_full <- file.path(
        outputdir,
        paste0("hla", hla_region, "full.out")
    )
    return(quilt_hla_full)
}


file_quilt_hla_snpformatalleles <- function(outputdir, hla_region) {
    return(file.path(
        outputdir,
        paste0("hla", hla_region, "snpformatalleles.out")
    ))
}


file_quilt_final_RData_output_file <- function(outputdir, region) {
    suffix <- paste0("quilt.output.hla", region, ".RData")
    if (outputdir == "") {
        return(suffix)
    } else {
        return(file.path(outputdir, suffix))
    }
}

file_quilt_hla_phase_step_1 <- function(outputdir, region) {
    file.path(
        outputdir,
        paste("hla", region, "newphased.out")
    )
    ## paste("hla",region,"newphased.out",sep="")
}


file_quilt_hla_phased_haplotypes <- function(outputdir, region) {
    file.path(
        outputdir,
        paste("hla", region, "haptypes.out")
    )
}
