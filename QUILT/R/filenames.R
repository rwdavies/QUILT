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