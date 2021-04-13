validate_should_be_integer <- function(x, name) {
    if (is.na(x)) {
        stop(paste0(name, " must be an integer greater than 0"))
    }
    if ((class(x) != "integer") & (class(x) != "numeric")) {
        stop(paste0(name, " must be of class numeric or integer"))
    }
    if (round(x) != x) {
        stop(paste0(name, " must be an integer greater than 0 but you have input:", x))
    }
    if (x < 1) {
        stop(paste0(name,  " must be an integer greater than 0 but you have input:", x))
    }
}

validate_minimum_number_of_sample_reads <- function(minimum_number_of_sample_reads) {
    if (is.na(minimum_number_of_sample_reads)) {
        stop(paste0("minimum_number_of_sample_reads must be an integer greater than 0"))
    }
    if (round(minimum_number_of_sample_reads) != minimum_number_of_sample_reads) {
        stop(paste0("minimum_number_of_sample_reads must be an integer greater than 0 but you have input:", minimum_number_of_sample_reads))
    }
    return(NULL)
}

validate_nMaxDH <- function(nMaxDH) {
    if (is.na(nMaxDH)) {
        stop(paste0("nMaxDH must be an integer greater than 0"))
    }
    if (round(nMaxDH) != nMaxDH) {
        stop(paste0("nMaxDH must be an integer greater than 0 but you have input:", nMaxDH))
    }
    if (nMaxDH < 1) {
        stop(paste0("nMaxDH must be an integer greater than 0 but you have input:", nMaxDH))
    }
    return(NULL)
}


validate_panel_size <- function(panel_size) {
    if (is.na(panel_size)) {
        return(NULL)
    }
    if (!(class(panel_size) %in% c("integer", "numeric"))) {
        stop(paste0("Panel size should be a number but you have input:", panel_size, " with class ", class(panel_size)))
    }
    if (round(panel_size) != panel_size) {
        stop(paste0("Panel size should be an integer value"))
    }
    if (panel_size < 2) {
        stop("Panel size should be at least 2")
    }
    return(NULL)
}

validate_quilt_use_of_region_variables <- function(
    regionStart,
    new_regionStart,
    regionEnd,
    new_regionEnd,
    buffer,
    new_buffer
) {
    ## copy and paste these
    good <- (is.na(regionStart) & is.na(new_regionStart)) | (!is.na(regionStart) & (regionStart == new_regionStart))
    if (!good) {
        stop(paste0("From quilt-prepare-reference, you have selected regionStart = ", regionStart, " but here you have selected regionStart = ", new_regionStart, ". These need to be consistent"))
    }
    ## 
    good <- (is.na(regionEnd) & is.na(new_regionEnd)) | (!is.na(regionEnd) & (regionEnd == new_regionEnd))
    if (!good) {
        stop(paste0("From quilt-prepare-reference, you have selected regionEnd = ", regionEnd, " but here you have selected regionEnd = ", new_regionEnd, ". These need to be consistent"))
    }
    ##
    good <- (is.na(buffer) & is.na(new_buffer)) | (!is.na(buffer) & (buffer == new_buffer))
    if (!good) {
        stop(paste0("From quilt-prepare-reference, you have selected buffer = ", buffer, " but here you have selected buffer = ", new_buffer, ". These need to be consistent"))
    }
    return(NULL)
}


validate_niterations_and_block_gibbs <- function(
    block_gibbs_iterations,
    n_gibbs_burn_in_its
) {
    validate_should_be_integer(n_gibbs_burn_in_its, "n_gibbs_burn_in_its")
    ## need integers, with all entries between 1 and (n_gibbs_burn_in_its - 1)
    ## can have length 0 or be NULL
    if (length(block_gibbs_iterations) > 0) {
        for(i in 1:length(block_gibbs_iterations)) {
            validate_should_be_integer(
                block_gibbs_iterations[i],
                paste0("block_gibbs_iterations entry ", i)
            )
            if ((block_gibbs_iterations[i]) > (n_gibbs_burn_in_its - 1)) {
                stop(paste0("block_gibbs_iteration entries must be between 1 and n_gibbs_burn_in_its - 1 but in entry ", i, " you have selected ", block_gibbs_iterations[i]))
            }
        }
    }
    NULL
}
