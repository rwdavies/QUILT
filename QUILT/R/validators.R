validate_should_be_integer <- function(x, name, min_val = 1) {
    if (is.na(x)) {
        stop(paste0(name, " must be an integer greater than 0"))
    }
    if ((class(x) != "integer") & (class(x) != "numeric")) {
        stop(paste0(name, " must be of class numeric or integer"))
    }
    if (round(x) != x) {
        stop(paste0(name, " must be an integer greater than 0 but you have input:", x))
    }
    if (x < min_val) {
        stop(paste0(name,  " must be an integer greater than or equal to ", min_val, " but you have input:", x))
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


validate_niterations_and_small_ref_panel_block_gibbs <- function(
    small_ref_panel_block_gibbs_iterations,
    small_ref_panel_gibbs_iterations
) {
    validate_should_be_integer(small_ref_panel_gibbs_iterations, "small_ref_panel_gibbs_iterations")
    ## need integers, with all entries between 1 and (small_ref_panel_gibbs_iterations - 1)
    ## can have length 0 or be NULL
    if (length(small_ref_panel_block_gibbs_iterations) > 0) {
        for(i in 1:length(small_ref_panel_block_gibbs_iterations)) {
            validate_should_be_integer(
                small_ref_panel_block_gibbs_iterations[i],
                paste0("small_ref_panel_block_gibbs_iterations entry ", i)
            )
            if ((small_ref_panel_block_gibbs_iterations[i]) > (small_ref_panel_gibbs_iterations - 1)) {
                stop(paste0("small_ref_panel_block_gibbs_iteration entries must be between 1 and small_ref_panel_gibbs_iterations - 1 but in entry ", i, " you have selected ", small_ref_panel_block_gibbs_iterations[i]))
            }
        }
    }
    NULL
}


validate_n_seek_its_and_n_burn_in_seek_its <- function(
    n_seek_its,
    n_burn_in_seek_its
) {
    validate_should_be_integer(n_seek_its, "n_seek_its")
    validate_should_be_integer(n_burn_in_seek_its, "n_burn_in_seek_its", min_val = 0)
    if (n_burn_in_seek_its >= n_seek_its) {
        stop(paste0("n_burn_in_seek_its = ", n_burn_in_seek_its, " should be stricly less than n_seek_its = ", n_seek_its))
    }
    NULL
}
