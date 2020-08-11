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
