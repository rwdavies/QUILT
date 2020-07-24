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
