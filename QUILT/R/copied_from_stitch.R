print_message <- function(x, include_mem = FALSE) {
    if (include_mem) {
        mem <- system("ps auxww | grep 'scripts/profile.R' | grep slave | grep -v 'grep' | awk -v OFS='\t' '$1=$1' | cut -f6", intern = TRUE)
        if (length(mem) > 0) {
            mem <- paste0(paste0(round(as.integer(mem) / 2 ** 20, 3), ollapse = ", "), " - ")
        } else {
            mem <- ""
        }
    } else {
        mem <- ""
    }
    message(
        paste0(
            "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", mem, x
        )
    )
}
