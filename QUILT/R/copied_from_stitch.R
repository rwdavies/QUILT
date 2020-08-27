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


check_mclapply_OK <- function(out, stop_message = "An error occured during QUILT. The first such error is above") {
    te <- sapply(out, class) == "try-error"
    if (sum(te) > 0) {
        print_message(out[[which(te)[1]]]) # print first error
        stop(stop_message)
    }
    return(NULL)
}


## for this in the future use from STITCH
## am too lazy to export and re-version STITCH right now
quilt_get_chromosome_length <- function(iBam, bam_files, cram_files, chr) {
    if (length(bam_files) > 0) {
        file <- bam_files[iBam]
    } else if (length(cram_files) > 0) {
        file <- cram_files[iBam]
    }
    header <- system(paste0("samtools view -H ", file, " | grep ^@SQ"), intern = TRUE)
    seq_spots <- lapply(header, strsplit, split = "\t")
    if (length(seq_spots) == 0)
        stop(paste0("There is no @SQ tag (with sample name) for:", file))
    seq <- t(sapply(seq_spots, function(x) {
        y <- x[[1]]
        chrName <- substr(y[substr(y, 1, 3) == "SN:"], 4, 1000)
        chrLength <- substr(y[substr(y, 1, 3) == "LN:"], 4, 1000)
        return(c(chrName, chrLength))
    }))
    seq <- seq[seq[, 1] == chr, , drop = FALSE]
    if (nrow(seq) == 0)
        stop(paste0("Could not find chromosome length for file:", file))
    return(seq[1, 2])
}


