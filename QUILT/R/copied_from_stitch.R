print_message <- function(x, include_mem = FALSE) {
    if (include_mem) {
        ## mem <- system("ps auxww | grep 'scripts/profile.R' | grep slave | grep -v 'grep' | awk -v OFS='\t' '$1=$1' | cut -f6", intern = TRUE)
        a <- system("ps aux | grep '/usr/lib64/R/bin/exec/R -f scratch/ukbb_gel.R'", intern = TRUE)
        a <- a[-grep("grep", a)]
        b <- unlist(strsplit(a, " "))
        mem <- as.integer(b[b!= ""][6])
        if (length(mem) > 0) {
            mem <- paste0(paste0(round(as.integer(mem) / 2 ** 20, 3), collapse = ", "), " - ")
            ## mem <- paste0(round(as.integer(mem) / 2 ** 20, 3), collapse = ", ")
        } else {
            mem <- ""
        }
        
    } else {
        mem <- ""
    }
    
    print(
        paste0(
            "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", mem, x
        )
    )
}




check_mclapply_OK <- function(out, stop_message = "An error occured during QUILT. The first such error is above") {
    te <- (sapply(out, class) == "try-error") | sapply(out, is.null)
    if (sum(te) > 0) {
        print_message(out[[which(te)[1]]]) # print first error
        stop(stop_message)
    }
    return(NULL)
}

check_system_OK <- function(out, stop_message = "An error occured during QUILT. The first such error is above") {
    status <- attr(out, "status")
    if (length(status) > 0) {
        if (status > 0) {
            stop(stop_message)
        }
    }
    NULL
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




## get_object_sizes(ls())
get_object_sizes <- function(x) {
    sort( sapply(x,function(x){object.size(get(x))}), decreasing = TRUE) / 1024 / 1024 / 1024
}
