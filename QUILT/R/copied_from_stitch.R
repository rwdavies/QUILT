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
    message(
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




## modified from STITCH one for now
## possibly, long term, modify that one directly
validate_phase_header <- function(phasefile) {
    first_row <- as.character(unlist(read.table(phasefile,sep="\t",nrows=1)))
    ## check none are 0|0, 0|1, 1|0, 1|1
    m <- match(first_row, c("0|0", "0|1", "1|0", "1|1"))
    if (sum(is.na(m) == FALSE) > 0) {
        m2 <- which.max(is.na(m) == FALSE)
        stop(paste0(
            "The header for the phasefile is either invalid or missing. ",
            "The first invalid entry is in position ", m2, " and is entry ", first_row[m2]
        ))
    }
}

validate_phase_col <- function(col, i_samp, method = "diploid") {
    x <- sapply(col, length)
    n <- c( diploid = 2, nipt = 3)[method]
    if (sum(x != n) > 0) {
        m <- which.max(x)
        stop(paste0(
            "Unable to split column ", i_samp, " of phasefile at position ", m,
            " with entry ", col[m], " due to lack of field separator |"
        ))
    }
}

get_and_validate_phase <- function(
    phasefile,
    method = "diploid"
) {
    if (phasefile == "") {
        return(NULL)
    }
    phaseX <- read.table(phasefile, header = TRUE, stringsAsFactors = FALSE)
    validate_phase_header(phasefile)
    n <- c( diploid = 2, nipt = 3)[method]    
    phase <- array(0, c(nrow(phaseX), ncol(phaseX), n))
    for(i_samp in 1:ncol(phase)) {
        col <- strsplit(phaseX[, i_samp], "|", fixed = TRUE)
        validate_phase_col(col, i_samp, method)
        for(j in 1:n) {
            x <- sapply(col, function(x) x[j])
            ## only allow NA, 0, or 1
            w <- x %in% c(0, 1, NA, "0", "1", "NA") ## any of these are fine
            if (sum(!w) > 0) {
                ## report error
                i_row <- which.max(!w)
                p2 <- phaseX[i_row, i_samp]
                stop(paste0(
                    "The phasefile contains entries other than 0, 1 or NA. ",
                    "One such entry is in column ", i_samp, " and row ", i_row, " ",
                    " with value ", paste(p2, collapse = "|")
                ))
            }
            phase[, i_samp, j] <- suppressWarnings(as.integer(x))
        }
    }
    if (length(colnames(phaseX)) == 1) {
        dimnames(phase)[[2]] <- list(colnames(phaseX))
    } else {
        dimnames(phase)[[2]] <- colnames(phaseX)
    }
    if (length(dim(phase)) != 3) {
        stop("The phasefile does not have the right number of dimensions")
    }
    return(phase)
}


alpha_col <- function(col, alpha) {
    x <- col2rgb(col) / 255
    return(rgb(x["red", 1], x["green", 1], x["blue", 1], alpha = alpha)    )
}

## get_object_sizes(ls())
get_object_sizes <- function(x) {
    sort( sapply(x,function(x){object.size(get(x))}), decreasing = TRUE) / 1024 / 1024 / 1024
}
