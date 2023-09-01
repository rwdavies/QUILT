get_1000_genomes_samples <- function(
    pop = "CEU",
    n = 10
) {


    ##
    ## normal samples
    ##
    file <- "1000G_2504_high_coverage.sequence.index"
    if (!file.exists(file)) {
        system("wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index")
    }
    head <- unlist(strsplit(system(paste0("grep '^#ENA_FILE_PATH' ", file), intern = TRUE)[[1]], "\t"))
    length(head)
    index <- read.table(file, sep = "\t")
    colnames(index) <- head


    ##
    ## offspring
    ##
    file2 <- "1000G_698_related_high_coverage.sequence.index"
    if (!file.exists(file2)) {
        system("wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index")
    }
    head <- unlist(strsplit(system(paste0("grep '^#ENA_FILE_PATH' ", file2), intern = TRUE)[[1]], "\t"))
    index2 <- read.table(file2, sep = "\t")
    ## colSums(is.na(index2)) == nrow(index2)
    index2 <- index2[, -12]    ## weird error?
    colnames(index2) <- head

    ##
    ## known ped, weird format
    ##
    file3 <- "20130606_g1k.ped"
    if (!file.exists(file3)) {
        system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
    }
    ped <- read.table(file3, sep = "\t", header = TRUE)

    ## OK find samples where we have kid and parents
    kids <- ped[which((ped[, "Paternal.ID"] != 0) & (ped[, "Maternal.ID"] != 0)), ]
    ## now find which of those have all three sequenced
    ## keep it simple, do the normal way
    i1 <- (kids[, "Individual.ID"] %in% index2[, "SAMPLE_NAME"]) ## | (kids[, "Individual.ID"] %in% index2[, "SAMPLE_NAME"])
    i2 <- (kids[, "Paternal.ID"] %in% index[, "SAMPLE_NAME"]) ##| (kids[, "Paternal.ID"] %in% index2[, "SAMPLE_NAME"])
    i3 <- (kids[, "Maternal.ID"] %in% index[, "SAMPLE_NAME"]) ##| (kids[, "Maternal.ID"] %in% index2[, "SAMPLE_NAME"])        
    kids <- kids[which(i1 & i2 & i3), ]

    ## now - want a dataset, with unique samples, with
    ## pop
    ## sample names: kid, father, mother
    ## links: kid, father, mother
    i1 <- match(kids[, "Individual.ID"], index2[, "SAMPLE_NAME"])
    i2 <- match(kids[, "Paternal.ID"], index[, "SAMPLE_NAME"])
    i3 <- match(kids[, "Maternal.ID"], index[, "SAMPLE_NAME"])
    
    kids <- kids[which(i1 & i2 & i3), ]
    data <- cbind(
        POP = index2[i1, "POPULATION"],
        KID_NAME = index2[i1, "SAMPLE_NAME"],
        KID_LINK = index2[i1, "#ENA_FILE_PATH"],
        FATHER_NAME = index[i2, "SAMPLE_NAME"],
        FATHER_LINK = index[i2, "#ENA_FILE_PATH"],
        MOTHER_NAME = index[i3, "SAMPLE_NAME"],
        MOTHER_LINK = index[i3, "#ENA_FILE_PATH"]
    )

    ## go back, keep appropriate part of ped
    ##keep <- ped[, "Individual.ID"] %in% c(data[, "KID_NAME"], data[, "FATHER_NAME"], data[, "MOTHER_NAME"])
    ##ped <- ped[keep, ]
    ##a <- table(ped[, "Family.ID"])
    ##a[a != 3]
    
    ## also build actual ped bits
    to_out <- NULL
    for(i in 1:nrow(data)) {
        dataL <- data[i, ]
        fid <- rep(paste0("FAM", i), 3)
        iid <- dataL[c("KID_NAME", "FATHER_NAME", "MOTHER_NAME")]
        to_out <- rbind(
            to_out,
            cbind(fid = fid, iid = iid, parent1 = c(iid[2], 0, 0), parent2 = c(iid[3], 0, 0))
        )
    }
    ped <- to_out


    ## now, select
    ## take 10 CEU for now
    w <- which(data[, "POP"] == pop)[1:n]
    data <- data[w, ]

    ## to impute
    samples <- data[, "KID_NAME"]
    ## save(samples, file = paste0(output_date, "samples.RData"))

    ## to genotype (and later to help phase)
    samples_needed_for_phasing <- unlist(c(t(data[, c("KID_NAME", "FATHER_NAME", "MOTHER_NAME")])))
    links_needed_for_phasing  <- unlist(c(t(data[, c("KID_LINK", "FATHER_LINK", "MOTHER_LINK")])))

    ## re-size ped
    ped <- ped[ped[, 2] %in% samples_needed_for_phasing, ]
    
    return(
        list(
            samples = samples,
            data = data,
            ped = ped,
            samples_needed_for_phasing = samples_needed_for_phasing,
            links_needed_for_phasing = links_needed_for_phasing
        )
    )

}
    
                             




load_and_convert_vcf <- function(vcf_file) {
    ## truth_vcf <- paste0("bams/all.HC.chr20.19500000.22500000.vcf.gz")
    truth_vcf <- vcf_file
    truth <- fread(cmd = paste0("gunzip -c ", truth_vcf, " | grep -v '^##'"), data.table = FALSE)
    truthG <- matrix(0L, nrow = nrow(truth), ncol = ncol(truth) - 9)
    for(i in 1:(ncol(truth) - 9)) {
        suppressWarnings(
            truthG[, i]  <- as.integer(substr(truth[, i + 9], 1, 1)) + as.integer(substr(truth[, i + 9], 3, 3))
        )
    }
    colnames(truthG) <- colnames(truth)[-c(1:9)]
    rownames(truthG) <- paste0(truth[, 1], truth[, 2], truth[, 4], truth[, 5])    
    return(
        list(
            pos = truth[, c(1, 2, 4, 5)],
            G = truthG
        )
    )
}
