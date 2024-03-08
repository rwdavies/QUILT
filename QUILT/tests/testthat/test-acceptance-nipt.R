if ( 1 == 0 ) {

    library("testthat")
    library("QUILT")
    dir <- "~/proj/QUILT/"
    setwd(paste0(dir, "/QUILT/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)


}


set.seed(1911919)

n_snps <- 400
K <- 6
n_big_haps <- 100 ## 1000
chr <- 10


## make a spectrum or rare and common SNPs
phasemaster1 <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
choice_of_k <- sample(1:K, n_big_haps, replace = TRUE)
phasemaster2 <- phasemaster1[, choice_of_k]

## make about 20% common, and do not touch
n_common <- round(n_snps / 5)
common <- sample(1:n_snps, n_common)
## make about 90% rare, and make about 1% freq
phasemaster2[-common, ] <- 0
## make remaining frequency about 1%
phasemaster2[-common, ] <- sample(c(0, 1), (n_snps - n_common) * n_big_haps, prob = c(0.99, 0.01), replace = TRUE)

reads_span_n_snps <- 3
## want about 4X here
n_reads <- round(4 * n_snps / reads_span_n_snps)
## make mother father
f <- function(x) {
    which.max(choice_of_k == x)
}
haps_to_sample_all <- rbind(
    c(f(1), f(2)),
    c(f(1), f(3)),
    c(f(4), f(5)),
    c(f(4), f(6))
)
data_package <- STITCH::make_acceptance_test_data_package(
    reads_span_n_snps = reads_span_n_snps,
    n_samples = 4,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 2,
    chr = chr,
    K = K,
    phasemaster = phasemaster2,
    haps_to_sample_all = haps_to_sample_all
)

refpack <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 500,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster2
)
set.seed(010)
chr <- data_package$chr




test_that("QUILT can impute a few samples with NIPT, w/wo IRC, w/wo regular vs mspbwt vs zilong", {

    nipt <- TRUE
    impute_rare_common <- FALSE
    i_method <- 1
    
    regionStart <- 11
    regionEnd <- n_snps - 10
    buffer <- 5
    nCores <- 1
    rare_af_threshold <- 0.01

    for(nipt in c(TRUE, FALSE)) {

        if (nipt == TRUE) {
            
            ## just merge them together
            bams <- read.table(data_package$bamlist)
            stopifnot((nrow(bams) %% 2) == 0) ## assume even number
            nSamples <- nrow(bams) / 2

            ## 
            ff <- seq(0.1, 0.3, length.out = nSamples)
            fflist <- tempfile()
            cat(ff, file = fflist, sep = "\n")

            ## make new bams, merge them, in proportion
            newbams <- sapply(1:nSamples, function(iSample) {
                ##
                bam1 <- bams[2 * iSample - 1, 1]
                bam2 <- bams[2 * iSample - 0, 1]
                ## downsample each using samtools
                newsam <- tempfile(fileext = ".sam")
                newbam <- gsub(".sam", ".bam", newsam)
                system(paste0("samtools view -s ", (1 - ff[iSample]), " ", bam1, " > ", newsam))
                system(paste0("samtools view -s ", ff[iSample], " ", bam2, " >> ", newsam))
                system(paste0("samtools view -H ", bam1, " > ", newsam, "2"))
                system(paste0("sort -n -t$'\t' -k4 ", newsam, " >> ", newsam, "2"))
                ## sort
                system(paste0("samtools view -bS ", newsam, "2 > ", newbam))
                system(paste0("samtools index ", newbam))
                unlink(newsam)
                unlink(paste0(newsam, "2"))
                newbam
            })
            bamlist <- tempfile(fileext = ".txt")
            cat(newbams, file = bamlist, sep = "\n")
            
            ## merge together a few bams with some proportions
            method <- "nipt"

            ## write phase file here
            phase <- data_package$phase
            stopifnot(nSamples == 2) ## meh
            ## this relates to how haps_to_sample_all was made above
            m <- cbind(
                paste0(phase[, 1, 1], "|", phase[, 1, 2], "|", phase[, 2, 2]),
                paste0(phase[, 3, 1], "|", phase[, 3, 2], "|", phase[, 4, 2])
            )
            phasefile <- tempfile()
            colnames(m) <- dimnames(phase)[[2]][c(1, 3)] ## yuck but OK
            write.table(
                m,
                file = phasefile,
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t",
                quote = FALSE
            )
            
        } else {

            bamlist <- data_package$bamlist
            phasefile <- data_package$phasefile
            fflist <- ""
            method <- "diploid"
            
        }
        
        for(impute_rare_common in c(TRUE, FALSE)) {
            
            for(i_method in 1:2) {

                if (i_method == 1) {
                    use_mspbwt <- FALSE
                } else if (i_method == 2) {
                    use_mspbwt <- TRUE
                } else if (i_method == 3) {
                    use_mspbwt <- FALSE
                }

                print(paste0("nipt = ", nipt, ", impute_rare_common = ", impute_rare_common, ", i_method = ", i_method, ", ", date()))

                outputdir <- STITCH::make_unique_tempdir()

                QUILT(
                    outputdir = outputdir,
                    chr = chr,
                    regionStart = regionStart,
                    regionEnd = regionEnd,
                    buffer = buffer,
                    bamlist = bamlist,
                    fflist = fflist,
                    reference_vcf_file= refpack$reference_vcf_file,
                    reference_haplotype_file = refpack$reference_haplotype_file,
                    reference_legend_file = refpack$reference_legend_file,
                    genetic_map_file = refpack$reference_genetic_map_file,
                    posfile = data_package$posfile,
                    phasefile = phasefile,
                    nGibbsSamples = 3,
                    n_seek_its = 3,
                    nCores = nCores,
                    nGen = 100,
                    use_mspbwt = use_mspbwt,
                    impute_rare_common = impute_rare_common,
                    mspbwt_nindices = 1,
                    method = method
                )

                ## check it made stuff?
                regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
                which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)
                
                check_quilt_output(
                    file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
                    data_package = data_package,
                    which_snps = which_snps,
                    tol = 0.1,
                    min_info = 0.01,
                    check_info_only = TRUE
                )
                

            }
        }
    }
    


})
