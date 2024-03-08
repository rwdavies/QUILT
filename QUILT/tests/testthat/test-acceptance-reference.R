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

## n_snps <- 500
n_snps <- 100
chr <- 10
K <- 6
set.seed(919)
phasemaster <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
reads_span_n_snps <- 3
## want about 4X here
n_reads <- round(4 * n_snps / reads_span_n_snps)
data_package <- STITCH::make_acceptance_test_data_package(
    reads_span_n_snps = reads_span_n_snps,
    n_samples = 3,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 2,
    chr = chr,
    K = K,
    phasemaster = phasemaster
)

## 
L <- 1:n_snps
refs <- rep("A", n_snps)
alts <- rep("G", n_snps)
## make not a SNP
refs[round(n_snps / 10)] <- "AG"
alts[round(n_snps / 10) + 3] <- "AG"
## make not bi-allelic
alts[round(n_snps / 10) + 6] <- "A,G"
## make same position
i <- round(n_snps / 10) + 6
L[i] <- L[i - 1]
L[i + 3] <- L[i + 2]
L[i + 4] <- L[i + 3]

## clean means all SNPs ready to impute
refpack_clean <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 100,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster
)

## dirty means some SNPs should be skipped
refpack_dirty <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 100,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster,
    L = L,
    refs = refs,
    alts = alts
)

refpack <- refpack_clean

test_that("QUILT can remove samples from the reference", {

    set.seed(0105)

    regionStart <- 11
    regionEnd <- 80
    buffer <- 5
    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
    i <- 1
    use_reference_vcf <- TRUE
    use_zilong <- TRUE
    selection_method <- "sample"

    for(use_reference_vcf in c(FALSE, TRUE)) {

        if (use_reference_vcf) {
            refpack <- refpack_dirty
            reference_vcf_file <- refpack$reference_vcf_file
            reference_haplotype_file <- ""
            reference_legend_file <- ""
            use_zilong_options <- c(FALSE, TRUE)
        } else {
            refpack <- refpack_clean
            reference_vcf_file <- ""
            reference_haplotype_file <- refpack$reference_haplotype_file
            reference_legend_file <- refpack$reference_legend_file
            use_zilong_options <- c(FALSE)
        }

        ## make file with samples to exclude
        reference_samples <- refpack$reference_samples[, 1]
        refpack$reference_exclude_samplelist_file <- tempfile()
        refpack$reference_samples_to_exclude <- reference_samples[sample(length(reference_samples), 10)]
        refpack$reference_samples_to_keep <- setdiff(reference_samples, refpack$reference_samples_to_exclude)
        write.table(
            matrix(refpack$reference_samples_to_exclude, ncol = 1),
            file = refpack$reference_exclude_samplelist_file,
            row.names = FALSE,
            col.names = FALSE,
            sep = "",
            quote = FALSE
        )
        rm(reference_samples)

        for(use_zilong in use_zilong_options) {
            
            for(selection_method in c("sample", "pop")) {

                ## print(paste0("use_reference_vcf = ", use_reference_vcf, ", use_zilong = ", use_zilong, ", selection_method = ", selection_method))
                
                if (selection_method == "sample") {
                    reference_exclude_samplelist_file <- refpack$reference_exclude_samplelist_file
                    reference_populations <- NA
                } else {
                    reference_exclude_samplelist_file <- ""
                    reference_populations <- "CEU"
                }
                
                outputdir <- STITCH::make_unique_tempdir()

                ## print(paste0("use_reference_vcf = ", use_reference_vcf, ", selection_method = ", selection_method))
                
                QUILT_prepare_reference(
                    outputdir = outputdir,
                    chr = data_package$chr,
                    regionStart = regionStart,
                    regionEnd = regionEnd,
                    buffer = buffer,
                    reference_vcf_file = reference_vcf_file,
                    reference_haplotype_file = reference_haplotype_file,
                    reference_legend_file = reference_legend_file,
                    reference_sample_file = refpack$reference_sample_file,
                    reference_exclude_samplelist_file = reference_exclude_samplelist_file,
                    reference_populations = reference_populations,
                    nGen = 100
                )
                
                ## check prepared stuff
                load(file_quilt_prepared_reference(outputdir, regionName))

                if (selection_method == "sample") {
                    expect_equal(length(refpack$reference_samples_to_keep) * 2, nrow(hapMatcherR))
                    expect_equal(reference_samples[, 1], refpack$reference_samples_to_keep)
                } else {
                    w <- refpack$reference_samples[, 2] == "CEU"                
                    expect_equal(reference_samples[, "POP"], rep("CEU", sum(w)))
                    expect_equal(nrow(hapMatcherR), sum(w) * 2)
                    expect_equal(reference_samples[, 1], refpack$reference_samples[w, 1])
                }

            }
        }
    }

})



test_that("QUILT can use reference vcf that is dirty with impute_rare_common and use_zilong", {

    set.seed(0105)

    regionStart <- 11
    regionEnd <- 80
    buffer <- 5
    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)

    outputdir <- STITCH::make_unique_tempdir()
    
    QUILT(
        outputdir = outputdir,
        chr = data_package$chr,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        reference_vcf_file = refpack_dirty$reference_vcf_file,
        nGen = 100,
        impute_rare_common = TRUE,
        bamlist = data_package$bamlist        
    )

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
    
})

