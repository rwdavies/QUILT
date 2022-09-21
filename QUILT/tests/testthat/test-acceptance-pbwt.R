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

n_snps <- 50
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
refpack <- STITCH::make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 500,
    reference_populations = c("CEU", "GBR"),
    chr = chr,
    phasemaster = phasemaster
)
set.seed(010)


test_that("QUILT can impute a few samples using robbie mspbwt", {

    outputdir <- STITCH::make_unique_tempdir()

    regionStart <- 11
    regionEnd <- 40
    buffer <- 5
    QUILT_prepare_reference(
        outputdir = outputdir,
        chr = data_package$chr,
        nGen = 100,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        genetic_map_file = refpack$reference_genetic_map_file,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        use_mspbwt = TRUE
    )
    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
    expect_true(file.exists(file_quilt_prepared_reference(outputdir, regionName)))
    i <- 1

    ##phasefile = data_package$phasefile,
    phasefile <- ""

    QUILT(
        outputdir = outputdir,
        chr = data_package$chr,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        phasefile = phasefile,
        use_mspbwt = TRUE
    )

    which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)

    ## now evaluate versus truth!
    check_quilt_output(
        file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
        data_package = data_package,
        which_snps = which_snps,
        tol = 0.1,
        min_info = 0.9
    )

})



test_that("QUILT can impute a few samples using zilong pbwt", {

    outputdir <- STITCH::make_unique_tempdir()

    print("Fix me - figure out how best to make this work when reference files are larger than region to be imputed")

    skip("work in progress")
    
    regionStart <- 1
    regionEnd <- 50
    buffer <- 0
    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
    i <- 1
    phasefile <- ""
    QUILT_prepare_reference(
        outputdir = outputdir,
        chr = data_package$chr,
        nGen = 100,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        genetic_map_file = refpack$reference_genetic_map_file,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        use_mspbwt = TRUE
    )


    prefix <- gsub(".txt.gz", "", refpack$reference_haplotype_file)

    x <- read.table(refpack$reference_legend_file, header = TRUE)
    x[, 1] <- paste0(chr, ":", x[, 2], "_", x[, 3], "_", x[, 4])
    new_ref_legend_file <- tempfile(fileext = ".vcf")
    write.table(x, new_ref_legend_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    system(paste0("gzip -1 ", new_ref_legend_file))
    new_ref_legend_file <- paste0(new_ref_legend_file, ".gz")

    fake_ref_vcf <- tempfile(fileext=".vcf.gz")
    cmd <- paste0(
        "bcftools convert ",
        "--haplegendsample2vcf ",
        shQuote(refpack$reference_haplotype_file), ",",
        shQuote(new_ref_legend_file), ",",
        shQuote(refpack$reference_sample_file), " ",
        "--output-type z ",
        "--output ", fake_ref_vcf
    )
    system(cmd)
    system(paste("tabix ", fake_ref_vcf))


    QUILT(
        outputdir = outputdir,
        chr = data_package$chr,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        bamlist = data_package$bamlist,
        zilong = TRUE,
        vcf = fake_ref_vcf,
        reference_sample_file = refpack$reference_sample_file
    )

    which_snps <- (regionStart <= data_package$L) & (data_package$L <= regionEnd)

    ## now evaluate versus truth!
    check_quilt_output(
        file = file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz")),
        data_package = data_package,
        which_snps = which_snps,
        tol = 0.1,
        min_info = 0.9
    )

})



