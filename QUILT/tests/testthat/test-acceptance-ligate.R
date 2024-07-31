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



test_that("can concatenate using external software", {

  
    n_snps <- 220
    chr <- 10
    K <- 6
    set.seed(919)
    phasemaster <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
    reads_span_n_snps <- 3
    ## want about 1X here
    n_reads <- round(1.0 * n_snps / reads_span_n_snps)
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


    outputdir <- STITCH::make_unique_tempdir()
    i_region <- 1

    to_ligate <- c()
    to_ligate2 <- c()
    to_ligate3 <- c()        
    
    for(i_region in 1:5) {

        regionStartInterior <-  10 + 30 * i_region + 1
        regionStart <- regionStartInterior + - 10
        regionEndInterior <- 40 + 30 * i_region
        regionEnd <- regionEndInterior + 10

        regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
        buffer <- 0

        set.seed(i_region * 2)
        QUILT(
            outputdir = outputdir,
            chr = data_package$chr,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            bamlist = data_package$bamlist,
            nGibbsSamples = 3,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            genetic_map_file = refpack$reference_genetic_map_file,
            nGen = 100,
            phasefile = data_package$phasefile,
            posfile = data_package$posfile
        )

        ## OK so can build new VCF using bcftools hack
        vcffile <- file.path(outputdir, paste0("quilt.", regionName, ".vcf.gz"))
        newfile <- file.path(outputdir, paste0("quilt.", regionName, ".phased.vcf.gz"))        
        system(paste0("bcftools view -h ", shQuote(vcffile), " > ", shQuote(gsub(".gz", "", newfile))))
        ##
        system(paste0("bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\tGT\t[%GT\\t]\\n' ", shQuote(vcffile), " >> ", shQuote(gsub(".gz", "", newfile))))
        system(paste0("bgzip -f ", shQuote(gsub(".gz", "", newfile))))
        system(paste0("tabix -f ", shQuote(newfile)))

        ## either do the fancy one, or just this one
        to_ligate <- c(to_ligate, newfile)
        to_ligate2 <- c(to_ligate2, vcffile)

        ## argh - subset for second version
        newfile2 <- file.path(outputdir, paste0("quilt.", regionName, ".exact.vcf.gz"))
        system(paste0("bcftools view ", shQuote(vcffile), " ", data_package$chr, ":", regionStartInterior, "-", regionEndInterior, " | bgzip > ", shQuote(newfile2)))
        to_ligate3 <- c(to_ligate3, newfile2)
        
    }


    ## yes ligation option
    newfile <- file.path(outputdir, paste0("quilt.phased.vcf.gz"))            
    system(paste0("bcftools concat --output-type z --ligate ", paste0(shQuote(to_ligate2), collapse = " "), " > ", shQuote(newfile))) #, ignore.stderr = TRUE)
    system(paste0("tabix -f ", shQuote(newfile)))

    ## no ligation option
    newfile2 <- file.path(outputdir, paste0("quilt.phasedwrong.vcf.gz"))            
    system(paste0("bcftools concat --output-type z ", paste0(shQuote(to_ligate3), collapse = " "), " > ", shQuote(newfile2)), ignore.stderr = TRUE)
    system(paste0("tabix -f ", shQuote(newfile2)))    

    ## check it is OK! especially phase
    check_orientation <- function(newfile) {
        vcf <- read.table(newfile)
        sapply(1:(ncol(vcf) - 9), function(iSample) {
            ## 
            hap1 <- as.integer(substr(vcf[, 9 + iSample], 1, 1))
            hap2 <- as.integer(substr(vcf[, 9 + iSample], 3, 3))
            truth <- data_package$phase[vcf[, 2], iSample, ]
            ## check agreement
            orient1 <- sum(truth[, 1] != hap1) + sum(truth[, 2] != hap2)
            orient2 <- sum(truth[, 1] != hap2) + sum(truth[, 2] != hap1)
            ## OK - phasing error left!
            (c(orient1, orient2))
        })
    }

    ## check no phasing errors in individual files
    for(file in to_ligate) {
        expect_true(sum(colSums(check_orientation(file) == 0) != 1) == 0)
    }
    
    ## check that each of them has a 0 i.e. is phased across entire length correctly
    expect_true(sum(colSums(check_orientation(newfile) == 0) != 1) == 0)
    ## check that by contrast that WITHOUT ligate this is not true
    expect_false(sum(colSums(check_orientation(newfile2) == 0) != 1) == 0)    
        
    ## check_orientation(newfile2)
    ## OK! phasing error in second file! HMM, why?
    ## ## ## OK looks good
    
    ## ## OK, see what happened here
    ## iSample <- 3
    ## truth <- data_package$phase[vcf[, 2], iSample, ]
    ## hap1 <- as.integer(substr(vcf[, 9 + iSample], 1, 1))
    ## hap2 <- as.integer(substr(vcf[, 9 + iSample], 3, 3))
    
    ## m <- vcf[1:110, c(2, 12)]
    ## m <- cbind(m, NA, NA, NA, NA)
    ## m[1:50, 3:4] <- read.table("quilt.10.31.80.phased.vcf.gz")[, c(2, 12)]
    ## m[31:80, 5:6] <- read.table("quilt.10.61.110.phased.vcf.gz")[, c(2, 12)]
    ## m <- cbind(m, NA, truth[1:110, ])


    ## OK - to confirm
    ## this code checks below that bcftools concat ligate will take the dosage etc from the middle outwards
    ## e.g. if two files with overlap of 20 SNPs (not sure if on SNPs or physical position - probably physical)
    ## then the first half will come from first file, second half from second file
    if (1 == 0) {

        ## phased
        f2 <- function(x) {
            as.numeric(sapply(strsplit(as.character(x[, 2]), ":"), function(x) {
                if (length(x) == 1)
                    return(NA)
                x[[3]]
            }))
        }
        f <- function(x) {
            ## overlap is 60-80
            f2(x[match(51:90, x[, 1]), ])
        }
        cbind(
            f(read.table(newfile)[, c(2, 10)]),
            f(read.table(to_ligate2[1])[, c(2, 10)]),
            f(read.table(to_ligate2[2])[, c(2, 10)])
        )
        
        
    }

    
})



