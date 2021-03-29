prepare_example_using_1000G <- function(
    package_output_date = "2021_01_14",
    chr = "chr20",
    regionStart = 2000001,
    regionEnd = 2100000,
    buffer = 10000
) {

    ##
    ## specify region and paths
    ##
    options(scipen= 999)
    out_dir <- paste0("package_", package_output_date, "/")
    dir.create(out_dir)
    regionName <- paste0(chr, ".", regionStart, ".", regionEnd)

    ##
    ## get bi-allelic 1000G VCF at region of interest
    ##
    vcf_path <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr20_GRCh38.genotypes.20170504.vcf.gz"
    local_file <- file.path(
        out_dir,
        paste0(
            "ALL.chr20_GRCh38.genotypes.20170504.",
            regionName, ".",
            "vcf.gz"
        )
    )
    cmd <- paste0(
        "bcftools view ",
        "--output-file ", local_file, ".temp ",
        "--output-type z ",
        "--min-alleles 2 --max-alleles 2 --types snps ",
        vcf_path, " ",
        gsub("chr", "", chr), ":", regionStart - buffer, "-", regionEnd + buffer
    )
    system(cmd)
    system(paste0("tabix -f ", local_file, ".temp"))
    
    ## argh - need to fix chromosome labelling here    
    chr_file <- tempfile()
    cat("20\tchr20\n", file = chr_file)
    cmd <- paste0(
        "bcftools annotate ",
        "--rename-chrs ", chr_file, " ",
        "--output ", local_file, " ",
        "--output-type z ",
        local_file, ".temp"
    )
    system(cmd)
    system(paste0("tabix -f ", local_file))




    ##
    ## extract and convert reference panel with no NA12878
    ##
    ref_prefix <- file.path(
        out_dir,
        paste0(
            "ALL.chr20_GRCh38.genotypes.20170504.",
            regionName
        )
    )
    for(i in 1:2) {
        if (i == 1) {
            suffix <- ".noNA12878.vcf.gz "
            samples <- "^NA12878"
        } else {
            suffix <- ".onlyNA12878.vcf.gz "
            samples <- "NA12878"
        }
        cmd <- paste0(
            "bcftools view ",
            "--output-file ", ref_prefix, suffix, 
            "--output-type z ",
            "--samples ", samples, " ", 
            local_file, " ",
            chr, ":", regionStart - buffer, "-", regionEnd + buffer
        )
        system(cmd)
        system(paste0("tabix -f ", ref_prefix, suffix))
    }

    
    ##
    ## convert ref panel to proper input format
    ##
    cmd <- paste0(
        "bcftools convert", " ",
        "--haplegendsample ",
        ref_prefix, ".noNA12878 ", 
        ref_prefix, ".noNA12878.vcf.gz"
    )
    system(cmd)


    ##
    ## genetic map
    ##
    system(paste0("rsync -av /well/davies/shared/recomb/CEU/CEU-chr20-final.b38.txt.gz ", out_dir))


    ##
    ## subset ONT, HT, NYGC files
    ##
    for(i_depth in 1:2) {
        if (i_depth == 1) {
            depth <- "1.0"
        } else {
            depth <- "0.1"
        }
        for(i in 1:3) {
            if (i == 1) {
                inbam <- "/well/davies/shared/haplotagged_v2/GM-only.sorted.reheadered.bam"
                outbam <- paste0(out_dir, "NA12878.haplotagged.", depth, ".bam")
                s <- 1 / 3.4
            } else if (i == 2) {
                inbam <- "/well/davies/shared/ont_bowden/na12878.GRCh38.sorted.bam"
                outbam <- paste0(out_dir, "NA12878.ont.", depth, ".bam")
                s <- 1 / 7.06          
            } else {
                inbam <- "/well/davies/shared/1000G/NA12878.final.bam"
                outbam <- paste0(out_dir, "NA12878.illumina.", depth, ".bam")
                s <- 1 / 32.9
            }
            s <- s * as.numeric(depth)
            system(paste0(
                "samtools view ",
                "-s ", s, " ", 
                "-o ", outbam, " ", 
                inbam, " ",
                chr, ":", regionStart - buffer, "-", regionEnd + buffer
            ))
            system(paste0("samtools index ", outbam))
            if (i == 1) {
                ## re-header, change name
                t <- tempfile()
                system(paste0("samtools view -H ", outbam, " > ", t))
                system(paste0("sed 's/SM:NA12878/SM:NA12878HT/g' ", t, " > ", t, "2"))
                system(paste0("samtools reheader ", t, "2 ", outbam, " > ", outbam, ".temp"))
                system(paste0("mv ", outbam, ".temp ", outbam))
                system(paste0("samtools index -@ 4 ", outbam))
            }
        }
        ##
        ## make bamlist
        ##
        bamlist <- file.path(out_dir, paste0("bamlist.", depth, ".txt"))
        write.table(
            matrix(c(
                paste0(out_dir, "NA12878.haplotagged.", depth, ".bam"),
                paste0(out_dir, "NA12878.ont.", depth, ".bam"),
                paste0(out_dir, "NA12878.illumina.", depth, ".bam")
            ), ncol = 1),
            file = bamlist,
            row.names = FALSE,
            col.names = FALSE,
            sep = "",
            quote = FALSE
        )
    }


    
    ##
    ## make posfile and phasefile
    ##
    vcf <- read.table(paste0(ref_prefix, ".onlyNA12878.vcf.gz"))
    write.table(
        vcf[, c(1, 2, 4, 5)],
        file = paste0(ref_prefix, ".posfile.txt"),
        col.names = FALSE,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    temp <- vcf[, c(10, 10, 10)]
    colnames(temp) <- c("NA12878HT", "NA12878ONT", "NA12878")
    write.table(
        temp,
        file = paste0(ref_prefix, ".phasefile.txt"),
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )

   
        
    
    ##
    ## build ref panel
    ##
    outputdir <- paste0("package_test_", package_output_date)
    QUILT_prepare_reference(
        outputdir = outputdir,
        chr = chr,
        nGen = 100,
        reference_haplotype_file = paste0(ref_prefix, ".noNA12878.hap.gz"),
        reference_legend_file = paste0(ref_prefix, ".noNA12878.legend.gz"),
        genetic_map_file = file.path(out_dir, "CEU-chr20-final.b38.txt.gz"),
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )


    ##
    ## check imputes OK
    ## 
    QUILT(
        outputdir = outputdir,
        chr = chr,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        bamlist = file.path(out_dir, "bamlist.1.0.txt"),
        posfile = paste0(ref_prefix, ".posfile.txt"),        
        phasefile = paste0(ref_prefix, ".phasefile.txt"),        
        bqFilter = 10,
        nCores = 1,
        make_plots = TRUE
    )

    ##  --- haplotagged
    ## [2021-01-14 20:26:04] Final imputation dosage accuracy for sample NA12878HT, r2:0.973
    ## [2021-01-14 20:26:04] Final phasing accuracy for sample NA12878HT, pse:0.2, disc(%):2.4%
    ##  --- ont
    ## [2021-01-14 20:26:25] Final imputation dosage accuracy for sample NA12878ONT, r2:0.926
    ## [2021-01-14 20:26:25] Final phasing accuracy for sample NA12878ONT, pse:0.2, disc(%):15.6%
    ##  --- nygc
    ## [2021-01-14 20:26:54] Final imputation dosage accuracy for sample NA12878, r2:0.971
    ## [2021-01-14 20:26:54] Final phasing accuracy for sample NA12878, pse:0.2, disc(%):1.5%
    
    ##
    ## make it a tarball
    ##
    system(paste0("tar -czvf QUILT_example_", package_output_date, ".tgz package_", package_output_date, "/"))
    ##system(paste0("rsync -av --progress QUILT_example_", package_output_date, ".tgz smew:~/pub_html/"))
    ##system(paste0("ssh smew 'chmod 755 ~/pub_html/QUILT_example_", package_output_date, ".tgz'"))

}


prepare_example_using_hrc <- function(
    package_output_date = "2020_08_25",
    main_run_date = "2020_08_25",
    chr = "chr20",
    regionStart = 2000001,
    regionEnd = 4000000,
    buffer = 500000
) {

    ##
    ## specify region and paths
    ##
    options(scipen= 999)
    out_dir <- paste0("package_", package_output_date, "/")
    dir.create(out_dir)
    in_ref_panel_dir <- file.path(main_run_date, "ref_panels")
    regionName <- paste0(chr, ".", regionStart, ".", regionEnd)


    ##
    ## make smaller .hap .legend .sample file
    ##
    library("data.table")
    inhap <- file.path(in_ref_panel_dir, paste0("hrc.", chr, ".hap.gz"))
    inlegend <- file.path(in_ref_panel_dir, paste0("hrc.", chr, ".legend.gz"))
    insample <- file.path(in_ref_panel_dir, paste0("hrc.", chr, ".samples"))
    outhap <- file.path(out_dir, paste0("ref.", regionName, ".hap.clean.gz"))
    outhap2 <- file.path(out_dir, paste0("ref.", regionName, ".hap.clean.example.gz"))
    outsample <- file.path(out_dir, paste0("ref.", chr, ".forRegion.samples"))
    outlegend <- file.path(out_dir, paste0("ref.", regionName, ".legend.clean.example.gz"))
    phasefile <- file.path(out_dir, paste0("phase.", regionName, ".txt"))
    posfile <- file.path(out_dir, paste0("pos.", regionName, ".txt"))
    outhap2_1000G <- file.path(out_dir, paste0("ref.", regionName, ".hap.clean.example.1000Gonly.gz"))
    outsample_1000G <- file.path(out_dir, paste0("ref.", chr, ".forRegion.1000Gonly.samples"))
    ##
    ##
    ##
    hrc_legend <- fread(cmd = paste0("gunzip -c ", inlegend, ""), data.table = FALSE)
    colnames(hrc_legend) <- c("CHR", "POS", "REF", "ALT")
    stopifnot(length(unique(nchar(hrc_legend[, "REF"]))) == 1)
    stopifnot(length(unique(nchar(hrc_legend[, "ALT"]))) == 1)
    ## want unique and in region
    keep <- which(((regionStart - buffer) <= hrc_legend[, "POS"]) & (hrc_legend[, "POS"] <= (regionEnd + buffer)))
    keep <- match(unique(hrc_legend[keep, "POS"]), hrc_legend[, "POS"])

    ## keep <- match(unique(hrc_legend[matches, "POS"]), hrc_legend[, "POS"])

    for(i in 1:2) {
        if (i == 1) {
            keepX <- matrix(c(1, 1 + keep), ncol = 1)    
            infile <- inlegend
            outfile <- outlegend
        } else if (i == 2) {
            keepX <- matrix(keep, ncol = 1)    
            infile <- inhap
            outfile <- outhap
        }
        message(paste0("Work on infile:", infile))
        message(paste0("Write to outfile:", outfile))
        message(date())
        ## 
        file1 <- tempfile()
        write.table(
            keepX,
            file = file1,
            row.names = FALSE,
            col.names = FALSE,
            sep = "",
            quote = FALSE
        )
        ## now rebuild
        cmd <- paste0(
            "awk '{if(NR == FNR) {a[$1]=$1; next} if (FNR in a) {print $0}} ' ",
            file1, " ",
            "<(gunzip -c ", infile, ") | sed 's/\t/ /g' | gzip -1 > ",
            outfile
        )
        file2 <- tempfile()
        cat(cmd, file = file2)
        system(paste0("bash ", file2))
    }


    ##
    ## after the fact, blank .hap file entries for non-1000G entries, also NA12878 (as opposed to removing)
    ##
    samples <- read.table(insample, header = TRUE, stringsAsFactors = FALSE)
    hap <- fread(cmd = paste0("gunzip -c ", outhap), data.table = FALSE)
    ##
    ## write posfile and phasefile for NA12878 from the bams NA12878 and NA12878ONT
    ##
    a <- 2 * rep(which(samples[, "SAMPLE"] == "NA12878"), 2)
    a[1] <- a[1] - 1
    phase <- cbind(paste0(hap[, a[1]], "|", hap[, a[2]]), paste0(hap[, a[1]], "|", hap[, a[2]]), paste0(hap[, a[1]], "|", hap[, a[2]]))
    colnames(phase) <- c("NA12878ONT", "NA12878HT_tagged", "NA12878HT_untagged")
    write.table(phase, file = phasefile, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    legend <- read.table(outlegend, header =TRUE, stringsAsFactors = FALSE)
    write.table(cbind(chr, legend[, 2], legend[, 3], legend[, 4]), file = posfile, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    hap[, a] <- 0
    n <- 2175
    ## argh, blank out the rest of them for release. not sure why not working matrix-wise
    for(i in (2*n + 1):ncol(hap)) {
        if ((i %% 100) == 0) {print(paste0(i, " / ", ncol(hap), ", ", date()))}
        hap[, i] <- 0
    }

    ##
    ## make full version output
    ##
    w <- (1:n)
    samples[-w, ] <- cbind(
        paste0("SAMPLE", 1:(nrow(samples) - n)),
        "REDACTED",
        "REDACTED",
        2
    )
    write.table(
        samples,
        file = outsample,
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    fwrite(hap, file = gsub(".gz", "", outhap2), quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
    system(paste0("gzip -f -1 ", gsub(".gz", "", outhap2)))

    ##
    ## also make version with only 1000 Genomes entries for faster run time
    ##
    samples <- samples[1:n, ]
    write.table(
        samples,
        file = outsample_1000G,
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    fwrite(hap[, 1:(2 * n)], file = gsub(".gz", "", outhap2_1000G), quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
    system(paste0("gzip -f -1 ", gsub(".gz", "", outhap2_1000G)))


    ##
    ## now, remove subsetted hap file
    ##
    unlink(outhap)


    ##
    ## genetic map
    ##
    system(paste0("rsync -av /well/davies/shared/recomb/CEU/CEU-chr20-final.b38.txt.gz ", out_dir))


    ##
    ##  make bamlist
    ##  get haplotagged bam keep tags, haplotagged bam remove tags, and ONT bam
    ##
    bamdir <- file.path(main_run_date, "bams")
    ##
    dir.create(file.path(out_dir, "bams"))
    for(cov in c("0.25", "1.0")) {
        s <- regionStart - 2 * buffer ## make bigger
        e <- regionEnd + 2 * buffer
        r <- paste0(chr, ":", s, "-", e)
        o <- NULL
        for(i in 1:3) {
            ## in order
            ## 1 = haplotagged, keep bx tag
            ## 2 = haplotagged, remove bx tag
            ## 3 = ONT
            if (i == 1) { f1 <- paste0("NA12878.NA12878HT.cov.", cov, ".chr20.downTagged.bam"); f2 <- "NA12878HT_tagged"; sm <- "SM:NA12878"}
            if (i == 2) { f1 <- paste0("NA12878.NA12878HT.cov.", cov, ".chr20.downTagged.bam"); f2 <- "NA12878HT_untagged"; sm <- "SM:NA12878"}
            if (i == 3) { f1 <- paste0("NA12878.NA12878ONT.cov.", cov, ".chr20.downOnly.bam"); f2 <- "NA12878ONT"; sm <- "SM:NA12878ONT"}        
            outbam <- paste0(out_dir, "/bams/", f2, ".cov.", cov, ".chr20.", s, ".", e, ".bam")
            if (i == 2) {
                additional_command <- "| sed 's/\tBX:Z:*//' "
            } else {
                additional_command <- ""
            }
            command <- paste0("samtools view -h ", bamdir, "/", f1, " ", r, " ", additional_command, " | samtools view -b > ", outbam)
            system(command)
            o <- c(o, paste0(paste0("package_", package_output_date), "/bams/", f2, ".cov.", cov, ".chr20.", s, ".", e, ".bam"))
            ## re-name in header
            t <- tempfile()
            system(paste0("samtools view -H ", outbam, " > ", t))
            system(paste0("sed 's/", sm, "/SM:", f2, "/g' ", t, " > ", t, "2"))
            system(paste0("samtools reheader ", t, "2 ", outbam, " > ", outbam, ".temp"))
            ##
            system(paste0("mv ", outbam, ".temp ", outbam))
            system(paste0("samtools index -@ 4 ", outbam))
        }
        write.table(matrix(o, ncol = 1), file = file.path(out_dir, paste0("bamlist.", cov, ".txt")), row.names = FALSE, col.names =FALSE, sep = "\t", quote = FALSE)
    }


    ##
    ## make it a tarball
    ##
    system(paste0("tar -czvf QUILT_example_", package_output_date, ".tgz package_", package_output_date, "/"))
    system(paste0("rsync -av --progress QUILT_example_", package_output_date, ".tgz smew:~/pub_html/"))
    system(paste0("ssh smew 'chmod 755 ~/pub_html/QUILT_example_", package_output_date, ".tgz'"))

}
