##
## not for general use
## use to prepare example package to demonstrate usage
## to be run from rescomp in single_imp results directory on well
## e.g. /well/davies/users/user-name/single_imp/ && R -f ~/proj/QUILT/scripts/prepare_example.R
## 


package_output_date <- "2020_08_25" ## how to name this output
main_run_date <- "2020_08_25" ## what to base the run on


##
## specify region and paths
##
options(scipen= 999)
out_dir <- paste0("package_", package_output_date, "/")
dir.create(out_dir)
in_ref_panel_dir <- file.path(main_run_date, "ref_panels")
chr <- "chr20"
regionStart <- 2000001
regionEnd <- 4000000
buffer <- 500000
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
