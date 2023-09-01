print(commandArgs(trailingOnly = TRUE))
print(paste0(commandArgs(trailingOnly = TRUE), collapse = "', '"))
print(system("hostname"))
args <- commandArgs(trailingOnly = TRUE)

if (1 == 0) {

    setwd("/well/davies/users/dcc832/single_imp/2020_06_10/")
    args <- c('genotypes/hc.ont.chr6.vcf.gz', 'genotypes/phase.ont.chr6', 'PR33_PR33_PR33_PR22_PR22_PR22_BB15_BB15_BB15_VN068_VN068_VN068_GB52_GB52_GB52_SL16_SL16_SL16_PK43_PK43_PK43', 'HG01109_HG01107_HG01108_HG01243_HG01241_HG01242_HG02055_HG02053_HG02054_HG02080_HG02082_HG02081_HG02723_HG02721_HG02722_HG03098_HG03096_HG03097_HG03492_HG03490_HG03491', 'HG01107_0_0_HG01241_0_0_HG02053_0_0_HG02082_0_0_HG02721_0_0_HG03096_0_0_HG03490_0_0', 'HG01108_0_0_HG01242_0_0_HG02054_0_0_HG02081_0_0_HG02722_0_0_HG03097_0_0_HG03491_0_0', '/well/davies/shared/recomb/CEU/CEU-chr6-final.b38.txt.gz', 'chr6', 'ref_panels/hrc.chr6.clean.bcf', '1')
    args <- c('genotypes/vcf.NA12878.chr6.30000001.32000000.vcf.gz', 'genotypes/phase.NA12878.chr6.30000001.32000000', 'fam1_fam1_fam1', 'NA12878_NA12891_NA12892', 'NA12891_0_0', 'NA12892_0_0', '/well/davies/shared/recomb/CEU/CEU-chr6-final.b38.txt.gz', 'chr6', 'ref_panels/hrc.chr6.clean.bcf', '1')
    setwd("/data/smew1/rdavies/ukbb_gel_2023_01_26")
    args <- c('bams/specific_clean_everyone.HC.20.19000000.31000000.vcf.gz', 'bams/everyone.HC.20.19000000.31000000', 'FAM51_FAM51_FAM51_FAM52_FAM52_FAM52_FAM53_FAM53_FAM53_FAM54_FAM54_FAM54_FAM55_FAM55_FAM55_FAM56_FAM56_FAM56_FAM57_FAM57_FAM57_FAM58_FAM58_FAM58_FAM59_FAM59_FAM59_FAM60_FAM60_FAM60', 'NA12329_NA06984_NA06989_NA12344_NA12347_NA12348_NA12335_NA12340_NA12341_NA07029_NA06994_NA07000_NA10852_NA12045_NA12046_NA10857_NA12043_NA12044_NA10855_NA11831_NA11832_NA10856_NA11829_NA11830_NA12376_NA12546_NA12489_NA12386_NA12399_NA12400', 'NA06984_0_0_NA12347_0_0_NA12340_0_0_NA06994_0_0_NA12045_0_0_NA12043_0_0_NA11831_0_0_NA11829_0_0_NA12546_0_0_NA12399_0_0', 'NA06989_0_0_NA12348_0_0_NA12341_0_0_NA07000_0_0_NA12046_0_0_NA12044_0_0_NA11832_0_0_NA11830_0_0_NA12489_0_0_NA12400_0_0', '/data/smew1/rdavies/ukbb_gel_2023_01_26/CEU/CEU-chr20-final.b38.txt.gz', 'chr20:19000000-31000000', '/data/smew1/rdavies/hrc/hg38_liftover_clean_chr20.vcf.gz', '16')
    
}



vcffile <- args[1]
outstem <- args[2]
fid <- strsplit(args[3], "_")[[1]]
iid <- strsplit(args[4], "_")[[1]]
parent1 <- strsplit(args[5], "_")[[1]]
parent2 <- strsplit(args[6], "_")[[1]]
recomb_file <- args[7]
region <- args[8]
refpanel <- args[9]
nthreads <- args[10]

## argh - hack baby
if (length(grep("FamilyD", iid)) > 0) {
    iid <- paste0(iid[seq(1, length(iid), 2)], "_", iid[seq(2, length(iid), 2)])
    parent1 <- paste0(parent1[seq(1, length(iid), 2)], "_", parent1[seq(2, length(iid), 2)])
    parent2 <- paste0(parent2[seq(1, length(iid), 2)], "_", parent2[seq(2, length(iid), 2)])
}


library("data.table")
vcf <- fread(cmd = paste0("gunzip -c ", vcffile, "| grep -v '##'"), data.table = FALSE)
phased_vcf <- vcf
phased_vcf[, -c(1:9)] <- "./."
for(i in 10:ncol(vcf)) {
    x <- vcf[, i]
    g <- substr(x, 1, 3)
    vcf[, i] <- g
    w <- g == "0/0"
    phased_vcf[w, i] <- "0|0"
    w <- g == "1/1"
    phased_vcf[w, i] <- "1|1"
}


change <- rbind(
    c(0,2,1,0,1),
    c(0,1,1,0,1),
    c(1,0,1,1,0),
    c(1,1,1,NA,NA),
    c(1,2,1,0,1),
    c(2,0,1,1,0),
    c(2,1,1,1,0)
)



print(colnames(vcf))
a <- match(iid, colnames(vcf))
if (sum(is.na(a)) > 0) {
    print("cannot find something")
    print(a)
    print(iid)
    stop("problem!")
}

##
## now, want to fix all sites that are not triple het for each family
##
for(family in unique(fid)) {
    ## proband is the one with parents
    fidL <- fid[fid == family]
    iidL <- iid[fid == family]
    p1L <- parent1[fid == family]
    p2L <- parent2[fid == family]
    i_pros <- which(p1L != "0")
    for(i_pro in i_pros) {
        i_p1 <- setdiff(1:length(fidL), i_pro)[1]
        i_p2 <- setdiff(1:length(fidL), i_pro)[2]
        stopifnot(length(i_pro) == 1)
        ##
        ## manually phase hets here
        ##
        is_het <- vcf[, iidL[i_pro]] == "0/1"
        pro <- vcf[is_het, iidL[i_pro]]
        par1 <- vcf[is_het, iidL[i_p1]]
        par2 <- vcf[is_het, iidL[i_p2]]
        is_het <- which(vcf[, iidL[i_pro]] == "0/1")
        nHetSNPs <- length(is_het)
        ##    
        for(i in c(1,2,3,5,6,7)) {
            x <- change[i, ]
            pro_gen <- c("0/0", "0/1", "1/1")[x[3] + 1]
            p1_gen <- c("0/0", "0/1", "1/1")[x[1] + 1]
            p2_gen <- c("0/0", "0/1", "1/1")[x[2] + 1]      
            phase1 <- x[4]
            phase2 <- x[5]
            w <- (par1 == p1_gen) & (par2 == p2_gen) & (pro == pro_gen)
            print(paste0("There are ", sum(w), " type ", i, " SNPs"))
            v <- paste0(phase1, "|", phase2)
            phased_vcf[is_het[w], iidL[i_pro]] <- v
        }
        ##
        ##
    }
}


## write out to disk here
system(paste0("gunzip -c ", vcffile, " | head -n 10000 | grep '##'  > ", outstem, ".scaffold.vcf"))
fwrite(phased_vcf, file = paste0(outstem, ".scaffold.vcf"), row.names = FALSE, col.names = TRUE, sep = "\t", append = TRUE)
system(paste0("bgzip -f ", outstem, ".scaffold.vcf"))
system(paste0("bcftools index -f ", outstem, ".scaffold.vcf.gz"))

## make bcf
## system(paste0("bcftools view ", paste0(outstem, ".scaffold.vcf"), " -Ob -o ", paste0(outstem, ".scaffold.bcf")))
##system(paste0("bcftools view ", vcffile, " -Ob -o ", paste0(outstem, ".originalvcf.bcf")))
##system(paste0("bcftools index ", outstem, ".scaffold.bcf"))
##system(paste0("bcftools index ", outstem, ".originalvcf.bcf"))


## in shell
## . ~/.conda_activate
## conda create --name shapeit4
## conda install shapeit4=4.0 -c bioconda

## Holy fuck boost for shapeit4
## easier to conda it
## system(    ". ~/.conda_activate && conda activate shapeit4env && shapeit4 --help")
##     ". ~/.conda_activate && conda activate shapeit4_v4.0 && ",
input <- paste0(outstem, ".originalvcf.vcf")
cmd <- paste0(
    ". ~/.conda_activate && conda activate shapeit4 && ",
    "shapeit4 ",
    "--input ", vcffile, " ", 
    "--scaffold ", outstem, ".scaffold.vcf.gz", " ", 
    "--map ", recomb_file, " ", 
    "--region ", region, " ", 
    "--reference ", refpanel, " ",
    "--log ", outstem, ".log ",
    "--output ", outstem, ".phased.vcf.gz ",
    "--thread ", nthreads
)
print(paste0("running cmd:", cmd))
system(cmd)

## system(paste0("bcftools view ", outstem, ".phased.bcf -o ", outstem, ".phased.vcf.gz -O z"))


quit()

##    "--sequencing ",

## compare against original
after_phased_vcf <- fread(cmd = paste0("gunzip -c ", outstem, ".phased.vcf.gz| grep -v '##'"), data.table = FALSE)
## compare

table(after_phased_vcf[, 10], phased_vcf[, 10])

haps <- fread("/gpfs3/well/davies/users/dcc832/single_imp/2020_06_05/genotypes/phase.experimental.chr20.txt", data.table = FALSE)
## so pretty good, some disagreements, not a huge number, and they are always me being ref/ref
table(paste0(haps[, 5], haps[, 6]), after_phased_vcf[, 10])
2
