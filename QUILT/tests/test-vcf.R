library(QUILT)

vcffile <- "/Users/zilong/Project/zll/vcfpp/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
N <- 16325
M <- 2504 * 2

pbwt <- pbwt_build(vcffile = vcffile , samples = "-",  region = "22:19500000-20000000")

str(pbwt)
