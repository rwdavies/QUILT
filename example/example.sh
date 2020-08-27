#!/usr/bin/env bash

set -e

## script_dir=`dirname "$0"`
## cd "${script_dir}"/../


##
## example of running QUILT
## using either approx 5000 or 50000 haplotypes
## step 0 - download some example data
## step 1 - prepare a reference package using a haplotype reference file
## step 2 - impute samples using the prepared reference package



##
## step 0 - download some example data
##
os=`uname -a | cut -f 1 -d " "`
if [ "${os}" == "Darwin" ]
then
   download_command="curl -O"
elif [ "${os}" == "Linux" ]
then
    download_command="wget"
else
    echo "cannot figure out operating system for selecting download command, please run file interactively"
    exit 1
fi
## download package - this should download the files from this website including for example "bamlist.1.0.txt"
if [ ! -f QUILT_example_2020_08_25.tgz ]
then
    ${download_command} http://www.stats.ox.ac.uk/~rdavies/QUILT_example_2020_08_25.tgz ## either wget or curl -O or otherwise as appropriate
fi
tar -xzvf QUILT_example_2020_08_25.tgz





##
## step 1 - prepare a reference package using a haplotype reference file
##
rm -r -f quilt_output
./QUILT_prepare_reference.R \
    --outputdir=quilt_output \
    --chr=chr20 \
    --nGen=100 \
    --reference_haplotype_file=package_2020_08_25/ref.chr20.2000001.4000000.hap.clean.example.1000Gonly.gz \
    --reference_legend_file=package_2020_08_25/ref.chr20.2000001.4000000.legend.clean.example.gz \
    --genetic_map_file=package_2020_08_25/CEU-chr20-final.b38.txt.gz \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000
## NOTES
## 1) succesful completion of this run results in a file at quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData, which is used by QUILT itself during the actual imputation in the next step
## 2) there are two haplotype reference file options included in this package
##   ref.chr20.2000001.4000000.hap.clean.example.1000Gonly.gz and
##   ref.chr20.2000001.4000000.hap.clean.example.gz
##   the former is just 1000 Genomes samples
##   the later is the full size of the HRC but all non-1000G sample genoytpes have been set to 0 and sample names redatcted
##   the site list for both is the HRC site list which has been released publicly
##   as such the larger haplotype reference panel just demonstrates performance run times of QUILT rather than generating more accurate imputation
##   the choice of file can be modified using the reference_haplotype_file option above



## 
## step 2 - impute samples using the prepared reference package
## 
./QUILT.R \
    --outputdir=quilt_output \
    --chr=chr20 \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000 \
    --bamlist=package_2020_08_25/bamlist.1.0.txt \
    --posfile=package_2020_08_25/pos.chr20.2000001.4000000.txt \
    --phasefile=package_2020_08_25/phase.chr20.2000001.4000000.txt \
    --bqFilter=10 \
    --nCores=1
## NOTES
## 1) this imputes NA12878 on this region, using a haplotagged Illumina example and ONT
##   succesful completion of this run results in a VCF at quilt_output/quilt.chr20.2000001.4000000.vcf.gz
## 2) per-sample per-SNP output includes hard-called (integer) genotype, genotype posteriors, diploid dosage, and haploid dosages (i.e. haplotypes, can be turned into standard 0 or 1 haplotpes using rounding)
## 3) you can try bamlist.0.25.txt to try 0.25X bams
## 4) try ./QUILT.R --help for more options
## 5) for normal operation when you do not have high quality phased truth genotypes, you can omit posfile and phasefile, and sites will be imputed 

