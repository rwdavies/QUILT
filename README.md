QUILT
=====
**__Current Version: pre-release__**
Release date: pre-release

For details of past changes please see [CHANGELOG](CHANGELOG.md).

QUILT is an R and C++ program for rapid genotype imputation from low-coverage sequence using a large reference panel.

# Table of contents
[1. Introduction](#paragraph-introduction)
[2. Installation](#paragraph-installation)
    [1. github](#paragraph-installation-github)
    [2. conda](#paragraph-installation-conda)
[3. Quick start](#paragraph-quickstart)

## Introduction <a name="paragraph-introduction"></a>
Some introduction text



## installation <a name="paragraph-installation"></a>

QUILT is available to download either through this github repository, or through conda.

### github <a name="paragraph-installation-github"></a>

First, install STITCH, installed in a similar way to QUILT, as specified on the STITCH website [here](https://github.com/rwdavies/STITCH). Next, install QUILT, as follows

```
git clone --recursive https://github.com/rwdavies/QUILT.git
cd QUILT
./scripts/install-dependencies.sh
cd releases
wget https://github.com/rwdavies/quilt/releases/download/1.0.0/QUILT_1.0.0.tar.gz ## or curl -O
R CMD INSTALL QUILT_1.0.0.tar.gz
```

### conda <a name="paragraph-installation-conda"></a>
QUILT (as r-quilt) can be installed using [conda](https://conda.io/miniconda.html). Full tutorials can be found elsewhere, but briefly, something like this should work
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install r-quilt -c defaults -c bioconda -c conda-forge
source activate
R -e 'library("QUILT")'
```
Note that currently the command like `QUILT.R` is not included with the bioconda installation, so from the command line, you can either run something like `R -e 'library("QUILT"); QUILT(chr="chr19", etc)'`, or clone the repo to get `QUILT.R`. 



## Quick start <a name="paragraph-quickstart"></a>

A quick start to ensure QUILT is properly installed and working can be performed using the following

Download package
```
wget http://www.stats.ox.ac.uk/~rdavies/QUILT_example_2020_08_25.tgz ## or curl -O
tar -xzvf QUILT_example_2020_08_25.tgz
```

First, to re-format the reference panel
```
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
```

Second, to perform imputation
```
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
```







## QUILT HLA

### Preparing files

#### Preparing IPD-IGMT files

Download IPD_IGMT files
```
## To download IPD-IGMT version 3.39, for example
wget https://github.com/ANHIG/IMGTHLA/blob/032815608e6312b595b4aaf9904d5b4c189dd6dc/Alignments_Rel_3390.zip?raw=true
mv Alignments_Rel_3390.zip?raw=true Alignments_Rel_3390.zip
```
Prepare supplementary information file (or use provided one, if using above release, and GRCh38)

Prepare mapping related files
```
./QUILT_HLA_prepare_reference.R \
--outputdir=/well/davies/users/dcc832/single_imp/HLA_TEST_2021_01_06/ \
--ipd_igmt_alignments_zip_file=Alignments_Rel_3390.zip \
--quilt_hla_supplementary_info_file=quilt_hla_supplementary_info.txt
```

#### Preparing haplotype files
```
HLA_GENE="A"
regionStart=29942554
regionEnd=29945741
    
./QUILT_prepare_reference.R \
--outputdir=/well/davies/users/dcc832/single_imp/HLA_TEST_2021_01_06/ \
--nGen=100 \
--reference_haplotype_file=/well/davies/users/dcc832/single_imp/2020_06_25/ref_panels/hrc.chr6.hap.clean.gz \
--reference_legend_file=/well/davies/users/dcc832/single_imp/2020_06_25/ref_panels/hrc.chr6.legend.clean.gz \
--reference_sample_file=/well/davies/users/dcc832/single_imp/2020_06_25/ref_panels/hrc.chr6.samples.reheadered2 \
--chr=chr6 \
--regionStart=${regionStart} \
--regionEnd=${regionEnd} \
--buffer=500000 \
--genetic_map_file=/well/davies/shared/recomb/CEU/CEU-chr6-final.b38.txt.gz \
--reference_exclude_samplelist_file=/well/davies/shared/1000G/robbie_files/hlauntyped${HLA_GENE}.excludefivepop.txt \
--output_file=quilt.hrc.chr6.hla.${HLA_GENE}.haplotypes.RData \
--region_exclude_file=/well/davies/shared/1000G/robbie_files/hlagenes.txt \
--minRate=0.01
```
