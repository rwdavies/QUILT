QUILT2: genotype calling from low coverage reads
=====
<!-- badges: start -->
![Build Status](https://github.com/rwdavies/QUILT/workflows/CI/badge.svg)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/r-quilt/README.html)
![version](https://anaconda.org/bioconda/r-quilt/badges/version.svg)
<!-- badges: end -->

**__Current Version: 2.0.0__**. Release date: July 06, 2024

We are excited to release QUILT2, a major upgrade over QUILT (QUILT1). For details of past changes to QUILT2 and QUILT1, please see the [CHANGELOG](CHANGELOG.md).

QUILT2 is an R and C++ program for fast genotype imputation from low-coverage sequence using a large reference panel. QUILT2 is accurate and versatile, able to handle imputation from ***short read, long read, ancient DNA and cell-free DNA from NIPT***.

Please use this README for general information about QUILT2. For more information about QUILT1 and QUILT-HLA, please see the [QUILT1 README](README_QUILT1.md) and the [QUILT-HLA README](README_QUILT-HLA.md).

# Table of contents
1. [Introduction](#paragraph-introduction)
2. [Installation](#paragraph-installation)
3. [Quick start run](#paragraph-quickstartrun)
4. [READMEs and tutorials](#paragraph-readme-tutorial)
5. [NIPT](#paragraph-nipt)
6. [License](#paragraph-license)
7. [Citation](#paragraph-citation)
8. [Testing](#paragraph-testing)
9. [Bug reports](#paragraph-bugreports)
10. [Contacts](#paragraph-contacts)


## Introduction <a name="paragraph-introduction"></a>

QUILT2 is a fast and memory-efficient method for imputation from low coverage sequence. Statistically, QUILT2 operates on a per-read basis, and is base quality aware, meaning it can accurately impute from diverse inputs, including short read (e.g. Illumina), long read sequencing (that might be noisy) (e.g. Oxford Nanopore Technologies), barcoded Illumina sequencing (e.g. Haplotagging) and ancient DNA. In addition, QUILT2 can impute both the mother and fetal genome using cfDNA NIPT data. Methodologically, QUILT2 introduces [mspbwt](https://github.com/rwdavies/mspbwt), and a two-stage imputation strategy for rare and common variants, to facilitate analysis using haplotype reference panels derived from hundreds of thousands or millions of whole genome sequenced haplotypes. Accuracy using QUILT2 and lc-WGS meets or exceeds other methods for imputation, particularly for high diversity regions or genomes (e.g. MHC, or non-human species). Relative to DNA genotyping microarrays, QUILT2 offers improved accuracy at reduced cost, particularly for diverse populations, with the potential to nearly double accurate at rare SNPs. Links to published references with more details and detailed evaluations are available in the [Citation](#paragraph-citation).

## Installation <a name="paragraph-installation"></a>

QUILT2 is available on bioconda which can be installed by 

```
conda create -c conda-forge -c defaults -c bioconda -n quilt2 'r-quilt>=2.0.0'
```

Also, QUILT2 is available to download and install through this GitHub repository, which depends on [STITCH>=1.7.0](https://github.com/rwdavies/STITCH) and [mspbwt>=0.1.0](https://github.com/rwdavies/mspbwt).

```
git clone --recursive https://github.com/rwdavies/QUILT.git
cd QUILT
bash ./scripts/install-dependencies.sh ## skip this if STITCH>=1.7.0 and mspbwt>=0.1.0 installed
wget https://github.com/rwdavies/QUILT/releases/download/2.0.0/QUILT_2.0.0.tar.gz ## or curl -OL
R CMD INSTALL QUILT_2.0.0.tar.gz
```

## Quick start run <a name="paragraph-quickstartrun"></a>

A quick start to ensure QUILT2 is properly installed and working can be performed using the following. **Note QUILT2 only takes VCF file as input for reference panel.** 

Download example data package, containing 1000 Genomes haplotypes, and NA12878 bams

```
wget https://zenodo.org/records/12786681/files/QUILT2_example_2024.tar.xz  ## or curl -OL
tar --xz -xf QUILT2_example_2024.tar.xz
```

Perform imputation as below. Note that reference panel data can be processed separately to speed up repeated imputation of the same region in independent jobs. For a detailed tutorial on preparing reference, options, input/output files, see [READMEs and tutorials](#paragraph-readme-tutorial)


```
outdir=quilt2_output && rm -rf $outdir
./QUILT2.R \
--outputdir=$outdir \
--chr=chr20 \
--regionStart=2000001 \
--regionEnd=4000000 \
--buffer=500000 \
--nGen=100 \
--bamlist=package/bamlist.1.0.txt \
--genetic_map_file=package/CEU-chr20-final.b38.txt.gz \
--reference_vcf_file=package/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.noNA12878.vcf.gz \
--save_prepared_reference=TRUE
```

## READMEs and tutorials <a name="paragraph-readme-tutorial"></a>

- **[QUILT2 README](README.md)**. Please see the **[QUILT2 Tutorial](README_QUILT2.org)** for a guide and example code for running QUILT2 efficiently. 
- **[QUILT1 README](README_QUILT1.md)**. For tutorials, please see **[Examples](README_QUILT1.md#paragraph-examples)** section of the QUILT1 README.
- **[QUILT-HLA README](README_QUILT-HLA.md)**. For tutorials, please see **[Example making and running a reference panel](QUILT/example/QUILT_hla_reference_panel_construction.Md)**.

Note that in the future, the QUILT1 and QUILT-HLA READMEs will not be further updated. All of `QUILT2.R`, `QUILT.R` and `QUILT_HLA.R` can be accessed through the main QUILT github page, which serves all three. Note that `QUILT2.R` and `QUILT.R` differ only in that in `QUILT2.R`, msPBWT is turned on by default `use_mspbwt=TRUE`, as is the two stage strategy of first using common variants, and then using all variants `impute_rare_common=TRUE`.


## NIPT <a name="paragraph-nipt"></a>

QUILT2 can impute using cfDNA NIPT. To impute using cfDNA NIPT, you need to set `method=nipt` and specify a path to `fflist`, which is a textfile, with one entry per line, with a one to one correspondence to `bamlist`, with the fetal fraction of that sample  (i.e. the first row of fflist is the fetal fraction for the first entry of bamlist, the second row is for the second entry of bamlist, etc.). For more details, including output details, please see [here](README_QUILT2.org#perform-nipt-imputation). 

## License <a name="paragraph-license"></a>

QUILT2, QUILT and QUILT_HLA, and the code in this repo is available under a GPL3 license. For more information please see the [LICENSE](LICENSE).

## Citation <a name="paragraph-citation"></a>

- For QUILT2, please cite Li Z., Albrechtsen A., Davies R. W. Rapid and accurate genotype imputation from low coverage short read, long read, and cell free DNA sequence. bioRxiv. [https://doi.org/10.1101/2024.07.18.604149](https://doi.org/10.1101/2024.07.18.604149)
- For QUILT1, please cite: Davies, R. W., Kucka M., Su D., Shi S., Flanagan M., Cunniff C. M., Chan Y. F. , Myers S. Rapid genotype imputation from sequence with reference panels. *Nat. Genet.* 53, 1104–1111 (2021)

## Testing <a name="paragraph-testing"></a>

Tests in QUILT2 are split into unit or acceptance run using ```./scripts/test-unit.sh``` and ```./scripts/test-acceptance.sh```. To run all tests use ```./scripts/all-tests.sh```, which also builds and installs a release version of QUILT2. To make compilation go faster do something like ```export MAKE="make -j 8"```.

## Bug reports <a name="paragraph-bugreports"></a>

The best way to get help is to submit a bug report on GitHub in the [Issues](https://github.com/rwdavies/QUILT/issues) section. Please also use the Issues section if you have a more general question, such issues will be left open for others to see. Similarly, please check the issues before posting to see if your issue has already been addressed

## Contacts <a name="paragraph-contacts"></a>

For more detailed questions or other concerns please contact Robert Davies robertwilliamdavies@gmail.com and Zilong Li zilong.dk@gmail.com
 
