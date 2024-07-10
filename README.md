QUILT2
=====
**__Current Version: 2.0.0__**
Release date: July 06, 2024

![Build Status](https://github.com/rwdavies/QUILT/workflows/CI/badge.svg)

We are excited about the upgraded QUILT, that is QUILT2. For details of past changes please see [CHANGELOG](CHANGELOG.md).

QUILT2 is an R and C++ program for fast genotype imputation from low-coverage sequence using a large reference panel, which is also accurate and versatile for various data, such as ***short reads, long reads, ancient DNA and cell-free DNA from NIPT***.

Please use this README for general information about QUILT2. Also, specific information about original QUILT1 and QUILT-HLA, please see the [QUILT1 README](README_QUILT1.md) and [QUILT-HLA README](README_QUILT-HLA.md) for specific details.

# Table of contents
1. [Introduction](#paragraph-introduction)
2. [Installation](#paragraph-installation)
3. [Quick start run](#paragraph-quickstartrun)
4. [Tutorials](#paragraph-tutorial)
5. [License](#paragraph-license)
6. [Citation](#paragraph-citation)
7. [Testing](#paragraph-testing)
8. [Bug reports](#paragraph-bugreports)


## Introduction <a name="paragraph-introduction"></a>

First of all, QUILT2 is a major upgrade of QUILT1, and is fast and memory-efficient for millions reference haplotypes. Methodologically, in addition to diploid method, QUILT2 can impute the mother and fetus simultaneously using cell-free DNA from NIPT. Computationally, QUILT2 introduces [mspbwt](https://github.com/rwdavies/mspbwt) and a two-stage imputation strategy for rare and common variants. Statistically, QUILT2 models directly on a per-read basis, and is base quality aware, meaning it can accurately impute from diverse inputs, including noisy long read sequencing (e.g. Oxford Nanopore Technologies), barcoded Illumina sequencing (e.g. Haplotagging) and ancient DNA. Accuracy using QUILT2 and lc-WGS meets or exceeds other methods for lc-WGS imputation, particularly for high diversity regions or genomes (e.g. MHC, or non-human species). Relative to DNA genotyping microarrays, QUILT2 offers improved accuracy at reduced cost, particularly for diverse populations, with the potential for accuracy to nearly double at rare SNPs. Further details and detailed evaluations are available in the [QUILT paper](README.md#paragraph-citation).

## Installation <a name="paragraph-installation"></a>

QUILT2 is available to download and install through this Github repository.
QUILT2 depends on [STITCH](https://github.com/rwdavies/STITCH) and [mspbwt](https://github.com/rwdavies/mspbwt).

```
git clone --recursive https://github.com/rwdavies/QUILT.git
cd QUILT
bash ./scripts/install-dependencies.sh ## skip this if STITCH and mspbwt already installed
Rscript ./scripts/build-and-install.R
```

## Quick start run <a name="paragraph-quickstartrun"></a>

A quick start to ensure QUILT is properly installed and working can be performed using the following.

Download example data package, containing 1000 Genomes haplotypes, and NA12878 bams

```
wget https://zenodo.org/records/12697284/files/QUILT2_example_2024.tar.xz ## or curl -O
tar --xz -xf QUILT2_example_2024.tar.xz
```

Perform imputation as below. Note that reference panel data can be processed separately to speed up repeated imputation of the same region in independent jobs, for a detailed tutorial on preparing reference, options, input/output files, see [Tutorials](#paragraph-tutorial)


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

## Tutorials <a name="paragraph-tutorial"></a>

- Please use the **[QUILT2 Tutorial](README_QUILT2.md)** for guide of running QUILT2 efficiently. 
- Please use the **[QUILT1 Tutorial](README_QUILT1.md)** for guide of running QUILT1 efficiently. 
- Please use the **[QUILT-HLA Tutorial](README_QUILT-HLA.md)** for guide of running QUILT-HLA efficiently. 

## License <a name="paragraph-license"></a>

QUILT and the code in this repo is available under a GPL3 license. For more information please see the [LICENSE](LICENSE).

## Citation <a name="paragraph-citation"></a>

Davies, R. W., Kucka M., Su D., Shi S., Flanagan M., Cunniff C. M., Chan Y. F. , Myers S. Rapid genotype imputation from sequence with reference panels. In press, Nature Genetics

## Testing <a name="paragraph-testing"></a>

Tests in QUILT are split into unit or acceptance run using ```./scripts/test-unit.sh``` and ```./scripts/test-acceptance.sh```. To run all tests use ```./scripts/all-tests.sh```, which also builds and installs a release version of QUILT. To make compilation go faster do something like ```export MAKE="make -j 8"```.

## Bug reports <a name="paragraph-bugreports"></a>

The best way to get help is to submit a bug report on GitHub in the [Issues](https://github.com/rwdavies/STITCH/issues) section. Please also use the Issues section if you have a more general question, such issues will be left open for others to see. Similarly, please check the issues before posting to see if your issue has already been addressed

For more detailed questions or other concerns please contact Robert Davies robertwilliamdavies@gmail.com
 
