QUILT
=====
**__Current Version: pre-release__**
Release date: pre-release

QUILT is an R and C++ program for rapid genotype imputation from low-coverage sequence using a large reference panel.

QUILT-HLA is an R and C++ program for rapid HLA imputation from low-coverage sequence using a labelled reference panel.

Please use this README for general information about QUILT and QUILT-HLA, and specific information about QUILT. Please see the [QUILT-HLA README](README_QUILT-HLA.md) for specific details about QUILT-HLA.

For details of past changes please see [CHANGELOG](CHANGELOG.md).

# Table of contents
1. [Introduction](#paragraph-introduction)
2. [Installation](#paragraph-installation)
    1. [github](#paragraph-installation-github)
    2. [conda](#paragraph-installation-conda)
3. [Quick start](#paragraph-quickstart)
4. [Input and output formats](#paragraph-io)
    1. [Input](#paragraph-io-input)
    2. [Output](#paragraph-io-output)
5. [Help, options and parameters](#paragraph-helpoptionsparams)
6. [Important parameters that influence run time and accuracy](#paragraph-paramsimportant)
7. [Examples and plots](#paragraph-examples)
8. [License](#paragraph-license)
9. [Citation](#paragraph-citation)
10. [Testing](#paragraph-testing)
11. [Bug reports](#paragraph-bugreports)


## Introduction <a name="paragraph-introduction"></a>

Forthcoming. 

## Installation <a name="paragraph-installation"></a>

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

Download data package, containing 1000 Genomes haplotypes, and NA12878 bams
```
wget http://www.stats.ox.ac.uk/~rdavies/QUILT_example_2021_01_15A.tgz ## or curl -O
tar -xzvf QUILT_example_2021_01_15A.tgz
```

First, to re-format the reference panel
```
rm -r -f quilt_output
./QUILT_prepare_reference.R \
--outputdir=quilt_output \
--chr=chr20 \
--nGen=100 \
--reference_haplotype_file=package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.noNA12878.hap.gz \
--reference_legend_file=package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.noNA12878.legend.gz \
--genetic_map_file=package_2021_01_15A/CEU-chr20-final.b38.txt.gz \
--regionStart=2000001 \
--regionEnd=2100000 \
--buffer=10000
```

Second, to perform imputation
```
./QUILT.R \
--outputdir=quilt_output \
--chr=chr20 \
--regionStart=2000001 \
--regionEnd=2100000 \
--buffer=10000 \
--bamlist=package_2021_01_15A/bamlist.1.0.txt \
--posfile=package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.posfile.txt \
--phasefile=package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.phasefile.txt
```
Succesful completion of this run results in a VCF at `quilt_output/quilt.chr20.2000001.2100000.vcf.gz`. For a slightly longer version of this example, see [Examples](#paragraph-examples)



## Input and output formats <a name="paragraph-io"></a>

### Input <a name="paragraph-io-input"></a>

For all of these, it can be useful to take a look at the example files provided as part of the quick start example above.

- Reference panel. IMPUTE format hap and legend format files with reference haplotypes. These can be made from haplotype VCFs using `bcftools convert --haplegendsample`. Alternatively, they can be made manually. The haplotype file is a gzipped file with no header and no rownames, with one row per SNP, with one column per reference haplotype, space separated, and values of 0 (ref) and 1 (alt). The legend file is a gzipped file with no rownames, a header file including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele. An optional sample file and file with samples to exclude can be useful for changing who is used in the reference panel.
- Genetic map. File with genetic map information, with 3 white-space delimited columns giving position (1-based), genetic rate map in cM/Mbp, and genetic map in cM
- Bams. Given as a bamlist (i.e. a file with one row per sample, the path to the bam)
- (Optional) Truth data. phasefile and posfile. Useful for understanding performance. Phasefile has a header row with a name for each sample, matching what is found in the bam file. File is tab separated, one subject per column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header. For posfile, this is a file with positions of where to impute, lining up one-to-one with the SNPs of phasefile. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab>

### Output <a name="paragraph-io-output"></a>

- VCF. Details are given in the VCF header, which is copied here
```
##INFO=<ID=EAF,Number=.,Type=Float,Description="Estimated allele frequency">
##INFO=<ID=HWE,Number=.,Type=Float,Description="Hardy-Weinberg p-value">
##INFO=<ID=ERC,Number=.,Type=Float,Description="Estimated number of copies of the reference allele from the pileup">
##INFO=<ID=EAC,Number=.,Type=Float,Description="Estimated number of copies of the alternate allele from the pileup">
##INFO=<ID=PAF,Number=.,Type=Float,Description="Estimated allele frequency using the pileup of reference and alternate alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Most likely genotype, given posterior probability of at least 0.90">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Posterior genotype probability of 0/0, 0/1, and 1/1">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Diploid dosage">
##FORMAT=<ID=HD,Number=2,Type=Float,Description="Haploid dosages">
```



## Help, options and parameters <a name="paragraph-helpoptionsparams"></a>

For a full list of options, query `?QUILT::QUILT`, or alternatively, type 
```
./QUILT.R --help
```



## Important parameters that influence run time and accuracy <a name="paragraph-paramsimportant"></a>

These parameters are most likely to influence run time and accuracy

- nGibbsSamples. Number of Gibbs samples performed. Run time is approximately linear in this, and increasing it will increase accuracy, but with diminishing returns. Generally less important for higher coverage.
- n_seek_its. Number of iterations between Gibbs sampling on small reference panel, and imputing using the full reference panel to find good matches for the small reference panel. Run time is approximately linear in this. Setting it higher will increase accuracy but with diminishing returns. 
- Ksubset. Size of the small reference panel. Higher values might increase accuracy but will increase run time.



## Examples and plots <a name="paragraph-examples"></a>

In this directory you will find [example/README.Md](example/), with for now just one example, a larger version of the quick start example. This can be run using either  `./example/examples.sh`, or can be run interactively line by line from the [example/README.Md](example/) file. It also contains an explanation of what the plots are, and why they might be useful.

## License <a name="paragraph-license"></a>

QUILT and the code in this repo is available under a GPL3 license. For more information please see the [LICENSE](LICENSE).

## Citation <a name="paragraph-citation"></a>

Forthcoming

## Testing <a name="paragraph-testing"></a>

Tests in QUILT are split into unit or acceptance run using ```./scripts/test-unit.sh``` and ```./scripts/test-acceptance.sh```. To run all tests use ```./scripts/all-tests.sh```, which also builds and installs a release version of QUILT. To make compilation go faster do something like ```export MAKE="make -j 8"```.

## Bug reports <a name="paragraph-bugreports"></a>

The best way to get help is to submit a bug report on GitHub in the [Issues](https://github.com/rwdavies/STITCH/issues) section. Please also use the Issues section if you have a more general question, such issues will be left open for others to see. Similarly, please check the issues before posting to see if your issue has already been addressed

For more detailed questions or other concerns please contact Robert Davies robertwilliamdavies@gmail.com
