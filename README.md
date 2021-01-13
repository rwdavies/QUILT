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
4. [Help, options and parameters](#paragraph-helpoptionsparams)
5. [Important parameters that influence run time, memory usage and accuracy](#paragraph-paramsimportant)
6. [Examples](#paragraph-examples)
7. [License](#paragraph-license)
8. [Citation](#paragraph-citation)
9. [Testing](#paragraph-testing)
10. [Bug reports](#paragraph-bugreports)



## Introduction <a name="paragraph-introduction"></a>
Some introduction text

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

Download data package
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
```
Succesful completion of this run results in a VCF at `quilt_output/quilt.chr20.2000001.4000000.vcf.gz`. For more details on this and other examples, see [Examples](#paragraph-examples)

## Help, options and parameters <a name="paragraph-helpoptionsparams"></a>

For a full list of options, query `?QUILT::QUILT`, or alternatively, type 
```
./QUILT.R --help
```



## Important parameters that influence run time, memory usage and accuracy <a name="paragraph-paramsimportant"></a>

Text here


## Examples <a name="paragraph-examples"></a>

In this directory you will find [example/examples.sh](example/examples.sh), a script with examples, including the above. They can either be run in their entirety using `./example/examples.sh`, or can be run interactively line by line.

## License <a name="paragraph-license"></a>

Text here

## Citation <a name="paragraph-citation"></a>

Text here

## Testing <a name="paragraph-testing"></a>

Text here

## Bug reports <a name="paragraph-bugreports"></a>

Text here









