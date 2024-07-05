QUILT
=====
**__Current Version: 1.0.5__**
Release date: Sept 11, 2023

![Build Status](https://github.com/rwdavies/QUILT/workflows/CI/badge.svg)

Changes in latest version

1. Be able to work with cram files

For details of past changes please see [CHANGELOG](CHANGELOG.md).

QUILT is an R and C++ program for rapid genotype imputation from low-coverage sequence using a large reference panel.

QUILT-HLA is an R and C++ program for rapid HLA imputation from low-coverage sequence using a labelled reference panel.

Please use this README for general information about QUILT and QUILT-HLA, and specific information about QUILT. Please see the [QUILT-HLA README](README_QUILT-HLA.md) for specific details about QUILT-HLA.

# Table of contents
1. [Introduction](#paragraph-introduction)
2. [Installation](#paragraph-installation)
    1. [github](#paragraph-installation-github)
    2. [conda](#paragraph-installation-conda)
3. [Quick start run](#paragraph-quickstartrun)
4. [Input and output formats](#paragraph-io)
    1. [Input](#paragraph-io-input)
    2. [Output](#paragraph-io-output)
5. [Help, options and parameters](#paragraph-helpoptionsparams)
6. [Separating reference panel processing from imputation](#paragraph-separate)
7. [Important parameters that influence run time and accuracy](#paragraph-paramsimportant)
8. [Examples](#paragraph-examples)
9. [License](#paragraph-license)
10. [Citation](#paragraph-citation)
11. [Testing](#paragraph-testing)
12. [Bug reports](#paragraph-bugreports)


## Introduction <a name="paragraph-introduction"></a>

QUILT is a program for rapid diploid genotype imputation from low-coverage sequence using a large reference panel. Statistically, the QUILT model works on a per-read basis, and is base quality aware, meaning it can accurately impute from diverse inputs, including noisy long read sequencing (e.g. Oxford Nanopore Technologies), and barcoded Illumina sequencing (e.g. Haplotagging). Accuracy using QUILT and lc-WGS meets or exceeds other methods for lc-WGS imputation, particularly for high diversity regions or genomes (e.g. MHC, or non-human species). Relative to DNA genotyping microarrays, QUILT offers improved accuracy at reduced cost, particularly for diverse populations, with the potential for accuracy to nearly double at rare SNPs (e.g. 2.0X lc-WGS vs microarrays for SNPs at 0.1% frequency). Further details and detailed evaluations are available in the [QUILT paper](README.md#paragraph-citation).

## Installation <a name="paragraph-installation"></a>

QUILT is available to download either through this github repository, or through conda.

### github <a name="paragraph-installation-github"></a>

First, install STITCH, installed in a similar way to QUILT, as specified on the STITCH website [here](https://github.com/rwdavies/STITCH). Next, install QUILT, as follows

```
git clone --recursive https://github.com/rwdavies/QUILT.git
cd QUILT
./scripts/install-dependencies.sh
cd releases
wget https://github.com/rwdavies/quilt/releases/download/1.0.5/QUILT_1.0.5.tar.gz ## or curl -O
R CMD INSTALL QUILT_1.0.5.tar.gz
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



## Quick start run <a name="paragraph-quickstartrun"></a>

A quick start to ensure QUILT is properly installed and working can be performed using the following

Download example data package, containing 1000 Genomes haplotypes, and NA12878 bams
```
wget http://www.stats.ox.ac.uk/~rdavies/QUILT_example_2021_01_15A.tgz ## or curl -O
tar -xzvf QUILT_example_2021_01_15A.tgz
```

Perform imputation. Note that reference panel data can be processed separately to speed up repeated imputation of the same region in independent jobs, see [Separating reference panel processing from imputation](#paragraph-separate).
```
rm -r -f quilt_output
./QUILT.R \
--outputdir=quilt_output \
--chr=chr20 \
--regionStart=2000001 \
--regionEnd=2100000 \
--buffer=10000 \
--bamlist=package_2021_01_15A/bamlist.1.0.txt \
--posfile=package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.posfile.txt \
--phasefile=package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.phasefile.txt \
--reference_haplotype_file=package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.noNA12878.hap.gz \
--reference_legend_file=package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.noNA12878.legend.gz \
--genetic_map_file=package_2021_01_15A/CEU-chr20-final.b38.txt.gz \
--nGen=100 \
--save_prepared_reference=TRUE
```
Succesful completion of this run results in a VCF at `quilt_output/quilt.chr20.2000001.2100000.vcf.gz`. For a slightly longer version of this example, see [Examples](#paragraph-examples)



## Input and output formats <a name="paragraph-io"></a>

### Input <a name="paragraph-io-input"></a>

For all of these, it can be useful to take a look at the example files provided as part of the quick start example above.

- Reference panel. IMPUTE format hap and legend format files with reference haplotypes. These can be made from haplotype VCFs using `bcftools convert --haplegendsample`. Alternatively, they can be made manually. The haplotype file is a gzipped file with no header and no rownames, with one row per SNP, with one column per reference haplotype, space separated, and values of 0 (ref) and 1 (alt). The legend file is a gzipped file with no rownames, a header file including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele. An optional sample file and file with samples to exclude can be useful for changing who is used in the reference panel.
- Genetic map. File with genetic map information, with 3 white-space delimited columns giving position (1-based), genetic rate map in cM/Mbp, and genetic map in cM
- Bams. Given as a bamlist (i.e. a file with one row per sample, the path to the bam)
- (Optional) Truth data. phasefile and posfile. Useful for understanding performance. Phasefile has a header row with a name for each sample, matching what is found in the bam file. File is tab separated, one subject per column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header. For posfile, this is a file with positions of where to impute, lining up one-to-one with the SNPs of phasefile. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1\<tab\>1000\<tab\>A\<tab\>G\<tab\>

### Output <a name="paragraph-io-output"></a>

- VCF with both SNP annotation information (see below) and per-sample genotype information. Per-sample genotype information includes the following entries

- GT **Phased genotypes** Phased genotype, where each allele is the rounded per-haplotype posterior probability (HD below)
- GP **Genotype posteriors** Posterior probabilities of the three genotypes given the data
- DS **Diploid dosage** Posterior expectation of the diploid genotype i.e. the expected number of copies of the alternate allele
- HD **Haploid dosages** Per-haplotype posterior probability of an alternate allele

Note that in QUILT, genotype posteriors (GP) and dosages (DS) are taken from the main Gibbs sampling, while the phasing results (GT and HD) are taken from an additional special phasing Gibbs sample. As such, phasing results (GT and HD) might not be consistent with genotype information (GP and DS). If consistency is necessary, note that you can create a consistent GP and DS from HD.

Per-SNP annotation is available as follows
```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">,
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Posterior genotype probability of 0/0, 0/1, and 1/1">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Diploid dosage">
##FORMAT=<ID=HD,Number=2,Type=Float,Description="Haploid dosages">
```

SNP annotation information
```
##INFO=<ID=EAF,Number=.,Type=Float,Description="Estimated allele frequency">
##INFO=<ID=HWE,Number=.,Type=Float,Description="Hardy-Weinberg p-value">
##INFO=<ID=ERC,Number=.,Type=Float,Description="Estimated number of copies of the reference allele from the pileup">
##INFO=<ID=EAC,Number=.,Type=Float,Description="Estimated number of copies of the alternate allele from the pileup">
##INFO=<ID=PAF,Number=.,Type=Float,Description="Estimated allele frequency using the pileup of reference and alternate alleles">
```


## Separating reference panel processing from imputation <a name="paragraph-separate"></a>

For large reference panels, and for many jobs involving imputing few samples, it can be computationally efficient to pre-process the reference panel and save the output, and use this output for multiple independent runs. Here is an example for how we would do this, for the case of the quick start example. Note that any parameters available jointly in `QUILT` and `QUILT_prepare_reference` that inform how the reference panel is processed must be set in `QUILT_prepare_reference` (for example, `maxRate` bounds the recombination rate, and must be set when running`QUILT_prepare_reference` as the recombination rate is processed in this step).


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
Note that when running multiple versions of QUILT against the same reference data, it is useful to set `output_filename` to change the default filename for each job, and to keep the temporary directories used independent (which is the behaviour for default `tempdir`).

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



## Examples <a name="paragraph-examples"></a>

- **[example/QUILT_usage.Md](example/QUILT_usage.Md)** A larger version of the quick start example. This can be run using either  `./example/run_example.sh example/QUILT_usage.Md`, or can be run interactively line by line from the [example/QUILT_usage.Md](example/QUILT_usage.Md) file.
- **[example/ligation.Md](example/ligation.Md)** An example of how to run QUILT in chunks and ligate results together to get correctly oriented phased results across VCFs.

## License <a name="paragraph-license"></a>

QUILT and the code in this repo is available under a GPL3 license. For more information please see the [LICENSE](LICENSE).

## Citation <a name="paragraph-citation"></a>

Davies, R. W., Kucka M., Su D., Shi S., Flanagan M., Cunniff C. M., Chan Y. F. , Myers S. Rapid genotype imputation from sequence with reference panels. In press, Nature Genetics

## Testing <a name="paragraph-testing"></a>

Tests in QUILT are split into unit or acceptance run using ```./scripts/test-unit.sh``` and ```./scripts/test-acceptance.sh```. To run all tests use ```./scripts/all-tests.sh```, which also builds and installs a release version of QUILT. To make compilation go faster do something like ```export MAKE="make -j 8"```.

## Bug reports <a name="paragraph-bugreports"></a>

The best way to get help is to submit a bug report on GitHub in the [Issues](https://github.com/rwdavies/STITCH/issues) section. Please also use the Issues section if you have a more general question, such issues will be left open for others to see. Similarly, please check the issues before posting to see if your issue has already been addressed

For more detailed questions or other concerns please contact Robert Davies robertwilliamdavies@gmail.com
 
