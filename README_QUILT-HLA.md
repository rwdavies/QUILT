QUILT-HLA
=========

For general details including installation, version, and changelog, see the main QUILT [README](https://github.com/rwdavies/QUILT).

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
6. [Preparing a reference package](#paragraph-preparing)
    1. [Preparing ancillary files](#paragraph-ancillary-files)
    2. [Preparing IPD-IGMT files](#paragraph-preparing-ipdigmt)
    3. [Preparing haplotype files](#paragraph-preparing-haplotypes)


## Introduction <a name="paragraph-introduction"></a>

Forthcoming. 

## Installation <a name="paragraph-installation"></a>

QUILT-HLA is installed through the installation of QUILT. No additional installation is necessary.

## Quick start run <a name="paragraph-quickstart"></a>

A quick start to ensure QUILT-HLA is properly installed and working can be performed using the following

First, download some prepared reference panel data. This reference panel package was prepared using [example/QUILT_hla_reference_panel_construction.Md](example/QUILT_hla_reference_panel_construction.Md), and uses data from IPD-IGMT version 3.39, 1000 Genomes Project haplotypes (20201028), and 1000 Genomes Project HLA types (20181129).

```
wget something
# or curl -O on mac

put into here quilt_hla_reference_data
ACTUALLY
name package, once ready
```

Next, download data package, containing NA12878 bam
```
Incorporate this into package
echo -e ${inputs_dir}"NA12878.mhc.2.0.bam" > bamlist.txt
echo -e ${inputs_dir}"NA18566.mhc.2.0.bam" >> bamlist.txt

wget something
# or curl -O on mac
```

Now, HLA imputation for a particular region (here A) can be done as follows
```

HLA_GENE="A"
./QUILT_HLA.R \
--outputdir=quilt_output \
--bamlist=bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=quilt_hla_reference_data \
--quilt_hla_haplotype_panelfile=quilt_hla_reference_data/quilt.hrc.hla.${HLA_GENE}.haplotypes.RData \
--hla_gene_region_file=hla_ancillary_files/hlagenes.txt \
--dict_file=hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict
```

## Input and output formats <a name="paragraph-io"></a>

### Input <a name="paragraph-io-input"></a>

For all of these, it can be useful to take a look at the example files provided as part of the quick start example above.

- Reference package. 
- Bams. Given as a bamlist (i.e. a file with one row per sample, the path to the bam)

### Output <a name="paragraph-io-output"></a>

- Text files, explained here
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

For a full list of options, query `?QUILT::QUILT_HLA`, or alternatively, type 
```
./QUILT_HLA.R --help
```

## Preparing a reference package <a name="paragraph-preparing"></a>

An example of this is presented in detail in [example/QUILT_hla_reference_panel_construction.Md](example/QUILT_hla_reference_panel_construction.Md), which was used to make the reference panel package from 1000 Genomes Project data presented above.



