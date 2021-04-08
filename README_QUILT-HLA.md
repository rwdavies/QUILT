QUILT-HLA
=========

For general details including installation of QUILT, citation, versions, and changelog, see the main QUILT [README](https://github.com/rwdavies/QUILT).

# Table of contents
1. [Introduction](#paragraph-introduction)
2. [Installation](#paragraph-installation)
3. [Quick start run](#paragraph-quickstartrun)
4. [Input and output formats](#paragraph-io)
    1. [Input](#paragraph-io-input)
    2. [Output](#paragraph-io-output)
5. [Help, options and parameters](#paragraph-helpoptionsparams)
6. [Reference packages](#paragraph-reference-package)
7. [Preparing a reference package](#paragraph-preparing)


## Introduction <a name="paragraph-introduction"></a>

QUILT-HLA is program for rapid HLA imputation from low-coverage sequence. QUILT-HLA uses reads inside an HLA locus for direct read mapping, and uses the remaining reads for imputation using a labelled reference panel (Online Methods). QUILT-HLA is highly accurate across populations and depth of coverage, including coverages as low as 0.1X, though higher coverage increases the accuracy of imputation through the read mapping component of QUILT-HLA. Imputation of HLA types from lc-WGS and QUILT-HLA generally outperforms imputation from genotyping microarray input. Further details and detailed evaluations are available in the QUILT paper.

## Installation <a name="paragraph-installation"></a>

QUILT-HLA is installed through the installation of QUILT. No additional installation is necessary.

## Quick start run <a name="paragraph-quickstartrun"></a>

A quick start to ensure QUILT-HLA is properly installed and working can be performed using the following.

First, download some prepared reference panel data. This reference panel package was prepared as described in [Preparing a reference package](#paragraph-preparing), and uses data from IPD-IGMT version 3.39, 1000 Genomes Project haplotypes (20201028), and 1000 Genomes Project HLA types (20181129).

<span style="color:red">IMPORTANT:</span> This example uses a reference panel data package with a few example samples removed (e.g. NA12878). For normal use, use a panel without samples removed. See [Reference packages](#reference-packages) for more details.

```
wget something
tar -xzvf 
```

Next, download some example bam data
```
wget something
tar -xzvf 
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

This imputes three samples, with true HLA types of X. Here we take a look at the primary output file from the above run
```

```
Which should show that two samples (X and Y) are imputed correctly, while for the other, one allele is imputed correctly (X:), while the other allele is the second most likely.

## Input and output formats <a name="paragraph-io"></a>

### Input <a name="paragraph-io-input"></a>

- Bams. Given as a bamlist (i.e. a file with one row per sample, the path to the bam)
- Reference package. For more detail, see [Preparing haplotype files](#paragraph-preparing-haplotypes). For making a reference package, for any questions, please email or file a bug report.

### Output <a name="paragraph-io-output"></a>

- Text files, explained here
```
```

## Help, options and parameters <a name="paragraph-helpoptionsparams"></a>

For a full list of options, query `?QUILT::QUILT_HLA`, or alternatively, type 
```
./QUILT_HLA.R --help
```


## Reference packages <a name="reference-packages"></a>


## Preparing a reference package <a name="paragraph-preparing"></a>

An example of this is presented in detail in [example/QUILT_hla_reference_panel_construction.Md](example/QUILT_hla_reference_panel_construction.Md), which was used to make the reference panel package from 1000 Genomes Project data presented above.
