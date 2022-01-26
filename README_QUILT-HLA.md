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

QUILT-HLA is program for rapid HLA imputation from low-coverage sequence. QUILT-HLA uses information from reads within an HLA locus through direct read mapping, and uses information from imputation using a labelled haplotype reference panel using reads outside an HLA locus. QUILT-HLA is highly accurate, including at low coverage, and has been tested with coverages as low as 0.1X, with higher coverage slightly increasing accuracy through more reads being available to inform direct read mapping. Imputation results are most accurate for samples most closely matching the reference panel, but are generally high across populations. Imputation results using low-coverage sequence and QUILT-HLA generally outperform imputation from genotyping microarray input. Further details and detailed evaluations are available in the [QUILT paper](README.md#paragraph-citation).

## Installation <a name="paragraph-installation"></a>

QUILT-HLA is installed through the installation of QUILT. No additional installation is necessary.

## Quick start run <a name="paragraph-quickstartrun"></a>

A quick start to ensure QUILT-HLA is properly installed and working can be performed using the following.

First, download some prepared reference panel data. This reference panel package was prepared as described in [Preparing a reference package](#paragraph-preparing), and uses data from IPD-IGMT version 3.39, 1000 Genomes Project haplotypes (20201028), and 1000 Genomes Project HLA types (20181129).

<a><img src="important.png"/></a>
This example uses a reference panel data package with some of the reference samples excluded(e.g. NA12878). For normal use, use a panel without samples removed. See [Reference packages](#reference-packages) for more details.

```
wget http://www.stats.ox.ac.uk/~rdavies/QUILT_HLA_reference_package_demonstration_3.43_2021_12_26.tar ## or curl -O
tar -xvf QUILT_HLA_reference_package_demonstration_3.43_2021_12_26.tar
```

Download some example bams
```
wget http://www.stats.ox.ac.uk/~rdavies/QUILT_HLA_example_bams_2021_12_28.tar
tar -xvf QUILT_HLA_example_bams_2021_12_28.tar
```

HLA imputation for a particular region (here A) can be done as follows
```
HLA_GENE="A"
REF_DIR="quilt_hla_reference_panel_files_2021_12_26_demonstration_3.43"
./QUILT_HLA.R \
--outputdir=quilt_output \
--bamlist=bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=${REF_DIR} \
--quilt_hla_haplotype_panelfile=${REF_DIR}/quilt.hrc.hla.${HLA_GENE}.haplotypes.RData \
--dict_file=hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict
```

The above imputes three samples, with true HLA types as follows
```
Population Sample.ID HLA.A.1 HLA.A.2
CEU   NA12878   01:01   11:01
ASW   NA19625   02:01   23:17
ASW   NA19700   03:01   30:01
```

The output of QUILT-HLA will vary slightly run to run, as it uses random sampling, but should look approximately like the following (from file `quilt_output/quilt.hla.output.combined.all.txt`)
```
sample_number	sample_name	bestallele1	bestallele2	post_prob	sums
1	NA12878	A*01:01	A*11:01	0.999996657444772	0.999996657444772
2	NA19625	A*02:01	A*23:01	0.673515530204319	0.673515530204319
2	NA19625	A*02:01	A*23:17	0.326258885270228	0.999774415474547
3	NA19700	A*03:01	A*30:01	0.999999942895535	0.999999942895535
```

Here we see that NA12878 is correctly imputed, as is NA19700, while for NA19625, the first allele is imputed correctly (A\*02:01), while the second one is imputed less confidently, with the most confident allele (A\*23:01, confidence 0.67) is incorrect, while the second allele (A\*23:17, confidence 0.3262) is correct. More details about the output formats are given in [Output](#paragraph-io-output).

## Input and output formats <a name="paragraph-io"></a>

### Input <a name="paragraph-io-input"></a>

- Bams. Given as a bamlist (i.e. a file with one row per sample, the path to the bam)
- Reference package and ancillary files. For available pre-made options, see [Reference packages](#reference-packages). For more detail about building your own, see [Preparing haplotype files](#paragraph-preparing-haplotypes). For help in making a reference package, or any other questions or bugs, please feel free to email or file an issue on github.

### Output <a name="paragraph-io-output"></a>

Output is given as text files, with default names `quilt.hla.output.<combined/onlystates>.<all/topresult>.txt`. Here we begin by revisiting `quilt.hla.output.combined.all.txt`, explaining the output, and then explain the variations in the other files
```
sample_number	sample_name	bestallele1	bestallele2	lhoods	sums
1	NA12878	A*01:01	A*11:01	0.999996657444772	0.999996657444772
2	NA19625	A*02:01	A*23:01	0.673515530204319	0.673515530204319
2	NA19625	A*02:01	A*23:17	0.326258885270228	0.999774415474547
3	NA19700	A*03:01	A*30:01	0.999999942895535	0.999999942895535
```
The first two columns are straightforward, of `sample_number`, being the 1-based integer index of the sample from the original bamlist, and `sample_name`, being the sample name of that sample as taken from the BAM header. Next, we have the pair of imputed alleles, given in columns `bestallele1` and `bestallele2`. Finally, `post_prob` gives the posterior probability of the combination of alleles, and `sums` gives the sum of posterior probabilities of successive alleles (i.e. cumulative sum of post_prob across pairs of alleles). Pairs of alleles are outputted until the `sums` argument exceeds the QUILT-HLA parameter `summary_best_alleles_threshold` with default value `0.99`.

The difference between `all/topresult` is that `topresult` only outputs the single most likely pair of alleles for each sample, while `all` continues until the `sums` value exceeds the threshold as explained just above.

The difference between `combined/onlystates` is that `combined` uses information from all reads using both read mapping and imputation using a labelled haplotype reference panel, while `onlystates` uses only the imputation and not direct read mapping. As such you are recommended to use `combined` as your default file to use.

## Help, options and parameters <a name="paragraph-helpoptionsparams"></a>

For a full list of options, query `?QUILT::QUILT_HLA`, or alternatively, type 
```
./QUILT_HLA.R --help
```

## Reference packages <a name="reference-packages"></a>

Reference packages built

Max N=5132 haplotyes built using QUILT 0.1.5, IPD-IGMT version 3.39, 1000 Genomes Project haplotypes (20201028), and 1000 Genomes Project HLA types (20181129). 
```
http://www.stats.ox.ac.uk/~rdavies/QUILT_HLA_reference_package_2021_04_09.tar
```

Max N=5132 haplotyes built using QUILT 1.0.1, IPD-IGMT version 3.43, 1000 Genomes Project haplotypes (20201028), and 1000 Genomes Project HLA types (20181129). 
```
http://www.stats.ox.ac.uk/~rdavies/QUILT_HLA_reference_package_2021_12_24.tar
```


## Preparing a reference package <a name="paragraph-preparing"></a>

An example of this is presented in detail in [example/QUILT_hla_reference_panel_construction.Md](example/QUILT_hla_reference_panel_construction.Md), which was used to make the reference panel package from 1000 Genomes Project data presented above. This file can also be run non-interactively, using `bash example/run_example.sh example/QUILT_hla_reference_panel_construction.Md`, and run multiple times using `scripts/hla_prepare_workflow.sh`.
