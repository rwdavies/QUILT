QUILT-HLA
=========

For general details including installation, version, and changelog, see the main QUILT [README](https://github.com/rwdavies/QUILT).

# Table of contents
1. [Quick start run](#paragraph-quickstart)
2. [Options and parameters](#paragraph-optionsparams)
3. [Preparing a reference package](#paragraph-preparing)
    1. [Preparing ancillary files](#paragraph-ancillary-files)
    2. [Preparing IPD-IGMT files](#paragraph-preparing-ipdigmt)
    3. [Preparing haplotype files](#paragraph-preparing-haplotypes)


## Quick start run <a name="paragraph-quickstart"></a>

Example here, is this still the right way to do it?
```
echo -e ${inputs_dir}"NA12878.mhc.2.0.bam" > bamlist.txt
echo -e ${inputs_dir}"NA18566.mhc.2.0.bam" >> bamlist.txt
HLA_GENE="A"
~/proj/QUILT/QUILT_HLA.R \
--outputdir=${test_dir} \
--bamlist=bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=${test_dir} \
--quilt_hla_haplotype_panelfile=${test_dir}quilt.hrc.hla.${HLA_GENE}.haplotypes.RData \
--hla_gene_region_file=hla_ancillary_files/hlagenes.txt \
--dict_file=hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict
```

## Options and parameters <a name="paragraph-optionsparams"></a>

For a full list of options, query `?QUILT::QUILT_HLA`, or alternatively, type 
```
./QUILT_HLA.R --help
```

## Preparing a reference package <a name="paragraph-preparing"></a>

Some preamble

### Preparing ancillary files <a name="paragraph-ancillary-files"></a>


List here and explain any other random file

Provided files are for GRCh38
quilt_hla_supplementary_info.txt
hlagenes.txt 


#### Dictionary file

Reference sequence dictionary. Described elsewhere. Made using for example Picard CreateSequenceDictionary. Version here from 
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict
```

### Preparing IPD-IGMT files <a name="paragraph-preparing-ipdigmt"></a>

Download IPD_IGMT files
```
## To download IPD-IGMT version 3.39, for example
wget https://github.com/ANHIG/IMGTHLA/blob/032815608e6312b595b4aaf9904d5b4c189dd6dc/Alignments_Rel_3390.zip?raw=true
mv Alignments_Rel_3390.zip?raw=true Alignments_Rel_3390.zip
```
Prepare mapping related files
```
./QUILT_HLA_prepare_reference.R \
--outputdir=/well/davies/users/dcc832/single_imp/HLA_TEST_2021_01_06/ \
--ipd_igmt_alignments_zip_file=Alignments_Rel_3390.zip \
--quilt_hla_supplementary_info_file=quilt_hla_supplementary_info.txt
```


### Preparing haplotype files <a name="paragraph-preparing-haplotypes"></a>
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

