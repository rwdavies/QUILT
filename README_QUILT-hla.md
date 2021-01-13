QUILT-HLA
=========

For general details including installation, version, and changelog, see the main QUILT [README](README.md).

## Preparing files

### Preparing IPD-IGMT files

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

### Preparing haplotype files
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
