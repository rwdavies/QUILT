## QUILT HLA reference panel construction example

Here code to generate the provided reference panel data package is provided and walked through.

### Specify working directories

Where you want to store the files as you work on them, and where the output files should go

```
inputs_dir=/data/smew1/rdavies/quilt_hla_2021_04_08/
test_dir=/data/smew1/rdavies/quilt_hla_2021_04_08/
reference_package_dir=/data/smew1/rdavies/quilt_hla_reference_panel_files_2021_04_08/
mkdir -p ${test_dir}
mkdir -p ${inputs_dir}
mkdir -p ${reference_package_dir}
```

### Reference haplotype data including HLA alleles database

Here we show how we get reference haplotype data, in this case 1000 Genomes project haplotypes, and re-format them
```
cd ${inputs_dir}
oneKG_vcf_name=CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.filtered.shapeit2-duohmm-phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${oneKG_vcf_name}
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${oneKG_vcf_name}.tbi
bcftools view --output-file oneKG.temp.vcf.gz --output-type z --min-alleles 2 --max-alleles 2 --types snps ${oneKG_vcf_name} chr6:25000000-34000000
tabix oneKG.temp.vcf.gz
bcftools convert --haplegendsample oneKG oneKG.temp.vcf.gz
reference_haplotype_file=${inputs_dir}oneKG.hap.gz
reference_legend_file=${inputs_dir}oneKG.legend.gz
reference_sample_file=${inputs_dir}oneKG.samples
## convert to uppercase, slight wording change
sed -i 's/sample population group sex/SAMPLE POP GROUP SEX/g' ${reference_sample_file}
```

Here we download the IPD-IGMT data
```
## To download IPD-IGMT version 3.39, for example
cd ${test_dir}
wget https://github.com/ANHIG/IMGTHLA/blob/032815608e6312b595b4aaf9904d5b4c189dd6dc/Alignments_Rel_3390.zip?raw=true
mv Alignments_Rel_3390.zip?raw=true Alignments_Rel_3390.zip
```

Here we download the database of HLA alleles
```
cd ${test_dir}
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt
```

Here we can optionally make an exclusion file of samples we want to remove from constructing the reference panel (for example if using them to test performance)
```
# exclude NA12878 and two ASW samples for example usage below
echo NA12878 > ${test_dir}exclude_ref_samples_for_testing.txt
echo NA19625 >> ${test_dir}exclude_ref_samples_for_testing.txt
echo NA19700 >> ${test_dir}exclude_ref_samples_for_testing.txt
```


### Recombination rate

Here we use a CEU recombination rate, but any population can be used, or a cross-population recombination rate

```
cd ${inputs_dir}
wget ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/CEU_omni_recombination_20130507.tar
tar -xvf CEU_omni_recombination_20130507.tar
```

In this particular case, we need to use liftOver, to lift over the recombination rate from build 37, to build 38. We do that using a helper R script available in the QUILT repository.

```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
R -f ~/proj/QUILT/scripts/make_b38_recomb_map.R --args ${inputs_dir} CEU 6

```

### Make the reference data package

Here we run a single large self-contained function to build the reference panel. There are two other dependencies used below not explained below (`hla_gene_region_file`, `quilt_hla_supplementary_info_file`), for which an explanation is given with the QUILT help i.e. by typing `./QUILT_HLA_prepare_reference.R --help`. You can get current versions of these that work for GRCh38 by cloning the QUILT repository, if not already done so, or downloading the individual files and changing the paths below as appropriate.

```
cd ~/proj/QUILT/ ## change to the directory where you've cloned QUILT
./QUILT_HLA_prepare_reference.R \
--outputdir=${reference_package_dir} \
--nGen=100 \
--hla_gene_region_file=hla_ancillary_files/hlagenes.txt \
--hla_types_panel=${test_dir}/20181129_HLA_types_full_1000_Genomes_Project_panel.txt \
--ipd_igmt_alignments_zip_file=${test_dir}Alignments_Rel_3390.zip \
--quilt_hla_supplementary_info_file=hla_ancillary_files/quilt_hla_supplementary_info.txt \
--full_regionStart=25587319 \
--full_regionEnd=33629686 \
--buffer=500000 \
--region_exclude_file=hla_ancillary_files/hlagenes.txt \
--genetic_map_file=${inputs_dir}CEU/CEU-chr6-final.b38.txt.gz \
--reference_haplotype_file=${reference_haplotype_file} \
--reference_legend_file=${reference_legend_file} \
--reference_sample_file=${reference_sample_file} \
--reference_exclude_samplelist_file=${test_dir}exclude_ref_samples_for_testing.txt \
--reference_exclude_samples_for_initial_phasing=FALSE \
--hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')" \
--nCores=6

```

### Test the reference data package

Note in the below we need the reference sequence dictionary for `GRCh38_full_analysis_set_plus_decoy_hla`, which is described elsewhere. It can be made using for example Picard CreateSequenceDictionary, of which the current version was used
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict
```

Here we download some example 1000 Genomes data, using links given in an index file. Since we're actually downloading CRAM files, we need the reference genome to decode and turn into bams
```
cd ${inputs_dir}
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
for sample in NA12878 NA19625 NA19700
do
    link=`grep ${sample} 1000G_2504_high_coverage.sequence.index | cut -f1`
    wget ${link}.crai
    ## here first we download the sequence, then subsample to approximately 2X (assuming 30X input)
    samtools view -b -T GRCh38_full_analysis_set_plus_decoy_hla.fa ${link} chr6:26000000-34000000 | samtools view -b -s 0.0667 > ${sample}.mhc.2.0X.bam
    samtools index ${sample}.mhc.2.0X.bam
done

```


Here we test it out by running QUILT_HLA
```
echo -e ${inputs_dir}"NA12878.mhc.2.0X.bam" > ${test_dir}bamlist.txt
echo -e ${inputs_dir}"NA19625.mhc.2.0X.bam" >> ${test_dir}bamlist.txt
echo -e ${inputs_dir}"NA19700.mhc.2.0X.bam" >> ${test_dir}bamlist.txt

HLA_GENE="A"
cd ~/proj/QUILT/
./QUILT_HLA.R \
--outputdir=${test_dir} \
--bamlist=${test_dir}bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=${reference_package_dir} \
--quilt_hla_haplotype_panelfile=${reference_package_dir}quilt.hrc.hla.${HLA_GENE}.haplotypes.RData \
--hla_gene_region_file=hla_ancillary_files/hlagenes.txt \
--dict_file=hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict
```

Here we test against known truth. Here two samples that are imputed correctly, while the third matches one allele correctly, while the other is the second most likely result
```
echo Truth
file=${test_dir}20181129_HLA_types_full_1000_Genomes_Project_panel.txt
cut -f1-5 ${file} | head -n1
cut -f1-5 ${file} | grep NA12878
cut -f1-5 ${file} | grep NA19625
cut -f1-5 ${file} | grep NA19700
echo QUILT-HLA imputed
cat ${test_dir}quilt.hla.output.combined.all.txt
```
