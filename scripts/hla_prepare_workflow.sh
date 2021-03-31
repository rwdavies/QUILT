set -e

cd ~/proj/QUILT/


##
## Code to test minimal functionality of QUILT_HLA
##
inputs_dir=/data/smew1/rdavies/quilt_hla_finalize_with_simon/inputs_2021_03_29/
test_dir=/data/smew1/rdavies/quilt_hla_finalize_with_simon/HLA_TEST_2021_03_29/
mkdir -p ${test_dir}
mkdir -p ${inputs_dir}

##
## other dependencies, now captured
##
## todo, explain in README
ls -lth ~/proj/QUILT/hla_ancillary_files/quilt_hla_supplementary_info.txt
ls -lth ~/proj/QUILT/hla_ancillary_files/hlagenes.txt
ls -lth ~/proj/QUILT/hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict 


##
## normal dependencies
##
## recombination
WELL_RECOMB_DIR=/well/davies/shared/recomb/CEU/
rsync -av rescompNew2:${WELL_RECOMB_DIR}CEU-chr6-final.b38.txt.gz ${inputs_dir}
## bam file
rsync -av rescompNew2:/well/davies/shared/1000G/mhc_hla/NA12878.mhc.2.0.bam* ${inputs_dir}
rsync -av rescompNew2:/well/davies/shared/1000G/mhc_hla/NA18566.mhc.2.0.bam* ${inputs_dir}

##
## 1000 Genomes based version
##
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



## HRC version
if [ 1 == 0 ]
then
    ## make input HRC smaller, easier that way
    ## do pretty hack, this is fine for now.
    WELL_HRC_DIR=/well/davies/users/dcc832/single_imp/2020_06_25/ref_panels/
    rsync -av rescompNew2:${WELL_HRC_DIR}hrc.chr6.hap.clean.gz ${inputs_dir}
    rsync -av rescompNew2:${WELL_HRC_DIR}hrc.chr6.legend.clean.gz ${inputs_dir}
    rsync -av rescompNew2:${WELL_HRC_DIR}hrc.chr6.samples.reheadered2 ${inputs_dir}
    ##gunzip -c hrc.chr6.legend.clean.gz | awk 'NR % 100000 == 0' 
    ##gunzip -c hrc.chr6.legend.clean.gz | grep -n 23141138
    ##gunzip -c hrc.chr6.legend.clean.gz | grep -n 36253842
    ## unnecessary but should make some next steps faster until re-structuring done
    cd ${inputs_dir}
    if [ ! -f hrc.chr6.legend.clean.gz ]
    then
	echo "Shrinking input files"
	gunzip -c hrc.chr6.legend.clean.gz | awk '{if((NR == 1) || (NR >= 350000 && NR <= 550000)) {print $0}}' | gzip > hrc.chr6.legend.clean.small.gz
	gunzip -c hrc.chr6.hap.clean.gz | awk '{if((NR >= 349999 && NR <= 549000)) {print $0}}' | gzip > hrc.chr6.hap.clean.small.gz
	echo "Done shrinking input files"
    fi
    reference_haplotype_file=${inputs_dir}hrc.chr6.hap.clean.small.gz
    reference_legend_file=${inputs_dir}hrc.chr6.legend.clean.small.gz
    reference_sample_file=${inputs_dir}hrc.chr6.samples.reheadered2
fi



##
## Download HLA database
##
## To download IPD-IGMT version 3.39, for example
cd ${test_dir}
wget https://github.com/ANHIG/IMGTHLA/blob/032815608e6312b595b4aaf9904d5b4c189dd6dc/Alignments_Rel_3390.zip?raw=true
mv Alignments_Rel_3390.zip?raw=true Alignments_Rel_3390.zip
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt


##
## make exclusion file with some samples that we want to test
## 
awk '{if ((($2 == "ASW") || ($2 == "CEU") || ($2 == "CHB") || ($2 == "PJL") || ($2 == "PUR"))) {print $0}}'  ${test_dir}/20181129_HLA_types_full_1000_Genomes_Project_panel.txt | cut -f3 > ${test_dir}exclude_ref_samples_for_testing.txt




##
## Do all prep here, single function
## Slow, but self-contained
##
cd ~/proj/QUILT/
./QUILT_HLA_prepare_reference.R \
--outputdir=${test_dir} \
--nGen=100 \
--hla_gene_region_file=hla_ancillary_files/hlagenes.txt \
--hla_types_panel=${test_dir}/20181129_HLA_types_full_1000_Genomes_Project_panel.txt \
--ipd_igmt_alignments_zip_file=${test_dir}Alignments_Rel_3390.zip \
--quilt_hla_supplementary_info_file=hla_ancillary_files/quilt_hla_supplementary_info.txt \
--full_regionStart=25587319 \
--full_regionEnd=33629686 \
--buffer=500000 \
--region_exclude_file=hla_ancillary_files/hlagenes.txt \
--genetic_map_file=${inputs_dir}CEU-chr6-final.b38.txt.gz \
--reference_haplotype_file=${reference_haplotype_file} \
--reference_legend_file=${reference_legend_file} \
--reference_sample_file=${reference_sample_file} \
--reference_exclude_samplelist_file=${test_dir}exclude_ref_samples_for_testing.txt \
--reference_exclude_samples_for_initial_phasing=FALSE \
--hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')" \
--nCores=6


##
## Try running QUILT_HLA here
##
#echo -e "/well/davies/shared/1000G/mhc_hla/NA12878.mhc.2.0.bam\n" > bamlist.txt
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

## check versus previous values
## seems to have worked? A bit less confident I think, but OK?
