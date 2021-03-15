##
## Code to test minimal functionality of QUILT_HLA
##
inputs_dir=/data/smew1/rdavies/quilt_hla_finalize_with_simon/inputs/
test_dir=/data/smew1/rdavies/quilt_hla_finalize_with_simon/HLA_TEST_2021_03_15/
mkdir -p ${test_dir}

## To download IPD-IGMT version 3.39, for example
cd ${test_dir}
wget https://github.com/ANHIG/IMGTHLA/blob/032815608e6312b595b4aaf9904d5b4c189dd6dc/Alignments_Rel_3390.zip?raw=true
mv Alignments_Rel_3390.zip?raw=true Alignments_Rel_3390.zip

## 
## big unknown dependency 
## 
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/hla*haptypesexcludefivepops.out ${test_dir}

##
## other dependencies, to clean up
##
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/hlageneboundaries.out  ${inputs_dir} ## can be made a text file, included in repo
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/refseq.txt ${inputs_dir}
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/hlagenes.txt ${inputs_dir}
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/hlauntyped*.excludefivepop.txt ${inputs_dir}


##
## normal dependencies
##
WELL_HRC_DIR=/well/davies/users/dcc832/single_imp/2020_06_25/ref_panels/
rsync -av rescompNew2:${WELL_HRC_DIR}hrc.chr6.hap.clean.gz ${inputs_dir}
rsync -av rescompNew2:${WELL_HRC_DIR}hrc.chr6.legend.clean.gz ${inputs_dir}
rsync -av rescompNew2:${WELL_HRC_DIR}hrc.chr6.samples.reheadered2 ${inputs_dir}
## recombination
WELL_RECOMB_DIR=/well/davies/shared/recomb/CEU/
rsync -av rescompNew2:${WELL_RECOMB_DIR}CEU-chr6-final.b38.txt.gz ${inputs_dir}
## bam file
rsync -av rescompNew2:/well/davies/shared/1000G/mhc_hla/NA12878.mhc.2.0.bam ${inputs_dir}



##
## Prepare mapping related files
## Decently slow, think some of Simon's code is inefficient, but meh, can probably release
##
rsync -av ~/proj/QUILT/quilt_hla_supplementary_info.txt . ## what is this - check it out
cd ~/proj/QUILT/
./QUILT_HLA_prepare_reference.R \
--outputdir=${test_dir} \
--ipd_igmt_alignments_zip_file=Alignments_Rel_3390.zip \
--quilt_hla_supplementary_info_file=quilt_hla_supplementary_info.txt


##
## Preparing haplotype files
## Also slow-ish, less than above, files are massive. Can be made faster using cleverness. Also faster on 1000G
##
HLA_GENE="A"
regionStart=29942554
regionEnd=29945741
./QUILT_prepare_reference.R \
--outputdir=${test_dir} \
--nGen=100 \
--reference_haplotype_file=${inputs_dir}hrc.chr6.hap.clean.gz \
--reference_legend_file=${inputs_dir}hrc.chr6.legend.clean.gz \
--reference_sample_file=${inputs_dir}hrc.chr6.samples.reheadered2 \
--chr=chr6 \
--regionStart=${regionStart} \
--regionEnd=${regionEnd} \
--buffer=500000 \
--genetic_map_file=${inputs_dir}CEU-chr6-final.b38.txt.gz \
--reference_exclude_samplelist_file=${inputs_dir}/hlauntyped${HLA_GENE}.excludefivepop.txt \
--output_file=quilt.hrc.chr6.hla.${HLA_GENE}.haplotypes.RData \
--region_exclude_file=${test_dir}hlagenes.txt \
--minRate=0.01

##
## Try running QUILT_HLA here
##
#echo -e "/well/davies/shared/1000G/mhc_hla/NA12878.mhc.2.0.bam\n" > bamlist.txt
cd ${test_dir}
echo -e ${inputs_dir}"NA12878.mhc.2.0.bam\n" > bamlist.txt
~/proj/QUILT/QUILT_HLA.R \
--bamlist=bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=${test_dir} \
--quilt_hla_haplotype_panelfile=${test_dir}quilt.hrc.chr6.hla.${HLA_GENE}.haplotypes.RData
