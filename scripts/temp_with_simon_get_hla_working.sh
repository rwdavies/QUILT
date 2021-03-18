##
## TODO
## 1) Re-do from scratch, make sure it works
## 2) Capture other dependencies more properly so can rsync nothing into testdir
## 3) Re-facor making reference panel to be more efficient
## 

inputs_dir=/data/smew1/rdavies/quilt_hla_finalize_with_simon/inputs/
test_dir=/data/smew1/rdavies/quilt_hla_finalize_with_simon/HLA_TEST_2021_03_15B/
mkdir -p ${test_dir}

set -e

cd ~/proj/QUILT/


##
## Code to test minimal functionality of QUILT_HLA
##
inputs_dir=/data/smew1/rdavies/quilt_hla_finalize_with_simon/inputs_2021_03_18/
test_dir=/data/smew1/rdavies/quilt_hla_finalize_with_simon/HLA_TEST_2021_03_18/
mkdir -p ${test_dir}




## 
## not yet properly captured dependencies
## 
## now made!
## rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/hla*haptypesexcludefivepops.out ${test_dir}

##
## other dependencies, to clean up
##
## not 100% sure, why isn't this in repo
rsync -av rescompNew2:~/proj/QUILT/quilt_hla_supplementary_info.txt ~/proj/QUILT/
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/hlageneboundaries.out  ${test_dir} ## can be made a text file, included in repo
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/refseq.txt ${test_dir}
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/hlagenes.txt ${test_dir}
rsync -av rescompNew2:/well/davies/shared/1000G/robbie_files/hlauntyped*.excludefivepop.txt ${test_dir}
rsync -av /datas/muscovy/not-backed-up/shi/covid/hlauntyped.exclude.txt ${test_dir}



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
rsync -av rescompNew2:/well/davies/shared/1000G/mhc_hla/NA12878.mhc.2.0.bam* ${inputs_dir}

## make input HRC smaller, easier that way
## do pretty hack, this is fine for now.
##gunzip -c hrc.chr6.legend.clean.gz | awk 'NR % 100000 == 0' 
##gunzip -c hrc.chr6.legend.clean.gz | grep -n 23141138
##gunzip -c hrc.chr6.legend.clean.gz | grep -n 36253842

## unnecessary but should make some next steps faster until re-structuring done
gunzip -c hrc.chr6.legend.clean.gz | awk '{if((NR == 1) || (NR >= 350000 && NR <= 550000)) {print $0}}' | gzip > hrc.chr6.legend.clean.small.gz
gunzip -c hrc.chr6.hap.clean.gz | awk '{if((NR >= 349999 && NR <= 549000)) {print $0}}' | gzip > hrc.chr6.hap.clean.small.gz




##
## Download HLA database
##
## To download IPD-IGMT version 3.39, for example
cd ${test_dir}
wget https://github.com/ANHIG/IMGTHLA/blob/032815608e6312b595b4aaf9904d5b4c189dd6dc/Alignments_Rel_3390.zip?raw=true
mv Alignments_Rel_3390.zip?raw=true Alignments_Rel_3390.zip
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt








##
## Prepare panels
##
cd ~/proj/QUILT
regionStart=25587319
regionEnd=33629686
./QUILT_prepare_reference.R \
--outputdir=${test_dir} \
--nGen=100 \
--reference_haplotype_file=${inputs_dir}hrc.chr6.hap.clean.small.gz \
--reference_legend_file=${inputs_dir}hrc.chr6.legend.clean.small.gz \
--reference_sample_file=${inputs_dir}hrc.chr6.samples.reheadered2 \
--chr=chr6 \
--regionStart=${regionStart} \
--regionEnd=${regionEnd} \
--buffer=500000 \
--genetic_map_file=${inputs_dir}CEU-chr6-final.b38.txt.gz \
--reference_exclude_samplelist_file=${test_dir}/hlauntyped.exclude.txt \
--output_file=${tet_dir}quilt.hrc.chr6.hla.all.haplotypes.RData \
--minRate=0.01 &
##
## Preparing haplotype files, need to do all of them it seems?
## Clearly need to fix this, easily, internally. Should be do-able with above
##
for i in $(seq 0 5)
do
    if [ $i == 0 ]
    then
	HLA_GENE="A"; regionStart=29942554; regionEnd=29945741
    elif [ $i == 1 ]
    then
	HLA_GENE="B"; regionStart=31353367; regionEnd=31357155
    elif [ $i == 2 ]
    then
	HLA_GENE="C"; regionStart=31268257; regionEnd=31353367
    elif [ $i == 3 ]
    then
	HLA_GENE="DRB1"; regionStart=32578780; regionEnd=32589729
    elif [ $i == 4 ]
    then
	HLA_GENE="DQB1"; regionStart=32660035; regionEnd=32666603
    elif [ $i == 5 ]
    then
	HLA_GENE="DQA1"; regionStart=32637480; regionEnd=32643199
    fi && 
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
	--output_file=${test_dir}quilt.hrc.chr6.hla.${HLA_GENE}.haplotypes.RData \
	--region_exclude_file=${test_dir}hlagenes.txt \
	--minRate=0.01 &
done
wait





##
## Prepare mapping related files
## Decently slow, think some of Simon's code is inefficient, but meh, can probably release
##
cd ~/proj/QUILT/
./QUILT_HLA_prepare_reference.R \
--outputdir=${test_dir} \
--ipd_igmt_alignments_zip_file=${test_dir}Alignments_Rel_3390.zip \
--quilt_hla_supplementary_info_file=quilt_hla_supplementary_info.txt \
--full_reference_hap_file=${test_dir}quilt.hrc.chr6.hla.all.haplotypes.RData \
--local_reference_hap_file=${test_dir}quilt.hrc.chr6.hla.*.haplotypes.RData \
--hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')"




##
## Try running QUILT_HLA here
##
#echo -e "/well/davies/shared/1000G/mhc_hla/NA12878.mhc.2.0.bam\n" > bamlist.txt

echo -e ${inputs_dir}"NA12878.mhc.2.0.bam" > bamlist.txt
HLA_GENE="A"
~/proj/QUILT/QUILT_HLA.R \
--outputdir=${test_dir} \
--bamlist=bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=${test_dir} \
--quilt_hla_haplotype_panelfile=${test_dir}quilt.hrc.chr6.hla.${HLA_GENE}.haplotypes.RData

## seems to have worked? A bit less confident I think, but OK?
