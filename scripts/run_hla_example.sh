set -e

if [ ! -f QUILT_HLA_reference_package_samples_excluded_2021_04_08.tar.gz ]
then
    wget https://www.stats.ox.ac.uk/~rdavies/QUILT_HLA_reference_package_samples_excluded_2021_04_08.tar.gz ## or curl -O
fi
tar -xvf QUILT_HLA_reference_package_samples_excluded_2021_04_08.tar.gz

if [ ! -f QUILT_HLA_example_bams_2021_12_28.tar ]
then
    wget https://www.stats.ox.ac.uk/~rdavies/QUILT_HLA_example_bams_2021_12_28.tar
fi
tar -xvf QUILT_HLA_example_bams_2021_12_28.tar

HLA_GENE="A"
./QUILT_HLA.R \
--outputdir=quilt_output \
--bamlist=bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=quilt_hla_reference_panel_files_2021_04_08 \
--quilt_hla_haplotype_panelfile=quilt_hla_reference_panel_files_2021_04_08/quilt.hrc.hla.${HLA_GENE}.haplotypes.RData \
--hla_gene_region_file=hla_ancillary_files/hlagenes.txt \
--dict_file=hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict


