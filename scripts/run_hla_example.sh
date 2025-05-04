set -e

if [ ! -f QUILT_HLA_reference_package_2021_04_09.tar.gz ]
then
    wget https://zenodo.org/records/14891167/files/QUILT_HLA_reference_package_2021_04_09.tar.gz ## or curl -O
fi
tar -xvf QUILT_HLA_reference_package_2021_04_09.tar.gz

if [ ! -f QUILT_HLA_example_bams_2021_12_28.tar ]
then
    wget https://zenodo.org/records/14891167/files/QUILT_HLA_example_bams_2021_12_28.tar
fi
tar -xvf QUILT_HLA_example_bams_2021_12_28.tar

mkdir -p quilt_output
HLA_GENE="A"
./QUILT_HLA.R \
--outputdir=quilt_output \
--bamlist=bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=quilt_hla_reference_panel_files_2021_04_09 \
--quilt_hla_haplotype_panelfile=quilt_hla_reference_panel_files_2021_04_09/quilt.hrc.hla.${HLA_GENE}.haplotypes.RData \
--hla_gene_region_file=hla_ancillary_files/hlagenes.txt \
--dict_file=hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict


