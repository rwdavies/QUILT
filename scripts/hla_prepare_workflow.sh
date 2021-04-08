set -e

cd ~/proj/QUILT/

## specify directory
./example/run_example.sh example/QUILT_hla_reference_panel_construction.Md 

exit

## code to make larger exclusion list for testing 
awk '{if ((($2 == "ASW") || ($2 == "CEU") || ($2 == "CHB") || ($2 == "PJL") || ($2 == "PUR"))) {print $0}}'  ${test_dir}/20181129_HLA_types_full_1000_Genomes_Project_panel.txt | cut -f3 > ${test_dir}exclude_ref_samples_for_testing.txt



