set -e

cd ~/proj/QUILT/

set -e

## Run twice, to make the two packages


##
## WITH exclusion of samples
##
script=example/reference_panel_with_exclusion.sh
rm -f ${script}
./example/run_example.sh example/QUILT_hla_reference_panel_construction.Md ${script}
bash ${script}
## Once done, make tar-ball of required outputs
## grab directory specified in above

exit


/data/smew1/rdavies/quilt_hla_reference_panel_build_2021_04_08/

exit




## WITHOUT exclusion of samples
script=example/reference_panel_no_exclusion.sh
rm -f ${script}
./example/run_example.sh example/QUILT_hla_reference_panel_construction.Md ${script}
## remove lines about exclusion, replace with empty file
where=`grep -n "exclude NA12878 and two ASW samples for example usage below" ${script} | cut -f1 --delimiter=":"`
cat ${script} | 
    awk '{if(NR=='${where}') {print "echo  > ${test_dir}exclude_ref_samples_for_testing.txt"} else {print $0}}' | 
    awk '{if(NR=='${where}' + 1) {print ""}else {print $0}}' | 
    awk '{if(NR=='${where}' + 2) {print ""}else {print $0}}' | 
    awk '{if(NR=='${where}' + 3) {print ""}else {print $0}}' > ${script}.temp
mv ${script}.temp ${script}

## remove exclusion here
example/reference_panel_no_exclusion.sh
## Make tar-ball of required outputs


exit

## code to make larger exclusion list for testing 
awk '{if ((($2 == "ASW") || ($2 == "CEU") || ($2 == "CHB") || ($2 == "PJL") || ($2 == "PUR"))) {print $0}}'  ${test_dir}/20181129_HLA_types_full_1000_Genomes_Project_panel.txt | cut -f3 > ${test_dir}exclude_ref_samples_for_testing.txt



