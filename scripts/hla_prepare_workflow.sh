set -e

cd ~/proj/QUILT/

set -e

##
## WITH exclusion of samples
##
output_date=2021_04_08B
script=example/reference_panel_with_exclusion.sh
rm -f ${script}
./example/run_example.sh example/QUILT_hla_reference_panel_construction.Md ${script}
bash ${script}
## Once done, make tar-ball of required outputs
## grab directory specified in above
reference_package_dir=`cat ${script} | grep reference_package_dir= | sed 's/reference_package_dir=//g'`
current_dir=`pwd`
cd ${reference_package_dir}
rm -f quilt.hrc.hla.all.haplotypes.RData
cd ..
temp=`basename ${reference_package_dir}`
temp2=`echo "${temp}/*.RData"`
tar -cvf QUILT_HLA_reference_package_samples_excluded_${output_date}.tar ${temp2}
chmod 755 QUILT_HLA_reference_package_samples_excluded_${output_date}.tar
rsync -av QUILT_HLA_reference_package_samples_excluded_${output_date}.tar ~/pub_html/
cd ${current_dir}

## also bams
inputs_dir=`cat ${script} | grep inputs_dir= | sed 's/inputs_dir=//g'`
current_dir=`pwd`
cd ${inputs_dir}
cat bamlist.txt | xargs -l basename > bamlist2.txt
mv bamlist2.txt bamlist.txt
tar -cvf QUILT_HLA_example_bams_${output_date}.tar *2.0X*bam* bamlist.txt
chmod 755 QUILT_HLA_example_bams_${output_date}.tar
rsync -av QUILT_HLA_example_bams_${output_date}.tar ~/pub_html/




##
## WITHOUT exclusion of samples
##
output_date=2021_04_09
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
## run here
bash example/reference_panel_no_exclusion.sh

reference_package_dir=`cat ${script} | grep reference_package_dir= | sed 's/reference_package_dir=//g'`
current_dir=`pwd`
cd ${reference_package_dir}
rm -f quilt.hrc.hla.all.haplotypes.RData
cd ..
temp=`basename ${reference_package_dir}`
temp2=`echo "${temp}/*.RData"`
tar -cvf QUILT_HLA_reference_package_${output_date}.tar ${temp2}
chmod 755 QUILT_HLA_reference_package_${output_date}.tar
rsync -av QUILT_HLA_reference_package_${output_date}.tar ~/pub_html/
cd ${current_dir}









exit








