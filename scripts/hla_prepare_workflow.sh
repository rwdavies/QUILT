set -e

QUILT_DIR=~/proj/QUILT/
cd "${QUILT_DIR}"

SAVE_DIR="/data/smew1/rdavies/quilt_hla_packages/"
    
set -e

output_date=2021_12_28

ipdigmt_version=3.39
ipdigmt_link=https://github.com/ANHIG/IMGTHLA/blob/032815608e6312b595b4aaf9904d5b4c189dd6dc/Alignments_Rel_3390.zip?raw=true

ipdigmt_version=3.43
ipdigmt_link=https://github.com/ANHIG/IMGTHLA/blob/3430/Alignments_Rel_3430.zip?raw=true


##
## change some things here, potentially
##

for i in $(seq 2 3)
do

    echo ========================= i=${i} ======================== 

    cd ${QUILT_DIR}
    script=example/reference_panel_builder.sh
    rm -f ${script}
    ./example/run_example.sh example/QUILT_hla_reference_panel_construction.Md ${script}

    if [ ${i} == 1 ]
    then
	##
	## WITH exclusion of 3 samples, for demonstration purposes
	## the default behaviour of the example script
	##
	exclude_set="demonstration"
    elif  [ ${i} == 2 ]
    then
	##
	## WITH exclusion of many samples, for benchmarking
	##
	exclude_set="benchmarking"
	where=`grep -n "exclude NA12878 and two ASW samples for example usage below" ${script} | cut -f1 --delimiter=":"`
	cat ${script} | 
	    awk '{if(NR=='${where}' + 3) {print ""}else {print $0}}' | 
	    awk '{if(NR=='${where}' + 4) {print ""}else {print $0}}' | 
	    awk '{if(NR=='${where}' + 5) {print ""}else {print $0}}' > ${script}.temp
	mv ${script}.temp ${script}
	## now replace the file too
	where=`grep -n "^exclude_sample_list" ${script} | cut -f1 --delimiter=":"`
	exclude_sample_list="~/proj/QUILT/hla_ancillary_files/hla_samples_to_exclude.txt"
	what=`echo ${where} i exclude_sample_list=${exclude_sample_list}`
	awk '{if(FNR=='${where}') {print "exclude_sample_list='${exclude_sample_list}'"} else {print $0}}' ${script} > ${script}.temp
	mv ${script}.temp ${script}
    elif  [ ${i} == 3 ]
    then
	##
	## NO exclusion of samples, for released use
	##
	exclude_set="full"
	## same as above, but no insertion of new file, use blank one
	where=`grep -n "exclude NA12878 and two ASW samples for example usage below" ${script} | cut -f1 --delimiter=":"`
	cat ${script} | 
	    awk '{if(NR=='${where}' + 3) {print ""}else {print $0}}' | 
	    awk '{if(NR=='${where}' + 4) {print ""}else {print $0}}' | 
	    awk '{if(NR=='${where}' + 5) {print ""}else {print $0}}' > ${script}.temp
	mv ${script}.temp ${script}
    fi

    ## change output date
    where=`grep -n "^output_date" ${script} | cut -f1 --delimiter=":"`
    what=`echo ${where} i output_date=${output_date}`    
    awk '{if(FNR=='${where}') {print "output_date='${output_date}'"} else {print $0}}' ${script} > ${script}.temp
    mv ${script}.temp ${script}
    ## change download link
    where=`grep -n "^ipdigmt_link" ${script} | cut -f1 --delimiter=":"`
    what=`echo ${where} i ipdigmt_link=${ipdigmt_link}`    
    awk '{if(FNR=='${where}') {print "ipdigmt_link='${ipdigmt_link}'"} else {print $0}}' ${script} > ${script}.temp
    mv ${script}.temp ${script}

    ## run!
    source ${script} ## AM HERE - DID THIS WORK? FINISH OFF THE REST!
    description=${exclude_set}_${ipdigmt_version}
    
    ## Once done, make tar-ball of required outputs
    ## grab directory specified in above
    cd `dirname ${reference_package_dir}`
    temp=`basename ${reference_package_dir}`
    temp2=`echo "${temp}/*.RData"`
    tar -cvf QUILT_HLA_reference_package_${description}_${output_date}.tar ${temp2}
    chmod 755 QUILT_HLA_reference_package_${description}_${output_date}.tar
    mkdir -p /data/smew1/rdavies/quilt_hla_packages/
    rsync -av QUILT_HLA_reference_package_${description}_${output_date}.tar ${SAVE_DIR}
    cd "${QUILT_DIR}"

    if [ ${i} == 1 ]
    then
	## also bams
	## inputs_dir=`cat ${script} | grep inputs_dir= | sed 's/inputs_dir=//g'`
	cd ${inputs_dir}
	cat bamlist.txt | xargs -l basename > bamlist2.txt
	mv bamlist2.txt bamlist.txt
	tar -cvf QUILT_HLA_example_bams_${output_date}.tar *2.0X*bam* bamlist.txt
v	chmod 755 QUILT_HLA_example_bams_${output_date}.tar
	rsync -av QUILT_HLA_example_bams_${output_date}.tar ${SAVE_DIR}
	cd "${QUILT_DIR}"	
    fi

    ## nuke directory that was being used
    rm -f ${script}
    rm -r -f "${inputs_dir}"
    
    
done

    







exit








