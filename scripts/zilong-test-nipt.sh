## local run

gtime=/usr/bin/time
workdir=$1
rdata=$2
output=$3

cd $workdir

## mspbwt
${gtime} -v QUILT.R \
         --output_filename=${output} \
         --prepared_reference_filename=${rdata} \
         --bamlist="bamlist.nipt.txt" \
         --method="nipt" \
         --fflist="fflist.nipt.txt" \
         --posfile="pos.txt" \
         --phasefile="phasefile.txt" \
         --use_hapMatcherR=TRUE \
         --impute_rare_common=TRUE \
         --chr="chr20" \
         --regionStart=15281340 \
         --regionEnd=20023517 \
         --buffer=500000 \
         --zilong=FALSE \
         --use_mspbwt=TRUE \
         --mspbwtM=1 \
         --mspbwtL=10 \
         --Knew=400 \
         --Ksubset=400 \
         --nGibbsSamples=7 \
         --n_seek_its=3  &> ${output}.log
cd -
exit

${gtime} -v QUILT.R \
         --output_filename=${output} \
         --prepared_reference_filename=${rdata} \
         --bamlist="bamlist.nipt.txt" \
         --method="nipt" \
         --fflist="fflist.nipt.txt" \
         --posfile="pos.txt" \
         --phasefile="phasefile.txt" \
         --chr="chr20" \
         --regionStart=15281340 \
         --regionEnd=20023517 \
         --buffer=500000 \
         --nGibbsSamples=7 \
         --n_seek_its=3  &> ${output}.log

