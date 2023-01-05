cd ~/proj/QUILT/

## version: QUILT1 vs QUILT2
## method: full vs mspbwt
## objects: regular vs efficient

## version:QUILT1 method:full objects:regular
## version:QUILT2 method:full objects:regular
## version:QUILT2 method:full objects:efficient
## version:QUILT2 method:mspbwt objects:regular
## version:QUILT2 method:mspbwt objects:efficient
## version:QUILT2 method:zilong objects:regular
## version:QUILT2 method:zilong objects:efficient

versions=(QUILT1 QUILT2 QUILT2 QUILT2 QUILT2 QUILT2 QUILT2)
methods=(full full full mspbwt mspbwt zilong zilong)
objects=(regular regular efficient regular efficient regular efficient)

export REBUILD=TRUE

## argh - have to do installations manually

for i_method in $(seq 3 6) ## change as appropriate
do
    export VERSION=${versions[${i_method}]}
    export METHOD=${methods[${i_method}]}
    export OBJECT=${objects[${i_method}]}    
    export OUTPUT_LOG1=~/proj/QUILT/scratch/build_${VERSION}_${METHOD}_${OBJECT}.txt
    export OUTPUT_LOG2=~/proj/QUILT/scratch/run_${VERSION}_${METHOD}_${OBJECT}.txt    
    
    # if [ "${REBUILD}" == "TRUE" ] || [ ! -e ${OUTPUT_LOG1} ]
    # then
    # 	export WHAT="BUILD"
    # 	/usr/bin/time -vvv ./scripts/profile2.R 2>&1 | tee ${OUTPUT_LOG1}
    # fi

    if [ "${REBUILD}" == "TRUE" ] || [ ! -e ${OUTPUT_LOG2} ]
    then
	export WHAT="RUN"
	/usr/bin/time -vvv ./scripts/profile2.R 2>&1 | tee ${OUTPUT_LOG2}	
    fi

done

exit

cd ~/proj/QUILT/scratch/

echo --- Build time ---
grep "Elapsed" build*
echo --- Build RAM ---
grep "Maximum resident set size" build*
echo --- Run time ---
grep "Elapsed" run*
echo --- Run RAM ---
grep "Maximum resident set size" run*

