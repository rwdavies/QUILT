#!/usr/bin/env bash

set -e

script_dir=`dirname "$0"`
cd "${script_dir}"/../

MARKDOWN_FILE="${1:-example/QUILT_usage.Md}"

## get start and end bits
lines_temp_file=$(mktemp)
grep -n "\`\`\`" example/README.Md > ${lines_temp_file}

n=`wc -l ${lines_temp_file} | cut -f1 --delimiter=" "`
n=`echo $((${n} / 2))`

script_temp_file=$(mktemp)
for i in $(seq 1 ${n})
do
    a1=`echo $((2*${i} - 1))`
    a2=`echo $((2*${i} - 0))`
    b1=`head -n ${a1} ${lines_temp_file}| tail -n 1 | cut -f1 --delimiter=":"`
    b2=`head -n ${a2} ${lines_temp_file}| tail -n 1 | cut -f1 --delimiter=":"`    
    awk '{if((NR > '${b1}') && (NR < '${b2}')) {print $0}}' ${MARKDOWN_FILE} >> ${script_temp_file}
done

## hack, get rid of extra print out things
grep -v "^\[" ${script_temp_file} > ${script_temp_file}.temp
mv ${script_temp_file}.temp ${script_temp_file}

bash -e ${script_temp_file}

## clean up
rm ${script_temp_file}
rm ${lines_temp_file}
