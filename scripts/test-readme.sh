#!/bin/bash

set -e

get_url () {
    url=${1}
    one_if_curl_installed=`which curl | wc -l`
    one_if_wget_installed=`which wget | wc -l`
    if [ ${one_if_curl_installed} == 1 ]
    then
	    curl -L -O ${url}
    elif [ ${one_if_wget_installed} == 1 ]
    then
	    wget ${url}
    fi
}

get_url "https://zenodo.org/records/12786681/files/QUILT2_example_2024.tar.xz"

tar --xz -xf QUILT2_example_2024.tar.xz

outdir=quilt2_output && rm -rf $outdir

./QUILT2.R \
    --outputdir=$outdir \
    --chr=chr20 \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000 \
    --nGen=100 \
    --bamlist=package/bamlist.1.0.txt \
    --genetic_map_file=package/CEU-chr20-final.b38.txt.gz \
    --reference_vcf_file=package/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.noNA12878.vcf.gz \
    --save_prepared_reference=TRUE


## Nipt

./QUILT2.R \
    --prepared_reference_filename=quilt2_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData \
    --bamlist=package/bamlist.1.0.txt \
    --fflist=package/fflist.1.0.txt \
    --method=nipt \
    --chr=chr20 \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000 \
    --nGen=100 \
    --output_filename=quilt2_output/quilt2.nipt.vcf.gz
