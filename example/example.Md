## for now there is only one example, a slightly larger version of the quick start example
## imputation is performed in a more reasonable chunk size for this example
## you should be able to run this example non-interactively (./examples/examples.sh)
## this is also how the example plots were generated
## see scripts/prepare_example.R for code for how it was made
## note this includes not directly avaialble public data (the bams, though these all should be available for download through links not used here in their creation)


## structure of remaining example
## step 0 - download some example data
## step 1 - prepare a reference package using a haplotype reference file
## step 2 - impute samples using the prepared reference package


##
## step 0 - download some example data
##
```
os=`uname -a | cut -f 1 -d " "`
if [ "${os}" == "Darwin" ]
then
   download_command="curl -O"
elif [ "${os}" == "Linux" ]
then
    download_command="wget"
else
    echo "cannot figure out operating system for selecting download command, please run file interactively"
    exit 1
fi
## download package - this should download the files from this website including for example "bamlist.1.0.txt"
if [ ! -f QUILT_example_2021_01_14B.tgz ]
then
    ${download_command} http://www.stats.ox.ac.uk/~rdavies/QUILT_example_2021_01_14B.tgz ## either wget or curl -O or otherwise as appropriate
fi
tar -xzvf QUILT_example_2021_01_14B.tgz
```




##
## step 1 - prepare a reference package using a haplotype reference file
##
```
rm -r -f quilt_output
./QUILT_prepare_reference.R \
    --outputdir=quilt_output \
    --chr=chr20 \
    --nGen=100 \
    --reference_haplotype_file=package_2021_01_14B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.noNA12878.hap.gz \
    --reference_legend_file=package_2021_01_14B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.noNA12878.legend.gz \
    --genetic_map_file=package_2021_01_14B/CEU-chr20-final.b38.txt.gz \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000
```
## note that succesful completion of this run results in a file at quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData, which is used by QUILT itself during the actual imputation that follows in the next step
## also note that making a reference panel in .hap and .legend format from a haplotype VCF is facilitated by using bcftools convert --haplegendsample




## 
## step 2 - impute samples using the prepared reference package
##
```
./QUILT.R \
    --outputdir=quilt_output \
    --chr=chr20 \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000 \
    --bamlist=package_2021_01_14B/bamlist.1.0.txt \
    --posfile=package_2021_01_14B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.posfile.txt \
    --phasefile=package_2021_01_14B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.phasefile.txt
```    
## note the following
##  this imputes NA12878 on this region, using a haplotagged Illumina, an ONT, and an unlinked Illumina example
##  succesful completion of this run results in a VCF at quilt_output/quilt.chr20.2000001.4000000.vcf.gz (this location can be changed with --output_filename)
##  output format is described on the main README
##  you can try bamlist.0.1.txt for 0.1X versions of the above bams
##  for normal operation when you do not have high quality phased truth genotypes, you can omit posfile and phasefile, and sites will be imputed 

Plots here