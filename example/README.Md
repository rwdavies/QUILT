## Example

For now there is only one example, a slightly larger version of the quick start example, with imputation performed over a more reasonable chunk size.

Here we first walk through downloading some data, preparing the reference package using a haplotype reference file, and imputing samples using the prepared reference package. 


## Downloading example data
Here we download some example data to use. Note that the code below is just a fancy way of ensuring it should download on either mac or linux.
```
path=http://www.stats.ox.ac.uk/~rdavies/QUILT_example_2021_01_14B.tgz
wget ${path} || true
curl -O || true
if [ ! -f basename ${path} ]
then
    exit 1
fi
```




## Prepare reference package using a haplotype reference file
In QUILT, the step of loading in and converting the reference data into a binary form is done separately from the imputation, to avoid performing this step multiple times when imputing a small number of samples. The principal output of this is a file at quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData, which is used by QUILT itself during the actual imputation that follows in the next step.

Note that QUILT requires a reference panel in the .hap and .legend format. One easy way to make these from a haplotype VCF is by using the `bcftools convert --haplegendsample` command.
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