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
In QUILT, the step of loading in and converting the reference data into a binary form is done separately from the imputation, to avoid performing this step multiple times when imputing a small number of samples. The principal output of this is a file at `quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData`, which is used by QUILT itself during the actual imputation that follows in the next step.

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




## Impute samples using the prepared reference package
This imputes NA12878 on this region, using three different bams, a haplotagged Illumina, an ONT, and an unlinked Illumina example, all downsampled to about 1.0X. You can try imputing from lower coverage using the other included bamlist in this example.

If this run completes succesfully, you will get a VCF at quilt_output/quilt.chr20.2000001.4000000.vcf.gz. Note that this can be changed with `--output_filename`, which is useful if you're running multiple jobs off of the same reference dataset built above (they can all be pointed at the same `outputdir` that has the reference data).

For normal operation, when you do not have high quality phased truth genotypes, you can omit posfile and phasefile, and normal imputation will be performed.

Here we have turned on `make_plots_block_gibbs`, which is very slow, but does make some fun graphs when you have truth data, which will be discussed below. You should disable this for non-explaratory jobs.

```
./QUILT.R \
    --outputdir=quilt_output \
    --chr=chr20 \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000 \
    --bamlist=package_2021_01_14B/bamlist.1.0.txt \
    --posfile=package_2021_01_14B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.posfile.txt \
    --phasefile=package_2021_01_14B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.phasefile.txt \
    --make_plots=TRUE \
    --make_plots_block_gibbs=TRUE
```    

 plot_attempt_to_reblock_snps
So make plots is the one I want (probably)
Need to re-figure out what is plot attempt to reblock SNPs


## Plots