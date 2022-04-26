## Example

This example is a slightly larger version of the quick start example, with imputation performed over a more reasonable chunk size.

Here we first walk through downloading some data, preparing the reference package using a haplotype reference file, and imputing samples using the prepared reference package. 


## Downloading example data
Here we download some example data to use. Note that the code below is just a fancy way of ensuring it should download on either mac or linux, otherwise throw an error.
```
path=https://www.stats.ox.ac.uk/~rdavies/QUILT_example_2021_01_15B.tgz
file=`basename ${path}`
if [ ! -f ${file} ]
then
    wget ${path} || true
    curl -O ${path} || true
    tar -xzvf ${file}
fi
if [ ! -f ${file} ]
then
    exit 1
fi
```




## Prepare reference package using a haplotype reference file
QUILT can be run in two ways, either jointly converting the haplotype reference data and performing imputation, or splitting those into two distinct steps. The later option is useful to avoid performing this conversion step multiple times when imputing samples in separate runs, especially if the reference panel is very large (tens of thousands of haplotypes, for example).

Here we consider the two step procedure. The principal output of the first step is a file at `quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData`, which is used by QUILT itself during the actual imputation that follows in the next step.

Note that QUILT requires a reference panel in the IMPUTE .hap and .legend format. One easy way to make these from a haplotype VCF is by using the `bcftools convert --haplegendsample` command.
```
rm -r -f quilt_output
./QUILT_prepare_reference.R \
    --outputdir=quilt_output \
    --chr=chr20 \
    --nGen=100 \
    --reference_haplotype_file=package_2021_01_15B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.noNA12878.hap.gz \
    --reference_legend_file=package_2021_01_15B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.noNA12878.legend.gz \
    --genetic_map_file=package_2021_01_15B/CEU-chr20-final.b38.txt.gz \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000
```




## Impute samples using the prepared reference package
This imputes NA12878 on this region, using three different bams, a haplotagged Illumina, an ONT, and an unlinked Illumina example, all downsampled to about 1.0X. You can try imputing from lower coverage using the other included bamlist in this example.

If this run completes succesfully, you will get a VCF at quilt_output/quilt.chr20.2000001.4000000.vcf.gz. Note that this can be changed with `--output_filename`, which is useful if you're running multiple jobs off of the same reference dataset built above (they can all be pointed at the same `outputdir` that has the reference data).

For normal operation, when you do not have high quality phased truth genotypes, you would omit posfile and phasefile.

```
./QUILT.R \
    --outputdir=quilt_output \
    --chr=chr20 \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000 \
    --bamlist=package_2021_01_15B/bamlist.1.0.txt \
    --posfile=package_2021_01_15B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.posfile.txt \
    --phasefile=package_2021_01_15B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.phasefile.txt
```

While this runs, because the phasefile (with truth phased haplotypes) and posfile (companion file, see QUILT.R --help for details), accuracy measures are printed during the run. They look like the following
```
[2021-01-15 18:52:59] Final imputation dosage accuracy for sample NA12878HT, r2:0.974
[2021-01-15 18:52:59] Final phasing accuracy for sample NA12878HT, pse:0, disc(%):1.7%
```
The r2 is the squared correlation between the imputed dosages and truth genotypes, in the central region (without buffer). The phasing accuracy lists two measures. First, at truth sites that are heterozygous, it lists the proportion of imputed sites that are not heterozygous (`disc` i.e. disrepency), and second, among sites both truth and imputed heterozygous, it lists the phase switch error (`pse`, the number of phase switch errors) (note that consecutive phase switch errors are not counted (these are assumed to be artefacts of the truth data)).

