## Example

For now there is only one example, a slightly larger version of the quick start example, with imputation performed over a more reasonable chunk size.

Here we first walk through downloading some data, preparing the reference package using a haplotype reference file, and imputing samples using the prepared reference package. 


## Downloading example data
Here we download some example data to use. Note that the code below is just a fancy way of ensuring it should download on either mac or linux.
```
path=http://www.stats.ox.ac.uk/~rdavies/QUILT_example_2021_01_15B.tgz
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
In QUILT, the step of loading in and converting the reference data into a binary form is done separately from the imputation, to avoid performing this step multiple times when imputing samples in separate runs. The principal output of this is a file at `quilt_output/RData/QUILT_prepared_reference.chr20.2000001.4000000.RData`, which is used by QUILT itself during the actual imputation that follows in the next step.

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

Here we have turned on `make_plots_block_gibbs`, which is very slow, but does make some fun graphs when you have truth data, which will be discussed below. You should disable this for non-exploratory jobs.

```
./QUILT.R \
    --outputdir=quilt_output \
    --chr=chr20 \
    --regionStart=2000001 \
    --regionEnd=4000000 \
    --buffer=500000 \
    --bamlist=package_2021_01_15B/bamlist.1.0.txt \
    --posfile=package_2021_01_15B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.posfile.txt \
    --phasefile=package_2021_01_15B/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.4000000.phasefile.txt \
    --make_plots=TRUE
```

While this runs, because the phasefile (with truth phased haplotypes) and posfile (companion file, see QUILT.R --help for details), accuracy measures are printed during the run. They look like the following
```
[2021-01-15 18:52:59] Final imputation dosage accuracy for sample NA12878HT, r2:0.974
[2021-01-15 18:52:59] Final phasing accuracy for sample NA12878HT, pse:0, disc(%):1.7%
```
The r2 is the squared correlation between the imputed dosages and truth genotypes, in the central region (without buffer). The phasing accuracy lists two measures. First, at truth sites that are heterozygous, it lists the proportion of imputed sites that are not heterozygous (`disc` i.e. disrepency), and second, among sites both truth and imputed heterozygous, it lists the phase switch error (`pse`, the number of phase switch errors) (note that consecutive phase switch errors are not counted).

## Plot to understand per-sample imputation performance

Turning on `make_plots` will generate some plots that can be informative about how the method is doing. Here we take a look at a few of these plots, and interpret them. Note that since the random number generator in R is intialized differently each run, unless you set the seed in QUILT, re-running the same command will generate slightly different plots each time.

First, plots are named `haps.{sample_name}.{regionName}_igs.{iGibbsSample}.it{iSeekIt}.gibbs.png` where assuming you set `regionStart` and `regionEnd`, then `regionName={chr}.{regionStart}.{regionEnd}`. The `{iGibbsSample}` is what number of Gibbs sample this is (the method runs `nGibbsSamples` number of Gibbs sample iterations, and a final additional one where phasing is performed, so `iGibbsSample=nGibbsSamples + 1` is the phasing iteration). The `{iSeekIt}` is how many times the iterative process of first performing Gibbs sampling using a reduced reference panel, and then getting a new version of that reduced reference panel using full haploid imputation, is performed.

Here is an example of such a plot. Below we'll walk through what the plot contains, then we'll interpret it, alongside a later plot.

TODO, add image

The plot is divided into sections. All sections are aligned vertically, by physical position.

The first two sections, "Imputed haplotype 1" and "Imputed haplotype 2", show the same information, but for the two imputed haplotypes. Within these subsections are first two plots with red and green lines, and then a large colourful rectangular region. The upper red and green lines represent a rolling difference between the dosages of that imputed haplotype and the two truth haplotypes, one at a time, where a flat line represents perfect agreement between the truth and imputed haplotype, and higher values off the line indicate more differences. To the left is an r2 correlation between the truth and imputed haplotype dosage, in the region excluding the buffer. Below these, the large colourful rectangle gives, for any given physical position, a stacked bar chart of what fraction of the posterior state probability is taken up by a given haplotype in the small haplotype reference panel used by the Gibbs sampler.

The next section gives read assignments. The top of the two regions is for reads assigned to haplotype 1, the second to haplotype 2. The reads are coloured by their truth labelling of blue = truth haplotype 1, and red = truth haplotype 2. The reads are represented as continuous horizontal lines from their start to their end, ignoring gaps, e.g. a paired end read would be a horizontal line starting at the start of the first read in the pair, and ending at the end of the second read in the pair. Reads with uncertain truth status (posterior probability for both truth haplotypes of more than 0.33 each) are coloured black.

Finally there is the Imputed genotypes section, which, similar to before, shows a running measure of how accurate the imputed dosage is compated to the truth dosage, with a correlation r2 on the left indicating squared correlation between imputed and truth dosage. Finally there is a scaled smoothed version of the recombination rate map, where values higher on the y axis indicate higher recombination rate. 

Putting this all together allows us to interpret results. 
