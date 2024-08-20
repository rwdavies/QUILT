The genetic map file and recombination rate
=====

[Here](hg38)  we provide the genetic map files of GRCh38 coordinates using CEU recombination rate, but any population can be used, or a cross-population recombination rate. Also, we give instruction on how one can download the recombination rate of interested population and create the genetic map file that QUILT needs.


Download the genetic map file of build 37 for CEU population from the 1000 Genome Project.

```
inputs_dir="."
cd ${inputs_dir}
wget ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/CEU_omni_recombination_20130507.tar
tar -xvf CEU_omni_recombination_20130507.tar
```

In this particular case, we can use liftOver, to lift over the recombination rate from build 37, to build 38. We do that using a helper R [script](https://github.com/rwdavies/QUILT/blob/master/scripts/make_b38_recomb_map.R) available in the QUILT repository.

```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
R -f QUILT/scripts/make_b38_recomb_map.R --args ${inputs_dir} CEU 6
```
