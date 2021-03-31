QUILT-HLA
=========

For general details including installation, version, and changelog, see the main QUILT [README](https://github.com/rwdavies/QUILT).

# Table of contents
1. [Quick start run](#paragraph-quickstart)
2. [Options and parameters](#paragraph-optionsparams)
3. [Preparing a reference package](#paragraph-preparing)
    1. [Preparing ancillary files](#paragraph-ancillary-files)
    2. [Preparing IPD-IGMT files](#paragraph-preparing-ipdigmt)
    3. [Preparing haplotype files](#paragraph-preparing-haplotypes)


## Quick start run <a name="paragraph-quickstart"></a>

Example here, is this still the right way to do it?
```
echo -e ${inputs_dir}"NA12878.mhc.2.0.bam" > bamlist.txt
echo -e ${inputs_dir}"NA18566.mhc.2.0.bam" >> bamlist.txt
HLA_GENE="A"
~/proj/QUILT/QUILT_HLA.R \
--outputdir=${test_dir} \
--bamlist=bamlist.txt \
--region=${HLA_GENE} \
--prepared_hla_reference_dir=${test_dir} \
--quilt_hla_haplotype_panelfile=${test_dir}quilt.hrc.hla.${HLA_GENE}.haplotypes.RData \
--hla_gene_region_file=hla_ancillary_files/hlagenes.txt \
--dict_file=hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict
```

## Options and parameters <a name="paragraph-optionsparams"></a>

For a full list of options, query `?QUILT::QUILT_HLA`, or alternatively, type 
```
./QUILT_HLA.R --help
```

## Preparing a reference package <a name="paragraph-preparing"></a>

An example of this is presented in detail in [example/QUILT_hla_reference_panel_construction.Md](example/QUILT_hla_reference_panel_construction.Md), which was used to make the reference panel package from 1000 Genomes Project data presented above.



