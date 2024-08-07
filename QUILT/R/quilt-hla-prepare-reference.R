#' @title QUILT_HLA_prepare_reference
#' @param outputdir What output directory to use
#' @param nGen Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure, where K is number of haplotypes.
#' @param hla_types_panel Path to file with 1000 Genomes formatted HLA types (see example for format details)
#' @param ipd_igmt_alignments_zip_file Path to zip file with alignments from IPD-IGMT (see README and example for more details)
#' @param ref_fasta Path to reference genome fasta 
#' @param refseq_table_file Path to file with UCSC refseq gene information (see README and example for more details)
#' @param full_regionStart When building HLA full reference panel, start of maximal region spanning all HLA genes. The 1-based position x is kept if regionStart <= x <= regionEnd
#' @param full_regionEnd As above, but end of maximal region spanning all HLA genes
#' @param buffer Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd
#' @param region_exclude_file File with regions to exclude from constructing the reference panel. Particularly useful for QUILT_HLA, where you want to exclude SNPs in the HLA genes themselves, so that reads contribute either to the read mapping or state inference. This file is space separated with a header of Name, Chr, Start and End, with Name being the HLA gene name (e.g. HLA-A), Chr being the chromosome (e.g. chr6), and Start and End are the 1-based starts and ends of the genes (i.e. where we don't want to consider SNPs for the Gibbs sampling state inference)
#' @param genetic_map_file Path to file with genetic map information, a file with 3 white-space delimited entries giving position (1-based), genetic rate map in cM/Mbp, and genetic map in cM. If no file included, rate is based on physical distance and expected rate (expRate)
#' @param reference_haplotype_file Path to reference haplotype file in IMPUTE format (file with no header and no rownames, one row per SNP, one column per reference haplotype, space separated, values must be 0 or 1)
#' @param reference_legend_file Path to reference haplotype legend file in IMPUTE format (file with one row per SNP, and a header including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele)
#' @param reference_sample_file Path to reference sample file (file with header, one must be POP, corresponding to populations that can be specified using reference_populations)
#' @param reference_exclude_samplelist_file File with one column of samples to exclude from final reference panels. These samples can be removed at one of two points, depending on reference_exclude_samples_for_initial_phasing. If reference_exclude_samples_for_initial_phasing is FALSE, then these samples, if present, are used for the initial phasing and allele assignment. If TRUE, then they are removed immediately
#' @param reference_exclude_samples_for_initial_phasing See help for reference_exclude_samplelist_file
#' @param all_hla_regions Character vector of all HLA genes for which IPD-IMGT files are available for 
#' @param hla_regions_to_prepare Character vector of HLA genes to prepare for imputation
#' @param chr What chromosome to run (probably chr6 or similar)
#' @param minRate Minimum recomb rate cM/Mb
#' @param nCores How many cores to use
#' @return Results in properly formatted version
#' @author Robert Davies
#' @export
QUILT_HLA_prepare_reference <- function(
    outputdir,
    nGen,
    hla_types_panel,
    ipd_igmt_alignments_zip_file,
    ref_fasta,
    refseq_table_file,
    full_regionStart,
    full_regionEnd,
    buffer,
    region_exclude_file,
    genetic_map_file = "",
    reference_haplotype_file,
    reference_legend_file,
    reference_sample_file,
    reference_exclude_samplelist_file = "",
    reference_exclude_samples_for_initial_phasing = FALSE,
    all_hla_regions = c('A','B','C','DMA','DMB','DOA','DOB','DPA1','DPA2','DPB1','DPB2','DQA1','DQA2','DQB1','DRA','DRB1','DRB3','DRB4','DRB5','E','F','G','HFE','H','J','K','L','MICA','MICB','N','P','S','TAP1','TAP2','T','U','V','W','Y'),
    hla_regions_to_prepare = c('A','B','C','DQA1','DQB1','DRB1'),
    chr = 'chr6',
    minRate = 0.1,
    nCores = 1
) {


    ## DEPRECATED
    ## ' @param hla_gene_region_file Path to file with gene boundaries. 4 columns, named Name Chr Start End, with respectively gene name (e.g. HLA-A), chromsome (e.g. chr6), and 1 based start and end positions of gene
    ## ' @param quilt_hla_supplementary_info_file Path to file with supplementary information about the genes, necessary for proper converstion. File is tab separated with header, with 3 columns. First (allele) is the allele that matches the reference genome. Second (genome_pos) is the position of this allele in the reference genome, finally the strand (strand) (options 1 or -1)

    ## INTEGRATED
    ## ref_fasta <- "/data/smew1/rdavies/quilt_hla_2021_12_24_3430//GRCh38_full_analysis_set_plus_decoy_hla.fa"
    ## refseq_table_file <- "hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz"
    
    x <- as.list(environment())
    command_line <- paste0(
        "QUILT_HLA_prepare_reference(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))

    system(paste0("rsync -av ", ipd_igmt_alignments_zip_file, " ", outputdir))
    system(paste0("cd ", outputdir, " && unzip -u ", basename(ipd_igmt_alignments_zip_file)))

    hla_gene_information <- get_hla_gene_information(
        table_file = refseq_table_file,
        all_hla_regions = all_hla_regions,
        chr = chr,
        what = "refseq"
    )

   
    make_and_save_hla_all_alleles_kmers(
        outputdir = outputdir,
        all_hla_regions = all_hla_regions,
        hla_gene_information = hla_gene_information,
        nCores = nCores
    )

    make_and_save_hla_files_for_imputation(
        outputdir = outputdir,
        hla_regions_to_prepare = hla_regions_to_prepare,
        hla_gene_information = hla_gene_information,
        ref_fasta = ref_fasta,
        nCores = nCores
    )

    ## this is somewhat distinct in flavour, uses above, but also HRC haplotypes
    ## could remove steps and put in here conceivably, or wrap up above
    regions <- hla_regions_to_prepare
    phase_hla_haplotypes(    
        outputdir = outputdir,
        chr = chr,
        hla_gene_information = hla_gene_information,
        full_regionStart = full_regionStart,
        full_regionEnd = full_regionEnd,
        buffer = buffer,
        regions = regions,
        region_exclude_file = region_exclude_file,
        reference_haplotype_file = reference_haplotype_file,
        reference_legend_file = reference_legend_file,
        reference_sample_file = reference_sample_file,
        reference_exclude_samplelist_file = reference_exclude_samplelist_file,
        reference_exclude_samples_for_initial_phasing = reference_exclude_samples_for_initial_phasing,
        hla_types_panel = hla_types_panel,
        genetic_map_file = genetic_map_file,
        minRate = minRate,
        nGen = nGen,
        nCores = nCores
    )

    print_message("Done preparing HLA reference")
    
}






if (1 == 0) {

    
    
    ##
    ## input files (for first set of functions)
    ##
    ## input files from IPD-IMGT version 3.39.0 needed, one per region AND for other HLA regions that are to be excluded
    ## A_gen.txt, B_gen.txt, C_gen.txt, DQA1_gen.txt,DRB1_gen.txt
    ##
    ## list of hla regions ("ourfiles") hlareglist.out (now a character vector)
    ## the above has entries for each HLA gene i.e. A,B,C,....DRB, also DQB, etc that is going to be removed as a possible confusion source for kmers from reads


    ##
    ## formerly hlapos.out now a tab delimited text file
    ##
    ## files hlaxxpos.out giving for each _gen file in the list "regs", the position of the first base relative to the reference, and which sequence matches the reference allele
    ## hlaxxpos.out stores in R format a 1x3 matrix matches
    ##  - whose rowname is the matching allele (actually - seems irrelevant?)
    ##  - first (now second) entry the corresponding row in the _gen.txt file
    ##  - second (now third) entry the position of the first column in the reference genome
    ##  - finally (now fourth) the strand e.g. DRB1*15:01:01:01   96 32589742   -1
    ## above can be constructed using the genome sequence but not included here. Included files will then only work for the current build of HRC and regions

    ## "~/proj/QUILT/quilt_hla_supplementary_info.txt"    
    
    ##
    ## output files from sequence
    ## 
    ## output files, first line is complete data and then with imputation of missing vals
    ## essential overall file: HLAallalleleskmers.out (now file_quilt_hla_all_alleles_kmers)
    ## interim regional files: hla*snpformatalleles.out (now file_quilt_hla_snpformatalleles)
    ## essential regional files: HLAAfullallelesfilledin.out hlaAfull.out (now file_quilt_hla_full_alleles_filled_in and file_quilt_hla_full)

    ##
    ## output files from haplotypes
    ## 
    ## essential phasing output files, some individuals end up unphased (important??): hlaAphased.out hlaAhaptypes.out  hlaAhaptypesexcludefivepops.out  
    ## is below from Robbie or me? Maybe Robbie?: hlauntypedA.excludefivepop.txt  hlauntypedA.exclude.txt
    ## possibly below are also given from Robbie: hrc.chr6.hla.A.hlatypedexcludefivepop.RData  hrc.chr6.hla.A.hlatyped.RData

    
    ## regions
    
    ## all kmers that exist
    ## record all kmers that exist

    ## load in regions (ourfiles)
    ## load("hlareglist.out")

    
    ## convert Simons previous thing into single file
    x <- t(sapply(c("A", "B", "C", "DQA1", "DQB1", "DRB1"), function(hla_gene) {
        load(paste0("hla", hla_gene, "pos.out"))
        return(matches)
    }))
    
    write.table(
        x,
        file = "~/proj/QUILT/quilt_hla_supplementary_info.txt",
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    
}


