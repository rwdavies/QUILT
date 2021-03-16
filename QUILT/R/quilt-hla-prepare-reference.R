#' @title QUILT_HLA_prepare_reference
#' @param outputdir What output directory to use
#' @param ipd_igmt_alignments_zip_file Path to zip file with alignments from IPD-IGMT (see README and example for more details)
#' @param quilt_hla_supplementary_info_file Path to file with supplementary information about the genes, necessary for proper converstion. File is tab separated without header, with 4 columns, with the following entries. First is the HLA gene name (e.g. A for HLA-A). Second is the correponding row in the _gen.txt IPD-IMGT file. Third is the position of the first column in the reference genome. Fourth is the strand (options 1 or -1).
#' @param all_hla_regions Character vector of all HLA genes for which IPD-IMGT files are available for 
#' @param hla_regions_to_prepare Character vector of HLA genes to prepare for imputation
#' @param full_reference_hap_file Path to file with full haplotype reference file 
#' @param local_reference_hap_file Path to file with haplotype reference file for HLA genes with asterisk for gene name
#' @return Results in properly formatted version
#' @author Robert Davies
#' @export
QUILT_HLA_prepare_reference <- function(
    outputdir,
    ipd_igmt_alignments_zip_file,
    quilt_hla_supplementary_info_file,
    all_hla_regions = c('A','B','C','DMA','DMB','DOA','DOB','DPA1','DPA2','DPB1','DPB2','DQA1','DQA2','DQB1','DRA','DRB1','DRB3','DRB4','DRB5','E','F','G','HFE','H','J','K','L','MICA','MICB','N','P','S','TAP1','TAP2','T','U','V','W','Y'),
    hla_regions_to_prepare = c('A','B','C','DQA1','DQB1','DRB1'),
    full_reference_hap_file = "",
    local_reference_hap_file = ""
) {

    x <- as.list(environment())
    command_line <- paste0(
        "QUILT_HLA_prepare_reference(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))

    
    ## outputdir = "~/proj/QUILT/hla-data/",
    ## hla_zip_file = "~/proj/QUILT/HLA_input_2021_01_04.zip"
    
    system(paste0("rsync -av ", ipd_igmt_alignments_zip_file, " ", outputdir))
    system(paste0("cd ", outputdir, " && unzip -u ", basename(ipd_igmt_alignments_zip_file)))




    
    ##
    ## this is the checking and naming for the phasing which happens below
    ##
    regions <- hla_regions_to_prepare
    ## conceptually more separate
    if (full_reference_hap_file == "") {
        full_reference_hap_file <- file.path(outputdir, "quilt.hrc.chr6.hla.all.haplotypes.RData")
        if (!file.exists(full_reference_hap_file)) {
            stop(paste0("Cannot find full haplotype reference file:", full_reference_hap_file))
        }
    }
    if (local_reference_hap_file == "") {
        local_reference_hap_file <- file.path(outputdir, paste0("quilt.hrc.chr6.hla.*.haplotypes.RData")) ## OK?
        for(region in regions) {
            x <- gsub("*", region, local_reference_hap_file, fixed = TRUE)
            if (!file.exists(x)) {
                stop(paste0("Cannot find specific haplotype reference file:", x))
            }
        }
    }
    hla_types_panel <- file.path(outputdir, "20181129_HLA_types_full_1000_Genomes_Project_panel.txt")



    
    
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

    if (1 == 0) {
        
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

    ## "~/proj/QUILT/quilt_hla_supplementary_info.txt"    
    supplementary_gene_info <- read.table(quilt_hla_supplementary_info_file)
    colnames(supplementary_gene_info) <- c("gene", "first_row", "genome_pos", "strand")
    
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

    ## I think this is what Simon means here
    hla_regions <- all_hla_regions
    regs <- hla_regions_to_prepare
    ourfiles <- hla_regions

    make_and_save_hla_all_alleles_kmers(
        outputdir = outputdir,
        ourfiles = ourfiles
    )
    
    make_and_save_hla_snpformatalleles(
        outputdir = outputdir,
        regs = regs,
        supplementary_gene_info = supplementary_gene_info
    )

    make_and_save_hla_full_alleles_filled_in (
        outputdir = outputdir,
        regs = regs,
        supplementary_gene_info = supplementary_gene_info
    )

    make_quilt_hla_full(
        outputdir = outputdir,
        regs = regs
    )

    ## this is distinct in flavour, uses above, but also HRC haplotypes
    phase_hla_haplotypes(    
        outputdir = outputdir,
        regions = regions,
        full_reference_hap_file = full_reference_hap_file,
        local_reference_hap_file = local_reference_hap_file,
        hla_types_panel = hla_types_panel
    )
        
    print_message("Done preparing HLA reference")
    
}




