#!/usr/bin/env Rscript

if (!suppressPackageStartupMessages(require("optparse")))
    install.packages("optparse", repos="http://cran.rstudio.com/")

option_list <- list(
    make_option(
        "--outputdir",
        type = "character",
        help = "What output directory to use"
    ), 
    make_option(
        "--ipd_igmt_alignments_zip_file",
        type = "character",
        help = "Path to zip file with alignments from IPD-IGMT (see README and example for more details)"
    ), 
    make_option(
        "--quilt_hla_supplementary_info_file",
        type = "character",
        help = "Path to file with supplementary information about the genes, necessary for proper converstion. File is tab separated without header, with 4 columns, with the following entries. First is the HLA gene name (e.g. A for HLA-A). Second is the correponding row in the _gen.txt IPD-IMGT file. Third is the position of the first column in the reference genome. Fourth is the strand (options 1 or -1)."
    ), 
    make_option(
        "--all_hla_regions",
        type = "character",
        help = "Character vector of all HLA genes for which IPD-IMGT files are available for [default c('A','B','C','DMA','DMB','DOA','DOB','DPA1','DPA2','DPB1','DPB2','DQA1','DQA2','DQB1','DRA','DRB1','DRB3','DRB4','DRB5','E','F','G','HFE','H','J','K','L','MICA','MICB','N','P','S','TAP1','TAP2','T','U','V','W','Y')] ",
        default = "c('A','B','C','DMA','DMB','DOA','DOB','DPA1','DPA2','DPB1','DPB2','DQA1','DQA2','DQB1','DRA','DRB1','DRB3','DRB4','DRB5','E','F','G','HFE','H','J','K','L','MICA','MICB','N','P','S','TAP1','TAP2','T','U','V','W','Y')"
    ), 
    make_option(
        "--hla_regions_to_prepare",
        type = "character",
        help = "Character vector of HLA genes to prepare for imputation [default c('A','B','C','DQA1','DQB1','DRB1')] ",
        default = "c('A','B','C','DQA1','DQB1','DRB1')"
    ), 
    make_option(
        "--full_reference_hap_file",
        type = "character",
        help = "Path to file with full haplotype reference file [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--local_reference_hap_file",
        type = "character",
        help = "Path to file with haplotype reference file for HLA genes with asterisk for gene name [default \"\"] ",
        default = ""
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(QUILT))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
QUILT_HLA_prepare_reference(
    outputdir = opt$outputdir,
    ipd_igmt_alignments_zip_file = opt$ipd_igmt_alignments_zip_file,
    quilt_hla_supplementary_info_file = opt$quilt_hla_supplementary_info_file,
    all_hla_regions = eval(parse(text=opt$all_hla_regions)),
    hla_regions_to_prepare = eval(parse(text=opt$hla_regions_to_prepare)),
    full_reference_hap_file = opt$full_reference_hap_file,
    local_reference_hap_file = opt$local_reference_hap_file
)
