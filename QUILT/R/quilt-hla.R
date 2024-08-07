#' @title QUILT_HLA
#' @param bamlist Path to file with bam file locations. File is one row per entry, path to bam files. Bam index files should exist in same directory as for each bam, suffixed either .bam.bai or .bai
#' @param region HLA region to be analyzed, for example A for HLA-A
#' @param dict_file Path to dictionary file for reference build
#' @param hla_gene_region_file For reference packages built after QUILT 1.0.2, this is not used. For older reference packages, this is needed, and is a path to file with gene boundaries. 4 columns, named Name Chr Start End, with respectively gene name (e.g. HLA-A), chromsome (e.g. chr6), and 1 based start and end positions of gene
#' @param outputdir What output directory to use. Otherwise defaults to current directory
#' @param summary_output_file_prefix Prefix for output text summary files
#' @param nCores How many cores to use
#' @param prepared_hla_reference_dir Output directory containing HLA reference material necessary for QUILT HLA
#' @param quilt_hla_haplotype_panelfile Prepared HLA haplotype reference panel file
#' @param final_output_RData_file Final output RData file path, if desired
#' @param write_summary_text_files Whether to write out final summary text files or not
#' @param nGibbsSamples How many QUILT Gibbs samples to perform
#' @param n_seek_iterations How many seek iterations to use in QUILT Gibbs sampling
#' @param quilt_seed When running QUILT Gibbs sampling, what seed to use, optionally
#' @param chr What chromosome, probably chr6 or maybe 6
#' @param quilt_buffer For QUILT Gibbs sampling, what buffer to include around center of gene
#' @param quilt_bqFilter For QUILT Gibbs sampling, do not consider sequence information if the base quality is below this threshold
#' @param summary_best_alleles_threshold When reporting results, give results until posterior probability exceeds this value
#' @param downsampleToCov For imputing states specifically using QUILT, what coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage. This is not used in the direct read mapping
#' @return Results in properly formatted version
#' @author Robert Davies
#' @export
QUILT_HLA <- function(
    bamlist,
    region,
    dict_file,
    hla_gene_region_file = NULL,
    outputdir = "",
    summary_output_file_prefix = 'quilt.hla.output',
    nCores = 1,
    prepared_hla_reference_dir = "",
    quilt_hla_haplotype_panelfile = "",
    final_output_RData_file = NA,
    write_summary_text_files = TRUE,
    overrideoutput = FALSE,
    nGibbsSamples = 15,
    n_seek_iterations = 3,
    quilt_seed = NA,
    chr = 'chr6',
    quilt_buffer = 500000,
    quilt_bqFilter = 10,
    summary_best_alleles_threshold = 0.99,
    downsampleToCov = 30
) {

    x <- as.list(environment())
    command_line <- paste0(
        "QUILT_HLA(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))
    print_message("Begin QUILT-HLA")    

    ##
    ## rwd: simon notes from his version of the code
    ## note that many of these variables are now obsolete
    ##
    ## NOTES:
    ## 1000 genomes types from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/README_20181129_HLA_types_full_1000_Genomes_Project_panel.txt
    ## CITE: Immune diversity sheds light on missing variation in worldwide genetic diversity panels. Abi-Rached L, Gouret P, Yeh J-H, Di Cristofaro J, Pontarotti P, Picard C, and Paganini J.
    
    ## function to take a likelihood matrix and output the possible ordered combinations (by decreasing likelihood) to total a threshhold thresh (default 99 percent coverage overall)
    
    ## function to run a single HLA region for a single individual bam file for region "A","B","C","DQB1","DRB1"
    
    ## outputs a lot of stuff but the key for maximally informed inference is hla_results$processed_output$combinedscaledlikelihoodmat
    
    ## however also stored are results from QUILT only, and from READS (in the HLA gene), only
    ## it would be a red flag if these are strongly discordant
    ## there is also likelihood at higher resolution potentially available e.g. 6-digit or 8-digit resolution from the reads, which likely is useful mainly at higher coverage
    ## recommended to give full path for the finaloutputfile file to avoid unexpected behaviour, otherwise it will save into outputdir if no path given, or into the starting directory using the flag "overrideoutput=T"
    ## bamfile is name of bamfile to be analysed (now bamlist)
    ## region is HLA region, must be "A","B","C","DQB1", or "DRB1"
    ## outputdir is where output files (intermediate and final) will be written to
    ## finaloutputfile is name of main output file
    ## overrideoutput is a flag and defines if finaloutputfile should be placed into directory from which this code is run if no full path is given
    

    ##
    ## validation
    ##
    if (!file.exists(bamlist)) {
        stop(paste0("Cannot find file:", bamlist))
    }
    

    
    ##
    ## initialization stuff
    ##
    check_samtools() ## needs samtools 1.10 or greater in PATH
    startdir <- getwd()
    if (is.na(final_output_RData_file)) {
        final_output_RData_file <- file_quilt_final_RData_output_file(outputdir, region)
    }

    ##
    ## other
    ##
    ##refseq_file <- file.path(prepared_hla_reference_dir, "refseq.txt")
    refseq_file <- dict_file
    if (!file.exists(refseq_file)) {
        stop(paste0("Cannot find file with sequence information refseq_file:", refseq_file))
    }

    
    ##
    ## Load various pre-made files
    ##
    print_message("Load input files")

    ##
    ## inputs from sequences
    ##
    load(file_quilt_hla_full(prepared_hla_reference_dir, region))
    load(file_quilt_hla_all_alleles_kmers(prepared_hla_reference_dir))
    load(file_quilt_hla_full_alleles_filled_in(prepared_hla_reference_dir, region))
    load(file_quilt_hla_full(prepared_hla_reference_dir, region))
    ## load(file.path(ancillary_file_dir, "HLAallalleleskmers.out")) ## from simon file, now incorporated
    ##load(file.path(ancillary_file_dir, paste0("HLA", region, "fullallelesfilledin.out"))) ## from simon file, now incorporated
    ## load(file.path(ancillary_file_dir, paste0("hla", region, "full.out"))) ## from simon file, now incorporated
    
    if (!exists("hla_gene_information")) {
        if (!file.exists(hla_gene_region_file)) {
            stop(paste0("An older reference package is being used, so you need to supply the file with gene boundaries. Cannot find file with HLA gene boundaries (hla_gene_region_file):", hla_gene_region_file))
        }
        ourpos2 <- read.table(hla_gene_region_file, header = TRUE)
        regstart <- ourpos2[ourpos2[, "Name"] == paste0("HLA-", region), "Start"]
        regend <- ourpos2[ourpos2[, "Name"] == paste0("HLA-", region), "End"]
        regmid <- (regstart + regend) / 2
    } else {
        regstart <- hla_gene_information[hla_gene_information[, "Name"] == paste0("HLA-", region), "Start"]
        regend <- hla_gene_information[hla_gene_information[, "Name"] == paste0("HLA-", region), "End"]        
        regmid <- (regstart + regend) / 2
    }
    




    ##
    ## things that depend on the haplotypes
    ##
    ## phaseallelesfile <- file.path(
    ##     prepared_hla_reference_dir,
    ##     paste0("hla", region, "haptypesexcludefivepops.out")
    ## )
    phaseallelesfile <- file_quilt_hla_phased_haplotypes(prepared_hla_reference_dir, region)    
    if (!file.exists(phaseallelesfile)) {
        stop(paste0("Cannot find phase alleles file ", phaseallelesfile))
    }
    ## load(file.path(ancillary_file_dir, phaseallelesfile)) ##
    load(phaseallelesfile) ##     
    print_message(paste0("The number of HLA reference haplotypes is:", length(hlahaptypes)))

    if (quilt_hla_haplotype_panelfile == "") {
        quilt_hla_haplotype_panelfile <- file_quilt_hla_panelfile(prepared_hla_reference_dir, region) 
    }
    if (!file.exists(quilt_hla_haplotype_panelfile)) {
        stop(paste0("Cannot find prepared HLA haplotype reference file, expecting:", quilt_hla_haplotype_panelfile))
    }






    ##
    ##
    ##
    ## regs <- c("A","B","C","DMA" , "DMB" , "DOA"  ,"DOB",  "DPA1", "DPA2", "DPB1" ,"DPB2" ,"DQA1" ,"DQA2" ,"DQB1" ,"DRA",  "DRB1", "DRB5", "E", "F" ,   "G",    "H", "J" ,   "K",    "L",   "MICA", "MICB", "P", "S",    "T" ,   "TAP1" ,"TAP2", "U",    "V",    "W")

    kk <- unique(kmers)
    names(kk) <- kk
    for(i in 1:length(kk)) {
        kk[i] <- paste(positions[kmers==names(kk)[i]],collapse=",")
    }

    ##
    ## do the normal QUILT part first. don't need to save this per-se
    ##
    dir.create(tempdir(), showWarnings = FALSE) ## shouldn't need to create?
    if (!dir.exists(tempdir())) {
        stop("R Cannot create a useable tempdir using tempdir(). Please investigate this")
    }

    
    outfile1 <- tempfile()
    
    QUILT(
        outputdir = tempdir(),
        chr = chr,
        regionStart = regstart,
        regionEnd = regend,
        buffer = quilt_buffer,
        bamlist = bamlist,
        prepared_reference_filename = quilt_hla_haplotype_panelfile,
        RData_objects_to_save = c("sampleNames", "final_set_of_results"),
        output_RData_filename = outfile1,
        nCores = nCores,
        record_interim_dosages = FALSE,
        bqFilter = quilt_bqFilter,
        nGibbsSamples = nGibbsSamples,
        n_seek_its = n_seek_iterations,
        gamma_physically_closest_to = regmid,
        seed = quilt_seed,
        hla_run = TRUE,
        verbose = FALSE,
        downsampleToCov = downsampleToCov
    )
    load(outfile1)
    unlink(outfile1)
    ## will return final_set_of_results
    ## kind of cumbersome, but meh
    

    ##
    ## run many samples through
    ##
    bamfiles <- scan(bamlist, what = "char", quiet = TRUE)
    all_results <- mclapply(
        1:length(bamfiles),
        mc.cores = nCores,
        FUN = quilt_hla_one_sample,
        final_set_of_results = final_set_of_results,
        bamfiles = bamfiles,
        chr = chr,
        region = region,
        regstart = regstart,
        regend = regend,
        regmid = regmid,
        refseq_file = refseq_file,
        newkmers = newkmers,
        lookup = lookup,
        revlookup = revlookup,
        fullalleles = fullalleles,
        kk = kk,
        quilt_hla_haplotype_panelfile = quilt_hla_haplotype_panelfile,
        quilt_seed = quilt_seed,
        quilt_buffer = quilt_buffer,
        quilt_bqFilter = quilt_bqFilter,
        nGibbsSamples = nGibbsSamples,
        hlahaptypes = hlahaptypes,
        summary_best_alleles_threshold = summary_best_alleles_threshold
    )

    check_mclapply_OK(all_results)

    ## save everything to a giant RData file
    save(
        all_results,
        file = final_output_RData_file
    )

    ##
    ## make table summaries of the most important information
    ##
    ## OK, can have QUILT only, or with mapping
    ## then, can either be based on single, final phased run
    ## all of them (Gibbs)
    ## or somehow un-phased? leave this for now
    ## either for single one (the final, phased one), or proper Gibbs
    ## OK, so should have 
    
    ## so per-sample, write out top results, for either option
    ## and/or all results, so that sums >= 0.99

    ##
    ## just QUILT, joint results, i.e. across all the Gibbs samples
    ##
    if (write_summary_text_files) {
        for(i in 1:4) {
            if (i == 1) {
                only_take_top_result <- FALSE
                what <- "joint_quiltresults"
                summary_suffix <- "onlystates.all.txt"
            } else if (i == 2) {
                only_take_top_result <- TRUE
                what <- "joint_quiltresults"
                summary_suffix <- "onlystates.topresult.txt"
            } else if (i == 3) {
                only_take_top_result <- FALSE
                what <- "joint_combinedresults"
                summary_suffix <- "combined.all.txt"
            } else if (i == 4) {
                only_take_top_result <- TRUE
                what <- "joint_combinedresults"
                summary_suffix <- "combined.topresult.txt"
            }
            summary_quilt_joint_results <- summarize_all_results_and_write_to_disk(
                all_results = all_results,
                sampleNames = sampleNames,
                what = what,
                outputdir = outputdir,
                summary_prefix = summary_output_file_prefix,
                summary_suffix = summary_suffix,
                only_take_top_result = only_take_top_result
            )
        }
    }

    ##
    ## also, write out text tables of most common output?
    ## and/or, explain output?
    ##

    print_message("Done QUILT-HLA")

}



