#' @title QUILT_HLA
#' @param outputdir What output directory to use
#' @param bamfile Path to bamfile to analyze
#' @param region HLA region to be analyzed, for example A for HLA-A
#' @param rundir Where some files are, to be deprecated
#' @param ancillary_file_dir To be deprecated likely
#' @param finaloutputfile Final output file path
#' @param nGibbsSamples How many QUILT Gibbs samples to perform
#' @param n_seek_iterations How many seek iterations to use in QUILT Gibbs sampling
#' @param quilt_seed When running QUILT Gibbs sampling, what seed to use, optionally
#' @param chr What chromosome, probably chr6 or maybe 6
#' @param quilt_buffer For QUILT Gibbs sampling, what buffer to include around center of gene
#' @param quilt_bqFilter For QUILT Gibbs sampling, do not consider sequence information if the base quality is below this threshold
#' @return Results in properly formatted version
#' @author Robert Davies
#' @export
QUILT_HLA <- function(
    outputdir,
    bamfile,
    region,
    rundir,
    ancillary_file_dir,
    finaloutputfile = NA,
    overrideoutput = FALSE,
    nGibbsSamples = 15,
    n_seek_iterations = 3,
    quilt_seed = NA,
    chr = 'chr6',
    quilt_buffer = 500000,
    quilt_bqFilter = 10
) {

    x <- as.list(environment())
    command_line <- paste0(
        "QUILT_HLA(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))

    ## to be deprecated
    excludefivepops <- TRUE
    
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
    ## bamfile is name of bamfile to be analysed
    ## region is HLA region, must be "A","B","C","DQB1", or "DRB1"
    ## excludefivepops should be set to T for testing but F in general
    ## rundir is the location of the script and other files needed for running the code, including pre-processed data files
    ## outputdir is where output files (intermediate and final) will be written to
    ## finaloutputfile is name of main output file
    ## overrideoutput is a flag and defines if finaloutputfile should be placed into directory from which this code is run if no full path is given
    

    ##
    ## validation
    ##
    if (!file.exists(bamfile)) {
        stop(paste0("Cannot find file:", bamfile))
    }
    

    
    ##
    ## initialization stuff
    ##
    check_samtools() ## needs samtools 1.10 or greater in PATH
    startdir <- getwd()
    if (outputdir == "output") {
        outputdir <- file.path(rundir, outputdir)
    }
    if (is.na(finaloutputfile)) {
        finaloutputfile=paste(outputdir,"/","quiltandreadcombinedanalysis",bamfile,"hla",region,".out",sep="")
    }
    ##by default, write into the output directory, but can override to write into directory run script from 
    ## or if full path given will use that, this is safest
    if(!length(grep("/",finaloutputfile)) & overrideoutput==FALSE) {
        finaloutputfile=paste(outputdir,"/",finaloutputfile,sep="")
    }
    if(!length(grep("/",finaloutputfile)) & overrideoutput==TRUE) {
        finaloutputfile=paste(outputdir,"/",finaloutputfile,sep="")
    }
    setwd(rundir)
    load(file.path(ancillary_file_dir, paste("hla",region,"full.out",sep="")))
    load(file.path(ancillary_file_dir, "hlageneboundaries.out"))
    regstart <- ourpos[region,1]
    regend <- ourpos[region,2]
    regmid <- (regstart+regend)/2

    ##ids
    phaseallelesfile <- paste("hla",region,"haptypes.out",sep="")
    if (excludefivepops) {
        phaseallelesfile <- paste("hla",region,"haptypesexcludefivepops.out",sep="")
    }
    load(file.path(ancillary_file_dir, phaseallelesfile))
    print_message(paste0("The number of HLA reference haplotypes is:", length(hlahaptypes)))

    load(file.path(ancillary_file_dir, "HLAallalleleskmers.out"))
    load(file.path(ancillary_file_dir, paste("HLA",region,"fullallelesfilledin.out",sep="")))
    load(file.path(ancillary_file_dir, paste("hla",region,"full.out",sep="")))
    ## make scope so all functions can use
    ## lookup <<- lookup ## Robbie: what is this sorcery!
    ## revlookup <<- revlookup
    ## fullalleles <<- fullalleles
    ## kmers <<- kmers 
    ## ourpos<<-ourpos     
    ## positions <<-positions   
    ## xx <<-xx


    ## regs <- c("A","B","C","DMA" , "DMB" , "DOA"  ,"DOB",  "DPA1", "DPA2", "DPB1" ,"DPB2" ,"DQA1" ,"DQA2" ,"DQB1" ,"DRA",  "DRB1", "DRB5", "E", "F" ,   "G",    "H", "J" ,   "K",    "L",   "MICA", "MICB", "P", "S",    "T" ,   "TAP1" ,"TAP2", "U",    "V",    "W")
    kk <- unique(kmers)
    names(kk) <- kk
    for(i in 1:length(kk)) {
        kk[i] <- paste(positions[kmers==names(kk)[i]],collapse=",")
    }




    ##now all the reads that map alternatively to HLA
    ##refseq.txt contains all the names of alleles that can be mapped to
    this <- read.delim("refseq.txt",header=F)


    ##
    ## this appears to get reads that map to the region of interest on chr6
    ##
    that <- get_that(
        bamfile = bamfile,
        chr = chr,
        regstart = regstart,
        regend = regend
    )
    
    ##
    ## this appears to get reads that map to any of the HLA regions from the ref
    ##
    that2 <- get_that2(
        this = this,
        region = region,
        that = that,
        bamfile = bamfile
    )


    that <- filter_that(
        that = that,
        chr = chr,
        regstart = regstart,
        regend = regend
    )

    that2 <- filter_that2(
        that2 = that2,
        chr = chr,
        regstart = regstart,
        regend = regend
    )


    ## now have filtered matrices
    ## now more sophisticated filtering to remove reads mapping to alternative HLA alleles
    ## uses the set of ALL kmers found in all HLA alleles, to determine whether kmers can be found ordered in correct vs. other alleles
    out <- further_filter_that_and_that2(
        that = that,
        that2 = that2,
        newkmers = newkmers,
        region = region
    )
    that <- out[["that"]]
    that2 <- out[["that2"]]




    ## print("#########could worry have no reads; perhaps a dummy??")
    ## have at least two reads, both that and that2?
    ## still need to be careful later!
    ## now have both types of read, filter and process
    ## deal with any problems of conversion,w now, just make sure both have at least two reads; duplicate reads are not overcounted in the end
    if(nrow(that)==1) that <- rbind(that,that)
    if(nrow(that2)==1) that2 <- rbind(that2,that2)
    if(nrow(that)==0 & nrow(that2)>0) that <- that2[1:2,]
    if(nrow(that2)==0 & nrow(that)>0) that2 <- that[1:2,]


    ##
    ## don't entirely know what this code is doing, the first part of Simon's read mapping part
    ##
    out <- do_simon_read_stuff_with_that_and_that2(
        that = that,
        that2 = that2,
        lookup = lookup,
        revlookup = revlookup,
        fullalleles = fullalleles,
        regstart = regstart,
        regend = regend,
        region = region,
        kk = kk
    )
    readlikelihoodmat <- out[["readlikelihoodmat"]]
    readset1 <- out[["readset1"]]
    readset2 <- out[["readset2"]]
    pairedscores <- out[["pairedscores"]]
    readscaledlikelihoodmat <- out[["readscaledlikelihoodmat"]]
    fourdigitreadscaledlikelihoodmat <- out[["fourdigitreadscaledlikelihoodmat"]]
    intersectfourdigitreadscaledlikelihoodmat <- out[["intersectfourdigitreadscaledlikelihoodmat"]]
    intersectquiltscaledlikelihoodmat <- out[["intersectquiltscaledlikelihoodmat"]]
    intersectreadscaledlikelihoodmat <- out[["intersectreadscaledlikelihoodmat"]]
    combinedscaledlikelihoodmat <- out[["combinedscaledlikelihoodmat"]]
    combinedresults <- out[["combinedresults"]]
    mappingonlyresults <- out[["mappingonlyresults"]]
    overall <- out[["overall"]]
            

    

    final_set_of_results <- run_QUILT_as_part_of_QUILT_HLA(
        outputdir = outputdir,
        bamfile = bamfile,
        region = region,
        rundir = rundir,
        excludefivepops = excludefivepops,
        regstart = regstart,
        regend = regend,
        regmid = regmid,
        quilt_seed = quilt_seed,
        chr = chr,
        quilt_buffer = quilt_buffer,
        quilt_bqFilter = quilt_bqFilter
    )     

    ##
    ## really strong Robbie hack because I don't know nature of how below works
    ## basically, do 1:20 are Gibbs samples
    ##
    i_gibbs_sample <- 1
    
    for(i_gibbs_sample in 1:(nGibbsSamples + 1)) {

        print_message(paste0("Re-shaping QUILT output ", i_gibbs_sample, " / ", nGibbsSamples))
        if (i_gibbs_sample == (nGibbsSamples + 1)) {
            use_final_phase <- TRUE
            use_averaging <- TRUE
        } else {
            use_final_phase <- FALSE
            use_averaging <- FALSE
        }
        
        ##resset is the results from QUILT
        ##can use this to make likelihood for different 4-digit codes, see below processing
        ##QUILT already takes into account allele frequency in the population
        ##there are some 4-digit codes that are non-overlapping, and even some weight on unknown alleles
        ##now take reads in our region and process
        resset <- reshape_QUILT_output(
            final_set_of_results = final_set_of_results,
            hlahaptypes = hlahaptypes,
            use_final_phase = use_final_phase,
            i_gibbs_sample = i_gibbs_sample
        )
        ##
        ## store a list of raw data things
        ##
        quiltprobs <- resset
        colnames(quiltprobs) <- c("Allele_1_best_prob","Allele_2_best_prob","Summed_means")
        raw_output <- list(
            quiltprobs = quiltprobs,
            readlikelihoodmat = readlikelihoodmat,
            readset1 = readset1,
            readset2 = readset2,
            pairedscores = pairedscores,
            ndistinctfragments = nrow(pairedscores)
        )
        
        out <- reshape_and_filter_resset(
            resset = resset,
            region = region,
            use_averaging = use_averaging
        )
        newresset <- out[["newresset"]]        
        newresset2 <- out[["newresset2"]]
        newquiltprobs <- out[["newquiltprobs"]]
    
        newphasedquiltprobs <- newresset2
        quiltscaledlikelihoodmat <- matrix(1,nrow=nrow(newresset),ncol=nrow(newresset))
        rownames(quiltscaledlikelihoodmat) <- rownames(newphasedquiltprobs)
        colnames(quiltscaledlikelihoodmat) <- rownames(newphasedquiltprobs)
        quiltscaledlikelihoodmat <- quiltscaledlikelihoodmat*newphasedquiltprobs[,1]
        quiltscaledlikelihoodmat <- t(t(quiltscaledlikelihoodmat)*newphasedquiltprobs[,2])
        quiltscaledlikelihoodmat <- 0.5*(quiltscaledlikelihoodmat+t(quiltscaledlikelihoodmat))


        ## make a list of processed output from reads
        ## now if we have some interesting results, process read based inferences
        if(length(that) | length(that2)){
            if (i_gibbs_sample == 1) {
                ## this should depend on quilt only through alleles which are constant
                output_fourdigitreadscaledlikelihoodmat <- get_fourdigitreadscaledlikelihoodmat(
                    overall = overall,
                    newphasedquiltprobs = newphasedquiltprobs
                )
            }
            fourdigitreadscaledlikelihoodmat <- output_fourdigitreadscaledlikelihoodmat[["fourdigitreadscaledlikelihoodmat"]]
            mappingonlyresults <- getbestalleles(fourdigitreadscaledlikelihoodmat)
            vv2 <- output_fourdigitreadscaledlikelihoodmat[["vv2"]]
            readscaledlikelihoodmat <- output_fourdigitreadscaledlikelihoodmat[["readscaledlikelihoodmat"]]
            intersectreadscaledlikelihoodmat <- output_fourdigitreadscaledlikelihoodmat[["intersectreadscaledlikelihoodmat"]]
            ##
            ##intersection of this with quilt four digit inferences, summing the intersection, above (should sum to 1)
            ## 
            intersectfourdigitreadscaledlikelihoodmat=fourdigitreadscaledlikelihoodmat[names(vv2),names(vv2)]
            intersectfourdigitreadscaledlikelihoodmat=intersectfourdigitreadscaledlikelihoodmat/sum(intersectfourdigitreadscaledlikelihoodmat)
            ##
            ##intersection of quilt matrix with four digit inferences (rescaled to sum to one, ordered to match above)
            ## 
            intersectquiltscaledlikelihoodmat=quiltscaledlikelihoodmat[rownames(intersectfourdigitreadscaledlikelihoodmat),colnames(intersectfourdigitreadscaledlikelihoodmat)]
            intersectquiltscaledlikelihoodmat=intersectquiltscaledlikelihoodmat/sum(intersectquiltscaledlikelihoodmat)
            ##product of intersection matrices gives overall likelihood, rescaled to sum to one
            combinedscaledlikelihoodmat=intersectfourdigitreadscaledlikelihoodmat*intersectquiltscaledlikelihoodmat
            combinedscaledlikelihoodmat=combinedscaledlikelihoodmat/sum(combinedscaledlikelihoodmat)
            combinedresults <- getbestalleles(combinedscaledlikelihoodmat,0.99)
            ##for a matrix, make a little function to output allele pair probabilities, the answer
        }

        quiltresults <- getbestalleles(quiltscaledlikelihoodmat,0.99)
        processed_output <- list(
            newquiltprobs = newquiltprobs,
            newphasedquiltprobs = newphasedquiltprobs,
            quiltscaledlikelihoodmat = quiltscaledlikelihoodmat,
            readscaledlikelihoodmat = readscaledlikelihoodmat,
            intersectreadscaledlikelihoodmat = intersectreadscaledlikelihoodmat,
            fourdigitreadscaledlikelihoodmat = fourdigitreadscaledlikelihoodmat,
            intersectfourdigitreadscaledlikelihoodmat = intersectfourdigitreadscaledlikelihoodmat,
            intersectquiltscaledlikelihoodmat = intersectquiltscaledlikelihoodmat,
            combinedscaledlikelihoodmat = combinedscaledlikelihoodmat
        )

        ## now per-gibbs sample, save
        ## need this
        if (i_gibbs_sample == 1) {
            joint_quiltscaledlikelihoodmat <- quiltscaledlikelihoodmat
            if(length(that) | length(that2)){
                joint_combinedscaledlikelihoodmat <- combinedscaledlikelihoodmat
            }
        } else if (i_gibbs_sample <= nGibbsSamples) {
            joint_quiltscaledlikelihoodmat <-
                joint_quiltscaledlikelihoodmat + 
                quiltscaledlikelihoodmat
            if(length(that) | length(that2)){
                joint_combinedscaledlikelihoodmat <-
                    joint_combinedscaledlikelihoodmat + 
                    combinedscaledlikelihoodmat
            }
        }

    }

    ##
    ## normalize joint version
    ##
    joint_quiltscaledlikelihoodmat <- joint_quiltscaledlikelihoodmat / nGibbsSamples
    joint_quiltresults <- getbestalleles(joint_quiltscaledlikelihoodmat, 0.99)    
    if(length(that) | length(that2)){
        joint_combinedscaledlikelihoodmat <- joint_combinedscaledlikelihoodmat / nGibbsSamples
        joint_combinedresults <- getbestalleles(joint_combinedscaledlikelihoodmat, 0.99)    
    } else {
        joint_combinedresults <- NULL
    }




    
    ##
    ## robbie hacks
    ## try and do proper unphased version
    ## 
    ## first, for QUILT, this is easy. take two highest posterior probabilities
    quilt_unphased_probs <- newresset[, 3]
    y <- quilt_unphased_probs[order(-quilt_unphased_probs)]
    unphased_summary_quilt_only <- c(y[1:2], conf = sum(y[1:2]))
    ## now, for Simon's thing
    ## first, get intersection
    if (!is.null(intersectfourdigitreadscaledlikelihoodmat)) {
        ## work on intersection, and make Simon's thing unphased
        a <- rownames(intersectfourdigitreadscaledlikelihoodmat)
        ## great
        simon_map_phased_probs <- intersectfourdigitreadscaledlikelihoodmat[a %in% names(y), a %in% names(y)]
        simon_map_unphased_probs <- rowSums(simon_map_phased_probs) ## I think
        ## OK, now merge
        joint <- intersect(names(quilt_unphased_probs), names(simon_map_unphased_probs))
        simon_unphased_joint <- simon_map_unphased_probs[match(joint, names(simon_map_unphased_probs))]
        quilt_unphased_joint <- quilt_unphased_probs[match(joint, names(quilt_unphased_probs))]
        both_unphased <- simon_unphased_joint * quilt_unphased_joint
        both_unphased <- both_unphased / sum(both_unphased)
        ##
        y <- both_unphased[order(-both_unphased)]
        unphased_summary_both <- c(y[1:2], conf = sum(y[1:2]))
    } else {
        unphased_summary_both <- NULL
    }
    ## simon_unphased_joint[c("A*02:01", "A*23:01")]
    unphased_summary_quilt_only
    unphased_summary_both

    
    ##various of above terms are empty if there is no additional informartion from  reads in a gene (at low coverage)
    ##below is a function to predict hla combinations from the above
    ## save in finaloutfile
    ##if that and that2 are empty, what happens, is one question
    ##should just have empty (null) containers for the other things I think

    hla_results <- list(
        raw_output = raw_output,
        processed_output = processed_output,
        region = region,
        quiltresults = quiltresults,
        combinedresults = combinedresults,
        mappingonlyresults = mappingonlyresults,
        unphased_summary_quilt_only = unphased_summary_quilt_only,
        unphased_summary_both = unphased_summary_both,
        joint_quiltresults = joint_quiltresults,
        joint_combinedresults = joint_combinedresults
    )
    
    save(hla_results, file = finaloutputfile)


}

