quilt_hla_one_sample <- function(
    iSample,
    final_set_of_results,
    bamfiles,
    chr,
    region,
    regstart,
    regend,
    regmid,
    refseq_file,
    newkmers,
    lookup,
    revlookup,
    fullalleles,
    kk,
    quilt_hla_haplotype_panelfile,
    quilt_seed,
    quilt_buffer,
    quilt_bqFilter,
    nGibbsSamples,
    hlahaptypes,
    summary_best_alleles_threshold
) {

    ##
    ## print 
    ##
    n <- length(bamfiles)
    if (iSample %in% (match(1:10, ceiling(10 * (1:n / n))))) {
        print_message(paste0("Processing file ", iSample, " out of ", n))
    }
    
    bamfile <- bamfiles[iSample]

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
    ## refseq.txt contains all the names of alleles that can be mapped to
    this <- read.delim(refseq_file, header = FALSE)
    
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
        regend = regend,
        region = region
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
    ## print_message("Interpret mapped sequence data")
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
            

    ##
    ## really strong Robbie hack because I don't know nature of how below works
    ## basically, do 1:20 are Gibbs samples
    ##
    i_gibbs_sample <- 1
    
    for(i_gibbs_sample in 1:(nGibbsSamples + 1)) {

        ## print_message(paste0("Re-shaping QUILT output ", i_gibbs_sample, " / ", nGibbsSamples))
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
            iSample = iSample,
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
            mappingonlyresults <- getbestalleles(fourdigitreadscaledlikelihoodmat, summary_best_alleles_threshold)
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
            combinedresults <- getbestalleles(combinedscaledlikelihoodmat, summary_best_alleles_threshold)
            ##for a matrix, make a little function to output allele pair probabilities, the answer
        }

        quiltresults <- getbestalleles(quiltscaledlikelihoodmat, summary_best_alleles_threshold)
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
    joint_quiltresults <- getbestalleles(joint_quiltscaledlikelihoodmat, summary_best_alleles_threshold)
    if(length(that) | length(that2)){
        joint_combinedscaledlikelihoodmat <- joint_combinedscaledlikelihoodmat / nGibbsSamples
        joint_combinedresults <- getbestalleles(joint_combinedscaledlikelihoodmat, summary_best_alleles_threshold)
    } else {
        ##joint_combinedresults <- NULL
        ## change this, so if no reads, this just gives top QUILT results
        joint_combinedresults <- joint_quiltresults
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
    ## unphased_summary_quilt_only
    ##unphased_summary_both

    
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

    return(hla_results)
    
}





## final_set_of_results <- run_QUILT_as_part_of_QUILT_HLA(
##     bamfile = bamfile,
##     quilt_hla_haplotype_panelfile = quilt_hla_haplotype_panelfile,
##     regstart = regstart,
##     regend = regend,
##     regmid = regmid,
##     quilt_seed = quilt_seed,
##     chr = chr,
##     quilt_buffer = quilt_buffer,
##     quilt_bqFilter = quilt_bqFilter
## )     


## run_QUILT_as_part_of_QUILT_HLA <- function(
##     bamfile,
##     quilt_hla_haplotype_panelfile,
##     regstart,
##     regend,
##     regmid,
##     quilt_seed = NA,
##     chr = "chr6",
##     quilt_buffer = 500000,
##     quilt_bqFilter = 10
## ) {
##     ## imputation results from QUILT
##     ## outfile1 <- paste(outputdir,"/quiltoutput.",bamfile,".hla.",region,".RData",sep="")
##     outfile1 <- tempfile()
##     ## imputation results from reads and database
##     bamlist <- tempfile()
##     cat(bamfile, file = bamlist)
##     ## run QUILT
##     ## quilt_cmd <- paste0(
##     ##     "cd ",bamdir,"; ",
##     ##     "R -f ", quilt_imp_file, "  --args chr6 ",regstart," ",regend," 500000 NA12878 ,",bamfile," ",quiltpanelfile," ",outfile1," 1 /well/davies/users/dcc832/single_imp/2020_06_25/genotypes/gen.NA12878.chr20.1.2000000.txt  1.0 simon_is_awesome FALSE 10 FALSE ", nGibbsSamples, " ", n_seek_iterations, " full ",regmid," ", quilt_seed
##     ## )
##     ## genfile <- "/well/davies/users/dcc832/single_imp/2020_06_25/genotypes/gen.NA12878.chr20.1.2000000.txt"
##     QUILT(
##         outputdir = tempdir(),
##         chr = chr,
##         regionStart = regstart,
##         regionEnd = regend,
##         buffer = quilt_buffer,
##         bamlist = bamlist,
##         prepared_reference_filename = quilt_hla_haplotype_panelfile,
##         RData_objects_to_save = "final_set_of_results",
##         output_RData_filename = outfile1,
##         nCores = 1,
##         record_interim_dosages = FALSE,
##         bqFilter = quilt_bqFilter,
##         nGibbsSamples = nGibbsSamples,
##         n_seek_its = n_seek_iterations,
##         gamma_physically_closest_to = regmid,
##         seed = quilt_seed,
##         hla_run = TRUE
##     )
##     ##
##     load(outfile1)
##     unlink(outfile1)
##     ## 
##     return(final_set_of_results)
## }




check_samtools <- function() {
    ## /data/smew2/myers/1000G/samtools-1.10/
    v <- strsplit(system("samtools --version", intern = TRUE)[1], " ", fixed = TRUE)[[1]][2]
    vs <- as.numeric(strsplit(v, ".", fixed = TRUE)[[1]])
    ok <- vs[1] >= 2 | (vs[1] == 1 & vs[2] >= 10)
    if (!ok) {
        stop("Need samtools >=1.10 in PATH")
    }
    return(NULL)
}





reshape_QUILT_output <- function(
    final_set_of_results,
    iSample,
    hlahaptypes,
    use_final_phase,
    i_gibbs_sample
) {
    ## OK here's where Simons starting on this
    if (use_final_phase) {
        g1 <- final_set_of_results[[iSample]]$gamma1
        g2 <- final_set_of_results[[iSample]]$gamma2
    } else {
        g1 <- final_set_of_results[[iSample]]$list_of_gammas[[i_gibbs_sample]][[1]]
        g2 <- final_set_of_results[[iSample]]$list_of_gammas[[i_gibbs_sample]][[2]]
    }
    g3 <- final_set_of_results[[iSample]]$gamma_total    
    names(g1) <- hlahaptypes
    names(g2) <- hlahaptypes
    names(g3) <- hlahaptypes
    ##tabulate
    uniques=unique(names(g1))
    resset=matrix(0,nrow=length(uniques),ncol=3)
    rownames(resset)=uniques
    for(i in 1:length(g1)) resset[names(g1)[i],1]=as.double(resset[names(g1)[i],1])+g1[i]
    for(i in 1:length(g2)) resset[names(g2)[i],2]=as.double(resset[names(g2)[i],2])+g2[i]
    for(i in 1:length(g3)) resset[names(g3)[i],3]=as.double(resset[names(g3)[i],3])+g3[i]
    ## print(resset[order(resset[,1],decreasing=T),][1:10,])
    ## print(resset[order(resset[,2],decreasing=T),][1:10,])
    ## print(resset[order(resset[,3],decreasing=T),][1:10,])
    return(resset)
}



get_that <- function(
    bamfile,
    chr,
    regstart,
    regend
) {
    file1 <- tempfile()
    command <- paste0(
        "samtools view ",
        bamfile, " ",
        chr, ":", regstart, "-", regend, " > ",
        file1
    )
    system(command)
    if (!file.exists(file1)) {
        Sys.sleep(4)
    }
    if (!file.exists(file1)) {
        print(command)
        stop(paste0("Error making temporary file to hold samtools output", file1))
    }
    ##read in reads mapping to our chosen region
    if (file.info(file1)["size"] == 0) {
        that <- matrix(nrow = 0, ncol = 21)
    } else {
        that <- read.delim(file1, header = FALSE)
    }
    that <- as.matrix(that)
    unlink(file1)
    if (is.vector(that)) {
        that <- matrix(that,nrow=1)
    }
    return(that)
}

get_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

filter_that <- function(
    that,
    chr,
    regstart,
    regend
) {
    ## now have reads for our region
    ## now process and score
    ## remove anything mapping to wrong  chromosome (should not be any)
    if (nrow(that) > 0){
        ## above is all that is used so can filter these initially
        ## This filters out all reads that are not the desired read length, which is assumed to be the mode of all reads
        rl = get_mode(nchar(that[, 11]))
        that <- that[nchar(that[, 11]) == rl, ]

        check <- as.vector(that[,7])
        check2 <- as.vector(that[,3])
        check3 <- as.vector(that[,12])
        check4 <- as.double(as.vector(that[,8]))
        check5 <- nchar(as.vector(that[,10]))
        keep <- 1:length(check)*0+1
        keep[check2 != chr]=0
        ## where mate maps to a different chromosome, remove
        keep[check %in% paste("chr",c(1:5,7:22,"X","Y","M"),sep="")]=0
        ## where read alternately maps to a different chromosome, remove
        crosscheck=paste("chr",c(1:5,7:22,"X","Y","M"),sep="")
        for(i in 1:length(crosscheck)){
            h=grep(crosscheck[i],check3)
            if(length(h)) keep[h]=0
        }
        ##where mate maps elsewhere on chromosome 6, remove
        rule <- check == chr | (check == "=" & check2 == chr)
        if(sum(rule)) keep[rule & check4+check5+1000<=regstart]=0
        if(sum(rule)) keep[rule & check4-1000>=regend]=0
        ##table(keep)
        ##check if alternately maps elsewhere on chromosome 6
        r=rep(0,nrow(that))
        for(i in 1:nrow(that)){
            zz=unlist(strsplit(check3[i],","))
            zz=unlist(strsplit(zz,":"))
            if(length(zz)){
                h=which(zz == chr)
                r[i]=length(h)
            }
        }
        keep[r>0]=0
        that=that[keep==1,,drop=F]
    }
    return(that)
}


## robbie made function from Simon spaghetti code
get_that2 <- function(
    this,
    region,
    that,
    bamdir,
    bamfile
) {
    ## find all the hla-a cases to pull reads in
    ## only for the canonical six regions is below helpful!!
    w <- grep(paste("HLA-", region,sep=""), this[,2])
    ##temp <- as.vector(this[w, 2:3, drop = FALSE])
    temp <- this[w, 2:3, drop = FALSE]
    that2 <- that
    check <- length(w) > 0
    if (is.logical(check) && length(check) == 1 && !is.na(check) && check) {
        tryCatch(
            expr = {
                temp[, 1] <- gsub("SN:","",temp[,1])
            },
            error = function(e){
                print("---this check---")
                print(class(this))
                print(head(this))
                print("---intermediate check---")
                print(w)
                print(class(head(this[w, 2:3, drop = FALSE])))
                print(class(as.vector(head(this[w, 2:3, drop = FALSE]))))
                print("temp check")
                print(temp)
                print(class(temp))
                print(length(temp))
                print(dim(temp))
                stop("there was an error with the dimension of temp")
            }
        )
        temp[, 2] <- gsub("LN:","",temp[,2])
        temp2 <- paste(temp[, 1],collapse=" ")
        temp2 <- gsub("[*]",'\\\\*',temp2)
        ## needs samtools 1.10 or greater
        ## file2 <- file.path(outputdir, paste0(bamfile,"hla",region,"othermapping.txt"))
        file2 <- tempfile()
        cmd <- paste0(
            "samtools view ", bamfile, " ",
            temp2,
            " > ",
            file2
        )
        system(cmd, intern = TRUE, wait = TRUE)
        if (!file.exists(file2)) {
            Sys.sleep(4)
        }
        ## reads mapping to our chosen region
        ## read in, first check if empty
        aa <- system(paste0("ls -las ",file2), intern = TRUE)
        aa <- substring(aa, 1, 1)
        if(aa==0) {
            that2 <- matrix(nrow=0,ncol=ncol(that))
        }
        if(aa!=0){
            that2 <- read.delim(file2,header = FALSE)
            that2 <- as.matrix(that2)
            if(is.vector(that2)) {
                that2 <- matrix(that2, nrow = 1)
            }
        }
        unlink(file2)
    }
    that2
}


filter_that2 <- function(
    that2,
    chr,
    regstart,
    regend,
    region
) {
    ##now the hla defined read mapping filters
    ##above is all that is used so can filter these initially
    if(nrow(that2)>0){
        ## this also filters these to be their mode length
        rl = get_mode(nchar(that2[, 11]))
        that2 <- that2[nchar(that2[, 11]) == rl, ]
        
        check=as.vector(that2[,7])
        check2=as.vector(that2[,3])
        check3=as.vector(that2[,12])
        check4=as.double(as.vector(that2[,8]))
        check5=nchar(as.vector(that2[,10]))
        keep=1:length(check)*0+0
        keep[grep(paste("HLA-",region,sep=""),check2)]=1
        ##where mate maps to a different chromosome, remove
        keep[check %in% paste("chr",c(1:5,7:22,"X","Y","M"),sep="")]=0
        ##where read alternately maps to a different chromosome, remove
        crosscheck=paste("chr",c(1:5,7:22,"X","Y","M"),sep="")
        for(i in 1:length(crosscheck)){
            h=grep(crosscheck[i],check3)
            if(length(h)) keep[h]=0
        }
        ##table(keep)
        ##where mate maps elsewhere on chromosome 6, remove
        rule= check == chr | (check== "=" & check2 == chr)
        if(sum(rule)) keep[rule & check4+check5+1000<=regstart]=0
        if(sum(rule)) keep[rule & check4-1000>=regend]=0
        ##table(keep)
        ##check if alternately maps elsewhere on chromosome 6
        r=rep(0,length(check3))
        for(i in 1:length(r)){
            zz=unlist(strsplit(check3[i],","))
            zz=unlist(strsplit(zz,":"))
            if(length(zz)){
                h <- which(zz == chr)
                pos=abs(as.double(zz[h+1]))
                if(sum(pos<regstart-check5[i] | pos>regend)) r[i]=1
            }
        }
        keep[r==1]=0
        that2=that2[keep==1,,drop=F]
    }
    return(that2)
}



further_filter_that_and_that2 <- function(
    that,
    that2,
    newkmers,
    region
) {
    if(nrow(that) > 0){
        temp <- filterbyregion(
            seqs = as.vector(that[, 10]),
            seqnames = as.vector(that[, 1]),
            newkmers = newkmers,
            regname = region
        )
        if(sum(temp$keepc) > sum(temp$keep)) {
            that <- that[temp$keepc == 1, , drop = FALSE]
        }
        if(sum(temp$keepc) <= sum(temp$keep)) {
            that <- that[temp$keep == 1, , drop = FALSE]
        }
    }
    ## now more sophisticated filtering to remove reads mapping to alternative HLA alleles
    if(nrow(that2) > 0){
        if(is.vector(that2)) {
            that2 <- matrix(that2, nrow = 1)
        }
        temp2 <- filterbyregion(as.vector(that2[,10]),as.vector(that2[,1]),newkmers,region)
        if(sum(temp2$keepc)>sum(temp2$keep)) that2=that2[temp2$keepc==1,,drop=F]
        if(sum(temp2$keepc)<=sum(temp2$keep)) that2=that2[temp2$keep==1,,drop=F]
        if(is.vector(that2)) {
            that2 <- matrix(that2, nrow = 1)
        }
    }
    return(
        list(
            that = that,
            that2 = that2
        )
    )
}

## process quilt output
## process resset to avoid the forward slashes
reshape_and_filter_resset <- function(resset, region, use_averaging = TRUE) {
    newnames <- matrix(nrow=0,ncol=2)
    dd <- grep("/",rownames(resset))
    cc <- rownames(resset)
    ee <- cbind(cc,1:nrow(resset),rep(1,nrow(resset)))
    if (length(dd)){
        ee <- ee[-dd,]
	for(j in dd){
            check=cc[j]
            nn=unlist(strsplit(unlist(strsplit(check,":")),"/"))
            nn2=paste(nn[1],nn[2:length(nn)],sep=":")
            ee=rbind(ee,cbind(nn2,j,1/length(nn2)))
	}
    }
    newresset=matrix(0,nrow=length(unique(ee[,1])),ncol=3)
    rownames(newresset)=unique(ee[,1])
    for(i in 1:nrow(ee)) {
        newresset[ee[i,1],]=newresset[ee[i,1],]+as.double(ee[i,3])*resset[as.double(ee[i,2]),]
    }
    colnames(newresset)=colnames(resset)
    rownames(newresset)=paste(region,"*",rownames(newresset),sep="")
    newquiltprobs=newresset
    ## get some approximate phasing information for resset
    newresset2 <- newresset
    if (use_averaging) {
        for(i in 1:nrow(newresset2)){
            newresset2[i,1] <- newresset[i,1] / (newresset[i,1] + newresset[i,2]) * newresset[i,3] * 2
            newresset2[i,2] <- newresset[i,2] / (newresset[i,1] + newresset[i,2]) * newresset[i,3] * 2
        }
        newresset2[,1]=newresset2[,1]/sum(newresset2[,1])
        newresset2[,2]=newresset2[,2]/sum(newresset2[,2])
    }
    return(
        list(
            newresset = newresset,            
            newresset2 = newresset2,
            newquiltprobs = newquiltprobs
        )
    )
}






get_fourdigitreadscaledlikelihoodmat <- function(overall, newphasedquiltprobs) {
    ##first scale (by number of types for each four-digit code) and convert overall to likelihood scale, scaled to sum to 1
    ## only keep alleles typed to 4-digit accuracy or above
    tt=rownames(overall)
    ##count colons
    colons=1:length(tt)*0
    for(i in 1:length(tt)) colons[i]=sum(substring(tt[i],1:nchar(tt[i]),1:nchar(tt[i]))==":")
    ##keep alleles typed at four digit accuracy or above
    ##could go to six or even eight digit accuracy
    keep=tt[colons>=1]
    overall2=overall[keep,keep]
    overall2=overall2-max(overall2)
    overall2=exp(overall2)
    ##so this is raw likelihood for all 4-digit or above
    ##now extract 4-digit codes
    ##same total weight for each
    keep2=1:length(keep)
    for(i in 1:length(keep2)) {
        keep2[i]=paste(unlist(strsplit(keep[i],":"))[1:2],collapse=":")
    }
    vv=table(keep2)
    weights=1:nrow(overall2)*0
    for(i in 1:length(weights)) weights[i]=vv[keep2[i]]
    ##scale by rows
    overall2=overall2/weights
    ##scale by columns
    overall2=t(t(overall2)/weights)
    ##so weight on a particular allele is 1/number of alleles seen, corresponding to an equal likelihood of all four-digit HLA types (this is better for combining with quilt results which do weight by observed four-digit types)
    overall2=overall2/sum(overall2)
    readscaledlikelihoodmat=overall2
    ##intersection of this with quilt four digit inferences, scaled to sum to 1
    fourdigitsseen=keep2
    cond <- fourdigitsseen %in% rownames(newphasedquiltprobs)
    intersectreadscaledlikelihoodmat=readscaledlikelihoodmat[cond,cond]
    intersectreadscaledlikelihoodmat=intersectreadscaledlikelihoodmat/sum(intersectreadscaledlikelihoodmat)
    keep3=keep2[cond]
    vv2=table(keep3)
    ##
    ##
    ##
    ## four digit inferences, summing the above (should sum to 1)
    fourdigitreadscaledlikelihoodmat=matrix(0,nrow=length(vv),ncol=length(vv))
    rownames(fourdigitreadscaledlikelihoodmat)=names(vv)
    colnames(fourdigitreadscaledlikelihoodmat)=names(vv)
    ## fourdigit intersection with quilt (1000G) alleles, summing the above
    rows=match(keep2,names(vv))
    cols=rows
    for(i in 1:length(keep2)) {
        for(j in 1:length(keep2)) {
            fourdigitreadscaledlikelihoodmat[rows[i],cols[j]] <-
                fourdigitreadscaledlikelihoodmat[rows[i],cols[j]]+readscaledlikelihoodmat[i,j]
        }
    }
    return(
        list(
            fourdigitreadscaledlikelihoodmat = fourdigitreadscaledlikelihoodmat,
            vv2 = vv2,
            readscaledlikelihoodmat = readscaledlikelihoodmat,
            intersectreadscaledlikelihoodmat = intersectreadscaledlikelihoodmat
        )
    )
}


#######from below, is calling pipeline functions and then code, for the read-based calling of HLA type using only reads within each gene



getsubpos <-function(
    lookup,
    revlookup,                     
    offset = 10,
    column,
    readlength = 150
){
    ## position in bases for each allele (can be 0)
    aa <- lookup[,column]
    ## obtain position in columns for that base minus offset
    aa=aa-offset
    ## default to starting at base 1
    aa[aa<=0]=1
    aa[aa+readlength-1>ncol(revlookup)]=1
    ##end at base xx
    bb=aa+readlength-1
    bb[bb>ncol(revlookup)]=ncol(revlookup)
    temp <- revlookup[cbind(1:length(aa),bb)]
    aa[temp==0]=1
    return(
        cbind(
            revlookup[cbind(1:length(aa),aa)],
            revlookup[cbind(1:length(aa),aa+readlength-1)]
        )
    )
}


getsubmat <- function(
    lookup,
    revlookup,
    offset = 10,
    column,
    readlength = 150,
    fullalleles,
    subset = NA
) {
    if (is.na(subset[1])) {
        subset <- 1:nrow(fullalleles)
    }
    ## position in bases for each allele (can be 0)
    aa <- lookup[subset, column]
    ## obtain position in columns for that base minus offset
    aa <- aa - offset
    ## default to starting at base 1
    aa[aa<=0]=1
    aa[aa+readlength-1>ncol(revlookup)]=1
    ##end at base xx
    bb=aa+readlength-1
    ##bb[bb>ncol(revlookup)]=ncol(revlookup)
    temp <- revlookup[cbind(subset,bb)]
    aa[temp==0]=1
    retmat=matrix(nrow=length(subset),ncol=readlength)
    for(i in 1:readlength) {
        retmat[,i] <- fullalleles[cbind(subset,revlookup[cbind(subset,aa+i-1)])]
    }
    return(
        list(
            subset=subset,
            retmat=retmat,
            epos=revlookup[cbind(subset,aa+readlength-1)],
            spos=revlookup[cbind(subset,aa)]
        )
    )
}

## readposrow is a set of different positions of kmers so possible columns
## readposrow has e.g. 8 columns and offset is 1/2 length of readposrow
## readposrow identifies alignment positions based on exact kmer matches (could be adapted)

getalleles <- function(
    readposrow,
    lookup,
    revlookup,
    fullalleles,    
    offsets = c(10, 20, 120, 130),
    readlength = 150,
    thresh = 2
){

    if(sum(!is.na(readposrow))<thresh) return(list())

    offsets=offsets[!is.na(readposrow)]
    readposrow=readposrow[!is.na(readposrow)]
    newreadposrow=vector(length=0)
    newoffsets=vector(length=0)
    
    for(k in 1:length(readposrow)){
        
        temp=as.double(unlist(strsplit(readposrow[k],",")))
        newreadposrow=c(newreadposrow,temp)
        newoffsets=c(newoffsets,rep(offsets[k],length(temp)))
        
    }
    offsets=newoffsets
    readposrow=newreadposrow
    
    init <- getsubmat(
        lookup = lookup,
        revlookup = revlookup,
        fullalleles = fullalleles,
        offset = offsets[1],
        column = as.double(readposrow[1]),
        readlength = readlength
    )

    newseqs=init$retmat
    lefts=init$spos
    rights=init$epos
    allelelist=init$subset
    colshad=lefts
    if(length(offsets)>1){
	for(k in 2:length(offsets)){
            check <- getsubpos(
                lookup = lookup,
                revlookup = revlookup,
                offset = offsets[k],
                readposrow[k],
                readlength
            )
            colshad=cbind(colshad,check[,1])
            keep=rep(1,length(check[,1]))
            for(l in 1:(k-1)) keep[colshad[,k]==colshad[,l]]=0
            if(sum(keep)){
		subset=which(keep==1)
                next2 <- getsubmat(
                    lookup = lookup,
                    revlookup = revlookup,
                    fullalleles = fullalleles,
                    offset = offsets[k],
                    column = readposrow[k],
                    readlength = readlength,
                    subset = subset
                )
		lefts=c(lefts,next2$spos)
		rights=c(rights,next2$epos)
		newseqs=rbind(newseqs,next2$retmat)
		allelelist=c(allelelist,next2$subset)
            }
	}
    }
    
    return(list(lefts=lefts,rights=rights,newseqs=newseqs,allelelist=allelelist))
}

## allelelist is a vector of indexes for alignments
## newseqs is a matrix of correct length for targetseq
## targetseq is a vector of bases
## targetqual is a vector of qualities
getscores <- function(allelelist,newseqs,targetseq,targetqual,nalleles=NA,lefts=NA){
    if(is.na(nalleles)) {
        nalleles <- max(allelelist)
    }
    mm=1/3*10^(-targetqual/10)
    match=1-3*mm
    mm=log(mm)
    match=log(match)
    scoremat=matrix(match,nrow=length(targetseq),ncol=length(allelelist))
    mmmat=matrix(mm,nrow=length(targetseq),ncol=length(allelelist))
    scoremat[t(newseqs)!=targetseq]=mmmat[t(newseqs)!=targetseq]
    scores=colSums(scoremat)
    finalscores=rep(-1e9,nalleles)
    finalpos=rep(0,nalleles)
    for(i in 1:length(allelelist)) {
        if(scores[i]>finalscores[allelelist[i]]) {
            finalscores[allelelist[i]]=scores[i]
            finalpos[allelelist[i]]=lefts[i]
        }
    }
    ## best score for each allele, and best position
    return(cbind(finalscores,finalpos))
}


##use set of reads that and positions to filter reads out, returning a flag of which reads to filter
filter <- function(that,start,end,suffix="HLA-A"){
    flag=rep(0,nrow(that))
    c1=as.vector(that[,7])
    c1[c1=="="]=as.vector(that[c1=="=",3])
    c2=as.vector(that[,3])
    pos=as.double(that[,8])
    alt=as.vector(that[,12])
    typeam=rep(0,length(alt))
    typeam[grep("SA",substring(alt,1,2))]=1
    typeam[grep("XA",substring(alt,1,2))]=2
    ##remove anything whose mate-pair is wrong place on chromosome 6 (so if HLA-A or "alt", is OK)
    flag[c1=="chr6" & (pos<=start-1000+nchar(as.vector(that[,10]))-1 | pos>=end+1000)]=1
    zz=grep("HLA",c1)
    ww=grep(suffix,c1)
    if(sum(!zz %in%ww)) flag[zz[!zz%in%ww]]=1
    ## make a vector ok if both mate pairs localise specifically to this HLA allele.
    okl=flag*0
    okl[c2=="chr6" | 1:length(c2) %in% grep(suffix,c2)]=1
    okl[flag==1]=0
    okr=flag*0
    okr[c1=="chr6" | 1:length(c1) %in% grep(suffix,c1)]=1
    okr[flag==1]=0
    ok=okl*okr	
    ## print(sum(flag))
    ##remove anything with a mapping on chromosome 6 elsewhere so is "best" mapping
    for(k in 1:length(alt)){
	cc=unlist(strsplit(alt[k],";"))
        ##other alleles in region are OK
	cc=cc[-grep(suffix,cc)]
        ##alternative assemblies are OK
	if(length(cc)) cc=cc[-grep("alt",cc)]
        ##positions on chromosome 6 might not be OK and other alleles are definitely not OK
	if(length(cc)) {
            ## print(c(k,cc))
            ##if at least one possible different HLA  region interpretation, flag
            if(ok[k]==0) flag[k]=1
	}
    }
    ## print(sum(flag))
    names(flag)=that[,1]
    return(flag)
}




## function checking for alternative mappings
## now have both types of read, filter and process

## start with some sequences
## need ourfiles
## return for each read whether it should be kept, and if so which region
## newkmers is from a file indexing all HLA 10mers
## regname is target
## seqs is a string of reads
## seqnames is read IDs
filterbyregion=function(seqs,seqnames,newkmers,regname){

    readmat=matrix(nrow=length(seqs),ncol=nchar(seqs)[1])
    #*****
    rl = nchar(seqs)[1]
    #*****
    for(i in 1:nchar(seqs[1])) readmat[,i]=substring(seqs,i,i)

    compmat=readmat[,ncol(readmat):1,drop=F]
    compmat[compmat=="A"]="t"
    compmat[compmat=="C"]="g"
    compmat[compmat=="G"]="c"
    compmat[compmat=="T"]="a"
    compmat=toupper(compmat)
    compreadmat=compmat
    compseqs=vector(length=nrow(compmat))
    
    for( i in 1:length(compseqs)) compseqs[i]=paste(compmat[i,],collapse="")

######given a read it is now trivial to get info:

    readpos=cbind(newkmers[match(substring(seqs,11,20),newkmers[,1]),2],newkmers[match(substring(seqs,21,30),newkmers[,1]),2],newkmers[match(substring(seqs,rl-29,rl-20),newkmers[,1]),2],newkmers[match(substring(seqs,rl-19,rl-10),newkmers[,1]),2])
    readreg=cbind(newkmers[match(substring(seqs,11,20),newkmers[,1]),3],newkmers[match(substring(seqs,21,30),newkmers[,1]),3],newkmers[match(substring(seqs,rl-29,rl-20),newkmers[,1]),3],newkmers[match(substring(seqs,rl-19,rl-10),newkmers[,1]),3])


    readposc=cbind(newkmers[match(substring(compseqs,11,20),newkmers[,1]),2],newkmers[match(substring(compseqs,21,30),newkmers[,1]),2],newkmers[match(substring(compseqs,rl-29,rl-20),newkmers[,1]),2],newkmers[match(substring(compseqs,rl-19,rl-10),newkmers[,1]),2])
    readregc=cbind(newkmers[match(substring(compseqs,11,20),newkmers[,1]),3],newkmers[match(substring(compseqs,21,30),newkmers[,1]),3],newkmers[match(substring(compseqs,rl-29,rl-20),newkmers[,1]),3],newkmers[match(substring(compseqs,rl-19,rl-10),newkmers[,1]),3])

    ## readpos=cbind(newkmers[match(substring(seqs,11,20),newkmers[,1]),2],newkmers[match(substring(seqs,21,30),newkmers[,1]),2],newkmers[match(substring(seqs,121,130),newkmers[,1]),2],newkmers[match(substring(seqs,131,140),newkmers[,1]),2])
    ## readreg=cbind(newkmers[match(substring(seqs,11,20),newkmers[,1]),3],newkmers[match(substring(seqs,21,30),newkmers[,1]),3],newkmers[match(substring(seqs,121,130),newkmers[,1]),3],newkmers[match(substring(seqs,131,140),newkmers[,1]),3])


    ## readposc=cbind(newkmers[match(substring(compseqs,11,20),newkmers[,1]),2],newkmers[match(substring(compseqs,21,30),newkmers[,1]),2],newkmers[match(substring(compseqs,121,130),newkmers[,1]),2],newkmers[match(substring(compseqs,131,140),newkmers[,1]),2])
    ## readregc=cbind(newkmers[match(substring(compseqs,11,20),newkmers[,1]),3],newkmers[match(substring(compseqs,21,30),newkmers[,1]),3],newkmers[match(substring(compseqs,121,130),newkmers[,1]),3],newkmers[match(substring(compseqs,131,140),newkmers[,1]),3])


######turn this into a matrix of which regions and positions we might have

####all regions
    allreg=unlist(strsplit(readreg[1,1],","))
    if (nrow(readpos) >= 2) {
        for(i in 2:nrow(readpos)){
            allreg=c(allreg,unlist(strsplit(readreg[i,1],",")))
        }
    }
    
    for(i in 1:nrow(readpos)){
        allreg=c(allreg,unlist(strsplit(readreg[i,2],",")))
    }
    for(i in 1:nrow(readpos)){
        allreg=c(allreg,unlist(strsplit(readreg[i,3],",")))

    }
    for(i in 1:nrow(readpos)){
        allreg=c(allreg,unlist(strsplit(readreg[i,4],",")))

    }
    for(i in 1:nrow(readpos)){
        allreg=c(allreg,unlist(strsplit(readregc[i,1],",")))

    }
    for(i in 1:nrow(readpos)){
        allreg=c(allreg,unlist(strsplit(readregc[i,2],",")))

    }
    for(i in 1:nrow(readpos)){
        allreg=c(allreg,unlist(strsplit(readregc[i,3],",")))

    }
    for(i in 1:nrow(readpos)){
        allreg=c(allreg,unlist(strsplit(readregc[i,4],",")))

    }


    regnames=sort(unique(allreg))


########for a given read get pos
    thresh=2
    regions=matrix(0,nrow=nrow(readpos),ncol=length(regnames))
    colnames(regions)=regnames
    for(i in 1:nrow(readpos)){

####if no possibilities
        if(sum(!is.na(readpos[i,])>=thresh)){
            ww=unique(unlist(strsplit(readreg[i,1],",")))
            ww=c(ww,unique(unlist(strsplit(readreg[i,2],","))))
            ww=c(ww,unique(unlist(strsplit(readreg[i,3],","))))
            ww=c(ww,unique(unlist(strsplit(readreg[i,4],","))))

            ww=ww[!is.na(ww)]
            vv=table(ww)
            regions[i,names(vv)]=vv
        }
    }

    regionsc=matrix(0,nrow=nrow(readpos),ncol=length(regnames))
    colnames(regionsc)=regnames
    for(i in 1:nrow(readpos)){

####if no possibilities, don't do anything
        if(sum(!is.na(readpos[i,])>=thresh)){
            ww=unique(unlist(strsplit(readregc[i,1],",")))
            ww=c(ww,unique(unlist(strsplit(readregc[i,2],","))))
            ww=c(ww,unique(unlist(strsplit(readregc[i,3],","))))
            ww=c(ww,unique(unlist(strsplit(readregc[i,4],","))))

            ww=ww[!is.na(ww)]
            vv=table(ww)
            regionsc[i,names(vv)]=vv
        }
    }

    correctgaps=regions*0
###relative positions
    pos=c(11,21,rl-29,rl-19)
    # pos=c(11,21,121,131)
    k=10

    for(i in 1:nrow(regions)){
	zz=which(regions[i,]>=thresh)
	## print(i)
	if(length(zz)){

            for(k in zz){

		nn=unlist(strsplit(readreg[i,1],","))
		mm=unlist(strsplit(readpos[i,1],","))
		mm1=mm[nn==colnames(regions)[k]]
		nn=unlist(strsplit(readreg[i,2],","))
		mm=unlist(strsplit(readpos[i,2],","))
		mm2=mm[nn==colnames(regions)[k]]
		nn=unlist(strsplit(readreg[i,3],","))
		mm=unlist(strsplit(readpos[i,3],","))
		mm3=mm[nn==colnames(regions)[k]]
		nn=unlist(strsplit(readreg[i,4],","))
		mm=unlist(strsplit(readpos[i,4],","))
		mm4=mm[nn==colnames(regions)[k]]
		tt=list(mm1=mm1,mm2=mm2,mm3=mm3,mm4=mm4)
		##print(tt)
		for(p1 in 1:3) for(p2 in (p1+1):4){
                                   d1=as.double(tt[[p1]])
                                   d2=as.double(tt[[p2]])
                                   if(length(d1) & length(d2)){
                                       right=0
                                       val=1
                                       if(p2>=3 & p1<=2) val=4
                                       for(h in 1:length(d1)){
                                           if(sum(d2-d1[h]==pos[p2]-pos[p1]) & !is.na(sum(d2-d1[h]==pos[p2]-pos[p1]))) 						right=val
                                       }
                                       correctgaps[i,k]=correctgaps[i,k]+right
                                   }
                               }
            }
	}
    }

    correctgapsc=regionsc*0
###relative positions
    pos=c(11,21,rl-19,rl-29)
    # pos=c(11,21,121,131)

    for(i in 1:nrow(regionsc)){
	zz=which(regionsc[i,]>=thresh)
	##print(i)
	if(length(zz)){

            for(k in zz){

		nn=unlist(strsplit(readregc[i,1],","))
		mm=unlist(strsplit(readposc[i,1],","))
		mm1=mm[nn==colnames(regionsc)[k]]
		nn=unlist(strsplit(readregc[i,2],","))
		mm=unlist(strsplit(readposc[i,2],","))
		mm2=mm[nn==colnames(regionsc)[k]]
		nn=unlist(strsplit(readregc[i,3],","))
		mm=unlist(strsplit(readposc[i,3],","))
		mm3=mm[nn==colnames(regionsc)[k]]
		nn=unlist(strsplit(readregc[i,4],","))
		mm=unlist(strsplit(readposc[i,4],","))
		mm4=mm[nn==colnames(regionsc)[k]]
		tt=list(mm1=mm1,mm2=mm2,mm3=mm3,mm4=mm4)
		##print(tt)
		for(p1 in 1:3) for(p2 in (p1+1):4){
                                   d1=as.double(tt[[p1]])
                                   d2=as.double(tt[[p2]])
                                   if(length(d1) & length(d2)){
                                       right=0
                                       val=1
                                       if(p2>=3 & p1<=2) val=4
                                       for(h in 1:length(d1)){
                                           if(sum(d2-d1[h]==pos[p2]-pos[p1]) & !is.na(sum(d2-d1[h]==pos[p2]-pos[p1]))) 						right=val
                                       }
                                       correctgapsc[i,k]=correctgapsc[i,k]+right
                                   }
                               }
            }
	}
    }

######gaps etc.

####pick regions
    regbest=rep("None",nrow(correctgaps))
    for(i in which(rowSums(correctgaps)>0)){
	regbest[i]=paste(colnames(correctgaps)[which(correctgaps[i,]==max(correctgaps[i,]))],collapse=",")

    }

####pick regions
    regbestc=rep("None",nrow(correctgapsc))
    for(i in which(rowSums(correctgapsc)>0)){
	regbestc[i]=paste(colnames(correctgapsc)[which(correctgapsc[i,]==max(correctgapsc[i,]))],collapse=",")

    }

#####mate pair info on best mappings
    zz=unique(as.vector(seqnames))
    zzmat=matrix(nrow=length(zz),ncol=2)
    for(i in 1:length(zz)) zzmat[i,]=c(regbest[seqnames==zz[i]],"No_pair")[1:2]

####complement
    zzmatc=matrix(nrow=length(zz),ncol=2)
    for(i in 1:length(zz)) zzmatc[i,]=c(regbestc[seqnames==zz[i]],"No_pair")[1:2]

    keep=rep(0,length(regbest))
#########now back to original reads
    mat2=cbind(regbest,regbest)
    for(i in 1:nrow(mat2)) mat2[i,]=zzmat[zz==seqnames[i]]

    for(i in 1:nrow(mat2)) if(mat2[i,2] %in% c("No_pair","None") & mat2[i,1]==regname) keep[i]=1
    for(i in 1:nrow(mat2)) if(mat2[i,1]%in% c("No_pair","None") & mat2[i,2]==regname) keep[i]=1

    for(i in 1:nrow(mat2)) if(length(grep(regname,mat2[i,]))==2 & mat2[i,1]==regname) keep[i]=1

    for(i in 1:nrow(mat2)) if(length(grep(regname,mat2[i,]))==2 & mat2[i,2]==regname) keep[i]=1

#########now back to original reads
    mat2c=cbind(regbestc,regbestc)
    for(i in 1:nrow(mat2)) mat2c[i,]=zzmatc[zz==seqnames[i]]

    keepc=keep*0
    for(i in 1:nrow(mat2c)) if(mat2c[i,2] %in% c("No_pair","None") & mat2c[i,1]==regname) keepc[i]=1
    for(i in 1:nrow(mat2c)) if(mat2c[i,1]%in% c("No_pair","None") & mat2c[i,2]==regname) keepc[i]=1

    for(i in 1:nrow(mat2c)) if(length(grep(regname,mat2c[i,]))==2 & mat2c[i,1]==regname) keepc[i]=1

    for(i in 1:nrow(mat2c)) if(length(grep(regname,mat2c[i,]))==2 & mat2c[i,2]==regname) keepc[i]=1


######keep regions with this particular label
######keep regions with a best match to the appropriate region, in terms of total score, where we weight matches across length of a read higher than one end or other
######need at least 2 of 4 tested kmers to be found
######only one problematic (and odd!) paired end read for DRB1 after doing this
######largely filters out problematic reads
#######these map to DRB3 in particular and occasionally other places
#######keeps the large majority of "good" reads also
#######almost all things found this way end up being good and "kept" read pairs, so it seems a very good filter



    return(list(seqnames=seqnames,seqs=seqs,compseqs=compseqs,keep=keep,keepc=keepc,matc=mat2c,mat=mat2,regbest=regbest,regbestc=regbestc,correctgaps=correctgaps,correctgapsc=correctgapsc))
}

#####command example: temp=filterbyregion(as.vector(that[,10]),as.vector(that[,1]),newkmers,"DRB1")




##for a given region, now we have to read in the data
getbestalleles <- function(matrix,thresh=0.99){
    ##make diagonal
    for(i in 2:nrow(matrix)) {
        matrix[i,1:(i-1)]=0
    }
    diag(matrix)=diag(matrix)/2
    matrix=matrix/sum(matrix)
    bestallele1=rep(rownames(matrix),nrow(matrix))[order(matrix,decreasing=T)]
    bestallele2=rep(rownames(matrix),nrow(matrix))[order(t(matrix),decreasing=T)]
    lhoods=sort(matrix,decreasing=T)
    sums=cumsum(lhoods)
    results=cbind(bestallele1,bestallele2,lhoods,sums)
    row2=min(which(sums>=thresh))
    return(results[1:row2, , drop = FALSE])
}



do_simon_read_stuff_with_that_and_that2 <- function(
    that,
    that2,
    lookup,
    revlookup,
    fullalleles ,
    regstart,
    regend,
    region,
    kk
) {
    ## below is now true unless BOTH have no data
    if(nrow(that) > 0 & nrow(that2) > 0){
        seqs <- as.vector(that[,10])
        seqs2 <- as.vector(that2[,10])
        #**
        rl = nchar(seqs[1])
        #**
        readmat=matrix(nrow=length(seqs),ncol=nchar(seqs[1]))
        for(i in 1:nchar(seqs[1])) readmat[,i]=substring(seqs,i,i)
        readmat2=matrix(nrow=length(seqs2),ncol=nchar(seqs2[1]))
        for(i in 1:nchar(seqs2[1])) readmat2[,i]=substring(seqs2,i,i)
        ## 
        compmat=readmat[,ncol(readmat):1]
        compmat[compmat=="A"]="t"
        compmat[compmat=="C"]="g"
        compmat[compmat=="G"]="c"
        compmat[compmat=="T"]="a"
        compmat=toupper(compmat)
        compreadmat=compmat
        compseqs=vector(length=nrow(compmat))
        for( i in 1:length(compseqs)) compseqs[i]=paste(compmat[i,],collapse="")
        ##
        compmat=readmat2[,ncol(readmat):1]
        compmat[compmat=="A"]="t"
        compmat[compmat=="C"]="g"
        compmat[compmat=="G"]="c"
        compmat[compmat=="T"]="a"
        compmat=toupper(compmat)
        compseqs2=vector(length=nrow(compmat))
        for( i in 1:length(compseqs2)) compseqs2[i]=paste(compmat[i,],collapse="")
        compreadmat2=compmat
        ## 
        quals=as.vector(that[,11])
        qualmat=matrix(nrow=length(seqs),ncol=nchar(seqs[1]))
##       qualmat=matrix(nrow=length(seqs),ncol=nchar(seqs)[1])
        for(i in 1:nchar(seqs[1])) qualmat[,i]=substring(quals,i,i)
        for(i in 1:nrow(qualmat)) for(j in 1:ncol(qualmat)) qualmat[i,j]=utf8ToInt(qualmat[i,j])
        quals2=as.vector(that2[,11])
        qualmat2=matrix(nrow=length(seqs2),ncol=nchar(seqs2[1]))
##        qualmat2=matrix(nrow=length(seqs2),ncol=nchar(seqs2)[1])
        for(i in 1:nchar(seqs2[1])) qualmat2[,i]=substring(quals2,i,i)
        for(i in 1:nrow(qualmat2)) {
            for(j in 1:ncol(qualmat2)) {
                ## a handful of these might have ""
                ## I think this is that the reads are too short?
                qualmat2[i, j] <- utf8ToInt(qualmat2[i, j])
            }
        }
        qualmat=matrix(as.double(qualmat)-33,ncol=ncol(qualmat))
        qualmat2=matrix(as.double(qualmat2)-33,ncol=ncol(qualmat2))
        ## given a read it is now trivial to get info:
        readpos=cbind(kk[substring(seqs,11,20)],kk[substring(seqs,21,30)],kk[substring(seqs,rl-29,rl-20)],kk[substring(seqs,rl-19,rl-10)])
        readpos2=cbind(kk[substring(seqs2,11,20)],kk[substring(seqs2,21,30)],kk[substring(seqs2,rl-29,rl-20)],kk[substring(seqs2,rl-19,rl-10)])
        readposc=cbind(kk[substring(compseqs,11,20)],kk[substring(compseqs,21,30)],kk[substring(compseqs,rl-29,rl-20)],kk[substring(compseqs,rl-19,rl-10)])
        readpos2c=cbind(kk[substring(compseqs2,11,20)],kk[substring(compseqs2,21,30)],kk[substring(compseqs2,rl-29,rl-20)],kk[substring(compseqs2,rl-19,rl-10)])
        ## readpos=cbind(kk[substring(seqs,11,20)],kk[substring(seqs,21,30)],kk[substring(seqs,121,130)],kk[substring(seqs,131,140)])
        ## readpos2=cbind(kk[substring(seqs2,11,20)],kk[substring(seqs2,21,30)],kk[substring(seqs2,121,130)],kk[substring(seqs2,131,140)])
        ## readposc=cbind(kk[substring(compseqs,11,20)],kk[substring(compseqs,21,30)],kk[substring(compseqs,121,130)],kk[substring(compseqs,131,140)])
        ## readpos2c=cbind(kk[substring(compseqs2,11,20)],kk[substring(compseqs2,21,30)],kk[substring(compseqs2,121,130)],kk[substring(compseqs2,131,140)])
        ## print("##find alignments")
        vv <- rep(-1000,nrow(readpos))
        scoresmat=matrix(nrow=nrow(readpos),ncol=nrow(fullalleles))
        posmat=scoresmat
        compmat=scoresmat
        for(i in 1:nrow(readpos)) {
            ##print(i)
            temp <- getalleles(
                readpos[i,],
                lookup = lookup,
                revlookup = revlookup,
                fullalleles = fullalleles,
                readlength = rl
            )
            ##comp
            temp2 <- getalleles(
                readposc[i,],
                lookup = lookup,
                revlookup = revlookup,
                fullalleles = fullalleles,
                readlength = rl
            )
            ##score them
            if(length(temp)){
                ourscores=getscores(temp$allelelist,temp	$newseqs,readmat[i,],qualmat[i,],nrow(fullalleles),temp$lefts)
                scoresmat[i,]=ourscores[,1]
                posmat[i,]=ourscores[,2]
                ## print(max(ourscores[,1]))
                compmat[i,]=0
            }
            ##score them
            if(length(temp2)){
                ourscores2=getscores(temp2$allelelist,temp2	$newseqs,compreadmat[i,],qualmat[i,ncol(qualmat):1],nrow(fullalleles),temp2$lefts)
                if(!length(temp)){
                    scoresmat[i,]=ourscores2[,1]
                    posmat[i,]=ourscores2[,2]
                    compmat[i,]=1
                    ##print(max(scoresmat[i,]))
                } else {
                    if(max(ourscores2[,1])>max(ourscores[,1])){
                        scoresmat[i,]=ourscores2[,1]
                        posmat[i,]=ourscores2[,2]
                        compmat[i,]=1
                        ##print(max(scoresmat[i,]))
                    }
                }
            }
        }
        ## print("#####find alignments")
        scoresmat2=matrix(nrow=nrow(readpos2),ncol=nrow(fullalleles))
        posmat2=scoresmat2
        compmat2=scoresmat2
        for(i in 1:nrow(readpos2)) {
            ## print(i)
            temp <- getalleles(
                readpos2[i,],
                lookup = lookup,
                revlookup = revlookup,
                fullalleles = fullalleles,
                readlength = rl
            )
            temp2 <- getalleles(
                readpos2c[i,],
                lookup = lookup,
                revlookup = revlookup,
                fullalleles = fullalleles ,
                readlength = rl               
            )
            ##score them
            if(length(temp)){
                ourscores=getscores(temp$allelelist,temp	$newseqs,readmat2[i,],qualmat2[i,],nrow(fullalleles),temp$lefts)
                scoresmat2[i,]=ourscores[,1]
                posmat2[i,]=ourscores[,2]
                compmat2[i,]=0
                ##print(max(ourscores[,1]))
            }
            ##score them
            if(length(temp2)){
                ourscores2=getscores(temp2$allelelist,temp2	$newseqs,compreadmat2[i,],qualmat2[i,ncol(qualmat2):1],nrow(fullalleles),temp2$lefts)
                if(!length(temp)){
                    scoresmat2[i,]=ourscores2[,1]
                    posmat2[i,]=ourscores2[,2]
                    compmat2[i,]=1
                    ##print(max(scoresmat2[i,]))
                }
                if(length(temp)) if(max(ourscores2[,1])>max(ourscores[,1])){
                                     scoresmat2[i,]=ourscores2[,1]
                                     posmat2[i,]=ourscores2[,2]
                                     compmat2[i,]=1
                                     ##print(max(scoresmat2[i,]))
                                 }
            }
        }
        ## print("#####scoring and mismatches, first have a look")
        maxmismatch=5
        minscore=maxmismatch*log(0.001/3)
        maxes2=1:nrow(scoresmat2)*0
        maxpos2=maxes2
        for(i in 1:length(maxes2)){
            maxes2[i]=max(scoresmat2[i,])
            maxpos2[i]=mean(posmat2[i,which(scoresmat2[i,]==maxes2[i])])
        }
        maxes=1:nrow(scoresmat)*0
        maxpos=maxes
        for(i in 1:length(maxes)){
            maxes[i]=max(scoresmat[i,])
            maxpos[i]=mean(posmat[i,which(scoresmat[i,]==maxes[i])])
        }
        ##plot(maxpos[maxes>=minscore],maxes[maxes>=minscore])
        ##plot(maxpos2[maxes2>=minscore],maxes2[maxes2>=minscore])
        maxes3=c(maxes,maxes2)
        maxpos3=c(maxpos,maxpos2)
        ##plot(maxpos3[maxes3>=minscore],maxes3[maxes3>=minscore])
        ##print("#######scoring - we need to name reads etc.")
        rownames(scoresmat)=that[,1]
        rownames(scoresmat2)=that2[,1]
        rownames(posmat)=that[,1]
        rownames(posmat2)=that2[,1]
        rownames(compmat)=that[,1]
        rownames(compmat2)=that2[,1]
        ##filters to get rid of reads with mapping problems
        zz=filter(that,start=regstart,end=regend,suffix=paste("HLA-",region,sep=""))
        zz2=filter(that2,start=regstart,end=regend,suffix=paste("HLA-",region,sep=""))
        ##print("########indexing which of reads is it, read 1 or read2")
        readind=as.integer(as.numeric(that[,2])/64)%%4
        readind2=as.integer(as.numeric(that2[,2])/64)%%4
        cond=!is.na(maxes) & maxes>=minscore & zz==0
        cond2=!is.na(maxes2) & maxes2>=minscore & zz2==0
        allreadind=c(readind[cond],readind2[cond2])
        allscores=rbind(scoresmat[cond,],scoresmat2[cond2,])
        allpos=rbind(posmat[cond,],posmat2[cond2,])
        ##print("########paired reads, combine")
        pairedscores=matrix(ncol=ncol(allscores),nrow=length(unique(rownames(allscores))))
        rownames(pairedscores)=unique(rownames(allscores))
        if (nrow(pairedscores) > 0) {
            for(i in 1:nrow(pairedscores)){
                t1=allscores[rownames(allscores)==rownames(pairedscores)[i],]
                t2=allreadind[rownames(allscores)==rownames(pairedscores)[i]]
                pairedscores[i,]=0
                if(length(t2)==1) pairedscores[i,]=t1
                if(length(t2)>1){
                    if(sum(t2==1)){
                        pairedscores[i,]=pairedscores[i,]+t1[which(t2==1)[1],]
                    }
                    if(sum(t2==2)){
                        pairedscores[i,]=pairedscores[i,]+t1[which(t2==2)[1],]
                    }
                }
            }
        }
        ##print_message("Get score for pairs of alleles, and again filter so now we only allow e.g. equivalent of 5 mismatches among 300bp")
        ## print_message("First remove read pairs mismatching all alleles strongly")
        if (nrow(pairedscores) > 0) {
            newmaxes=1:nrow(pairedscores)*0
            for(i in 1:length(newmaxes)) {
                newmaxes[i]=max(pairedscores[i,])
            }
            ##remove
            pairedscores=pairedscores[newmaxes>=minscore, , drop = FALSE]
            newmaxes=newmaxes[newmaxes>=minscore]
            ##get scores relative to best allele
            temp=pairedscores-newmaxes
            ##for very bad scores, add a small amount
            ##or could come from elsewhere in the genome so it is well worth being conservative here
            temp=0.5*exp(temp)+1e-100
        }
        overall=matrix(0,nrow=ncol(pairedscores),ncol=ncol(pairedscores))
        colnames(overall)=rownames(fullalleles)
        rownames(overall)=rownames(fullalleles)
        qq=overall*0
        ##print_message("Percent final scoring complete:")
        if (nrow(pairedscores) > 0) {
            for(i in 1:nrow(pairedscores)){
                ##print(i/nrow(pairedscores)*100)
                qq=qq*0
                m1=qq+temp[i,]
                m1=m1+t(m1)
                overall=overall+log(m1)
            }
        }
        readlikelihoodmat <- overall
        readset1 <- that
        readset2 <- that2
        readscaledlikelihoodmat <- NULL ## overwritten later it looks like
        fourdigitreadscaledlikelihoodmat <- NULL ## overwritten later it looks like
        intersectfourdigitreadscaledlikelihoodmat <- NULL
        intersectquiltscaledlikelihoodmat <- NULL
        intersectreadscaledlikelihoodmat <- NULL
        combinedscaledlikelihoodmat <- NULL
        combinedresults <- NULL
        mappingonlyresults <- NULL
    } else {
        readlikelihoodmat <- NULL
        readset1 <- NULL
        readset2 <- NULL
        pairedscores <- NULL
        readscaledlikelihoodmat <- NULL
        fourdigitreadscaledlikelihoodmat <- NULL
        intersectfourdigitreadscaledlikelihoodmat <- NULL
        intersectquiltscaledlikelihoodmat <- NULL
        intersectreadscaledlikelihoodmat <- NULL
        combinedscaledlikelihoodmat <- NULL
        combinedresults <- NULL
        mappingonlyresults <- NULL
        overall <- NULL
    }
    return(
        list(
            readlikelihoodmat = readlikelihoodmat,
            readset1 = readset1,
            readset2 = readset2,
            pairedscores = pairedscores,
            readscaledlikelihoodmat = readscaledlikelihoodmat,
            fourdigitreadscaledlikelihoodmat = fourdigitreadscaledlikelihoodmat,
            intersectfourdigitreadscaledlikelihoodmat = intersectfourdigitreadscaledlikelihoodmat,
            intersectquiltscaledlikelihoodmat = intersectquiltscaledlikelihoodmat,
            intersectreadscaledlikelihoodmat = intersectreadscaledlikelihoodmat,
            combinedscaledlikelihoodmat = combinedscaledlikelihoodmat,
            combinedresults = combinedresults,
            mappingonlyresults = mappingonlyresults,
            overall = overall
        )
    )
}







summarize_all_results_and_write_to_disk <- function(
    all_results,
    sampleNames,
    what,
    outputdir,
    summary_prefix,
    summary_suffix,
    only_take_top_result = FALSE
) {
    result <- data.frame(rbindlist(lapply(1:length(all_results), function(iSample) {
        x <- all_results[[iSample]]
        y <- data.frame(
            sample_number = iSample,
            sample_name = sampleNames[iSample],
            x[[what]]
        )
        ## re-name some of the output
        colnames(y)[colnames(y) == "lhoods"] <- "post_prob"
        if (only_take_top_result) {
            y <- y[1, , drop = FALSE]
        }
        y
    })))
    file <- paste0(summary_prefix, ".", summary_suffix)    
    if (outputdir != "") {
        file <- file.path(outputdir, file)
    }
    write.table(
        result,
        file = file,
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )
    result
}
    
