phase_hla_haplotypes <- function(
    outputdir,
    chr,
    hla_gene_information,
    full_regionStart,
    full_regionEnd,
    buffer,
    regions,    
    region_exclude_file,
    reference_haplotype_file,
    reference_legend_file,
    reference_sample_file,
    reference_exclude_samplelist_file,
    reference_exclude_samples_for_initial_phasing,
    hla_types_panel,
    genetic_map_file,    
    minRate,
    nGen,
    nCores
) {


    if (1 == 0 )  {

        ## not needed per-se
        inputs_dir <- "/data/smew1/rdavies/quilt_hla_finalize_with_simon/inputs_2021_03_25/"
        outputdir <- "/data/smew1/rdavies/quilt_hla_finalize_with_simon/HLA_TEST_2021_03_25"
        chr <- "chr6"
        hla_gene_region_file <- "~/proj/QUILT/hla_ancillary_files/hlagenes.txt"
        full_regionStart <- 25587319
        full_regionEnd <- 33629686
        buffer <- 500000
        regions <- c("A", "B", "C", "DQB1", "DRB1")
        region_exclude_file <- "~/proj/QUILT/hla_ancillary_files/hlagenes.txt"
        reference_haplotype_file <- file.path(inputs_dir, "hrc.chr6.hap.clean.small.gz")
        reference_legend_file <- file.path(inputs_dir, "hrc.chr6.legend.clean.small.gz")
        reference_sample_file <- file.path(inputs_dir, "hrc.chr6.samples.reheadered2")
        reference_exclude_samplelist_file <- ""
        reference_exclude_samples_for_initial_phasing <- FALSE
        hla_types_panel <- file.path(outputdir, "20181129_HLA_types_full_1000_Genomes_Project_panel.txt")
        genetic_map_file <- file.path(inputs_dir, "CEU-chr6-final.b38.txt.gz")
        minRate <- 0.1
        nGen <- 100
        nCores <- 6
        ##reference_exclude_samplelist_file <- file.path(outputdir, "hlauntyped.exclude.txt")
        ##removeinds_file <- "" ## need to do this!
        ##cat(c("NA12878", "NA18566"), file = file.path(outputdir, "exclude.test.txt"), sep = "\n")
        ##removeinds_file <-  file.path(outputdir, "exclude.test.txt")
        ## new Robbie prepared input file list
        ##full_reference_hap_file <- "hrc.chr6.25587319.33629686.pretendhlatyped.RData" ## re-name later
        
    }


    ##
    ## determine people we don't want at all to extract
    ##
    print_message("Determine who to remove from full panel")
    hlatypes <- read.table(hla_types_panel, header = TRUE, skip = 0, sep = "\t")
    ref_samples <- read.table(reference_sample_file, header = TRUE)
    ## keep only well defined samples
    if (!("SAMPLE" %in% colnames(ref_samples))) {
        msg <- paste0("Expecting columns named 'SAMPLE POP GROUP SEX' in reference_sample_file=", reference_sample_file)
        stop(msg)
    }
    remove <- !(ref_samples[, "SAMPLE"] %in% hlatypes[, "Sample.ID"])
    samples_to_remove <- ref_samples[remove, "SAMPLE"]
    if (reference_exclude_samples_for_initial_phasing) {
        more_samples_to_remove <- as.character(read.table(reference_exclude_samplelist_file, stringsAsFactors = FALSE)[, 1])
        samples_to_remove <- c(samples_to_remove, more_samples_to_remove)
    }
    if (length(samples_to_remove) == 0) {
        full_reference_exclude_samplelist_file <- ""
    } else {
        full_reference_exclude_samplelist_file <- file.path(outputdir, "hlauntyped.exclude.txt")
        write.table(
            matrix(samples_to_remove, ncol = 1),
            file = full_reference_exclude_samplelist_file,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
        )
    }
    print_message("Done determining who to remove from full panel")
    


    
    ##
    ## do first extraction here (maybe?)
    ##
    QUILT_prepare_reference(
        outputdir = outputdir,
        nGen = nGen,
        chr = chr,
        regionStart = full_regionStart,
        regionEnd = full_regionEnd,
        buffer = buffer,
        reference_haplotype_file = reference_haplotype_file,
        reference_legend_file = reference_legend_file,
        reference_sample_file = reference_sample_file,
        genetic_map_file = genetic_map_file,
        reference_exclude_samplelist_file = full_reference_exclude_samplelist_file,
        output_file = file_quilt_hla_all_haplotypes(outputdir),
        minRate = minRate
    )




    
    ##
    ## this performs the initial phasing with a large set, using all SNPs
    ## later on, we remove individuals, and keep only the good set, and good SNPs
    ##
    print_message("Begin assigning HLA types to haplotypes")
    hla_perform_step_1_phasing(
        outputdir = outputdir,
        regions = regions,
        hla_types_panel = hla_types_panel
    )    

    ##
    ## determine who we want to keep
    ##
    print_message("Determine who to remove post HLA phasing")    
    hla_phasing_determine_who_to_remove(
        outputdir = outputdir,
        hla_types_panel = hla_types_panel,
        reference_sample_file = reference_sample_file,
        reference_exclude_samplelist_file = reference_exclude_samplelist_file,
        regions = regions
     )

    
    ##
    ## do second extraction here (maybe?)
    ## TODO - make this more efficient based on above already done data
    ##
    print_message("Extract new subsetted panels")
    ## what actually changes?
    ## exclude individuals, SNPs?
    ## for now, just operationalize, then fix later?
    out <- parallel::mclapply(regions, mc.cores = nCores, function(region) {
        i <- match(paste0("HLA-", region), hla_gene_information[, "Name"])
        regionStart <- hla_gene_information[i, "Start"]
        regionEnd <- hla_gene_information[i, "End"]
        ## if (region == "A")    {regionStart=29942554; regionEnd=29945741}
        ## if (region == "B")    {regionStart=31353367; regionEnd=31357155}
        ## if (region == "C")    {regionStart=31268257; regionEnd=31353367}
        ## if (region == "DRB1") {regionStart=32578780; regionEnd=32589729}
        ## if (region == "DQB1") {regionStart=32660035; regionEnd=32666603}
        ## if (region == "DQA1") {regionStart=32637480; regionEnd=32643199}
        QUILT_prepare_reference(
            outputdir = outputdir,
            nGen = nGen,
            reference_haplotype_file = reference_haplotype_file,
            reference_legend_file = reference_legend_file,
            reference_sample_file = reference_sample_file,
            chr = chr,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            genetic_map_file = genetic_map_file,
            reference_exclude_samplelist_file = file_quilt_hla_after_step_1_who_to_remove(outputdir, region),
            output_file = file_quilt_hla_specific_haplotypes(outputdir, region),
            region_exclude_file = region_exclude_file,
            minRate = minRate
        )
        return("Done")
    })
    check_mclapply_OK(out)

    ##
    ## get HLA types for each individual in a panel
    ##
    print_message("Get HLA types for each individual in a panel")
    hla_perform_step_2_phasing(
        outputdir = outputdir,
        regions = regions,
        hla_types_panel = hla_types_panel
    )
    
    print_message("Done assigning HLA types to haplotypes")

    NULL
}




if (1 == 0) {

    ## rwd: simon previous preamble comments
    
    ##
    ## input file description
    ##

    ## downloaded file of HLA types from publication on 1000G:
    ## 20181129_HLA_types_full_1000_Genomes_Project_panel.txt
    ## From: Immune diversity sheds light on missing variation in worldwide genetic diversity panels.
    ## Abi-Rached L, Gouret P, Yeh J-H, Di Cristofaro J, Pontarotti P, Picard C, and Paganini J.
    ## Published: October 26, 2018 DOI: 10.1371/journal.pone.0206512
    ## https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0206512

    ## reference panel in HLA regions, file hrc.chr6.25587319.33629686.pretendhlatyped.RData
    ## made with command (old code)
    ##   cd /data/smew1/rdavies/single_imp/2020_06_25/
    ##   R -f /datas/muscovy/not-backed-##up/shi/covid/prepare_region_hrc.R --args \
    ##   ref_panels/hrc.chr6.hap.clean.gz \
    ##   ref_panels/hrc.chr6.legend.clean.gz \
    ##   chr6 \
    ##   25587319 \
    ##   33629686 \
    ##   500000 \
    ##   /datas/muscovy/not-backed-up/shi/covid/hrc.chr6.25587319.33629686.RData \
    ##   /datas/muscovy/not-backed-up/shi/covid/hrc.chr6.25587319.33629686.sites.vcf.gz \
    ##   ref_panels/hrc.chr6.samples \
    ##   /data/smew1/rdavies/external/recomb/CEU/CEU-chr6-final.b38.txt.gz \
    ##   /datas/muscovy/not-backed-up/shi/covid/hlauntyped.exclude.txt \
    ##   /datas/muscovy/not-backed-up/shi/covid/hlapretendgenes.txt \
    ##   /datas/muscovy/not-backed-up/shi/covid/hrc.chr6.25587319.33629686.pretendhlatyped.RData


    ## other input files
    ## hla*snpformatalleles.out
    ## should already be made
    ## hla*fullallelesfilledin.out
    ## should already be made

    ## input file
    ## hrc.chr6.hla*.hlatyped.RData

    
    ##
    ## output file description
    ##

    ## temporary output file
    ## hla*newphased.out
    ## output file 
    ## hla*haptypes.out

}





hla_perform_step_1_phasing <- function(
    outputdir,
    regions,
    hla_types_panel
) {

    ## reorder the HLA types to match the reference panel information
    ## at this point we have two tables, "ref_samples" which gives a set of individuals in (in our case) HRC with reference panel data on SNPs, and "hlatypes2" which gives a set of unphased HLA types for these individuals in each region
    ## now, we need to phase the HLA types to combine with reference panel information, and to do this we need SNP format HLA data



    ## then can intersect prior to phasing
    ## hlatypes2 <- hlatypes[match(ref_samples[,1],hlatypes[,3]),]
    load(file_quilt_hla_all_haplotypes(outputdir))
    hlatypes <- read.table(hla_types_panel,header=T,skip=0,sep="\t")
    hlatypes2 <- hlatypes[match(reference_samples[,1],hlatypes[,3]),]
    
    for(region in regions){

        print_message(paste0("Work on phasing in region:", region))
        ## alternative
        ## load(paste("hla",region,"snpformatalleles.out",sep=""))
        load(file_quilt_hla_snpformatalleles(outputdir, region))
        load(file_quilt_hla_full_alleles_filled_in(outputdir, region))
        ## load(paste("hla",region,"fullallelesfilledin.out",sep=""))

        ##keep only biallelic SNPs
        keep=knownvarsfiltered[,2] %in% c("A","C","G","T") & knownvarsfiltered[,3] %in% c("A","C","G","T")

        ids=paste("chr6:",knownvarsfiltered[keep,1],sep="")
        oldsnpinfo=cbind(ids,knownvarsfiltered[keep,])
        ##haps=t(resmat[,keep])
        samples=rownames(resmat)
        colnames(oldsnpinfo)=c("id","position","a0","a1")

        ss=fullalleles[,match(oldsnpinfo[,2],ourpos)]
        oldsnpinfo=oldsnpinfo[colSums(is.na(ss))==0,]
        ss=ss[,colSums(is.na(ss))==0]

        ss=t(ss)
        ss[ss==oldsnpinfo[,3]]=0
        ss[ss==oldsnpinfo[,4]]=1
        qq=as.double(oldsnpinfo[,2])
        zz=1:nrow(ss)*0;for(i in 1:length(zz)) zz[i]=sum(ourpos==qq[i])

        zz2=zz*0
        zz2=rowSums(!(ss==0 | ss==1))
        zz3=zz*0
        for(i in 1:length(zz3)) zz3[i]=sum(qq==qq[i])
        zz4=rowSums(!(ss==0 | ss==1 | ss=="."))
        cond=(zz3==1 & zz==1 & zz2<ncol(ss)*0.1 & zz2<0.5*rowSums(ss==1) & zz4==0)


        ##keep sites uniquely mapping, not overlapping another SNP, with at most 10% gaps and 2-fold more non-ancestral than gaps
        ##print(c(sum(cond),sum(!cond)))

        snpinfo=oldsnpinfo[cond,]
        ss=ss[cond,]
        ss[ss=="."]=0
        haps=ss

        ##load(paste("hla",region,"snpformatalleles.out",sep=""))
        ##keep only biallelic SNPs
        ##keep=knownvarsfiltered[,2] %in% c("A","C","G","T") & ##knownvarsfiltered[,3] %in% c("A","C","G","T")
        ##ids=paste("chr6:",knownvarsfiltered[keep,1],sep="")
        ##snpinfo=cbind(ids,knownvarsfiltered[keep,])
        ##haps=t(resmat[,keep])
        ##samples=rownames(resmat)
        ##colnames(snpinfo)=c("id","position","a0","a1")

        ##check levels of missingness dealt with

        ##will check alleles match also (might want to check unique positions too!)

        ##knownvarsfiltered[as.double(knownvarsfiltered[,1])>=32579017 & ##knownvarsfiltered[,2] %in% c("A","C","G","T") & ##knownvarsfiltered[,3] %in% c("A","C","G","T"),]
        ##zz=1:nrow(ss)*0;for(i in 1:length(zz)) ##zz[i]=sum(ourpos==oldsnpinfo[i,2])

        ##snpinfo=read.table("HLA A .snpinfo",as.is=T)
        ##haps=read.table(paste("HLA ",region," .haps",sep=""), as.is=T)
        ##snpinfo=read.table(paste("HLA ",region," .snpinfo",sep=""),as.is=T,header=T)
        ##samples=read.table(paste("HLA ",region," .samples",sep=""),as.is=T,header=T)

        ##haps=haps[,seq(1,ncol(haps),2)]
        ##new
        samples=matrix(samples,ncol=1)
        c1=strsplit(samples[,1],"[*]")
        vv=matrix(nrow=nrow(samples),ncol=3)
        ww=1:length(c1)
        for(i in 1:length(ww)) ww[i]=unlist(c1[[i]][2])
        for(i in 1:length(c1)) vv[i,1]=unlist(c1[[i]][1])
        c1=strsplit(ww,"[:]")
        for(i in 1:length(c1)) vv[i,2]=unlist(c1[[i]][1])
        for(i in 1:length(c1)) vv[i,3]=unlist(c1[[i]][2])

        fourdigit=paste(vv[,2],":",vv[,3],sep="")

        ufourdigit=unique(fourdigit)

        haps=matrix(as.double(haps),ncol=ncol(haps))

        newhaps=matrix(ncol=length(ufourdigit),nrow=nrow(haps))
        colnames(newhaps)=ufourdigit

        for(i in 1:ncol(newhaps)) {
            ## print(i)
            tt=which(fourdigit==ufourdigit[i])
            if(length(tt)<2) tt=c(tt,tt)
            newhaps[,i]=rowMeans(haps[,tt])
        }

        v1=paste(pos[,2],pos[,3],pos[,4])
        v2=paste(snpinfo[,2],snpinfo[,3],snpinfo[,4])
        v2f=paste(snpinfo[,2],snpinfo[,4],snpinfo[,3])


        ##need to be on same strand, same ref/alt
        same=match(v1,v2)
        flipped=match(v1,v2f)

        locs=which(!is.na(same) | !is.na(flipped))
        ##print(length(locs))
        if(sum(!is.na(flipped))) same[!is.na(flipped)]=flipped[!is.na(flipped)]

        same=same[locs]
        flipped=flipped[locs]

        newhaps2=newhaps[same,]

        if(sum(!is.na(flipped))) newhaps2[(!is.na(flipped)),]=1-newhaps2[(!is.na(flipped)),]


        ##only a subset of SNPs
        translate=function(matchlist){
            qq=matchlist
            signs=rep(0,length(matchlist))
            signs[qq<0]=1
            qq[signs==1]=qq[signs==1]+2^31
            res=matrix(nrow=length(matchlist),ncol=0)
            for(i in 1:31){
                tempres=qq%%2
                qq=(qq-tempres)/2
                res=cbind(res,tempres)
            }
            res=cbind(res,signs)
            res[res==1]=0.999
            res[res==0]=0.001
            return(res)
        }
        ##now find types

        cols=grep(paste("HLA.",region,".",sep=""),colnames(hlatypes2))

        ourtypes1=as.vector(hlatypes2[,cols[1]])
        ourtypes1=unlist(gsub("[*]","",ourtypes1))

        ##if multiple possible types, assign to lowest number
        vv=strsplit(ourtypes1,"[/]")
        for(i in 1:length(ourtypes1)) ourtypes1[i]=vv[[i]][1]

        ourtypes2=as.vector(hlatypes2[,cols[2]])
        ourtypes2=unlist(gsub("[*]","",ourtypes2))

        vv=strsplit(ourtypes2,"[/]")
        for(i in 1:length(ourtypes2)) ourtypes2[i]=vv[[i]][1]


        ##now SNP info
        startandend=range(locs)

        start32=floor(startandend[1]/32)+1
        end32=floor(startandend[2]/32)+1

        start=(start32-1)*32+1
        end=end32*32

        tempIE=distinctHapsIE[,start:end]
        tempB=distinctHapsB[,start32:end32]

        temp=rhb_t[,start32:end32]
        temp2=matrix(nrow=nrow(temp),ncol=ncol(tempIE))
        ##NA case interpretation
        cc=c(rep(0.001,31),0.999)
        for(i in 1:ncol(temp)){
            temp2[,(1+(i-1)*32):(i*32)]=tempIE[match(temp[,i],tempB[,i]),(1+(i-1)*32):(i*32)]
            cc2=which(is.na(temp[,i]))
            if(length(cc2)) {
                cc3=matrix(cc,nrow=length(cc2),ncol=length(cc),byrow=T)
                temp2[cc2,(1+(i-1)*32):(i*32)]=cc3
            }
            ##final ones
            needed=which(is.na(temp2[,i*32]))
            if(length(needed)) temp2[needed,(1+(i-1)*32):(i*32)]=translate(temp[needed,i])
        }

        newlocs=locs-start+1
        hrchapstomatch=temp2[,newlocs]
        hrcfirstalleles=hrchapstomatch[seq(1,nrow(hrchapstomatch),2),]
        hrcsecondalleles=hrchapstomatch[seq(2,nrow(hrchapstomatch),2),]
        temp=1:length(ourtypes1)*0-1
        temp[ourtypes1%in%colnames(newhaps2)]=1
        predfirstalleles=matrix(ncol=nrow(newhaps2),nrow=length(temp))
        predfirstalleles[temp==1,]=t(newhaps2[,ourtypes1[temp==1]])
        temp=1:length(ourtypes2)*0-1
        temp[ourtypes2%in%colnames(newhaps2)]=1
        predsecondalleles=matrix(ncol=nrow(newhaps2),nrow=length(temp))
        predsecondalleles[temp==1,]=t(newhaps2[,ourtypes2[temp==1]])
        dist11=rowSums(abs((hrcfirstalleles)-predfirstalleles))
        dist12=rowSums(abs((hrcfirstalleles)-predsecondalleles))
        dist21=rowSums(abs((hrcsecondalleles)-predfirstalleles))
        dist22=rowSums(abs((hrcsecondalleles)-predsecondalleles))
        w=1:length(dist11)*0
        w[dist11>dist12]=1
        w[dist11<dist12]= -1
        ##table(w)
        w1=w
        w=1:length(dist21)*0
        w[dist22<dist21]= 1
        w[dist22>dist21]= -1
        w2=w
        ##table(w1,w2)


        obsgen=hrcfirstalleles+hrcsecondalleles
        predgen=predfirstalleles+predsecondalleles
        ## plot(colMeans(obsgen),colMeans(predgen,na.rm=T))
        corr=1:ncol(obsgen)
        for(i in 1:length(corr)) {
            corr[i] <- suppressWarnings(cor(obsgen[,i],predgen[,i],use="pair"))
        }

        dist11=rowSums((abs((hrcfirstalleles)-predfirstalleles))[,corr>0.8 & !is.na(corr)])
        dist12=rowSums((abs((hrcfirstalleles)-predsecondalleles))[,corr>0.8 & !is.na(corr)])
        dist21=rowSums((abs((hrcsecondalleles)-predfirstalleles))[,corr>0.8 & !is.na(corr)])
        dist22=rowSums((abs((hrcsecondalleles)-predsecondalleles))[,corr>0.8 & !is.na(corr)])

        w=1:length(dist11)*0
        w[dist11<dist12]=1
        w[dist11>dist12]= -1
        table(w)
        w1=w
        w=1:length(dist21)*0
        w[dist22<dist21]= 1
        w[dist22>dist21]= -1
        w2=w
        table(w1,w2)

        ##cbind(ourtypes1,ourtypes2,dist11,dist12,dist21,dist22)[which(w2==1 & w1==-1),]

        phase1=dist11+dist22
        phase2=dist12+dist21

        ## cbind(ourtypes1,ourtypes2,dist11,dist12,dist21,dist22)[!is.na(dist11) & !is.na(phase1) & phase1>4 & phase2>4,]
        d11=dist11
        d21=dist21
        d12=dist12
        d22=dist22

        phased= (!is.na(phase1) & !is.na(phase2) & phase1<4 & phase2>4) |
            (!is.na(phase1) & !is.na(phase2) & phase1>4 & phase2<4) | (!is.na(phase1) & !is.na(phase2) & phase1-phase2>2 ) | (!is.na(phase1) & !is.na(phase2) & phase2-phase1>2) | (!is.na(ourtypes1) & !is.na(ourtypes2) & ourtypes1==ourtypes2) | (is.na(d21) & !is.na(d12) & d22-d12>2) | (is.na(d21) & !is.na(d12) & d12-d22>2) | (is.na(d12) & !is.na(d21) & d11-d21>2) | (is.na(d12) & !is.na(d21) & d21-d11>2)

        ##print("Remaining to be phased")
        ##print(table(phased))

        ##106 left to phase

        ##temp=cbind(ourtypes1,ourtypes2,dist11,dist12,dist21,dist22,phase1,phase2)[!phased,][1:20,]

        ##2388 phased alleles
        ##to phase rest, try using known phases to remake predictions,  using all sites

        ##so in genic region, overlapping HRC and database sites, phase if there are <4 mismatches with one phasing and not the other, or one phasing improves the other by at least 2 mismatches, or if an allele is not predicted or not in the database, the other allele predicts phasing with a difference of  at least 2 mismatches
        ##some phasing is trivial also
        ##other cases to be phased later

        ##non-flipped
        phased1=(!is.na(phase1) & !is.na(phase2) & phase1<4 & phase2>4) |
            (!is.na(phase1) & !is.na(phase2) & phase2-phase1>2 & phase1<4) | (!is.na(ourtypes1) & !is.na(ourtypes2) & ourtypes1==ourtypes2) | 
            (is.na(d21) & !is.na(d12) & d12-d22>2 & d22<2) |  (is.na(d12) & !is.na(d21) & d21-d11>2 & d11<2)

        ##flipped
        phased2=(!is.na(phase1) & !is.na(phase2) & phase1>4 & phase2<4) |
            (!is.na(phase1) & !is.na(phase2) & phase1-phase2>2 & phase2<4) | 
            (is.na(d21) & !is.na(d12) & d22-d12>2 & d12<2) |  (is.na(d12) & !is.na(d21) & d11-d21>2 & d21<2)

        ##first and second alleles

        allele1=vector(length=length(ourtypes1))
        allele2=allele1
        allele1[phased1==T]=ourtypes1[phased1==T]
        allele2[phased1==T]=ourtypes2[phased1==T]
        allele1[phased2==T]=ourtypes2[phased2==T]
        allele2[phased2==T]=ourtypes1[phased2==T]



        ##now use current phasing to make haplotypes for each possible allele

        ##to start with
        phased1old=phased1
        phased2old=phased2
        oldphase1=phase1
        oldphase2=phase2
        phased1=phased1old
        phased2=phased2old
        for(extension in seq(50,1000,50)){
            plot(oldphase1,oldphase2,col=1+as.double(phased1)+2*as.double(phased2))
            startandend=range(locs)
            startandend[1]=startandend[1]-extension
            startandend[2]=startandend[2]+extension
            start32=floor(startandend[1]/32)+1
            end32=floor(startandend[2]/32)+1
            start=(start32-1)*32+1
            end=end32*32
            tempIE=distinctHapsIE[,start:end]
            tempB=distinctHapsB[,start32:end32]
            temp=rhb_t[,start32:end32]
            temp2=matrix(nrow=nrow(temp),ncol=ncol(tempIE))
            cc=c(rep(0.001,31),0.999)
            for(i in 1:ncol(temp)){
                temp2[,(1+(i-1)*32):(i*32)]=tempIE[match(temp[,i],tempB[,i]),(1+(i-1)*32):(i*32)]
                cc2=which(is.na(temp[,i]))
                if(length(cc2)) {
                    cc3=matrix(cc,nrow=length(cc2),ncol=length(cc),byrow=T)
                    temp2[cc2,(1+(i-1)*32):(i*32)]=cc3
                }
                needed=which(is.na(temp2[,i*32]))
                if(length(needed)) temp2[needed,(1+(i-1)*32):(i*32)]=translate(temp[needed,i])
            }
            temp2=temp2[,(startandend[1]-start+1):(startandend[2]-start+1)]
            alleles=vector(length=nrow(temp2))
            alleles[seq(1,length(alleles),2)]=allele1
            alleles[seq(2,length(alleles),2)]=allele2
            nameset=unique(alleles)
            nameset=nameset[!is.na(nameset)]
            predmat=matrix(nrow=length(nameset),ncol=ncol(temp2))
            for(i in 1:nrow(predmat)){
                tt=temp2[!is.na(alleles) & alleles==nameset[i],]
                tt=matrix(tt,ncol=ncol(temp2))
                predmat[i,]=colMeans(tt)
            }
            rownames(predmat)=nameset
            predmatallele1=matrix(nrow=length(ourtypes1),ncol=ncol(temp2))
            predmatallele2=matrix(nrow=length(ourtypes1),ncol=ncol(temp2))
            predmatallele1[ourtypes1 %in%nameset]=predmat[ourtypes1[ourtypes1 %in%nameset],]
            predmatallele2[ourtypes2 %in%nameset]=predmat[ourtypes2[ourtypes2 %in%nameset],]
            obsmatallele1=temp2[seq(1,nrow(temp2),2),]
            obsmatallele2=temp2[seq(2,nrow(temp2),2),]
            dist11=rowSums((abs((obsmatallele1)-predmatallele1))>0.9)
            dist12=rowSums((abs((obsmatallele1)-predmatallele2))>0.9)
            dist21=rowSums((abs((obsmatallele2)-predmatallele1))>0.9)
            dist22=rowSums((abs((obsmatallele2)-predmatallele2))>0.9)
            w=1:length(dist11)*0
            w[dist11<dist12]=1
            w[dist11>dist12]= -1
            table(w)
            w1=w
            w=1:length(dist21)*0
            w[dist22<dist21]= 1
            w[dist22>dist21]= -1
            w2=w
            ##print(table(w1,w2))
            phase1b=dist11+dist22
            phase2b=dist12+dist21
            d11=dist11
            d21=dist21
            d12=dist12
            d22=dist22
            ##define phased as simply the nearest, or both alleles the same
            ##non-flipped
            phased1b=(!is.na(phase1b) & !is.na(phase2b) & phase1b<phase2b) |
                (!is.na(ourtypes1) & !is.na(ourtypes2) & ourtypes1==ourtypes2) | 
                (is.na(d21) & !is.na(d12) & d12-d22>2) |  (is.na(d12) & !is.na(d21) & d21-d11>2)
            ##flipped
            phased2b=(!is.na(phase1b) & !is.na(phase2b) & phase1b>phase2b) |
                (is.na(d21) & !is.na(d12) & d22-d12>2) |  (is.na(d12) & !is.na(d21) & d11-d21>2)
            ##does it extend, does it always agree?
            ##print(table((phased1 | phased2), (phased1b | phased2b)))
            ##print(table(phased2b[phased1 | phased2],phased2[phased1 | phased2]))
            ##
            ##
            ##cbind(ourtypes1,ourtypes2,dist11,dist12,dist21,dist22,phase1b,phase2b)[!phased1b & !phased2b,]
            update=!phased1 & !phased2
            phased1[update]=phased1b[update]
            phased2[update]=phased2b[update]
        }
        
        plot(oldphase1,oldphase2,col=1+as.double(phased1)+2*as.double(phased2))

        save(
            phased1,
            phased2,
            ourtypes1,
            ourtypes2,
            file = file_quilt_hla_phase_step_1(outputdir, region)
        )

    }


    NULL
}


hla_phasing_determine_who_to_remove <- function(
    outputdir,
    hla_types_panel,
    reference_sample_file,
    reference_exclude_samplelist_file,
    regions
) {


    ## new code, to process these!!!
    ## used here to define the reference samples essentially, above to identify hla types and phase
    ## load(full_reference_hap_file)
    load(file_quilt_hla_all_haplotypes(outputdir))
    ## this needs to have the HLA types in it
    hlatypes <- read.table(hla_types_panel,header=T,skip=0,sep="\t")
    ##subset of inds in the reference panels
    hlatypes2=hlatypes[match(reference_samples[,1],hlatypes[,3]),]

    ## all individuals that could be used to make a panel, i.e. full set of reference panel individuals
    ## full <- read.table("hrc.chr6.samples",header=T)
    full <- read.table(reference_sample_file,header=T)
    ##full <- read.table("/data/smew1/rdavies/quilt_hla_finalize_with_simon/inputs_2021_03_18/hrc.chr6.samples.reheadered2", header=T)

    ## exclude individuals without phasing information
    for(region in regions){
        ##region specific

        ## load(paste("hla",region,"newphased.out",sep=""))
        load(file_quilt_hla_phase_step_1(outputdir, region))
        
        ## make an input file for getting imputation panels
        ## remove untyped cases (including no HLA type but not "None")
        keep=hlatypes2[phased1==T | phased2==T,3]
        ##pos=match(full[,1],hlatypes[,3])
        newremoves <- full[!full[,1]%in%keep,1]
        newremoves <- as.vector(newremoves)
        
        ## check if additional individuals are to be removed
        if(reference_exclude_samplelist_file != "") {
            removes <- scan(reference_exclude_samplelist_file, what = 'char')
            newremoves2 <- full[
                !full[, 1] %in% keep |
                 full[, 1] %in% removes,
                1
            ]
            newremoves2 <- as.vector(newremoves2)
            newremoves <- newremoves2
        }
        print(length(newremoves))
        cat(
            newremoves,
            file = file_quilt_hla_after_step_1_who_to_remove(outputdir, region),
            sep = "\n"
        )
    }

}



hla_perform_step_2_phasing <- function(
    outputdir,
    regions,
    hla_types_panel
) {


    load(file_quilt_hla_all_haplotypes(outputdir))
    hlatypes <- read.table(hla_types_panel, header=T,skip=0,sep="\t")
    hlatypes2 <- hlatypes[match(reference_samples[,1],hlatypes[,3]),]
    
    
    ## not entirely sure where this goes

    
    ##below is commented out as deals with excluding five groups for testing purposes
    ##for(region in c("A","B","C","DQB1","DRB1")){
    ##load(paste("hrc.chr6.hla.",region,".hlatypedexcludefivepop.RData",sep=""))
    ##sample_names=ref_samples
    ##finds=match(sample_names[,2],hlatypes2[,3])
    ##load(paste("hla",region,"newphased.out",sep=""))
    ##hlatypes3=hlatypes2[finds,]
    ##phased1=phased1[finds]
    ##phased2=phased2[finds]
    ##al1=as.vector(hlatypes3[,paste("HLA.",region,".1",sep="")])
    ##al2=as.vector(hlatypes3[,paste("HLA.",region,".2",sep="")])
    ##allele1=as.vector(al1)
    ##allele2=as.vector(al2)
    ##allele1[phased1==F]=al2[phased1==F]
    ##allele2[phased1==F]=al1[phased1==F]
    ##names(allele1)=paste(as.vector(hlatypes3[,3]),"_a1",sep="")
    ##names(allele2)=paste(as.vector(hlatypes3[,3]),"_a2",sep="")
    ##alleles=c(allele1,allele2)
    ##nhaps=dim(sample_names)[1]*2
    ##vv=1:nhaps
    ##vv[seq(1,length(vv),2)]=paste(sample_names[,2],"_a1",sep="")
    ##vv[seq(2,length(vv),2)]=paste(sample_names[,2],"_a2",sep="")
    ##copya1=alleles[vv]
    ##hlahaptypes=copya1
    ##save(hlahaptypes,file=paste("hla",region,"haptypesexcludefivepops.out",sep=""))
    ##}
    

    

    for(region in regions){

        print_message(paste0("Integrate results from region:", region))        
        ##panel hrc.chr6.hla. is the RData file of the reference panel for this region, must contain the individuals whose HLA types we seek
        load(file_quilt_hla_specific_haplotypes(outputdir, region))
        ## load(gsub("*", region, local_reference_hap_file, fixed = TRUE))
        ## load(paste("hrc.chr6.hla.",region,".hlatyped.RData",sep=""))
        
        sample_names <- reference_samples

        ## identify who corresponds to our HLA types
        finds <- match(sample_names[,2],hlatypes2[,3])

        ## load in our phasing inference for the HLA types, matching the reference panel phasing
        load(file = file_quilt_hla_phase_step_1(outputdir, region))
        ## load(paste("hla",region,"newphased.out",sep=""))

        ## reorder to match reference panel
        hlatypes3=hlatypes2[finds,]
        phased1=phased1[finds]
        phased2=phased2[finds]

        ##two alleles for each case
        al1=as.vector(hlatypes3[,paste("HLA.",region,".1",sep="")])
        al2=as.vector(hlatypes3[,paste("HLA.",region,".2",sep="")])
        allele1=as.vector(al1)
        allele2=as.vector(al2)
        allele1[phased1==F]=al2[phased1==F]
        allele2[phased1==F]=al1[phased1==F]
        names(allele1)=paste(as.vector(hlatypes3[,3]),"_a1",sep="")
        names(allele2)=paste(as.vector(hlatypes3[,3]),"_a2",sep="")
        alleles=c(allele1,allele2)
        nhaps=dim(sample_names)[1]*2
        vv=1:nhaps
        vv[seq(1,length(vv),2)]=paste(sample_names[,2],"_a1",sep="")
        vv[seq(2,length(vv),2)]=paste(sample_names[,2],"_a2",sep="")
        copya1=alleles[vv]

        ##note that the below will have values inherited from the data table itself - in particular this could mean ambiguous values, None, or missing data
        ## however these might be excluded at the ref panel construction stage

        hlahaptypes <- copya1
        save(
            hlahaptypes,
            file = file_quilt_hla_phased_haplotypes(outputdir, region)
        )

    }


    NULL

}
