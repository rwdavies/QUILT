getnames <- function(
    startpos,
    endpos,
    ourname = ":",
    this
){
    ourstuff=this[startpos:endpos]
    ourstuff=this[startpos:endpos]
    seqstarts=grep(ourname,ourstuff)
    ttt=gsub(ourname,"",ourstuff[seqstarts])
    check=substring(ttt,1,1)
    seqstarts=seqstarts[check %in% 0:9]
    lnames=ourstuff[seqstarts]
    return(lnames)
}

getseqs <- function(
    startpos,
    endpos,
    ourname = ":",
    this
){
    ourstuff=this[startpos:endpos]
    seqstarts=grep(ourname,ourstuff)
    ttt=gsub(ourname,"",ourstuff[seqstarts])
    check=substring(ttt,1,1)
    seqstarts=seqstarts[check %in% 0:9]+1
    ourseq=rep("",length(seqstarts))
    names(ourseq)=ourstuff[seqstarts[1:length(seqstarts)]-1]
    seqstarts=c(seqstarts,length(ourstuff)+2)
    for(j in 1:(length(seqstarts)-1)) ourseq[j]=paste(ourstuff[seqstarts[j]:(seqstarts[j+1]-2)],collapse="")
    return(ourseq)
}


getvars2 <- function(refrow=temp[ourrow,],testrow,ourpos){
    res=matrix(nrow=0,ncol=3)
    ##first insertions and snps
    p1=unique(ourpos[testrow!=refrow])
    if(length(p1)){

        for(i in p1){
            vv=paste(refrow[ourpos==i],collapse="")
            ww=paste(testrow[ourpos==i],collapse="")
            if(ww!=vv & !length(grep("[*]",ww))) res=rbind(res,c(i,vv,ww))
        }
        ##return(res)
        ##}
        ##print(res[1:10,])
#####deletions have e.g. CCC vs C.. so need to collapse
        j=1
        res2=matrix(nrow=0,ncol=3)
        ccc=gsub("[.]","",res[,3])
        ddd=gsub("[.]","",res[,2])
        rm=rep(0,nrow(res))
        while(j <nrow(res)){
            
            if(ccc[j]=="" & ddd[j]!=""){
                ##print("DEL")
#####have a deletion in new
                thispos=res[j,1]
                all1=res[j,2]
                all2=res[j,3]
                rm[j]=1;
                j=j+1
                while(j <=nrow(res) & ccc[j]==""){
                    all1=paste(all1,res[j,2],sep="")
                    all2=paste(all2,res[j,3],sep="")
                    rm[j]=1
                    j=j+1
                    ##print(j)
                }
                res2=rbind(res2,c(thispos,all1,all2))
                ##print(res2)
            } 
            j=j+1
        }
        res=matrix(rbind(res[rm==0,],res2),ncol=3)
        if(nrow(res)){
            ##print(nrow(res))
            res[,2]=gsub("[.]","",res[,2])
            res[,3]=gsub("[.]","",res[,3])
            res=res[res[,2]!=res[,3],]
            res=matrix(res,ncol=3)
            res=res[order(as.double(res[,1])),]
        }
    }
    
    return(res)
}


## test2=getvars2(temp[1,],temp[616,],ourpos)
## now slightly modify
process <- function(varset){
    changed=1
    while(changed==1){
        changed=0
        test=substring(varset[,2],1,1)
        test2=substring(varset[,3],1,1)
        changed=as.double(sum(test==test2)>0)
        if(changed){
            alter=which(test==test2)
            ##print(length(alter))
            for(k in 1:length(alter)){
		varset[alter[k],1]=as.double(varset[alter[k],1])+1
		varset[alter[k],2]=substring(varset[alter[k],2],2,nchar(varset[alter[k],2]))
		varset[alter[k],3]=substring(varset[alter[k],3],2,nchar(varset[alter[k],3]))
                
            }
        }
    }
    
######other end
    changed=1
    while(changed==1){
        changed=0
        test=substring(varset[,2],nchar(varset[,2]),nchar(varset[,2]))
        test2=substring(varset[,3],nchar(varset[,3]),nchar(varset[,3]))
        changed=as.double(sum(test==test2)>0)
        if(changed){
            alter=which(test==test2)
            ##print(length(alter))
            for(k in 1:length(alter)){
                varset[alter[k],2]=substring(varset[alter[k],2],1,(nchar(varset[alter[k],2])-1))
		varset[alter[k],3]=substring(varset[alter[k],3],1,(nchar(varset[alter[k],3])-1))
                
            }
        }
    }
    return(varset)
}





get_kmers_from_one_file <- function(
    i_file,
    ourfiles,
    outputdir
) {

    zztop <- i_file ##  dangit Simon
    n <- length(ourfiles)
    ##if (zztop %in% (match(1:10, ceiling(10 * (1:n / n))))) {
    print_message(paste0("Processing file ", zztop, " out of ", length(ourfiles)))
    ##}
    
    curkmers=matrix(nrow=0,ncol=3)

    this <- scan(
        file.path(
            outputdir,
            "alignments",
            paste0(ourfiles[zztop], "_gen.txt")
        ),
        what = 'char',
        quiet = TRUE
    )
    
    temp=grep("Please",this)
    this=this[1:(temp-1)]
    
    starts=grep("gDNA",this)
    ll=getseqs(starts[1]+2,starts[2]-1,paste(ourfiles[zztop],"[*]",sep=""), this = this)
    starts=c(starts,length(this)+2)
    
    ## amount to trim for first codon
    ## offset=as.double(this[starts[1]+1])
    
    for(k in 2:(length(starts)-1)) {
        ll=paste(
            ll,
            getseqs(starts[k]+2,starts[k+1]-1,paste(ourfiles[zztop],"[*]",sep=""), this = this),
            sep=""
        )
        ##print(k)
    }
    names(ll)=getnames(starts[1]+2,starts[2]-1,paste(ourfiles[zztop],"[*]",sep=""), this = this)
    
    ## find a match
    
    temp=matrix(nrow=length(ll),ncol=nchar(ll[1]))
    
    for(i in 1:ncol(temp)) temp[,i]=substring(ll,i,i)
    for(i in 1:ncol(temp)) temp[temp[,i]=="-",i]=temp[1,i]
    
    spos=which(temp[1,]=="|")[1]
    temp=temp[,(spos+1):ncol(temp)]
    temp=temp[,temp[1,]!="|"]
    
    newseqs=vector(length=nrow(temp))
    for(i in 1:length(newseqs)){ newseqs[i]=paste(temp[i,],collapse="")
        newseqs[i]=unlist(gsub("\\.","",newseqs[i]))
    }
    
    for(curpos in 1:(max(nchar(newseqs))-9)){
        
        vv=unique(substring(newseqs,curpos,curpos+9))
        qq=grep("\\*",vv)
        if(length(qq) )
            vv=vv[-qq]
        vv=vv[nchar(vv)==10]
        if(length(vv))curkmers=rbind(curkmers,cbind(vv,curpos,ourfiles[zztop]))
    }

    curkmers

}


make_and_save_hla_all_alleles_kmers <- function(
    outputdir,
    all_hla_regions,
    hla_gene_information,
    nCores
) {

    ourfiles <- all_hla_regions
    
    print_message("Begin making HLA all alleles kmers file")
    ## code to make HLAallalleleskmers.out

    out <- mclapply(
        1:length(ourfiles),
        FUN = get_kmers_from_one_file,
        ourfiles = ourfiles,
        outputdir = outputdir,
        mc.cores = nCores
    )
    check_mclapply_OK(out)
    
    kmers <- out[[1]]
    for(i in 2:length(out)) {
        kmers <- rbind(kmers,out[[i]])
    }
    ## kmers <- get_kmers_for_hla_alleles_all_of_them(ourfiles)
    
    newnames=unique(kmers[,1])
    newkmers=matrix(nrow=length(newnames),ncol=3)
    newkmers[,1]=newnames
    rownames(newkmers)=newnames
    
    temps=kmers[order(kmers[,1]),]
    starts=c(1,which(temps[2:nrow(kmers),1]!=temps[1:(nrow(kmers)-1),1])+1)
    ends=c(starts[2:length(starts)]-1,nrow(temps))
    
    newnames=unique(temps[,1])
    newkmers=matrix(nrow=length(newnames),ncol=3)
    newkmers[,1]=newnames
    rownames(newkmers)=newnames
    for(i in 1:length(newnames)) newkmers[i,2]=paste(temps[starts[i]:ends[i],2],collapse=",")
    for(i in 1:length(newnames)) newkmers[i,3]=paste(temps[starts[i]:ends[i],3],collapse=",")

    ## easy enough to save hla_gene_information here
    save(
        kmers,
        newkmers,
        hla_gene_information,
        file = file_quilt_hla_all_alleles_kmers(outputdir)
    )

    print_message("Done making HLA all alleles kmers file")
    ##
    ## end of code to make HLAallalleleskmers.out
    ##

    return(NULL)
}









make_single_snpformatalleles <- function(
    outputdir,
    hla_region,
    temp,
    ll,
    ourpos,
    ourrow
) {
        
    ##varset will store an initial set of variant positions
    ##print("Calling initial variants, sequences analysed:")
    varset=matrix(nrow=0,ncol=4)
    for(i in 1:nrow(temp)) {
        ##if(!i%%50) print(i)
        cc=getvars2(temp[ourrow,],temp[i,],ourpos)
        if(length(cc)) varset=rbind(varset,cbind(matrix(cc,ncol=3),i))
    }
    ##print("...done: filtering variants")
    
    cccc=paste(varset[,1],rep("W",nrow(varset)),varset[,2],rep("W",nrow(varset)),varset[,3],rep("W",nrow(varset)),sep="")
    length(unique(cccc))
    cccc=unique(cccc)
    cccc2=cccc
    cccc=strsplit(cccc,"W")
    newset=matrix(nrow=length(cccc),ncol=3)
    for(i in 1:nrow(newset)) newset[i,]=cccc[i][[1]]
    newset[1:5,]
    dim(newset)
    
    ##modify to simpler way of calling mutations
    varset2=process(newset)
    knownvars=varset2
    knownvars=knownvars[order(as.double(knownvars[,1])),]
    ##this has removed identical bases
    ##$filter variants where can co-occur depending on calling protocol
    
    
    newvars=matrix(nrow=0,ncol=ncol(knownvars))
    all1=knownvars[,2]
    all2=knownvars[,3]
    for(i in 1:length(all1)) if(nchar(all1[i])<nchar(all2[i])){
                                 
                                 all1[i]=paste(paste(all1[i],paste(rep(".",nchar(all2[i])-nchar(all1[i])),collapse=""),sep=""),collapse="")
                             }
    for(i in 1:length(all1)) if(nchar(all2[i])<nchar(all1[i])){
                                 
                                 all2[i]=paste(paste(all2[i],paste(rep(".",nchar(all1[i])-nchar(all2[i])),collapse=""),sep=""),collapse="")
                             }
    spos=as.double(knownvars[,1])
    epos=spos+nchar(all1)-1
    
    rm=rep(0,nrow(knownvars))
    toto=rm
    for(i in 1:nrow(knownvars)){
        
        if(nchar(all1[i])>1000 | nchar(all2[i])>1000) rm[i]=1
        else{
            
            overlaps=which(spos<=spos[i] & epos>=epos[i])
            toto[i]=length(overlaps)
            
######first just take subset case, only worry about shorter indels
            
            if(length(overlaps)>1){
                overlaps=overlaps[overlaps!=i]
                for(k in 1:length(overlaps)){
                    if(length(grep(all1[i],all1[overlaps[k]])) & length(grep(all2[i],all2[overlaps[k]])) & (nchar(all1[overlaps[k]])-nchar(all1[i])<=25 | nchar(knownvars[i,2])!=nchar(knownvars[i,3]))){
                        
                        cccc=substring(all1[overlaps[k]],1:(nchar(all1[overlaps[k]])-nchar(all1[i])+1),nchar(all1[i]):nchar(all1[overlaps[k]]))
                        cccc2=substring(all2[overlaps[k]],1:(nchar(all1[overlaps[k]])-nchar(all1[i])+1),nchar(all1[i]):nchar(all1[overlaps[k]]))
                        
                        if(sum(all1[i]==cccc & all2[i]==cccc2)) rm[i]=1
                        
                    }
                    if(!length(grep(all1[i],all2[overlaps[k]])) & length(grep(all2[i],all1[overlaps[k]])) & 
                       (nchar(all1[overlaps[k]])-nchar(all1[i])<=25 | nchar(knownvars[i,2])!=nchar(knownvars[i,3]))){
                        
                        cccc=substring(all1[overlaps[k]],1:(nchar(all1[overlaps[k]])-nchar(all1[i])+1),nchar(all1[i]):nchar(all1[overlaps[k]]))
                        cccc2=substring(all2[overlaps[k]],1:(nchar(all1[overlaps[k]])-nchar(all1[i])+1),nchar(all1[i]):nchar(all1[overlaps[k]]))
                        
                        if(sum(all2[i]==cccc & all1[i]==cccc2)) rm[i]=1
                        
                        
                    }
                }
            }
            ##partial overlap case
            overlaps=which(spos[i]<spos & epos[i]>=spos & epos[i]<epos)
            if(length(overlaps)) rm[i]=1
            overlaps=which(spos[i]>spos & spos[i]<=epos & epos[i]>epos)
            if(length(overlaps)) rm[i]=1
            
            
        }
        ##print(length(overlaps))
    }
    
    knownvarsfiltered=knownvars[rm==0,]
    ##print("....done:Making variant calls, percentage complete:")
    ##this is a final set of variant positions after filtering
    
    ##use to make calls
    tempmat=matrix(nrow=0,ncol=7)
    ##resmat now gives the called type at each of the variant positions
    resmat=matrix(0,nrow=nrow(temp),ncol=nrow(knownvarsfiltered))
    
    for(q in 1:nrow(knownvarsfiltered)){
        ## if(!q%%50) print(q/nrow(knownvarsfiltered)*100)
        a1=knownvarsfiltered[q,2]
        
        a2=knownvarsfiltered[q,3]
        pos=as.double(knownvarsfiltered[q,1])
        start=pos
        end=pos+nchar(a1)-1
        flag=0
        if(end<start){ start=end
            flag=1
        }
        s1=which(ourpos>=start & ourpos<=end)
        if(flag==1){
            s1=s1[2:length(s1)]
        }
        checkmat=temp[,s1]
        checkmat[checkmat=="."]=""
        checkmat=matrix(checkmat,ncol=length(s1))
        types=rep("",nrow(checkmat))
        for(i in 1:ncol(checkmat)){
            types=paste(types,checkmat[,i],sep="")
        }
        resmat[types==a2,q]=1
        rs=grep("[*]",types)
        if(length(rs)) resmat[rs,q]=-1
        ##print(c(a1,a2))
        ##print(table(types))
        ##print(c(q,sum(resmat[,q])))
        tempmat=rbind(tempmat,c(a1,a2,length(unique(types)),unique(types)[1:2],q,sum(resmat[,q])))
    }
    ##print("...done")
    ##print("Final clean-up...")
    ##can be an issue where some variants are never seen, remove these (it is because there is local similarity and more than one correct local aligment so may be inconsisently called)
    
    knownvarsfiltered=knownvarsfiltered[colSums(resmat>0)>0,]
    resmat=resmat[,colSums(resmat>0)>0]
    rownames(resmat)=names(ll)
    ##print("...done, saving to file")
    save(
        resmat,
        knownvarsfiltered,
        file = file_quilt_hla_snpformatalleles(outputdir, hla_region)
    )

    1
    
}

make_single_hla_full_alleles_filled_in <- function(
    outputdir,
    hla_region,
    temp,
    ll,
    ourpos,
    resmat,
    knownvarsfiltered
) {

    temp2=list(ol=knownvarsfiltered,ourallelemat=resmat)
    
    sample=cbind(rownames(temp2$ourallelemat),rownames(temp2$ourallelemat),0)
    colnames(sample)=c("ID_1","ID_2","missing")
    
    pos=cbind(paste("chr6:",temp2$ol[,1],sep=""),temp2$ol[,1],temp2$ol[,2],temp2$ol[,3])
    colnames(pos)=c("id","position","a0","a1")
    
    ##print("Filtering variants using missing data...")
    
    haps=t(temp2$ourallelemat)
    cond=pos[,3]%in%c("A","C","G","T") & pos[,4]%in%c("A","C","G","T") & rowMeans(haps==-1)<0.1
    haps=haps[cond ,]
    pos=pos[cond,]
    
    ##print("...done")
    
    qq=range(as.double(pos[,2]))
    qq[1]=qq[1]-10
    qq[2]=qq[2]+10
    
    temp=temp[,ourpos>=qq[1] & ourpos<=qq[2]]
    ourpos=ourpos[ourpos>=qq[1] & ourpos<=qq[2]]
    
    ##fill in
    
    ##space separated
    oldvalue=nrow(temp)*ncol(temp)
    newvalue=sum(temp=="*")
    
    while(oldvalue>newvalue){
        
        hapstore=haps
        tempstore=temp
        
        ##guess missing values
        ##print("Total missing still:")
        ##print(sum(temp=="*"))
        ##print("Total columns: ")
        ##print(ncol(haps))
        ##print("Processing column to find nearest neighbour:")
        
        for(i in 1:ncol(haps)){
            
            if(sum(temp[i,]=="*")){
                
                ## if(!i%%20) {print(c(i,sum(temp[i,]=="*")))}
####use haplotypes to find match
                cond2=(temp[i,]!="*")
                cond=haps[,i]!=-1
                c1=hapstore[cond,i]
                c2=hapstore[cond,]-hapstore[cond,i]
                dist=0;dist2=0;
                if(sum(c1>0)>1) dist=colSums(c2[c1>0,]== -1)
                if(sum(c1==0)>1) dist2=colSums(c2[c1==0,]==1)
                if(sum(c1>0)==1) dist=as.double(c2[c1>0,]== -1)
                if(sum(c1==0)==1) dist2=as.double(c2[c1==0,]==1)
                dist=dist+dist2
                newalleles=matrix(tempstore[dist==min(dist[-i]),!cond2],ncol=sum(!cond2))
                ##newhaps=matrix(hapstore[!cond,dist==min(dist[-i])],nrow=sum(!cond))
                allelemat=matrix(nrow=6,ncol=sum(!cond2))
                rownames(allelemat)=c("A","C","G","T","-","*")
                for(k in 1:nrow(allelemat)) allelemat[k,]=colSums(newalleles==rownames(allelemat)[k])
                best=1:ncol(allelemat);for(k in 1:ncol(allelemat)) best[k]=rownames(allelemat)[1:5][which(allelemat[1:5,k]==max(allelemat[1:5,k]))[1]]
                best[colSums(allelemat[1:5,, drop = FALSE])==0]="*"
                temp[i,!cond2]=best	
                
                
            }
        }
        oldvalue=newvalue
        newvalue=sum(temp=="*")
        
    }
    ##print("....done")
    
    ##guess missing values
    ## print_message("Inferring remaining missing values")
    
    ##print(sum(haps==-1))
    hapstore=haps
    tempstore=temp
    ##print("Processing column:")
    for(i in 1:ncol(haps)){
        
        if(sum(temp[i,]=="*")){
            
            ##print(i)
            cond=haps[,i]!=-1
            cond2=(temp[i,]=="*")
            alter=which(cond2==T)
            c1=hapstore[cond,i]
            c2=hapstore[cond,]-hapstore[cond,i]
            dist=0;dist2=0;
            if(sum(c1>0)>1) dist=colSums(c2[c1>0,]== -1)
            if(sum(c1==0)>1) dist2=colSums(c2[c1==0,]==1)
            if(sum(c1>0)==1) dist=as.double(c2[c1>0,]== -1)
            if(sum(c1==0)==1) dist2=as.double(c2[c1==0,]==1)
            dist=dist+dist2
            for(j in alter){
                t1=temp[,j]
                dd=dist[t1!="*"]
                ##print(min(dd))
                possalleles=t1[t1!="*"][dd==min(dd)]
                cc=1:5*0;names(cc)=c("A","C","G","T","-")
                for(k in 1:5) cc[k]=sum(possalleles==names(cc)[k])
                temp[i,j]=names(cc)[which(cc==max(cc))[1]]
                
            }
        }
    }
    ##print("...done")
    
    ##print("Saving...")
    fullalleles=temp
    rownames(fullalleles)=names(ll)
    
    save(
        ourpos,
        fullalleles,
        file = file_quilt_hla_full_alleles_filled_in(outputdir, hla_region)
    )

}



get_and_reformat_gen_alignments_for_hla_region <- function(
    outputdir,
    hla_region,
    hla_strand
) {

    this <- scan(
        file.path(
            outputdir,
            "alignments",
            paste0(hla_region, "_gen.txt")
        ),
        what = 'char',
        quiet = TRUE
    )
    
    temp=grep("Please",this)
    this=this[1:(temp-1)]
    
    starts=grep("gDNA",this)
    ll=getseqs(
        starts[1]+2,
        starts[2]-1,
        paste(hla_region,"[*]",sep=""),
        this = this
    )
    starts=c(starts,length(this)+2)
    
    ##amount to trim for first codon
    ##offset=as.double(this[starts[1]+1])
    
    for(k in 2:(length(starts)-1)) {
        ll <- paste0(
            ll,
            getseqs(
                starts[k]+2,
                starts[k+1]-1,
                paste(hla_region,"[*]",sep=""),
                this = this
            )
        )
        ##print(k)
    }
    names(ll)=getnames(starts[1]+2,starts[2]-1,paste(hla_region,"[*]",sep=""), this = this)
    
    
    ##temp is matrix of alignments, stranded appropriately
    
    temp=matrix(nrow=length(ll),ncol=nchar(ll[1]))
    
    for(i in 1:ncol(temp)) temp[,i]=substring(ll,i,i)
    for(i in 1:ncol(temp)) temp[temp[,i]=="-",i]=temp[1,i]

    ## OK looks like the logic was OK here (Simon original)
    ## not sure what the numbering is then if it is sometimes wrong!
    ## in any case, easy to calculate
    spos1 <- which(temp[1,]=="|")[1]
    offset <- as.numeric(this[grep("gDNA", this)[1] + 1]) * -1
    ## want the offset + 1 base after removing periods
    spos2 <- which(cumsum(temp[1, ] != ".") == (offset + 1))
    
    n <- c(
        spos1 = spos1,
        spos2 = spos2,
        nperiod1 = sum(temp[1, 1:spos1] == "."),        
        nperiod2 = sum(temp[1, 1:spos2] == "."),
        nbar1 = sum(temp[1, 1:spos1] == "|"),        
        nbar2 = sum(temp[1, 1:spos2] == "|")
    )
        
        
    ## remove before the start of CDS
    temp <- temp[,(spos1 + 1):ncol(temp)]
    temp <- temp[, temp[1,] != "|"]
    
    if(hla_strand != 1){
        temp[temp=="A"]="t"
        temp[temp=="C"]="g"
        temp[temp=="G"]="c"
        temp[temp=="T"]="a"
        temp[temp=="a"]="A"
        temp[temp=="c"]="C"
        temp[temp=="g"]="G"
        temp[temp=="t"]="T"
        temp=temp[,ncol(temp):1]
    }

    return(
        list(
            ll = ll,
            temp = temp,
            n = n
        )
    )
}


make_single_hla_full <- function(
    outputdir,
    hla_region,
    fullalleles
) {

    kmers=vector(length=0)
    positions=vector(length=0)
    xx=vector(length=0)

    ##print("Making kmers by allele:")
    for(i in 1:nrow(fullalleles)){
        ##if(!i%%20) print(i)
        zz=1:(ncol(fullalleles)-9)
        ww=fullalleles[i,]
        zz=zz[ww!="."]
        ww=ww[ww!="."]
        dd=paste(ww,collapse="")
        ee=substring(dd,1:(nchar(dd)-9),10:nchar(dd))
        ##unique only
        hh=paste(ee,zz)
        cond=!(hh %in% xx)
        if(sum(cond)){
            ee=ee[cond]
            zz=zz[cond]
            hh=hh[cond]
            kmers=c(kmers,ee)
            positions=c(positions,zz)
            xx=c(xx,hh)
        }
    }
    xx=xx[!is.na(kmers)]
    positions=positions[!is.na(kmers)]
    kmers=kmers[!is.na(kmers)]
    ##can save this
    ##print("..done")

    ##uniquify
    kk=unique(kmers)
    names(kk)=kk
    for(i in 1:length(kk)) kk[i]=paste(positions[kmers==names(kk)[i]],collapse=",")

    ##for each allele, what base is column i in alignment (0 means gap)
    lookup=matrix(0,nrow=nrow(fullalleles),ncol=ncol(fullalleles))
    ##for each allele, where is base i in columns of alignment
    revlookup=lookup

    ##print("Reverse-lookup calculation")
    for(i in 1:nrow(fullalleles)){
        t1=1:ncol(fullalleles)
        t2=fullalleles[i,]
        t1=t1[t2!="."]
        lookup[i,t1]=1:length(t1)
        revlookup[i,1:length(t1)]=t1
    }
    ##print("....done")
    ##print("Saving kmers...")
    save(
        kmers,
        positions,
        xx,
        revlookup,
        lookup,
        fullalleles,
        file = file_quilt_hla_full(outputdir, hla_region) 
    )
    ##print("...done")

    0
    
}




make_and_save_hla_files_for_imputation <- function(
    outputdir,
    hla_regions_to_prepare,
    hla_gene_information,
    ref_fasta,
    nCores
) {

    ## code to make hla*snpformatalleles.out
    print_message("Begin making HLA snp format alleles file")

    out <- mclapply(1:length(hla_regions_to_prepare), mc.cores = nCores, function(i_region){

        hla_region <- hla_regions_to_prepare[i_region]

        print_message(paste0("Working on region:", hla_region))

        m <- match(paste0("HLA-", hla_region), hla_gene_information[, "Name"])
        if (is.na(m)) {
            msg <- paste0(
                "Cannot match ", shQuote(paste0("HLA-", hla_region)), " to any of the available hla genes:",
                paste0(hla_gene_information[, "Name"], collapse= ", ")
            )
            stop(msg)
        }

        hla_strand <- hla_gene_information[m, "Strand"]
        if (hla_strand == 1) {
            genome_pos <- hla_gene_information[m, "Start"]
        } else {
            genome_pos <- hla_gene_information[m, "End"]
        }

        out <- get_and_reformat_gen_alignments_for_hla_region(
            outputdir = outputdir,
            hla_region = hla_region,
            hla_strand = hla_strand
        )
        ll <- out[["ll"]]
        temp <- out[["temp"]]

        aligned <- sum(temp[1,] %in% c("A","C","G","T"))
        
        if(hla_strand == 1) {
            start <- genome_pos
        } else {
            start <- genome_pos - (aligned - 1)
        }
        end <- start + aligned - 1

        reference_allele <- determine_reference_genome_hla_allele(
            ref_fasta = ref_fasta,
            chr = hla_gene_information[1, "Chr"],
            start = start,
            end = end,
            ll = ll,
            temp = temp,
            hla_strand = hla_strand
        )

        print_message(paste0("Determined automatically that the reference genome contains allele:", reference_allele))

        ##
        ## I don't fully understand this in simon's code
        ## but now need positions (somehow) with respect to this allele
        ## does not align obviously to genome, but works later!
        ##
        first_row <- which(names(ll) == reference_allele)
        aligned <- sum(temp[first_row,] %in% c("A","C","G","T"))
        
        if(hla_strand == 1) {
            start <- genome_pos
        } else {
            start <- genome_pos - (aligned - 1)
        }

        ourpos <- rep(0,ncol(temp))
        ourpos[1] <- start
        
        curpos <- ourpos[1]
        ourrow <- first_row
        for(i in 2:length(ourpos)){
            if(temp[ourrow,i] %in% c("A","C","G","T","*")) {
                curpos <- curpos + 1
            }
            ourpos[i] <- curpos
        }
        ## ourpos gives positions relative to the genome reference sequence

        ##
        ## now make the snpformatalleles and the other one
        ##
        print_message(paste0("Make SNP format alleles file for:", hla_region))        
        out <- make_single_snpformatalleles(
            outputdir = outputdir,
            hla_region = hla_region,
            temp = temp,
            ll = ll,
            ourpos = ourpos,
            ourrow = ourrow
        )
        load(file_quilt_hla_snpformatalleles(outputdir, hla_region))
        
        print_message(paste0("Make full alleles filled in file for:", hla_region))                
        make_single_hla_full_alleles_filled_in(
            outputdir = outputdir,
            hla_region = hla_region,
            temp = temp,
            ll = ll,
            ourpos = ourpos,
            resmat = resmat,
            knownvarsfiltered = knownvarsfiltered
        )

        load(file_quilt_hla_full_alleles_filled_in(outputdir, hla_region))

        print_message(paste0("Make final files:", hla_region))        
        make_single_hla_full(
            outputdir = outputdir,
            hla_region = hla_region,
            fullalleles = fullalleles
        )
        
        print_message(paste0("Done working on region:", hla_region))

        0
        
    })

    check_mclapply_OK(out)
    
    print_message("Done making HLA snp format alleles file")    

}





















per_entry_in_temp_check_match <- function(temp, refseq) {
    t(sapply(1:nrow(temp), function(i) {
        x <- temp[i, ] != "."
        a <- temp[i, ][x]
        b <- refseq[1:sum(x)]
        c <- a != b
        c(sum(c, na.rm = TRUE), sum(!is.na(c)))
    }))
}


determine_reference_genome_hla_allele <- function(
    ref_fasta,
    chr,
    start,
    end,
    ll,
    temp,
    hla_strand,
    error_check = TRUE
) {

    if (1 == 0) {
        ref_fasta <- "/data/smew1/rdavies/quilt_hla_2021_12_24_3430//GRCh38_full_analysis_set_plus_decoy_hla.fa"
    }
        
    command <- paste0(
        "samtools faidx ",
        ref_fasta, " ",
        chr, ":", start, "-", end
    )
    seq <- system(command, intern = TRUE)
    refseq <- unlist(strsplit(seq[-1], ""))

    if (hla_strand == -1) {
        refseq <- refseq[length(refseq):1]
        temp <- temp[, ncol(temp):1]
    }
    
    check_entries <- per_entry_in_temp_check_match(temp, refseq)

    x <- which.min(check_entries[, 1])
    if (sum(check_entries[, 1] == 0) == 0) {
        if (error_check) {
            f <- stop
        } else {
            f <- warning
        }
        f(paste0("Error in automatically determining the HLA allele carried by the reference sequence. No allele is a perfect match to the reference sequence. The nearest match is ", names(ll)[x], " which has ", check_entries[x, 1], " disagreements out of ", check_entries[x, 2], " positions checked"))
    }

    names(ll)[which.min(check_entries[, 1])]
    
}



get_hla_gene_information <- function(
    table_file,
    all_hla_regions,
    chr,
    what = "refseq"
) {
    ## table_file <- "hla_ancillary_files/GENCODEv38.hg38.chr6.26000000.34000000.txt.gz"
    ## table_file <- "hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz"
    ## genes <- read.table(gencode_table_file, header = TRUE, comment.char = "")
    if (what == "refseq") {
        genename <- "name2"
        chrom <- "chrom"
        start <- "cdsStart"
        end <- "cdsEnd"
    } else {
        genename <- "geneName"
        chrom <- "X.chrom"
        start <- "thickStart"
        end <- "thickEnd"
    }
    genes <- read.table(table_file, header = TRUE, comment.char = "")    
    ##
    x1 <- match(paste0("HLA-", all_hla_regions), genes[, genename])    
    x2 <- match(paste0(all_hla_regions), genes[, genename])
    w <- is.na(x1) & !is.na(x2)
    x1[w] <- x2[w]
    ## 
    strand <- c(1, -1)[match(genes[x1[!is.na(x1)], "strand"], c("+", "-"))]
    hla_gene_regions_new <- data.frame(
        Name = paste0("HLA-", all_hla_regions[!is.na(x1)]),
        Chr = genes[x1[!is.na(x1)], chrom],
        Start = genes[x1[!is.na(x1)], start] + 1,
        End = genes[x1[!is.na(x1)], end],
        Strand = strand
    )
    hla_gene_regions_new
}


if (1 == 0) {

    ## table_file <- 
    ## 
    
    get_hla_gene_information(
        table_file = "hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        all_hla_regions,
        chr,
        what = "refseq"
    ) [1:3, ]

    get_hla_gene_information(
        table_file = "hla_ancillary_files/GENCODEv38.hg38.chr6.26000000.34000000.txt.gz",
        all_hla_regions,
        chr,
        what = "gen"
    ) [1:3, ]
    
    
}
