library("STITCH")

## input
recomb_dir <- commandArgs(trailingOnly = TRUE)[1]
panel <- commandArgs(trailingOnly = TRUE)[2]
chr <- commandArgs(trailingOnly = TRUE)[3]

if (1 == 0) {

    ## e.g.
    recomb_dir <- "/data/smew1/rdavies/external/"
    panel <- "CEU"
    chr <- "6"

}

dir.create(recomb_dir)
setwd(recomb_dir)

liftOver <- file.path(recomb_dir, "liftOver")
chain <- file.path(recomb_dir, "hg19ToHg38.over.chain.gz")


##
## step 1 - in R, load in rates, output for liftOver
##      
ratesOri <- read.table(
    paste0(recomb_dir, "/", panel, "/", panel, "-", chr, "-final.txt.gz"),
    header = TRUE
)
n <- dim(ratesOri)[1]

## so for example, row 1 = 30, 0.05, 1, row2 = 40, X, 1+(40-30)*0.05/1e6
## interval 1 is 30-39, 30-40 bedFormat
## rate for each interval is 
out <- cbind(
    paste("chr",chr,sep=""),
    ratesOri[-n,1],
    ratesOri[-1,1],
    ratesOri[-n,2],"X","+"
)
##
reformated_file <- paste0(recomb_dir, "/", panel, "/", panel,"-", chr, "_cleaned_cm.hapmapFormat.forLiftOver.txt")
write.table(
    out,
    file = reformated_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)
out_stem <- paste0(recomb_dir, "/", panel, "/", panel, "-", chr, "_b38")

##
## step 2 - run liftOver
##
system(
    paste0(
        liftOver, " ",
        reformated_file, " ",
        chain, " ",
        out_stem, ".txt ",
        out_stem, ".unmapped.txt "
    )
)
system(paste0("gzip -f -1 ", out_stem, ".txt"))
system(paste0("gzip -f -1 ", out_stem, ".unmapped.txt"))

##
## step 3 - read in, fix
##
new <- read.table(paste0(out_stem, ".txt.gz"))
unmapped_file <- paste0(out_stem, ".unmapped.txt.gz")
a <- as.numeric(system(paste("gunzip -c ", unmapped_file, "| wc -l ",sep=""),intern=TRUE))
if (a > 0) {
    un <- read.table(unmapped_file)
} else { 
    un <- NULL
}

## out is the output we sent in mm9 - now resize to "new" mm10
if(is.null(un)==FALSE)  {
    rates <- out[-match(un[,2],out[,2]),]
}
## rates, the original, is now the same size as new
## remove interval if it's changed size
oldIntSize <- as.integer(rates[,3])-as.integer(rates[,2]    )
newIntSize <- as.integer(new[,3])-as.integer(new[,2]    )
## now work with a new matrix - with new mm10 coordinates, has interval
ratesCur <- new[oldIntSize==newIntSize & new[,1]==paste("chr",chr,sep=""),]
## now order by their location. remove if subsequent intervals overlap
ratesCur <- ratesCur[order(ratesCur[,2]),]
n <- dim(ratesCur)[1]    
x <- ratesCur[-n,3]-ratesCur[-1,2]
ratesCur <- ratesCur[c(TRUE,x<=0),]
                                        # now - we have a matrix, ratesCur
                                        # we have the rates in column 4. strand shouldnt matter
                                        # there may be "gaps" - fill in
n=dim(ratesCur)[1]    
x=ratesCur[-n,3]-ratesCur[-1,2]
whichX=(1:length(x))[x<0]    
                                        # now - x<0 implies there is a jump
                                        # for each of these, get average over 50 kb before, after

##
## where this is a gap in the new matrix, fill in using average 25 + and 25 behind
##
addOn <- cbind(as.character(ratesCur[1,1]),ratesCur[c(x<0,FALSE),3],ratesCur[c(FALSE,x<0),2],0.5,"X","+")
##     
## now fill in using 50000 bp before, after
for(i in 1:dim(addOn)[1]) {
    startOfInterval=ratesCur[whichX[1],2]
    endOfInterval=ratesCur[whichX[1],3]
    ## get rate before
    ## we will do this with a while loop - keep adding until done
    j=whichX[i]-1 # interval before
    rateSum=0
    toAdd=25000 # want over 50,000
    while(toAdd>0) {
        if(j==1) {
            averageBefore=rateSum/(25000-toAdd)
            toAdd <- -1
        } else if (j == 0) {
            averageBefore <- 1
            toAdd <- -1
        } else {
            len=ratesCur[j,3]-ratesCur[j,2]
            if((toAdd-len)>0) # add whole thing
            {
                rateSum=rateSum+len*ratesCur[j,4]
                toAdd=toAdd-len
            } else {
                rateSum=rateSum+(25000-toAdd)*ratesCur[j,4]
                ## add fraction
                toAdd=-1
                averageBefore=rateSum/25000
            }
        }
        j=j-1
    }
    ## get rate after
    j=whichX[i]+1 # interval before
    rateSum=0
    toAdd=25000 # want over 50,000
    while(toAdd>0) {
        len=ratesCur[j,3]-ratesCur[j,2]
        if((toAdd-len)>0) # add whole thing
        {
            rateSum=rateSum+len*ratesCur[j,4]
            toAdd=toAdd-len
        } else {
            rateSum=rateSum+(25000-toAdd)*ratesCur[j,4]
            ## add fraction
            toAdd=-1
            averageAfter=rateSum/25000
        }
        j=j+1
        if(j>dim(ratesCur)[1])
        {
            averageAfter=rateSum/(25000-toAdd)
            toAdd=-1
        }
    }
    val <- 0.5*averageBefore +0.5*averageAfter
    if (is.na(val)) {
        print("had to reset a value")
        val <- 0.1
    }
    addOn[i,4]=val
}

## so we want to add end of x<0, start of next x<0
ratesCur=rbind(ratesCur,addOn)
for(i in 2:4) {
    ratesCur[,i]=as.numeric(as.character(ratesCur[,i]))
}
## now check - there should be no discontinuities
ratesCur=ratesCur[order(ratesCur[,2]),]
n=dim(ratesCur)[1]    
x=ratesCur[-n,3]-ratesCur[-1,2]
print(c(chr,"proper build",table(x)))
## now build normal rate map
n=dim(ratesCur)[1]    
recomb=array(0,c(n+1,3))
colnames(recomb)=c("position","COMBINED_rate.cM.Mb.","Genetic_Map.cM.")

## FIXED bug in v1.1
## rates were off by 1 but cumulative sum were correct
recomb[,1]=c(ratesCur[,2],ratesCur[n,3])
recomb[,2]=c(ratesCur[,4],0)
## now make cumulative sum
recomb[-1,3]= recomb[-n,2] * (recomb[-1,1] - recomb[-n,1])/1000000
recomb[,3]=cumsum(recomb[,3])
##f=function(i,recomb)
##{
##jpeg(paste("/data/outbredmice/imputation/rates/rate.mm",i,".chr",chr,".jpg",sep=""))
##plot(recomb[,1],recomb[,3],type="l")
##dev.off()
##}
## plot and inspect visually
##f("9",ratesOri)
##f("10",recomb)

if (sum(is.na(recomb)) > 0) {
    stop("There are NA in the recomb!")
}

##
## new - smooth!
##
library("STITCH")
smoothed_rate <- c(rcpp_make_smoothed_rate(
    sigma_rate = recomb[, "COMBINED_rate.cM.Mb."],
    L_grid = recomb[, "position"],
    shuffle_bin_radius = 2000
), 0)
recomb[, "COMBINED_rate.cM.Mb."] <- smoothed_rate
fill_in_genetic_map_cm_column <- function(genetic_map) {
    genetic_map[1, "Genetic_Map.cM."] <- 0
    for(iRow in 2:nrow(genetic_map)) {
        distance <-
            genetic_map[iRow, "position"] -
            genetic_map[iRow - 1, "position"]
        rate <- genetic_map[iRow - 1, "COMBINED_rate.cM.Mb."]
        genetic_map[iRow, "Genetic_Map.cM."] <-
            genetic_map[iRow - 1, "Genetic_Map.cM."] + 
            (distance / 1e6) * rate
    }
    return(genetic_map)
}
recomb <- fill_in_genetic_map_cm_column(recomb)

## can now output!
## write to table
recomb[,1] <- as.integer(recomb[,1])
options("scipen"=100, "digits"=4)
out_file <- paste0(recomb_dir, "/", panel, "/", panel, "-chr", chr, "-final.b38.txt")
write.table(
    recomb,
    file = out_file,
    row.names=FALSE,col.names=TRUE,quote=FALSE
)


## gzip
system(paste0("gzip -1 -f ", out_file))
## print output as well
print(paste("chr",chr,"pos range - ori",range(ratesOri[,1])[1],range(ratesOri[,1])[2],"new",range(recomb[,1])[1],range(recomb[,1])[2]))
print(paste("chr",chr,"rate range - ori",range(ratesOri[,2])[1],range(ratesOri[,2])[2],"new",range(recomb[,2])[1],range(recomb[,2])[2]))
print(paste("chr",chr,"cumulative range - ori",range(ratesOri[,3])[1],range(ratesOri[,3])[2],"new",range(recomb[,3])[1],range(recomb[,3])[2]))

unlink(paste0(out_stem, ".txt.gz"))
unlink(paste0(out_stem, ".unmappped.txt.gz"))
unlink(reformated_file)
