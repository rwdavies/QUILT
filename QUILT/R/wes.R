#' @export
run_haplotype_selection <- function(
    WES,
    rhb_t,
    pos,
    ref_alleleCount,
    nSNPs,
    depth_threshold,
    region_divide,
    toolForSelection
) {
  
    WES <- refine_WES_data(WES,depth_threshold)
    pos.MAF <- merge_pos_MAF(pos,ref_alleleCount)

    if (toolForSelection == 'AF') {
        haps_for_subset <- select_K_haps_by_rare_alleles(
            WES, 
            rhb_t, 
            pos.MAF,
            pos, 
            nSNPs,
            depth_threshold,
            region_divide
        )
    } else {
      haps_for_subset <- select_K_haps_by_match_length(
          WES, 
          rhb_t, 
          pos, 
          nSNPs,
          depth_threshold,
          region_divide
      )
  }
    
    return(haps_for_subset)
    
}





Extract.Genotypes <- function(vcf) {
  vcf.separated <- separate(data = vcf, col = V10, into = c("GT","PL","DP"), sep="\\:")
  genotyped <- separate(data = vcf.separated, col = GT, into = c("H1", "H2"), sep="\\/")
  genotyped$GEN <- as.numeric(genotyped$H1) + as.numeric(genotyped$H2)
  genotyped$pos.chr.ref.alt <- paste(genotyped$V2, ":", genotyped$V1,":", genotyped$V4,":", genotyped$V5 )
  
  genotyped <- separate(data = genotyped, col = V8, into = c('DP','MQ0F','AC','AN','DP4','MQ'), sep="\\;")
  genotyped <- separate(data = genotyped, col = DP, into = c('Key','Depth'), sep="\\=")
  
  
  
  pos.gen.depth <- data.frame(genotyped$pos.chr.ref.alt, genotyped$GEN, genotyped$Depth)
  colnames(pos.gen.depth) <- c('pos.chr.ref.alt', 'GEN', 'Depth')
  
  na <- which(is.na(pos.gen.depth$GEN))
  pos.gen.depth.na.removed <- pos.gen.depth[-na,]
  
  return(pos.gen.depth.na.removed)
}

# Filter the genotyped regions to only sites with sufficient depth of reads.
# hap_data to input is in the form gen.depth

filter_to_exomic_regions <- function(hap_data, depth.threshold){
  filtered <- hap_data[as.numeric(hap_data$Depth)>depth.threshold,]
  return(filtered)
}


# Input the indices of the exonic sites 
find_exomic_regions_indices <- function(exonic_sites){
  
  no_exon_sites <- length(exonic_sites)
  differences   <- exonic_sites[1:no_exon_sites]-c(0,exonic_sites[1:(no_exon_sites-1)])
  exonic_breaks  <- which(differences>5)
  return(exonic_breaks)
  
}

inflate_haps_to_check <- function(haps_to_check, rhb_t, nSNPs) {
  #print(nSNPs)
  rhb <- t(rhb_t)
  inflated_haps <- inflate_fhb(rhb, haps_to_check, nSNPs)
  colnames(inflated_haps) <- haps_to_check
  return(inflated_haps)
}

get_genotyped_with_refs <- function(genotyped, pos.refs) {
  gen.refs <- merge(genotyped, pos.refs, by='pos.chr.ref.alt') 
  return(gen.refs)
}


check_if_match <- function(WES, ref_hap) {
  
  
  if (WES==0 && ref_hap==0)  {
    match <- 1
  }
  else if (WES==1 && (ref_hap==1 || ref_hap==0)){
    match <- 1
  }
  else if (WES==2 && ref_hap==1) {
    match <- 1
  }
  else {
    match <- 0
  }
  
  return(match)
}

##
find_het_SNPs_in_sample <- function(sample_gen_data) {
  het.pos <- sample_gen_data[which(sample_gen_data$GEN==1),c('pos.chr.ref.alt','MAF')]
  colnames(het.pos) <- c('pos.chr.ref.alt', 'MAF')
  rare_alleles_in_sample <- het.pos[order(het.pos$MAF),]
  return(rare_alleles_in_sample)
}
##


find_het_sites <- function(WES_region){
  het.sites <- which(WES_region$GEN==1)
  return(het.sites)
}

## het.sites <- find_het_sites(WES)



find_rarest_SNP_in_region<- function(WES_region, allele_freqs_in_region) {
  het.sites <- find_het_sites(WES_region)
  print('HET SITES')
  print(het.sites)
  if (length(het.sites)==0) {
    print('NO HET SITES')
    rarest_SNP_in_region <- NA
    return(list(NA,NA))
  }
  else {
    print(het.sites)
    het_allele_freqs <- allele_freqs_in_region[het.sites, c('pos.chr.ref.alt','MAF')]
    #print(het_allele_freqs)
    rarest_SNP_in_region <- which.min(het_allele_freqs$MAF)
    #print(rarest_SNP_in_region)
    #print(het_allele_freqs$MAF[rarest_SNP_in_region])
    return(list(rarest_SNP_in_region,het_allele_freqs$MAF[rarest_SNP_in_region]))
  }
}

find_het_sites_per_region <- function(WES, region_divide) {
  no_regions <- length(region_divide)-1
  het.sites.per.region <- vector(mode='list', length=no_regions)
  
  for (j in (1:(no_regions))) {
    REGION_START <- region_divide[j]
    REGION_END   <- region_divide[j+1]-1
    het.sites <- find_het_sites(WES[REGION_START:REGION_END,])
    het.sites.per.region[[j]] <- het.sites
  }
  return(het.sites.per.region)
}

## find_het_sites_per_region(WES.refined, region_divide)

find_rarest_SNP_per_region <- function(WES, allele_freqs, region_divide) {
  no_regions <- length(region_divide)-1
  rarest.SNP.index.per.region <- data.frame(REGION=c(1:no_regions), SNP.index=rep(NA,no_regions))
  for (j in (1:(no_regions))) {
    REGION_START <- region_divide[j]
    REGION_END   <- region_divide[j+1]-1
    index.rarest.SNP <- find_rarest_SNP_in_region(WES[REGION_START:REGION_END,],
                                                  allele_freqs[REGION_START:REGION_END,])
    rarest.SNP.index.per.region$SNP.index[j] <- index.rarest.SNP
  }
  return(rarest.SNP.index.per.region)
}


find_rarest_SNP_per_region <- function(WES, allele_freqs, region_divide) {
  no_regions <- length(region_divide)-1
  rarest.SNP.index.per.region <- data.frame(REGION=c(1:no_regions), SNP.index=rep(NA,no_regions), MAF=rep(NA,no_regions))
  for (j in (1:(no_regions))) {
    REGION_START <- region_divide[j]
    REGION_END   <- region_divide[j+1]-1
    print('ALLELE FREQS IN REGION')
    print(allele_freqs[REGION_START:REGION_END,])
    index.rarest.SNP <- find_rarest_SNP_in_region(WES[REGION_START:REGION_END,],
                                                  allele_freqs[REGION_START:REGION_END,])[[1]]
    MAF.rarest.SNP <- find_rarest_SNP_in_region(WES[REGION_START:REGION_END,],
                                                allele_freqs[REGION_START:REGION_END,])[[2]]
    
    rarest.SNP.index.per.region$SNP.index[j] <- index.rarest.SNP
    rarest.SNP.index.per.region$MAF[j] <- MAF.rarest.SNP
  }
  return(rarest.SNP.index.per.region)
}

check_if_het_sites <- function(region) {
  if (length(region)==0) {
    return(NA)
  }
  else {
    return(1)
  }
}

find_hap_not_added <- function(haps_to_consider, haps_added){
  UNIQUE_HAP_FOUND <- FALSE 
  index <- 1
  while (UNIQUE_HAP_FOUND==FALSE) {
    hap <- haps_to_consider[[index]]
    
    if (!(hap %in% haps_added)) {
      UNIQUE_HAP <- hap
      UNIQUE_HAP_FOUND = TRUE
    }
    
    else {
      index <- index + 1
    }
    
  }
  return(UNIQUE_HAP)
}

order_by_match <- function(df) {
  df <- df[order(df$MATCH_LENGTH, decreasing=TRUE),] 
  return(df)
}

refine_WES_data <- function(WES, depth_threshold){
  pos.gen.depths <- Extract.Genotypes(WES)
  pos.gen.depths.exon <- filter_to_exomic_regions(pos.gen.depths, depth_threshold)
  
}



merge_pos_MAF <- function(pos, ref_alleleCount) {
  
  pos$pos.chr.ref.alt <- paste(pos$POS, ":", pos$CHR,":", pos$REF,":", pos$ALT )
  ref_alleleCount <- as.data.frame(ref_alleleCount)
  ref_alleleCount$MAF <- as.numeric(ref_alleleCount[,3])
  w <- (ref_alleleCount[, 3] > 0.5)
  ref_alleleCount[w, "MAF"] <- 1 - ref_alleleCount[w, "MAF"]
  
  pos.MAF <- as.data.frame(cbind(pos$pos.chr.ref.alt,ref_alleleCount$MAF))
  colnames(pos.MAF) <- c('pos.chr.ref.alt', 'MAF')
  return(pos.MAF)
  
}



select_K_haps_by_match_length <- function(pos.gen.depths.exon, 
                                          rhb_t, 
                                          pos, 
                                          nSNPs,
                                          depth_threshold,
                                          region_divide){
  
  no_haps_in_reference <- nrow(rhb_t)
  no_regions <- length(region_divide)-1
  pos$pos.chr.ref.alt <- paste(pos$POS, ":", pos$CHR,":", pos$REF,":", pos$ALT )
  
  pos_gen_data <- pos.gen.depths.exon[,1]
  ref_positions_to_keep <- which(pos$pos.chr.ref.alt %in% pos.gen.depths.exon$pos.chr.ref.alt)
  pos <-pos[ref_positions_to_keep,]
  WES     <-pos.gen.depths.exon[which(pos.gen.depths.exon$pos.chr.ref.alt %in% pos$pos.chr.ref.alt),]

  
  region.list <- vector(mode='list', length = no_regions)
  for (i in (1:no_regions)) {
    region.list[[i]] <- data.frame(REF_HAP=c(1:no_haps_in_reference),MATCH_LENGTH=rep(NA,no_haps_in_reference))
  }
  HAPS_FOR_SUBSET <- rep(NA,no_regions)
  
  for (i in 0:(no_haps_in_reference-1)) {
    hap_to_evaluate <- inflate_haps_to_check(i, rhb_t, nSNPs)
    hap_refined <- hap_to_evaluate[ref_positions_to_keep]
    #print(length(hap_refined))
    #print(length(WES[,1]))
    #colnames(pos.ref)[1] <-'pos.chr.ref.alt'
    # Merges genotyped data with reference hap by position
    #pos.WES.hap <- get_genotyped_with_refs(pos.gen.depths.exon[,1:2], pos.ref)
    
    for (j in (1:(no_regions))) {
      REGION_START <- region_divide[j]
      REGION_END   <- region_divide[j+1]-1
      #MATCHES_THIS_REGION <- mapply(check_if_match,pos.WES.hap[REGION_START:REGION_END,2],
      #                              pos.WES.hap[REGION_START:REGION_END,3])
      MATCHES_THIS_REGION <- mapply(check_if_match,WES[REGION_START:REGION_END,2],
                                    hap_refined[REGION_START:REGION_END])
      MATCH_LENGTHS <- ave(MATCHES_THIS_REGION, rleid(MATCHES_THIS_REGION), FUN = seq_along)*MATCHES_THIS_REGION
      MAX_MATCH <- which.max(MATCH_LENGTHS)
      region.list[[j]][i,2] <- MAX_MATCH
      
    }
  }
  ordered_region_list <- lapply(region.list,order_by_match)
  for (i in (1:no_regions)) {
    potential_haps_for_region <- ordered_region_list[[i]]$REF_HAP
    HAPS_FOR_SUBSET[i] <- find_hap_not_added(potential_haps_for_region,HAPS_FOR_SUBSET)
    #HAPS_FOR_SUBSET[i] <- ordered_region_list[[i]][1,1]
  } 
  return(ordered_region_list[[4]])
}



select_K_haps_by_rare_alleles <- function(pos.gen.depths.exon, 
                                          rhb_t, 
                                          pos.MAF,
                                          pos, 
                                          nSNPs,
                                          depth_threshold,
                                          region_divide){
  
  no_haps_in_reference <- nrow(rhb_t)
  no_regions <- length(region_divide)-1
  
  ref_positions_to_keep <- which(pos.MAF$pos.chr.ref.alt %in% pos.gen.depths.exon$pos.chr.ref.alt)
  pos.MAF <-pos.MAF[ref_positions_to_keep,]
  WES     <-pos.gen.depths.exon[which(pos.gen.depths.exon$pos.chr.ref.alt %in% pos.MAF$pos.chr.ref.alt),]
  
  #print(length(pos.MAF[,1]))
  ##print(length(WES[,1]))

  
  het.sites.per.region <- find_het_sites_per_region(WES, region_divide)
  
  regions_with_het_sites <- lapply(het.sites.per.region, check_if_het_sites)
  no_het_sites <- which(is.na(regions_with_het_sites))
  
  #print(het.sites.per.region)
  rare_SNP_indices <- find_rarest_SNP_per_region(WES, pos.MAF, region_divide)
  print(rare_SNP_indices)
  hap.per.region <- data.frame(REGION=c(1:no_regions), HAP_SELECTED=rep(NA,no_regions), MAF=rep(NA,no_regions))
  hap.per.region$HAP_SELECTED[no_het_sites] <- 'NO HET SITES'
  
  
  for (i in 0:(no_haps_in_reference-1)) {
      if (length(which(is.na(hap.per.region$HAP_SELECTED)))==0) {
        print('HAPS ALL SELECTED')
        feasible_haps <- as.numeric(hap.per.region$HAP_SELECTED)[-which(is.na(as.numeric(hap.per.region$HAP_SELECTED)))]
        return(feasible_haps)
      }
      
      hap_to_evaluate <- inflate_haps_to_check(i, rhb_t, nSNPs)
      hap_refined <- hap_to_evaluate[ref_positions_to_keep]
      regions_to_consider <- which(is.na(hap.per.region$HAP_SELECTED))
      for (region in (regions_to_consider)) {
          if (length(regions_to_consider)>0) {
            REGION_START <- region_divide[region]
            REGION_END   <- region_divide[region+1]-1
            hap_region_to_consider <- hap_refined[REGION_START:REGION_END]
            print('HET SITES PER REGION')
            print(het.sites.per.region[[region]])
            hap_refined_to_het <- hap_region_to_consider[het.sites.per.region[[region]]]
            if (is.na(hap_refined_to_het[rare_SNP_indices$SNP.index[region]])){
              print('NO HET SITES IN REGION')
            }
            
            else if (hap_refined_to_het[rare_SNP_indices$SNP.index[region]]==1) {
              if (i %in% hap.per.region$HAP_SELECTED) {
                print('HAP ALREADY SELECTED')
              }
              else {
             hap.per.region$HAP_SELECTED[region] <- i
             hap.per.region$MAF[region] <- rare_SNP_indices$MAF[region]
              }
           }
        }
    }
  }
  print(hap.per.region)
  feasible_haps <- as.numeric(hap.per.region$HAP_SELECTED)[-which(is.na(as.numeric(hap.per.region$HAP_SELECTED)))]
  
  return(feasible_haps)
}

