% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quilt-hla-prepare-reference.R
\name{QUILT_HLA_prepare_reference}
\alias{QUILT_HLA_prepare_reference}
\title{QUILT_HLA_prepare_reference}
\usage{
QUILT_HLA_prepare_reference(
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
  all_hla_regions = c("A", "B", "C", "DMA", "DMB", "DOA", "DOB", "DPA1", "DPA2", "DPB1",
    "DPB2", "DQA1", "DQA2", "DQB1", "DRA", "DRB1", "DRB3", "DRB4", "DRB5", "E", "F", "G",
    "HFE", "H", "J", "K", "L", "MICA", "MICB", "N", "P", "S", "TAP1", "TAP2", "T", "U",
    "V", "W", "Y"),
  hla_regions_to_prepare = c("A", "B", "C", "DQA1", "DQB1", "DRB1"),
  chr = "chr6",
  minRate = 0.1,
  nCores = 1
)
}
\arguments{
\item{outputdir}{What output directory to use}

\item{nGen}{Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure, where K is number of haplotypes.}

\item{hla_types_panel}{Path to file with 1000 Genomes formatted HLA types (see example for format details)}

\item{ipd_igmt_alignments_zip_file}{Path to zip file with alignments from IPD-IGMT (see README and example for more details)}

\item{ref_fasta}{Path to reference genome fasta}

\item{refseq_table_file}{Path to file with UCSC refseq gene information (see README and example for more details)}

\item{full_regionStart}{When building HLA full reference panel, start of maximal region spanning all HLA genes. The 1-based position x is kept if regionStart <= x <= regionEnd}

\item{full_regionEnd}{As above, but end of maximal region spanning all HLA genes}

\item{buffer}{Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd}

\item{region_exclude_file}{File with regions to exclude from constructing the reference panel. Particularly useful for QUILT_HLA, where you want to exclude SNPs in the HLA genes themselves, so that reads contribute either to the read mapping or state inference. This file is space separated with a header of Name, Chr, Start and End, with Name being the HLA gene name (e.g. HLA-A), Chr being the chromosome (e.g. chr6), and Start and End are the 1-based starts and ends of the genes (i.e. where we don't want to consider SNPs for the Gibbs sampling state inference)}

\item{genetic_map_file}{Path to file with genetic map information, a file with 3 white-space delimited entries giving position (1-based), genetic rate map in cM/Mbp, and genetic map in cM. If no file included, rate is based on physical distance and expected rate (expRate)}

\item{reference_haplotype_file}{Path to reference haplotype file in IMPUTE format (file with no header and no rownames, one row per SNP, one column per reference haplotype, space separated, values must be 0 or 1)}

\item{reference_legend_file}{Path to reference haplotype legend file in IMPUTE format (file with one row per SNP, and a header including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele)}

\item{reference_sample_file}{Path to reference sample file (file with header, one must be POP, corresponding to populations that can be specified using reference_populations)}

\item{reference_exclude_samplelist_file}{File with one column of samples to exclude from final reference panels. These samples can be removed at one of two points, depending on reference_exclude_samples_for_initial_phasing. If reference_exclude_samples_for_initial_phasing is FALSE, then these samples, if present, are used for the initial phasing and allele assignment. If TRUE, then they are removed immediately}

\item{reference_exclude_samples_for_initial_phasing}{See help for reference_exclude_samplelist_file}

\item{all_hla_regions}{Character vector of all HLA genes for which IPD-IMGT files are available for}

\item{hla_regions_to_prepare}{Character vector of HLA genes to prepare for imputation}

\item{chr}{What chromosome to run (probably chr6 or similar)}

\item{minRate}{Minimum recomb rate cM/Mb}

\item{nCores}{How many cores to use}
}
\value{
Results in properly formatted version
}
\description{
QUILT_HLA_prepare_reference
}
\author{
Robert Davies
}
