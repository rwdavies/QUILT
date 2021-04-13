#!/usr/bin/env Rscript

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))

library("testthat")
library("parallel")
library("STITCH") ## for make_STITCH_cli

## testthat doesn't do what I want outside of package form
## so don't bother wrappping, just fail

cli_function_build <- Sys.getenv("CLI_FUNCTION_BUILD")
if (cli_function_build != "") {
    print(paste0("Using ", cli_function_build))
    dir <- tempdir()
    system(paste0("cp ", cli_function_build, " ", dir, "/"))
    system(paste0("(cd ", dir, " && tar -zxvf ", dir, "/*tar.gz QUILT/R/functions.R)"))
    function_file <- file.path(dir, "STITCH/R/functions.R")
} else {
    function_file <- "QUILT/R/quilt.R"
}

cli_output_file <- "QUILT.R"
STITCH::make_STITCH_cli(
    function_file = "QUILT/R/quilt.R",
    cli_output_file = cli_output_file,
    other_character_params = c("phasefile", "output_RData_filename", "RData_objects_to_save", "prepared_reference_filename", "reference_exclude_samplelist_file", "output_sites_filename", "region_exclude_file"),
    integer_vectors = c("block_gibbs_iterations"),
    character_vectors = c("RData_objects_to_save"),
    other_logical_params = c("make_plots", "verbose", "record_read_label_usage", "record_interim_dosages", "use_bx_tag", "addOptimalHapsToVCF", "estimate_bq_using_truth_read_labels", "make_plots_block_gibbs", "override_default_params_for_small_ref_panel", "hla_run", "make_fake_vcf_with_sites_list", "print_extra_timing_information", "plot_per_sample_likelihoods"),
    other_integer_params = c("nGibbsSamples", "n_seek_its", "Ksubset", "Knew", "K_top_matches", "panel_size", "bxTagUpperLimit", "seed", "gamma_physically_closest_to", "nMaxDH", "minimum_number_of_sample_reads", "block_gibbs_iterations", "n_gibbs_burn_in_its"),
    other_double_params = c("heuristic_match_thin", "minGLValue"),
    function_name = "QUILT",
    library_name = "QUILT"
)
system(paste0("chmod +x ", cli_output_file))


cli_output_file <- "QUILT_prepare_reference.R"
STITCH::make_STITCH_cli(
    function_file = "QUILT/R/quilt-prepare-reference.R",
    cli_output_file = cli_output_file,
    other_character_params = c("output_file", "reference_exclude_samplelist_file", "output_sites_filename", "region_exclude_file"),
    other_logical_params = c("make_fake_vcf_with_sites_list"),
    other_integer_params = c("nMaxDH"),
    function_name = "QUILT_prepare_reference",
    library_name = "QUILT"
)
system(paste0("chmod +x ", cli_output_file))



cli_output_file <- "QUILT_HLA.R"
STITCH::make_STITCH_cli(
    function_file = "QUILT/R/quilt-hla.R",
    cli_output_file = cli_output_file,
    other_character_params = c("bamlist", "region", "finaloutputfile", "chr", "quilt_hla_haplotype_panelfile", "prepared_hla_reference_dir", "summary_output_file_prefix", "final_output_RData_file", "hla_gene_region_file", "dict_file"),
    other_logical_params = c("overrideoutput", "write_summary_text_files"),
    other_integer_params = c("nGibbsSamples", "n_seek_iterations", "quilt_seed", "quilt_buffer", "quilt_bqFilter"),
    other_double_params = c("summary_best_alleles_threshold"),
    function_name = "QUILT_HLA",
    library_name = "QUILT"
)
system(paste0("chmod +x ", cli_output_file))

cli_output_file <- "QUILT_HLA_prepare_reference.R"
STITCH::make_STITCH_cli(
    function_file = "QUILT/R/quilt-hla-prepare-reference.R",
    cli_output_file = cli_output_file,
    other_character_params = c("outputdir", "ipd_igmt_alignments_zip_file", "quilt_hla_supplementary_info_file", "all_hla_regions", "hla_regions_to_prepare", "reference_exclude_samplelist_file", "region_exclude_file", "hla_gene_region_file", "hla_types_panel"),
    other_double_params = c("full_regionStart", "full_regionEnd", "full_buffer"),
    other_logical_params = c("reference_exclude_samples_for_initial_phasing"),
    character_vectors = c("all_hla_regions", "hla_regions_to_prepare"),
    function_name = "QUILT_HLA_prepare_reference",
    library_name = "QUILT"
)

system(paste0("chmod +x ", cli_output_file))



for(cli_output_file in c("QUILT.R", "QUILT_prepare_reference.R", "QUILT_HLA.R", "QUILT_HLA_prepare_reference.R")) {
    message(paste0("test that ", cli_output_file, " CLI produces help message"))
    ## behaviour of optparse changed!
    ## now exits code 0 as one would hope on --help
    stdout_file <- tempfile()
    stderr_file <- tempfile()
    out <- system2(
        cli_output_file,
        args = c(
            "--help"
        ),
        stdout = stdout_file, stderr = stderr_file
    )
    stderr <- system(paste0("cat ", stderr_file), intern = TRUE)
    stdout <- system(paste0("cat ", stdout_file), intern = TRUE)
    if (out > 0) {
        message("---stderr---")
        print(stderr)
        message("---stdout---")
        print(stdout)
    }
    expect_equal(0, out)
}



##
## todo, add more CLI tests
##
