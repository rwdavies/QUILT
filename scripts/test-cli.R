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
make_STITCH_cli(
    function_file = "QUILT/R/quilt.R",
    cli_output_file = cli_output_file,
    integer_vectors = c("unwindIterations", "sample_unwindIterations"),
    other_character_params = c("phasefile", "unwindIterations", "sample_unwindIterations"),
    other_logical_params = c("outputHaplotypeProbabilities", "very_verbose"),
    function_name = "QUILT",
    library_name = "QUILT"
)
system(paste0("chmod +x ", cli_output_file))

message("test that QUILT CLI produces help message")
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



##
## todo, add more CLI tests
##
