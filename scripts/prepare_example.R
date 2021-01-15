##
## not for general use
## use to prepare example package to demonstrate usage
## to be run from rescomp in single_imp results directory on well
## e.g. /well/davies/users/user-name/single_imp/ && R -f ~/proj/QUILT/scripts/prepare_example.R
## 

library("QUILT")
setwd("/well/davies/users/dcc832/single_imp/")

source("~/proj/QUILT/scripts/prepare_example_functions.R")

## smaller package for quick test 
prepare_example_using_1000G(
    package_output_date = "2021_01_15A",
    chr = "chr20",
    regionStart = 2000001,
    regionEnd = 2100000,
    buffer = 10000
)

## larger package for example
prepare_example_using_1000G(
    package_output_date = "2021_01_15B",
    chr = "chr20",
    regionStart = 2000001,
    regionEnd = 4000000,
    buffer = 500000
)


1 




##prepare_example_using_hrc(
##    package_output_date = "2020_08_25",
##    main_run_date = "2020_08_25",
##    chr = "chr20",
##    regionStart = 2000001,
##    regionEnd = 4000000,
##    buffer = 500000
## )

