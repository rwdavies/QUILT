test_that("can detect failed jobs where processes died", {

    y <- suppressWarnings(parallel::mclapply(1:2, FUN = function(x) x))
    expect_null(check_mclapply_OK(y))
    
    y <- suppressWarnings(parallel::mclapply(1:2, FUN = function(x) if (x == 1) stop("no") else x))
    expect_error(check_mclapply_OK(y))
    
    y <- suppressWarnings(parallel::mclapply(1:2, FUN = function(x) if (x == 1) quit("no") else x))
    expect_error(check_mclapply_OK(y))
    
})

