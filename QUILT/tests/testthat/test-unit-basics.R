test_that("dummy unit test", {

    expect_equal(1, 1)
    
})

test_that("dummy cpp unit test", {

    expect_equal(Rcpp_quilt_test_doubler(2), 4)
    
})
