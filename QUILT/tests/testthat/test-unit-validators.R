test_that("can validation niterations and number of block gibbs iterations", {

    expect_null(validate_niterations_and_block_gibbs(
        block_gibbs_iterations = c(3, 6, 9),
        n_gibbs_burn_in_its = 20
    ))

    ## now check errors
    expect_error(validate_niterations_and_block_gibbs(
        block_gibbs_iterations = c(3, 6, 9),
        n_gibbs_burn_in_its = -1
    ))

    expect_error(validate_niterations_and_block_gibbs(
        block_gibbs_iterations = c(3, 6, 9),
        n_gibbs_burn_in_its = "wer"
    ))

    expect_error(validate_niterations_and_block_gibbs(
        block_gibbs_iterations = c(3, 6, 9),
        n_gibbs_burn_in_its = 10.1
    ))

    expect_error(validate_niterations_and_block_gibbs(
        block_gibbs_iterations = c(3, -1, 9),
        n_gibbs_burn_in_its = 20
    ))

    expect_error(validate_niterations_and_block_gibbs(
        block_gibbs_iterations = c(3, "wer", 9),
        n_gibbs_burn_in_its = 20
    ))

    expect_error(validate_niterations_and_block_gibbs(
        block_gibbs_iterations = c(21, 6, 9),
        n_gibbs_burn_in_its = 20
    ))
    
})
