test_that("can validation niterations and number of small_ref_panel_block gibbs iterations", {

    expect_null(validate_niterations_and_small_ref_panel_block_gibbs(
        small_ref_panel_block_gibbs_iterations = c(3, 6, 9),
        small_ref_panel_gibbs_iterations = 20
    ))

    ## now check errors
    expect_error(validate_niterations_and_small_ref_panel_block_gibbs(
        small_ref_panel_block_gibbs_iterations = c(3, 6, 9),
        small_ref_panel_gibbs_iterations = -1
    ))

    expect_error(validate_niterations_and_small_ref_panel_block_gibbs(
        small_ref_panel_block_gibbs_iterations = c(3, 6, 9),
        small_ref_panel_gibbs_iterations = "wer"
    ))

    expect_error(validate_niterations_and_small_ref_panel_block_gibbs(
        small_ref_panel_block_gibbs_iterations = c(3, 6, 9),
        small_ref_panel_gibbs_iterations = 10.1
    ))

    expect_error(validate_niterations_and_small_ref_panel_block_gibbs(
        small_ref_panel_block_gibbs_iterations = c(3, -1, 9),
        small_ref_panel_gibbs_iterations = 20
    ))

    expect_error(validate_niterations_and_small_ref_panel_block_gibbs(
        small_ref_panel_block_gibbs_iterations = c(3, "wer", 9),
        small_ref_panel_gibbs_iterations = 20
    ))

    expect_error(validate_niterations_and_small_ref_panel_block_gibbs(
        small_ref_panel_block_gibbs_iterations = c(21, 6, 9),
        small_ref_panel_gibbs_iterations = 20
    ))
    
})


test_that("can validate n_seek_its and n_burn_in_seek_its", {


    expect_null(validate_n_seek_its_and_n_burn_in_seek_its(
        n_seek_its = 3,
        n_burn_in_seek_its = 2
    ))

    expect_null(validate_n_seek_its_and_n_burn_in_seek_its(
        n_seek_its = 4,
        n_burn_in_seek_its = 2
    ))

    expect_null(validate_n_seek_its_and_n_burn_in_seek_its(
        n_seek_its = 4,
        n_burn_in_seek_its = 0
    ))

    expect_error(validate_n_seek_its_and_n_burn_in_seek_its(
        n_seek_its = 4,
        n_burn_in_seek_its = 4
    ))

    expect_error(validate_n_seek_its_and_n_burn_in_seek_its(
        n_seek_its = 4,
        n_burn_in_seek_its = 5
    ))

    expect_error(validate_n_seek_its_and_n_burn_in_seek_its(
        n_seek_its = 4,
        n_burn_in_seek_its = -1
    ))
    
    expect_error(validate_n_seek_its_and_n_burn_in_seek_its(
        n_seek_its = 4,
        n_burn_in_seek_its = 2.2
    ))
    
    
    
})
