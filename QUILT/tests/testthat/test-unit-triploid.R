test_that("triploid mode uses equal three-haplotype priors", {
    expect_equal(
        QUILT:::get_three_haplotype_prior_probs(ff = 0, sample_is_triploid = TRUE),
        rep(1 / 3, 3)
    )
    expect_equal(
        QUILT:::get_three_haplotype_prior_probs(ff = 0.8, sample_is_triploid = TRUE),
        rep(1 / 3, 3)
    )
    expect_equal(
        QUILT:::get_three_haplotype_prior_probs(ff = 0.2, sample_is_triploid = FALSE),
        c(0.5, 0.4, 0.1)
    )
})

test_that("triploid read label initialization uses equal haplotype weights", {
    nReads <- 30000
    set.seed(1)
    H_triploid <- QUILT:::random_gibbs_nipt_read_labels(
        nReads = nReads,
        ff = 0,
        sample_is_triploid = TRUE
    )
    observed_triploid <- as.numeric(prop.table(table(factor(H_triploid, levels = 1:3))))
    expect_true(max(abs(observed_triploid - 1 / 3)) < 0.02)

    set.seed(1)
    H_nipt <- QUILT:::random_gibbs_nipt_read_labels(
        nReads = nReads,
        ff = 0,
        sample_is_triploid = FALSE
    )
    observed_nipt <- as.numeric(prop.table(table(factor(H_nipt, levels = 1:3))))
    expect_true(max(abs(observed_nipt - c(0.5, 0.5, 0))) < 0.02)
})

test_that("triploid genotype probabilities combine three haploid probabilities", {
    hap_probs_t <- rbind(
        c(0.1, 0.5),
        c(0.2, 0.5),
        c(0.3, 0.5)
    )
    gen_probs_t <- QUILT:::triploid_genotype_probs_from_hap_probs(hap_probs_t)
    expect_equal(colSums(gen_probs_t), c(1, 1))
    expect_equal(gen_probs_t[, 2], c(0.125, 0.375, 0.375, 0.125))
    expect_equal(
        gen_probs_t[, 1],
        c(0.504, 0.398, 0.092, 0.006)
    )
})

test_that("triploid genotype probabilities preserve dosage expectations", {
    hap_probs_t <- rbind(
        c(0, 0.1, 0.4, 1),
        c(0.2, 0.3, 0.5, 1),
        c(0.7, 0.8, 0.6, 1)
    )
    gen_probs_t <- QUILT:::triploid_genotype_probs_from_hap_probs(hap_probs_t)
    dosage <- colSums(gen_probs_t * 0:3)

    expect_equal(colSums(gen_probs_t), rep(1, ncol(hap_probs_t)))
    expect_equal(dosage, colSums(hap_probs_t))
})

test_that("triploid likelihood calculations use supplied equal priors", {
    c1 <- c(0.5, 0.25, 0.125)
    c2 <- c(0.2, 0.4, 0.8)
    c3 <- c(0.1, 0.8, 0.5)
    H <- as.integer(c(1, 1, 2, 3, 3, 3))
    prior_probs <- rep(1 / 3, 3)

    out <- calculate_likelihoods_values(
        c1 = c1,
        c2 = c2,
        c3 = c3,
        H = H,
        nGrids = length(c1),
        prior_probs = prior_probs,
        ff = 0
    )
    read_counts <- as.numeric(table(factor(H, levels = 1:3)))

    expect_equal(out[1], -sum(log(c1)))
    expect_equal(out[2], -sum(log(c2)))
    expect_equal(out[3], -sum(log(c3)))
    expect_equal(out[5], sum(log(prior_probs[H])))
    expect_equal(out[7], as.numeric(dmultinom(read_counts, prob = prior_probs, log = TRUE)))
})

