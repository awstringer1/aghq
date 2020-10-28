context("Normalization")

test_that("Normalization works",{
    expect_length(norm_sparse_3,3)
    expect_length(norm_sparse_5,3)
    expect_length(norm_sparse_7,3)

    expect_length(norm_trust_3,3)
    expect_length(norm_trust_5,3)
    expect_length(norm_trust_7,3)

    expect_length(norm_trust_3,3)
    expect_length(norm_trust_5,3)
    expect_length(norm_trust_7,3)

    expect_equal(nrow(norm_sparse_3$nodesandweights),3)
    expect_equal(nrow(norm_sparse_5$nodesandweights),5)
    expect_equal(nrow(norm_sparse_7$nodesandweights),7)

    expect_equal(nrow(norm_trust_3$nodesandweights),3)
    expect_equal(nrow(norm_trust_5$nodesandweights),5)
    expect_equal(nrow(norm_trust_7$nodesandweights),7)

    expect_equal(nrow(norm_trust_3$nodesandweights),3)
    expect_equal(nrow(norm_trust_5$nodesandweights),5)
    expect_equal(nrow(norm_trust_7$nodesandweights),7)

    expect_false(is.infinite(norm_sparse_3$lognormconst))
    expect_false(is.infinite(norm_sparse_5$lognormconst))
    expect_false(is.infinite(norm_sparse_7$lognormconst))

    expect_false(is.infinite(norm_trust_3$lognormconst))
    expect_false(is.infinite(norm_trust_5$lognormconst))
    expect_false(is.infinite(norm_trust_7$lognormconst))

    expect_false(is.infinite(norm_trust_3$lognormconst))
    expect_false(is.infinite(norm_trust_5$lognormconst))
    expect_false(is.infinite(norm_trust_7$lognormconst))

    expect_error(normalize_logpost(opt_sparsetrust,2.5))
})