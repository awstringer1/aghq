context("Normalization")

test_that("Normalization works",{
    # Return correct object
    expect_length(norm_sparse_3,3)
    expect_length(norm_sparse_5,3)
    expect_length(norm_sparse_7,3)

    expect_length(norm_trust_3,3)
    expect_length(norm_trust_5,3)
    expect_length(norm_trust_7,3)

    expect_length(norm_bfgs_3,3)
    expect_length(norm_bfgs_5,3)
    expect_length(norm_bfgs_7,3)

    expect_length(norm_sr1_3,3)
    expect_length(norm_sr1_5,3)
    expect_length(norm_sr1_7,3)

    expect_length(norm_sparse_2d_3,3)
    expect_length(norm_sparse_2d_5,3)
    expect_length(norm_sparse_2d_7,3)

    expect_length(norm_trust_2d_3,3)
    expect_length(norm_trust_2d_5,3)
    expect_length(norm_trust_2d_7,3)

    expect_length(norm_bfgs_2d_3,3)
    expect_length(norm_bfgs_2d_5,3)
    expect_length(norm_bfgs_2d_7,3)

    # Correct number of nodes
    expect_equal(nrow(norm_sparse_3$nodesandweights),3)
    expect_equal(nrow(norm_sparse_5$nodesandweights),5)
    expect_equal(nrow(norm_sparse_7$nodesandweights),7)

    expect_equal(nrow(norm_sr1_3$nodesandweights),3)
    expect_equal(nrow(norm_sr1_5$nodesandweights),5)
    expect_equal(nrow(norm_sr1_7$nodesandweights),7)

    expect_equal(nrow(norm_trust_3$nodesandweights),3)
    expect_equal(nrow(norm_trust_5$nodesandweights),5)
    expect_equal(nrow(norm_trust_7$nodesandweights),7)

    expect_equal(nrow(norm_bfgs_3$nodesandweights),3)
    expect_equal(nrow(norm_bfgs_5$nodesandweights),5)
    expect_equal(nrow(norm_bfgs_7$nodesandweights),7)

    expect_equal(nrow(norm_sparse_2d_3$nodesandweights),3^2)
    expect_equal(nrow(norm_sparse_2d_5$nodesandweights),5^2)
    expect_equal(nrow(norm_sparse_2d_7$nodesandweights),7^2)

    expect_equal(nrow(norm_trust_2d_3$nodesandweights),3^2)
    expect_equal(nrow(norm_trust_2d_5$nodesandweights),5^2)
    expect_equal(nrow(norm_trust_2d_7$nodesandweights),7^2)

    expect_equal(nrow(norm_bfgs_2d_3$nodesandweights),3^2)
    expect_equal(nrow(norm_bfgs_2d_5$nodesandweights),5^2)
    expect_equal(nrow(norm_bfgs_2d_7$nodesandweights),7^2)

    # Normconst is not infinite
    expect_false(is.infinite(norm_sparse_3$lognormconst))
    expect_false(is.infinite(norm_sparse_5$lognormconst))
    expect_false(is.infinite(norm_sparse_7$lognormconst))

    expect_false(is.infinite(norm_sr1_3$lognormconst))
    expect_false(is.infinite(norm_sr1_5$lognormconst))
    expect_false(is.infinite(norm_sr1_7$lognormconst))

    expect_false(is.infinite(norm_trust_3$lognormconst))
    expect_false(is.infinite(norm_trust_5$lognormconst))
    expect_false(is.infinite(norm_trust_7$lognormconst))

    expect_false(is.infinite(norm_bfgs_3$lognormconst))
    expect_false(is.infinite(norm_bfgs_5$lognormconst))
    expect_false(is.infinite(norm_bfgs_7$lognormconst))

    expect_false(is.infinite(norm_sparse_2d_3$lognormconst))
    expect_false(is.infinite(norm_sparse_2d_5$lognormconst))
    expect_false(is.infinite(norm_sparse_2d_7$lognormconst))

    expect_false(is.infinite(norm_trust_2d_3$lognormconst))
    expect_false(is.infinite(norm_trust_2d_5$lognormconst))
    expect_false(is.infinite(norm_trust_2d_7$lognormconst))

    expect_false(is.infinite(norm_bfgs_2d_3$lognormconst))
    expect_false(is.infinite(norm_bfgs_2d_5$lognormconst))
    expect_false(is.infinite(norm_bfgs_2d_7$lognormconst))



    # Require integer number of quad points
    expect_error(normalize_logpost(opt_sparsetrust,2.5))

    # Normconst actually equals the correct thing
    expect_equal(round(norm_sparse_3$lognormconst,2),round(truelognormconst,2))
    expect_equal(round(norm_sparse_5$lognormconst,3),round(truelognormconst,3))
    expect_equal(round(norm_sparse_7$lognormconst,3),round(truelognormconst,3))

    expect_equal(round(norm_sr1_3$lognormconst,2),round(truelognormconst,2))
    expect_equal(round(norm_sr1_5$lognormconst,3),round(truelognormconst,3))
    expect_equal(round(norm_sr1_7$lognormconst,3),round(truelognormconst,3))

    expect_equal(round(norm_trust_3$lognormconst,2),round(truelognormconst,2))
    expect_equal(round(norm_trust_5$lognormconst,3),round(truelognormconst,3))
    expect_equal(round(norm_trust_7$lognormconst,3),round(truelognormconst,3))

    expect_equal(round(norm_bfgs_3$lognormconst,2),round(truelognormconst,2))
    expect_equal(round(norm_bfgs_5$lognormconst,3),round(truelognormconst,3))
    expect_equal(round(norm_bfgs_7$lognormconst,3),round(truelognormconst,3))

    expect_equal(round(norm_sparse_2d_3$lognormconst,1),round(truelognormconst2d,1))
    expect_equal(round(norm_sparse_2d_5$lognormconst,2),round(truelognormconst2d,2))
    expect_equal(round(norm_sparse_2d_7$lognormconst,2),round(truelognormconst2d,2))

    expect_equal(round(norm_trust_2d_3$lognormconst,1),round(truelognormconst2d,1))
    expect_equal(round(norm_trust_2d_5$lognormconst,2),round(truelognormconst2d,2))
    expect_equal(round(norm_trust_2d_7$lognormconst,2),round(truelognormconst2d,2))

    expect_equal(round(norm_bfgs_2d_3$lognormconst,1),round(truelognormconst2d,1))
    expect_equal(round(norm_bfgs_2d_5$lognormconst,2),round(truelognormconst2d,2))
    expect_equal(round(norm_bfgs_2d_7$lognormconst,2),round(truelognormconst2d,2))

    # Parameter vector reordering works
    # Note: some of these tests only work because the parameters are independent in this example!
    expect_equal(unique(round(norm_sparse_2d_3$nodesandweights$theta1,3)),unique(round(norm_sparse_2d_reorder_3$nodesandweights$theta1,3)))
    expect_equal(unique(round(norm_sparse_2d_3$nodesandweights$theta2,3)),unique(round(norm_sparse_2d_reorder_3$nodesandweights$theta2,3)))

    expect_equal(sort(norm_sparse_2d_3$nodesandweights$logpost),sort(norm_sparse_2d_reorder_3$nodesandweights$logpost))
    expect_equal(sort(norm_sparse_2d_3$nodesandweights$logpost_normalized),sort(norm_sparse_2d_reorder_3$nodesandweights$logpost_normalized))

    expect_equal(norm_sparse_2d_3$grid$features$m,norm_sparse_2d_reorder_3$grid$features$m[c(2,1)])
    expect_equal(norm_sparse_2d_3$grid$features$C,norm_sparse_2d_reorder_3$grid$features$C[c(2,1),c(2,1)])

    expect_equal(norm_sparse_2d_3$lognormconst,norm_sparse_2d_reorder_3$lognormconst)

})