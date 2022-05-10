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

    expect_length(norm_sr1_2d_3,3)
    expect_length(norm_sr1_2d_5,3)
    expect_length(norm_sr1_2d_7,3)

    expect_length(norm_sparse_3d_3,3)
    expect_length(norm_sparse_3d_5,3)
    expect_length(norm_sparse_3d_7,3)

    expect_length(norm_trust_3d_3,3)
    expect_length(norm_trust_3d_5,3)
    expect_length(norm_trust_3d_7,3)

    expect_length(norm_bfgs_3d_3,3)
    expect_length(norm_bfgs_3d_5,3)
    expect_length(norm_bfgs_3d_7,3)

    expect_length(norm_sr1_3d_3,3)
    expect_length(norm_sr1_3d_5,3)
    expect_length(norm_sr1_3d_7,3)

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

    expect_equal(nrow(norm_sr1_2d_3$nodesandweights),3^2)
    expect_equal(nrow(norm_sr1_2d_5$nodesandweights),5^2)
    expect_equal(nrow(norm_sr1_2d_7$nodesandweights),7^2)

    expect_equal(nrow(norm_sparse_3d_3$nodesandweights),3^3)
    expect_equal(nrow(norm_sparse_3d_5$nodesandweights),5^3)
    expect_equal(nrow(norm_sparse_3d_7$nodesandweights),7^3)

    expect_equal(nrow(norm_trust_3d_3$nodesandweights),3^3)
    expect_equal(nrow(norm_trust_3d_5$nodesandweights),5^3)
    expect_equal(nrow(norm_trust_3d_7$nodesandweights),7^3)

    expect_equal(nrow(norm_bfgs_3d_3$nodesandweights),3^3)
    expect_equal(nrow(norm_bfgs_3d_5$nodesandweights),5^3)
    expect_equal(nrow(norm_bfgs_3d_7$nodesandweights),7^3)

    expect_equal(nrow(norm_sr1_3d_3$nodesandweights),3^3)
    expect_equal(nrow(norm_sr1_3d_5$nodesandweights),5^3)
    expect_equal(nrow(norm_sr1_3d_7$nodesandweights),7^3)

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

    expect_false(is.infinite(norm_sr1_2d_3$lognormconst))
    expect_false(is.infinite(norm_sr1_2d_5$lognormconst))
    expect_false(is.infinite(norm_sr1_2d_7$lognormconst))

    expect_false(is.infinite(norm_sparse_3d_3$lognormconst))
    expect_false(is.infinite(norm_sparse_3d_5$lognormconst))
    expect_false(is.infinite(norm_sparse_3d_7$lognormconst))

    expect_false(is.infinite(norm_trust_3d_3$lognormconst))
    expect_false(is.infinite(norm_trust_3d_5$lognormconst))
    expect_false(is.infinite(norm_trust_3d_7$lognormconst))

    expect_false(is.infinite(norm_bfgs_3d_3$lognormconst))
    expect_false(is.infinite(norm_bfgs_3d_5$lognormconst))
    expect_false(is.infinite(norm_bfgs_3d_7$lognormconst))

    expect_false(is.infinite(norm_sr1_3d_3$lognormconst))
    expect_false(is.infinite(norm_sr1_3d_5$lognormconst))
    expect_false(is.infinite(norm_sr1_3d_7$lognormconst))

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

    expect_equal(round(norm_sr1_2d_3$lognormconst,1),round(truelognormconst2d,1))
    expect_equal(round(norm_sr1_2d_5$lognormconst,2),round(truelognormconst2d,2))
    expect_equal(round(norm_sr1_2d_7$lognormconst,2),round(truelognormconst2d,2))

    expect_equal(round(norm_sparse_3d_3$lognormconst,0),round(truelognormconst3d,0)) # Less accurate with 3 points in 3d
    expect_equal(round(norm_sparse_3d_5$lognormconst,2),round(truelognormconst3d,2))
    expect_equal(round(norm_sparse_3d_7$lognormconst,2),round(truelognormconst3d,2))

    expect_equal(round(norm_trust_3d_3$lognormconst,0),round(truelognormconst3d,0))
    expect_equal(round(norm_trust_3d_5$lognormconst,2),round(truelognormconst3d,2))
    expect_equal(round(norm_trust_3d_7$lognormconst,2),round(truelognormconst3d,2))

    expect_equal(round(norm_bfgs_3d_3$lognormconst,0),round(truelognormconst3d,0))
    expect_equal(round(norm_bfgs_3d_5$lognormconst,2),round(truelognormconst3d,2))
    expect_equal(round(norm_bfgs_3d_7$lognormconst,2),round(truelognormconst3d,2))

    expect_equal(round(norm_sr1_3d_3$lognormconst,0),round(truelognormconst3d,0))
    expect_equal(round(norm_sr1_3d_5$lognormconst,2),round(truelognormconst3d,2))
    expect_equal(round(norm_sr1_3d_7$lognormconst,2),round(truelognormconst3d,2))

    # Parameter vector reordering works
    # Note: some of these tests only work because the parameters are independent in this example!
    # 2d
    expect_equal(unique(round(norm_sparse_2d_3$nodesandweights$theta1,3)),unique(round(norm_sparse_2d_reorder_3$nodesandweights$theta1,3)))
    expect_equal(unique(round(norm_sparse_2d_3$nodesandweights$theta2,3)),unique(round(norm_sparse_2d_reorder_3$nodesandweights$theta2,3)))

    expect_equal(sort(norm_sparse_2d_3$nodesandweights$logpost),sort(norm_sparse_2d_reorder_3$nodesandweights$logpost))
    expect_equal(sort(norm_sparse_2d_3$nodesandweights$logpost_normalized),sort(norm_sparse_2d_reorder_3$nodesandweights$logpost_normalized))

    expect_equal(norm_sparse_2d_3$grid$features$m,norm_sparse_2d_reorder_3$grid$features$m[c(2,1)])
    expect_equal(norm_sparse_2d_3$grid$features$C,norm_sparse_2d_reorder_3$grid$features$C[c(2,1),c(2,1)])

    expect_equal(norm_sparse_2d_3$lognormconst,norm_sparse_2d_reorder_3$lognormconst)

    # 3d

    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta1,3)),unique(round(norm_sparse_3d_reorder_1$nodesandweights$theta1,3)))
    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta2,3)),unique(round(norm_sparse_3d_reorder_1$nodesandweights$theta2,3)))
    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta3,3)),unique(round(norm_sparse_3d_reorder_1$nodesandweights$theta3,3)))

    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta1,3)),unique(round(norm_sparse_3d_reorder_2$nodesandweights$theta1,3)))
    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta2,3)),unique(round(norm_sparse_3d_reorder_2$nodesandweights$theta2,3)))
    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta3,3)),unique(round(norm_sparse_3d_reorder_2$nodesandweights$theta3,3)))

    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta1,3)),unique(round(norm_sparse_3d_reorder_3$nodesandweights$theta1,3)))
    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta2,3)),unique(round(norm_sparse_3d_reorder_3$nodesandweights$theta2,3)))
    expect_equal(unique(round(norm_sparse_3d_3$nodesandweights$theta3,3)),unique(round(norm_sparse_3d_reorder_3$nodesandweights$theta3,3)))


    expect_equal(sort(norm_sparse_3d_3$nodesandweights$logpost),sort(norm_sparse_3d_reorder_1$nodesandweights$logpost))
    expect_equal(sort(norm_sparse_3d_3$nodesandweights$logpost_normalized),sort(norm_sparse_3d_reorder_1$nodesandweights$logpost_normalized))

    expect_equal(sort(norm_sparse_3d_3$nodesandweights$logpost),sort(norm_sparse_3d_reorder_2$nodesandweights$logpost))
    expect_equal(sort(norm_sparse_3d_3$nodesandweights$logpost_normalized),sort(norm_sparse_3d_reorder_2$nodesandweights$logpost_normalized))

    expect_equal(sort(norm_sparse_3d_3$nodesandweights$logpost),sort(norm_sparse_3d_reorder_3$nodesandweights$logpost))
    expect_equal(sort(norm_sparse_3d_3$nodesandweights$logpost_normalized),sort(norm_sparse_3d_reorder_3$nodesandweights$logpost_normalized))

    expect_equal(norm_sparse_3d_3$grid$features$m,norm_sparse_3d_reorder_1$grid$features$m[c(1,2,3)])
    expect_equal(norm_sparse_3d_3$grid$features$C,norm_sparse_3d_reorder_1$grid$features$C[c(1,2,3),c(1,2,3)])

    expect_equal(norm_sparse_3d_3$grid$features$m,norm_sparse_3d_reorder_2$grid$features$m[c(2,1,3)])
    expect_equal(norm_sparse_3d_3$grid$features$C,norm_sparse_3d_reorder_2$grid$features$C[c(2,1,3),c(2,1,3)])

    expect_equal(norm_sparse_3d_3$grid$features$m,norm_sparse_3d_reorder_3$grid$features$m[c(2,3,1)])
    expect_equal(norm_sparse_3d_3$grid$features$C,norm_sparse_3d_reorder_3$grid$features$C[c(2,3,1),c(2,3,1)])


    expect_equal(norm_sparse_3d_3$lognormconst,norm_sparse_3d_reorder_1$lognormconst)
    expect_equal(norm_sparse_3d_3$lognormconst,norm_sparse_3d_reorder_2$lognormconst)
    expect_equal(norm_sparse_3d_3$lognormconst,norm_sparse_3d_reorder_3$lognormconst)


    # Extraction of normalizing constants
    expect_equal(normconst1,thequadrature$normalized_posterior$lognormconst)
    expect_equal(normconst2,thelaplace$lognormconst)
    expect_equal(normconst3,themarginallaplace$normalized_posterior$lognormconst)

    # Extraction of number of quadrature points
    expect_equal(get_numquadpoints(thequadrature),3)
    expect_equal(get_numquadpoints(thequadrature_k7),25)
    expect_equal(get_numquadpoints(thequadrature3d),3)



    ## Nested, possibly sparse, quadrature ##
    expect_is(quadtable_p1_k1_prod,'data.frame')
    expect_is(quadtable_p2_k1_prod,'data.frame')
    expect_is(quadtable_p3_k1_prod,'data.frame')
    expect_is(quadtable_p4_k1_prod,'data.frame')
    expect_is(quadtable_p5_k1_prod,'data.frame')

    expect_is(quadtable_p1_k3_prod,'data.frame')
    expect_is(quadtable_p2_k3_prod,'data.frame')
    expect_is(quadtable_p3_k3_prod,'data.frame')
    expect_is(quadtable_p4_k3_prod,'data.frame')
    expect_is(quadtable_p5_k3_prod,'data.frame')

    expect_is(quadtable_p1_k5_prod,'data.frame')
    expect_is(quadtable_p2_k5_prod,'data.frame')
    expect_is(quadtable_p3_k5_prod,'data.frame')
    expect_is(quadtable_p4_k5_prod,'data.frame')
    expect_is(quadtable_p5_k5_prod,'data.frame')

    expect_equal(dim(quadtable_p1_k1_prod),c(1,2))
    expect_equal(dim(quadtable_p2_k1_prod),c(1,3))
    expect_equal(dim(quadtable_p3_k1_prod),c(1,4))
    expect_equal(dim(quadtable_p4_k1_prod),c(1,5))
    expect_equal(dim(quadtable_p5_k1_prod),c(1,6))

    expect_equal(dim(quadtable_p1_k3_prod),c(3,3))
    expect_equal(dim(quadtable_p2_k3_prod),c(9,4))
    expect_equal(dim(quadtable_p3_k3_prod),c(27,5))
    expect_equal(dim(quadtable_p4_k3_prod),c(81,6))
    expect_equal(dim(quadtable_p5_k3_prod),c(243,7))

    expect_equal(dim(quadtable_p1_k5_prod),c(9,4))
    expect_equal(dim(quadtable_p2_k5_prod),c(81,5))
    expect_equal(dim(quadtable_p3_k5_prod),c(729,6))
    expect_equal(dim(quadtable_p4_k5_prod),c(6561,7))
    expect_equal(dim(quadtable_p5_k5_prod),c(59049,8))

    expect_is(quadtable_p1_k1_sparse,'data.frame')
    expect_is(quadtable_p2_k1_sparse,'data.frame')
    expect_is(quadtable_p3_k1_sparse,'data.frame')
    expect_is(quadtable_p4_k1_sparse,'data.frame')
    expect_is(quadtable_p5_k1_sparse,'data.frame')

    expect_is(quadtable_p1_k3_sparse,'data.frame')
    expect_is(quadtable_p2_k3_sparse,'data.frame')
    expect_is(quadtable_p3_k3_sparse,'data.frame')
    expect_is(quadtable_p4_k3_sparse,'data.frame')
    expect_is(quadtable_p5_k3_sparse,'data.frame')

    expect_is(quadtable_p1_k5_sparse,'data.frame')
    expect_is(quadtable_p2_k5_sparse,'data.frame')
    expect_is(quadtable_p3_k5_sparse,'data.frame')
    expect_is(quadtable_p4_k5_sparse,'data.frame')
    expect_is(quadtable_p5_k5_sparse,'data.frame')

    expect_equal(dim(quadtable_p1_k1_sparse),c(1,2))
    expect_equal(dim(quadtable_p2_k1_sparse),c(1,3))
    expect_equal(dim(quadtable_p3_k1_sparse),c(1,4))
    expect_equal(dim(quadtable_p4_k1_sparse),c(1,5))
    expect_equal(dim(quadtable_p5_k1_sparse),c(1,6))

    expect_equal(dim(quadtable_p1_k3_sparse),c(3,3))
    expect_equal(dim(quadtable_p2_k3_sparse),c(9,4))
    expect_equal(dim(quadtable_p3_k3_sparse),c(19,5))
    expect_equal(dim(quadtable_p4_k3_sparse),c(33,6))
    expect_equal(dim(quadtable_p5_k3_sparse),c(51,7))

    expect_equal(dim(quadtable_p1_k5_sparse),c(9,4))
    expect_equal(dim(quadtable_p2_k5_sparse),c(37,5))
    expect_equal(dim(quadtable_p3_k5_sparse),c(93,6))
    expect_equal(dim(quadtable_p4_k5_sparse),c(201,7))
    expect_equal(dim(quadtable_p5_k5_sparse),c(433,8))
    # Not adapted
    expect_equal(exp(nq_p1_k1_prod),c('w1'=1))
    expect_equal(exp(nq_p2_k1_prod),c('w1'=1))
    expect_equal(exp(nq_p3_k1_prod),c('w1'=1))
    expect_equal(exp(nq_p4_k1_prod),c('w1'=1))
    expect_equal(exp(nq_p5_k1_prod),c('w1'=1))

    expect_equal(exp(nq_p1_k3_prod),c('w3'=1,'w1'=1))
    expect_equal(exp(nq_p2_k3_prod),c('w3'=1,'w1'=1))
    expect_equal(exp(nq_p3_k3_prod),c('w3'=1,'w1'=1))
    expect_equal(exp(nq_p4_k3_prod),c('w3'=1,'w1'=1))
    expect_equal(exp(nq_p5_k3_prod),c('w3'=1,'w1'=1))

    expect_equal(exp(nq_p1_k5_prod),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(nq_p2_k5_prod),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(nq_p3_k5_prod),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(nq_p4_k5_prod),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(nq_p5_k5_prod),c('w5'=1,'w3'=1,'w1'=1))

    expect_equal(exp(nq_p1_k1_sparse),c('w1'=1))
    expect_equal(exp(nq_p2_k1_sparse),c('w1'=1))
    expect_equal(exp(nq_p3_k1_sparse),c('w1'=1))
    expect_equal(exp(nq_p4_k1_sparse),c('w1'=1))
    expect_equal(exp(nq_p5_k1_sparse),c('w1'=1))

    expect_equal(exp(nq_p1_k3_sparse),c('w3'=1,'w1'=1))
    expect_equal(exp(nq_p2_k3_sparse),c('w3'=1,'w1'=1))
    expect_equal(exp(nq_p3_k3_sparse),c('w3'=1,'w1'=1))
    expect_equal(exp(nq_p4_k3_sparse),c('w3'=1,'w1'=1))
    expect_equal(exp(nq_p5_k3_sparse),c('w3'=1,'w1'=1))

    expect_equal(exp(nq_p1_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(nq_p2_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(nq_p3_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(nq_p4_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(nq_p5_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    # Adapted
    expect_equal(exp(anq_p1_k1_prod),c('w1'=1))
    expect_equal(exp(anq_p2_k1_prod),c('w1'=1))
    expect_equal(exp(anq_p3_k1_prod),c('w1'=1))
    expect_equal(exp(anq_p4_k1_prod),c('w1'=1))
    expect_equal(exp(anq_p5_k1_prod),c('w1'=1))

    expect_equal(exp(anq_p1_k3_prod),c('w3'=1,'w1'=1))
    expect_equal(exp(anq_p2_k3_prod),c('w3'=1,'w1'=1))
    expect_equal(exp(anq_p3_k3_prod),c('w3'=1,'w1'=1))
    expect_equal(exp(anq_p4_k3_prod),c('w3'=1,'w1'=1))
    expect_equal(exp(anq_p5_k3_prod),c('w3'=1,'w1'=1))

    expect_equal(exp(anq_p1_k5_prod),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(anq_p2_k5_prod),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(anq_p3_k5_prod),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(anq_p4_k5_prod),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(anq_p5_k5_prod),c('w5'=1,'w3'=1,'w1'=1))

    expect_equal(exp(anq_p1_k1_sparse),c('w1'=1))
    expect_equal(exp(anq_p2_k1_sparse),c('w1'=1))
    expect_equal(exp(anq_p3_k1_sparse),c('w1'=1))
    expect_equal(exp(anq_p4_k1_sparse),c('w1'=1))
    expect_equal(exp(anq_p5_k1_sparse),c('w1'=1))

    expect_equal(exp(anq_p1_k3_sparse),c('w3'=1,'w1'=1))
    expect_equal(exp(anq_p2_k3_sparse),c('w3'=1,'w1'=1))
    expect_equal(exp(anq_p3_k3_sparse),c('w3'=1,'w1'=1))
    expect_equal(exp(anq_p4_k3_sparse),c('w3'=1,'w1'=1))
    expect_equal(exp(anq_p5_k3_sparse),c('w3'=1,'w1'=1))

    expect_equal(exp(anq_p1_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(anq_p2_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(anq_p3_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(anq_p4_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))
    expect_equal(exp(anq_p5_k5_sparse),c('w5'=1,'w3'=1,'w1'=1))

})