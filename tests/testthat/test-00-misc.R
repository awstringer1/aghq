context("Misc")

test_that("Misc functions work", {
  tt <- c(1,2,4,5)
  expect_equal(splice(tt,3,1),c(3,1,2,4,5))
  expect_equal(splice(tt,3,2),c(1,3,2,4,5))
  expect_equal(splice(tt,3,3),c(1,2,3,4,5))
  expect_equal(splice(tt,3,4),c(1,2,4,3,5))
  expect_equal(splice(tt,3,5),c(1,2,4,5,3))

})
