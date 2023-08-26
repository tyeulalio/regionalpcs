## - check functionality
test_that("get_sig_pcs() returns correct dimension", {
   expect_equal(get_sig_pcs(t(mtcars))$est_dim, 3)
})

## -- checking parameters
test_that("get_sig_pcs() checks method parameter", {
  expect_error(get_sig_pcs(t(mtcars), "method"))
})

test_that("get_sig_pcs() checks x parameter for matrix", {
  expect_equal(get_sig_pcs(as.matrix(t(mtcars)), "gd")$est_dim, 3)
})

test_that("get_sig_pcs() checks x parameter for dataframe", {
  expect_equal(get_sig_pcs(as.data.frame(t(mtcars)), "gd")$est_dim, 3)
})

test_that("get_sig_pcs() checks x parameter for wrong input", {
  expect_error(get_sig_pcs("x", "gd"))
})
