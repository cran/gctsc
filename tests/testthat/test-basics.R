
test_that("basic functions exist", {
  expect_true(exists("poisson.marg", where = asNamespace("gctsc")))
  expect_true(exists("arma.cormat",  where = asNamespace("gctsc")))
})
    