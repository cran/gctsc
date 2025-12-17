test_that("poisson.marg returns object of correct class", {
  m <- poisson.marg()
  expect_true(is.list(m))
  expect_equal(class(m), "marginal.gctsc")
})

test_that("poisson.marg has expected components", {
  m <- poisson.marg()
  expect_true("start" %in% names(m))
  expect_true("npar"  %in% names(m))
  expect_true("bounds" %in% names(m))
})

test_that("poisson.marg npar returns number of predictors", {
  m <- poisson.marg()
  X <- matrix(1, 10, 2)
  expect_equal(m$npar(X), 2)
})

test_that("poisson.marg bounds works for small example", {
  m <- poisson.marg()
  y <- c(0, 1)
  X <- matrix(1, 2, 1)
  lambda <- 1
  ab <- m$bounds(y, X, lambda)
  expect_equal(dim(ab), c(2,2))
})

