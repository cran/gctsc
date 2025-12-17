test_that("arma.cormat returns correct structure", {
  obj <- arma.cormat(p = 1, q = 1)
  
  expect_true(is.list(obj))
  expect_true("npar" %in% names(obj))
  expect_true("start" %in% names(obj))
  expect_true("od" %in% names(obj))
  expect_equal(class(obj), c("arma.gctsc", "cormat.gctsc"))
})
test_that("arma.cormat npar is p + q", {
  obj <- arma.cormat(p = 2, q = 1)
  expect_equal(obj$npar, 3)
})
test_that("arma.cormat stores ARMA order correctly", {
  obj <- arma.cormat(p = 1, q = 2)
  expect_equal(obj$od, c(1, 2))
})
test_that("arma.cormat start() returns named parameters", {
  set.seed(1)
  
  y <- rnorm(20)
  obj <- arma.cormat(p = 1, q = 1)
  
  tau <- obj$start(y)
  
  expect_equal(length(tau), 2)
  expect_equal(names(tau), c("ar1", "ma1"))
})
test_that("arma.cormat start() attaches bounds if given", {
  y <- rnorm(20)
  obj <- arma.cormat(p = 1, q = 0, tau.lower = -0.9, tau.upper = 0.9)
  
  tau <- obj$start(y)
  
  expect_equal(attr(tau, "lower"), -0.9)
  expect_equal(attr(tau, "upper"), 0.9)
})
