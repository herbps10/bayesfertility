test_that("resolve_exposure() resolves NULL, column name, and numeric-vector cases", {
  newdata <- data.frame(exposure = c(100, 200, 300), other_exposure = c(1, 2, 3))
  object <- list(exposure_col = "exposure")

  expect_equal(resolve_exposure(NULL, newdata, object), c(100, 200, 300))
  expect_equal(resolve_exposure("other_exposure", newdata, object), c(1, 2, 3))
  expect_equal(resolve_exposure(c(5, 6, 7), newdata, object), c(5, 6, 7))
})

test_that("resolve_exposure() errors on an invalid exposure argument type", {
  newdata <- data.frame(exposure = c(100, 200, 300))
  object <- list(exposure_col = "exposure")
  expect_error(resolve_exposure(list(1, 2, 3), newdata, object))
})

test_that("posterior_predict.pk_fit() returns an integer n_draws x n_obs matrix", {
  fake <- make_fake_pk_fit(n_draws = 10)
  yrep <- posterior_predict(fake)

  expect_true(is.matrix(yrep))
  expect_equal(dim(yrep), c(10, nrow(fake$data)))
  expect_true(all(yrep == round(yrep)))
  expect_true(all(yrep >= 0))
})
