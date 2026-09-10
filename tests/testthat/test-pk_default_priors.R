param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")

test_that("pk_default_prior_centers() returns the expected structure for every parameter", {
  for (p in param_names) {
    centers <- pk_default_prior_centers(p)
    expect_setequal(names(centers), c("(Intercept)", "parametric", "smooth"))
    expect_true(is.numeric(centers[["(Intercept)"]]))
    expect_equal(centers[["parametric"]], 0)
    expect_equal(centers[["smooth"]], 0)
  }
})

test_that("pk_default_prior_scales() returns the expected structure for every parameter", {
  for (p in param_names) {
    scales <- pk_default_prior_scales(p)
    expect_setequal(names(scales), c("(Intercept)", "parametric", "smooth"))
    expect_true(is.numeric(scales[["(Intercept)"]]))
    expect_gt(scales[["(Intercept)"]], 0)
    expect_true(is.numeric(scales[["parametric"]]))
    expect_null(scales[["smooth"]])
  }
})

test_that("pk_default_tau_scales() returns a single positive number for every parameter", {
  for (p in param_names) {
    scale <- pk_default_tau_scales(p)
    expect_length(scale, 1)
    expect_gt(scale, 0)
  }
})
