test_that("pk_effects() requires one-sided formulas", {
  expect_error(pk_effects(c1 = births ~ age), "one-sided")
})

test_that("pk_effects() rejects non-formula arguments", {
  expect_error(pk_effects(c1 = "not a formula"))
})

test_that("pk_effects() defaults to intercept-only formulas for all six parameters", {
  eff <- pk_effects()
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")
  for (p in param_names) {
    expect_identical(as.character(eff[[p]]), as.character(~1))
  }
  expect_s3_class(eff, "pk_effects")
})

test_that("pk_effects() stores user-supplied formulas", {
  eff <- pk_effects(c1 = ~1 + race, mu2 = ~1 + province)
  expect_identical(as.character(eff$c1), as.character(~1 + race))
  expect_identical(as.character(eff$mu2), as.character(~1 + province))
})

test_that("pk_effects() rejects unknown parameter names in prior_scales", {
  expect_error(
    pk_effects(prior_scales = list(not_a_param = list(`(Intercept)` = 1))),
    "Unknown schedule parameter"
  )
})

test_that("pk_effects() rejects unknown parameter names in tau_priors", {
  expect_error(
    pk_effects(tau_priors = list(not_a_param = list("s(x)" = list(scale = 1)))),
    "Unknown schedule parameter"
  )
})

test_that("pk_effects() validates prior_scales values are non-negative numbers", {
  expect_error(
    pk_effects(prior_scales = list(c1 = list(`(Intercept)` = -1)))
  )
})

test_that("pk_effects() fills in empty lists for unspecified prior overrides", {
  eff <- pk_effects()
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")
  expect_identical(names(eff$prior_scales), param_names)
  expect_identical(names(eff$prior_centers), param_names)
  expect_identical(names(eff$tau_priors), param_names)
  for (p in param_names) {
    expect_length(eff$prior_scales[[p]], 0)
  }
})