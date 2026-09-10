test_that("pk_fit() fits a real (tiny) model end-to-end", {
  skip_if_not(
    identical(Sys.getenv("BAYESFERTILITY_DEV"), "TRUE") &&
      requireNamespace("cmdstanr", quietly = TRUE) &&
      !inherits(try(cmdstanr::cmdstan_path(), silent = TRUE), "try-error"),
    "cmdstan not available in dev mode"
  )

  data <- make_synthetic_data()
  fit <- pk_fit(
    data,
    "births",
    "exposure",
    "age",
    effects = pk_effects(),
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    refresh = 0,
    sample_prior = "no"
  )

  expect_s3_class(fit, "pk_fit")
  expect_false(is.null(fit$fit))
  expect_equal(nrow(fit$beta_draws$c1), 50)
})
