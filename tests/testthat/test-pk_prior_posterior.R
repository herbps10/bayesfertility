test_that("pk_prior_posterior() type = 'coefficients' returns a ggplot", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data, with_prior = TRUE)
  p <- pk_prior_posterior(fake, type = "coefficients")
  expect_s3_class(p, "ggplot")
})

test_that("pk_prior_posterior() type = 'schedule_params' returns a ggplot", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data, with_prior = TRUE)
  p <- pk_prior_posterior(fake, type = "schedule_params")
  expect_s3_class(p, "ggplot")
})

test_that("pk_prior_posterior() type = 'tfr' returns a ggplot", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data, with_prior = TRUE)
  p <- pk_prior_posterior(fake, type = "tfr")
  expect_s3_class(p, "ggplot")
})

test_that("pk_prior_posterior() type = 'tau' returns a ggplot for a model with penalized terms", {
  data <- data.frame(province = factor(rep(paste0("p", 1:6), each = 5)))
  data$age <- rep(seq(15, 45, by = 10), length.out = nrow(data))
  data$exposure <- 1000
  set.seed(1)
  data$births <- rpois(nrow(data), lambda = 50)

  fake <- make_fake_pk_fit(
    effects = pk_effects(c1 = ~1 + s(province, bs = "re")),
    data = data,
    with_prior = TRUE
  )
  p <- pk_prior_posterior(fake, type = "tau")
  expect_s3_class(p, "ggplot")
})

test_that("pk_prior_posterior() errors when no prior draws are available", {
  fake <- make_fake_pk_fit(with_prior = FALSE)
  expect_error(pk_prior_posterior(fake, type = "coefficients"), "prior")
})
