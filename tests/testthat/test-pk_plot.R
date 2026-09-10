test_that("plot.pk_fit() returns a ggplot object", {
  fake <- make_fake_pk_fit()
  p <- plot(fake)
  expect_s3_class(p, "ggplot")
})

test_that("plot.pk_fit() facets by covariate when effects reference one, and not otherwise", {
  data <- make_synthetic_data()

  fake_no_covariate <- make_fake_pk_fit(effects = pk_effects(), data = data)
  p1 <- plot(fake_no_covariate)
  expect_s3_class(p1$facet, "FacetNull")

  fake_with_covariate <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data)
  p2 <- plot(fake_with_covariate)
  expect_s3_class(p2$facet, "FacetWrap")
})

test_that("plot.pk_fit() errors when prior draws are requested but unavailable", {
  fake <- make_fake_pk_fit(with_prior = FALSE)
  expect_error(plot(fake, source = "prior"))
})

test_that("plot_schedule_curve() adds an observed-data layer only when obs_df is supplied", {
  curve_df <- data.frame(age = c(20, 30), .fitted = c(0.1, 0.2), .lower = c(0.05, 0.15), .upper = c(0.15, 0.25))

  p1 <- plot_schedule_curve(curve_df)
  expect_s3_class(p1, "ggplot")
  expect_length(p1$layers, 2) # ribbon + line

  obs_df <- data.frame(age = c(20, 30), asfr = c(0.11, 0.19), .lower = c(0.08, 0.16), .upper = c(0.14, 0.22))
  p2 <- plot_schedule_curve(curve_df, obs_df = obs_df)
  expect_length(p2$layers, 4) # ribbon + line + errorbar + point
})
