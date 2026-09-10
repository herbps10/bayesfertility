test_that("pk_coef() returns a tibble when summarize = TRUE", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data)

  result <- pk_coef(fake, parameter = "c1", term = "group")
  expect_s3_class(result, "tbl_df")
  expect_setequal(names(result), c("level", "mean", "median", "lower", "upper"))
})

test_that("pk_coef() returns a matrix when summarize = FALSE", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data, n_draws = 15)

  result <- pk_coef(fake, parameter = "c1", term = "group", summarize = FALSE)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 15)
})

test_that("pk_coef() errors on an unknown parameter or term", {
  fake <- make_fake_pk_fit()
  expect_error(pk_coef(fake, parameter = "not_a_param", term = "(Intercept)"))
  expect_error(pk_coef(fake, parameter = "c1", term = "not_a_term"))
})
