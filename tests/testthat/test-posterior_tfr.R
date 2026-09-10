test_that("posterior_tfr() returns a tibble when summarize = TRUE", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data)

  result <- posterior_tfr(fake)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c(".fitted", ".median", ".lower", ".upper") %in% names(result)))
  expect_equal(nrow(result), length(unique(data$group)))
})

test_that("posterior_tfr() returns a matrix when summarize = FALSE", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data, n_draws = 15)

  result <- posterior_tfr(fake, summarize = FALSE)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 15)
  expect_equal(ncol(result), length(unique(data$group)))
})
