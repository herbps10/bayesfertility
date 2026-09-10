test_that("eta_to_params() applies the correct working-scale transforms", {
  eta <- list(
    c1 = matrix(log(2), 1, 1),
    c2 = matrix(log(3), 1, 1),
    mu1 = matrix(20, 1, 1),
    mu2 = matrix(log(5), 1, 1),
    sigma1 = matrix(log(4), 1, 1),
    sigma2 = matrix(log(6), 1, 1)
  )
  params <- eta_to_params(eta)

  expect_equal(params$c1[1, 1], 2)
  expect_equal(params$c2[1, 1], 3)
  expect_equal(params$mu1[1, 1], 20)
  expect_equal(params$mu2[1, 1], 25)
  expect_equal(params$sigma1[1, 1], 4)
  expect_equal(params$sigma2[1, 1], 6)
})

test_that("evaluate_schedule_grid() evaluates the PK schedule at each age", {
  params <- list(
    c1 = matrix(1, 1, 1), c2 = matrix(0, 1, 1),
    mu1 = matrix(25, 1, 1), mu2 = matrix(35, 1, 1),
    sigma1 = matrix(5, 1, 1), sigma2 = matrix(5, 1, 1)
  )
  result <- evaluate_schedule_grid(params, age_grid = c(25, 30))

  expect_equal(dim(result), c(1, 1, 2))
  expect_equal(result[1, 1, 1], 1)
  expect_equal(result[1, 1, 2], exp(-1))
})

test_that("resolve_age() resolves NULL, column name, and numeric-vector cases", {
  newdata <- data.frame(age = c(20, 30, 40), other_age = c(1, 2, 3))
  object <- list(age_col = "age")

  expect_equal(resolve_age(NULL, newdata, object), c(20, 30, 40))
  expect_equal(resolve_age("other_age", newdata, object), c(1, 2, 3))
  expect_equal(resolve_age(c(5, 6, 7), newdata, object), c(5, 6, 7))
})

test_that("resolve_age() errors on an invalid age argument type", {
  newdata <- data.frame(age = c(20, 30, 40))
  object <- list(age_col = "age")
  expect_error(resolve_age(list(1, 2, 3), newdata, object))
})

test_that("predict.pk_fit() type = 'asfr' returns correct shapes", {
  fake <- make_fake_pk_fit()

  summarized <- predict(fake, type = "asfr")
  expect_s3_class(summarized, "tbl_df")
  expect_equal(nrow(summarized), nrow(fake$data))

  raw <- predict(fake, type = "asfr", summarize = FALSE)
  expect_true(is.matrix(raw))
  expect_equal(ncol(raw), nrow(fake$data))
})

test_that("predict.pk_fit() type = 'linear_predictor'/'schedule_params' return correct shapes", {
  fake <- make_fake_pk_fit()

  for (type in c("linear_predictor", "schedule_params")) {
    summarized <- predict(fake, type = type)
    expect_s3_class(summarized, "tbl_df")
    expect_equal(nrow(summarized), nrow(fake$data))

    raw <- predict(fake, type = type, summarize = FALSE)
    expect_type(raw, "list")
    expect_setequal(names(raw), c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2"))
  }
})

test_that("predict.pk_fit() type = 'schedule_curve' returns correct shapes", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data)
  cell_data <- data.frame(group = factor(c("A", "B"), levels = levels(data$group)))
  age_grid <- seq(15, 45, by = 5)

  summarized <- predict(fake, newdata = cell_data, type = "schedule_curve", age_grid = age_grid)
  expect_s3_class(summarized, "tbl_df")
  expect_equal(nrow(summarized), length(age_grid) * nrow(cell_data))

  raw <- predict(fake, newdata = cell_data, type = "schedule_curve", age_grid = age_grid, summarize = FALSE)
  expect_equal(dim(raw), c(25, nrow(cell_data), length(age_grid)))
})

test_that("predict.pk_fit() errors when prior draws are requested but unavailable", {
  fake <- make_fake_pk_fit(with_prior = FALSE)
  expect_error(predict(fake, type = "asfr", source = "prior"))
})
