test_that("maybe_subsample_draws() is a no-op when ndraws is NULL or too large", {
  draws <- matrix(1:20, nrow = 5, ncol = 4)
  expect_identical(maybe_subsample_draws(draws, NULL), draws)
  expect_identical(maybe_subsample_draws(draws, 10), draws)
})

test_that("maybe_subsample_draws() subsamples to the requested number of rows", {
  set.seed(1)
  draws <- matrix(1:20, nrow = 5, ncol = 4)
  sub <- maybe_subsample_draws(draws, 2)
  expect_equal(nrow(sub), 2)
  expect_equal(ncol(sub), 4)
})

test_that("maybe_subsample_draws_list() applies the same subsample to every matrix", {
  set.seed(1)
  a <- matrix(rep(1:5, 2), nrow = 5, ncol = 2)
  b <- matrix(rep(1:5, 2), nrow = 5, ncol = 2)
  sub <- maybe_subsample_draws_list(list(a = a, b = b), 2)

  expect_equal(nrow(sub$a), 2)
  expect_equal(nrow(sub$b), 2)
  # same rows were selected from both matrices
  expect_equal(sub$a[, 1], sub$b[, 1])
})

test_that("maybe_subsample_draws_list() is a no-op when ndraws is NULL", {
  a <- matrix(1:10, nrow = 5, ncol = 2)
  out <- maybe_subsample_draws_list(list(a = a), NULL)
  expect_identical(out$a, a)
})

test_that("summarize_draws_matrix() computes per-column point estimates and intervals", {
  draws <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  result <- summarize_draws_matrix(draws, conf.level = 1)

  expect_equal(result$.fitted, c(2, 5))
  expect_equal(result$.median, c(2, 5))
  expect_equal(result$.lower, c(1, 4))
  expect_equal(result$.upper, c(3, 6))
})

test_that("summarize_draws_list() prefixes columns by list name", {
  draws_list <- list(
    c1 = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2),
    c2 = matrix(c(10, 20, 30, 40, 50, 60), nrow = 3, ncol = 2)
  )
  result <- summarize_draws_list(draws_list, conf.level = 1)

  expect_setequal(
    names(result),
    c("c1.fitted", "c1.median", "c1.lower", "c1.upper", "c2.fitted", "c2.median", "c2.lower", "c2.upper")
  )
  expect_equal(result$c1.fitted, c(2, 5))
  expect_equal(result$c2.fitted, c(20, 50))
})

test_that("summarize_curve_draws() summarizes an (draws x cells x ages) array into long format", {
  curve_draws <- array(0, dim = c(3, 2, 2))
  curve_draws[, 1, 1] <- c(1, 2, 3)
  curve_draws[, 1, 2] <- c(4, 5, 6)
  curve_draws[, 2, 1] <- c(7, 8, 9)
  curve_draws[, 2, 2] <- c(10, 11, 12)

  newdata <- data.frame(group = c("A", "B"))
  age_grid <- c(20, 30)

  result <- summarize_curve_draws(curve_draws, newdata, age_grid, conf.level = 1)

  expect_equal(nrow(result), 4)
  expect_equal(result$age, rep(age_grid, times = 2))
  expect_equal(result$group, rep(c("A", "B"), each = 2))
  expect_equal(result$.fitted, c(2, 5, 8, 11))
  expect_equal(result$.lower, c(1, 4, 7, 10))
  expect_equal(result$.upper, c(3, 6, 9, 12))
})

test_that("schedule_cell_vars() extracts covariates from schedule formulas, excluding age", {
  eff <- pk_effects(c1 = ~1 + group, mu1 = ~1 + group)
  fake <- list(effects = eff, age_col = "age")
  expect_setequal(schedule_cell_vars(fake), "group")
})

test_that("schedule_cell_vars() returns character(0) for an intercept-only model", {
  fake <- list(effects = pk_effects(), age_col = "age")
  expect_length(schedule_cell_vars(fake), 0)
})

test_that("schedule_cell_vars() ignores prior_scales/prior_centers/tau_priors overrides (regression)", {
  # Prior to the plot.pk_fit() refactor, unique_cells() called all.vars() over
  # every element of object$effects, including the non-formula prior_scales/
  # prior_centers/tau_priors list entries. Guard against that regressing.
  eff <- pk_effects(
    c1 = ~1 + group,
    prior_scales = list(c1 = list(group = 0.5)),
    tau_priors = list(c1 = list("s(x)" = list(scale = 0.5)))
  )
  fake <- list(effects = eff, age_col = "age")
  expect_no_error(vars <- schedule_cell_vars(fake))
  expect_setequal(vars, "group")
})

test_that("unique_cells() deduplicates training data by schedule covariates", {
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data)
  cells <- unique_cells(fake)

  expect_setequal(names(cells), "group")
  expect_setequal(as.character(cells$group), levels(data$group))
  expect_equal(nrow(cells), length(unique(data$group)))
})
