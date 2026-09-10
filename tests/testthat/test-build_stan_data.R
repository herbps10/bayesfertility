test_that("build_stan_data() builds correct top-level fields", {
  data <- make_synthetic_data()
  processed <- process_effects(pk_effects(c1 = ~1 + group), data)
  sd <- build_stan_data(processed, data, "births", "exposure", "age")

  expect_equal(sd$N, nrow(data))
  expect_equal(sd$A, length(unique(data$age)))
  expect_equal(sd$B, data$births)
  expect_equal(sd$E, data$exposure)
  expect_equal(sd$age, match(data$age, sort(unique(data$age))))
  expect_equal(sd$ages, sort(unique(data$age)))
})

test_that("build_stan_data() builds per-parameter design matrix fields", {
  data <- make_synthetic_data()
  processed <- process_effects(pk_effects(c1 = ~1 + group), data)
  sd <- build_stan_data(processed, data, "births", "exposure", "age")

  expect_equal(sd$P_c1, processed$c1$n_coef)
  expect_equal(dim(sd$X_c1), c(nrow(data), processed$c1$n_coef))
  expect_equal(sd$n_terms_c1, length(processed$c1$terms_info))
  expect_length(sd$term_starts_c1, sd$n_terms_c1)
  expect_length(sd$term_sizes_c1, sd$n_terms_c1)

  # intercept-only parameters have a single term
  expect_equal(sd$P_mu1, 1)
  expect_equal(sd$n_terms_mu1, 1)
})

test_that("build_stan_data() fills in empty penalty fields when there are no smooths", {
  data <- make_synthetic_data()
  processed <- process_effects(pk_effects(c1 = ~1 + group), data)
  sd <- build_stan_data(processed, data, "births", "exposure", "age")

  expect_equal(sd$n_penalties_c1, 0)
  expect_length(sd$penalty_sizes_c1, 0)
  expect_equal(dim(sd$S_packed_c1), c(0, 0))
})

test_that("build_stan_data() builds penalty bookkeeping for a smooth term", {
  data <- data.frame(
    age = rep(seq(15, 45, by = 5), 6),
    province = factor(rep(paste0("p", 1:6), each = 7)),
    exposure = 1000
  )
  set.seed(1)
  data$births <- rpois(nrow(data), lambda = 50)

  processed <- process_effects(pk_effects(c1 = ~1 + s(province, bs = "re")), data)
  sd <- build_stan_data(processed, data, "births", "exposure", "age")

  expect_gt(sd$n_penalties_c1, 0)
  expect_length(sd$penalty_sizes_c1, sd$n_penalties_c1)
  expect_equal(sum(sd$penalty_sizes_c1), nrow(sd$S_packed_c1))
})
