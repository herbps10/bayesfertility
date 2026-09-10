test_that("formula_has_smooth() detects smooth terms", {
  expect_false(formula_has_smooth(~1))
  expect_false(formula_has_smooth(~1 + group))
  expect_true(formula_has_smooth(~1 + s(x)))
  expect_true(formula_has_smooth(~1 + te(x, z)))
})

test_that("process_formula() handles an intercept-only formula", {
  data <- make_synthetic_data()
  info <- process_formula(~1, data)

  expect_equal(info$n_coef, 1)
  expect_equal(info$n_smooths, 0L)
  expect_null(info$pregam)
  expect_length(info$penalties, 0)
  expect_equal(info$terms_info[[1]]$label, "(Intercept)")
  expect_equal(info$terms_info[[1]]$type, "parametric")
  expect_identical(info$terms_info[[1]]$col_idx, 1L)
})

test_that("process_formula() groups parametric columns by term for a factor covariate", {
  data <- make_synthetic_data()
  info <- process_formula(~1 + group, data)

  labels <- vapply(info$terms_info, function(t) t$label, character(1))
  expect_setequal(labels, c("(Intercept)", "group"))
  expect_equal(info$n_coef, ncol(info$X))
  expect_null(info$pregam)
})

test_that("process_formula() builds penalties for a smooth/random-effect term", {
  data <- data.frame(
    province = factor(rep(paste0("p", 1:6), each = 5)),
    x = rnorm(30)
  )
  info <- process_formula(~1 + s(province, bs = "re"), data)

  expect_equal(info$n_smooths, 1L)
  expect_true(!is.null(info$pregam))
  expect_gt(length(info$penalties), 0)

  types <- vapply(info$terms_info, function(t) t$type, character(1))
  expect_true("smooth" %in% types)
})

test_that("build_terms_info() combines parametric and smooth terms in order", {
  out <- build_terms_info(
    parametric_terms = list("(Intercept)" = 1L, race = 2:3),
    smooth_labels = "s(province)",
    smooth_col_idx = list(4:6),
    smooth_penalty_idx = list(1L)
  )

  expect_length(out, 3)
  expect_equal(vapply(out, function(t) t$label, character(1)), c("(Intercept)", "race", "s(province)"))
  expect_equal(out[[2]]$type, "parametric")
  expect_equal(out[[3]]$type, "smooth")
  expect_identical(out[[3]]$col_idx, 4:6)
  expect_identical(out[[3]]$penalty_idx, 1L)
})

test_that("resolve_priors() applies package defaults with no user overrides", {
  data <- make_synthetic_data()
  info <- process_formula(~1 + group, data)$terms_info
  resolved <- resolve_priors(info, user_scales = list(), user_centers = list(), param_name = "c1")

  expect_length(resolved$scales, length(info))
  expect_length(resolved$centers, length(info))

  intercept_pos <- which(vapply(info, function(t) t$label, character(1)) == "(Intercept)")
  expect_equal(resolved$centers[intercept_pos], pk_default_prior_centers("c1")[["(Intercept)"]])
})

test_that("resolve_priors() honors user overrides", {
  data <- make_synthetic_data()
  info <- process_formula(~1 + group, data)$terms_info
  resolved <- resolve_priors(info, user_scales = list(group = 2.5), user_centers = list(), param_name = "c1")

  group_pos <- which(vapply(info, function(t) t$label, character(1)) == "group")
  expect_equal(resolved$scales[group_pos], 2.5)
})

test_that("resolve_priors() errors on an unknown term name", {
  data <- make_synthetic_data()
  info <- process_formula(~1 + group, data)$terms_info
  expect_error(
    resolve_priors(info, user_scales = list(not_a_term = 1), user_centers = list(), param_name = "c1"),
    "Unknown term name"
  )
})

test_that("resolve_tau_priors() defaults to the package tau scale and honors overrides", {
  data <- data.frame(
    province = factor(rep(paste0("p", 1:6), each = 5)),
    x = rnorm(30)
  )
  info <- process_formula(~1 + s(province, bs = "re"), data)
  n_pen <- length(info$penalties)

  default_resolved <- resolve_tau_priors(info$terms_info, n_pen, list(), "c1")
  expect_length(default_resolved, n_pen)
  expect_true(all(default_resolved == pk_default_tau_scales("c1")))

  overridden <- resolve_tau_priors(
    info$terms_info,
    n_pen,
    list("s(province)" = list(scale = 3)),
    "c1"
  )
  expect_true(all(overridden == 3))
})

test_that("resolve_tau_priors() errors on an unknown smooth label", {
  data <- data.frame(
    province = factor(rep(paste0("p", 1:6), each = 5)),
    x = rnorm(30)
  )
  info <- process_formula(~1 + s(province, bs = "re"), data)
  n_pen <- length(info$penalties)

  expect_error(
    resolve_tau_priors(info$terms_info, n_pen, list(not_a_smooth = list(scale = 1)), "c1"),
    "Unknown smooth label"
  )
})

test_that("process_effects() processes all six schedule parameters", {
  data <- make_synthetic_data()
  effects <- pk_effects(c1 = ~1 + group)
  processed <- process_effects(effects, data)

  expect_setequal(names(processed), c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2"))
  for (p in names(processed)) {
    expect_length(processed[[p]]$prior_scales, length(processed[[p]]$terms_info))
    expect_length(processed[[p]]$prior_centers, length(processed[[p]]$terms_info))
  }
  expect_equal(processed$c1$n_coef, 2) # (Intercept) + group
})
