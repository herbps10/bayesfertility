test_that("pp_check.pk_fit() returns a ggplot for dens_overlay and ribbon", {
  fake <- make_fake_pk_fit(n_draws = 10)

  p1 <- pp_check(fake, type = "dens_overlay")
  expect_s3_class(p1, "ggplot")

  p2 <- pp_check(fake, type = "ribbon")
  expect_s3_class(p2, "ggplot")
})

test_that("pp_check.pk_fit() type = 'stat_grouped' runs without error (regression test)", {
  # pp_check.pk_fit() previously called pp_check_stat_grouped(y_obs, yrep, group,
  # stat, ...) - missing the `newdata` argument, which shifted `group` into
  # newdata's slot and `stat` into group's. Confirm the fixed call works.
  data <- make_synthetic_data()
  fake <- make_fake_pk_fit(effects = pk_effects(c1 = ~1 + group), data = data, n_draws = 10)

  p <- pp_check(fake, type = "stat_grouped", group = "group", stat = "mean")
  expect_s3_class(p, "ggplot")
})

test_that("pp_check.pk_fit() type = 'stat_grouped' errors when group is missing", {
  fake <- make_fake_pk_fit(n_draws = 10)
  expect_error(pp_check(fake, type = "stat_grouped"))
})
