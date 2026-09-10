test_that("print.pk_fit() prints key summary sections", {
  fake <- make_fake_pk_fit()
  expect_message(print(fake), "Hierarchical Peristera-Kostaki")
  expect_message(print(fake), "Chains")
  expect_message(print(fake), "Max R-hat")
})

test_that("print.pk_fit() returns its input invisibly", {
  fake <- make_fake_pk_fit()
  out <- utils::capture.output(result <- withVisible(print(fake)))
  expect_false(result$visible)
  expect_identical(result$value, fake)
})

test_that("print.pk_fit() only prints the six schedule-parameter formulas (regression)", {
  # Previously computed param_names <- names(x$effects), which has 9 elements
  # (the 6 schedule formulas plus prior_scales/prior_centers/tau_priors), so
  # it tried to deparse those list-valued elements as formulas too.
  fake <- make_fake_pk_fit()
  msgs <- capture_messages(print(fake))
  output <- paste(msgs, collapse = "\n")
  expect_false(grepl("prior_scales", output, fixed = TRUE))
  expect_false(grepl("prior_centers", output, fixed = TRUE))
  expect_false(grepl("tau_priors", output, fixed = TRUE))
})
