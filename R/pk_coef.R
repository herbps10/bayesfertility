#' Extract posterior draws or summaries for a specific term
#'
#' @param fit a pk_fit object
#' @param parameter one of "c1", "c2", "mu1", "mu2", "sigma1", "sigma2"
#' @param term the term label
#' @param summarize if TRUE, return a tibble of posterior summaries;
#'   if FALSE, return the n_draws x n_coef matrix of draws
#' @param source "posterior" (default) or "prior"
#' @param conf.level credible interval width
#'
#' @export
pk_coef <- function(
  fit,
  parameter,
  term,
  summarize = TRUE,
  source = "posterior",
  conf.level = 0.95
) {
  checkmate::assert_class(fit, "pk_fit")
  checkmate::assert_choice(
    parameter,
    choices = c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")
  )

  terms_info <- fit$processed[[parameter]]$terms_info
  avail <- vapply(terms_info, function(t) t$label, character(1))
  checkmate::assert_choice(term, choices = avail)

  term_info <- Find(function(t) t$label == term, terms_info)
  col_idx <- term_info$col_idx

  beta <- switch(
    source,
    posterior = fit$beta_draws[[parameter]],
    prior = fit$prior_beta_draws[[parameter]]
  )

  if (is.null(beta)) {
    cli::cli_abort(c(
      "No {source} draws available.",
      "i" = if (source == "prior") "Refit with sample_prior = \"yes\"."
    ))
  }

  draws <- beta[, col_idx, drop = FALSE]
  level_names <- colnames(fit$processed[[parameter]]$X)[col_idx]
  if (is.null(level_names)) {
    level_names <- paste0(term, "[", seq_along(col_idx), "]")
  }
  colnames(draws) <- level_names

  if (!summarize) {
    return(draws)
  }

  alpha <- (1 - conf.level) / 2

  tibble::tibble(
    level = level_names,
    mean = matrixStats::colMeans2(draws),
    median = matrixStats::colMedians(draws),
    lower = matrixStats::colQuantiles(draws, probs = alpha),
    upper = matrixStats::colQuantiles(draws, probs = 1 - alpha)
  )
}
