#' Extract posterior draws or summaries for a specific term
#'
#' @description
#' Extracts the working-scale regression coefficients for a single term of a
#' single schedule parameter's formula (e.g. the `race` term of `c1`'s
#' formula), either as raw draws or summarized with a point estimate and
#' credible interval. For a term with multiple levels or basis columns (a
#' factor, or a smooth/random-effect term), one column is returned per level.
#'
#' @param fit A fitted `pk_fit` object, as returned by \code{pk_fit()}.
#' @param parameter One of `"c1"`, `"c2"`, `"mu1"`, `"mu2"`, `"sigma1"`,
#'   `"sigma2"` — which schedule parameter's formula the term belongs to.
#' @param term The term label to extract, as it appears in that parameter's
#'   formula (e.g. `"(Intercept)"`, `"race"`, `"s(province)"`).
#' @param summarize If `TRUE` (default), return a tibble of posterior
#'   summaries; if `FALSE`, return the `n_draws` x `n_coef` matrix of draws.
#' @param source `"posterior"` (default) or `"prior"`; `"prior"` requires the
#'   fit to have been created with `sample_prior = "yes"` or `"only"`.
#' @param conf.level Credible interval width used when `summarize = TRUE`.
#'
#' @return If `summarize`, a tibble with one row per level/column of `term`
#'   and columns `level`, `mean`, `median`, `lower`, `upper`. If not, an
#'   `n_draws` x `n_levels` matrix of draws, with columns named by level.
#'
#' @examples
#' \dontrun{
#' pk_coef(fit, parameter = "c1", term = "race")
#' pk_coef(fit, parameter = "c1", term = "(Intercept)", summarize = FALSE)
#' }
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
