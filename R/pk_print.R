#' Print a fitted PK model
#'
#' @description
#' Prints a human-readable summary of a \code{pk_fit} object: the per-parameter
#' effects specification, the number of training observations, sampling
#' metadata (chains, iterations, runtime), and posterior convergence
#' diagnostics (max R-hat, min effective sample size, and number of divergent
#' transitions), each flagged as passing or failing standard thresholds.
#'
#' @param x A fitted `pk_fit` object, as returned by \code{pk_fit()}.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns `x`. Called for its side effect of printing to
#'   the console.
#'
#' @importFrom stringr str_replace
#' @importFrom prettyunits pretty_sec
#' @export
print.pk_fit <- function(x, ...) {
  cli::cli_h1("Hierarchical Peristera-Kostaki Model 2")

  # Effects spec
  cli::cli_h3("Effects specification")
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")
  for (p in param_names) {
    formula_str <- deparse(x$effects[[p]], width.cutoff = 500L)
    formula_str <- str_replace(formula_str, "^~\\s*", "")
    cli::cli_li("{.field {p}} ~ {formula_str}", )
  }

  # Data summary
  cli::cli_h3("Data")
  n_obs <- nrow(x$data)
  cli::cli_text("Observations: {.val {n_obs}}")

  # Sampling info
  cli::cli_h3("Sampling")
  meta <- x$fit$metadata()
  n_chains <- meta$num_chains
  n_iter <- meta$iter_sampling
  n_warmup <- meta$iter_warmup
  runtime <- x$fit$time()$total

  cli::cli_text("Chains: {.val {n_chains}}")
  cli::cli_text("Iterations: {.val {n_iter}} ({.val {n_warmup}} warmup)")
  cli::cli_text("Runtime: {prettyunits::pretty_sec(runtime)}")

  # Convergence diagnostics
  cli::cli_h3("Convergence")
  diag <- x$diagnostics

  rhat_ok <- !is.na(diag$max_rhat) && diag$max_rhat < 1.01
  ess_ok <- !is.na(diag$min_ess) && diag$min_ess > 400
  div_ok <- !is.na(diag$n_divergent) && diag$n_divergent == 0

  rhat_str <- sprintf("Max R-hat: %0.3f", diag$max_rhat)
  ess_str <- sprintf("Min ESS: %d", as.integer(diag$min_ess))
  div_str <- sprintf(
    "%d divergent transition%s",
    diag$n_divergent,
    ifelse(diag$n_divergent == 1, "", "s")
  )

  if (rhat_ok) {
    cli::cli_alert_success(rhat_str)
  } else {
    cli::cli_alert_warning(rhat_str)
  }
  if (ess_ok) {
    cli::cli_alert_success(ess_str)
  } else {
    cli::cli_alert_warning(ess_str)
  }
  if (div_ok) {
    cli::cli_alert_success(div_str)
  } else {
    cli::cli_alert_warning(div_str)
  }

  cli::cli_rule()
  cli::cli_text(c(
    "Use {.code summary()} for a detailed summary, ",
    "or {.code plot()} to visualize."
  ))

  invisible(x)
}
