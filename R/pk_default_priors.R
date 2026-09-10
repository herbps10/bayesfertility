#' Default prior centers for a schedule parameter's terms
#'
#' @param param_name one of "c1", "c2", "mu1", "mu2", "sigma1", "sigma2"
#'
#' @return a named list with elements `"(Intercept)"` (a demographically
#'   anchored default for `param_name`), `parametric` (default for
#'   non-intercept parametric terms, `0`), and `smooth` (default for
#'   smooth/random-effect terms, `0`)
#' @noRd
pk_default_prior_centers <- function(param_name) {
  list(
    "(Intercept)" = switch(
      param_name,
      c1 = -3.0,
      #c1 = -1.5,
      c2 = -2.0,
      mu1 = 20.0,
      mu2 = log(8),
      sigma1 = log(3.5),
      sigma2 = log(9)
    ),
    parametric = 0,
    smooth = 0
  )
}

#' Default prior scales for a schedule parameter's terms
#'
#' @param param_name one of "c1", "c2", "mu1", "mu2", "sigma1", "sigma2"
#'
#' @return a named list with elements `"(Intercept)"` and `parametric`
#'   (default Normal prior SDs for `param_name`) and `smooth` (`NULL`,
#'   meaning smooth/random-effect terms use the penalty-based prior only)
#' @noRd
pk_default_prior_scales <- function(param_name) {
  list(
    "(Intercept)" = switch(
      param_name,
      c1 = 0.5,
      #c1 = 1,
      c2 = 0.5,
      mu1 = 0.7,
      mu2 = 0.2,
      sigma1 = 0.15,
      sigma2 = 0.4
    ),
    parametric = switch(
      param_name,
      c1 = 1.0,
      c2 = 1.0,
      mu1 = 0.5,
      mu2 = 0.5,
      sigma1 = 0.3,
      sigma2 = 0.3
    ),
    smooth = NULL
  )
}

#' Default smoothness-hyperparameter (tau) prior scale for a schedule parameter
#'
#' @param param_name one of "c1", "c2", "mu1", "mu2", "sigma1", "sigma2"
#'
#' @return a single number: the default scale of the half-Student-t(3, 0,
#'   scale) prior on `tau` for `param_name`'s smooth/random-effect terms
#' @noRd
pk_default_tau_scales <- function(param_name) {
  switch(
    param_name,
    c1 = 0.5,
    c2 = 0.5,
    mu1 = 0.2,
    mu2 = 0.3,
    sigma1 = 0.2,
    sigma2 = 0.2
  )
}
