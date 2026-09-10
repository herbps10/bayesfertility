#' bayesfertility: Fit Bayesian fertility schedule models
#' @name bayesfertility
#' @description Fit Bayesian fertility schedule models, including the Peristera-Kostaki model
#' @family help
#' @importFrom instantiate stan_package_model
#' @importFrom stats predict setNames
#' @importFrom utils head
NULL

`%||%` <- function(a, b) if (is.null(a)) b else a

# Column names referenced via NSE (ggplot2::aes(), tibble::tibble()'s `:=`)
# that R CMD check's static analysis can't otherwise resolve.
utils::globalVariables(c(".fitted", ".lower", ".upper", "age", "asfr", "value", ":="))
