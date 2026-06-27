#' Specify per-parameter formulas and prior customization for a PK model
#'
#' @param c1,c2,mu1,mu2,sigma1,sigma2 one-sided formulas specifying the
#'   linear predictor for each Peristera-Kostaki schedule parameter on its
#'   working scale (log for c1/c2/sigma1/sigma2, identify for mu1,
#'   log of the gap mu2-mu1 for mu2).
#' @param prior_scales optional named list of prior scale overrides, keyed
#'   first by schedule parameter (e.g., "c1") and then by term name
#'   (e.g., "(Intercept)", "race", "s(province)"). Use `NULL` for a term
#'   to indicate "use the penalty-based prior only" (only meaningful for
#'   smooth/random-effect terms). Unspecified terms fall back to
#'   package defaults.
#' @param prior_centers optional named list of prior center overrides,
#'   same structure as `prior_scales`. Non-intercept terms default to 0;
#'   intercepts default to demographically anchored values.
#' @param tau_priors optional named list of smoothness-hyperparameter
#'   (tau) prior overrides, keyed first by schedule parameter (e.g., "c1")
#'   and then by smooth label (e.g., "s(province)"). Each entry is
#'   itself a list of parameters; currently only `scale` is supported,
#'   controlling the scale of the half-Student-t(3, 0, scale) prior
#'   on tau. Unspecified smooths fall back to per-parameter defaults.
#' @export
pk_effects <- function(
  c1 = ~1,
  c2 = ~1,
  mu1 = ~1,
  mu2 = ~1,
  sigma1 = ~1,
  sigma2 = ~1,
  prior_scales = NULL,
  prior_centers = NULL,
  tau_priors = NULL
) {
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")

  formulas <- list(
    c1 = c1,
    c2 = c2,
    mu1 = mu1,
    mu2 = mu2,
    sigma1 = sigma1,
    sigma2 = sigma2
  )

  for (p in param_names) {
    checkmate::assert_formula(formulas[[p]], .var.name = p)
    if (length(formulas[[p]]) != 2L) {
      cli::cli_abort(c(
        "Formula for {.field {p}} must be one-sided.",
        "x" = "You supplied {.code {format(formulas[[p]])}}."
      ))
    }
  }

  prior_scales <- validate_prior_spec(
    prior_scales,
    param_names,
    what = "prior_scales"
  )
  prior_centers <- validate_prior_spec(
    prior_centers,
    param_names,
    what = "prior_centers"
  )
  tau_priors <- validate_tau_priors(
    tau_priors,
    param_names
  )

  structure(
    c(
      formulas,
      list(
        prior_scales = prior_scales,
        prior_centers = prior_centers,
        tau_priors = tau_priors
      )
    ),
    class = "pk_effects"
  )
}

#' Validate the structure of prior_scalecs or prior_centers
#' @noRd
validate_prior_spec <- function(spec, param_names, what) {
  # Return empty lists so downstream code can rely on consistent shape
  if (is.null(spec)) {
    out <- stats::setNames(
      replicate(length(param_names), list(), simplify = FALSE),
      param_names
    )
    return(out)
  }

  checkmate::assert_list(spec, names = "named", .var.name = what)

  # Top-level names must be valid schedule parameters
  unknown <- setdiff(names(spec), param_names)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown schedule parameter(s) in {.arg {what}}: {.val {unknown}}",
      "i" = "Valid parameters are: {.val {param_names}}"
    ))
  }

  # Each nested entry must be a named list of NULL or numeric values
  for (p in names(spec)) {
    entry <- spec[[p]]
    checkmate::assert_list(
      entry,
      names = "named",
      .var.name = paste0(what, "$", p)
    )

    for (term_name in names(entry)) {
      val <- entry[[term_name]]

      if (is.null(val)) {
        next
      }

      if (what == "prior_scales") {
        checkmate::assert_number(
          val,
          lower = 0,
          finite = TRUE,
          .var.name = paste0(what, "$", p, "$`", term_name, "`")
        )
      } else {
        checkmate::assert_number(
          val,
          finite = TRUE,
          .var.name = paste0(what, "$", p, "$`", term_name, "`")
        )
      }
    }
  }

  for (p in param_names) {
    if (is.null(spec[[p]])) spec[[p]] <- list()
  }

  spec[param_names]
}

#' Validate the structure of tau_priors
#'
#' @noRd
validate_tau_priors <- function(spec, param_names) {
  supported_params <- c("scale")

  if (is.null(spec)) {
    out <- stats::setNames(
      replicate(length(param_names), list(), simplify = FALSE),
      param_names
    )
    return(out)
  }

  checkmate::assert_list(spec, names = "named", .var.name = "tau_priors")

  unknown <- setdiff(names(spec), param_names)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown schedule parameter(s) in {.arg tau_priors}: {.val {unknown}}",
      "i" = "Valid parameters are: {.val {param_names}}"
    ))
  }

  for (p in names(spec)) {
    entry <- spec[[p]]
    checkmate::assert_list(
      entry,
      names = "named",
      .var.name = paste0("tau_priors$", p)
    )

    for (smooth_label in names(entry)) {
      smooth_spec <- entry[[smooth_label]]
      var_name <- paste0("tau_priors$", p, "$`", smooth_label, "`")

      checkmate::assert_list(
        smooth_spec,
        names = "named",
        .var.name = var_name
      )

      unknown_params <- setdiff(names(smooth_spec), supported_params)
      if (length(unknown_params) > 0) {
        cli::cli_abort(c(
          "Unknown parameter(s) in {.code {var_name}}: {.val {unknown_params}}",
          "i" = "Currently supported parameters are: {.val {supported_params}}"
        ))
      }

      if ("scale" %in% names(smooth_spec)) {
        checkmate::assert_number(
          smooth_spec$scale,
          lower = 0,
          finite = TRUE,
          .var.name = paste0(var_name, "$scale")
        )
      }
    }
  }

  for (p in param_names) {
    if (is.null(spec[[p]])) spec[[p]] <- list()
  }

  spec[param_names]
}


#' Print a pk_effects object
#' @export
print.pk_effects <- function(x, ...) {
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")

  cli::cli_h3("PK effects specification")

  max_name_width <- max(nchar(param_names))
  for (p in param_names) {
    formula_str <- deparse(x[[p]], width.cutoff = 500L)
    formula_str <- str_replace(formula_str, "^~\\s*", "")
    padded <- format(p, width = max_name_width)
    cli::cli_li("{.field {padded}} ~ {formula_str}")
  }

  has_scale_overrides <- any(lengths(x$prior_scales) > 0)
  has_center_overrides <- any(lengths(x$prior_centers) > 0)
  has_tau_overrides <- any(lengths(x$tau_priors) > 0)

  if (has_scale_overrides || has_center_overrides || has_tau_overrides) {
    cli::cli_h3("Prior customizations")
    for (p in param_names) {
      scales <- x$prior_scales[[p]]
      centers <- x$prior_centers[[p]]
      taus <- x$tau_priors[[p]]

      if (length(scales) == 0 && length(centers) == 0 && length(taus) == 0) {
        next
      }

      cli::cli_text("{.field {p}}:")
      for (term in union(names(scales), names(centers))) {
        parts <- character(0)
        if (term %in% names(scales)) {
          val <- scales[[term]]
          parts <- c(
            parts,
            paste0("scale = ", if (is.null(val)) "NULL (penalty only)" else val)
          )
        }
        if (term %in% names(centers)) {
          val <- centers[[term]]
          parts <- c(
            parts,
            paste0("scale = ", if (is.null(val)) "NULL (penalty only)" else val)
          )
        }

        cli::cli_text("  {.code {term}}: {paste(parts, collapse = ', ')}")
      }

      for (smooth_label in names(taus)) {
        tau_spec <- taus[[smooth_label]]
        parts <- character(0)
        if ("scale" %in% names(tau_spec)) {
          parts <- c(parts, paste0("tau scale = ", tau_spec$scale))
        }
        cli::cli_text(
          " {.code {smooth_label}} (tau): {paste(parts, collapse = ', ')}"
        )
      }
    }
  }
  invisible(x)
}
