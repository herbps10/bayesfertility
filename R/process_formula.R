#' Process formulas in a pk_effects specification
#'
#' @param effects output of pk_effects()
#' @param data data frame
#'
#' @return a named list (one entry per schedule parameter) of
#'   `process_formula()` output, each augmented with resolved
#'   `prior_scales`, `prior_centers` (from `resolve_priors()`), and
#'   `tau_scales` (from `resolve_tau_priors()`)
#'
#' @noRd
process_effects <- function(effects, data) {
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")

  processed <- lapply(param_names, function(p) {
    tryCatch(
      process_formula(effects[[p]], data),
      error = function(e) {
        stop(sprintf(
          "Error processing formula for parameter '%s': %s",
          p,
          conditionMessage(e)
        ))
      }
    )
  })
  names(processed) <- param_names

  user_scales <- effects$prior_scales %||%
    setNames(replicate(6, list(), simplify = FALSE), param_names)
  user_centers <- effects$prior_centers %||%
    setNames(replicate(6, list(), simplify = FALSE), param_names)
  user_tau_priors <- effects$tau_priors %||%
    setNames(replicate(6, list(), simplify = FALSE), param_names)

  for (p in param_names) {
    resolved <- resolve_priors(
      terms_info = processed[[p]]$terms_info,
      user_scales = user_scales[[p]] %||% list(),
      user_centers = user_centers[[p]] %||% list(),
      param_name = p
    )
    processed[[p]]$prior_scales <- resolved$scales
    processed[[p]]$prior_centers <- resolved$centers

    # Per-penalty tau priors
    processed[[p]]$tau_scales <- resolve_tau_priors(
      terms_info = processed[[p]]$terms_info,
      n_penalties = length(processed[[p]]$penalties),
      user_tau_priors = user_tau_priors[[p]] %||% list(),
      param_name = p
    )
  }

  processed
}

#' Build a list of parametric terms from a model matrix and formula
#'
#' @param X model matrix (with "assign" attribute) or design matrix from jagam
#' @param terms_obj the `terms()`/`pterms` object the columns of `X` came from
#'
#' @return named list: term label -> integer vector of column indices in X
#' @noRd
extract_parametric_terms <- function(X, terms_obj) {
  asgn <- attr(X, "assign")

  if (is.null(asgn)) {
    return(list("(Intercept)" = seq_len(ncol(X))))
  }

  labels <- attr(terms_obj, "term.labels")
  has_intercept <- attr(terms_obj, "intercept") == 1L

  # Map each column to its term label
  col_term <- character(length(asgn))
  for (i in seq_along(asgn)) {
    if (asgn[i] == 0L) {
      col_term[i] <- "(Intercept)"
    } else {
      col_term[i] <- labels[asgn[i]]
    }
  }

  # Group column indices by term
  split(seq_along(asgn), factor(col_term, levels = unique(col_term)))
}


#' Resolve prior scales and centers per term, given user overrides and defaults
#'
#' @param terms_info a parameter's `terms_info` list (from `build_terms_info()`)
#' @param user_scales named list of user-supplied per-term scale overrides
#'   for this schedule parameter (as in `pk_effects(prior_scales = ...)`)
#' @param user_centers named list of user-supplied per-term center overrides
#'   for this schedule parameter
#' @param param_name one of "c1", "c2", "mu1", "mu2", "sigma1", "sigma2"
#'
#' @return a list with elements `scales` and `centers`, each a numeric vector
#'   with one entry per element of `terms_info` (`NA` for smooth/random-effect
#'   terms left to use the penalty-based prior only)
#' @noRd
resolve_priors <- function(terms_info, user_scales, user_centers, param_name) {
  available_terms <- vapply(terms_info, function(t) t$label, character(1))

  unknown_scales <- setdiff(names(user_scales), available_terms)
  if (length(unknown_scales) > 0) {
    cli::cli_abort(c(
      "Unknown term name(s) in {.arg prior_scales} for {.field {param_name}}: {.val {unknown_scales}}",
      "i" = "Available terms: {.val {available_terms}}"
    ))
  }

  unknown_centers <- setdiff(names(user_centers), available_terms)
  if (length(unknown_scales) > 0) {
    cli::cli_abort(c(
      "Unknown term name(s) in {.arg prior_centers} for {.field {param_name}}: {.val {unknown_centers}}",
      "i" = "Available terms: {.val {available_terms}}"
    ))
  }

  defaults_s <- pk_default_prior_scales(param_name)
  defaults_c <- pk_default_prior_centers(param_name)

  scales <- numeric(length(terms_info))
  centers <- numeric(length(terms_info))

  for (i in seq_along(terms_info)) {
    t <- terms_info[[i]]
    lbl <- t$label

    # Resolve scale
    if (lbl %in% names(user_scales)) {
      val <- user_scales[[lbl]]
      scales[i] <- if (is.null(val)) NA_real_ else val
    } else if (lbl == "(Intercept)") {
      scales[i] <- defaults_s[["(Intercept)"]] %||% NA_real_
    } else if (t$type == "smooth") {
      scales[i] <- defaults_s[["smooth"]] %||% NA_real_
    } else {
      scales[i] <- defaults_s[["parametric"]] %||% NA_real_
    }

    # Resolve center
    if (lbl %in% names(user_centers)) {
      val <- user_centers[[lbl]]
      centers[i] <- if (is.null(val)) NA_real_ else val
    } else if (lbl == "(Intercept)") {
      centers[i] <- defaults_c[["(Intercept)"]] %||% NA_real_
    } else if (t$type == "smooth") {
      centers[i] <- defaults_c[["smooth"]] %||% NA_real_
    } else {
      centers[i] <- defaults_c[["parametric"]] %||% NA_real_
    }
  }

  list(scales = scales, centers = centers)
}

#' Resolve per-penalty tau prior scales from user overrides and defaults
#'
#' @param terms_info a parameter's `terms_info` list (from `build_terms_info()`)
#' @param n_penalties number of penalty terms for this schedule parameter
#' @param user_tau_priors named list of user-supplied per-smooth tau overrides
#'   for this schedule parameter (as in `pk_effects(tau_priors = ...)`)
#' @param param_name one of "c1", "c2", "mu1", "mu2", "sigma1", "sigma2"
#'
#' @return a numeric vector of length `n_penalties` giving the tau prior
#'   scale for each penalty
#' @noRd
resolve_tau_priors <- function(
  terms_info,
  n_penalties,
  user_tau_priors,
  param_name
) {
  # Validate user-supplied smooth labels against actual smooth terms
  smooth_labels_available <- vapply(
    Filter(function(t) t$type == "smooth", terms_info),
    function(t) t$label,
    character(1)
  )

  unknown <- setdiff(names(user_tau_priors), smooth_labels_available)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown smooth label(s) in {.arg tau_priors} for {.field {param_name}}: {.val {unknown}}",
      "i" = "Available smooth labels: {.val {smooth_labels_available}}"
    ))
  }

  default_scale <- pk_default_tau_scales(param_name)

  resolved_scales <- rep(default_scale, n_penalties)

  for (smooth_label in names(user_tau_priors)) {
    user_spec <- user_tau_priors[[smooth_label]]

    if ("scale" %in% names(user_spec)) {
      term <- Find(
        function(t) t$type == "smooth" && t$label == smooth_label,
        terms_info
      )
      penalty_idx <- term$penalty_idx

      if (length(penalty_idx) > 0) {
        resolved_scales[penalty_idx] <- user_spec$scale
      }
    }
  }

  resolved_scales
}


#' Process a single formula into design + penalty matrices
#'
#' @param formula one-side formula
#' @param data data frame
#'
#' @return a list with
#'   - X: design matrix (n_rows x n_coef)
#'   - S_list: list of penalty matrices, one per smooth term
#'   - S_col_idx: list of integer vectors, indicating the columns of
#'                X each penalty matrix applies to
#'   - n_coef: number of coefficients
#'   - n_smooths: number of smooth terms
#'   - smooth_labels: character vector naming each smooth
#'
#' @noRd
process_formula <- function(formula, data) {
  has_smooth <- formula_has_smooth(formula)

  if (!has_smooth) {
    X <- stats::model.matrix(formula, data = data)
    n_coef <- ncol(X)
    parametric_terms <- extract_parametric_terms(X, stats::terms(formula))

    terms_info <- build_terms_info(
      parametric_terms = parametric_terms,
      smooth_labels = character(0),
      smooth_col_idx = list(),
      smooth_penalty_idx = list()
    )

    return(list(
      X = X,
      penalties = list(),
      n_coef = n_coef,
      n_smooths = 0L,
      parametric_idx = unlist(parametric_terms, use.names = FALSE),
      terms_info = terms_info,
      pregam = NULL
    ))
  }

  # jagam needs a left-hand side, so we create a dummy outcome
  dummy_y <- rep(0, nrow(data))
  data_with_y <- data
  data_with_y[["..dummy_y.."]] <- dummy_y
  full_formula <- stats::update(formula, ..dummy_y.. ~ .)

  jagam_file <- tempfile(fileext = ".jags")
  on.exit(unlink(jagam_file), add = TRUE)

  jd <- mgcv::jagam(
    formula = full_formula,
    data = data_with_y,
    file = jagam_file,
    family = stats::gaussian(),
    diagonalize = FALSE
  )

  X <- jd$jags.data$X
  n_coef <- ncol(X)

  smooths <- jd$pregam$smooth
  n_smooths <- length(smooths)

  # Extract penalties
  penalties <- list()
  smooth_col_idx <- list()
  smooth_penalty_idx <- list()
  smooth_labels <- character(n_smooths)

  if (n_smooths > 0) {
    for (k in seq_len(n_smooths)) {
      sm <- smooths[[k]]

      # Coefficient indices for this smooth
      idx <- sm$first.para:sm$last.para

      smooth_col_idx[[k]] <- idx
      smooth_labels[k] <- sm$label
      pen_indices_for_smooth <- integer(0)

      if (!is.null(sm$S) && length(sm$S) > 0) {
        for (p in seq_along(sm$S)) {
          is_re <- inherits(sm, "random.effect")
          penalties[[length(penalties) + 1]] <- list(
            S = sm$S[[p]],
            col_idx = idx,
            smooth_idx = k,
            penalty_idx = p,
            smooth_label = sm$label,
            is_re = is_re
          )
          pen_indices_for_smooth <- c(pen_indices_for_smooth, length(penalties))
        }
      }
      smooth_penalty_idx[[k]] <- pen_indices_for_smooth
    }
  }

  # Extract parametric terms using pregam$pterms
  if (!is.null(jd$pregam$pterms)) {
    if (n_smooths > 0) {
      first_smooth_col <- min(vapply(
        smooths,
        function(sm) as.integer(sm$first.para),
        integer(1)
      ))
      n_param_cols <- first_smooth_col - 1L
    } else {
      n_param_cols <- n_coef
    }

    if (n_param_cols > 0) {
      X_param <- X[, seq_len(n_param_cols), drop = FALSE]
      xlevels_p <- jd$pregam$xlevels %||%
        stats::.getXlevels(jd$pregam$pterms, data_with_y)
      X_param_with_assign <- stats::model.matrix(
        jd$pregam$pterms,
        data = data_with_y,
        xlev = xlevels_p
      )
      parametric_terms <- extract_parametric_terms(
        X_param_with_assign,
        jd$pregam$pterms
      )
    } else {
      parametric_terms <- list()
    }
  } else {
    parametric_terms <- list()
  }

  terms_info <- build_terms_info(
    parametric_terms = parametric_terms,
    smooth_labels = smooth_labels,
    smooth_col_idx = smooth_col_idx,
    smooth_penalty_idx = smooth_penalty_idx
  )

  parametric_idx <- if (length(parametric_terms) > 0) {
    unlist(parametric_terms, use.names = FALSE)
  } else {
    integer(0)
  }

  list(
    X = X,
    penalties = penalties,
    n_coef = n_coef,
    n_smooths = n_smooths,
    parametric_idx = parametric_idx,
    terms_info = terms_info,
    pregam = jd$pregam
  )
}

#' Build a unified terms_info list combining parametric and smooth terms
#'
#' @param parametric_terms named list from `extract_parametric_terms()`:
#'   term label -> integer vector of column indices
#' @param smooth_labels character vector of smooth term labels
#' @param smooth_col_idx list of integer vectors, one per smooth, giving the
#'   columns of `X` belonging to that smooth
#' @param smooth_penalty_idx list of integer vectors, one per smooth, giving
#'   indices into the flat `penalties` list for that smooth
#'
#' @return a list with one element per logical term. Each element is a list:
#'   - label:       character (e.g., "(Intercept)", "race", "s(province)")
#'   - type:        "parametric" or "smooth"
#'   - col_idx:     integer vector of columns in X belonging to this term
#'   - penalty_idx: integer vector of indices into the flat `penalties` list
#'                  (empty for parametric terms)
#' @noRd
build_terms_info <- function(
  parametric_terms,
  smooth_labels,
  smooth_col_idx,
  smooth_penalty_idx
) {
  out <- list()

  for (lbl in names(parametric_terms)) {
    out[[length(out) + 1]] <- list(
      label = lbl,
      type = "parametric",
      col_idx = parametric_terms[[lbl]],
      penalty_idx = integer(0)
    )
  }

  for (k in seq_along(smooth_labels)) {
    out[[length(out) + 1]] <- list(
      label = smooth_labels[[k]],
      type = "smooth",
      col_idx = smooth_col_idx[[k]],
      penalty_idx = smooth_penalty_idx[[k]]
    )
  }

  out
}


#' Check whether a formula contains an mgcv-style smooth term
#'
#' @param formula a one-sided formula
#'
#' @return `TRUE` if any term label starts with `s(`, `te(`, `ti(`, or `t2(`
#' @noRd
formula_has_smooth <- function(formula) {
  term_labels <- attr(stats::terms(formula), "term.labels")
  if (length(term_labels) == 0) {
    return(FALSE)
  }

  # Check if any term starts with s(, te(, ti(, or t2(
  any(grepl("^(s|te|ti|t2)\\(", term_labels))
}
