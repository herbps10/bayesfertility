#' Plot prior vs. posterior for diagnostic comparison
#'
#' @param object a pk_fit object.
#'   Must have been fit with with \code{sample_prior = "yes"} or \code{sample_prior = "only"}.
#' @param type which quantity to compare: "coefficients" (individual beta_k[j]),
#'   "tau" (smoothness hyperparameters), "schedule_params" (constrained-scale
#'   schedule parameters at a cell), or "tfr" (TFR at a covariate grid)
#' @param parameter which schedule parameter(s) to include. NULL means all.
#' @param newdata for \code{type = "schedule_params"} or \code{type = "tfr"}: where
#'   to evaluate.
#' @param ndraws optional; subsample draws for speed.
#' @param ... reserved for future use.
#'
#' @return a ggplot object
#' @export
pk_prior_posterior <- function(
  object,
  type = c("coefficients", "tau", "schedule_params", "tfr"),
  parameter = NULL,
  terms = NULL,
  newdata = NULL,
  ndraws = NULL,
  ...
) {
  require_prior_samples(object)
  type <- match.arg(type)

  switch(
    type,
    coefficients = plot_pp_coefficients(object, parameter, terms, ndraws),
    schedule_params = plot_pp_schedule_params(
      object,
      newdata,
      parameter,
      ndraws
    ),
    tau = plot_pp_tau(object, parameter, ndraws),
    tfr = plot_pp_tfr(object, newdata, ndraws)
  )
}

#' @noRd
require_prior_samples <- function(object) {
  checkmate::assert_class(object, "pk_fit")
  if (is.null(object$prior_fit)) {
    cli::cli_abort(c(
      "No prior samples available in this {.cls pk_fit} object.",
      "i" = "Refit with {.code sample_prior = \"yes\"} or {.code \"only\"} to enable prior/posterior comparison."
    ))
  }

  invisible(object)
}

#' @noRd
pair_to_long <- function(prior_mat, post_mat, column_labels = NULL) {
  if (is.null(column_labels)) {
    column_labels <- if (!is.null(colnames(prior_mat))) {
      colnames(prior_mat)
    } else {
      paste0("V", seq_len(ncol(prior_mat)))
    }
  }

  to_long <- function(mat, source_label) {
    tibble::tibble(
      source = source_label,
      column = rep(column_labels, each = nrow(mat)),
      value = as.vector(mat)
    )
  }

  dplyr::bind_rows(
    to_long(prior_mat, "prior"),
    to_long(post_mat, "posterior")
  )
}

#' @noRd
plot_prior_posterior_density <- function(
  long_df,
  facet_var = "column",
  scales = "free",
  title = NULL
) {
  ggplot2::ggplot(long_df, ggplot2::aes(x = value, fill = source)) +
    ggplot2::geom_density(alpha = 0.4, color = NA) +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", facet_var)),
      scales = scales
    ) +
    ggplot2::scale_fill_manual(
      values = c(prior = "#9E9AC8", posterior = "#3F7EBA")
    ) +
    ggplot2::labs(title = title, x = "Value", y = "Density", fill = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}


#' @noRd
plot_pp_coefficients <- function(object, parameter, terms, ndraws) {
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")

  # Resolve which parameters to show
  if (is.null(parameter)) {
    parameters <- param_names
  } else {
    checkmate::assert_subset(parameter, param_names)
    parameters <- parameter
  }

  # Build a tibble of (parameter, term, coef_idx) for each coefficient
  rows <- list()
  for (p in parameters) {
    terms_info <- object$processed[[p]]$terms_info

    # Resolve which terms to show for this parameter
    avail_terms <- vapply(terms_info, function(t) t$label, character(1))
    if (is.null(terms)) {
      terms_to_show <- avail_terms
    } else {
      unknown <- setdiff(terms, avail_terms)
      if (length(unknown) > 0) {
        cli::cli_warn(c(
          "Some terms not found in {.field {p}}: {.val {unknown}}",
          "i" = "Available: {.val {avail_terms}}"
        ))
      }
      terms_to_show <- intersect(terms, avail_terms)
    }

    for (term_label in terms_to_show) {
      term <- Find(function(t) t$label == term_label, terms_info)
      for (j in seq_along(term$col_idx)) {
        coef_idx <- term$col_idx
        rows[[length(rows) + 1]] <- tibble::tibble(
          parameter = p,
          coef_idx = coef_idx[j],
          term = term_label,
          coef_label = if (length(coef_idx) == 1) {
            paste0(p, ": ", term_label)
          } else {
            paste0(p, ": ", term_label, "[", j, "]")
          }
        )
      }
    }
  }

  if (length(rows) == 0) {
    cli::cli_abort("No coefficients selected to plot.")
  }

  coef_map <- dplyr::bind_rows(rows)

  long_pieces <- lapply(seq_len(nrow(coef_map)), function(i) {
    p <- coef_map$parameter[i]
    idx <- coef_map$coef_idx[i]
    lbl <- coef_map$coef_label[i]

    prior_col <- object$prior_beta_draws[[p]][, idx, drop = TRUE]
    post_col <- object$beta_draws[[p]][, idx, drop = TRUE]

    if (!is.null(ndraws)) {
      prior_col <- prior_col[sample.int(
        length(prior_col),
        min(ndraws, length(prior_col))
      )]

      post_col <- post_col[sample.int(
        length(post_col),
        min(ndraws, length(post_col))
      )]
    }

    dplyr::bind_rows(
      tibble::tibble(coef_label = lbl, source = "prior", value = prior_col),
      tibble::tibble(coef_label = lbl, source = "posterior", value = post_col)
    )
  })

  long_df <- dplyr::bind_rows(long_pieces)

  long_df$coef_label <- factor(
    long_df$coef_label,
    levels = unique(coef_map$coef_label)
  )

  plot_prior_posterior_density(
    long_df,
    facet_var = "coef_label",
    title = "Prior vs. posterior: coefficients (working scale)"
  )
}


#' @noRd
plot_pp_schedule_params <- function(object, newdata, parameter, ndraws) {
  if (is.null(newdata)) {
    newdata <- unique_cells(object)
  }

  post_params <- predict(
    object,
    newdata,
    type = "schedule_params",
    source = "posterior",
    summarize = FALSE,
    ndraws = ndraws
  )
  prior_params <- predict(
    object,
    newdata,
    type = "schedule_params",
    source = "prior",
    summarize = FALSE,
    ndraws = ndraws
  )

  param_names <- if (is.null(parameter)) names(post_params) else parameter

  long_pieces <- lapply(param_names, function(p) {
    cell_labels <- paste0("cell_", seq_len(ncol(post_params[[p]])))

    dplyr::bind_rows(
      tibble::tibble(
        parameter = p,
        cell = rep(cell_labels, each = nrow(prior_params[[p]])),
        source = "prior",
        value = as.vector(prior_params[[p]])
      ),
      tibble::tibble(
        parameter = p,
        cell = rep(cell_labels, each = nrow(post_params[[p]])),
        source = "posterior",
        value = as.vector(post_params[[p]])
      )
    )
  })

  long_df <- dplyr::bind_rows(long_pieces)

  plot_prior_posterior_density(
    long_df,
    facet_var = "parameter",
    title = "Prior vs. posterior: schedule parameters (natural scale)"
  )
}

#' @noRd
plot_pp_tfr <- function(
  object,
  newdata,
  ndraws,
  age_range = c(15, 49),
  age_step = 0.5
) {
  # Default to unique training cells if newdata is not supplied
  if (is.null(newdata)) {
    newdata <- unique_cells(object)
  }

  checkmate::assert_data_frame(newdata)
  checkmate::assert_numeric(age_range, len = 2, finite = TRUE, sorted = TRUE)
  checkmate::assert_number(age_step, lower = 0, finite = TRUE)

  # Compute TFR draws from both sources
  post_tfr <- posterior_tfr(
    object,
    newdata = newdata,
    age_range = age_range,
    age_step = age_step,
    source = "posterior",
    summarize = FALSE,
    ndraws = ndraws
  )
  prior_tfr <- posterior_tfr(
    object,
    newdata = newdata,
    age_range = age_range,
    age_step = age_step,
    source = "prior",
    summarize = FALSE,
    ndraws = ndraws
  )

  cell_labels <- if (ncol(newdata) == 1) {
    as.character(newdata[[1]])
  } else {
    apply(newdata, 1, function(row) paste(row, collapse = ", "))
  }

  long_df <- dplyr::bind_rows(
    tibble::tibble(
      cell = rep(cell_labels, each = nrow(prior_tfr)),
      source = "prior",
      value = as.vector(prior_tfr)
    ),
    tibble::tibble(
      cell = rep(cell_labels, each = nrow(post_tfr)),
      source = "posterior",
      value = as.vector(post_tfr)
    )
  )

  long_df$cell <- factor(long_df$cell, levels = cell_labels)

  plot_prior_posterior_density(
    long_df,
    facet_var = "cell",
    title = "Prior vs. posterior: TFR by cell"
  )
}

#' @noRd
plot_pp_tau <- function(object, parameter, ndraws) {
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")

  if (is.null(parameter)) {
    parameters <- param_names
  } else {
    checkmate::assert_subset(parameter, param_names)
    parameters <- parameter
  }

  # Build a penalty -> label map, for parameters that have penalties
  rows <- list()
  for (p in parameters) {
    penalties <- object$processed[[p]]$penalties
    if (length(penalties) == 0) {
      next
    }

    for (k in seq_along(penalties)) {
      pen <- penalties[[k]]

      same_smooth_idx <- which(vapply(
        penalties,
        function(e) e$smooth_idx == pen$smooth_idx,
        logical(1)
      ))

      label <- if (length(same_smooth_idx) == 1) {
        paste0(p, ": ", pen$smooth_label)
      } else {
        within_smooth_idx <- which(same_smooth_idx == k)
        paste0(p, ": ", pen$smooth_label, "[", within_smooth_idx, "]")
      }

      rows[[length(rows) + 1]] <- tibble::tibble(
        parameter = p,
        penalty_idx = k,
        label = label
      )
    }
  }

  if (length(rows) == 0) {
    cli::cli_abort(c(
      "No penalties found for the selected parameters.",
      "i" = "Tau parameters only exist for terms with smooths or random effects."
    ))
  }

  pen_map <- dplyr::bind_rows(rows)

  # Build long-format tibble of prior and posterior draws.
  long_pieces <- lapply(seq_len(nrow(pen_map)), function(i) {
    p <- pen_map$parameter[i]
    k <- pen_map$penalty_idx[i]
    lbl <- pen_map$label[i]

    prior_col <- object$prior_tau_draws[[p]][, k, drop = TRUE]
    post_col <- object$tau_draws[[p]][, k, drop = TRUE]

    if (!is.null(ndraws)) {
      prior_col <- prior_col[sample.int(
        length(prior_col),
        min(ndraws, length(prior_col))
      )]

      post_col <- post_col[sample.int(
        length(post_col),
        min(ndraws, length(post_col))
      )]
    }

    dplyr::bind_rows(
      tibble::tibble(label = lbl, source = "prior", value = prior_col),
      tibble::tibble(label = lbl, source = "posterior", value = post_col)
    )
  })

  long_df <- dplyr::bind_rows(long_pieces)
  long_df$label <- factor(long_df$label, levels = pen_map$label)

  plot_prior_posterior_density(
    long_df,
    facet_var = "label",
    title = "Prior vs. posterior: smoothness hyperparameters (tau)"
  )
}
