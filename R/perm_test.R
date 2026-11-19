# Covariate-adjusted blocked residual permutation tests ----------------

#' Permutation tests for regression coefficients
#'
#' @description
#' A tidier implementation of the Freedman--Lane residual permutation procedure
#' for testing one or more coefficients of a Gaussian linear model while
#' adjusting for covariates and, optionally, restricting permutations to
#' within exchangeability blocks. Unlike permuting the response directly, this
#' procedure preserves the covariate structure of the data and supplies a
#' joint null distribution for all of the tested coefficients from a single
#' set of permutations.
#'
#' @param x A data frame containing all of the variables referenced by
#'   `full_formula`, `reduced_formula`, and `blocks`.
#' @param full_formula A two-sided [formula][stats::formula] giving the
#'   complete Gaussian linear model, including both the covariates that are
#'   retained under the null hypothesis and the terms whose coefficients are
#'   being tested.
#' @param reduced_formula A two-sided [formula][stats::formula] giving the
#'   model nested within `full_formula` that contains only the terms retained
#'   under the null hypothesis. The coefficients that distinguish
#'   `full_formula` from `reduced_formula` are the coefficients under test.
#' @param terms A character vector naming the coefficients (model-matrix
#'   column names) under test, in the order they should be reported. Must be
#'   exactly the set of coefficients in `full_formula` that are not in
#'   `reduced_formula`.
#' @param reps The number of random permutations to draw. Must be a single
#'   positive integer.
#' @param blocks An optional bare (unquoted) column name in `x`---or,
#'   equivalently, a single string naming such a column---giving the
#'   exchangeability blocks. Residuals are only permuted within blocks. When
#'   `NULL` (the default), the whole data set is treated as one block, the
#'   classical Freedman--Lane procedure.
#' @param direction A character string giving the direction of the test.
#'   Options are `"two-sided"` (default), `"greater"`, or `"less"`. Aliases
#'   `"left"`, `"right"`, `"both"`, `"two_sided"`, `"two sided"`, and
#'   `"two.sided"` are also accepted.
#'
#' @details
#' The procedure, for each of `reps` random permutations, is:
#'
#' 1. Fit the reduced model implied by `reduced_formula` and store its fitted
#'    values and residuals.
#' 2. Draw one permutation of the observation indices that only reorders
#'    observations within each level of `blocks` (a single global
#'    permutation shared by every tested coefficient) and build a permuted
#'    response as the reduced-model fitted values plus the permuted
#'    reduced-model residuals.
#' 3. Fit `full_formula` to the permuted response using infer's native
#'    Gaussian linear model fitter and record the studentized (t) statistics
#'    for the tested coefficients.
#'
#' The observed statistics are the studentized statistics from fitting
#' `full_formula` to the original response. Because the permuted response is
#' a real numeric vector that is passed to the native fitter, every fit
#' carries the same model-matrix structure as the observed fit, so the
#' coefficient names, formula, and sample correspondence are preserved.
#'
#' Raw p-values use the finite-sample correction that counts the observed
#' arrangement as one of the permutations:
#' \deqn{p_k = (1 + \#\{b: T_{bk} \ge T_k\}) / (reps + 1)}
#' with absolute statistics for two-sided tests. The step-down maxT adjusted
#' p-values (Westfall and Young, 1993) reuse the same permutations: at each
#' step the running maximum of the (absolute, for two-sided tests) null
#' statistics across all remaining coefficients supplies the joint null,
#' enforcing multiplicity control without a Bonferroni correction. The tail
#' comparisons are inclusive, so equality between a null and an observed
#' statistic counts against the null.
#'
#' Every permutation must produce a valid, full-rank fit with finite
#' studentized statistics; a failure in any permutation raises an error
#' rather than silently dropping the permutation.
#'
#' @return A [tibble][tibble::tibble] with one row per tested coefficient and
#' the columns:
#'
#' \itemize{
#'   \item `term`: The model-matrix column name of the tested coefficient.
#'   \item `estimate`: The observed coefficient estimate from the full model.
#'   \item `statistic`: The observed studentized (t) statistic for the
#'     coefficient from the full model.
#'   \item `p_value`: The raw, unadjusted permutation p-value.
#'   \item `p_value_adj`: The step-down maxT adjusted p-value.
#' }
#'
#' The number of permutations and the test direction are stored as the
#' `reps` and `direction` attributes of the result.
#'
#' @references
#' Freedman, D., & Lane, D. (1983). A Nonstochastic Interpretation of
#' Reported Significance Levels. *Journal of Business & Economic Statistics*,
#' 1(4), 292--298. \doi{10.1080/07350015.1983.10509354}
#'
#' Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
#' Testing: Examples and Methods for p-Value Adjustment*. Wiley.
#'
#' @examples
#' # test the partial effects of age and college on hours worked,
#' # adjusting for the intercept (the only term in the reduced model)
#' set.seed(1)
#' perm_test(
#'   gss,
#'   full_formula = hours ~ age + college,
#'   reduced_formula = hours ~ 1,
#'   terms = c("age", "collegedegree"),
#'   reps = 100
#' )
#'
#' # adjust for a covariate that is retained under the null and only
#' # permute residuals within levels of a blocking variable
#' set.seed(1)
#' perm_test(
#'   gss,
#'   full_formula = hours ~ age + college + sex,
#'   reduced_formula = hours ~ sex,
#'   terms = c("age", "collegedegree"),
#'   reps = 100,
#'   blocks = sex
#' )
#'
#' # one-sided test
#' set.seed(1)
#' perm_test(
#'   gss,
#'   full_formula = hours ~ age,
#'   reduced_formula = hours ~ 1,
#'   terms = "age",
#'   reps = 100,
#'   direction = "greater"
#' )
#'
#' @family wrapper functions
#' @export
perm_test <- function(
  x,
  full_formula,
  reduced_formula,
  terms,
  reps,
  blocks = NULL,
  direction = "two-sided"
) {
  call <- caller_env()

  check_type(x, is.data.frame, call = call)
  x <- standardize_variable_types(x)

  check_perm_formulas(full_formula, reduced_formula, call = call)
  check_perm_terms(terms, call = call)
  check_perm_reps(reps, call = call)
  check_direction(direction, call = call)
  direction <- norm_direction(direction)

  blocks_name <- check_perm_blocks(x, {{ blocks }}, call = call)

  data <- check_perm_data(
    x, full_formula, reduced_formula, terms, blocks_name, call = call
  )

  # Fixed design matrices: their structure, and thus coefficient names and
  # the correspondence of observations, never changes across permutations.
  full_mf <- stats::model.frame(full_formula, data = data, na.action = stats::na.fail)
  reduced_mf <- stats::model.frame(
    reduced_formula, data = data, na.action = stats::na.fail
  )
  X_full <- stats::model.matrix(full_formula, full_mf)
  X_reduced <- stats::model.matrix(reduced_formula, reduced_mf)
  y <- stats::model.response(full_mf)

  check_perm_design(
    X_full, X_reduced, y, full_formula, reduced_formula, terms, call = call
  )

  # The observed full-model fit, via infer's native Gaussian fitter.
  observed <- fit_gaussian_model(data, full_formula, call = call)
  term_idx <- match(terms, observed$coef$term)
  obs_t <- observed$coef$statistic[term_idx]

  # The reduced model supplies the fixed fitted component and the residuals
  # that are permuted within blocks.
  reduced_fit <- fit_gaussian_model(data, reduced_formula, call = call)
  if (any(!is.finite(reduced_fit$resid))) {
    cli_abort(
      "The reduced model produced non-finite residuals, so a valid \\
       permutation response cannot be constructed.",
      call = call
    )
  }

  block_id <- if (is.null(blocks_name)) {
    rep(1L, nrow(data))
  } else {
    as.integer(as.factor(data[[blocks_name]]))
  }
  block_rows <- split(seq_len(nrow(data)), block_id)
  if (!any(lengths(block_rows) > 1L)) {
    cli_abort(
      "Every block contains a single observation, so no residuals can be \\
       permuted. Supply a blocking variable with at least one block \\
       containing multiple observations.",
      call = call
    )
  }

  # One shared within-block permutation per replicate for every coefficient.
  null_t <- permutation_null_stats(
    data = data,
    response_name = as.character(f_lhs(full_formula)),
    full_formula = full_formula,
    fitted = reduced_fit$fitted,
    resid = reduced_fit$resid,
    block_rows = block_rows,
    term_idx = term_idx,
    terms = terms,
    reps = reps,
    call = call
  )

  raw_p <- raw_permutation_p(null_t, obs_t, direction, reps)
  adj_p <- stepdown_maxt_p(null_t, obs_t, direction, reps)

  res <- tibble::tibble(
    term = terms,
    estimate = observed$coef$estimate[term_idx],
    statistic = obs_t,
    p_value = raw_p,
    p_value_adj = adj_p
  )

  attr(res, "reps") <- reps
  attr(res, "direction") <- direction

  res
}

# Validation helpers -----------------------------------------------------

check_perm_formulas <- function(full_formula, reduced_formula, call) {
  if (!rlang::is_formula(full_formula) || length(full_formula) != 3L) {
    cli_abort(
      "`full_formula` must be a two-sided formula, e.g. `y ~ x + z`.",
      call = call
    )
  }
  if (!rlang::is_formula(reduced_formula) || length(reduced_formula) != 3L) {
    cli_abort(
      "`reduced_formula` must be a two-sided formula, e.g. `y ~ z`.",
      call = call
    )
  }

  full_response <- as.character(f_lhs(full_formula))
  reduced_response <- as.character(f_lhs(reduced_formula))
  if (!identical(full_response, reduced_response)) {
    cli_abort(
      "`full_formula` and `reduced_formula` must have the same response \\
       variable, but got {.val {full_response}} and {.val {reduced_response}}.",
      call = call
    )
  }

  if (length(attr(stats::terms(full_formula), "offset")) ||
        length(attr(stats::terms(reduced_formula), "offset"))) {
    cli_abort(
      "Offset terms are not supported by `perm_test()`.",
      call = call
    )
  }

  invisible(TRUE)
}

check_perm_terms <- function(terms, call) {
  check_type(terms, is.character, call = call)
  if (length(terms) == 0L || any(!nzchar(terms)) || anyNA(terms)) {
    cli_abort(
      "`terms` must be a non-empty character vector of coefficient names.",
      call = call
    )
  }
  if (anyDuplicated(terms)) {
    cli_abort("`terms` must not contain duplicate coefficient names.", call = call)
  }

  invisible(TRUE)
}

check_perm_reps <- function(reps, call) {
  if (!is.numeric(reps) || length(reps) != 1L || !is.finite(reps) ||
        reps < 1 || !rlang::is_scalar_integerish(reps)) {
    cli_abort(
      "`reps` must be a single positive integer, not {.obj_type_friendly {reps}}.",
      call = call
    )
  }

  invisible(TRUE)
}

check_perm_blocks <- function(x, blocks, call) {
  blocks_expr <- rlang::enquo(blocks)

  if (rlang::quo_is_null(blocks_expr)) {
    return(NULL)
  }

  blocks_expr_obj <- rlang::get_expr(blocks_expr)
  if (rlang::is_symbol(blocks_expr_obj)) {
    # A bare column name: use the name without evaluating it.
    blocks_name <- rlang::as_name(blocks_expr_obj)
  } else {
    blocks_value <- rlang::eval_tidy(blocks_expr, data = x)
    if (is.null(blocks_value)) {
      return(NULL)
    }
    if (rlang::is_string(blocks_value)) {
      blocks_name <- blocks_value
    } else {
      cli_abort(
        "`blocks` must be a single bare column name or string naming a \\
         column in `x`.",
        call = call
      )
    }
  }

  if (!blocks_name %in% names(x)) {
    cli_abort(
      "The blocking variable {.val {blocks_name}} is not a column in `x`.",
      call = call
    )
  }

  blocks_name
}

check_perm_data <- function(
  x, full_formula, reduced_formula, terms, blocks_name, call
) {
  full_vars <- all.vars(full_formula)
  reduced_vars <- all.vars(reduced_formula)
  used_vars <- unique(c(full_vars, reduced_vars, blocks_name))
  missing_vars <- setdiff(used_vars, names(x))
  if (length(missing_vars) > 0L) {
    cli_abort(
      "{qty(missing_vars)}The variable{?s} {.val {missing_vars}} referenced \\
       by the formulas or `blocks` {?is/are} not in `x`.",
      call = call
    )
  }

  data <- tibble::as_tibble(x)[used_vars]

  if (anyNA(data)) {
    cli_abort(
      "`perm_test()` does not support missing values; {sum(is.na(data))} \\
       missing value{?s} {?was/were} found in the variables used by the \\
       test. Complete or remove the relevant rows before calling the \\
       function.",
      call = call
    )
  }

  non_finite_cols <- names(data)[
    vapply(data, function(col) is.numeric(col) && any(!is.finite(col)), logical(1))
  ]
  if (length(non_finite_cols) > 0L) {
    cli_abort(
      "{qty(non_finite_cols)}The variable{?s} {.val {non_finite_cols}} \\
       contain{?s/} non-finite values (`NaN` or `Inf`), which are not valid \\
       for a Gaussian linear model.",
      call = call
    )
  }

  response_name <- full_vars[[1]]
  response <- data[[response_name]]
  if (!is.numeric(response)) {
    cli_abort(
      "Gaussian linear models require a numeric response, but \\
       {.val {response_name}} is {.obj_type_friendly {response}}.",
      call = call
    )
  }
  if (any(!is.finite(response))) {
    cli_abort(
      "The response {.val {response_name}} contains non-finite values.",
      call = call
    )
  }

  invisible(data)
}

check_perm_design <- function(
  X_full, X_reduced, y, full_formula, reduced_formula, terms, call
) {
  n <- length(y)
  if (nrow(X_full) != n || nrow(X_reduced) != n) {
    cli_abort(
      "The full and reduced models must be fit to the same observations.",
      call = call
    )
  }

  rank_full <- qr(X_full)$rank
  rank_reduced <- qr(X_reduced)$rank
  if (rank_full < ncol(X_full)) {
    cli_abort(
      "The full model's design matrix is rank deficient (rank {rank_full} \\
       with {ncol(X_full)} columns), so the tested coefficients are not \\
       uniquely defined. Remove redundant terms before calling \\
       `perm_test()`.",
      call = call
    )
  }
  if (rank_reduced < ncol(X_reduced)) {
    cli_abort(
      "The reduced model's design matrix is rank deficient (rank \\
       {rank_reduced} with {ncol(X_reduced)} columns). Remove redundant \\
       terms before calling `perm_test()`.",
      call = call
    )
  }

  # The reduced design must span a subspace of the full design.
  nesting_resid <- stats::lm.fit(X_full, X_reduced)$residuals
  if (max(abs(nesting_resid)) > sqrt(.Machine$double.eps)) {
    cli_abort(
      "`reduced_formula` must be nested within `full_formula`: every term \\
       in the reduced model must also appear in the full model.",
      call = call
    )
  }

  full_names <- colnames(X_full)
  reduced_names <- colnames(X_reduced)
  added <- setdiff(full_names, reduced_names)

  missing_terms <- setdiff(terms, full_names)
  if (length(missing_terms) > 0L) {
    cli_abort(
      "{qty(missing_terms)}The coefficient{?s} {.val {missing_terms}} \\
       supplied to `terms` {?is/are} not in the full model. The full \\
       model's coefficients are {.val {full_names}}.",
      call = call
    )
  }

  extra_terms <- setdiff(terms, added)
  if (length(extra_terms) > 0L) {
    cli_abort(
      "{qty(extra_terms)}The coefficient{?s} {.val {extra_terms}} \\
       supplied to `terms` {?is/are} also in the reduced model and so \\
       {?is/are} not under test.",
      call = call
    )
  }

  untested <- setdiff(added, terms)
  if (length(untested) > 0L) {
    cli_abort(
      "`terms` must name every coefficient that distinguishes the full \\
       model from the reduced model. Also include {.val {untested}}.",
      call = call
    )
  }

  invisible(TRUE)
}

# Model fitting and permutation ------------------------------------------

# Fit a Gaussian GLM with infer's native fitter and return the studentized
# statistics alongside the fitted values and residuals needed by the
# Freedman--Lane procedure.
fit_gaussian_model <- function(data, formula, call) {
  fit <- tryCatch(
    stats::glm(formula = formula, data = data, family = stats::gaussian()),
    error = function(cnd) {
      cli_abort(
        "Fitting {.code {deparse1(formula)}} failed.",
        parent = cnd,
        call = call
      )
    }
  )

  if (!fit$converged) {
    cli_abort(
      "Fitting {.code {deparse1(formula)}} did not converge.",
      call = call
    )
  }

  tidy <- broom::tidy(fit)
  if (any(!is.finite(tidy$statistic)) || any(!is.finite(tidy$estimate))) {
    cli_abort(
      "Fitting {.code {deparse1(formula)}} produced non-finite estimates \\
       or statistics.",
      call = call
    )
  }

  list(
    coef = tibble::tibble(
      term = tidy$term,
      estimate = tidy$estimate,
      statistic = tidy$statistic
    ),
    fitted = stats::fitted(fit),
    resid = stats::residuals(fit)
  )
}

permutation_null_stats <- function(
  data,
  response_name,
  full_formula,
  fitted,
  resid,
  block_rows,
  term_idx,
  terms,
  reps,
  call
) {
  n <- length(resid)
  null_t <- matrix(NA_real_, nrow = reps, ncol = length(terms))
  colnames(null_t) <- terms

  for (b in seq_len(reps)) {
    perm <- integer(n)
    for (rows in block_rows) {
      if (length(rows) == 1L) {
        perm[rows] <- rows
      } else {
        perm[rows] <- sample(rows, size = length(rows), replace = FALSE)
      }
    }

    # Real, numeric response: fitted component plus permuted residuals.
    perm_data <- data
    perm_data[[response_name]] <- fitted + resid[perm]

    perm_fit <- tryCatch(
      fit_gaussian_model(perm_data, full_formula, call = call),
      error = function(cnd) {
        cli_abort(
          "Permutation {b} of {reps} failed; permutation samples cannot be \\
           discarded. Check the data and block structure and rerun.",
          parent = cnd,
          call = call
        )
      }
    )

    if (!identical(perm_fit$coef$term[term_idx], terms)) {
      cli_abort(
        "Permutation {b} of {reps} produced different coefficient names \\
         than the observed fit, so the joint null distribution is not \\
         well-defined.",
        call = call
      )
    }

    null_t[b, ] <- perm_fit$coef$statistic[term_idx]
  }

  if (anyNA(null_t) || any(!is.finite(null_t))) {
    cli_abort(
      "At least one permutation produced non-finite studentized statistics. \\
       Permutation samples cannot be discarded; resolve the degenerate data \\
       and rerun.",
      call = call
    )
  }

  null_t
}

# P-value calculations ---------------------------------------------------

# Transform to the statistic whose large values are evidence against the
# null: absolute value for two-sided tests, the raw or negated statistic
# for the directional tests.
directed_stats <- function(stat, direction) {
  switch(
    direction,
    both = abs(stat),
    right = stat,
    left = -stat
  )
}

# Finite-sample-corrected raw p-values: the observed arrangement counts as
# one of the permutations and the tail is inclusive of equality.
raw_permutation_p <- function(null_t, obs_t, direction, reps) {
  null_abs <- directed_stats(null_t, direction)
  obs_abs <- directed_stats(obs_t, direction)

  vapply(
    seq_along(obs_abs),
    function(j) {
      (sum(null_abs[, j] >= obs_abs[j]) + 1) / (reps + 1)
    },
    numeric(1)
  )
}

# Westfall--Young free-combinations step-down maxT adjustment on a single
# shared set of permutations. Enforced to be monotone as the observed
# statistics decrease in magnitude.
stepdown_maxt_p <- function(null_t, obs_t, direction, reps) {
  null_abs <- directed_stats(null_t, direction)
  obs_abs <- directed_stats(obs_t, direction)

  running_max <- apply(null_abs, 1L, max)
  order_desc <- order(obs_abs, decreasing = TRUE)

  adj <- numeric(length(obs_abs))
  previous <- 0
  for (j in order_desc) {
    p <- (sum(running_max >= obs_abs[j]) + 1) / (reps + 1)
    previous <- max(previous, p)
    adj[j] <- previous
  }

  adj
}
