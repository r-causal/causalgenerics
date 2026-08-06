#' Methods for a pooled inverse probability weighted result
#'
#' @description
#' The methods the class [pool_ipw()] returns carries.
#'
#' * `print()` summarizes the pooling and tabulates the pooled effects.
#' * `coef()` returns the pooled estimates.
#' * `vcov()` returns their pooled covariance.
#' * `confint()` returns their confidence limits.
#' * `nobs()` returns the smallest number of observations any pooled analysis
#'   was estimated from.
#' * `as.data.frame()` reports the pooled effects as a tidier-shaped table.
#'
#' @details
#' These live here for the reason the `ipw` methods do. Two packages each
#' registering `print.ipw_pooled()` would collide in the shared S3 method table,
#' and a caller writing against a pooled result would then get whichever package
#' was installed last rather than the contract.
#'
#' Everything they read is a field [pool_ipw()] wrote. They recompute no
#' pooling, so a result read back out of a saved file answers exactly as one
#' just built does.
#'
#' `ipw_pooled` deliberately does not inherit from `ipw`, and these methods are
#' deliberately separate from the ones registered against that class. The two
#' share column names and differ in what the numbers under them mean: a pooled
#' result reports several analyses rather than one fit and refers its inference
#' to t rather than to z, so a pooled result reaching `confint.ipw()` would come
#' back with normal limits and nothing would say so.
#'
#' # The degrees of freedom
#'
#' Every pooled effect carries its own degrees of freedom, in the `df` column of
#' the `estimates` frame. They are the Barnard-Rubin adjusted count, which
#' depends on how much of that effect's variance came from between the
#' imputations, so two effects pooled from the same analyses ordinarily differ
#' in them. `confint()` and `as.data.frame()` refer each row to t on its own
#' count rather than to the normal or to a single count for the table; with few
#' imputations those counts can be in single figures, where the difference is
#' large.
#'
#' There is deliberately no `df.residual()` method. Residual degrees of freedom
#' are a property of one fit, and a pooled result has as many as it had
#' imputations. The per-effect count above is not a residual count and differs
#' from row to row, so it is reported in the column beside the interval and the
#' p-value it was used for, and `df.residual()` on a pooled result finds no
#' method and no field and gives `NULL`.
#'
#' # The confidence limits
#'
#' `confint()` returns the limits the result stores for any row reported at the
#' level asked for, and rebuilds the rest. That is the rule `confint()` on an
#' unpooled result keeps, row by row: a stored pair need not be the one
#' recomputing gives, since it may have been rounded on its way into the frame.
#'
#' `as.data.frame()` applies the same rule to the frame as a whole, which is
#' what `as.data.frame()` on an unpooled result does. The stored pair comes back
#' only when every row records the level asked for, since a table mixing stored
#' limits with rebuilt ones would report two intervals under one pair of column
#' headings. `conf.level = NULL`, its default, names the level the frame records,
#' which is the level [pool_ipw()] built its limits at; a frame whose rows
#' disagree records none, and the limits are rebuilt at `0.95`.
#'
#' # Exponentiating
#'
#' On a marginal table `exponentiate = TRUE` means what it means for an unpooled
#' result. The rows labelled `log(rr)` and `log(or)` are matched exactly, their
#' estimate and their bounds move to the natural scale, and the two terms are
#' relabelled `rr` and `or`. The interval is settled before the scale is, so a
#' rebuilt bound is a t half width on the log scale added to an estimate on the
#' log scale and exponentiated afterwards. The standard error, the statistic,
#' and the p-value describe the log scale and stay there, and the `ipw_vcov`
#' attribute is dropped rather than carried, since it would describe neither the
#' table it sits on nor anything else.
#'
#' A conditional table has no rows labelled as ratios to pick out, so the link
#' the outcome models were fitted with settles the question for the whole table:
#' a `logit` link puts every coefficient on the log odds scale and a `log` link
#' puts every coefficient on the log risk scale, and both are scales an
#' exponential undoes. Every estimate moves and no term is relabelled, since the
#' terms are coefficient names and a coefficient does not change its name with
#' the scale its estimate is reported on. Every other link raises an error of
#' class `causalgenerics_exponentiate_link`, and of the classes
#' `causalgenerics_invalid_argument_exponentiate` and
#' `causalgenerics_invalid_argument`, rather than exponentiating coefficients
#' that describe nothing once exponentiated.
#'
#' @param x An `ipw_pooled` object.
#' @param object An `ipw_pooled` object.
#' @param parm The rows to report an interval for, given either as the effect
#'   labels or as their positions. Missing means all of them.
#' @param level The confidence level. At the level a row stores, that row's own
#'   limits are returned; at any other level they are rebuilt from t on the
#'   row's degrees of freedom.
#' @param row.names A character vector of row names for the returned table, or
#'   `NULL` for the automatic ones.
#' @param optional Accepted for the [base::as.data.frame()] generic. Every
#'   column of the table is named, so there is nothing for it to make optional.
#' @param conf.int If `TRUE`, append `conf.low` and `conf.high` columns after
#'   the rest of the table. Default is `FALSE`.
#' @param conf.level The confidence level the bounds report. `NULL`, the
#'   default, uses the level the result records.
#' @param exponentiate If `TRUE`, move the estimates that are on a log scale to
#'   their natural scale, as the section above describes. Default is `FALSE`.
#' @param ... Further arguments. These methods ignore them.
#'
#' @return
#' `print()` returns its input invisibly.
#'
#' `coef()` returns a named numeric vector of pooled estimates, named by effect
#' label in the marginal reading and by coefficient name in the conditional one.
#'
#' `vcov()` returns the pooled covariance of those estimates, with the same
#' labels as dimnames on both margins. A result whose analyses did not all carry
#' a covariance records none, and raises an error of class
#' `causalgenerics_no_vcov_ipw_pooled`, and of the general class
#' `causalgenerics_no_vcov`, rather than returning one built from the standard
#' errors, which would report the effects as uncorrelated.
#'
#' `confint()` returns a matrix with one row per effect `parm` selects, in the
#' order `parm` gives them, and two columns holding the lower and upper limit,
#' named by effect label and by the two tail probabilities as percentages. A
#' character `parm` naming an effect the result does not report raises an error
#' of class `causalgenerics_invalid_argument`.
#'
#' `nobs()` returns a single integer.
#'
#' `as.data.frame()` returns a plain data frame with the columns `term`, then
#' `contrast` when the result names contrasts, then `estimate`, `std.error`,
#' `statistic`, `df`, and `p.value`, with `conf.low` and `conf.high` appended
#' when they are asked for. The pooled covariance travels on it under the
#' `ipw_vcov` attribute unless the table was exponentiated.
#'
#' @seealso [pool_ipw()], which produces these results, and [new_ipw()] for the
#'   unpooled result they are pooled from.
#'
#' @examples
#' dat <- data.frame(
#'   x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
#'   z = rep(c(0, 1), 10),
#'   y = rep(c(0, 1, 1, 0, 1), 4)
#' )
#'
#' imputed <- lapply(1:3, function(i) {
#'   completed <- dat
#'   completed$y[i] <- 1 - completed$y[i]
#'   completed
#' })
#'
#' estimate <- list(c(0.20, 0.55), c(0.30, 0.60), c(0.40, 0.65))
#' std_err <- list(c(0.10, 0.25), c(0.20, 0.30), c(0.30, 0.35))
#'
#' fits <- Map(
#'   function(completed, estimate, std.err) {
#'     new_ipw(
#'       estimand = "ate",
#'       wt_mod = glm(z ~ x, family = binomial(), data = completed),
#'       outcome_mod = glm(y ~ z, family = quasibinomial(), data = completed),
#'       estimates = data.frame(
#'         effect = c("rd", "log(rr)"),
#'         estimate = estimate,
#'         std.err = std.err,
#'         z = estimate / std.err,
#'         ci.lower = estimate - 1.96 * std.err,
#'         ci.upper = estimate + 1.96 * std.err,
#'         conf.level = 0.95,
#'         p.value = 2 * pnorm(-abs(estimate / std.err))
#'       ),
#'       se_method = "linearization",
#'       fit = NULL
#'     )
#'   },
#'   imputed,
#'   estimate,
#'   std_err
#' )
#'
#' pooled <- pool_ipw(fits)
#'
#' pooled
#'
#' coef(pooled)
#'
#' # Wider than the normal limits would be, since the degrees of freedom are
#' # what few imputations leave.
#' confint(pooled)
#'
#' nobs(pooled)
#'
#' # The tidier-shaped table, with the degrees of freedom the statistic beside
#' # it is referred to.
#' as.data.frame(pooled)
#'
#' # With an interval, and the ratio on its natural scale.
#' as.data.frame(pooled, conf.int = TRUE, exponentiate = TRUE)
#'
#' # Residual degrees of freedom belong to one fit, so a pooled result has no
#' # method for them; the per-effect count is in the table above.
#' df.residual(pooled)
#'
#' @name ipw-pooled-methods
NULL

#' @rdname ipw-pooled-methods
#' @export
print.ipw_pooled <- function(x, ...) {
  # The mode is read first, so that a result recording something that is not a
  # reading is refused before any of the summary reaches the console.
  effects <- ipw_effects(x)

  cat("Pooled Inverse Probability Weight Estimator\n")
  cat("Estimand:", toupper(x$estimand), "\n")
  cat("Effects:", ipw_effects_label(effects), "\n")
  # The two facts an unpooled result has no counterpart for. How many analyses
  # went in decides how much of each interval is between-imputation variance,
  # and the complete-data count is what the small-sample adjustment spent, so a
  # reader cannot judge the degrees of freedom in the table without both.
  cat("Imputations:", x$m, "\n")
  cat("Complete-data df:", format(x$dfcom), "\n")
  cat("\n")

  # The reading names itself in the heading as well as above, since the two
  # readings are different tables of different numbers and a reader handed the
  # table alone would otherwise take a pooled `(Intercept)` for a causal effect.
  if (effects == "conditional") {
    cat("Pooled conditional estimates (outcome model):\n")
  } else {
    cat("Pooled marginal estimates:\n")
  }
  print_effect_table(x$estimates)

  # `riv`, `lambda`, and `fmi` are deliberately absent. They describe how much
  # of the uncertainty came from the imputation rather than from the data, which
  # is a diagnostic rather than a result; they are in the `pooling` field for a
  # caller who wants them; and five more columns beside the estimates would push
  # the table past the width where `printCoefmat()` keeps a row on one line.

  invisible(x)
}

#' @rdname ipw-pooled-methods
#' @export
#' @importFrom stats coef setNames
coef.ipw_pooled <- function(object, ...) {
  # The same labels the unpooled result names its estimates with, read through
  # the same helper, so a caller who pools a set of results gets the names they
  # were reading before they pooled them. In the conditional reading the frame's
  # `effect` column holds the coefficient names and there is no contrast column,
  # so the helper gives those back unchanged.
  stats::setNames(
    object$estimates$estimate,
    ipw_effect_labels(object$estimates)
  )
}

#' @rdname ipw-pooled-methods
#' @export
#' @importFrom stats vcov
vcov.ipw_pooled <- function(object, ...) {
  # `pool_ipw()` attaches no covariance when any of the results it pooled
  # carried none, since the average of the within-imputation covariances needs
  # one from every imputation. There is nothing to fall back on: the standard
  # errors in the frame give the diagonal and say nothing about the off-diagonal
  # entries, which are far from zero for effects estimated from the same
  # weighted means.
  covariance <- attr(object$estimates, "ipw_vcov", exact = TRUE)
  if (is.null(covariance)) {
    stop_no_vcov("ipw_pooled")
  }
  covariance
}

#' @rdname ipw-pooled-methods
#' @export
#' @importFrom stats confint qt
confint.ipw_pooled <- function(object, parm, level = 0.95, ...) {
  estimates <- object$estimates
  labels <- ipw_effect_labels(estimates)
  rows <- if (missing(parm)) {
    seq_along(labels)
  } else {
    select_effects(parm, labels)
  }

  # One half width, added and subtracted, taken at each row's own pooled degrees
  # of freedom. `qt()` is not exactly antisymmetric, so taking the lower limit
  # from the lower tail instead would put the two bounds a bit apart from each
  # other, and a single count for the table would be the wrong one for every row
  # but at most one.
  half_width <- stats::qt(1 - (1 - level) / 2, estimates$df) * estimates$std.err
  lower <- estimates$estimate - half_width
  upper <- estimates$estimate + half_width

  # At the level a row was reported for, the reported limits are returned rather
  # than rebuilt, row by row, which is the rule `confint()` on an unpooled
  # result keeps. A stored pair need not be the one recomputing gives: one
  # rounded on its way into the frame is not that number. `which()` rather than
  # a logical index so that a missing `conf.level` selects nothing instead of
  # raising a subscript error.
  stored <- which(estimates$conf.level == level)
  lower[stored] <- estimates$ci.lower[stored]
  upper[stored] <- estimates$ci.upper[stored]

  limits <- cbind(lower[rows], upper[rows])
  dimnames(limits) <- list(labels[rows], percent_labels(level))
  limits
}

#' @rdname ipw-pooled-methods
#' @export
#' @importFrom stats nobs
nobs.ipw_pooled <- function(object, ...) {
  # The field rather than a computation. `nobs.default()` would find a list
  # element of this name and answer with it, which is the right number by
  # coincidence of two names agreeing; a registered method is what makes it the
  # contract, and what keeps the answer an integer.
  as.integer(object$nobs)
}

#' @rdname ipw-pooled-methods
#' @export
#' @importFrom stats qt
as.data.frame.ipw_pooled <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  conf.int = FALSE,
  conf.level = NULL,
  exponentiate = FALSE
) {
  # Checked on every call rather than on the branch that reads them, for the
  # reason `as.data.frame()` on an unpooled result checks them all: a level no
  # interval can be built at is wrong whichever way `conf.int` was set.
  check_flag(conf.int, "conf.int")
  if (!is.null(conf.level)) {
    check_conf_level(conf.level)
  }
  check_flag(exponentiate, "exponentiate")

  estimates <- x$estimates
  effects <- ipw_effects(x)

  # Before anything is built, so a table that cannot be reported on the scale
  # asked for is refused rather than half-built.
  if (exponentiate && effects == "conditional") {
    check_exponentiate_link(x$outcome_link)
  }

  # `NULL` names the level the frame records. A fixed default would rebuild the
  # limits of a result pooled at any other level, which is a table of numbers
  # the result already holds, computed a second time and reported as though it
  # had not been.
  level <- if (is.null(conf.level)) {
    pool_stored_level(estimates)
  } else {
    conf.level
  }
  if (is.null(level)) {
    level <- 0.95
  }

  # `term` first, then the column naming the contrast it qualifies when the
  # result reports one. `df` sits after the statistic, since that is what the
  # statistic beside it is referred to and a table without it leaves a reader
  # nothing to refer it to.
  contrast <- ipw_contrast_column(estimates)
  columns <- list(term = as.character(estimates$effect))
  if (!is.null(contrast)) {
    columns$contrast <- estimates[[contrast]]
  }
  columns$estimate <- estimates$estimate
  columns$std.error <- estimates$std.err
  columns$statistic <- estimates$t
  columns$df <- estimates$df
  columns$p.value <- estimates$p.value

  # Built before the scale is changed, so that a rebuilt bound is a half width
  # on the log scale added to an estimate on the log scale.
  bounds <- if (conf.int) pooled_interval_bounds(estimates, level) else NULL

  if (exponentiate) {
    # The marginal reading picks its rows out by label, matched exactly so that
    # a table whose ratios are already on the natural scale is left alone rather
    # than exponentiated a second time. The conditional reading has no such
    # labels: the link checked above settles it for every row at once.
    is_log_rr <- columns$term == "log(rr)"
    is_log_or <- columns$term == "log(or)"
    ratios <- if (effects == "conditional") {
      rep(TRUE, length(columns$term))
    } else {
      is_log_rr | is_log_or
    }

    # Only the point estimate and the bounds move. `std.error`, `statistic`, and
    # `p.value` stay on the log scale, which is where the inference is done.
    columns$estimate[ratios] <- exp(columns$estimate[ratios])
    if (!is.null(bounds)) {
      bounds$lower[ratios] <- exp(bounds$lower[ratios])
      bounds$upper[ratios] <- exp(bounds$upper[ratios])
    }

    # The label names the scale, so it moves with the value it labels. A
    # coefficient name does not: it names the term rather than the scale, and a
    # conditional table relabelled here would report a coefficient the outcome
    # model never had.
    if (effects != "conditional") {
      columns$term[is_log_rr] <- "rr"
      columns$term[is_log_or] <- "or"
    }
  }

  # Last, after the columns the table always carries, so that asking for an
  # interval adds to the table rather than rearranging it.
  if (!is.null(bounds)) {
    columns$conf.low <- bounds$lower
    columns$conf.high <- bounds$upper
  }

  df <- data.frame(columns, row.names = row.names, stringsAsFactors = FALSE)

  # The covariance belongs to the estimates rather than to a column of them, so
  # it travels on the table under the attribute the contract names it under. A
  # result that carries none carries none here: assigning `NULL` sets no
  # attribute, which is the right answer rather than an accident.
  if (!exponentiate) {
    attr(df, "ipw_vcov") <- attr(estimates, "ipw_vcov", exact = TRUE)
  }

  df
}

#' The confidence bounds `as.data.frame()` reports for a pooled result
#'
#' The frame-level rule `interval_bounds()` keeps for an unpooled one, with the
#' limits rebuilt from t on each row's own pooled degrees of freedom rather than
#' from the normal. The stored pair comes back only when every row of the frame
#' was reported at the level asked for, since a table mixing stored bounds with
#' rebuilt ones would report two intervals under one pair of column headings.
#'
#' This is deliberately not the row-by-row rule `confint()` keeps. A matrix of
#' limits is read a row at a time and says which level it is at in its column
#' names; a tidier-shaped table is read as a table, and the level is an argument
#' to the call rather than a column of it, so a table whose rows came from two
#' different rules would have nothing on it to say so.
#'
#' @param estimates The `estimates` component of an `ipw_pooled` object.
#' @param conf.level The level the bounds report.
#'
#' @return A list of two numeric vectors, `lower` and `upper`.
#'
#' @noRd
#' @importFrom stats qt
pooled_interval_bounds <- function(estimates, conf.level) {
  stored <- estimates$conf.level
  if (!is.null(stored) && isTRUE(all(stored == conf.level))) {
    return(list(lower = estimates$ci.lower, upper = estimates$ci.upper))
  }

  half_width <- stats::qt(1 - (1 - conf.level) / 2, estimates$df) *
    estimates$std.err
  list(
    lower = estimates$estimate - half_width,
    upper = estimates$estimate + half_width
  )
}

#' Refuse to exponentiate coefficients that are not on a log scale
#'
#' The conditional reading reports the outcome model's coefficients, and there
#' are no rows labelled as ratios among them to pick out. The link the models
#' were fitted with is what says whether there is anything for an exponential to
#' undo: a logit link puts every coefficient on the log odds scale and a log
#' link puts every coefficient on the log risk scale, and a coefficient on any
#' other scale exponentiates to a number describing nothing.
#'
#' There is no subset of rows to move instead, so the call is refused rather
#' than answered in part.
#'
#' @param link The result's `outcome_link`.
#' @param call The call to report the error against, which is the method's
#'   rather than this helper's.
#'
#' @return `link`, invisibly, when it is one an exponential undoes.
#'
#' @noRd
check_exponentiate_link <- function(link, call = sys.call(-1)) {
  exponentiable <- c("logit", "log")

  if (!link %in% exponentiable) {
    stop_exponentiate_link(link, exponentiable, call = call)
  }

  invisible(link)
}
