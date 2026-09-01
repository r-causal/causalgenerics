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
#' * `estimand()` returns the estimand the pooled analyses targeted.
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
#' # The reading these methods report
#'
#' A pooled result carries both readings whenever [pool_ipw()] could compute
#' both: the one it presents and, under `alternate`, the one it does not.
#' `coef()`, `vcov()`, `confint()`, and `as.data.frame()` report the reading the
#' result presents, and take an `effects` argument naming one for a single call.
#' Naming a reading reports it and leaves the result as it was; [as_marginal()]
#' and [as_conditional()] are what change which reading a result presents.
#'
#' The reading is settled before anything else the call asks for. `confint()`
#' selects the rows `parm` names from the reading it was asked for, and
#' `as.data.frame()` decides what `exponentiate = TRUE` moves and what it
#' relabels from that reading too, so a conditional table is gated on the outcome
#' model link whichever reading the result itself records. A reading the pooling
#' could not compute is refused rather than answered with the one the result
#' holds, and a value naming neither reading raises an error of class
#' `causalgenerics_invalid_argument_effects`, and of the general class
#' `causalgenerics_invalid_argument`, as it does everywhere else in this package.
#'
#' `print()` reports the reading the result presents and takes no such argument.
#' It names that reading in its heading, and a caller who wants the other one
#' moves the result to it first. `nobs()` and `estimand()` describe the pooled
#' analyses rather than a reading of them, so they answer the same way in either
#' reading and take no `effects` argument.
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
#' `confint()` and `as.data.frame()` return the limits the result stores when
#' every row of the frame records the level asked for, and rebuild every row
#' otherwise. That is the rule both keep on an unpooled result, and it is a rule
#' about the frame rather than about a row: a stored pair need not be the one
#' recomputing gives, since it may have been rounded on its way into the frame,
#' so a surface mixing stored limits with rebuilt ones would report two kinds of
#' interval under one pair of column headings with nothing on it to say which is
#' which.
#'
#' `conf.level = NULL`, the `as.data.frame()` default, names the level the frame
#' records, which is the level [pool_ipw()] built its limits at; a frame whose
#' rows disagree records none, and the limits are rebuilt at `0.95`.
#'
#' # Exponentiating
#'
#' On a marginal table `exponentiate = TRUE` means what it means for an unpooled
#' result. The rows labeled `log(rr)` and `log(or)` are matched exactly, their
#' estimate and their bounds move to the natural scale, and the two terms are
#' relabeled `rr` and `or`. The interval is settled before the scale is, so a
#' rebuilt bound is a t half width on the log scale added to an estimate on the
#' log scale and exponentiated afterwards. The standard error, the statistic,
#' and the p-value describe the log scale and stay there, and the `ipw_vcov`
#' attribute is dropped rather than carried, since it would describe neither the
#' table it sits on nor anything else.
#'
#' A conditional table has no rows labeled as ratios to pick out, so the link
#' the outcome models were fitted with settles the question for the whole table:
#' a `logit` link puts every coefficient on the log odds scale and a `log` link
#' puts every coefficient on the log risk scale, and both are scales an
#' exponential undoes. Every estimate moves and no term is relabeled, since the
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
#' @param level The confidence level. At the level every row of the result
#'   stores, the stored limits are returned; at any other level, and for a frame
#'   whose rows do not agree on one, every row is rebuilt from t on its own
#'   degrees of freedom.
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
#' @param effects The reading to report, either `"marginal"` or
#'   `"conditional"`. `NULL`, the default, reports the reading the result
#'   records; any other value overrides it for the one call and leaves the
#'   result as it is. A reading the pooling could not compute is refused with an
#'   error of class `causalgenerics_pool_missing_surface`, which is where
#'   [as_marginal()] and [as_conditional()] refuse it, and a result pooled
#'   before both readings were kept refuses a request for the other reading with
#'   the same classes, having nothing recorded to report.
#' @param call The call to report a refusal against. A tidier that builds its
#'   table by calling `as.data.frame()` passes the call a user wrote, so the
#'   refusal names the function they typed rather than the delegation behind it.
#'   The default is the method's own call, which is what a direct call reports.
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
#' `contrast` when the result names contrasts, then `group` when it names
#' subgroups, then `estimate`, `std.error`, `statistic`, `df`, and `p.value`,
#' with `conf.low` and `conf.high` appended when they are asked for. The pooled
#' covariance travels on it under the `ipw_vcov` attribute unless the table was
#' exponentiated.
#'
#' `estimand()` returns the estimand the pooled results agreed on, which is the
#' one the weights their estimates were computed under targeted.
#'
#' @seealso [pool_ipw()], which produces these results, [as_marginal()] and
#'   [as_conditional()], which move one between the readings it carries, and
#'   [new_ipw()] for the unpooled result they are pooled from.
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
coef.ipw_pooled <- function(object, ..., effects = NULL) {
  object <- pooled_surface(object, effects)

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
vcov.ipw_pooled <- function(object, ..., effects = NULL) {
  object <- pooled_surface(object, effects)

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
confint.ipw_pooled <- function(
  object,
  parm,
  level = 0.95,
  ...,
  effects = NULL
) {
  # Before `parm` is matched against anything, so the rows are selected from the
  # reading that was asked for rather than from the one the result records.
  object <- pooled_surface(object, effects)

  estimates <- object$estimates
  labels <- ipw_effect_labels(estimates)
  rows <- if (missing(parm)) {
    seq_along(labels)
  } else {
    select_effects(parm, labels)
  }

  # The limits the pooled table reports, built by the helper that table is built
  # by, which is the rule `confint()` on an unpooled result keeps: the stored
  # pair comes back only when every row of the frame records the level asked
  # for, and every row is rebuilt from t on its own pooled degrees of freedom
  # otherwise. The level names itself in the column headings and nowhere else,
  # so a matrix mixing stored limits with rebuilt ones would report two kinds of
  # interval under one heading.
  bounds <- pooled_interval_bounds(estimates, level)

  limits <- cbind(bounds$lower[rows], bounds$upper[rows])
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
  exponentiate = FALSE,
  effects = NULL,
  call = sys.call()
) {
  # Checked on every call rather than on the branch that reads them, for the
  # reason `as.data.frame()` on an unpooled result checks them all: a level no
  # interval can be built at is wrong whichever way `conf.int` was set. Each of
  # them reports against `call`, so a tidier that builds its table here names
  # one entry point across every refusal on the path.
  check_flag(conf.int, "conf.int", call = call)
  if (!is.null(conf.level)) {
    check_conf_level(conf.level, call = call)
  }
  check_flag(exponentiate, "exponentiate", call = call)

  # The reading is settled before anything is read off the result, so the gate
  # below and every column built after it belong to the reading that was asked
  # for. Once the result presents that reading, the reading is what it records.
  x <- pooled_surface(x, effects, call = call)

  estimates <- x$estimates
  effects <- ipw_effects(x, call = call)

  # Before anything is built, so a table that cannot be reported on the scale
  # asked for is refused rather than half-built.
  if (exponentiate && effects == "conditional") {
    check_exponentiate_link(x$outcome_link, call = call)
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

  # The columns that name a row come from the helper the unpooled table is keyed
  # by, so a caller who pools a set of results reads the same table under the
  # same headings. `df` sits after the statistic, since that is what the
  # statistic beside it is referred to and a table without it leaves a reader
  # nothing to refer it to.
  columns <- ipw_term_columns(estimates)
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
    # conditional table relabeled here would report a coefficient the outcome
    # model never had.
    if (effects != "conditional") {
      columns$term[is_log_rr] <- "rr"
      columns$term[is_log_or] <- "or"
    }
  }

  ipw_tidy_frame(
    columns,
    bounds,
    # The covariance describes the estimates on the scale they were estimated
    # on, so an exponentiated table carries none.
    covariance = if (!exponentiate) {
      attr(estimates, "ipw_vcov", exact = TRUE)
    } else {
      NULL
    },
    row.names = row.names
  )
}

#' @rdname ipw-pooled-methods
#' @export
estimand.ipw_pooled <- function(x, ...) {
  # The field rather than a computation, the way `estimand.ipw()` reads the
  # unpooled one. The estimand is one of the things the results had to agree on
  # before they could be pooled, so a pooled result records exactly one; it
  # describes the analyses rather than a reading of them, which is why this
  # method takes no `effects` argument and answers the same way either way.
  #
  # A method of its own rather than the unpooled one registered a second time,
  # since `ipw_pooled` does not inherit from `ipw` and a shared method would
  # carry a later change to the unpooled reading into the pooled result with
  # nothing to say so.
  x$estimand
}

#' The reading a pooled accessor reports
#'
#' An accessor that takes an `effects` argument reports the reading the caller
#' names, and the reading the result records when the caller names none. A
#' pooled result holds both readings, so the accessor reads the reading it was
#' asked for off the result the flip gives rather than reaching into the
#' `alternate` field itself.
#'
#' Routing it through the flip is what keeps one answer to the question of what
#' a reading a pooling could not compute does: the flip raises the recorded
#' reason, and every accessor raises it in the same words without repeating the
#' branch that decides it. The flip answers a request for the reading the result
#' already presents with the result itself, so declining to name one costs
#' nothing and reports what it reported before the argument existed.
#'
#' @param object An `ipw_pooled` object.
#' @param effects The `effects` argument as the caller supplied it, or `NULL`.
#' @param call The call to report the error against, which is the accessor's
#'   rather than this helper's.
#'
#' @return An `ipw_pooled` object presenting the reading asked for.
#'
#' @noRd
pooled_surface <- function(object, effects, call = sys.call(-1)) {
  flip_ipw_pooled(
    object,
    resolve_ipw_effects(object, effects, call = call),
    call = call
  )
}

#' The confidence bounds a pooled result reports
#'
#' The frame-level rule `interval_bounds()` keeps for an unpooled one, with the
#' limits rebuilt from t on each row's own pooled degrees of freedom rather than
#' from the normal. The stored pair comes back only when every row of the frame
#' was reported at the level asked for, since a surface mixing stored bounds
#' with rebuilt ones would report two kinds of interval under one pair of column
#' headings.
#'
#' Both surfaces the pooled result presents read the bounds through this helper,
#' as the unpooled ones read theirs through `interval_bounds()`. The level is
#' named in the column headings of a matrix and in the argument that built a
#' table, and neither says which rows were stored and which rebuilt, so the two
#' report the same numbers for the same result at the same level.
#'
#' @param estimates The `estimates` component of an `ipw_pooled` object.
#' @param conf.level The level the bounds report.
#'
#' @return A list of two numeric vectors, `lower` and `upper`.
#'
#' @noRd
#' @importFrom stats qt
pooled_interval_bounds <- function(estimates, conf.level) {
  # By exact name, for the reason `interval_bounds()` reads it that way: the
  # column is optional, and a frame whose `$` is stricter than a plain data
  # frame's warns about a column the contract never required.
  stored <- estimates[["conf.level"]]
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
#' are no rows labeled as ratios among them to pick out. The link the models
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
