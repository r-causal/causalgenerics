#' Construct an inverse probability weighted result
#'
#' @description
#' `new_ipw()` is the low-level constructor for the object that every [ipw()]
#' method returns. It is intended for package developers writing an `ipw()`
#' method, not for end users, and it assumes its arguments are already validated.
#'
#' @details
#' The result layer is shared so that an IPW estimate reads the same way
#' whichever package produced it. A package supplying an `ipw()` method builds
#' its return here and inherits the `print()` and `as.data.frame()` methods
#' registered against the class rather than writing its own. Two packages each
#' defining `print.ipw()` would collide in the shared S3 method table, which is
#' the situation this package exists to prevent.
#'
#' The field names and their order are part of the contract, since callers read
#' fields by name and print the object positionally. `fit` is present on every
#' path, including the ones that have no fitted variance object to report, and
#' `effects` is present whether or not the method that built the result named a
#' mode.
#'
#' `print()` writes the estimand and the call of each component model, then the
#' table of the surface the result's presentation mode names. The section below
#' describes the two modes and what each one tabulates.
#'
#' `as.data.frame()` reports the effect estimates as a tidier-shaped table
#' rather than as a copy of the `estimates` component. Its columns are `term`,
#' then `contrast` when the result names contrasts, then `estimate`,
#' `std.error`, `statistic`, and `p.value`. Those are the names the tidier
#' convention uses, so a fitting package's `tidy()` method is this table read as
#' a tibble and nothing more. The `estimates` component itself is unchanged by
#' any of what follows.
#'
#' `conf.int = TRUE` appends `conf.low` and `conf.high` after the other columns,
#' and `conf.level` names the level they report. The level is an argument rather
#' than a column, since a column would repeat one number down every row and be
#' read as part of the table rather than as the level the two bounds beside it
#' were built at. The bounds `estimates` stores are returned when every row of
#' the frame records the level asked for. They need not be the normal pair: a
#' bootstrap or profile interval is asymmetric about the estimate, and even a
#' normal one rounded on its way into the frame is not the number recomputing
#' gives. At any other level, and for a frame whose rows disagree about the
#' level or record none, the bounds are the normal approximation built from the
#' estimate and its standard error.
#'
#' With `exponentiate = TRUE` the `log(rr)` and `log(or)` rows move to their
#' natural scale, exponentiating the point estimate and the confidence bounds and
#' relabelling the two terms `"rr"` and `"or"`. Standard errors, statistics, and
#' p-values stay on the log scale, where the inference is done, and the interval
#' is settled before the scale is: bounds recomputed at another level are built
#' on the log scale and exponentiated afterwards. The covariance described below
#' travels on the returned table under the same `ipw_vcov` attribute while the
#' rows are on the scale they were estimated on, and is dropped when
#' `exponentiate = TRUE`, since a matrix left attached there would describe
#' neither the table it sits on nor anything else.
#'
#' # The effect labels
#'
#' Every row of `estimates` has a label, and it is the label rather than the
#' position that `print()` writes down the side of its table and that
#' [`coef()`][ipw-accessors], [`vcov()`][ipw-accessors], and
#' [`confint()`][ipw-accessors] name their results with. The label is the
#' `effect` column on its own when there is no `contrast` column, and `effect`
#' and `contrast` pasted together, such as `"rd b vs a"`, when there is. A
#' categorical exposure repeats each effect measure across its contrasts, so
#' `effect` alone would name several rows the same thing.
#'
#' # The presentation mode
#'
#' A result reports its effects in one of two readings, recorded in the
#' `effects` field. The `"marginal"` reading shows the causal contrast
#' estimates the method targeted; the `"conditional"` reading presents the
#' outcome model's coefficient surface. Both surfaces always exist on the
#' object, so the field says which one the result presents rather than which
#' one it holds. [as_marginal()] and [as_conditional()] are how a caller moves a
#' result between them.
#'
#' A printed result names its mode twice, since the two readings are different
#' tables of different numbers: once on an `Effects:` line beside the estimand,
#' and once in the heading of the table itself. The marginal reading tabulates
#' the effect estimates the result stores, under `Marginal estimates:`. The
#' conditional reading tabulates the outcome model's coefficients, under
#' `Conditional estimates (outcome model):`, with the standard errors implied by
#' the corrected covariance a fitting package attaches through [new_ipw_model()].
#' Which coefficient each entry of that block belongs to is what its labels say
#' rather than its row order, so a block attached in another order still prints
#' each standard error beside the coefficient it belongs to.
#'
#' An outcome model a fitting package never wrapped is still printed: the
#' coefficients are written on their own, followed by a note saying that no
#' covariance from the joint estimation is recorded, rather than beside the
#' standard errors the model computed for itself. An outcome model that carries
#' the [new_ipw_model()] class with no covariance behind it is refused instead,
#' with an error of class `causalgenerics_no_vcov_ipw_model`, and one whose block
#' cannot be paired with its coefficients with an error of class
#' `causalgenerics_conditional_vcov_mismatch`. The note answers a package that
#' has not adopted the contract by telling the reader to wrap the model, and that
#' advice has already been taken in both of those cases: the object is what is
#' wrong, not the package that produced it.
#'
#' An outcome model that reports no coefficients has no rows to tabulate under
#' that heading, and the printed form says so in place of the table.
#'
#' # The covariance of the effects
#'
#' A method that can compute the covariance of the effects it reports attaches
#' it to the `estimates` data frame as an attribute named `ipw_vcov`. The value
#' is a square numeric matrix whose row order is the row order of `estimates`
#' and whose dimnames on both margins are the effect labels above.
#' [`vcov()`][ipw-accessors] reads that attribute and raises an error when it is
#' absent, so a method that has no covariance to report attaches none rather
#' than a substitute: the standard errors in `estimates` give the diagonal of the
#' matrix and say nothing about the off-diagonal entries, which are not zero for
#' effects estimated from the same weighted means.
#'
#' @param estimand The causal estimand the method targeted, such as `"ate"` or
#'   `"att"`.
#' @param wt_mod The weighting object: the fitted model that produced the
#'   weights.
#' @param outcome_mod The fitted weighted outcome model.
#' @param estimates A data frame of effect estimates, in the shape the return
#'   value describes.
#' @param se_method The standard error method that ran, such as `"mestimation"`
#'   or `"linearization"`.
#' @param fit The fitted variance object, or `NULL` when the method has none.
#' @param effects The presentation mode the result reports its effects in,
#'   either `"marginal"` or `"conditional"`. A method that names no mode reports
#'   marginal effects.
#'
#' @return `new_ipw()` returns an S3 object of class `ipw`: a list of the
#'   following seven components, in this order.
#' \describe{
#'   \item{`estimand`}{The causal estimand, such as `"ate"` or `"att"`.}
#'   \item{`wt_mod`}{The weighting object: the fitted model that produced the
#'     weights.}
#'   \item{`outcome_mod`}{The fitted outcome model.}
#'   \item{`estimates`}{A data frame with one row per effect measure and the
#'     following columns: `effect` (the measure name), `estimate` (point
#'     estimate), `std.err` (standard error), `z` (z-statistic), `ci.lower` and
#'     `ci.upper` (confidence interval bounds), `conf.level`, and `p.value`. For
#'     a categorical exposure the data frame also has a `contrast` column,
#'     placed after `effect`, naming the non-reference level and reference level
#'     of each contrast. A frame stored against an earlier version of this
#'     contract names that column `comparison`. The older name is read as an
#'     alias for the canonical one wherever the column is read, so a result
#'     holding such a frame labels its rows and reports its table exactly as one
#'     holding a `contrast` column does. A method written now writes
#'     `contrast`.}
#'   \item{`se_method`}{The standard error method used, such as `"mestimation"`
#'     or `"linearization"`.}
#'   \item{`fit`}{The fitted object the variance estimator produced, or `NULL`.
#'     A method that stacks estimating equations records the M-estimator here;
#'     the linearization path has no such object and records `NULL`.}
#'   \item{`effects`}{The presentation mode, either `"marginal"` or
#'     `"conditional"`. The marginal reading shows the causal contrast
#'     estimates and the conditional reading presents the outcome model's
#'     coefficient surface; both surfaces exist on every result. See
#'     [as_marginal()] and [as_conditional()].}
#' }
#'
#'   `print()` returns its input invisibly. `as.data.frame()` returns a plain
#'   data frame of the effect estimates under the tidier column names described
#'   above, with the confidence bounds appended when they are asked for.
#'
#' @seealso [ipw()], the generic these results come from, and [as_marginal()]
#'   and [as_conditional()] for the presentation mode.
#'
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
#'   z = rep(c(0, 1), 10),
#'   y = rep(c(0, 1, 1, 0, 1), 4)
#' )
#'
#' # Written out literally, in the shape the return contract documents. These
#' # stand in for what an `ipw()` method would compute from the models below.
#' estimates <- data.frame(
#'   effect = c("rd", "log(rr)", "log(or)"),
#'   estimate = c(0.199882, 0.560414, 0.878313),
#'   std.err = c(0.092425, 0.273519, 0.418661),
#'   z = c(2.1626, 2.0489, 2.0979),
#'   ci.lower = c(0.018732, 0.024326, 0.057753),
#'   ci.upper = c(0.381032, 1.096502, 1.698873),
#'   conf.level = 0.95,
#'   p.value = c(0.030570, 0.040470, 0.035910)
#' )
#'
#' res <- new_ipw(
#'   estimand = "ate",
#'   wt_mod = glm(z ~ x, family = binomial(), data = dat),
#'   outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat),
#'   estimates = estimates,
#'   se_method = "linearization",
#'   fit = NULL
#' )
#'
#' res
#'
#' # The tidier-shaped table the result reports.
#' as.data.frame(res)
#'
#' # With an interval, and the ratios on their natural scale.
#' as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)
new_ipw <- function(
  estimand,
  wt_mod,
  outcome_mod,
  estimates,
  se_method,
  fit,
  effects = "marginal"
) {
  # The mode is the one field with a fixed set of values, and a misspelling
  # stored unchecked would sit in the result until something downstream branched
  # on it and took the branch neither reading names.
  check_ipw_effects(effects)

  structure(
    list(
      estimand = estimand,
      wt_mod = wt_mod,
      outcome_mod = outcome_mod,
      estimates = estimates,
      se_method = se_method,
      fit = fit,
      effects = effects
    ),
    class = "ipw"
  )
}

#' @param x An `ipw` object.
#' @param ... Further arguments. Neither method reads them.
#' @rdname new_ipw
#' @export
print.ipw <- function(x, ...) {
  # The mode decides which table is written, and it is read first so that a
  # result recording something that is not a reading is refused before any of
  # the summary reaches the console.
  effects <- ipw_effects(x)

  cat("Inverse Probability Weight Estimator\n")
  cat("Estimand:", toupper(x$estimand), "\n")
  # The mode is a fact about the result rather than about the table, so it is
  # reported here beside the estimand as well as in the table's own heading.
  cat("Effects:", ipw_effects_label(effects), "\n\n")

  cat("Weight Estimator:\n")
  cat("  Call:", format_model_call(x$wt_mod), "\n")
  cat("\n")

  cat("Outcome Model:\n")
  cat("  Call:", format_model_call(x$outcome_mod), "\n")

  cat("\n")

  if (effects == "conditional") {
    print_conditional_estimates(x$outcome_mod)
  } else {
    print_marginal_estimates(x$estimates)
  }

  invisible(x)
}

#' The presentation mode as the printed form names it
#'
#' The field holds the reading's name and nothing more, and a printed line
#' saying only "marginal" or "conditional" would say that the estimates are
#' marginal or conditional without saying what over. The parenthetical completes
#' each one.
#'
#' @param effects The mode the result records, either `"marginal"` or
#'   `"conditional"`.
#'
#' @return A single string.
#'
#' @noRd
ipw_effects_label <- function(effects) {
  switch(
    effects,
    marginal = "marginal (population-averaged)",
    conditional = "conditional (outcome model)"
  )
}

#' Write the table the marginal reading reports
#'
#' The effect estimates the result stores, keyed by the effect labels the
#' [new_ipw()] contract defines.
#'
#' @param estimates The `estimates` component of an `ipw` object.
#'
#' @return `NULL`, invisibly. Called for the table it writes.
#'
#' @noRd
#' @importFrom stats printCoefmat
print_marginal_estimates <- function(estimates) {
  cat("Marginal estimates:\n")

  # The rows are keyed by effect label, and the character columns the labels are
  # built from are dropped from the numeric matrix printCoefmat() formats. Both
  # spellings of the contrast column go, rather than whichever one
  # `ipw_contrast_column()` reads: the question here is which columns must not
  # reach `printCoefmat()`, which is every character column the frame might
  # carry, and not which one the labels are built from. A frame carrying both
  # would otherwise keep the unread one, and that does not error. It is
  # factor-coded by `data.matrix()` into the first numeric position, which shifts
  # every column along one and leaves `cs.ind` and `tst.ind` naming the wrong
  # ones, so the standard errors are formatted as a test statistic.
  numbers <- estimates[
    setdiff(names(estimates), c("effect", "contrast", "comparison"))
  ]
  rownames(numbers) <- ipw_effect_labels(estimates)
  stats::printCoefmat(numbers, has.Pvalue = TRUE, cs.ind = 1:2, tst.ind = 3)

  invisible(NULL)
}

#' Write the table the conditional reading reports
#'
#' The outcome model's coefficients, with the standard errors implied by the
#' covariance the joint estimation of the weights and the outcome gives, which a
#' fitting package attaches with [new_ipw_model()].
#'
#' A model a fitting package never wrapped is summarized rather than refused.
#' `print()` is the view of whatever a caller is holding, and refusing there
#' would leave a result that cannot be looked at at all. The coefficients are
#' written on their own in that case, with a note saying what is missing, rather
#' than beside the standard errors the model computed for itself: those treat the
#' estimated weights as fixed and report an uncertainty the coefficients do not
#' have, and a column of them under this heading would be read as the corrected
#' ones.
#'
#' A model carrying the `ipw_model` class with no covariance behind it is refused
#' instead, and the error travels out of `print()`. So is one whose block cannot
#' be paired with its coefficients. The note has nothing to tell either caller:
#' it answers a package that has not adopted the contract by saying to wrap the
#' model, and both of these models are wrapped. What is wrong is the object
#' itself, which no constructor here could have produced, and an object in that
#' state should fail the same way wherever it is met rather than reading as a
#' result with one part missing.
#'
#' A model that reports no coefficients has no rows to tabulate under this
#' heading whatever covariance it carries, so it is told in a sentence and the
#' covariance is never asked for.
#'
#' @param model The result's `outcome_mod`.
#'
#' @return `NULL`, invisibly. Called for the table it writes.
#'
#' @noRd
#' @importFrom stats coef pnorm printCoefmat
print_conditional_estimates <- function(model) {
  # Taken here rather than left to the helper's own default. The covariance is
  # asked for inside `tryCatch()`, and a call the helper resolved from there
  # would name one of `tryCatch()`'s internal frames instead of the method the
  # caller invoked. Reading it eagerly, before anything else runs, is what makes
  # a condition raised further down name `print.ipw()` the way the one an
  # accessor raises names the accessor.
  call <- sys.call(-1)

  cat("Conditional estimates (outcome model):\n")

  estimate <- stats::coef(model)

  # A model with no coefficients has no rows to tabulate, and the reason the
  # table is absent takes its place. The guard runs before the covariance is
  # looked up, since a table with no rows has no use for one and reaching for it
  # would refuse the result over a part of it nothing was going to read.
  if (length(estimate) == 0L) {
    cat(
      "The outcome model reports no coefficients, so there is no\n",
      "conditional table to print.\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  # Asked for through the helper the accessors read it with, so that what counts
  # as a covariance this reading can report is settled in one place. Two
  # conditions come back from it and only one is handled here. A model that was
  # never wrapped raises `causalgenerics_no_conditional_vcov`, which is the case
  # the note below answers. A model wrapped and then stripped raises
  # `causalgenerics_no_vcov_ipw_model`, which is left to travel past this handler
  # and out of `print()`.
  covariance <- tryCatch(
    conditional_vcov(model, call = call),
    causalgenerics_no_conditional_vcov = function(cnd) NULL
  )

  if (is.null(covariance)) {
    stats::printCoefmat(
      cbind(Estimate = estimate),
      cs.ind = 1L,
      tst.ind = integer(),
      has.Pvalue = FALSE,
      signif.stars = FALSE
    )
    cat("\n")
    # Printed rather than signalled, so that it travels with the output a caller
    # captures and stays beside the table it qualifies.
    cat(
      "Standard errors are not reported: this result's outcome model records\n",
      "no covariance from the joint estimation of the weights and the outcome.\n",
      "The package that produced it attaches one by wrapping the model with\n",
      "`new_ipw_model()`.\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  standard_error <- sqrt(diag(covariance))
  z <- estimate / standard_error
  stats::printCoefmat(
    cbind(
      Estimate = estimate,
      `Std. Error` = standard_error,
      `z value` = z,
      `Pr(>|z|)` = 2 * stats::pnorm(-abs(z))
    ),
    has.Pvalue = TRUE,
    cs.ind = 1:2,
    tst.ind = 3
  )

  invisible(NULL)
}

#' The label each row of an `ipw` result's estimates carries
#'
#' `print()` writes these down the side of its table, `coef()` names its vector
#' with them, `vcov()` uses them as dimnames, and `confint()` both labels its
#' rows and matches a character `parm` against them. One helper so that the
#' surfaces cannot drift: a caller who reads a covariance out by the name
#' `coef()` gave has to get the entry `print()` showed.
#'
#' A categorical exposure repeats each effect measure across its contrasts, so a
#' frame that names contrasts needs both columns to name a row uniquely. The
#' labels stay unique because the contrast labels are distinct.
#'
#' The label is the same string whichever column name a stored frame keeps its
#' contrasts under, since the column supplies the second half of the label and
#' not its own name. A result built against the earlier contract therefore
#' reports the labels it always reported.
#'
#' @param estimates The `estimates` component of an `ipw` object.
#'
#' @return A character vector, one element per row of `estimates`.
#'
#' @noRd
ipw_effect_labels <- function(estimates) {
  contrast <- ipw_contrast_column(estimates)
  if (is.null(contrast)) {
    return(as.character(estimates$effect))
  }
  paste(estimates$effect, estimates[[contrast]])
}

#' The column of an estimates frame that names its contrasts
#'
#' A categorical exposure names each row's contrast in a column of the
#' `estimates` frame, and the [new_ipw()] contract calls that column `contrast`.
#' A frame stored against an earlier version of the contract calls it
#' `comparison`, and a result holding one is still a result a caller has in hand,
#' so the older name is read as an alias for the canonical one rather than as a
#' shape this layer has no reading for. A frame carrying both is read under the
#' canonical name.
#'
#' The alias is resolved here and nowhere else, so that the three surfaces which
#' read the column cannot come to disagree about which frames name contrasts.
#' What it buys is on the way in only. `print()` labels its rows, the accessors
#' name their results, and `as.data.frame()` heads its column the canonical way
#' for a frame of either vintage, so nothing a caller reads says which one the
#' result was built from.
#'
#' @param estimates The `estimates` component of an `ipw` object.
#'
#' @return A single string naming the column, or `NULL` when the frame names no
#'   contrasts.
#'
#' @noRd
ipw_contrast_column <- function(estimates) {
  # `intersect()` keeps the order of its first argument, which is what puts the
  # canonical name ahead of the alias.
  columns <- intersect(c("contrast", "comparison"), names(estimates))
  if (length(columns) == 0L) {
    return(NULL)
  }
  columns[[1L]]
}

#' Format a model's originating call for the `ipw()` summary
#'
#' Objects that carry an accessible call, such as `glm` and `lm`, report the
#' deparsed call. `getCall()` reaches the call through `getElement()`, so an
#' object that records no call gives back `NULL` and an object that cannot be
#' subset raises a condition instead. Both fall back to a class label, which is
#' what lets `print()` work for a weighting object of any shape.
#'
#' @param mod A fitted model or weighting object.
#'
#' @return A single string.
#'
#' @importFrom stats getCall
#' @noRd
format_model_call <- function(mod) {
  call <- tryCatch(stats::getCall(mod), error = function(e) NULL)
  if (is.null(call)) {
    return(paste0("<", paste(class(mod), collapse = "/"), ">"))
  }
  paste(deparse(call), collapse = "\n")
}

#' @param row.names A character vector of row names for the returned table, or
#'   `NULL` for the automatic ones.
#' @param optional Accepted for the [base::as.data.frame()] generic. Every
#'   column of the table is named, so there is nothing for it to make optional.
#' @param conf.int If `TRUE`, append `conf.low` and `conf.high` columns after
#'   the rest of the table. Default is `FALSE`.
#' @param conf.level The confidence level the bounds report. At the level every
#'   row of `estimates` records, the stored bounds are returned; at any other
#'   level they are the normal approximation built from the estimate and its
#'   standard error. Default is `0.95`.
#' @param exponentiate If `TRUE`, exponentiate the log risk ratio and log odds
#'   ratio to produce risk ratios and odds ratios on their natural scale,
#'   relabelling the two terms `"rr"` and `"or"`. The confidence bounds move with
#'   them. Standard errors, statistics, and p-values remain on the log scale, and
#'   the `ipw_vcov` attribute is dropped rather than carried, since it describes
#'   the estimates on the scale they were estimated on. Default is `FALSE`.
#' @rdname new_ipw
#' @export
as.data.frame.ipw <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE
) {
  # All three are checked on every call rather than on the branch that reads
  # them. `conf.level` means the same thing whether or not the bounds are
  # reported, so a call naming a level no interval can be built at has said
  # something wrong whichever way `conf.int` was set, and checking it only where
  # it is read would accept that call now and refuse the same value the moment
  # the bounds were asked for.
  check_flag(conf.int, "conf.int")
  check_conf_level(conf.level)
  check_flag(exponentiate, "exponentiate")

  estimates <- x$estimates

  # `term` first, then the column naming the contrast it qualifies when the
  # result reports one. A binary or continuous exposure has a single contrast,
  # so a `contrast` column there would repeat one value down the table and read
  # as a contrast that was named. The heading is the canonical one whichever
  # name the stored frame used, since the heading is what a caller reads the
  # table by.
  contrast <- ipw_contrast_column(estimates)
  columns <- list(term = as.character(estimates$effect))
  if (!is.null(contrast)) {
    columns$contrast <- estimates[[contrast]]
  }
  columns$estimate <- estimates$estimate
  columns$std.error <- estimates$std.err
  columns$statistic <- estimates$z
  columns$p.value <- estimates$p.value

  # Built before the scale is changed, so that a recomputed bound is a half
  # width on the log scale added to an estimate on the log scale.
  bounds <- if (conf.int) interval_bounds(estimates, conf.level) else NULL

  if (exponentiate) {
    # The rows to move are the ones labelled `log(rr)` and `log(or)`, matched
    # exactly so that a frame whose ratios are already on the natural scale is
    # left alone rather than exponentiated a second time.
    is_log_rr <- columns$term == "log(rr)"
    is_log_or <- columns$term == "log(or)"
    ratios <- is_log_rr | is_log_or

    # Only the point estimate and the confidence bounds move to the natural
    # scale. `std.error`, `statistic`, and `p.value` stay on the log scale,
    # which is where the inference is done.
    columns$estimate[ratios] <- exp(columns$estimate[ratios])
    if (!is.null(bounds)) {
      bounds$lower[ratios] <- exp(bounds$lower[ratios])
      bounds$upper[ratios] <- exp(bounds$upper[ratios])
    }

    # The label names the scale, so it moves with the value it labels.
    columns$term[is_log_rr] <- "rr"
    columns$term[is_log_or] <- "or"
  }

  # Last, after the columns the table always carries, so that asking for an
  # interval adds to the table rather than rearranging it.
  if (!is.null(bounds)) {
    columns$conf.low <- bounds$lower
    columns$conf.high <- bounds$upper
  }

  df <- data.frame(columns, row.names = row.names, stringsAsFactors = FALSE)

  # The covariance belongs to the estimates rather than to a column of them, so
  # it travels on the table under the attribute the `new_ipw()` contract names
  # it under. A result whose method attached none carries none: assigning `NULL`
  # sets no attribute, which is the right answer rather than an accident.
  if (!exponentiate) {
    attr(df, "ipw_vcov") <- attr(estimates, "ipw_vcov", exact = TRUE)
  }

  df
}

#' The confidence bounds `as.data.frame()` reports
#'
#' The stored level is a property of the frame rather than of a row. The bounds
#' `estimates` records come back only when every row of it was reported at the
#' level asked for, since a table mixing stored bounds with recomputed ones
#' would report two intervals under one pair of column headings. A frame with no
#' `conf.level` column says nothing about what its bounds describe, so there is
#' no level for the stored pair to be returned at and the approximation is built
#' for whatever was asked for.
#'
#' The bounds are one half width added and subtracted, rather than a bound from
#' each tail: `qnorm()` is not exactly antisymmetric, and taking the lower bound
#' from the lower tail would put the two a little apart from each other.
#'
#' @param estimates The `estimates` component of an `ipw` object.
#' @param conf.level The level the bounds report.
#'
#' @return A list of two numeric vectors, `lower` and `upper`.
#'
#' @noRd
#' @importFrom stats qnorm
interval_bounds <- function(estimates, conf.level) {
  stored <- estimates$conf.level
  if (!is.null(stored) && isTRUE(all(stored == conf.level))) {
    return(list(lower = estimates$ci.lower, upper = estimates$ci.upper))
  }

  half_width <- stats::qnorm(1 - (1 - conf.level) / 2) * estimates$std.err
  list(
    lower = estimates$estimate - half_width,
    upper = estimates$estimate + half_width
  )
}

#' Refuse an argument that is not a flag
#'
#' A number and a string both take the true branch of a bare `if`, and a missing
#' logical errors from the `if` itself with a message about the condition rather
#' than about the argument the caller wrote. All three are refused here, so that
#' what an argument had to be is said once and in the same words wherever it is
#' said.
#'
#' @param value The argument as the caller supplied it.
#' @param arg The argument's name, which the message reports.
#' @param call The call to report the error against, which is the method's
#'   rather than this helper's.
#'
#' @return `value`, invisibly, when it is a flag.
#'
#' @noRd
check_flag <- function(value, arg, call = sys.call(-1)) {
  is_flag <- is.logical(value) && length(value) == 1L && !is.na(value)

  if (!is_flag) {
    stop_invalid_argument(
      arg,
      "be a single logical value, either `TRUE` or `FALSE`",
      call = call
    )
  }

  invisible(value)
}

#' Refuse a confidence level that is not a probability
#'
#' The bounds of the open interval are refused along with everything outside it:
#' the normal approximation gives no finite pair of numbers at level 0 or 1. A
#' vector of levels names several intervals where the table reports one, and a
#' missing level names none.
#'
#' @param conf.level The argument as the caller supplied it.
#' @param call The call to report the error against, which is the method's
#'   rather than this helper's.
#'
#' @return `conf.level`, invisibly, when it is a probability.
#'
#' @noRd
check_conf_level <- function(conf.level, call = sys.call(-1)) {
  is_level <- is.numeric(conf.level) &&
    length(conf.level) == 1L &&
    !is.na(conf.level) &&
    conf.level > 0 &&
    conf.level < 1

  if (!is_level) {
    stop_invalid_argument(
      "conf.level",
      "be a single number greater than 0 and less than 1",
      call = call
    )
  }

  invisible(conf.level)
}
