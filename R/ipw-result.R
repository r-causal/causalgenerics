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
#' `as.data.frame()` returns the `estimates` component. With
#' `exponentiate = TRUE` it moves the `log(rr)` and `log(or)` rows to their
#' natural scale, exponentiating the point estimate and the confidence limits and
#' relabelling the two effects `"rr"` and `"or"`. Standard errors, z statistics,
#' and p-values stay on the log scale, where the inference is done.
#'
#' # The effect labels
#'
#' Every row of `estimates` has a label, and it is the label rather than the
#' position that `print()` writes down the side of its table and that
#' [`coef()`][ipw-accessors], [`vcov()`][ipw-accessors], and
#' [`confint()`][ipw-accessors] name their results with. The label is the
#' `effect` column on its own when there is no `comparison` column, and `effect`
#' and `comparison` pasted together, such as `"rd b vs a"`, when there is. A
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
#'     a categorical exposure the data frame also has a `comparison` column,
#'     placed after `effect`, naming the non-reference level and reference level
#'     of each contrast.}
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
#'   `print()` returns its input invisibly. `as.data.frame()` returns the
#'   `estimates` component as a data frame.
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
#' # The ratios on their natural scale.
#' as.data.frame(res, exponentiate = TRUE)
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
#' @param ... Further arguments. `print()` ignores them; `as.data.frame()`
#'   passes them to [base::as.data.frame()].
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
  # built from are dropped from the numeric matrix printCoefmat() formats.
  numbers <- estimates[setdiff(names(estimates), c("effect", "comparison"))]
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
#' frame with a `comparison` column needs both columns to name a row uniquely.
#' The labels stay unique because the comparison labels are distinct.
#'
#' @param estimates The `estimates` component of an `ipw` object.
#'
#' @return A character vector, one element per row of `estimates`.
#'
#' @noRd
ipw_effect_labels <- function(estimates) {
  if ("comparison" %in% names(estimates)) {
    paste(estimates$effect, estimates$comparison)
  } else {
    as.character(estimates$effect)
  }
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

#' @param exponentiate If `TRUE`, exponentiate the log risk ratio and log odds
#'   ratio to produce risk ratios and odds ratios on their natural scale. The
#'   confidence interval bounds are also exponentiated. Standard errors, z
#'   statistics, and p-values remain on the log scale. Default is `FALSE`.
#' @param row.names,optional Passed to [base::as.data.frame()].
#' @rdname new_ipw
#' @export
as.data.frame.ipw <- function(
  x,
  row.names = NULL,
  optional = NULL,
  exponentiate = FALSE,
  ...
) {
  df <- as.data.frame(
    x$estimates,
    row.names = row.names,
    optional = optional,
    ...
  )

  if (!exponentiate) {
    # Return as-is.
    return(df)
  }

  # The rows to move are the ones labelled `log(rr)` and `log(or)`, matched
  # exactly so that a frame whose ratios are already on the natural scale is left
  # alone rather than exponentiated a second time.
  is_log_rr <- df$effect == "log(rr)"
  is_log_or <- df$effect == "log(or)"

  rows_to_expo <- is_log_rr | is_log_or

  # Only the point estimate and the confidence limits move to the natural scale.
  df$estimate[rows_to_expo] <- exp(df$estimate[rows_to_expo])
  df$ci.lower[rows_to_expo] <- exp(df$ci.lower[rows_to_expo])
  df$ci.upper[rows_to_expo] <- exp(df$ci.upper[rows_to_expo])

  # The labels name the scale, so they move with the values.
  df$effect[is_log_rr] <- "rr"
  df$effect[is_log_or] <- "or"

  # `std.err`, `z`, and `p.value` stay on the log scale, which is where the
  # significance testing is done.

  df
}
