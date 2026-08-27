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
#' `effects` and `readings` are present whether or not the method that built the
#' result named a mode or a set of readings.
#'
#' `print()` writes the estimand and the call of each component model, then the
#' table of the surface the result's presentation mode names. The section below
#' describes the two modes and what each one tabulates.
#'
#' `as.data.frame()` reports the result as a tidier-shaped table rather than as
#' a copy of the `estimates` component. In the marginal reading its columns are
#' `term`, then `contrast` when the result names contrasts, then `group` when it
#' names subgroups, then `estimate`, `std.error`, `statistic`, and `p.value`.
#' Those are the names the tidier convention uses, so a fitting package's
#' `tidy()` method is this table read as a tibble and nothing more. The
#' `estimates` component itself is unchanged by any of what follows.
#'
#' The conditional reading reports the outcome model's coefficients under those
#' same headings, one row per coefficient, with `term` naming the coefficient
#' and `std.error` the standard error the corrected covariance implies. It names
#' no contrasts and no subgroups: those key the rows of the effects table, and a
#' coefficient is named by itself. It is the table
#' [`as.data.frame()`][ipw-pooled-methods] reports for a pooled result in the
#' same reading, without the `df` column a pooled result has and an unpooled one
#' does not.
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
#' estimate and its standard error. The conditional reading has no stored pair
#' to prefer, since the bounds the frame holds belong to the effects the other
#' reading reports, so its bounds are the normal approximation at every level.
#'
#' With `exponentiate = TRUE` the `log(rr)` and `log(or)` rows of the marginal
#' table move to their natural scale, exponentiating the point estimate and the
#' confidence bounds and relabeling the two terms `"rr"` and `"or"`. Standard
#' errors, statistics, and p-values stay on the log scale, where the inference is
#' done, and the interval is settled before the scale is: bounds recomputed at
#' another level are built on the log scale and exponentiated afterwards. The
#' covariance described below travels on the returned table under the same
#' `ipw_vcov` attribute while the rows are on the scale they were estimated on,
#' and is dropped when `exponentiate = TRUE`, since a matrix left attached there
#' would describe neither the table it sits on nor anything else.
#'
#' A conditional table has no rows labeled as ratios to pick out, so the link the
#' outcome model was fitted with settles the question for the whole table: a
#' `logit` link puts every coefficient on the log odds scale and a `log` link
#' puts every coefficient on the log risk scale, and both are scales an
#' exponential undoes. Every estimate moves and no term is relabeled, since a
#' coefficient name names the term rather than the scale its estimate is
#' reported on. Every other link raises an error of class
#' `causalgenerics_exponentiate_link`, and of the classes
#' `causalgenerics_invalid_argument_exponentiate` and
#' `causalgenerics_invalid_argument`, rather than exponentiating coefficients
#' that describe nothing once exponentiated. That is the rule
#' [`as.data.frame()`][ipw-pooled-methods] keeps for a pooled result, in the same
#' words.
#'
#' # The effect labels
#'
#' Every row of `estimates` has a label, and it is the label rather than the
#' position that `print()` writes down the side of its table and that
#' [`coef()`][ipw-accessors], [`vcov()`][ipw-accessors], and
#' [`confint()`][ipw-accessors] name their results with. The label is the
#' columns that name a row, pasted together in the order they appear in the
#' frame's contract: the `effect` column, then `contrast` when the result names
#' contrasts, then `group` when it names subgroups. A binary exposure estimated
#' for the whole sample is labeled `"rd"`, a categorical one `"rd b vs a"`, and
#' one reported by subgroup `"rd b vs a sex = 0"`. A categorical exposure repeats
#' each effect measure across its contrasts and a subgroup analysis repeats each
#' one across its groups, so `effect` alone would name several rows the same
#' thing.
#'
#' # The presentation mode
#'
#' A result reports its effects in one of two readings, recorded in the
#' `effects` field. The `"marginal"` reading shows the causal contrast
#' estimates the method targeted; the `"conditional"` reading presents the
#' outcome model's coefficient surface. The field says which one the result
#' presents, and [as_marginal()] and [as_conditional()] are how a caller moves a
#' result between them.
#'
#' Which readings a result can present at all is the separate fact the
#' `readings` field records. Both surfaces exist on most results and not on all
#' of them: an exposure entering the outcome model through several columns has
#' no single coefficient to read as the conditional effect, and a package that
#' builds such a result records the reading it can answer for. A result records
#' both unless the method that built it said otherwise, and the mode it records
#' has to be one of them, so a result never presents a reading it does not
#' support. Asking one for the reading it does not support raises an error of
#' class `causalgenerics_unsupported_reading`, whether it is asked through the
#' mode generics, through [`coef()`][ipw-accessors],
#' [`vcov()`][ipw-accessors], or [`confint()`][ipw-accessors], or through the
#' `effects` argument of `as.data.frame()`. A result stored before the field
#' existed records no set and is read as supporting both.
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
#' @param effects For `new_ipw()`, the presentation mode the result reports its
#'   effects in, either `"marginal"` or `"conditional"`; a method that names no
#'   mode reports marginal effects. For `as.data.frame()`, the reading to
#'   report: `NULL`, the default, reports the reading the result records, and
#'   any other value overrides it for the one call and leaves the result as it
#'   is. A reading the result does not support is refused with an error of class
#'   `causalgenerics_unsupported_reading`.
#' @param readings The readings the result supports, one or both of
#'   `"marginal"` and `"conditional"`. The default is both, which is what every
#'   result supported before the field existed. A method whose analysis has no
#'   meaning under one of them names the other, and the mode the result records
#'   has to be one of the readings named.
#'
#' @return `new_ipw()` returns an S3 object of class `ipw`: a list of the
#'   following eight components, in this order.
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
#'     `contrast`. A result reported once per level of a grouping variable also
#'     has a `group` column, placed after the contrast column when one is
#'     present and after `effect` when the result names no contrasts, naming the
#'     subgroup each row was estimated in as a `"var = value"` string such as
#'     `"sex = 0"`. That column has one spelling and no alias. Both optional
#'     columns are absent rather than constant when the result reports one
#'     contrast or one group, since a column repeating a single value down the
#'     table would read as a contrast or a subgroup that was named. `group` must
#'     be character and must name a subgroup in every row; the constructor
#'     refuses anything else with an error of class
#'     `causalgenerics_invalid_argument_estimates`, and of the general class
#'     `causalgenerics_invalid_argument`, since a label pasted from such a column
#'     would relabel every row it keys without failing anywhere.}
#'   \item{`se_method`}{The standard error method used, such as `"mestimation"`
#'     or `"linearization"`.}
#'   \item{`fit`}{The fitted object the variance estimator produced, or `NULL`.
#'     A method that stacks estimating equations records the M-estimator here;
#'     the linearization path has no such object and records `NULL`.}
#'   \item{`effects`}{The presentation mode, either `"marginal"` or
#'     `"conditional"`. The marginal reading shows the causal contrast
#'     estimates and the conditional reading presents the outcome model's
#'     coefficient surface. See [as_marginal()] and [as_conditional()].}
#'   \item{`readings`}{The readings the result supports, one or both of
#'     `"marginal"` and `"conditional"`, and always including the mode above. A
#'     result stored before the field existed carries fewer components and is
#'     read as supporting both readings.}
#' }
#'
#'   `print()` returns its input invisibly. `as.data.frame()` returns a plain
#'   data frame under the tidier column names described above, with the
#'   confidence bounds appended when they are asked for: the effect estimates in
#'   the marginal reading, and the outcome model's coefficients, one row each,
#'   in the conditional one.
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
#'
#' # A result whose analysis has no conditional reading records the one it
#' # supports, and the other is refused rather than reported.
#' marginal_only <- new_ipw(
#'   estimand = "ate",
#'   wt_mod = glm(z ~ x, family = binomial(), data = dat),
#'   outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat),
#'   estimates = estimates,
#'   se_method = "linearization",
#'   fit = NULL,
#'   readings = "marginal"
#' )
#'
#' try(as.data.frame(marginal_only, effects = "conditional"))
new_ipw <- function(
  estimand,
  wt_mod,
  outcome_mod,
  estimates,
  se_method,
  fit,
  effects = "marginal",
  readings = c("marginal", "conditional")
) {
  # The mode is the one field with a fixed set of values, and a misspelling
  # stored unchecked would sit in the result until something downstream branched
  # on it and took the branch neither reading names.
  check_ipw_effects(effects)

  # The readings are checked against the mode as well as on their own, since a
  # result recording a reading it does not support is a result every method
  # downstream would have to have an answer for.
  check_ipw_readings(readings, effects)

  # The optional `group` column is the one part of the estimates frame that is
  # checked here. It completes a row's label, so a column of the wrong kind
  # relabels every row it keys without failing anywhere, which is the shape of
  # mistake a constructor is worth having.
  check_estimates_group(estimates)

  structure(
    list(
      estimand = estimand,
      wt_mod = wt_mod,
      outcome_mod = outcome_mod,
      estimates = estimates,
      se_method = se_method,
      fit = fit,
      effects = effects,
      readings = readings
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
  print_effect_table(estimates)
}

#' Write an estimates frame as a coefficient table
#'
#' The rows are keyed by effect label, and the identity columns the labels are
#' built from are dropped from the numeric matrix `printCoefmat()` formats. They
#' are asked for under every spelling rather than the ones the labels read, for
#' the reason `ipw_identity_columns()` documents: an identity column the frame
#' carries and does not read is still a character column, and one left in the
#' matrix does not error. It is factor-coded by `data.matrix()` into the first
#' numeric position, which shifts every column along one and leaves `cs.ind` and
#' `tst.ind` naming the wrong ones, so the standard errors are formatted as a
#' test statistic.
#'
#' The estimate and its standard error are the first two columns of every frame
#' this writes and the statistic is the third, whether that statistic is
#' referred to the normal or to t, so one set of column roles serves the pooled
#' table and the unpooled one alike.
#'
#' @param estimates An `estimates` frame, from an `ipw` or an `ipw_pooled`
#'   result.
#'
#' @return `NULL`, invisibly. Called for the table it writes.
#'
#' @noRd
#' @importFrom stats printCoefmat
print_effect_table <- function(estimates) {
  numbers <- estimates[
    setdiff(
      names(estimates),
      ipw_identity_columns(estimates, all_spellings = TRUE)
    )
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

  # The same four numbers the tidier-shaped table reports for this reading,
  # written here as a matrix keyed by coefficient name.
  coefficients <- conditional_coefficients(estimate, covariance)
  stats::printCoefmat(
    cbind(
      Estimate = coefficients$estimate,
      `Std. Error` = coefficients$std.error,
      `z value` = coefficients$statistic,
      `Pr(>|z|)` = coefficients$p.value
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
#' A categorical exposure repeats each effect measure across its contrasts, and
#' a result reported by subgroup repeats each measure across those, so a frame
#' that names either needs all of its identity columns to name a row uniquely.
#' The labels stay unique because the contrast labels and the group labels are
#' each distinct.
#'
#' The label is the identity columns pasted together in the order
#' `ipw_identity_columns()` gives them, so a label reads as the measure, then
#' the contrast it compares, then the subgroup it was estimated in:
#' `"rd b vs a sex = 0"`. The group comes last because it qualifies the whole
#' contrast rather than one side of it.
#'
#' The label is the same string whichever column name a stored frame keeps its
#' contrasts under, since the column supplies part of the label and not its own
#' name. A result built against the earlier contract therefore reports the
#' labels it always reported.
#'
#' @param estimates The `estimates` component of an `ipw` object.
#'
#' @return A character vector, one element per row of `estimates`.
#'
#' @noRd
ipw_effect_labels <- function(estimates) {
  columns <- ipw_identity_columns(estimates)
  parts <- lapply(columns, function(column) as.character(estimates[[column]]))

  # `paste()` handed one vector gives it back as a character vector, which is
  # what a frame keyed by the effect alone needs, so the single-column case is
  # the general one with nothing to join rather than a branch of its own. The
  # names are dropped because they name the roles the columns fill, and
  # `paste()` would read them as arguments of its own.
  do.call(paste, unname(parts))
}

#' The columns of an estimates frame that name a row
#'
#' A row's identity is the `effect` column, the column naming the contrast when
#' the frame reports contrasts, and the column naming the subgroup when it
#' reports subgroups, in that order. Both of the last two are optional: a binary
#' or continuous exposure has one contrast, and a result estimated for the whole
#' sample has one group, so a column of either kind would repeat one value down
#' the table and read as a contrast or a subgroup that was named.
#'
#' Four surfaces read that identity, and they have to read the same one: the
#' labels `ipw_effect_labels()` builds, the columns `as.data.frame()` emits, the
#' columns `print_effect_table()` keeps out of its numeric matrix, and the key
#' the pooling carries onto its frames. A caller who reads a covariance out by
#' the name `coef()` gave has to get the entry `print()` showed, and a pooled
#' frame has to be keyed the way the frames it pooled were, so the answer is
#' given here once rather than restated at each of them.
#'
#' The return is named by the role each column fills rather than by the name it
#' goes under in the frame, since a frame stored against an earlier version of
#' the contract keeps its contrasts under `comparison`. A caller building
#' something keyed canonically, such as the pooling, reads the names; a caller
#' reading values out of the frame reads the values.
#'
#' `all_spellings = TRUE` answers a different question: which columns of this
#' frame are identity columns under any name one may be stored under, rather
#' than which ones a label is built from. The two differ for a frame carrying
#' both spellings of the contrast column, which is read under the canonical one.
#' The unread column is still a character column, and `print_effect_table()`
#' needs it gone for a reason that is not about labels: a character column left
#' in the matrix `printCoefmat()` formats does not error, it is factor-coded
#' into a column of numbers beside the real estimates.
#'
#' @param estimates The `estimates` component of an `ipw` or an `ipw_pooled`
#'   object.
#' @param all_spellings If `TRUE`, every identity column the frame carries,
#'   under any name; otherwise the ones a row's identity is read from.
#'
#' @return A character vector of column names, named by the role each column
#'   fills, in identity order.
#'
#' @noRd
ipw_identity_columns <- function(estimates, all_spellings = FALSE) {
  if (all_spellings) {
    # `intersect()` keeps the order of its first argument, which is identity
    # order with the contrast's two spellings side by side.
    return(intersect(
      c("effect", "contrast", "comparison", "group"),
      names(estimates)
    ))
  }

  # A `NULL` from either resolver drops out of the vector rather than becoming
  # an element, which is what makes an absent column absent from the identity.
  c(
    effect = "effect",
    contrast = ipw_contrast_column(estimates),
    group = ipw_group_column(estimates)
  )
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

#' The column of an estimates frame that names its subgroups
#'
#' A result reported once per level of a grouping variable names each row's
#' subgroup in a column of the `estimates` frame, and the [new_ipw()] contract
#' calls that column `group`. There is one spelling and no alias: the column is
#' newer than the contract's first version, so no stored frame keeps its
#' subgroups anywhere else.
#'
#' It is still resolved through a function of its own rather than read inline,
#' the way `ipw_contrast_column()` is. What the surfaces ask is whether this
#' frame names subgroups, and asking it in one place is what keeps the label,
#' the reported table, the printed rows, and the pooled key agreeing about which
#' frames do.
#'
#' @param estimates The `estimates` component of an `ipw` object.
#'
#' @return A single string naming the column, or `NULL` when the frame names no
#'   subgroups.
#'
#' @noRd
ipw_group_column <- function(estimates) {
  if ("group" %in% names(estimates)) "group" else NULL
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
#' @param exponentiate If `TRUE`, move the estimates that are on a log scale to
#'   their natural scale. In the marginal reading those are the log risk ratio
#'   and the log odds ratio, and the two terms are relabeled `"rr"` and `"or"`;
#'   in the conditional reading the outcome model's link settles it for the whole
#'   table, as the section above describes. The confidence bounds move with the
#'   estimates. Standard errors, statistics, and p-values remain on the log
#'   scale, and the `ipw_vcov` attribute is dropped rather than carried, since it
#'   describes the estimates on the scale they were estimated on. Default is
#'   `FALSE`.
#' @rdname new_ipw
#' @export
as.data.frame.ipw <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  effects = NULL
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

  # The table is a presentation of the result rather than a copy of one
  # component of it, so the reading is settled before anything is read off the
  # result and everything built after this belongs to the reading asked for.
  if (resolve_ipw_effects(x, effects) == "conditional") {
    return(conditional_data_frame(
      x,
      row.names = row.names,
      conf.int = conf.int,
      conf.level = conf.level,
      exponentiate = exponentiate,
      call = sys.call()
    ))
  }

  estimates <- x$estimates

  # `term` first, then the column naming the contrast it qualifies when the
  # result reports one, then the column naming the subgroup it was estimated in
  # when the result reports those. That is the order the labels paste them in,
  # so the table names a row the way the printed form and the accessors do.
  columns <- ipw_term_columns(estimates)
  columns$estimate <- estimates$estimate
  columns$std.error <- estimates$std.err
  columns$statistic <- estimates$z
  columns$p.value <- estimates$p.value

  # Built before the scale is changed, so that a recomputed bound is a half
  # width on the log scale added to an estimate on the log scale.
  bounds <- if (conf.int) interval_bounds(estimates, conf.level) else NULL

  if (exponentiate) {
    # The rows to move are the ones labeled `log(rr)` and `log(or)`, matched
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

#' The table the conditional reading reports
#'
#' The outcome model's coefficients, one row each, under the headings the
#' marginal table uses for the effects. The standard errors are the ones the
#' corrected covariance implies, which is the block a fitting package attached
#' with [new_ipw_model()], and the statistic and the p-value are built from
#' them. There is no stored pair of bounds to prefer at any level: the ones the
#' estimates frame holds belong to the effects the other reading reports.
#'
#' A model that was never wrapped carries no such block, and the table is
#' refused rather than built from the standard errors the model computed for
#' itself, which treat the estimated weights as fixed. That is where this table
#' parts company with `print()`, which reports such a model with a note in place
#' of the errors: the printed form is the view of whatever a caller is holding,
#' and a table has no place to put a note.
#'
#' The rows are keyed by the `term` column alone. A conditional table names no
#' contrasts and no subgroups, since those key the rows of the effects table and
#' a coefficient is named by itself, and every other column is unnamed so that
#' the rows are not keyed twice over.
#'
#' @param x An `ipw` object.
#' @param row.names The row names for the returned table, or `NULL`.
#' @param conf.int Whether to append the confidence bounds.
#' @param conf.level The level those bounds report.
#' @param exponentiate Whether to move the estimates to their natural scale.
#' @param call The call to report a refusal against, which is the method's
#'   rather than this helper's.
#'
#' @return A plain data frame.
#'
#' @noRd
#' @importFrom stats coef qnorm
conditional_data_frame <- function(
  x,
  row.names,
  conf.int,
  conf.level,
  exponentiate,
  call = sys.call(-1)
) {
  # Before anything is built, so a table that cannot be reported on the scale
  # asked for is refused rather than half built. The reading has no rows labeled
  # as ratios to pick out, so the link settles it for the whole table at once.
  if (exponentiate) {
    check_exponentiate_link(ipw_outcome_link(x$outcome_mod), call = call)
  }

  covariance <- conditional_vcov(x$outcome_mod, call = call)
  estimate <- stats::coef(x$outcome_mod)
  coefficients <- conditional_coefficients(estimate, covariance)

  # The term is the coefficient's name, and the numbers are unnamed: the names
  # they carry are the same coefficient names, and a column keyed by them beside
  # a column holding them would key the rows twice over.
  columns <- list(term = names(estimate))
  columns$estimate <- unname(coefficients$estimate)
  columns$std.error <- unname(coefficients$std.error)
  columns$statistic <- unname(coefficients$statistic)
  columns$p.value <- unname(coefficients$p.value)

  # Built before the scale is changed, so that a bound is a half width on the
  # log scale added to an estimate on the log scale.
  bounds <- if (conf.int) {
    half_width <- stats::qnorm(1 - (1 - conf.level) / 2) * columns$std.error
    list(
      lower = columns$estimate - half_width,
      upper = columns$estimate + half_width
    )
  } else {
    NULL
  }

  if (exponentiate) {
    # Every row moves, since the link says every coefficient is on a scale an
    # exponential undoes. No term is relabeled: a coefficient name names the
    # term rather than the scale its estimate is reported on.
    columns$estimate <- exp(columns$estimate)
    if (!is.null(bounds)) {
      bounds$lower <- exp(bounds$lower)
      bounds$upper <- exp(bounds$upper)
    }
  }

  ipw_tidy_frame(
    columns,
    bounds,
    covariance = if (!exponentiate) covariance else NULL,
    row.names = row.names
  )
}

#' The coefficient surface the conditional reading reports
#'
#' The outcome model's coefficients and the three numbers the corrected
#' covariance implies for them: the standard error, the statistic each estimate
#' gives, and the p-value that statistic is referred to the normal at. The
#' printed table and the tidier-shaped one report the same surface in different
#' shapes, so the arithmetic is here rather than in each of them.
#'
#' The names the model gives its coefficients are kept, since the printed table
#' is keyed by them.
#'
#' @param estimate The outcome model's coefficients.
#' @param covariance Their corrected covariance, in coefficient order.
#'
#' @return A list of four numeric vectors, `estimate`, `std.error`,
#'   `statistic`, and `p.value`.
#'
#' @noRd
#' @importFrom stats pnorm
conditional_coefficients <- function(estimate, covariance) {
  standard_error <- sqrt(diag(covariance))
  statistic <- estimate / standard_error

  list(
    estimate = estimate,
    std.error = standard_error,
    statistic = statistic,
    p.value = 2 * stats::pnorm(-abs(statistic))
  )
}

#' The columns of a tidier-shaped table that name its rows
#'
#' `term` first, then the column naming the contrast it qualifies when the
#' result reports one, then the column naming the subgroup it was estimated in
#' when the result reports those. That is the order the effect labels paste them
#' in, so a table names a row the way the printed form and the accessors do.
#'
#' A binary or continuous exposure has a single contrast and an ungrouped result
#' has a single group, so a column of either kind there would repeat one value
#' down the table and read as a contrast or a subgroup that was named. The
#' headings are the canonical ones whichever names the stored frame used, since
#' the heading is what a caller reads the table by.
#'
#' The pooled table is keyed the same way, from a frame the pooling carried the
#' same key onto, so both methods build these columns here.
#'
#' @param estimates The `estimates` component of an `ipw` or an `ipw_pooled`
#'   object.
#'
#' @return A named list of one to three character columns.
#'
#' @noRd
ipw_term_columns <- function(estimates) {
  contrast <- ipw_contrast_column(estimates)
  group <- ipw_group_column(estimates)

  columns <- list(term = as.character(estimates$effect))
  if (!is.null(contrast)) {
    columns$contrast <- estimates[[contrast]]
  }
  if (!is.null(group)) {
    columns$group <- estimates[[group]]
  }
  columns
}

#' Assemble a tidier-shaped table from its columns
#'
#' The last few steps every table this package reports has in common: the
#' confidence bounds go on after the columns the table always carries, so that
#' asking for an interval adds to the table rather than rearranging it, and the
#' covariance travels on the finished frame rather than in a column of it, under
#' the attribute the [new_ipw()] contract names it under.
#'
#' A `NULL` covariance sets no attribute, which is what a result that carries
#' none and a table reported on another scale both want. It is the right answer
#' rather than an accident of assignment, so the callers pass `NULL` deliberately
#' rather than skipping the assignment.
#'
#' @param columns The named list of columns, in order.
#' @param bounds The confidence bounds, or `NULL` when none were asked for.
#' @param covariance The covariance to attach, or `NULL`.
#' @param row.names The row names for the returned table, or `NULL` for the
#'   automatic ones.
#'
#' @return A plain data frame.
#'
#' @noRd
ipw_tidy_frame <- function(columns, bounds, covariance, row.names = NULL) {
  if (!is.null(bounds)) {
    columns$conf.low <- bounds$lower
    columns$conf.high <- bounds$upper
  }

  df <- data.frame(columns, row.names = row.names, stringsAsFactors = FALSE)
  attr(df, "ipw_vcov") <- covariance
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

#' Refuse an estimates frame whose subgroups are not named
#'
#' The `group` column completes a row's label, and a label is a string built by
#' pasting the identity columns together, so anything at all pastes into
#' something and nothing downstream errors. What comes out is a label the
#' package that produced the result did not write, which is why the column is
#' checked where the result is constructed rather than where a label is built.
#'
#' The values are `"var = value"` strings, so the column has to be character. A
#' factor pastes as its levels rather than as what the frame holds, and it
#' carries levels the frame may have no rows for; a number pastes as a formatted
#' number, which names a value without saying which variable took it.
#'
#' Every row it groups has to name a group. A missing one pastes into a label
#' reading `"rd NA"`, which describes no subgroup, and two of them paste into the
#' same label, which leaves `coef()` naming two elements the same thing and
#' `printCoefmat()` with row names it cannot be handed.
#'
#' A frame that names no subgroups is the frame every fitting package writes
#' today, and there is nothing here for it to fail.
#'
#' @param estimates The `estimates` argument as the caller supplied it.
#' @param call The call to report the error against, which is the
#'   constructor's rather than this helper's.
#'
#' @return `estimates`, invisibly, when its subgroups are named.
#'
#' @noRd
check_estimates_group <- function(estimates, call = sys.call(-1)) {
  column <- ipw_group_column(estimates)
  if (is.null(column)) {
    return(invisible(estimates))
  }

  group <- estimates[[column]]

  if (!is.character(group)) {
    stop_invalid_argument(
      "estimates",
      paste0(
        "name its subgroups in a character `group` column, since a group is ",
        "written as a \"var = value\" string, but the column this frame ",
        "carries is <",
        class(group)[[1L]],
        ">"
      ),
      call = call
    )
  }

  unnamed <- sum(is.na(group))
  if (unnamed > 0L) {
    stop_invalid_argument(
      "estimates",
      paste0(
        "name a subgroup in every row of its `group` column, but ",
        unnamed,
        if (unnamed == 1L) " row records none" else " rows record none"
      ),
      call = call
    )
  }

  invisible(estimates)
}

#' Refuse a confidence level that is not a probability
#'
#' The bounds of the open interval are refused along with everything outside it:
#' the normal approximation gives no finite pair of numbers at level 0 or 1. A
#' vector of levels names several intervals where the table reports one, and a
#' missing level names none.
#'
#' @param conf.level The argument as the caller supplied it.
#' @param arg The name the caller's function gives the argument, which the
#'   message and the error's keyed class report. Every surface that reports
#'   bounds refuses the same values, and the one that differs between them is
#'   what the caller wrote at the call site: a message about `conf.level` for a
#'   caller who wrote `conf_level` names an argument the call does not have.
#' @param call The call to report the error against, which is the method's
#'   rather than this helper's.
#'
#' @return `conf.level`, invisibly, when it is a probability.
#'
#' @noRd
check_conf_level <- function(
  conf.level,
  arg = "conf.level",
  call = sys.call(-1)
) {
  is_level <- is.numeric(conf.level) &&
    length(conf.level) == 1L &&
    !is.na(conf.level) &&
    conf.level > 0 &&
    conf.level < 1

  if (!is_level) {
    stop_invalid_argument(
      arg,
      "be a single number greater than 0 and less than 1",
      call = call
    )
  }

  invisible(conf.level)
}
