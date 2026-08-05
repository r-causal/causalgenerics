#' Model accessors for an inverse probability weighted result
#'
#' @description
#' The accessors a fitted model answers to, defined for the class [new_ipw()]
#' constructs.
#'
#' * `coef()` returns the effect estimates.
#' * `vcov()` returns the covariance of those estimates.
#' * `confint()` returns their confidence limits.
#' * `nobs()` returns the number of observations the outcome model was fitted
#'   on.
#' * `df.residual()` returns the residual degrees of freedom the fitted variance
#'   object reports.
#' * `weights()` returns the weights the outcome model was fitted with.
#' * `model.frame()` returns the outcome model's model frame.
#' * `estimand()` returns the estimand the weights targeted.
#'
#' @details
#' These methods live here for the reason `print()` does. Two packages each
#' registering `coef.ipw()` would collide in the shared S3 method table, and a
#' caller writing against a result would then get whichever package was
#' installed last rather than the contract.
#'
#' Everything they read is part of the [new_ipw()] contract: the `estimates`
#' frame, the `ipw_vcov` attribute attached to it, `estimand`, `outcome_mod`,
#' `fit`, and the `effects` field the section below describes.
#' They never branch on `se_method`, since the fields they read already hold the
#' result of whatever computation that names, and they reach into `fit` only
#' through ordinary S3 dispatch. A package whose variance object is a bare list
#' rather than a fitted model therefore uses them unchanged.
#'
#' Rows are named by the effect labels the [new_ipw()] contract defines, so the
#' name `coef()` gives an estimate is the one `vcov()` and `confint()` use for
#' it and the one `print()` labels its row with.
#'
#' # The reading these methods report
#'
#' `coef()`, `vcov()`, and `confint()` report the surface the result's `effects`
#' field names, and take an `effects` argument that names one for a single call.
#' The marginal reading is the causal contrast estimates, which is what these
#' methods reported before the field existed and what a result that records no
#' mode still reports. The conditional reading is the outcome model's
#' coefficient surface: `coef()` gives the model's coefficients, `vcov()` gives
#' their covariance, and `confint()` bounds them.
#'
#' The covariance the conditional reading reports is the block of the joint
#' estimation that the fitting package attached with [new_ipw_model()], never
#' the one the outcome model computed for itself. A model fitted as though its
#' weights were fixed understates its uncertainty, because the weights were
#' estimated from the same data, so there is nothing to fall back on when no
#' corrected block is there. `vcov()` and `confint()` refuse the reading instead,
#' and name which of the two ways the block can be missing they met. An outcome
#' model that was never wrapped raises an error of class
#' `causalgenerics_no_conditional_vcov`, which the package that produced the
#' result answers by wrapping the outcome model with [new_ipw_model()] before it
#' builds the result. An outcome model that carries the wrapper class with no
#' covariance behind it raises an error of class
#' `causalgenerics_no_vcov_ipw_model`, since a model that is already wrapped has
#' nothing to gain from being wrapped again and the object itself is what is
#' wrong. All of these carry the general class `causalgenerics_no_vcov`. `coef()`
#' needs no such block and reports the coefficients either way.
#'
#' The block and the coefficients are paired by name, and the block is reported
#' in coefficient order whatever order it was attached in, so the variance read
#' beside a coefficient is that coefficient's. A block whose labels cannot be
#' paired with the coefficients raises an error of class
#' `causalgenerics_conditional_vcov_mismatch`. A block of another size, one
#' labelled with the parameter names of a stacked system, and a model whose
#' coefficients carry no names are the three ways the pairing fails. Reading such
#' a block by position instead would report the covariance of other parameters
#' under this model's coefficient names.
#'
#' `nobs()`, `df.residual()`, `weights()`, `model.frame()`, and `estimand()`
#' describe the fit rather than a surface of it, so they answer the same way in
#' either reading and take no `effects` argument.
#'
#' @param object An `ipw` object.
#' @param formula An `ipw` object. The `model.frame()` generic in \pkg{stats}
#'   names its first argument `formula`, and the method matches it.
#' @param x An `ipw` object. The `estimand()` generic names its first argument
#'   `x`, and the method matches it.
#' @param parm The rows to report an interval for, given either as the labels
#'   of the reported surface or as their positions in it. Missing means all of
#'   them.
#' @param level The confidence level. At the level the result stores, its own
#'   limits are returned; at any other level they are recomputed.
#' @param effects The reading to report, either `"marginal"` or
#'   `"conditional"`. `NULL`, the default, reports the reading the result
#'   records; any other value overrides it for the one call and leaves the
#'   result as it is.
#' @param ... Further arguments. These methods ignore them.
#'
#' @return
#' `coef()` returns a named numeric vector of effect estimates, one element per
#' row of the `estimates` frame, named by effect label. In the conditional
#' reading it returns the outcome model's coefficients as that model reports
#' them, names included.
#'
#' `vcov()` returns the square numeric covariance matrix of those estimates,
#' with the effect labels as dimnames on both margins, exactly as the fitting
#' package attached it. A result that records no such matrix raises an error of
#' class `causalgenerics_no_vcov` rather than returning one built from the
#' standard errors, which would report the effects as uncorrelated. In the
#' conditional reading it returns the corrected covariance of the outcome
#' model's coefficients, in the order `coef()` reports them whatever order it was
#' attached in. An outcome model that was never wrapped raises an error of class
#' `causalgenerics_no_conditional_vcov`, a wrapped model that no longer carries
#' the covariance raises an error of class `causalgenerics_no_vcov_ipw_model`,
#' and a covariance whose labels cannot be paired with the coefficients raises an
#' error of class `causalgenerics_conditional_vcov_mismatch`.
#'
#' `confint()` returns a matrix with one row per effect `parm` selects, in the
#' order `parm` gives them, and two columns holding the lower and upper limit.
#' The rows are named by effect label and the columns by the two tail
#' probabilities as percentages, the way the `confint()` methods in \pkg{stats}
#' name theirs. A character `parm` that names an effect the result does not
#' report raises an error of class `causalgenerics_invalid_argument`. In the
#' conditional reading the rows are the outcome model's coefficients, which
#' `parm` names and indexes in the same two ways, and the limits are the normal
#' ones built from the corrected covariance at every level, since the limits the
#' result stores belong to the effects the marginal reading reports.
#'
#' `coef()`, `vcov()`, and `confint()` raise an error of class
#' `causalgenerics_invalid_argument_effects` when `effects` names neither
#' reading, and when the result's own field does.
#'
#' `nobs()` returns a single integer, the number of observations the outcome
#' model was fitted on, which is the number the estimates were computed from and
#' not necessarily the number the weighting model saw.
#'
#' `df.residual()` returns the residual degrees of freedom the fitted variance
#' object reports, as an integer when that number is a whole one an integer can
#' hold and as the double the object gave when it is not. A fractional count is
#' what a penalized or smooth fit spends, and it comes back as it stands rather
#' than truncated to a count the fit did not spend. `NA_integer_` comes back
#' when the result records no fitted variance object or that object reports no
#' residual degrees of freedom.
#'
#' `weights()` returns the outcome model's weights as its model frame stores
#' them, so a concrete weight class such as `psw` comes back as itself, or
#' `NULL` when the outcome model was fitted unweighted.
#'
#' `model.frame()` returns the outcome model's model frame, which is the data
#' the reported estimates were computed from. The `(weights)` column a weighted
#' frame carries is deliberately dropped, since tooling that reads the columns
#' of a frame as the variables a model was fitted on, such as prediction and
#' averaging packages, would otherwise treat the estimation weights as one of
#' them; `weights()` is where those are reported. A model whose own method has
#' no frame to give raises that model's error rather than one of this package's.
#'
#' `estimand()` returns the estimand the result records, which is the one the
#' weights the estimates were computed under targeted. There is deliberately no
#' `estimand<-()` method for a result: assigning a new estimand would relabel
#' those estimates rather than recompute them, so a caller who writes the
#' assignment is told that no method exists.
#'
#' @seealso [new_ipw()] for the result class and the fields these methods read.
#'
#' @examples
#' dat <- data.frame(
#'   x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
#'   z = rep(c(0, 1), 10),
#'   y = rep(c(0, 1, 1, 0, 1), 4)
#' )
#'
#' ps <- fitted(glm(z ~ x, family = binomial(), data = dat))
#' wts <- ifelse(dat$z == 1, 1 / ps, 1 / (1 - ps))
#'
#' # Written out literally, in the shape the `ipw()` return contract documents.
#' estimates <- data.frame(
#'   effect = c("rd", "log(rr)"),
#'   estimate = c(0.199882, 0.560414),
#'   std.err = c(0.092425, 0.273519),
#'   z = c(2.1626, 2.0489),
#'   ci.lower = c(0.018732, 0.024326),
#'   ci.upper = c(0.381032, 1.096502),
#'   conf.level = 0.95,
#'   p.value = c(0.030570, 0.040470)
#' )
#'
#' # The covariance of the two effects, which a method attaches to the estimates
#' # it returns. Both are computed from the same weighted means, so the
#' # off-diagonal entry is far from zero.
#' attr(estimates, "ipw_vcov") <- matrix(
#'   c(0.008542, 0.022753, 0.022753, 0.074813),
#'   nrow = 2,
#'   dimnames = list(estimates$effect, estimates$effect)
#' )
#'
#' res <- new_ipw(
#'   estimand = "ate",
#'   wt_mod = glm(z ~ x, family = binomial(), data = dat),
#'   outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat, weights = wts),
#'   estimates = estimates,
#'   se_method = "linearization",
#'   fit = NULL
#' )
#'
#' coef(res)
#' vcov(res)
#'
#' # At the stored level the stored limits come back.
#' confint(res)
#'
#' # At any other level they are recomputed.
#' confint(res, parm = "rd", level = 0.9)
#'
#' nobs(res)
#'
#' # The linearization path records no fitted variance object.
#' df.residual(res)
#'
#' head(weights(res))
#'
#' # The data the estimates were computed from, without the `(weights)` column
#' # the fitting machinery stored beside the variables the formula named.
#' head(model.frame(res))
#'
#' estimand(res)
#'
#' # The conditional reading reports the outcome model's coefficient surface
#' # rather than the effects.
#' coef(res, effects = "conditional")
#'
#' # Its covariance is the corrected block a fitting package attaches to the
#' # outcome model. This one carries none, and the covariance the model computed
#' # for itself is not a substitute: it treats the estimated weights as fixed.
#' try(vcov(res, effects = "conditional"))
#'
#' outcome_mod <- glm(y ~ z, family = quasibinomial(), data = dat, weights = wts)
#'
#' # Standing in for the outcome block of a stacked sandwich, which a method
#' # computes from the estimating equations of both steps.
#' corrected <- vcov(outcome_mod) * 1.4
#'
#' conditional <- as_conditional(new_ipw(
#'   estimand = "ate",
#'   wt_mod = glm(z ~ x, family = binomial(), data = dat),
#'   outcome_mod = new_ipw_model(outcome_mod, corrected),
#'   estimates = estimates,
#'   se_method = "linearization",
#'   fit = NULL
#' ))
#'
#' # A result that records the conditional reading answers in it with nothing
#' # named at the call site.
#' coef(conditional)
#' vcov(conditional)
#' confint(conditional)
#'
#' @name ipw-accessors
NULL

#' @rdname ipw-accessors
#' @export
#' @importFrom stats coef setNames
coef.ipw <- function(object, ..., effects = NULL) {
  # The coefficient surface is the outcome model's own, so the wrapper a
  # corrected covariance travels on is beside the point here: `coef()` walks
  # past it either way.
  if (resolve_ipw_effects(object, effects) == "conditional") {
    return(stats::coef(object$outcome_mod))
  }

  stats::setNames(
    object$estimates$estimate,
    ipw_effect_labels(object$estimates)
  )
}

#' @rdname ipw-accessors
#' @export
#' @importFrom stats vcov
vcov.ipw <- function(object, ..., effects = NULL) {
  if (resolve_ipw_effects(object, effects) == "conditional") {
    return(conditional_vcov(object$outcome_mod))
  }

  # The attribute is what the contract names. A covariance the fitted variance
  # object happens to carry is of its own parameters, which are not the effects
  # reported here, so there is nothing in `fit` to fall back on.
  covariance <- attr(object$estimates, "ipw_vcov", exact = TRUE)
  if (is.null(covariance)) {
    stop_no_vcov("ipw")
  }
  covariance
}

#' @rdname ipw-accessors
#' @export
#' @importFrom stats confint qnorm
confint.ipw <- function(object, parm, level = 0.95, ..., effects = NULL) {
  if (resolve_ipw_effects(object, effects) == "conditional") {
    # The covariance first, so that a result with no corrected block is refused
    # before a `parm` is matched against a surface no interval can be built for.
    covariance <- conditional_vcov(object$outcome_mod)
    estimate <- stats::coef(object$outcome_mod)
    labels <- names(estimate)
    rows <- if (missing(parm)) {
      seq_along(labels)
    } else {
      select_effects(
        parm,
        labels,
        surface = "coefficients the outcome model reports"
      )
    }

    # No stored limits to prefer: the ones the frame holds belong to the effects
    # the other reading reports, so every level is computed here.
    half_width <- stats::qnorm(1 - (1 - level) / 2) * sqrt(diag(covariance))
    limits <- cbind(
      (estimate - half_width)[rows],
      (estimate + half_width)[rows]
    )
    dimnames(limits) <- list(labels[rows], percent_labels(level))
    return(limits)
  }

  estimates <- object$estimates
  labels <- ipw_effect_labels(estimates)
  rows <- if (missing(parm)) {
    seq_along(labels)
  } else {
    select_effects(parm, labels)
  }

  # One half width, added and subtracted. `qnorm()` is not exactly
  # antisymmetric, so taking the lower limit from the lower tail instead would
  # put the two bounds a bit apart from each other.
  half_width <- stats::qnorm(1 - (1 - level) / 2) * estimates$std.err
  lower <- estimates$estimate - half_width
  upper <- estimates$estimate + half_width

  # At the level a row was reported for, the reported limits are returned rather
  # than rebuilt. They need not be the normal-based pair at all: a bootstrap or
  # profile interval is asymmetric about the estimate, and even a normal one
  # that was rounded on its way into the frame is not the number recomputing
  # gives. `which()` rather than a logical index so that a missing `conf.level`
  # selects nothing instead of raising a subscript error.
  stored <- which(estimates$conf.level == level)
  lower[stored] <- estimates$ci.lower[stored]
  upper[stored] <- estimates$ci.upper[stored]

  limits <- cbind(lower[rows], upper[rows])
  dimnames(limits) <- list(labels[rows], percent_labels(level))
  limits
}

#' @rdname ipw-accessors
#' @export
#' @importFrom stats nobs
nobs.ipw <- function(object, ...) {
  as.integer(stats::nobs(object$outcome_mod))
}

#' @rdname ipw-accessors
#' @export
#' @importFrom stats df.residual
df.residual.ipw <- function(object, ...) {
  if (is.null(object$fit)) {
    return(NA_integer_)
  }

  # A variance object is whatever the fitting package records, so the generic
  # may have nothing for it. `df.residual.default()` gives `NULL` for a bare
  # list, and a registered method is free to answer with an empty vector or
  # `Inf` for a fit with no residual degrees of freedom left. Each of those has
  # to become `NA_integer_` rather than a warning from `as.integer()`.
  df <- stats::df.residual(object$fit)
  if (length(df) != 1L || !is.numeric(df) || !is.finite(df)) {
    return(NA_integer_)
  }

  # A whole number comes back as an integer, the type a count of parameters
  # spent one at a time has. Anything else is the fit's own answer: a fractional
  # count is what a penalized or smooth term spends, and truncating it would
  # report a fit that spent more parameters than it did. A whole number past the
  # largest integer stays a double for the reason `as.integer()` cannot take it,
  # which is that the conversion gives `NA` and a warning.
  if (df == trunc(df) && abs(df) <= .Machine$integer.max) as.integer(df) else df
}

#' @rdname ipw-accessors
#' @export
#' @importFrom stats model.frame model.weights weights
weights.ipw <- function(object, ...) {
  stats::model.weights(stats::model.frame(object$outcome_mod))
}

#' @rdname ipw-accessors
#' @export
#' @importFrom stats model.frame
model.frame.ipw <- function(formula, ...) {
  # Delegation and nothing else, so a model with no frame to give says so
  # itself. The `(weights)` column is the one thing taken off: the formula never
  # named it, and tooling that reads the columns of a frame as the variables a
  # model was fitted on would take the estimation weights for one of them.
  mf <- stats::model.frame(formula$outcome_mod)
  mf[, setdiff(colnames(mf), "(weights)"), drop = FALSE]
}

#' @rdname ipw-accessors
#' @export
estimand.ipw <- function(x, ...) {
  x$estimand
}

#' The covariance the conditional reading reports
#'
#' The conditional reading presents the outcome model's coefficients, and their
#' covariance is the block of the joint estimation that a fitting package
#' attached with [new_ipw_model()]. A model that was not wrapped carries no such
#' block, and its own covariance is not a substitute: it treats the weights as
#' fixed when they were estimated from the same data, so it reports an
#' uncertainty the coefficients do not have. The guard therefore runs before the
#' model is asked anything, since asking is the mistake.
#'
#' The block and the coefficients are paired by name, which is what
#' [new_ipw_model()] leaves to be done here: the constructor checks that both
#' margins are labelled rather than what the labels say, since only the fitting
#' package knows which block of its stacked system belongs to which model. A
#' block labelled with the same names in another order is reported in
#' coefficient order, so the variance a caller reads beside a coefficient is that
#' coefficient's. Anything the labels cannot be paired with is refused: a block
#' of another size, one naming other parameters, and a model whose coefficients
#' carry no names at all. Falling back on the positions in any of those cases
#' would report the covariance of something else under this model's coefficient
#' names, which is the one answer a caller has no way to check.
#'
#' @param model The result's `outcome_mod`.
#' @param call The call to report the error against, which is the accessor's
#'   rather than this helper's.
#'
#' @return The corrected covariance matrix, in coefficient order.
#'
#' @noRd
conditional_vcov <- function(model, call = sys.call(-1)) {
  if (!inherits(model, "ipw_model")) {
    stop_no_conditional_vcov(call = call)
  }
  covariance <- stats::vcov(model)

  labels <- names(stats::coef(model))
  block_labels <- rownames(covariance)

  # The size is checked alongside the names because a larger block can hold
  # every coefficient name and more. Indexing that one by name alone would take
  # the corner of it and answer with the covariance of two coefficients of a
  # different fit.
  paired <- !is.null(labels) &&
    identical(dim(covariance), c(length(labels), length(labels))) &&
    setequal(block_labels, labels) &&
    setequal(colnames(covariance), labels)

  if (!paired) {
    stop_conditional_vcov_mismatch(block_labels, labels, call = call)
  }

  # A block already in coefficient order comes back as the fitting package
  # attached it, rather than through an indexing that would rebuild the same
  # matrix.
  if (
    identical(block_labels, labels) && identical(colnames(covariance), labels)
  ) {
    return(covariance)
  }

  # Both margins, since a matrix reordered down one and not the other is no
  # longer the covariance of anything.
  covariance[labels, labels, drop = FALSE]
}

#' The rows of a reported surface a `parm` argument selects
#'
#' Character `parm` matches the row labels and numeric `parm` indexes them,
#' which is the split every `confint()` method in \pkg{stats} makes. Either way,
#' a selection that names a row the surface does not have is refused rather than
#' dropped: returning the rows that did match would answer a question the caller
#' did not ask, with a matrix of the wrong number of rows and nothing to say
#' why.
#'
#' The surface names itself in the message. A caller who asks the conditional
#' reading for an effect label, or the marginal one for a coefficient name, has
#' the labels of the surface they were reading in hand, and the error is where
#' they find out which surface they are addressing.
#'
#' @param parm The `parm` argument as the caller supplied it.
#' @param labels The reported surface's row labels.
#' @param surface What those labels name, completing the sentence "`parm` must
#'   name ..." and "`parm` must index ...".
#' @param call The call to report the error against, which is the accessor's
#'   rather than this helper's.
#'
#' @return An integer or numeric vector of row positions.
#'
#' @noRd
select_effects <- function(
  parm,
  labels,
  surface = "effects the result reports",
  call = sys.call(-1)
) {
  if (is.character(parm)) {
    rows <- match(parm, labels)
    if (anyNA(rows)) {
      stop_invalid_argument(
        "parm",
        paste0(
          "name ",
          surface,
          ", which are ",
          toString(encodeString(labels, quote = '"'))
        ),
        call = call
      )
    }
    return(rows)
  }

  # Indexing the positions rather than the labels keeps the ordinary meaning of
  # a numeric subscript, negative positions included, and turns anything outside
  # the surface into an `NA` to catch.
  rows <- seq_along(labels)[parm]
  if (anyNA(rows)) {
    stop_invalid_argument(
      "parm",
      paste0(
        "index ",
        surface,
        ", numbered 1 to ",
        length(labels)
      ),
      call = call
    )
  }
  rows
}

#' Label a confidence interval's columns the way \pkg{stats} does
#'
#' The `confint()` methods in \pkg{stats} label the two columns with the tail
#' probabilities as percentages, and a caller reading a limit out by column name
#' is entitled to the same labels here. Those methods share an internal helper
#' for the formatting, which is not ours to call, so the formatting is restated.
#'
#' @param level The confidence level.
#'
#' @return A character vector of length two.
#'
#' @noRd
percent_labels <- function(level) {
  tail <- (1 - level) / 2
  paste(
    format(
      100 * c(tail, 1 - tail),
      trim = TRUE,
      scientific = FALSE,
      digits = 3
    ),
    "%"
  )
}
