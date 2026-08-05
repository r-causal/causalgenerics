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
#' * `df.residual()` returns the residual degrees of freedom of the fitted
#'   variance object.
#' * `weights()` returns the weights the outcome model was fitted with.
#'
#' @details
#' These methods live here for the reason `print()` does. Two packages each
#' registering `coef.ipw()` would collide in the shared S3 method table, and a
#' caller writing against a result would then get whichever package was
#' installed last rather than the contract.
#'
#' Everything they read is part of the [new_ipw()] contract: the `estimates`
#' frame, the `ipw_vcov` attribute attached to it, `outcome_mod`, and `fit`.
#' They never branch on `se_method`, since the fields they read already hold the
#' result of whatever computation that names, and they reach into `fit` only
#' through ordinary S3 dispatch. A package whose variance object is a bare list
#' rather than a fitted model therefore uses them unchanged.
#'
#' Rows are named by the effect labels the [new_ipw()] contract defines, so the
#' name `coef()` gives an estimate is the one `vcov()` and `confint()` use for
#' it and the one `print()` labels its row with.
#'
#' @param object An `ipw` object.
#' @param parm The effects to report an interval for, given either as effect
#'   labels or as their positions in the result. Missing means all of them.
#' @param level The confidence level. At the level the result stores, its own
#'   limits are returned; at any other level they are recomputed.
#' @param ... Further arguments. These methods ignore them.
#'
#' @return
#' `coef()` returns a named numeric vector of effect estimates, one element per
#' row of the `estimates` frame, named by effect label.
#'
#' `vcov()` returns the square numeric covariance matrix of those estimates,
#' with the effect labels as dimnames on both margins, exactly as the fitting
#' package attached it. A result that records no such matrix raises an error of
#' class `causalgenerics_no_vcov` rather than returning one built from the
#' standard errors, which would report the effects as uncorrelated.
#'
#' `confint()` returns a matrix with one row per effect `parm` selects, in the
#' order `parm` gives them, and two columns holding the lower and upper limit.
#' The rows are named by effect label and the columns by the two tail
#' probabilities as percentages, the way the `confint()` methods in \pkg{stats}
#' name theirs. A character `parm` that names an effect the result does not
#' report raises an error of class `causalgenerics_invalid_argument`.
#'
#' `nobs()` returns a single integer, the number of observations the outcome
#' model was fitted on, which is the number the estimates were computed from and
#' not necessarily the number the weighting model saw.
#'
#' `df.residual()` returns a single integer, or `NA_integer_` when the result
#' records no fitted variance object or that object reports no residual degrees
#' of freedom.
#'
#' `weights()` returns the outcome model's weights as its model frame stores
#' them, so a concrete weight class such as `psw` comes back as itself, or
#' `NULL` when the outcome model was fitted unweighted.
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
#' @name ipw-accessors
NULL

#' @rdname ipw-accessors
#' @export
#' @importFrom stats coef setNames
coef.ipw <- function(object, ...) {
  stats::setNames(
    object$estimates$estimate,
    ipw_effect_labels(object$estimates)
  )
}

#' @rdname ipw-accessors
#' @export
#' @importFrom stats vcov
vcov.ipw <- function(object, ...) {
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
confint.ipw <- function(object, parm, level = 0.95, ...) {
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
  # list, and a registered method is free to answer with a double, an empty
  # vector, or `Inf` for a fit with no residual degrees of freedom left. Each of
  # those has to become `NA_integer_` rather than a warning from `as.integer()`.
  df <- stats::df.residual(object$fit)
  if (length(df) != 1L || !is.numeric(df) || !is.finite(df)) {
    return(NA_integer_)
  }
  as.integer(df)
}

#' @rdname ipw-accessors
#' @export
#' @importFrom stats model.frame model.weights weights
weights.ipw <- function(object, ...) {
  stats::model.weights(stats::model.frame(object$outcome_mod))
}

#' The rows of an `ipw` result a `parm` argument selects
#'
#' Character `parm` matches the effect labels and numeric `parm` indexes them,
#' which is the split every `confint()` method in \pkg{stats} makes. Either way,
#' a selection that names a row the result does not have is refused rather than
#' dropped: returning the rows that did match would answer a question the caller
#' did not ask, with a matrix of the wrong number of rows and nothing to say
#' why.
#'
#' @param parm The `parm` argument as the caller supplied it.
#' @param labels The result's effect labels.
#' @param call The call to report the error against, which is the accessor's
#'   rather than this helper's.
#'
#' @return An integer or numeric vector of row positions.
#'
#' @noRd
select_effects <- function(parm, labels, call = sys.call(-1)) {
  if (is.character(parm)) {
    rows <- match(parm, labels)
    if (anyNA(rows)) {
      stop_invalid_argument(
        "parm",
        paste0(
          "name effects the result reports, which are ",
          toString(encodeString(labels, quote = '"'))
        ),
        call = call
      )
    }
    return(rows)
  }

  # Indexing the positions rather than the labels keeps the ordinary meaning of
  # a numeric subscript, negative positions included, and turns anything outside
  # the result into an `NA` to catch.
  rows <- seq_along(labels)[parm]
  if (anyNA(rows)) {
    stop_invalid_argument(
      "parm",
      paste0(
        "index effects the result reports, numbered 1 to ",
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
