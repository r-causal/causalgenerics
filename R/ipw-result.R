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
#' path, including the ones that have no fitted variance object to report.
#'
#' `as.data.frame()` returns the `estimates` component. With
#' `exponentiate = TRUE` it moves the `log(rr)` and `log(or)` rows to their
#' natural scale, exponentiating the point estimate and the confidence limits and
#' relabelling the two effects `"rr"` and `"or"`. Standard errors, z statistics,
#' and p-values stay on the log scale, where the inference is done.
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
#'
#' @return `new_ipw()` returns an S3 object of class `ipw`: a list of the
#'   following six components, in this order.
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
#' }
#'
#'   `print()` returns its input invisibly. `as.data.frame()` returns the
#'   `estimates` component as a data frame.
#'
#' @seealso [ipw()], the generic these results come from.
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
new_ipw <- function(estimand, wt_mod, outcome_mod, estimates, se_method, fit) {
  structure(
    list(
      estimand = estimand,
      wt_mod = wt_mod,
      outcome_mod = outcome_mod,
      estimates = estimates,
      se_method = se_method,
      fit = fit
    ),
    class = "ipw"
  )
}

#' @param x An `ipw` object.
#' @param ... Further arguments. `print()` ignores them; `as.data.frame()`
#'   passes them to [base::as.data.frame()].
#' @rdname new_ipw
#' @export
#' @importFrom stats printCoefmat
print.ipw <- function(x, ...) {
  cat("Inverse Probability Weight Estimator\n")
  cat("Estimand:", toupper(x$estimand), "\n\n")

  cat("Weight Estimator:\n")
  cat("  Call:", format_model_call(x$wt_mod), "\n")
  cat("\n")

  cat("Outcome Model:\n")
  cat("  Call:", format_model_call(x$outcome_mod), "\n")

  cat("\n")

  cat("Estimates:\n")
  if ("comparison" %in% names(x$estimates)) {
    # A categorical exposure repeats effect labels across comparisons, so the
    # printed rows are keyed by effect and comparison together and the character
    # comparison column is dropped from the numeric matrix printCoefmat formats.
    estimates <- x$estimates[setdiff(
      names(x$estimates),
      c("effect", "comparison")
    )]
    rownames(estimates) <- paste(x$estimates$effect, x$estimates$comparison)
  } else {
    estimates <- x$estimates[-1]
    rownames(estimates) <- x$estimates$effect
  }
  stats::printCoefmat(estimates, has.Pvalue = TRUE, cs.ind = 1:2, tst.ind = 3)

  invisible(x)
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
