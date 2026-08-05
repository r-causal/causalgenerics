#' Attach a corrected covariance to a component model of an IPW fit
#'
#' @description
#' `new_ipw_model()` wraps one of the models an [ipw()] method was given so that
#' the model reports the covariance the two-step estimation implies rather than
#' the one it computed for itself. It is intended for package developers writing
#' an `ipw()` method, not for end users, and beyond the shape of `vcov` it
#' assumes its arguments are already validated.
#'
#' @details
#' A weighted outcome model fitted as though its weights were fixed understates
#' its own uncertainty, because the weights were estimated from the same data. A
#' method that solves the estimating equations of both steps as one stacked
#' system has the corrected block on hand, and wrapping the model is how that
#' block reaches a caller who writes `vcov(result$outcome_mod)`.
#'
#' The wrapper is two changes and no more: `"ipw_model"` is prepended to the
#' model's class vector, and the matrix is attached as the `ipw_vcov` attribute.
#' Nothing is registered for the class beyond `vcov()`, so `predict()`, `coef()`,
#' `model.frame()`, and every other generic walk past it to the model's own
#' methods. That is what lets the corrected covariance travel with the model
#' rather than beside it.
#'
#' The dimnames of `vcov` are the fitting package's to choose. They are
#' ordinarily `names(coef(model))`, but a method whose stacked system is
#' parameterized some other way labels the block with its own parameter names.
#' The constructor checks that both margins are labelled rather than what the
#' labels say, since only the fitting package knows which block of the sandwich
#' belongs to which model.
#'
#' @param model The fitted component model to wrap, such as the outcome model or
#'   the weighting model of an [ipw()] result.
#' @param vcov The corrected covariance of that model's parameters: a square
#'   numeric matrix with dimnames on both margins.
#' @param object An `ipw_model` object.
#' @param ... Further arguments. `vcov()` ignores them.
#'
#' @return `new_ipw_model()` returns `model` with `"ipw_model"` prepended to its
#'   class vector and `vcov` attached as the `ipw_vcov` attribute. Everything
#'   else about the model is untouched. `vcov()` returns that matrix.
#'
#' @seealso [new_ipw()] for the result these models are components of, and
#'   [ipw-accessors] for the accessors on the result itself.
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
#' ps <- fitted(glm(z ~ x, family = binomial(), data = dat))
#' wts <- ifelse(dat$z == 1, 1 / ps, 1 / (1 - ps))
#' outcome_mod <- glm(y ~ z, family = quasibinomial(), data = dat, weights = wts)
#'
#' # Standing in for the outcome block of a stacked sandwich. A real method
#' # computes this from the estimating equations of both steps; treating the
#' # weights as estimated rather than fixed ordinarily widens it.
#' corrected <- vcov(outcome_mod) * 1.4
#'
#' wrapped <- new_ipw_model(outcome_mod, corrected)
#'
#' class(wrapped)
#'
#' # The corrected covariance rather than the model's own.
#' vcov(wrapped)
#'
#' # Everything else about the model reaches its own methods.
#' coef(wrapped)
new_ipw_model <- function(model, vcov) {
  if (!is.matrix(vcov) || !is.numeric(vcov)) {
    stop_invalid_argument("vcov", "be a square numeric matrix")
  }
  if (nrow(vcov) != ncol(vcov)) {
    stop_invalid_argument(
      "vcov",
      paste0("be square, but it is ", nrow(vcov), " by ", ncol(vcov))
    )
  }
  if (is.null(rownames(vcov)) || is.null(colnames(vcov))) {
    stop_invalid_argument("vcov", "have dimnames on both margins")
  }

  attr(model, "ipw_vcov") <- vcov
  class(model) <- c("ipw_model", class(model))
  model
}

#' @rdname new_ipw_model
#' @export
#' @importFrom stats vcov
vcov.ipw_model <- function(object, ...) {
  attr(object, "ipw_vcov", exact = TRUE)
}
