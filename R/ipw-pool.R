#' Pool inverse probability weighted results over multiply imputed data
#'
#' @description
#' `pool_ipw()` combines the [ipw()] results fitted to each of a set of multiply
#' imputed datasets into a single result, by Rubin's rules. It takes a plain
#' list of results or the `mira` object `mice::with()` returns, and gives back
#' an object of class `ipw_pooled`.
#'
#' @details
#' Each result is estimated on one completed dataset and carries a standard
#' error that accounts for the weights having been estimated, and pooling those
#' by Rubin's rules adds the uncertainty the imputation itself contributed: the
#' order matters, since a standard error that treated the weights as fixed would
#' be too small before the pooling ever saw it, and pooling the imputed datasets
#' rather than the estimates would understate the uncertainty whatever the
#' standard errors were.
#'
#' For each effect the pooled estimate is the mean of the per-imputation
#' estimates, and the pooled variance is the mean of their squared standard
#' errors plus the between-imputation variance inflated by `1 + 1/m` (Rubin,
#' 1987). The degrees of freedom carry the Barnard-Rubin small-sample adjustment
#' (Barnard and Rubin, 1999), which is why the complete-data count matters: with
#' few imputations the pooled degrees of freedom can be far below what a single
#' analysis reports, and the interval and the p-value are referred to t on them
#' rather than to the normal.
#'
#' `mice` is not a dependency, and nothing here is imported from it. A `mira` is
#' recognized by its class and read through its `analyses` element, which is all
#' the object is, so a package that produces one by another route is answered
#' the same way.
#'
#' # The two readings
#'
#' Both readings of the results are pooled from one call whenever both can be
#' pooled. `effects` says which of them the returned result presents, and the
#' other is stored whole under `alternate`, so a caller moves a pooled result
#' between the two afterwards with [as_marginal()] and [as_conditional()] rather
#' than pooling again. Which reading was named therefore decides which pair of
#' frames is the active one rather than what was computed at all: pooling one
#' reading and moving to the other gives the frames pooling the other directly
#' gives.
#'
#' The reading the call names is pooled first and is not guarded, so a set that
#' cannot be pooled on it is refused as it always was. The other reading is
#' pooled under a guard, since a set that cannot be pooled on it is a result with
#' one reading rather than a failed pooling: nothing about the reading the caller
#' asked for is wrong. `alternate` then records that reading and the refusal it
#' raised, and asking the result for it later raises that refusal in the words it
#' used. The commonest case is a set of results whose outcome models carry no
#' corrected covariance, which is what the conditional reading is pooled from.
#'
#' The readings the results declare are read before either surface is, so a set
#' whose results support one reading pools that one and no other. Naming the
#' reading they do not support refuses the call with an error of class
#' `causalgenerics_unsupported_reading_marginal` or
#' `causalgenerics_unsupported_reading_conditional`, in the words a result raises
#' when the same reading is asked of it directly and a sentence saying which
#' result of the set that was, since what the caller asked for is a reading of
#' the analysis rather than a way of pooling it. Leaving `effects` at `NULL`
#' pools the reading the results record and puts that same refusal on
#' `alternate` as the reason the other reading has no frames, so
#' [as_marginal()] or [as_conditional()] on the pooled result afterwards raises
#' an error of class `causalgenerics_pool_missing_surface_marginal` or
#' `causalgenerics_pool_missing_surface_conditional` carrying it. Every result
#' has to support the reading being pooled, since the pooled estimate of an
#' effect is an average over all of them, so a set where one result alone lacks
#' it is refused for that result. The caller has to go back to that one result,
#' so the refusal names its position in `fits` in a sentence of its own and
#' carries the position under `position`.
#'
#' The components describing the analyses rather than a reading of them are
#' shared by both readings. The estimand, the standard error method, the number
#' of results, the complete-data degrees of freedom, the observation count, and
#' the outcome model link are settled once, from the results themselves, and are
#' what either reading reports.
#'
#' # What the results have to agree on
#'
#' The pooled estimate of an effect is an average of the per-imputation ones, so
#' the results have to be estimating the same things the same way. A set that
#' disagrees about its `estimand`, its `se_method`, the presentation mode it
#' records, the effects it reports, the confidence level it stores, or the link
#' its outcome model was fitted with is refused with an error of class
#' `causalgenerics_pool_mismatch`, and of a second class naming which of those
#' it was. The differing values travel on the condition under `values`.
#'
#' Those requirements are not all of one kind. The estimand, the standard error
#' method, and the outcome model link describe the results rather than a reading
#' of them, and so does the presentation mode when `effects` leaves it to be
#' read; a disagreement about any of those refuses the call. Which effects a
#' result reports and what level it reported them at are properties of one
#' reading of it, and those requirements bind the reading being pooled: a
#' disagreement confined to the reading the call did not name is recorded on
#' `alternate` as the reason that reading could not be pooled rather than
#' refused.
#'
#' The effects have to agree as an ordered vector rather than as a set. The
#' labels are what say which row is which, so two results reporting the same
#' contrasts in different orders would otherwise have the `b vs a` rows of one
#' averaged with the `c vs a` rows of the other.
#'
#' A label carries the subgroup a row was estimated in along with the effect and
#' the contrast, so a set whose results report different subgroups, or the same
#' subgroups in different orders, disagrees about its labels and is refused
#' through that requirement rather than through one of its own. A subgroup
#' observed in one imputation and not in another is the case this answers: there
#' is no row to average the missing one with, and pooling by position would
#' combine two subgroups that answer different questions.
#'
#' Two of those agreements are only required when the argument that would settle
#' the question is left at `NULL`. Naming `effects` says which surface to pool
#' and the stored mode is not read at all; naming `conf_level` says what the
#' bounds report and the stored level is not read at all.
#'
#' # The complete-data degrees of freedom
#'
#' `dfcom` is the number of observations the analysis had minus the number of
#' parameters it fitted, before any data went missing. It is looked for in three
#' places, in order:
#'
#' 1. The `dfcom` argument, when it is not `NULL`.
#' 2. `df.residual()` on each result, which reports what the fitted variance
#'    object records, when at least one of them reports a number.
#' 3. `df.residual()` on each result's outcome model, when at least one of those
#'    reports a number.
#'
#' The smallest count is taken rather than the first, since the pooled inference
#' is no stronger than the weakest analysis supports, and results that report
#' nothing are passed over rather than making the minimum missing. Wherever the
#' count comes from it is at least 1, since the adjustment has nothing to work
#' with at zero.
#'
#' When nothing reports a count, `dfcom` is `Inf` and a warning of class
#' `causalgenerics_pool_large_sample` says so. That is the honest answer for a
#' large sample and the widest degrees of freedom the adjustment can give, which
#' makes it the narrowest intervals, so it is said out loud rather than assumed
#' quietly.
#'
#' # The scale the effects are pooled on
#'
#' The estimates are pooled as they are stored. A result reports its ratio
#' effects on the log scale, under the labels `log(rr)` and `log(or)`, which is
#' the scale the inference is done on and the only one on which averaging the
#' per-imputation estimates means anything. Nothing here exponentiates, and the
#' pooled frame carries the same labels the results did. Moving the ratios to
#' their natural scale is presentation rather than pooling, and belongs where
#' the same choice is made for a single result.
#'
#' @param fits The results to pool: a list of `ipw` objects, one per imputed
#'   dataset, or a `mira`. At least two are needed, since the
#'   between-imputation variance is estimated from the spread across them.
#' @param ... These dots exist so that every argument after `fits` is matched by
#'   name. Passing anything through them is an error.
#' @param effects The reading the pooled result presents, either `"marginal"` or
#'   `"conditional"`. `NULL`, the default, presents the reading the results
#'   record. Both readings are pooled whenever both can be, so this says which
#'   one the result reports and which one it stores beside it. The marginal
#'   reading pools the causal contrast estimates; the conditional reading pools
#'   the outcome models' coefficients, with the standard errors implied by the
#'   corrected covariance each one carries.
#' @param dfcom The complete-data degrees of freedom. `NULL`, the default, reads
#'   them off the results as the section above describes.
#' @param conf_level The level the pooled bounds report. `NULL`, the default,
#'   uses the level the results stored, or `0.95` when they store none.
#'
#' @return An S3 object of class `ipw_pooled`: a list of the following ten
#'   components, in this order.
#' \describe{
#'   \item{`estimand`}{The causal estimand every pooled result targeted.}
#'   \item{`estimates`}{A data frame with one row per effect and the following
#'     columns: `effect` (the measure name), `contrast` after it when the
#'     results name contrasts, `group` after that when they name subgroups,
#'     `estimate` (the pooled point estimate),
#'     `std.err` (the pooled standard error), `t` (the test statistic),
#'     `df` (the pooled degrees of freedom), `ci.lower` and `ci.upper`,
#'     `conf.level`, and `p.value`. The statistic and the p-value are referred
#'     to t on `df` rather than to the normal. When every pooled result carried
#'     the covariance of its effects, the pooled covariance is attached to this
#'     frame as the `ipw_vcov` attribute, with the effect labels as dimnames on
#'     both margins; when any of them carried none, no attribute is attached,
#'     since a matrix built from a subset of the imputations would sit beside
#'     estimates built from all of them.}
#'   \item{`pooling`}{A data frame keyed by the same columns that name a row of
#'     the frame above, holding `ubar` (the within-imputation variance), `b` (the
#'     between-imputation variance), `riv` (the relative increase in variance),
#'     `lambda` (the proportion of the total variance due to missingness), and
#'     `fmi` (the fraction of missing information).}
#'   \item{`se_method`}{The standard error method every pooled result used.}
#'   \item{`effects`}{The reading that was pooled, either `"marginal"` or
#'     `"conditional"`.}
#'   \item{`m`}{The number of results pooled.}
#'   \item{`dfcom`}{The complete-data degrees of freedom the adjustment used.}
#'   \item{`nobs`}{The smallest number of observations any pooled result was
#'     estimated from.}
#'   \item{`outcome_link`}{The link every pooled result's outcome model was
#'     fitted with, which is the scale the effects are reported on.}
#'   \item{`alternate`}{The reading the result does not present, in one of two
#'     shapes. When that reading was pooled it is a list of `effects`, naming
#'     which reading it is, and `estimates` and `pooling`, the two frames
#'     described above built for it and built the same way, the pooled covariance
#'     on the estimates frame included. When it could not be pooled it is a list
#'     of `effects` and `reason`, the message of the refusal that reading raised.
#'     [as_marginal()] and [as_conditional()] read this component to move the
#'     result between the two readings, and the accessors that take an `effects`
#'     argument read it the same way.}
#' }
#'
#' @references
#' Barnard, J. and Rubin, D. B. (1999). Small sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#'
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*. New
#' York: John Wiley and Sons.
#'
#' @seealso [new_ipw()] for the results this pools and the fields it reads, and
#'   [as_marginal()] and [as_conditional()] for moving the pooled result between
#'   the two readings it carries.
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
#' # Three imputations of that dataset, differing in one filled-in cell each.
#' imputed <- lapply(1:3, function(i) {
#'   completed <- dat
#'   completed$y[i] <- 1 - completed$y[i]
#'   completed
#' })
#'
#' # The estimates each analysis produced, written out literally in the shape
#' # the `ipw()` return contract documents. What is being shown is the pooling,
#' # so these stand in for what a method would compute from the models beside
#' # them.
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
#' # The pooled effects. The standard errors are wider than the average of the
#' # per-imputation ones, which is the uncertainty the imputation contributed.
#' pooled$estimates
#'
#' # How much of that uncertainty came from the imputation rather than the data.
#' pooled$pooling
#'
#' # The complete-data count the adjustment used, read off the outcome models.
#' pooled$dfcom
#'
#' # Naming one settles it directly, and the pooled degrees of freedom move.
#' pool_ipw(fits, dfcom = 500)$estimates$df
#'
#' # The ratio effect stays on the log scale it was estimated on.
#' pooled$estimates$effect
#'
#' # Both readings are pooled from one call. These outcome models carry no
#' # corrected covariance, which is what the conditional reading is pooled from,
#' # so that reading records why it has none rather than refusing the call.
#' pooled$alternate$effects
#'
#' try(as_conditional(pooled))
pool_ipw <- function(
  fits,
  ...,
  effects = NULL,
  dfcom = NULL,
  conf_level = NULL
) {
  # Taken here rather than left to each helper's own default. Several of the
  # conditions below are raised from inside an `lapply()`, and a call resolved
  # from there would name `FUN(X[[i]], ...)` instead of the call the caller
  # wrote. Reading it eagerly is what makes every refusal name `pool_ipw()`.
  call <- sys.call()

  check_pool_dots(list(...), call = call)
  if (!is.null(effects)) {
    check_ipw_effects(effects, call = call)
  }
  if (!is.null(dfcom)) {
    check_dfcom(dfcom, call = call)
  }
  if (!is.null(conf_level)) {
    check_conf_level(conf_level, "conf_level", call = call)
  }

  fits <- pool_analyses(fits, call = call)
  m <- length(fits)

  estimand <- pool_common(fits, function(x) x$estimand, "estimand", call)
  se_method <- pool_common(fits, function(x) x$se_method, "se_method", call)
  link <- pool_common(
    fits,
    function(x) ipw_outcome_link(x$outcome_mod),
    "outcome_link",
    call
  )

  # A named mode says which surface to pool, so the stored field is not read at
  # all in that case and a set that disagrees about it is poolable by naming
  # one. This is how `resolve_ipw_effects()` reads an accessor's argument.
  mode <- if (is.null(effects)) {
    pool_mode(fits, call = call)
  } else {
    effects
  }

  # The reading the caller asked for is read first and is never guarded, so a
  # set that cannot be pooled on it is refused rather than answered with an
  # empty result. The order the refusals come in is part of the contract: what
  # the results report has to agree before the level they report it at does.
  read <- pool_mode_surfaces(fits, mode, call = call)
  level <- pool_mode_level(read$surfaces, conf_level, call = call)

  # Settled once, before either reading is pooled, since the complete-data count
  # is a property of the analyses rather than of a reading of them. Resolving it
  # per surface would raise the large-sample warning twice for one call, which a
  # reader would take as two separate assumptions.
  if (is.null(dfcom)) {
    dfcom <- pool_dfcom(fits)
    if (is.na(dfcom)) {
      warn_pool_large_sample(call = call)
      dfcom <- Inf
    }
  }
  dfcom <- max(dfcom, 1)

  frames <- pool_mode_frames(read, level, dfcom)

  new_ipw_pooled(
    estimand = estimand,
    estimates = frames$estimates,
    pooling = frames$pooling,
    se_method = se_method,
    effects = mode,
    m = m,
    dfcom = dfcom,
    nobs = min(vapply(fits, stats::nobs, integer(1))),
    outcome_link = link,
    alternate = pool_alternate(fits, mode, conf_level, dfcom, call = call)
  )
}

#' The mode the results agree on
#'
#' The stored mode of each result, read with the result's position in the set,
#' and refused where the results do not agree about it. It is a helper rather
#' than a `pool_common()` call that reads the field with `ipw_effects()`,
#' because that reader refuses two things of its own: a field holding neither
#' reading, and a mode the result records no reading of. Invoked from inside
#' `pool_common()` it would resolve its call to the `lapply()` that invoked it,
#' and both refusals would name that rather than [pool_ipw()].
#'
#' Pairing each result with its position is what names the result a cross-check
#' refusal is about, the way `pool_mode_surfaces()` names the one a surface
#' refusal is about. A caller holding several results cannot tell them apart
#' from the refusal otherwise.
#'
#' @param fits The results being pooled.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return A single string, either `"marginal"` or `"conditional"`.
#'
#' @noRd
pool_mode <- function(fits, call = sys.call(-1)) {
  modes <- Map(
    function(fit, position) {
      ipw_effects(fit, position = position, call = call)
    },
    fits,
    seq_along(fits)
  )

  pool_common(modes, identity, "effects", call)
}

#' The surfaces of one reading, and what keys their rows
#'
#' The first half of pooling a reading: every result's surface is read, what
#' they report has to agree, and the columns that name a row are taken from the
#' first of them. It is a helper rather than part of [pool_ipw()] because both
#' readings are pooled from one call, and two copies of it would be two places
#' for the agreements to drift apart.
#'
#' The key is the columns that name a row, which both returned frames are keyed
#' by, and it is the key of the first surface rather than one rebuilt here. The
#' surfaces have agreed on their labels by the time it is taken, and the labels
#' are the key columns pasted together, so the first surface's key is every
#' surface's. A binary or continuous exposure has one contrast and an ungrouped
#' result has one group, so a column of either kind would repeat one value down
#' the table and read as a contrast or a subgroup that was named. The
#' conditional reading names coefficients and has neither.
#'
#' @param fits The results being pooled.
#' @param effects The reading to read off them.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return A list of the surfaces, the effect labels they agree on, and the key.
#'
#' @noRd
pool_mode_surfaces <- function(fits, effects, call = sys.call(-1)) {
  # Paired with its position rather than read on its own, so that a result which
  # does not support the reading is named by where it sits in the set. A caller
  # holding several results cannot tell them apart from the refusal otherwise.
  surfaces <- Map(
    function(fit, position) {
      pool_surface(fit, effects = effects, position = position, call = call)
    },
    fits,
    seq_along(fits)
  )
  labels <- pool_common(surfaces, function(x) x$labels, "labels", call)

  list(surfaces = surfaces, labels = labels, key = surfaces[[1L]]$key)
}

#' The level one reading's bounds are reported at
#'
#' The same reading of a `NULL` argument that the mode gets. A level named at
#' the call site settles the bounds of every reading, and otherwise each reading
#' resolves its own from what its own surfaces stored: an estimates frame
#' records the level its bounds were reported at and a coefficient surface
#' records none. A set that records no level at all leaves nothing to disagree
#' about, and the level every other surface in this package defaults to is what
#' the bounds report.
#'
#' @param surfaces The surfaces of one reading.
#' @param conf_level The argument as the caller supplied it.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return A single number.
#'
#' @noRd
pool_mode_level <- function(surfaces, conf_level, call = sys.call(-1)) {
  level <- if (is.null(conf_level)) {
    pool_common(surfaces, function(x) x$conf.level, "conf_level", call)
  } else {
    conf_level
  }

  if (is.null(level)) 0.95 else level
}

#' The two frames one reading is pooled into
#'
#' The second half of pooling a reading: Rubin's rules over the surfaces, and
#' the estimates and diagnostics frames those give, keyed the same way so that a
#' row of one is found from a row of the other. The pooled covariance rides on
#' the estimates frame it describes and is named the way that frame is named.
#'
#' @param read The output of `pool_mode_surfaces()`.
#' @param level The level the bounds report.
#' @param dfcom The complete-data degrees of freedom, which both readings share.
#'
#' @return A list of the estimates frame and the diagnostics frame.
#'
#' @noRd
pool_mode_frames <- function(read, level, dfcom) {
  surfaces <- read$surfaces
  estimate <- do.call(rbind, lapply(surfaces, function(x) x$estimate))
  std_err <- do.call(rbind, lapply(surfaces, function(x) x$std.err))
  pooled <- rubin_rules(estimate, std_err, dfcom)

  estimates <- pool_estimates_frame(read$key, pooled, level)
  attr(estimates, "ipw_vcov") <- pool_vcov(surfaces, estimate, read$labels)

  list(
    estimates = estimates,
    pooling = pool_diagnostics_frame(read$key, pooled)
  )
}

#' The reading the caller did not ask for
#'
#' Both readings are pooled from one set of results, so which one the call named
#' decides which pair of frames is the active one rather than what gets computed
#' at all. The other pair is recorded here, which is what lets a caller move a
#' pooled result between the two readings afterwards without pooling again.
#'
#' This pass is guarded and the active one is not. A set whose other reading
#' cannot be pooled is a result with one reading rather than a failed pooling,
#' since nothing about the reading the caller asked for is wrong; what stopped
#' the other one is recorded so that the request it answers can raise it later,
#' in the words the refusal itself used. The guard catches errors alone, and the
#' pass raises nothing else: the complete-data count is resolved before it and
#' passed in, so the only warning [pool_ipw()] has is already behind it.
#'
#' @param fits The results being pooled.
#' @param effects The reading that was pooled actively.
#' @param conf_level The `conf_level` argument as the caller supplied it.
#' @param dfcom The complete-data degrees of freedom the active reading used.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return A list naming the other reading, holding either its two frames or the
#'   reason it has none.
#'
#' @noRd
pool_alternate <- function(
  fits,
  effects,
  conf_level,
  dfcom,
  call = sys.call(-1)
) {
  other <- setdiff(c("marginal", "conditional"), effects)

  frames <- tryCatch(
    {
      read <- pool_mode_surfaces(fits, other, call = call)
      level <- pool_mode_level(read$surfaces, conf_level, call = call)
      pool_mode_frames(read, level, dfcom)
    },
    error = function(cnd) cnd
  )

  if (inherits(frames, "error")) {
    return(list(effects = other, reason = conditionMessage(frames)))
  }

  list(
    effects = other,
    estimates = frames$estimates,
    pooling = frames$pooling
  )
}

#' Construct a pooled inverse probability weighted result
#'
#' The field names and their order are the contract, since callers read fields
#' by name and a printed form writes them positionally. The constructor exists
#' so that order is stated in one place rather than at whatever point in
#' [pool_ipw()] the last field happens to be computed.
#'
#' The class does not inherit from `ipw`. A pooled result reports several
#' analyses rather than one fit, holds no component models, and refers its
#' inference to t rather than to z, so the methods registered against `ipw`
#' would answer for fields it does not carry.
#'
#' `alternate` is last rather than beside the fields it mirrors. The nine
#' components before it are what the return contract has always documented, and
#' a caller reading a pooled result positionally, or comparing one against a
#' result stored before the other reading was kept, finds them where they were.
#'
#' @param estimand,estimates,pooling,se_method,effects,m,dfcom,nobs,outcome_link
#'   The components, as the [pool_ipw()] return contract describes them.
#' @param alternate The reading that was not pooled actively, as
#'   `pool_alternate()` records it: a list naming that reading and holding
#'   either its two frames or the reason it has none.
#'
#' @return An S3 object of class `ipw_pooled`.
#'
#' @noRd
new_ipw_pooled <- function(
  estimand,
  estimates,
  pooling,
  se_method,
  effects,
  m,
  dfcom,
  nobs,
  outcome_link,
  alternate
) {
  structure(
    list(
      estimand = estimand,
      estimates = estimates,
      pooling = pooling,
      se_method = se_method,
      effects = effects,
      m = m,
      dfcom = dfcom,
      nobs = nobs,
      outcome_link = outcome_link,
      alternate = alternate
    ),
    class = "ipw_pooled"
  )
}

#' The results a `fits` argument holds
#'
#' A `mira` is what `mice::with()` returns, and it is a list of fits under an
#' `analyses` element and nothing else. Reading it by class rather than by
#' asking \pkg{mice} anything is what lets this package pool one without
#' depending on it.
#'
#' A single result is refused with its own message rather than through the
#' element check below. An `ipw` object is itself a list, so the check would
#' otherwise report that all seven of its fields are not results, which is true
#' and tells the caller nothing about what they did.
#'
#' @param fits The argument as the caller supplied it.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return The list of results.
#'
#' @noRd
pool_analyses <- function(fits, call = sys.call(-1)) {
  if (inherits(fits, "mira")) {
    fits <- fits$analyses
  }

  if (inherits(fits, "ipw")) {
    stop_invalid_argument(
      "fits",
      "be a list of `ipw` results rather than a single result",
      call = call
    )
  }

  if (!is.list(fits)) {
    stop_invalid_argument("fits", "be a list of `ipw` results", call = call)
  }

  wrong <- which(!vapply(fits, inherits, logical(1), what = "ipw"))
  if (length(wrong) == 1L) {
    stop_invalid_argument(
      "fits",
      paste0(
        "hold `ipw` results throughout, but the element at position ",
        wrong,
        " is not one"
      ),
      call = call
    )
  }
  if (length(wrong) > 1L) {
    stop_invalid_argument(
      "fits",
      paste0(
        "hold `ipw` results throughout, but the elements at positions ",
        format_series(as.character(wrong)),
        " are not"
      ),
      call = call
    )
  }

  if (length(fits) < 2L) {
    stop_invalid_argument(
      "fits",
      paste0(
        "hold at least two results, since the between-imputation variance is ",
        "estimated from the spread across them, but it holds ",
        length(fits)
      ),
      call = call
    )
  }

  fits
}

#' Refuse arguments the function has no name for
#'
#' Every argument after `fits` is matched by name, so a positional argument and
#' a misspelled one both land in the dots. Ignoring them would run the pooling
#' at the default the caller was trying to change and return a result that looks
#' like the one that was asked for.
#'
#' @param dots The dots, as a list.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return `NULL`, invisibly, when there are none.
#'
#' @noRd
check_pool_dots <- function(dots, call = sys.call(-1)) {
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }

  names <- names(dots)
  written <- vapply(
    seq_along(dots),
    function(i) {
      if (is.null(names) || is.na(names[[i]]) || !nzchar(names[[i]])) {
        "an unnamed argument"
      } else {
        paste0("`", names[[i]], "`")
      }
    },
    character(1)
  )

  stop_invalid_argument(
    "...",
    paste0(
      "be empty, since every argument after `fits` is matched by name, but ",
      format_series(written),
      " was passed"
    ),
    call = call
  )
}

#' Refuse a complete-data degrees of freedom that is not a count
#'
#' `Inf` is accepted, since it is the large-sample assumption written out and
#' the value the fallback itself uses. A missing count is refused along with the
#' rest: it would travel through the adjustment into an `NA` degrees of freedom
#' and an `NA` p-value for every effect, with nothing in the result to say where
#' it came from.
#'
#' @param dfcom The argument as the caller supplied it.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return `dfcom`, invisibly, when it is a count.
#'
#' @noRd
check_dfcom <- function(dfcom, call = sys.call(-1)) {
  is_count <- is.numeric(dfcom) && length(dfcom) == 1L && !is.na(dfcom)

  if (!is_count) {
    stop_invalid_argument(
      "dfcom",
      "be a single number, or `Inf` to assume a large sample",
      call = call
    )
  }

  invisible(dfcom)
}

#' The one value a field takes across the results, or a refusal
#'
#' The pooled estimate of an effect is an average of the per-imputation ones, so
#' every field the average depends on has to be the same in all of them. The
#' values are gathered into a list rather than a vector because one of the
#' fields checked this way is the vector of effect labels a result reports.
#'
#' @param x The results, or the surfaces read off them.
#' @param get A function reading the field from one of those.
#' @param field The field's name, which keys the error's class.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return The common value, which is `NULL` when every result records `NULL`.
#'
#' @noRd
pool_common <- function(x, get, field, call = sys.call(-1)) {
  values <- unique(lapply(x, get))
  if (length(values) > 1L) {
    stop_pool_mismatch(field, values, call = call)
  }
  values[[1L]]
}

#' The surface of one result that is being pooled
#'
#' Both readings come back in the same shape, so the pooling itself is written
#' once. The marginal reading is the effect estimates the result stores, keyed
#' by the effect labels the [new_ipw()] contract defines, which is what makes a
#' frame stored under the older `comparison` spelling poolable with one stored
#' under `contrast`: the labels are read through the same helper the accessors
#' read them with.
#'
#' The conditional reading is the outcome model's coefficients and the corrected
#' covariance a fitting package attached with [new_ipw_model()]. That covariance
#' is asked for through the helper the accessors use, so a model that carries
#' none is refused here exactly as it is when a caller asks the result for it.
#' The covariance the model computed for itself is not a substitute, and its
#' diagonal is what the standard errors would otherwise be taken from.
#'
#' The key is the columns that name a row, read through
#' `ipw_identity_columns()` and named by the role each column fills rather than
#' by the name the stored frame keeps it under. That is what carries a result's
#' contrasts and its subgroups onto the pooled frames under the canonical
#' headings, whichever vintage of the contract each result was written against.
#'
#' The set of readings the result declares is read first, before either surface
#' is. A result that does not support the reading being pooled has no such
#' surface whatever its frames and its outcome model happen to hold, so the
#' refusal a caller gets says the analysis has no such reading rather than
#' reporting whichever fault reading the surface would have run into. Asking
#' here is what asks once per reading rather than once per surface: this is the
#' only place a reading meets a result, so both `pool_mode_surfaces()` and the
#' guarded pass `pool_alternate()` makes over the other reading are answered by
#' the one check.
#'
#' @param fit One result.
#' @param effects The reading being pooled.
#' @param position The result's position in the set, which the refusal names so
#'   that a caller holding several of them knows which one to go back to.
#' @param call The call to report the error against, which is [pool_ipw()]'s
#'   rather than this helper's.
#'
#' @return A list of the row labels, the columns that key them, the estimates,
#'   their standard errors, their covariance or `NULL`, and the level the result
#'   stored or `NULL`.
#'
#' @noRd
#' @importFrom stats coef
pool_surface <- function(fit, effects, position = NULL, call = sys.call(-1)) {
  check_ipw_reading(fit, effects, position = position, call = call)

  if (effects == "conditional") {
    covariance <- conditional_vcov(fit$outcome_mod, call = call)
    estimate <- stats::coef(fit$outcome_mod)
    return(list(
      labels = names(estimate),
      key = list(effect = names(estimate)),
      estimate = unname(estimate),
      std.err = unname(sqrt(diag(covariance))),
      vcov = covariance,
      conf.level = NULL
    ))
  }

  estimates <- fit$estimates
  columns <- ipw_identity_columns(estimates)
  list(
    labels = ipw_effect_labels(estimates),
    key = lapply(columns, function(column) as.character(estimates[[column]])),
    estimate = estimates$estimate,
    std.err = estimates$std.err,
    vcov = attr(estimates, "ipw_vcov", exact = TRUE),
    conf.level = pool_stored_level(estimates)
  )
}

#' The confidence level an estimates frame records
#'
#' The stored level is a property of the frame rather than of a row, which is
#' the same reading `interval_bounds()` and `pooled_interval_bounds()` both make
#' of it. A frame whose rows disagree about the level names no one level for the
#' bounds to be reported at, and neither does a frame with no such column, so
#' both give `NULL` and are answered by the default or by a level named at the
#' call site.
#'
#' Both sides of the pooling read the frame this way, so they read it through
#' one helper. [pool_ipw()] asks it what level the results it is pooling agree
#' on, and `as.data.frame()` on the pooled result asks it what level the frame it
#' is reporting was built at. Two copies would be two places for that reading to
#' drift.
#'
#' @param estimates The `estimates` component of an `ipw` or an `ipw_pooled`
#'   object.
#'
#' @return A single number, or `NULL`.
#'
#' @noRd
pool_stored_level <- function(estimates) {
  level <- unique(estimates$conf.level)
  if (length(level) == 1L) level else NULL
}

#' The link an outcome model was fitted with
#'
#' The link says what scale the effects are on, so it is part of what the
#' results have to agree about. A model that is not a generalized linear one
#' reports no link of its own and its coefficients are on the identity scale,
#' which is the reading the ecosystem's fitting packages already make of one.
#'
#' @param model A result's `outcome_mod`.
#'
#' @return A single string.
#'
#' @noRd
ipw_outcome_link <- function(model) {
  link <- if (inherits(model, "glm")) model$family$link else NULL
  if (is.character(link) && length(link) == 1L) link else "identity"
}

#' The complete-data degrees of freedom the results report
#'
#' Three places in order: the fitted variance object each result records, then
#' each result's outcome model. The smallest is taken rather than the first,
#' since the pooled inference is no stronger than the weakest analysis supports.
#'
#' Results that report nothing are passed over rather than making the minimum
#' missing. A set can hold both kinds at once, since the linearization path
#' records no variance object and the M-estimation path does, and a minimum
#' taken across the missing ones as well would be `NA`: that travels into an
#' `NaN` degrees of freedom and an `NaN` p-value for every effect, with nothing
#' in the result to say where it came from.
#'
#' @param fits The results.
#'
#' @return A single number, or `NA_real_` when nothing reports one.
#'
#' @noRd
#' @importFrom stats df.residual
pool_dfcom <- function(fits) {
  from_fit <- vapply(
    fits,
    function(fit) pool_count(stats::df.residual(fit)),
    numeric(1)
  )
  if (!all(is.na(from_fit))) {
    return(min(from_fit, na.rm = TRUE))
  }

  from_outcome <- vapply(
    fits,
    function(fit) {
      # A result's outcome model is whatever the fitting package was given, so
      # the generic may have no method for it and `df.residual.default()` may
      # find no element to report. Neither is an error here: the count is looked
      # for in one more place after this one.
      pool_count(tryCatch(
        stats::df.residual(fit$outcome_mod),
        error = function(e) NULL
      ))
    },
    numeric(1)
  )
  if (!all(is.na(from_outcome))) {
    return(min(from_outcome, na.rm = TRUE))
  }

  NA_real_
}

#' One reported count, or `NA` when it is not one
#'
#' `df.residual()` answers with `NULL` for an object it has nothing to say
#' about, and a registered method is free to answer with an empty vector or
#' `Inf`. Each of those has to become a missing count rather than a value the
#' minimum could take.
#'
#' @param df What `df.residual()` reported.
#'
#' @return A single number, or `NA_real_`.
#'
#' @noRd
pool_count <- function(df) {
  if (length(df) != 1L || !is.numeric(df) || !is.finite(df)) {
    return(NA_real_)
  }
  as.numeric(df)
}

#' Rubin's rules for a matrix of per-imputation estimates
#'
#' One column per effect and one row per imputation. The pooled estimate is the
#' mean of the estimates, and the total variance is the within-imputation
#' variance plus the between-imputation variance inflated by `1 + 1/m` (Rubin,
#' 1987).
#'
#' @param estimate The per-imputation estimates, as an m by k matrix.
#' @param std.err Their standard errors, in the same shape.
#' @param dfcom The complete-data degrees of freedom.
#'
#' @return A list of numeric vectors, each with one element per effect.
#'
#' @noRd
#' @importFrom stats var
rubin_rules <- function(estimate, std.err, dfcom) {
  m <- nrow(estimate)
  qbar <- colMeans(estimate)
  ubar <- colMeans(std.err^2)
  b <- apply(estimate, 2L, stats::var)
  total <- ubar + (1 + 1 / m) * b

  # One effect at a time, which is the shape the adjustment is written for and
  # the shape `mice::pool()` evaluates it in. Its arithmetic happens to vectorize
  # over the effects, but the contract it is documented under is the scalar one,
  # and reading the answer for an effect out of a call that was given only that
  # effect's variances is what keeps a later change to it from quietly reporting
  # one effect's degrees of freedom for all of them.
  df <- vapply(
    seq_along(b),
    function(j) barnard_rubin(m, b[[j]], total[[j]], dfcom),
    numeric(1)
  )

  riv <- (1 + 1 / m) * b / ubar

  list(
    estimate = qbar,
    std.err = sqrt(total),
    df = df,
    ubar = ubar,
    b = b,
    riv = riv,
    lambda = (1 + 1 / m) * b / total,
    fmi = (riv + 2 / (df + 3)) / (riv + 1)
  )
}

#' The Barnard-Rubin small-sample degrees of freedom
#'
#' The adjustment of Barnard and Rubin (1999), which combines the
#' between-imputation degrees of freedom with what the complete-data analysis
#' had left. With no complete-data count to spend there is nothing to adjust and
#' the answer is the unadjusted `(m - 1) / lambda^2`.
#'
#' This is one effect's worth of arithmetic. `b` and `t` are single numbers and
#' so is the answer, which is how `mice::pool()` uses the same formula: it
#' evaluates it per group of a grouped summary. Callers with a vector of effects
#' apply it to each in turn, and `rubin_rules()` says why.
#'
#' Transcribed from that formula with one line dropped. mice ends on
#' `ifelse(is.infinite(dfcom), dfold, df_br)`, which guards a vector `dfcom`
#' holding a mix of finite and infinite counts; nothing here builds one, the
#' early return above already answers a `dfcom` that is infinite, and the line
#' is worse than redundant kept: `ifelse()` returns a value the length of its
#' test, so a scalar `dfcom` reaching it would collapse a vector `b` to its
#' first element.
#'
#' @param m The number of imputations.
#' @param b The between-imputation variance of one effect.
#' @param t The total variance of that effect.
#' @param dfcom The complete-data degrees of freedom.
#'
#' @return A single number.
#'
#' @noRd
barnard_rubin <- function(m, b, t, dfcom = Inf) {
  lambda <- (1 + 1 / m) * b / t
  dfold <- (m - 1) / lambda^2
  if (length(dfcom) == 1L && is.infinite(dfcom)) {
    return(dfold)
  }
  tmp <- (1 - lambda) * (1 + dfcom) * dfcom
  (m - 1) * tmp / ((dfcom + 3) * (m - 1) + lambda^2 * tmp)
}

#' The pooled effect estimates, as the returned frame reports them
#'
#' The statistic and the p-value are referred to t on the pooled degrees of
#' freedom rather than to the normal, and so are the bounds: with few
#' imputations those degrees of freedom can be in single figures, where the two
#' distributions are far apart.
#'
#' The bounds are one half width added and subtracted rather than a bound from
#' each tail, for the reason `confint()` builds its limits that way: `qt()` is
#' not exactly antisymmetric, so taking the lower bound from the lower tail
#' would put the two a little apart from each other.
#'
#' @param key The columns that name a row.
#' @param pooled The output of `rubin_rules()`.
#' @param conf.level The level the bounds report.
#'
#' @return A data frame.
#'
#' @noRd
#' @importFrom stats pt qt
pool_estimates_frame <- function(key, pooled, conf.level) {
  statistic <- pooled$estimate / pooled$std.err
  half_width <- stats::qt(1 - (1 - conf.level) / 2, pooled$df) * pooled$std.err

  data.frame(
    c(
      key,
      list(
        estimate = pooled$estimate,
        std.err = pooled$std.err,
        t = statistic,
        df = pooled$df,
        ci.lower = pooled$estimate - half_width,
        ci.upper = pooled$estimate + half_width,
        conf.level = conf.level,
        p.value = 2 * stats::pt(-abs(statistic), pooled$df)
      )
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' The pooling diagnostics, as the returned frame reports them
#'
#' What the pooling itself produced, which says how much of the uncertainty came
#' from the imputation rather than from the data. They sit in a frame of their
#' own rather than beside the estimates, keyed the same way so that a row of one
#' is found from a row of the other.
#'
#' @param key The columns that name a row.
#' @param pooled The output of `rubin_rules()`.
#'
#' @return A data frame.
#'
#' @noRd
pool_diagnostics_frame <- function(key, pooled) {
  data.frame(
    c(
      key,
      list(
        ubar = pooled$ubar,
        b = pooled$b,
        riv = pooled$riv,
        lambda = pooled$lambda,
        fmi = pooled$fmi
      )
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' The pooled covariance of the effects
#'
#' The rule the standard errors follow, applied to the whole matrix: the average
#' of the within-imputation covariances plus the between-imputation one inflated
#' by `1 + 1/m`. Its diagonal is therefore the square of the pooled standard
#' errors, since the two are the same combination of the same numbers.
#'
#' A result that carries no covariance leaves nothing to average its entry with.
#' The alternative is a matrix built from the results that have one, which would
#' sit beside estimates built from all of them with nothing on the object to say
#' so, and the standard errors the frames hold give only the diagonal.
#'
#' @param surfaces The surfaces read off the results.
#' @param estimate The per-imputation estimates, as an m by k matrix.
#' @param labels The effect labels, which name both margins.
#'
#' @return A covariance matrix, or `NULL` when any result carries none.
#'
#' @noRd
#' @importFrom stats cov
pool_vcov <- function(surfaces, estimate, labels) {
  blocks <- lapply(surfaces, function(x) x$vcov)
  if (any(vapply(blocks, is.null, logical(1)))) {
    return(NULL)
  }

  m <- length(blocks)
  covariance <- Reduce(`+`, blocks) / m + (1 + 1 / m) * stats::cov(estimate)
  dimnames(covariance) <- list(labels, labels)
  covariance
}
