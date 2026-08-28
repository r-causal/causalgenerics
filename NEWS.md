# causalgenerics (development version)

* `new_ipw()` now checks the `estimates` component it is handed. It has to be a
  data frame, and it has to name each row's effect measure in an `effect` column
  of strings or of factor levels, with a measure named in every row; anything
  else is refused with an error of class
  `causalgenerics_invalid_argument_estimates`, and of the general class
  `causalgenerics_invalid_argument`. Every surface of a result reports from that
  frame and each of them begins by reading that column, so a component of
  another shape reached them as a subscript error about the column that failed
  rather than about the argument that was wrong, and which surface a caller met
  first depended on what they asked the result for. A missing measure was
  accepted where a missing subgroup was not: it pasted into a label reading
  `"NA"`, and two of them pasted into the same label, which left `coef()` naming
  two elements the same thing and `vcov()` repeating a name down both of its
  dimnames. A factor is accepted beside a character column, since the labels
  come off it as strings either way and a frame written when `data.frame()`
  coded its strings as factors is a frame a caller still has in hand.

* `new_ipw()` gains a `readings` field, the eighth and last of the result
  contract, recording which of the two readings a result supports. Both
  surfaces exist on most results and not on all of them: an exposure that
  enters the outcome model through several columns has no single coefficient to
  read as the conditional effect, so a package that builds such a result records
  the reading it can answer for, and the mode the result records has to be one
  of them. Asking a result for a reading it does not support is refused with an
  error of class `causalgenerics_unsupported_reading_marginal` or
  `causalgenerics_unsupported_reading_conditional`, and of the general class
  `causalgenerics_unsupported_reading`, which carries the reading asked for and
  the set the result records as fields. The refusal is raised wherever a reading
  is named: by `as_marginal()` and `as_conditional()`, which are no longer total
  on the class, and by the `effects` argument of `coef()`, `vcov()`,
  `confint()`, and `as.data.frame()`. `as.data.frame()` itself gains that
  `effects` argument and now presents the reading the result records rather than
  the effects table in either reading: the conditional reading reports the
  outcome model's coefficients, with the standard errors the corrected
  covariance implies, under the headings the pooled table uses for the same
  reading. A result stored before the field existed records no set and is read
  as supporting both, so an older result and a result built today behave the
  same way. `pool_ipw()` honors the field too: it reads the readings the
  results declare before either surface, so a set whose results support one
  reading pools that one alone. Naming the other refuses the call with the same
  error a result raises directly, and leaving `effects` at `NULL` pools the
  reading the results record and puts that refusal on `alternate` as the reason
  the other has no frames, so `as_marginal()` or `as_conditional()` on the
  pooled result raises `causalgenerics_pool_missing_surface_marginal` or
  `causalgenerics_pool_missing_surface_conditional` carrying it.

* New exposure-type machinery gives the ecosystem one reading of whether an
  exposure is binary, categorical, or continuous. `detect_exposure_type()`
  classifies a vector, `match_exposure_type()` resolves the type a caller
  declared against the types a function supports, `check_forced_type()` refuses
  a declaration the data cannot carry, and `observed_values()`,
  `has_two_levels()`, and `is_categorical()` are the readings those are built
  from. The classification was duplicated as an internal in four packages, and
  the copies had drifted: two of them counted a missing value as a level of its
  own, which reads a binary exposure recorded with any missingness as a
  three-level one and takes it off the binary path entirely. The reading here
  counts only the values an exposure is observed at, so a single-level factor is
  a degenerate binary exposure, a numeric vector is categorical when it takes
  fewer distinct values than 20 percent of its non-missing observations, and
  everything else is continuous. Detection announces what it read in a line
  naming the exposure argument and the type, which a caller silences for one
  call with `announce = FALSE` and a user silences for every call with
  `options(causalgenerics.quiet = TRUE)`. A declared type wins over the
  heuristics and is returned without the data being consulted at all, since a
  ten-level dose declared continuous is fitted as a dose; a type this package
  reads but the calling function does not support is refused with
  `causalgenerics_unsupported_exposure_type`, whether the caller named it or
  detection reached it, and a declaration no data of that shape could support is
  refused with `causalgenerics_forced_exposure_type`. Both carry the umbrella
  class `causalgenerics_error` and report the environment passed as `call`, so a
  consumer package names its own function in the refusal. These are the first
  conditions in the package built with `cli` rather than base
  `errorCondition()`, which is what keeps their wording identical to the
  refusals the consuming packages already show; the conditions raised elsewhere
  in the package are unchanged. `cli` and `rlang` are new imports.

* New `joint_exposure()` crosses two discrete treatments into one categorical
  exposure and records the declaration on the vector it returns: which two
  treatments were crossed, in what order, which levels each of them takes, and
  which cell the effects are reported against. `is_joint_exposure()`,
  `joint_components()`, and `joint_reference()` read that back, so the package
  that builds the weights, the package that checks balance over the cells, and
  the package that reports the effects agree on which cells exist without any of
  them owning the declaration. Conditioning on a variable and jointly
  intervening on two treatments are different causal questions, and the class
  serves the second. The result is a factor, with `"factor"` ahead of
  `"vctrs_vctr"` in its class vector so that `model.matrix()` gives exactly what
  a plain factor over the same cells gives, and the cells are ordered with the
  first component varying fastest, which is what puts the cell crossing the two
  components' own reference levels first. Construction refuses a crossing with a
  cell nothing falls in, with an error of class
  `causalgenerics_joint_exposure_empty_cell`: an empty cell is a positivity
  violation rather than a small sample, so the joint effect is not identified in
  those data and no weight, contrast, or balance check over the crossing can be
  computed. A component with a missing value is refused with
  `causalgenerics_joint_exposure_missing_value`, and so is a factor component
  that declares `NA` as one of its levels even when nothing currently takes it,
  since such a level gives its missing observations an ordinary code and the
  cell the crossing would write for them names an absence rather than a
  treatment. A component observed at one level, an unnamed component, and
  anything other than exactly two components are refused as well, each under the
  general class `causalgenerics_invalid_joint_exposure` alongside a class keyed
  to the fault. So are two components sharing a name, with
  `causalgenerics_joint_exposure_shared_name`: the names are what tell the two
  treatments apart in the cells, so one name used twice would label every cell
  with the same variable on either side of it and leave nothing downstream able
  to tell which component a cell varies. The declaration then survives every
  operation that keeps the full set of cells, including `[`,
  `vctrs::vec_slice()`, `sort()`, a round trip through a data frame, and
  combining two joint exposures that declare the same crossing. Operations that
  narrow or rewrite the level set give it up and say so, each warning under the
  general class `causalgenerics_joint_exposure_downgrade`: `droplevels()` and
  `x[i, drop = TRUE]` return a plain factor over the cells that remain,
  `levels<-` returns a plain factor with the new levels, and combining with a
  factor or with a joint exposure declaring a different crossing returns a plain
  factor over the levels of both; combining with a character vector returns a
  character vector. Casting into the class is refused outright rather than
  degraded, because a crossing cannot be recovered from labels that merely look
  like one. Base `c()` is routed through `vctrs::vec_c()` so that it takes those
  same rules, rather than reaching `c.factor()` and combining the codes into a
  bare factor without a word. Two routes escape all of it and degrade silently,
  because neither reaches a method this package can register: `unlist()`
  combines the underlying codes in base C code without dispatching at all, and
  `c()` dispatches on its first argument, so a combine that begins with a plain
  factor never reaches this class. Reach for `vctrs::vec_c()` where either could
  apply.

* The `estimates` frame of an `ipw` result may carry a `group` column, naming
  the subgroup each row was estimated in as a `"var = value"` string such as
  `"sex = 0"`. A row's identity is the `effect` column, then `contrast` when the
  result names contrasts, then `group` when it names subgroups, and every
  surface that reads that identity reads the new column with it. An effect label
  pastes the group on last, so a subgroup row is labeled `"rd b vs a sex = 0"`;
  `print()` writes those labels down the side of its table and keeps the column
  out of the numeric matrix it formats; `coef()` names its vector and `vcov()`
  its dimnames with them; `as.data.frame()` heads a `group` column after
  `contrast` and before `estimate`, on a result and on a pooled one alike; and
  `pool_ipw()` keys its pooled frames by the same columns, so a set of results
  reporting different subgroups, or the same subgroups in different orders,
  disagrees about its labels and is refused through that requirement rather than
  through one of its own. The column has one spelling and no alias, unlike
  `contrast`: it is newer than the first version of the contract, so no stored
  frame keeps its subgroups anywhere else. It is absent rather than constant
  when a result reports one group, since a column repeating one value down the
  table would read as a subgroup that was named. `new_ipw()` refuses a `group`
  column that is not character, and one that leaves any row without a subgroup,
  with an error of class `causalgenerics_invalid_argument_estimates`, and of the
  general class `causalgenerics_invalid_argument`. A label is the identity
  columns pasted together, so anything at all pastes into something and nothing
  downstream errors: a factor pastes as its levels rather than as what the frame
  holds, a number names a value without saying which variable took it, and a
  missing entry pastes into a label reading `"rd NA"`; two of them paste into
  the same label, naming two rows the same thing. What comes out either way is
  a label the package that produced the result did not write, which is why the
  column is checked where the result is constructed.

* `pool_ipw()` pools both readings of the results it is given whenever both can
  be pooled, and stores the one the call did not name whole, as a tenth
  `alternate` component holding that reading's pooled estimates and pooling
  diagnostics. Which reading `effects` names therefore decides which one the
  result presents rather than what was computed at all, and pooling one reading
  and moving to the other gives the frames pooling the other directly gives. The
  reading the call names is pooled first and is not guarded, so
  `pool_ipw(fits, effects = "conditional")` over results whose outcome models
  carry no corrected covariance is refused exactly as it was before. The other
  reading is pooled under a guard, since a set that cannot be pooled on it is a
  result with one reading rather than a failed pooling: `alternate` then names
  that reading and records the refusal it raised, in the words the refusal used.

* `as_marginal()` and `as_conditional()` move a pooled result between the two
  readings it carries, exchanging the pooled estimates, the pooling diagnostics,
  and the recorded mode with the reading stored beside them and leaving the
  components that describe the analyses where they are. `coef()`, `vcov()`,
  `confint()`, and `as.data.frame()` on a pooled result each gain an `effects`
  argument naming a reading for one call without changing the result, and
  `estimand()` now answers on a pooled result with the estimand its analyses
  agreed on. Asking a pooled result for a reading the pooling could not compute
  raises an error of class `causalgenerics_pool_missing_surface_marginal` or
  `causalgenerics_pool_missing_surface_conditional`, and of the general class
  `causalgenerics_pool_missing_surface`, carrying the recorded refusal as its
  reason; a result pooled before both readings were kept refuses the same way
  and says that pooling the results again gives one that carries both.

* A pooled result carries the methods a fitted model answers to: `print()`,
  `coef()`, `vcov()`, `confint()`, `nobs()`, and `as.data.frame()`. Every one of
  them refers the inference to t on each effect's own pooled degrees of freedom
  rather than to the normal, which is what the Barnard-Rubin adjustment leaves
  and can be in single figures when there are few imputations. `print()` reports
  the estimand, the reading, how many imputations went in, and the complete-data
  degrees of freedom, then tabulates the pooled effects; the pooling diagnostics
  stay in the `pooling` component rather than crowding that table.
  `as.data.frame()` reports the tidier-shaped table, with a `df` column after
  the statistic so the statistic beside it can be referred to something, and its
  `conf.level` defaults to the level the result records rather than to 0.95.
  `exponentiate = TRUE` moves the `log(rr)` and `log(or)` rows of a marginal
  table to their natural scale as it does for an unpooled result; on a
  conditional table it moves every coefficient when the outcome models were
  fitted with a `logit` or `log` link, relabeling nothing, and raises an error
  of class `causalgenerics_exponentiate_link` on any other link. There is
  deliberately no `df.residual()` method: residual degrees of freedom belong to
  one fit, and the per-effect pooled count is in the table.

* New `pool_ipw()` combines the `ipw` results fitted to each of a set of
  multiply imputed datasets into one result, by Rubin's rules, and returns it
  under the new class `ipw_pooled`. It takes a plain list of results or the
  `mira` object `mice::with()` returns, which it recognizes by class and reads
  through its `analyses` element, so mice is not a dependency. Each result
  already carries a standard error that accounts for the weights having been
  estimated, and pooling those adds the uncertainty the imputation contributed.
  The pooled degrees of freedom carry the Barnard-Rubin small-sample
  adjustment, so the statistic, the p-value, and the bounds are referred to t
  rather than to the normal; the complete-data degrees of freedom are read off
  the results when `dfcom` does not name them, and a large sample is assumed
  with a warning when nothing reports them. Results that disagree about their
  estimand, standard error method, presentation mode, reported effects, stored
  confidence level, or outcome model link are refused with an error of class
  `causalgenerics_pool_mismatch` rather than averaged. The ratio effects are
  pooled on the log scale they were estimated on.

* Breaking change: the column a categorical `ipw` result names its contrasts
  with is `contrast`, in the `estimates` frame `new_ipw()` takes and in the
  table `as.data.frame()` reports. It was `comparison`. Both things that read
  the table by column name spell it the new way: `mice::pool()` groups the rows
  it pools by a closed whitelist of names that holds `contrast` and not
  `comparison`, and marginaleffects heads the same column `contrast`, so a
  table under the old name lost the contrast it was keyed by and did not stack
  with one of theirs. A stored result whose frame carries the old name is still
  read: the older spelling is an alias wherever the column is read, so such a
  result labels its rows and reports its table exactly as one built now does. A
  package supplying an `ipw()` method writes `contrast`.

* Breaking change: `as.data.frame()` on an `ipw` result returns the
  tidier-shaped table of its effect estimates rather than a copy of the
  `estimates` component. The columns are `term`, `contrast` when the result
  names contrasts, `estimate`, `std.error`, `statistic`, and `p.value`, so the
  storage names `effect`, `std.err`, and `z` no longer appear and code reading
  the table by those names has to be updated. The interval is opt-in through the
  new `conf.int` argument, which appends `conf.low` and `conf.high` after the
  other columns, and the level it reports is the new `conf.level` argument
  rather than a `conf.level` column. These are the column names the tidier
  convention uses, so a fitting package's `tidy()` method for an `ipw` result is
  this table read as a tibble and nothing more, and what such a method returns
  is settled here rather than in every package that has one.

* `as.data.frame()` returns the bounds `estimates` stores only when every row of
  the frame records the level asked for, and otherwise builds the normal
  approximation from the estimate and its standard error. The stored level is a
  property of the frame rather than of a row: a table mixing stored bounds with
  recomputed ones would report two intervals under one pair of column headings.
  A frame that records no level has no level for its stored bounds to be
  returned at, so those are recomputed too.

* `as.data.frame()` refuses a `conf.int` or `exponentiate` that is not a single
  `TRUE` or `FALSE`, and a `conf.level` that is not a single number strictly
  between 0 and 1, with an error of class
  `causalgenerics_invalid_argument`. A number or a string takes the true branch
  of a bare `if`, so an unchecked `conf.int = "no"` would report an interval.
  `conf.level` is checked whether or not the bounds were asked for, since a call
  naming a level no interval can be built at has said something wrong either
  way.

* `as.data.frame(exponentiate = TRUE)` drops the `ipw_vcov` attribute rather
  than carrying it. The covariance describes the estimates on the scale they
  were estimated on, and once the ratio rows have moved to their natural scale
  it describes neither the table it is attached to nor anything else. There is
  no correct matrix to put in its place, since the delta method answer is not
  what the fitting package computed. On the log scale the attribute travels on
  the returned table, including when the bounds were recomputed.

* `df.residual()` on an `ipw` result reports a fractional degrees of freedom as
  the double the fit reports, rather than truncating it to an integer. A
  penalized or smooth fit spends a fractional count of parameters, and
  truncating it would report a fit that spent more of them than it did. A whole
  number still comes back as an integer.

* The conditional reading pairs the corrected covariance with the outcome
  model's coefficients by name. A block a fitting package attached in the order
  of its own stacked system is reported in coefficient order, so the variance
  read beside a coefficient is that coefficient's, and `print()` writes each
  standard error against the coefficient it belongs to. A block whose labels
  cannot be paired with the coefficients is refused with an error of class
  `causalgenerics_conditional_vcov_mismatch`. A block of another size, one
  labeled with the parameter names of a stacked system, and a model whose
  coefficients carry no names are the three ways the pairing fails, and reading
  such a block by position would report the covariance of other parameters under
  this model's coefficient names.

* `print()` on a conditional result whose outcome model reports no coefficients
  says so in place of the table, rather than writing an empty one. The
  covariance is not asked for on that path, since a table with no rows has no
  use for one.

* `vcov()` on a model wrapped by `new_ipw_model()` raises an error of class
  `causalgenerics_no_vcov_ipw_model` when the wrapper no longer carries a
  covariance, rather than returning `NULL`. The class and the attribute go on
  together, so a model carrying one without the other passed through code that
  dropped its attributes and kept its class. `NULL` traveled from there:
  standard errors taken from it are an empty vector, and the limits built from
  those come back as `NA` with nothing said about why. `print()` on a
  conditional result whose wrapper lost its covariance refuses for the same
  reason, rather than printing the coefficients under a note telling the reader
  to wrap a model that is already wrapped.

* The error a missing covariance raises now says that the package which produced
  the object attaches one when it supports the `ipw_vcov` contract, in place of
  telling the reader to refit with a current version. The same helper answers
  for a result and for a component model of one, and a package that never
  attached a covariance is not one an upgrade would fix.

* `model.frame()` on an `ipw` result returns the outcome model's model frame, so
  prediction and averaging tooling can recover the data an estimate was computed
  from out of the result itself. The `(weights)` column a weighted frame carries
  is dropped, since a package that reads the columns of a frame as the variables
  a model was fitted on would otherwise treat the estimation weights as one of
  them. `weights()` reports those.

* `estimand()` on an `ipw` result returns the estimand the weights it was
  computed under targeted. There is no replacement method for a result:
  assigning a new estimand would relabel the estimates rather than recompute
  them.

* `print()` on an `ipw` result names the reading it is showing, on an `Effects:`
  line beside the estimand and in the heading of the table. The marginal
  reading tabulates the effect estimates as before, under `Marginal estimates:`.
  The conditional reading tabulates the outcome model's coefficients with the
  standard errors the corrected covariance implies, and prints the coefficients
  on their own, with a note, when the outcome model records no such covariance.

* `coef()`, `vcov()`, and `confint()` on an `ipw` result report the reading the
  `effects` field records, and each gains an `effects` argument that names a
  reading for one call without changing the result. The conditional reading
  reports the outcome model's coefficient surface: its coefficients, the
  covariance the joint estimation of the weights and the outcome implies, and
  the normal limits built from that covariance. `nobs()`, `df.residual()`, and
  `weights()` describe the fit rather than a surface of it and answer the same
  way in either reading.

* `vcov()` and `confint()` refuse the conditional reading when the outcome model
  carries no corrected covariance, with an error of class
  `causalgenerics_no_conditional_vcov`. The covariance a weighted model computes
  for itself treats the estimated weights as fixed, so it is not a substitute; a
  fitting package supplies the corrected block by wrapping the outcome model
  with `new_ipw_model()`.

* `new_ipw()` gains an `effects` argument, stored as a seventh field, recording
  which reading a result presents: `"marginal"` for the causal contrast
  estimates or `"conditional"` for the outcome model's coefficient surface. A
  method that names no mode reports marginal effects, as every method did
  before.

* New generics `as_marginal()` and `as_conditional()` move a result between the
  two readings. Both surfaces exist on every result, so the methods on `ipw`
  set the mode and leave every other field alone.

* `ipw` results gain the accessors a fitted model answers to: `coef()`,
  `vcov()`, `confint()`, `nobs()`, `df.residual()`, and `weights()`. A method
  that can compute the covariance of the effects it reports now attaches it to
  its `estimates` data frame as the `ipw_vcov` attribute, which is what `vcov()`
  reads.

* New `new_ipw_model()` wraps a component model of an IPW fit with the
  covariance the two-step estimation implies, so that `vcov()` on the model
  accounts for the weights having been estimated rather than fixed.

# causalgenerics 0.1.0

* Initial CRAN submission.
