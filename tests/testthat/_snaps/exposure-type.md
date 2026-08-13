# detect_exposure_type() announces the type it read

    Code
      detect_exposure_type(cg_binary_numeric())
    Message
      i Treating `.exposure` as binary
    Output
      [1] "binary"

---

    Code
      detect_exposure_type(cg_three_level_factor())
    Message
      i Treating `.exposure` as categorical
    Output
      [1] "categorical"

---

    Code
      detect_exposure_type(cg_continuous_numeric())
    Message
      i Treating `.exposure` as continuous
    Output
      [1] "continuous"

# detect_exposure_type() announces the argument it was told to name

    Code
      detect_exposure_type(cg_binary_numeric(), arg = ".exposures")
    Message
      i Treating `.exposures` as binary
    Output
      [1] "binary"

# match_exposure_type() announces the detected type

    Code
      match_exposure_type("auto", cg_categorical_numeric())
    Message
      i Treating `.exposure` as categorical
    Output
      [1] "categorical"

# match_exposure_type() threads the argument name into detection

    Code
      match_exposure_type("auto", cg_binary_numeric(), arg = ".exposures")
    Message
      i Treating `.exposures` as binary
    Output
      [1] "binary"

# match_exposure_type() refuses a string that is not a type

    Code
      match_exposure_type("bimary", cg_binary_numeric())
    Condition
      Error:
      ! `exposure_type` must be one of "auto", "binary", "categorical", or "continuous", not "bimary".
      i Did you mean "binary"?

# match_exposure_type() refuses a known type the caller cannot fit

    Code
      match_exposure_type("continuous", cg_binary_numeric(), valid_types = c("auto",
        "binary", "categorical"))
    Condition
      Error:
      ! Exposure type "continuous" is not supported.
      i Supported exposure types: "binary" and "categorical".

# the unsupported refusal covers a detected type as well

    Code
      match_exposure_type("auto", cg_continuous_numeric(), valid_types = c("auto",
        "binary", "categorical"))
    Message
      i Treating `.exposure` as continuous
    Condition
      Error:
      ! Exposure type "continuous" is not supported.
      i Supported exposure types: "binary" and "categorical".

# the unsupported refusal does not offer auto as a type

    Code
      match_exposure_type("categorical", cg_binary_numeric(), valid_types = c("auto",
        "binary"))
    Condition
      Error:
      ! Exposure type "categorical" is not supported.
      i Supported exposure types: "binary".

# both refusals name the function the caller called

    Code
      cg_wrap_match("continuous", cg_binary_numeric(), c("auto", "binary"))
    Condition
      Error in `cg_wrap_match()`:
      ! Exposure type "continuous" is not supported.
      i Supported exposure types: "binary".

---

    Code
      cg_wrap_match("bimary", cg_binary_numeric(), c("auto", "binary"))
    Condition
      Error in `cg_wrap_match()`:
      ! `exposure_type` must be one of "auto" or "binary", not "bimary".
      i Did you mean "binary"?

# a caller narrowing valid_types narrows its own default too

    Code
      match_exposure_type(.exposure = cg_binary_numeric(), valid_types = c("auto",
        "binary"))
    Condition
      Error:
      ! `exposure_type` must be one of "auto" or "binary", not "auto".

# check_forced_type() refuses a binary claim over more than two levels

    Code
      check_forced_type("binary", cg_three_level_factor())
    Condition
      Error:
      ! `exposure_type` was set to "binary", but `.exposure` cannot be treated that way.
      x A "binary" exposure takes exactly two observed values, and `.exposure` takes 3.
      i Drop `exposure_type` to detect the type from the data, which reads `.exposure` as "categorical".

# check_forced_type() refuses a continuous claim over a factor

    Code
      check_forced_type("continuous", cg_two_level_character())
    Condition
      Error:
      ! `exposure_type` was set to "continuous", but `.exposure` cannot be treated that way.
      x A "continuous" exposure is numeric, and `.exposure` is a character vector.
      i Drop `exposure_type` to detect the type from the data, which reads `.exposure` as "binary".

# check_forced_type() reports what detection would have read

    Code
      check_forced_type("binary", rep(1, 10))
    Condition
      Error:
      ! `exposure_type` was set to "binary", but `.exposure` cannot be treated that way.
      x A "binary" exposure takes exactly two observed values, and `.exposure` takes 1.
      i Drop `exposure_type` to detect the type from the data, which reads `.exposure` as "categorical".

# check_forced_type() names the argument it was told to name

    Code
      check_forced_type("binary", cg_three_level_factor(), arg = ".exposures")
    Condition
      Error:
      ! `exposure_type` was set to "binary", but `.exposures` cannot be treated that way.
      x A "binary" exposure takes exactly two observed values, and `.exposures` takes 3.
      i Drop `exposure_type` to detect the type from the data, which reads `.exposures` as "categorical".

# check_forced_type() names the function the caller called

    Code
      cg_wrap_forced("binary", cg_three_level_factor())
    Condition
      Error in `cg_wrap_forced()`:
      ! `exposure_type` was set to "binary", but `.exposure` cannot be treated that way.
      x A "binary" exposure takes exactly two observed values, and `.exposure` takes 3.
      i Drop `exposure_type` to detect the type from the data, which reads `.exposure` as "categorical".

# check_forced_type() refuses a type it does not know

    Code
      check_forced_type("auto", cg_binary_numeric())
    Condition
      Error:
      ! `forced` must be one of "binary", "categorical", or "continuous", not "auto".

---

    Code
      check_forced_type("bimary", cg_binary_numeric())
    Condition
      Error:
      ! `forced` must be one of "binary", "categorical", or "continuous", not "bimary".
      i Did you mean "binary"?

