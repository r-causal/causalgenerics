# the two-component refusal states the contract

    Code
      joint_exposure(qsmk = smoking_qsmk())
    Condition
      Error in `joint_exposure()`:
      ! A joint exposure is the crossing of exactly two treatments, and 1 was supplied; pass one named vector per treatment.

# the unnamed-component refusal states the contract

    Code
      joint_exposure(smoking_qsmk(), exercise = smoking_exercise())
    Condition
      Error in `joint_exposure()`:
      ! Every component of a joint exposure must be named, because the names are the treatment names its cell labels are written from; pass each treatment as `name = value`.

# the shared-name refusal states the contract

    Code
      joint_exposure(qsmk = smoking_qsmk(), qsmk = smoking_exercise())
    Condition
      Error in `joint_exposure()`:
      ! Both components of a joint exposure are named `qsmk`, so its cells would name that one treatment twice and nothing reading the declaration could tell which component a cell varies; give the two treatments distinct names.

# the constant-component refusal names the component

    Code
      joint_exposure(qsmk = constant, exercise = smoking_exercise())
    Condition
      Error in `joint_exposure()`:
      ! `qsmk` takes one value in these data, so it is not a treatment the crossing can vary; a joint exposure needs both of its components observed at two or more levels.

# the missing-value refusal states the contract

    Code
      joint_exposure(qsmk = smoking_qsmk(), exercise = with_na)
    Condition
      Error in `joint_exposure()`:
      ! `exercise` admits a missing exposure, and an observation whose exposure is unknown falls in no cell of the crossing; drop or recover those observations, and any `NA` the component declares as a level, before declaring the joint exposure.

# the declared-NA-level refusal states the contract

    Code
      joint_exposure(qsmk = c(0, 1, 0, 1), exercise = declared_na)
    Condition
      Error in `joint_exposure()`:
      ! `exercise` admits a missing exposure, and an observation whose exposure is unknown falls in no cell of the crossing; drop or recover those observations, and any `NA` the component declares as a level, before declaring the joint exposure.

# the empty-cell refusal names the cell and the positivity failure

    Code
      joint_exposure(qsmk = never_both, exercise = exercise)
    Condition
      Error in `joint_exposure()`:
      ! Nothing in these data falls in the cell "qsmk = 1, exercise = yes", so the crossing violates positivity and the joint effect of the two treatments is not identified here; coarsen a component or restrict the analysis to the cells that are populated.

# the printed and formatted forms are recorded

    Code
      print(x)
    Output
      <joint_exposure[4]: qsmk, exercise>
      [1] qsmk = 0, exercise = no  qsmk = 1, exercise = no  qsmk = 0, exercise = yes
      [4] qsmk = 1, exercise = yes
      Reference: qsmk = 0, exercise = no

---

    Code
      format(x)
    Output
      [1] "qsmk = 0, exercise = no " "qsmk = 1, exercise = no "
      [3] "qsmk = 0, exercise = yes" "qsmk = 1, exercise = yes"

---

    Code
      vctrs::vec_ptype_full(x)
    Output
      [1] "joint_exposure<qsmk, exercise>"

---

    Code
      vctrs::vec_ptype_abbr(x)
    Output
      [1] "jnt_exp"

