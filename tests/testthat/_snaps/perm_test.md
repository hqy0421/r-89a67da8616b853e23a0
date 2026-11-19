# perm_test fails on rank deficiency

    Code
      perm_test(bad, full_f, red_f, tested, reps = 10, blocks = blk)
    Condition
      Error:
      ! The full model's design matrix is rank deficient (rank 3 with 4 columns), so the tested coefficients are not uniquely defined. Remove redundant terms before calling `perm_test()`.

# perm_test fails when the reduced model is not nested

    Code
      perm_test(d, y ~ x1 + z, y ~ x2, "x1", reps = 10)
    Condition
      Error:
      ! `reduced_formula` must be nested within `full_formula`: every term in the reduced model must also appear in the full model.

# perm_test fails when the formulas have different responses

    Code
      perm_test(d, y ~ x1, z ~ 1, "x1", reps = 10)
    Condition
      Error:
      ! `full_formula` and `reduced_formula` must have the same response variable, but got "y" and "z".

# perm_test fails on mismatched terms

    Code
      perm_test(d, full_f, red_f, "nope", reps = 10)
    Condition
      Error:
      ! The coefficient "nope" supplied to `terms` is not in the full model. The full model's coefficients are "(Intercept)", "z", "x1", and "x2".

---

    Code
      perm_test(d, full_f, red_f, c("x1", "z"), reps = 10)
    Condition
      Error:
      ! The coefficient "z" supplied to `terms` is also in the reduced model and so is not under test.

---

    Code
      perm_test(d, full_f, red_f, "x1", reps = 10)
    Condition
      Error:
      ! `terms` must name every coefficient that distinguishes the full model from the reduced model. Also include "x2".

# perm_test fails on invalid data

    Code
      perm_test(d_na, full_f, red_f, tested, reps = 10)
    Condition
      Error:
      ! `perm_test()` does not support missing values; 1 missing value was found in the variables used by the test. Complete or remove the relevant rows before calling the function.

---

    Code
      perm_test(d_inf, full_f, red_f, tested, reps = 10)
    Condition
      Error:
      ! The variable "z" contains non-finite values (`NaN` or `Inf`), which are not valid for a Gaussian linear model.

---

    Code
      perm_test(d_cat, full_f, red_f, tested, reps = 10)
    Condition
      Error:
      ! Gaussian linear models require a numeric response, but "y" is a <factor> object.

# perm_test fails on degenerate block structure

    Code
      perm_test(d_singleton, y ~ x1, y ~ 1, "x1", reps = 10, blocks = id)
    Condition
      Error:
      ! Every block contains a single observation, so no residuals can be permuted. Supply a blocking variable with at least one block containing multiple observations.

---

    Code
      perm_test(d, y ~ x1, y ~ 1, "x1", reps = 10, blocks = absent)
    Condition
      Error:
      ! The blocking variable "absent" is not a column in `x`.

# perm_test fails on invalid scalar arguments

    Code
      perm_test(d, full_f, red_f, tested, reps = 0)
    Condition
      Error:
      ! `reps` must be a single positive integer, not a number.

---

    Code
      perm_test(d, full_f, red_f, tested, reps = 3.5)
    Condition
      Error:
      ! `reps` must be a single positive integer, not a number.

---

    Code
      perm_test(d, full_f, red_f, tested, reps = 10, direction = "sideways")
    Condition
      Error:
      ! The provided value for `direction` is not appropriate. Possible values are "less", "greater", "two-sided", "left", "right", "both", "two_sided", "two sided", or "two.sided".

---

    Code
      perm_test(d, ~x1, red_f, tested, reps = 10)
    Condition
      Error:
      ! `full_formula` must be a two-sided formula, e.g. `y ~ x + z`.

---

    Code
      perm_test(list(a = 1), full_f, red_f, tested, reps = 10)
    Condition
      Error:
      ! `x` must be 'data.frame', not 'list'.

---

    Code
      perm_test(d, full_f, red_f, c("x1", "x1"), reps = 10)
    Condition
      Error:
      ! `terms` must not contain duplicate coefficient names.

