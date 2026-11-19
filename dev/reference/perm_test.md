# Permutation tests for regression coefficients

A tidier implementation of the Freedman–Lane residual permutation
procedure for testing one or more coefficients of a Gaussian linear
model while adjusting for covariates and, optionally, restricting
permutations to within exchangeability blocks. Unlike permuting the
response directly, this procedure preserves the covariate structure of
the data and supplies a joint null distribution for all of the tested
coefficients from a single set of permutations.

## Usage

``` r
perm_test(
  x,
  full_formula,
  reduced_formula,
  terms,
  reps,
  blocks = NULL,
  direction = "two-sided"
)
```

## Arguments

- x:

  A data frame containing all of the variables referenced by
  `full_formula`, `reduced_formula`, and `blocks`.

- full_formula:

  A two-sided [formula](https://rdrr.io/r/stats/formula.html) giving the
  complete Gaussian linear model, including both the covariates that are
  retained under the null hypothesis and the terms whose coefficients
  are being tested.

- reduced_formula:

  A two-sided [formula](https://rdrr.io/r/stats/formula.html) giving the
  model nested within `full_formula` that contains only the terms
  retained under the null hypothesis. The coefficients that distinguish
  `full_formula` from `reduced_formula` are the coefficients under test.

- terms:

  A character vector naming the coefficients (model-matrix column names)
  under test, in the order they should be reported. Must be exactly the
  set of coefficients in `full_formula` that are not in
  `reduced_formula`.

- reps:

  The number of random permutations to draw. Must be a single positive
  integer.

- blocks:

  An optional bare (unquoted) column name in `x`—or, equivalently, a
  single string naming such a column—giving the exchangeability blocks.
  Residuals are only permuted within blocks. When `NULL` (the default),
  the whole data set is treated as one block, the classical
  Freedman–Lane procedure.

- direction:

  A character string giving the direction of the test. Options are
  `"two-sided"` (default), `"greater"`, or `"less"`. Aliases `"left"`,
  `"right"`, `"both"`, `"two_sided"`, `"two sided"`, and `"two.sided"`
  are also accepted.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per tested coefficient and the columns:

- `term`: The model-matrix column name of the tested coefficient.

- `estimate`: The observed coefficient estimate from the full model.

- `statistic`: The observed studentized (t) statistic for the
  coefficient from the full model.

- `p_value`: The raw, unadjusted permutation p-value.

- `p_value_adj`: The step-down maxT adjusted p-value.

The number of permutations and the test direction are stored as the
`reps` and `direction` attributes of the result.

## Details

The procedure, for each of `reps` random permutations, is:

1.  Fit the reduced model implied by `reduced_formula` and store its
    fitted values and residuals.

2.  Draw one permutation of the observation indices that only reorders
    observations within each level of `blocks` (a single global
    permutation shared by every tested coefficient) and build a permuted
    response as the reduced-model fitted values plus the permuted
    reduced-model residuals.

3.  Fit `full_formula` to the permuted response using infer's native
    Gaussian linear model fitter and record the studentized (t)
    statistics for the tested coefficients.

The observed statistics are the studentized statistics from fitting
`full_formula` to the original response. Because the permuted response
is a real numeric vector that is passed to the native fitter, every fit
carries the same model-matrix structure as the observed fit, so the
coefficient names, formula, and sample correspondence are preserved.

Raw p-values use the finite-sample correction that counts the observed
arrangement as one of the permutations: \$\$p_k = (1 + \\\\b: T\_{bk}
\ge T_k\\) / (reps + 1)\$\$ with absolute statistics for two-sided
tests. The step-down maxT adjusted p-values (Westfall and Young, 1993)
reuse the same permutations: at each step the running maximum of the
(absolute, for two-sided tests) null statistics across all remaining
coefficients supplies the joint null, enforcing multiplicity control
without a Bonferroni correction. The tail comparisons are inclusive, so
equality between a null and an observed statistic counts against the
null.

Every permutation must produce a valid, full-rank fit with finite
studentized statistics; a failure in any permutation raises an error
rather than silently dropping the permutation.

## References

Freedman, D., & Lane, D. (1983). A Nonstochastic Interpretation of
Reported Significance Levels. *Journal of Business & Economic
Statistics*, 1(4), 292–298.
[doi:10.1080/07350015.1983.10509354](https://doi.org/10.1080/07350015.1983.10509354)

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
Testing: Examples and Methods for p-Value Adjustment*. Wiley.

## See also

Other wrapper functions:
[`chisq_stat()`](https://infer.tidymodels.org/dev/reference/chisq_stat.md),
[`chisq_test()`](https://infer.tidymodels.org/dev/reference/chisq_test.md),
[`observe()`](https://infer.tidymodels.org/dev/reference/observe.md),
[`prop_test()`](https://infer.tidymodels.org/dev/reference/prop_test.md),
[`t_stat()`](https://infer.tidymodels.org/dev/reference/t_stat.md),
[`t_test()`](https://infer.tidymodels.org/dev/reference/t_test.md)

## Examples

``` r
# test the partial effects of age and college on hours worked,
# adjusting for the intercept (the only term in the reduced model)
set.seed(1)
perm_test(
  gss,
  full_formula = hours ~ age + college,
  reduced_formula = hours ~ 1,
  terms = c("age", "collegedegree"),
  reps = 100
)
#> # A tibble: 2 × 5
#>   term          estimate statistic p_value p_value_adj
#>   <chr>            <dbl>     <dbl>   <dbl>       <dbl>
#> 1 age            0.00596     0.119   0.921       1    
#> 2 collegedegree  1.53        1.10    0.228       0.376

# adjust for a covariate that is retained under the null and only
# permute residuals within levels of a blocking variable
set.seed(1)
perm_test(
  gss,
  full_formula = hours ~ age + college + sex,
  reduced_formula = hours ~ sex,
  terms = c("age", "collegedegree"),
  reps = 100,
  blocks = sex
)
#> # A tibble: 2 × 5
#>   term          estimate statistic p_value p_value_adj
#>   <chr>            <dbl>     <dbl>   <dbl>       <dbl>
#> 1 age           -0.00108   -0.0222   0.990       1    
#> 2 collegedegree  1.57       1.16     0.228       0.475

# one-sided test
set.seed(1)
perm_test(
  gss,
  full_formula = hours ~ age,
  reduced_formula = hours ~ 1,
  terms = "age",
  reps = 100,
  direction = "greater"
)
#> # A tibble: 1 × 5
#>   term  estimate statistic p_value p_value_adj
#>   <chr>    <dbl>     <dbl>   <dbl>       <dbl>
#> 1 age    0.00781     0.157   0.495       0.495
```
