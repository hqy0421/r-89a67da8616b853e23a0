# Small synthetic data with a covariate and exchangeability blocks
make_perm_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z = rnorm(n),
    blk = rep(c("a", "b", "c"), length.out = n)
  )
  d$y <- 1 + 0.4 * d$x1 + 0 * d$x2 + 0.3 * d$z + rnorm(n, sd = 0.8)
  d
}

d <- make_perm_data()
full_f <- y ~ z + x1 + x2
red_f <- y ~ z
tested <- c("x1", "x2")

test_that("perm_test returns the expected structure", {
  set.seed(1)
  res <- perm_test(d, full_f, red_f, tested, reps = 49, blocks = blk)

  expect_s3_class(res, "tbl_df")
  expect_named(res, c("term", "estimate", "statistic", "p_value", "p_value_adj"))
  expect_equal(nrow(res), length(tested))
  expect_equal(res$term, tested)
  expect_true(all(vapply(res, is.numeric, logical(1))[-1]))
  expect_equal(attr(res, "reps"), 49)
  expect_equal(attr(res, "direction"), "both")
})

test_that("observed statistics come from the native Gaussian fit", {
  set.seed(1)
  res <- perm_test(d, full_f, red_f, tested, reps = 11, blocks = blk)

  native <- broom::tidy(stats::glm(full_f, data = d, family = stats::gaussian()))
  idx <- match(tested, native$term)

  expect_equal(res$estimate, native$estimate[idx])
  expect_equal(res$statistic, native$statistic[idx])
})

test_that("permutations are shared and restricted to within blocks", {
  reps <- 25L
  set.seed(7)
  res <- perm_test(d, full_f, red_f, tested, reps = reps, blocks = blk)

  # Recreate the package's RNG stream (same seed, same sample() ordering).
  set.seed(7)
  reduced <- stats::glm(red_f, data = d, family = stats::gaussian())
  fitted <- stats::fitted(reduced)
  resid <- stats::residuals(reduced)
  block_rows <- split(seq_len(nrow(d)), as.integer(as.factor(d$blk)))

  null_t <- matrix(NA_real_, reps, length(tested), dimnames = list(NULL, tested))
  for (b in seq_len(reps)) {
    perm <- integer(nrow(d))
    for (rows in block_rows) {
      perm[rows] <- sample(rows, size = length(rows), replace = FALSE)
    }

    # the permutation never moves an observation outside its block
    for (rows in block_rows) {
      expect_equal(sort(perm[rows]), sort(rows))
    }

    dd <- d
    dd$y <- fitted + resid[perm]
    tt <- broom::tidy(stats::glm(full_f, data = dd, family = stats::gaussian()))
    null_t[b, ] <- tt$statistic[match(tested, tt$term)]
  }

  obs <- broom::tidy(stats::glm(full_f, data = d, family = stats::gaussian()))
  obs_t <- obs$statistic[match(tested, obs$term)]

  raw <- (colSums(abs(null_t) >= rep(abs(obs_t), each = reps)) + 1) / (reps + 1)
  expect_equal(res$p_value, unname(raw))
})

test_that("finite-sample correction counts the observed arrangement", {
  set.seed(11)
  strong <- d
  strong$y <- 1 + 5 * d$x1 + 5 * d$x2 + rnorm(nrow(d), sd = 0.01)
  reps <- 39L

  res <- perm_test(strong, y ~ x1 + x2, y ~ 1, tested, reps = reps)

  # the smallest attainable p-value is 1 / (reps + 1), never zero
  expect_equal(min(res$p_value), 1 / (reps + 1))
  expect_equal(min(res$p_value_adj), 1 / (reps + 1))
  expect_true(all(res$p_value > 0))
})

test_that("step-down maxT is monotone and dominates the raw p-values", {
  set.seed(2)
  res <- perm_test(d, full_f, red_f, tested, reps = 79, blocks = blk)

  expect_true(all(res$p_value_adj >= res$p_value - 1e-12))
  expect_true(all(res$p_value_adj <= 1))

  # adjusted p-values are monotone in the observed statistic: the smaller the
  # |statistic|, the larger (or equal) the adjusted p-value
  ord <- order(abs(res$statistic), decreasing = TRUE)
  expect_true(all(diff(res$p_value_adj[ord]) >= -1e-12))
})

test_that("step-down maxT matches a hand-rolled implementation", {
  reps <- 99L
  set.seed(4)
  res <- perm_test(d, full_f, red_f, tested, reps = reps, blocks = blk)

  set.seed(4)
  reduced <- stats::glm(red_f, data = d, family = stats::gaussian())
  block_rows <- split(seq_len(nrow(d)), as.integer(as.factor(d$blk)))
  null_t <- matrix(NA_real_, reps, length(tested))
  for (b in seq_len(reps)) {
    perm <- integer(nrow(d))
    for (rows in block_rows) perm[rows] <- sample(rows, length(rows))
    dd <- d
    dd$y <- stats::fitted(reduced) + stats::residuals(reduced)[perm]
    tt <- broom::tidy(stats::glm(full_f, data = dd, family = stats::gaussian()))
    null_t[b, ] <- tt$statistic[match(tested, tt$term)]
  }
  obs <- broom::tidy(stats::glm(full_f, data = d, family = stats::gaussian()))
  obs_t <- obs$statistic[match(tested, obs$term)]

  tmax <- apply(abs(null_t), 1, max)
  o <- order(abs(obs_t), decreasing = TRUE)
  expected_adj <- numeric(length(tested))
  previous <- 0
  for (j in o) {
    previous <- max(previous, (sum(tmax >= abs(obs_t[j])) + 1) / (reps + 1))
    expected_adj[j] <- previous
  }

  expect_equal(res$p_value_adj, expected_adj)
})

test_that("directional tests use signed statistics", {
  set.seed(5)
  n <- 80
  ds <- data.frame(x = rnorm(n))
  ds$y <- 0.8 * ds$x + rnorm(n)

  greater <- perm_test(ds, y ~ x, y ~ 1, "x", reps = 99, direction = "greater")
  less <- perm_test(ds, y ~ x, y ~ 1, "x", reps = 99, direction = "less")
  two <- perm_test(ds, y ~ x, y ~ 1, "x", reps = 99)

  # a positive slope is strong evidence for "greater" but not "less"
  expect_lt(greater$p_value, less$p_value)
  # with one term there is nothing to adjust, so raw and adjusted coincide
  expect_equal(greater$p_value, greater$p_value_adj)
  expect_equal(two$p_value, two$p_value_adj)
  expect_equal(attr(less, "direction"), "left")
})

test_that("bare-name and string block arguments agree", {
  set.seed(9)
  bare <- perm_test(d, full_f, red_f, tested, reps = 40, blocks = blk)
  set.seed(9)
  string <- perm_test(d, full_f, red_f, tested, reps = 40, blocks = "blk")

  expect_equal(bare, string)
})

test_that("factor coefficient names are preserved", {
  set.seed(9)
  res <- perm_test(
    gss,
    hours ~ age + college + sex,
    hours ~ age,
    c("collegedegree", "sexfemale"),
    reps = 30,
    blocks = sex
  )

  expect_equal(res$term, c("collegedegree", "sexfemale"))
})

test_that("raw p-values are calibrated under the null", {
  set.seed(20261119)
  sims <- 40
  reps <- 99
  raw <- matrix(NA_real_, sims, length(tested))

  for (s in seq_len(sims)) {
    ds <- make_perm_data(seed = 1000 + s)
    ds$y <- 1 + 0.3 * ds$z + rnorm(nrow(ds)) # x1 and x2 have no effect
    r <- perm_test(ds, full_f, red_f, tested, reps = reps, blocks = blk)
    raw[s, ] <- r$p_value
  }

  # discrete permutation p-values contain ties, which ks.test warns about
  expect_gt(suppressWarnings(stats::ks.test(raw[, 1], "punif")$p.value), 0.05)
  expect_gt(suppressWarnings(stats::ks.test(raw[, 2], "punif")$p.value), 0.05)
})

test_that("perm_test fails on rank deficiency", {
  bad <- d
  bad$x2 <- bad$x1

  expect_snapshot(
    error = TRUE,
    perm_test(bad, full_f, red_f, tested, reps = 10, blocks = blk)
  )
})

test_that("perm_test fails when the reduced model is not nested", {
  expect_snapshot(
    error = TRUE,
    perm_test(d, y ~ x1 + z, y ~ x2, "x1", reps = 10)
  )
})

test_that("perm_test fails when the formulas have different responses", {
  expect_snapshot(
    error = TRUE,
    perm_test(d, y ~ x1, z ~ 1, "x1", reps = 10)
  )
})

test_that("perm_test fails on mismatched terms", {
  expect_snapshot(
    error = TRUE,
    perm_test(d, full_f, red_f, "nope", reps = 10)
  )
  expect_snapshot(
    error = TRUE,
    perm_test(d, full_f, red_f, c("x1", "z"), reps = 10)
  )
  expect_snapshot(
    error = TRUE,
    perm_test(d, full_f, red_f, "x1", reps = 10)
  )
})

test_that("perm_test fails on invalid data", {
  d_na <- d
  d_na$x1[3] <- NA
  expect_snapshot(
    error = TRUE,
    perm_test(d_na, full_f, red_f, tested, reps = 10)
  )

  d_inf <- d
  d_inf$z[2] <- Inf
  expect_snapshot(
    error = TRUE,
    perm_test(d_inf, full_f, red_f, tested, reps = 10)
  )

  d_cat <- d
  d_cat$y <- sample(c("a", "b"), nrow(d_cat), replace = TRUE)
  expect_snapshot(
    error = TRUE,
    perm_test(d_cat, full_f, red_f, tested, reps = 10)
  )
})

test_that("perm_test fails on degenerate block structure", {
  d_singleton <- d
  d_singleton$id <- as.character(seq_len(nrow(d_singleton)))

  expect_snapshot(
    error = TRUE,
    perm_test(d_singleton, y ~ x1, y ~ 1, "x1", reps = 10, blocks = id)
  )
  expect_snapshot(
    error = TRUE,
    perm_test(d, y ~ x1, y ~ 1, "x1", reps = 10, blocks = absent)
  )
})

test_that("perm_test fails on invalid scalar arguments", {
  expect_snapshot(
    error = TRUE,
    perm_test(d, full_f, red_f, tested, reps = 0)
  )
  expect_snapshot(
    error = TRUE,
    perm_test(d, full_f, red_f, tested, reps = 3.5)
  )
  expect_snapshot(
    error = TRUE,
    perm_test(d, full_f, red_f, tested, reps = 10, direction = "sideways")
  )
  expect_snapshot(
    error = TRUE,
    perm_test(d, ~ x1, red_f, tested, reps = 10)
  )
  expect_snapshot(
    error = TRUE,
    perm_test(list(a = 1), full_f, red_f, tested, reps = 10)
  )
  expect_snapshot(
    error = TRUE,
    perm_test(d, full_f, red_f, c("x1", "x1"), reps = 10)
  )
})
