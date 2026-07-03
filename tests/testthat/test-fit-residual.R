# Fit statistics of Andrich & Marais (2019): the log-residual fit
# statistic, its degrees of freedom, the ANOVA item fit, and the
# class-interval rule.

simd <- function(N, d, seed = 1) {
  set.seed(seed)
  th <- rnorm(N)
  X <- matrix(rbinom(N * length(d), 1, plogis(outer(th, d, "-"))),
              N, length(d))
  colnames(X) <- paste0("I", seq_along(d))
  list(X = X, th = th)
}

test_that("fit residual degrees of freedom reproduce the apportionment formula", {
  s <- simd(400, seq(-1.2, 1.2, length.out = 10), seed = 4)
  fit <- rasch(s$X)
  ne <- sum(!fit$person$extreme)
  I <- ncol(s$X)
  # complete dichotomous data: total df = (N-1)(I-1), apportioned per cell
  f_cell <- (ne - 1) * (I - 1) / (ne * I)
  expect_equal(fit$items$df_fit, rep(f_cell * ne, I), tolerance = 1e-10)
  expect_equal(fit$person$df_fit[!fit$person$extreme],
               rep(f_cell * I, ne), tolerance = 1e-10)
  # extreme persons carry no fit residual and no df
  if (any(fit$person$extreme)) {
    expect_true(all(is.na(fit$person$fit_resid[fit$person$extreme])))
    expect_true(all(is.na(fit$person$df_fit[fit$person$extreme])))
  }
})

test_that("fit residuals are roughly centred and scaled on model-true dichotomous data", {
  s <- simd(700, seq(-2, 2, length.out = 15), seed = 8)
  fit <- rasch(s$X)
  expect_lt(abs(fit$item_fit_summary$mean), 0.6)
  expect_gt(fit$item_fit_summary$sd, 0.3)
  expect_lt(fit$item_fit_summary$sd, 1.8)
  expect_lt(abs(fit$person_fit_summary$mean), 0.6)
  # natural residual present and ordered the same way
  expect_equal(order(fit$items$fit_resid), order(fit$items$natural_resid))
})

test_that("fit residual sign separates over- from under-discrimination", {
  set.seed(21)
  d <- seq(-1.5, 1.5, length.out = 10)
  N <- 600
  th <- rnorm(N)
  X <- matrix(rbinom(N * 10, 1, plogis(outer(th, d, "-"))), N, 10)
  X[, 5] <- rbinom(N, 1, plogis(3 * th - d[5]))  # over-discriminating, on target
  X[, 6] <- rbinom(N, 1, 0.5)                    # noise: under-discriminating
  colnames(X) <- paste0("I", 1:10)
  fit <- rasch(X)
  expect_lt(fit$items$fit_resid[5], -2.5)
  expect_gt(fit$items$fit_resid[6], 2.5)
  # ANOVA item fit flags the noise item most strongly (other items can show
  # compensating artificial misfit, so no specificity claim is made here)
  expect_lt(fit$items$p_anova[6], 0.01)
  expect_setequal(order(fit$items$F_anova, decreasing = TRUE)[1:2], c(5L, 6L))
})

test_that("a Guttman-pattern person is flagged deterministic, a reversed one erratic", {
  set.seed(5)
  d <- seq(-2.5, 2.5, length.out = 12)
  N <- 400
  X <- matrix(rbinom(N * 12, 1, plogis(outer(rnorm(N), d, "-"))), N, 12)
  X[1, ] <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)  # perfect Guttman, score 6
  X[2, ] <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)  # reversed, score 6
  colnames(X) <- paste0("I", 1:12)
  fit <- rasch(X)
  expect_lt(fit$person$fit_resid[1], 0)
  expect_gt(fit$person$fit_resid[2], 2.5)
  expect_lt(fit$person$fit_resid[1], fit$person$fit_resid[2])
})

test_that("the automatic class-interval rule sizes the intervals", {
  expect_equal(.default_n_groups(700), 10L)  # capped at 10
  expect_equal(.default_n_groups(520), 10L)
  expect_equal(.default_n_groups(320), 6L)   # floor(320/50)
  expect_equal(.default_n_groups(120), 2L)
  expect_equal(.default_n_groups(30), 2L)    # never below 2
  s <- simd(320, seq(-1, 1, length.out = 8), seed = 3)
  fit <- rasch(s$X)
  ne <- sum(!fit$person$extreme & !is.na(fit$person$theta))
  expect_equal(fit$n_groups, .default_n_groups(ne))
  expect_lte(max(fit$person$class_interval, na.rm = TRUE), fit$n_groups)
  # an explicit n_groups is honoured
  fit4 <- rasch(s$X, n_groups = 4)
  expect_equal(fit4$n_groups, 4)
  expect_equal(max(fit4$person$class_interval, na.rm = TRUE), 4)
})

test_that("item-trait chi-square df are per item and total df is their sum", {
  s <- simd(450, seq(-1.5, 1.5, length.out = 9), seed = 12)
  X <- s$X
  X[sample(length(X), 300)] <- NA   # missing data can starve intervals
  fit <- rasch(X)
  expect_true(all(fit$item_trait$df >= 1))
  expect_true(all(fit$item_trait$df <= fit$n_groups - 1))
  expect_equal(fit$total_df, sum(fit$item_trait$df))
  expect_true(all(c("p_bonf") %in% names(fit$item_trait)))
})

test_that("ANOVA item fit is calibrated under the model", {
  s <- simd(600, seq(-1.8, 1.8, length.out = 12), seed = 9)
  fit <- rasch(s$X)
  expect_gt(mean(fit$items$F_anova, na.rm = TRUE), 0.5)
  expect_lt(mean(fit$items$F_anova, na.rm = TRUE), 1.6)
  expect_true(all(c("df1", "df2", "p_adj", "p_bonf") %in% names(fit$item_anova)))
  expect_equal(fit$item_anova$df1, rep(fit$n_groups - 1L, ncol(s$X)))
})

test_that("polytomous and missing-data analyses carry coherent fit output", {
  set.seed(11)
  tau <- list(c(-1.5, -0.5, 0.5), c(-0.8, 0.8), c(-0.7, 0.5, 1.7),
              c(-0.5, 0.5), c(-1, 1), c(-1, 0))
  X <- sapply(tau, function(tt) vapply(rnorm(350), function(b) {
    p <- item_moments(b, tt)$P
    sample(0:(length(p) - 1), 1, prob = p)
  }, 0L))
  colnames(X) <- paste0("Q", 1:6)
  X[sample(length(X), 200)] <- NA
  fit <- rasch(X)
  ne <- !fit$person$extreme
  # df vary with each person's observed cells
  expect_gt(length(unique(na.omit(round(fit$person$df_fit, 6)))), 3)
  # item df equals cell df times observed non-extreme responses
  obs_i <- colSums(!is.na(X) & ne)
  ratio <- unname(fit$items$df_fit / obs_i)
  expect_equal(ratio, rep(ratio[1], 6), tolerance = 1e-10)
  expect_true(all(is.finite(fit$items$fit_resid)))
})

test_that("class intervals never split persons sharing a location", {
  s <- simd(430, seq(-1.8, 1.8, length.out = 9), seed = 21)
  fit <- rasch(s$X)
  th <- fit$person$theta; ci <- fit$person$class_interval
  ok <- !is.na(th) & !is.na(ci)
  # every distinct location maps to exactly one interval
  per_loc <- tapply(ci[ok], th[ok], function(g) length(unique(g)))
  expect_true(all(per_loc == 1))
  # intervals are contiguous in location and cover the requested count
  ord <- order(th[ok])
  expect_true(all(diff(ci[ok][ord]) >= 0))
  expect_equal(max(ci, na.rm = TRUE), fit$n_groups)
  # sizes are as equal as tie-preservation allows (no interval empty)
  expect_true(all(tabulate(ci[ok]) > 0))
})
