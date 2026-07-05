# Output tables: chi-square class-interval detail, score-table
# estimators and extreme-score extrapolation, the LR test block, and the
# summary distribution statistics.

simd2 <- function(N, d, seed = 1) {
  set.seed(seed)
  th <- rnorm(N)
  X <- matrix(rbinom(N * length(d), 1, plogis(outer(th, d, "-"))),
              N, length(d))
  colnames(X) <- paste0("I", seq_along(d))
  X
}

test_that("chisq_detail reproduces the item-trait chi-square from its components", {
  fit <- rasch(simd2(400, seq(-1.5, 1.5, length.out = 8), seed = 2))
  cd <- chisq_detail(fit, "I4")
  expect_equal(sum(cd$intervals$chisq[cd$intervals$used]),
               fit$item_trait$chisq[4], tolerance = 1e-10)
  expect_equal(cd$df, fit$item_trait$df[4])
  # ES is the residual with the interval-size factor removed
  used <- cd$intervals$used
  expect_equal(cd$intervals$residual[used] / sqrt(cd$intervals$n[used]),
               cd$intervals$es[used], tolerance = 1e-10)
  # per-interval category proportions sum to 1; dichotomous OBS.T = OBS.P(1)
  cats <- cd$categories
  for (g in unique(cats$interval)) {
    rows <- cats$interval == g
    if (any(is.na(cats$obs_p[rows]))) next
    expect_equal(sum(cats$obs_p[rows]), 1, tolerance = 1e-10)
    expect_equal(sum(cats$est_p[rows]), 1, tolerance = 1e-8)
    expect_equal(cats$obs_t[rows][2], cats$obs_p[rows][2], tolerance = 1e-10)
  }
  # interval means increase along the trait
  expect_true(all(diff(na.omit(cd$intervals$theta_mean)) > 0))
})

test_that("score_table extrapolates extremes by the geometric rule", {
  fit <- rasch(simd2(500, seq(-2, 2, length.out = 10), seed = 6))
  # MLE table: infinite extremes are NA, interior solves the score equation
  mle <- score_table(fit, method = "mle")
  expect_true(is.na(mle$theta[1]) && is.na(mle$theta[11]))
  r5 <- mle$theta[mle$score == 5]
  expect_equal(sum(vapply(fit$tau_list, function(tt)
    item_moments(r5, tt)$E, 0)), 5, tolerance = 1e-6)
  # WLE is finite everywhere and flatter than MLE (shrunk towards centre)
  wle <- score_table(fit)
  expect_true(all(is.finite(wle$theta)))
  expect_lt(wle$theta[10], mle$theta[10])
  # geometric extrapolation: d_next = b^2 / a on the preceding differences
  ex <- score_table(fit, method = "mle", extremes = "extrapolated")
  th <- mle$theta
  b <- th[10] - th[9]; a <- th[9] - th[8]
  expect_equal(ex$theta[11], th[10] + b^2 / a, tolerance = 1e-10)
  b0 <- th[3] - th[2]; a0 <- th[4] - th[3]
  expect_equal(ex$theta[1], th[2] - b0^2 / a0, tolerance = 1e-10)
  expect_true(all(ex$extrapolated[c(1, 11)]))
  # extrapolated SE exceeds the SE one score in
  expect_gt(ex$se[11], ex$se[10])
  # frequencies count complete responders
  expect_equal(sum(ex$freq), sum(stats::complete.cases(fit$X)))
  expect_equal(ex$cum_pct[11], 100)
})

test_that("lr_test prefers PCM when thresholds differ and RSM when they do not", {
  set.seed(9)
  simpoly <- function(taus) {
    X <- sapply(taus, function(tt) vapply(rnorm(400), function(b)
      sample(0:2, 1, prob = item_moments(b, tt)$P), 0L))
    colnames(X) <- paste0("Q", seq_along(taus)); X
  }
  d <- seq(-1, 1, length.out = 6)
  # common threshold structure: the adjusted test retains RSM; the raw
  # composite statistic is anticonservative (that is why it is adjusted)
  same <- lapply(d, function(dd) c(-0.8, 0.8) + dd)
  lr1 <- lr_test(rasch(simpoly(same)))
  expect_gt(lr1$p_adj, 0.01)
  expect_equal(lr1$df, 5)  # 6 items x 2 thresholds - 1 vs (6 - 1) + 1
  expect_lt(lr1$chisq_adj, lr1$chisq)
  # the Godambe eigenvalues carry the inflation of the composite statistic
  expect_gt(mean(lr1$lambda), 1)
  # item-specific spreads: PCM needed, decisively
  diff_ <- lapply(seq_along(d), function(i) c(-0.2, 0.2) * i + d[i])
  lr2 <- lr_test(rasch(simpoly(diff_)))
  expect_lt(lr2$p_adj, 0.001)
  expect_gt(lr2$chisq, 0)
  # guards
  expect_error(lr_test(rasch(simd2(150, c(-1, 0, 1)))), "dichotomous")
})

test_that("summary blocks carry the distribution statistics", {
  X <- simd2(400, seq(-1.5, 1.5, length.out = 8), seed = 13)
  fit <- rasch(X)
  for (blk in list(fit$item_fit_summary, fit$person_fit_summary,
                   fit$summary_stats$item_location,
                   fit$summary_stats$person_location_noext)) {
    expect_true(all(c("mean", "sd", "skewness", "kurtosis") %in% names(blk)))
    expect_true(is.finite(blk$sd))
  }
  expect_equal(fit$summary_stats$item_location$mean, 0, tolerance = 1e-8)
  expect_true(is.finite(fit$summary_stats$cor_item_fit_location))
  expect_true(fit$summary_stats$df_factor < 1 && fit$summary_stats$df_factor > 0.8)
  # item separation index behaves like a reliability
  expect_true(fit$isi$PSI >= 0 && fit$isi$PSI <= 1)
  # spread-out items with a large sample are well separated
  expect_gt(fit$isi$PSI, 0.9)
  # alpha applicability flag
  expect_true(fit$alpha$applicable)
  Xm <- X; Xm[1, 1] <- NA
  expect_false(rasch(Xm)$alpha$applicable)
})

test_that("report_html writes a complete self-contained report", {
  set.seed(3)
  d <- seq(-1.5, 1.5, length.out = 6)
  X <- matrix(rbinom(250 * 6, 1, plogis(outer(rnorm(250), d, "-"))), 250, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g = rep(c("a", "b"), each = 125)), factors = "g")
  out <- file.path(tempdir(), "rmt_report_test.html")
  on.exit(unlink(out), add = TRUE)
  report_html(fit, out, title = "Test report")
  expect_true(file.exists(out))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  # self-contained: plots embedded, no external references
  expect_gt(lengths(regmatches(html, gregexpr("data:image/png;base64", html))), 5)
  expect_false(grepl("src=\"http", html))
  for (sec in c("Summary", "Item statistics", "Thresholds", "Score to measure",
                "Dimensionality", "Local dependence", "Classical companions",
                "Differential item functioning", "Person estimates"))
    expect_true(grepl(paste0("<h2>", sec), html), label = sec)
  # the base-R base64 encoder matches the RFC 4648 test vectors
  enc <- function(txt) {
    p <- tempfile(); on.exit(unlink(p), add = TRUE)
    writeBin(charToRaw(txt), p)
    .b64(p)
  }
  expect_equal(enc("Man"), "TWFu")
  expect_equal(enc("Ma"), "TWE=")
  expect_equal(enc("M"), "TQ==")
  expect_equal(enc("foobar"), "Zm9vYmFy")
})

test_that("the fit and targeting summaries are complete tidy tables", {
  set.seed(1)
  d <- seq(-2, 2, length.out = 6)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  f <- rasch(X)
  ft <- fit_summary_table(f)
  tt <- targeting_table(f)
  expect_identical(names(ft), c("statistic", "value"))
  expect_identical(names(tt), c("statistic", "value"))
  expect_false(anyNA(ft$value))
  # spot-check against the fit object
  expect_equal(ft$value[ft$statistic == "Model"], f$model)
  expect_lt(abs(as.numeric(ft$value[ft$statistic == "Total item-trait chi-square"]) -
                f$total_chisq), 5e-4)
  expect_lt(abs(as.numeric(tt$value[tt$statistic == "PSI"]) - f$psi$PSI), 5e-4)
  expect_lt(abs(as.numeric(tt$value[tt$statistic == "Coefficient alpha"]) -
                f$alpha$alpha), 5e-4)
  # robust when alpha is not applicable (missing data)
  Xm <- X; Xm[1:150, 1] <- NA; Xm[151:300, 6] <- NA
  fm <- rasch(Xm)
  ttm <- targeting_table(fm)
  expect_true("NA" %in% ttm$value[ttm$statistic == "Coefficient alpha"] ||
              is.finite(as.numeric(ttm$value[ttm$statistic == "Coefficient alpha"])))
  expect_no_error(fit_summary_table(fm))
})
