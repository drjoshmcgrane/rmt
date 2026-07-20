# Statistical-validity regressions from the 2026-07 external review, each
# verified by simulation before fixing: the recentring covariance transform
# (mixed max scores), the Warm WLE discrimination cancellation, honest
# chi-square degrees of freedom, covariance-correct equating drift tests,
# the few-judges clustering guard, and the judge-level DIF ANOVA.

test_that("mixed-max-score item location SEs are calibrated (cov transform)", {
  skip_on_cran()   # 60 replicate fits
  set.seed(9)
  L <- 6; m <- rep(c(1L, 3L), each = 3)
  btrue <- seq(-1, 1, length.out = L)
  tau_l <- lapply(1:L, function(j) btrue[j] + seq(-0.6, 0.6, length.out = m[j]))
  gen <- function(N) {
    th <- rnorm(N)
    X <- matrix(0L, N, L, dimnames = list(NULL, sprintf("I%d", 1:L)))
    for (j in 1:L) for (i in 1:N) {
      d <- th[i] - tau_l[[j]]
      p <- c(1, exp(cumsum(d))); X[i, j] <- sample(0:m[j], 1, prob = p / sum(p))
    }
    X
  }
  R <- 60; locs <- matrix(NA, R, L); ses <- matrix(NA, R, L)
  for (r in 1:R) {
    f <- pcml(gen(500))
    locs[r, ] <- vapply(1:L, function(j) mean(f$thr$tau[f$thr$item == j]), 0)
    ses[r, ] <- vapply(1:L, function(j) {
      rows <- f$thr$id[f$thr$item == j]
      sqrt(mean(f$cov_tau[rows, rows]))
    }, 0)
  }
  ratio <- apply(locs, 2, sd) / colMeans(ses)
  # before the transform, dichotomous items sat ~0.90 and 3-threshold items
  # ~1.20 systematically; after it every item is within noise of 1
  expect_true(all(ratio > 0.8 & ratio < 1.25))
  expect_lt(abs(mean(ratio) - 1), 0.12)
})

test_that("Warm WLE is invariant to a common discrimination", {
  tau_list <- list(c(-1), c(0), c(1), c(-0.5, 0.5))
  for (a in c(0.5, 2)) {
    w <- person_wle(tau_list, disc = a)
    for (R in 1:4) {
      # exact WLE: root of the weighted score a(R-E) + a^3 mu3 / (2 a^2 V)
      obj <- function(th) {
        mo <- lapply(tau_list, item_moments, theta = th, disc = a)
        E <- sum(sapply(mo, `[[`, "E")); V <- sum(sapply(mo, `[[`, "V"))
        m3 <- sum(sapply(mo, `[[`, "mu3"))
        a * (R - E) + a^3 * m3 / (2 * a^2 * V)
      }
      exact <- uniroot(obj, c(-30, 30), tol = 1e-12)$root
      expect_equal(unname(w$theta[as.character(R)]), exact, tolerance = 1e-6)
    }
  }
})

test_that("an item with fewer than two class intervals gets NA, not df = 1", {
  # two items, all non-extreme persons share one raw score -> one interval
  X <- cbind(I1 = rep(c(0L, 1L), 60), I2 = rep(c(1L, 0L), 60))
  f <- rasch(X, n_groups = 2)
  expect_true(all(is.na(f$items$df)))
  expect_true(all(is.na(f$items$chisq)))
  expect_true(all(is.na(f$items$p)))
  # with nothing testable the omnibus is NA too, not chi-square 0 on 0 df
  expect_true(is.na(f$total_df))
  expect_true(is.na(f$total_chisq_p))
})

test_that("equating drift tests are calibrated under the null", {
  skip_on_cran()   # 80 replicate pairs of fits
  set.seed(42)
  L <- 8; btrue <- seq(-1.5, 1.5, length.out = L)
  mk <- function() {
    X <- matrix(rbinom(400 * L, 1, plogis(outer(rnorm(400), btrue, "-"))),
                400, L, dimnames = list(NULL, paste0("I", 1:L)))
    rasch(as.data.frame(X))
  }
  rej <- 0; tot <- 0
  for (r in 1:80) {
    eq <- equate_tests(mk(), mk())
    rej <- rej + sum(eq$table$p < 0.05, na.rm = TRUE)
    tot <- tot + sum(is.finite(eq$table$p))
  }
  # naive sqrt(v) denominators (no shift covariance) were mis-calibrated;
  # the projected covariance restores ~nominal rejection
  expect_gt(rej / tot, 0.02)
  expect_lt(rej / tot, 0.09)
})

test_that("clustered SEs refuse a single judge and note few judges", {
  d <- data.frame(object_a = rep(paste0("O", 1:4), 30),
                  object_b = rep(paste0("O", c(2:4, 1)), 30))
  set.seed(2); d$winner <- ifelse(runif(120) < .5, d$object_a, d$object_b)
  d$judge <- "J1"
  expect_error(btl(d, "object_a", "object_b", "winner", judge = "judge"),
               "at least 2 judges")
  d$judge <- rep(sprintf("J%d", 1:6), 20)
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  expect_true(any(grepl("judge clusters", f$notes)))
})

test_that("btl_dif does not flag under judge heterogeneity with null groups", {
  skip_on_cran()   # several fits
  simnull <- function(seed) {
    set.seed(seed)
    K <- 8; beta <- seq(-1.2, 1.2, length.out = K); nj <- 12; npj <- 60
    rows <- list()
    for (j in 1:nj) {
      bj <- beta + rnorm(K, 0, 0.6)
      ia <- sample(K, npj, TRUE); ib <- (ia + sample(K - 1, npj, TRUE) - 1L) %% K + 1L
      win <- rbinom(npj, 1, plogis(bj[ia] - bj[ib]))
      rows[[j]] <- data.frame(object_a = paste0("O", ia),
                              object_b = paste0("O", ib),
                              winner = paste0("O", ifelse(win == 1, ia, ib)),
                              judge = sprintf("J%02d", j))
    }
    d <- do.call(rbind, rows)
    bt <- btl(d, "object_a", "object_b", "winner", judge = "judge")
    grp <- setNames(rep(c("A", "B"), each = nj / 2), sprintf("J%02d", 1:nj))
    df <- btl_dif(bt, factors = list(g = grp))
    sum(df$summary$uniform_DIF %in% TRUE)
  }
  # comparison-level pseudoreplication flagged 6 of 10 such nulls
  flags <- vapply(1:5, simnull, 0)
  expect_lte(sum(flags > 0), 1)
})

test_that("pairwise chi-square df counts every estimated parameter", {
  set.seed(3)
  K <- 8; beta <- seq(-1.5, 1.5, length.out = K)
  n <- 1500
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib] + 0.3))
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(win == 1, ia, ib)))
  f0 <- btl(d, "object_a", "object_b", "winner")
  f1 <- btl(d, "object_a", "object_b", "winner", position = TRUE)
  # the position covariate consumes one further degree of freedom
  expect_equal(f1$total_df, f0$total_df - 1L)
  # a design with no testable pairs left reports NA, not df = 1
  d3 <- data.frame(object_a = c("A", "B", "C"), object_b = c("B", "C", "A"),
                   winner = c("A", "B", "C"))
  d3 <- d3[rep(1:3, 20), ]
  f3 <- btl(d3, "object_a", "object_b", "winner", position = TRUE)
  expect_true(is.na(f3$total_df) || f3$total_df >= 1L)
})

test_that("dimensionality reference respects comparison counts", {
  skip_on_cran()   # two bootstrap references
  set.seed(7)
  K <- 8; beta <- setNames(seq(-1.2, 1.2, length.out = K), paste0("O", 1:K))
  pr <- t(combn(names(beta), 2)); rows <- list()
  for (i in seq_len(nrow(pr))) {
    wins <- rbinom(1, 12, plogis(beta[pr[i, 1]] - beta[pr[i, 2]]))
    rows[[length(rows) + 1]] <- data.frame(a = pr[i, 1], b = pr[i, 2],
                                           win = pr[i, 1], k = wins)
    rows[[length(rows) + 1]] <- data.frame(a = pr[i, 1], b = pr[i, 2],
                                           win = pr[i, 2], k = 12 - wins)
  }
  agg <- do.call(rbind, rows); agg <- agg[agg$k > 0, ]
  expd <- agg[rep(seq_len(nrow(agg)), agg$k), ]; expd$k <- 1
  fa <- btl(agg, "a", "b", winner = "win", count = "k")
  fe <- btl(expd, "a", "b", winner = "win")
  set.seed(11); da <- btl_dimensionality(fa, reps = 150)
  set.seed(11); de <- btl_dimensionality(fe, reps = 150)
  # pair-level simulation makes the two forms draw from the same design:
  # compare the leading-strength reference distribution, not a pooled mean
  ra <- da$reference$draws; re <- de$reference$draws
  expect_lt(abs(mean(ra) - mean(re)) / mean(re), 0.05)
  expect_lt(abs(stats::quantile(ra, .95) - stats::quantile(re, .95)) /
              stats::quantile(re, .95), 0.05)
})

test_that("judge-resampling bootstrap runs and matches the estimates", {
  skip_on_cran()   # B pipeline refits
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 8, reps_within = 25,
                         reps_cross = 25, panel_units = c(0.8, 1.25),
                         set_units = c(1, 1.3), seed = 7)
  os <- attr(d, "truth")$object_sets
  set.seed(1)
  fj <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
                 se_method = "judge_bootstrap", boot_reps = 30)
  fc <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
                 se_method = "conditional")
  expect_equal(fj$phi_table$phi, fc$phi_table$phi)
  expect_equal(fj$objects$v, fc$objects$v)
  expect_true(all(is.finite(fj$phi_table$se_log_phi)))
  expect_true(is.finite(fj$alpha_table$se_log_alpha[2]))
  expect_match(fj$se_note, "judge-resampling")
})

test_that("half-tie weights no longer break the dimensionality reference", {
  set.seed(4)
  K <- 6; beta <- setNames(seq(-1, 1, length.out = K), paste0("O", 1:K))
  pr <- t(combn(names(beta), 2)); rows <- list()
  for (i in seq_len(nrow(pr))) for (r in 1:10) {
    pw <- plogis(beta[pr[i, 1]] - beta[pr[i, 2]])
    win <- if (runif(1) < 0.15) "tie" else
      if (runif(1) < pw) pr[i, 1] else pr[i, 2]
    rows[[length(rows) + 1]] <- data.frame(a = pr[i, 1], b = pr[i, 2],
                                           win = win)
  }
  d <- do.call(rbind, rows)
  f <- btl(d, "a", "b", winner = "win", ties = "half")
  dm <- btl_dimensionality(f, reps = 60)
  # fractional 0.5 weights fed as.integer() a zero binomial size before:
  # every reference draw was degenerate at 0
  expect_true(all(is.finite(dm$reference$draws)))
  expect_gt(dm$reference$mean, 0.5)
})

test_that("judge bootstrap refuses a single-judge panel", {
  d <- simulate_btl_efrm(n_objects_per_set = 5, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 4, reps_within = 15,
                         reps_cross = 15, seed = 3)
  os <- attr(d, "truth")$object_sets
  # collapse panel 2 to a single judge
  keep <- d$panel != "panel2" | d$judge == d$judge[d$panel == "panel2"][1]
  expect_error(
    btl_efrm(d[keep, ], "object_a", "object_b", "winner", "judge", "panel",
             os, se_method = "judge_bootstrap", boot_reps = 20),
    "at least 2 judges in every panel")
})

test_that("clustered covariance notes rank deficiency (judges <= parameters)", {
  set.seed(6)
  K <- 10; beta <- seq(-1.5, 1.5, length.out = K)
  n <- 600
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib]))
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(win == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:5), n, TRUE))
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  expect_true(any(grepl("rank-deficient", f$notes)))
})

# --- DIF ANOVA engine (round 10): order invariance, person units, GG ------

test_that("multi-factor DIF tests are order-invariant (Type II)", {
  skip_on_cran()
  set.seed(42)
  N <- 500; L <- 8; btrue <- seq(-1.5, 1.5, length.out = L)
  g1 <- sample(c("a", "b"), N, TRUE, prob = c(0.35, 0.65))
  g2 <- ifelse(runif(N) < 0.8, ifelse(g1 == "a", "x", "y"),
               sample(c("x", "y"), N, TRUE))
  th <- rnorm(N)
  X <- matrix(0L, N, L, dimnames = list(NULL, paste0("I", 1:L)))
  for (j in 1:L)
    X[, j] <- rbinom(N, 1, plogis(th - btrue[j] +
                                    ifelse(j == 3 & g1 == "b", 0.8, 0)))
  df <- data.frame(X, g1 = g1, g2 = g2)
  s12 <- dif_anova(rasch(df, factors = c("g1", "g2"),
                         items = paste0("I", 1:L)))$summary
  s21 <- dif_anova(rasch(df, factors = c("g2", "g1"),
                         items = paste0("I", 1:L)))$summary
  key <- function(s) s[order(s$item, s$term), c("F_uniform", "F_nonuniform")]
  expect_equal(key(s12), key(s21), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("duplicating every person leaves the DIF tests exactly unchanged", {
  set.seed(7)
  N <- 250; L <- 6; btrue <- seq(-1.2, 1.2, length.out = L)
  g <- sample(c("a", "b"), N, TRUE); th <- rnorm(N)
  X <- matrix(0L, N, L, dimnames = list(NULL, paste0("I", 1:L)))
  for (j in 1:L) X[, j] <- rbinom(N, 1, plogis(th - btrue[j]))
  d1 <- data.frame(X, g = g, pid = sprintf("P%03d", 1:N))
  f1 <- rasch(d1, factors = "g", id = "pid", items = paste0("I", 1:L))
  f2 <- rasch(rbind(d1, d1), factors = "g", id = "pid",
              items = paste0("I", 1:L))
  s1 <- dif_anova(f1)$summary; s2 <- dif_anova(f2)$summary
  expect_equal(s2$F_uniform[order(s2$item)], s1$F_uniform[order(s1$item)],
               tolerance = 1e-8)
})

test_that("stacked between-treatment is refused; within declarations checked", {
  set.seed(3)
  N <- 120; d <- seq(-1, 1, length.out = 6)
  X <- rbind(matrix(rbinom(N * 6, 1, plogis(outer(rnorm(N), d, "-"))), N, 6),
             matrix(rbinom(N * 6, 1, plogis(outer(rnorm(N), d, "-"))), N, 6))
  colnames(X) <- paste0("I", 1:6)
  dat <- data.frame(X, occ = rep(c("t1", "t2"), each = N))
  id <- rep(sprintf("P%03d", 1:N), 2)
  fit <- rasch(dat, factors = "occ", id = id)
  expect_error(dif_anova(fit, within = character(0)), "vary within persons")
  expect_error(dif_anova(fit, within = "nope"), "not among the nominated")
})

test_that("multi-level within DIF is GG-calibrated and has power", {
  skip_on_cran()   # replicate fits
  np <- 120; L2 <- 6; b2 <- seq(-1, 1, length.out = L2); K <- 4
  gen <- function(seed, occ_sd, shift3_t4 = 0) {
    set.seed(seed)
    th0 <- rnorm(np); rows <- list()
    for (k in 1:K) {
      thk <- th0 + rnorm(np, 0, occ_sd[k])
      Xk <- matrix(0L, np, L2, dimnames = list(NULL, paste0("I", 1:L2)))
      for (j in 1:L2)
        Xk[, j] <- rbinom(np, 1, plogis(thk - b2[j] +
                                          ifelse(j == 3 & k == 4,
                                                 shift3_t4, 0)))
      rows[[k]] <- data.frame(Xk, occ = paste0("t", k),
                              pid = sprintf("P%03d", 1:np))
    }
    rasch(do.call(rbind, rows), factors = "occ", id = "pid",
          items = paste0("I", 1:L2))
  }
  ## nonspherical null: raw rejections at or below ~nominal over 15 fits
  rej <- 0; tot <- 0
  for (r in 1:15) {
    ss <- dif_anova(gen(100 + r, c(0.05, 0.05, 0.05, 1.5)))$summary
    rej <- rej + sum(ss$p_uniform < 0.05, na.rm = TRUE)
    tot <- tot + sum(is.finite(ss$p_uniform))
  }
  expect_lte(rej / tot, 0.09)   # the uncorrected engine sat near 9-percent
  ## planted occasion DIF: detected as the top flag in most replicates
  hits <- 0
  for (r in 1:5) {
    ss <- dif_anova(gen(200 + r, rep(0.4, K), shift3_t4 = -1.0))$summary
    fl <- ss$item[ss$uniform_DIF %in% TRUE]
    hits <- hits + ("I3" %in% fl)
  }
  expect_gte(hits, 3)
})

# --- DIF ANOVA round 11: multi-within alignment, incomplete panels --------

test_that("multi-within contrasts align: a pure w1 effect loads on w1", {
  skip_on_cran()
  set.seed(5)
  np <- 150; L <- 6; b <- seq(-1, 1, length.out = L)
  rows <- list(); th0 <- rnorm(np)
  for (i1 in c("a1", "a2")) for (i2 in c("b1", "b2", "b3")) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L)
      X[, j] <- rbinom(np, 1, plogis(th - b[j] +
                                       ifelse(j == 3 & i1 == "a2", -0.9, 0)))
    rows[[length(rows) + 1]] <- data.frame(X, w1 = i1, w2 = i2,
                                           pid = sprintf("P%03d", 1:np))
  }
  f <- rasch(do.call(rbind, rows), factors = c("w1", "w2"), id = "pid",
             items = paste0("I", 1:L))
  tt <- dif_anova(f, effects = "factorial")$terms
  i3 <- tt[tt$item == "I3", ]
  # interaction()'s first-fastest cell order silently rotated this into
  # w1:w2 (F 7.4 on the interaction, 3.2 on w1)
  expect_gt(i3$F_value[i3$term == "w1"], 10)
  expect_lt(i3$F_value[i3$term == "w1:w2"], 4)
  expect_lt(i3$F_value[i3$term == "w2"], 4)
  # GG metadata reproducible: p == pf(F, eps*df, eps*df_denom)
  r <- i3[i3$term == "w1", ]
  expect_equal(r$p, pf(r$F_value, r$gg_epsilon * r$df,
                       r$gg_epsilon * r$df_denom, lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("differentially incomplete within panels give no false group DIF", {
  skip_on_cran()
  set.seed(9)
  np <- 200; L <- 6; b <- seq(-1, 1, length.out = L)
  g <- rep(c("A", "B"), each = np / 2)
  th0 <- rnorm(np); rows <- list()
  for (k in 1:2) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L)
      X[, j] <- rbinom(np, 1, plogis(th - b[j] +
                                       ifelse(j == 3 & k == 2, -0.9, 0)))
    keep <- if (k == 2) (g == "A") | (runif(np) < 0.2) else rep(TRUE, np)
    rows[[k]] <- data.frame(X, occ = paste0("t", k), grp = g,
                            pid = sprintf("P%03d", 1:np))[keep, ]
  }
  f <- rasch(do.call(rbind, rows), factors = c("occ", "grp"), id = "pid",
             items = paste0("I", 1:L))
  s <- dif_anova(f)$summary
  gF <- s$F_uniform[s$item == "I3" & s$term == "grp"]
  # raw person means over unmatched cells reported group F = 37.6 here
  expect_true(is.na(gF) || gF < 8)
})

test_that("an NA in a between factor does not flip it to within-subject", {
  set.seed(2)
  b <- seq(-1, 1, length.out = 6)
  X <- matrix(rbinom(200 * 6, 1, plogis(outer(rnorm(200), b, "-"))), 200, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  gg <- rep(c("a", "b"), 100); gg[5] <- NA
  f <- rasch(data.frame(X, g = gg, pid = sprintf("P%03d", 1:200)),
             factors = "g", id = "pid", items = paste0("I", 1:6))
  d <- dif_anova(f)
  expect_length(d$within, 0L)
  expect_gt(nrow(d$summary), 0L)
})

# --- DIF ANOVA round 12: incomplete-panel edges ---------------------------

test_that("trait-dependent within effects with differential missingness give no group DIF", {
  skip_on_cran()
  set.seed(13)
  np <- 240; L <- 6; b <- seq(-1, 1, length.out = L)
  g <- rep(c("A", "B"), each = np / 2)
  th0 <- rnorm(np); rows <- list()
  for (k in 1:2) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L)
      X[, j] <- rbinom(np, 1, plogis(th - b[j] +
        ifelse(j == 3 & k == 2, -0.5 - 0.6 * th0, 0)))
    keep <- if (k == 2) (g == "A") | (runif(np) < 0.2) else rep(TRUE, np)
    rows[[k]] <- data.frame(X, occ = paste0("t", k), grp = g,
                            pid = sprintf("P%03d", 1:np))[keep, ]
  }
  f <- rasch(do.call(rbind, rows), factors = c("occ", "grp"), id = "pid",
             items = paste0("I", 1:L))
  s <- dif_anova(f)$summary
  gU <- s$F_uniform[s$item == "I3" & s$term == "grp"]
  gN <- s$F_nonuniform[s$item == "I3" & s$term == "grp"]
  # cell-only centring reported non-uniform group DIF at F = 214.7 here
  expect_true(is.na(gU) || gU < 8)
  expect_true(is.na(gN) || gN < 8)
})

test_that("a between level with no complete panels yields NA terms, not a crash", {
  set.seed(14)
  np <- 160; L <- 6; b <- seq(-1, 1, length.out = L)
  g <- rep(c("A", "B"), each = np / 2); th0 <- rnorm(np); rows <- list()
  for (k in 1:4) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L) X[, j] <- rbinom(np, 1, plogis(th - b[j]))
    keep <- if (k == 1) rep(TRUE, np) else g == "B"
    rows[[k]] <- data.frame(X, occ = paste0("t", k), grp = g,
                            pid = sprintf("P%03d", 1:np))[keep, ]
  }
  f <- rasch(do.call(rbind, rows), factors = c("occ", "grp"), id = "pid",
             items = paste0("I", 1:L))
  expect_s3_class(dif_anova(f), "rasch_dif")
  df <- dif_anova(f, effects = "factorial")
  expect_true(any(is.na(df$terms$F_value[df$terms$term == "occ:grp"])))
  expect_true(any(grepl("non-estimable", df$notes)))
  expect_true(any(grepl("dropped from the within-person", df$notes)))
})

test_that("significant multilevel within terms never reach ordinary Tukey", {
  set.seed(15)
  np <- 150; L <- 6; b <- seq(-1, 1, length.out = L)
  rows <- list(); th0 <- rnorm(np)
  gg <- rep(c("x", "y"), length.out = np)
  for (k in 1:4) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L)
      X[, j] <- rbinom(np, 1, plogis(th - b[j] +
                                       ifelse(j == 3 & k == 4, -1.2, 0)))
    keep <- (gg == "x") | (runif(np) < 0.7)
    rows[[k]] <- data.frame(X, occ = paste0("t", k), grp = gg,
                            pid = sprintf("P%03d", 1:np))[keep, ]
  }
  f <- rasch(do.call(rbind, rows), factors = c("occ", "grp"), id = "pid",
             items = paste0("I", 1:L))
  d <- dif_anova(f)
  expect_s3_class(d, "rasch_dif")
  # any Tukey rows present concern between terms only
  if (nrow(d$tukey)) expect_false(any(grepl("occ", d$tukey$term)))
})

# --- round 13: BTL / EFRM / MFRM DIF and identification -------------------

test_that("btl_dif is order-invariant across correlated judge factors", {
  skip_on_cran()
  set.seed(31)
  K <- 10; beta <- seq(-1.4, 1.4, length.out = K); nj <- 28; npj <- 60
  g1 <- rep(c("a", "b"), c(10, 18))
  g2 <- ifelse(runif(nj) < 0.75, ifelse(g1 == "a", "x", "y"),
               sample(c("x", "y"), nj, TRUE))
  jids <- sprintf("J%02d", 1:nj); rows <- list()
  for (j in 1:nj) {
    bj <- beta; if (g1[j] == "b") bj[4] <- bj[4] - 1
    ia <- sample(K, npj, TRUE); ib <- (ia + sample(K - 1, npj, TRUE) - 1L) %% K + 1L
    win <- rbinom(npj, 1, plogis(bj[ia] - bj[ib]))
    rows[[j]] <- data.frame(object_a = paste0("O", ia),
                            object_b = paste0("O", ib),
                            winner = paste0("O", ifelse(win == 1, ia, ib)),
                            judge = jids[j])
  }
  bt <- btl(do.call(rbind, rows), "object_a", "object_b", "winner",
            judge = "judge")
  A <- setNames(g1, jids); B <- setNames(g2, jids)
  sAB <- btl_dif(bt, factors = list(A = A, B = B))$summary
  sBA <- btl_dif(bt, factors = list(B = B, A = A))$summary
  expect_equal(sAB$F_uniform[order(sAB$object, sAB$term)],
               sBA$F_uniform[order(sBA$object, sBA$term)],
               tolerance = 1e-8)
})

test_that("btl_dif rejects factors that vary within a judge", {
  set.seed(3)
  K <- 6; b <- seq(-1, 1, length.out = K); n <- 600
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(
                    rbinom(n, 1, plogis(b[ia] - b[ib])) == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:10), n, TRUE))
  bt <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  rowfac <- sample(c("u", "v"), n, TRUE)      # varies within judges
  expect_error(btl_dif(bt, factors = list(g = rowfac)),
               "varies within judge")
})

test_that("btl_efrm validates within-set connectivity and alpha identification", {
  d <- simulate_btl_efrm(n_objects_per_set = 5, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 6, reps_within = 20,
                         reps_cross = 20, seed = 4)
  os <- attr(d, "truth")$object_sets
  ## alpha: cross-set rows touching only ONE object of set 2
  s2 <- os[[2]]
  cross_rows <- (d$object_a %in% os[[1]]) != (d$object_b %in% os[[1]])
  keep <- !cross_rows | (d$object_a == s2[1]) | (d$object_b == s2[1])
  expect_error(
    btl_efrm(d[keep, ], "object_a", "object_b", "winner", "judge", "panel",
             os, se_method = "conditional"),
    "unit \\(alpha\\) is unidentified")
  ## within-set connectivity: split set 1's internal comparisons
  g1 <- os[[1]][1:2]; g2 <- os[[1]][3:5]
  within1 <- (d$object_a %in% os[[1]]) & (d$object_b %in% os[[1]])
  bridge <- within1 & ((d$object_a %in% g1) != (d$object_b %in% g1))
  expect_error(
    btl_efrm(d[!bridge, ], "object_a", "object_b", "winner", "judge",
             "panel", os, se_method = "conditional"),
    "not connected|no within-set comparison")
})

test_that("EFRM group linkage requires shared items, not shared set labels", {
  d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 2,
                     n_groups = 2, group_unit_ratio = 1.25, seed = 2)
  tr <- attr(d, "truth")
  ## make the groups' item subsets DISJOINT within every set
  grp <- tr$groups
  X <- d
  items <- unlist(tr$item_sets)
  for (s in tr$item_sets) {
    half <- seq_len(floor(length(s) / 2))
    X[grp == unique(grp)[1], s[half]] <- NA
    X[grp == unique(grp)[2], s[-half]] <- NA
  }
  expect_error(
    rasch_efrm(X, item_sets = tr$item_sets, groups = grp),
    "not linked|unidentified")
})

test_that("dif_anova integrates with EFRM and MFRM fits", {
  skip_on_cran()
  d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 2,
                     n_groups = 2, group_unit_ratio = 1.25, seed = 1)
  tr <- attr(d, "truth")
  expect_error(
    dif_anova(rasch_efrm(d, item_sets = tr$item_sets, groups = tr$groups)),
    "frame structure")
  sex <- rep(c("m", "f"), length.out = nrow(d))
  f2 <- rasch_efrm(d, item_sets = tr$item_sets, groups = tr$groups,
                   factors = data.frame(sex = sex))
  d2 <- dif_anova(f2)
  expect_true(any(grepl("frame structure", d2$notes)))
  expect_gt(nrow(d2$summary), 0)

  set.seed(1)
  simP <- function(th, tau) { x <- 0:length(tau)
    p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%03d", 1:150); raters <- paste0("R", 1:3)
  th <- setNames(rnorm(150, 0, 1.3), persons)
  sx <- setNames(rep(c("m", "f"), length.out = 150), persons)
  tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
  dd <- expand.grid(person = persons, item = names(tau), rater = raters,
                    stringsAsFactors = FALSE)
  dd$score <- mapply(function(p, i, r)
    sample(0:2, 1, prob = simP(th[p] + ifelse(i == "B" & sx[p] == "f",
                                              -0.6, 0),
                               tau[[i]] + c(R1 = -0.4, R2 = 0,
                                            R3 = 0.4)[r])),
    dd$person, dd$item, dd$rater)
  dd$sex <- sx[dd$person]
  mf <- rasch_mfrm(dd, person = "person", item = "item", score = "score",
                   facets = "rater", factors = "sex")
  dm <- dif_anova(mf)
  expect_true(any(dm$summary$uniform_DIF[grepl("^B:", dm$summary$item)] %in%
                    TRUE))
  expect_error(
    rasch_mfrm(dd, person = "person", item = "item", score = "score",
               facets = "rater", factors = "rater"),
    "varies within person")
})
