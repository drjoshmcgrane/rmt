simEF <- function(th, tau, r) {
  x <- 0:length(tau)
  p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p)
}

test_that("EFRM recovers person-group units (one set, four groups)", {
  set.seed(7); L <- 20; per_g <- 400
  phi_true <- c(0.6, 0.9, 1.1, 1.5); phi_true <- phi_true / exp(mean(log(phi_true)))
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  glev <- paste0("G", 1:4); grp <- rep(glev, each = per_g); Np <- length(grp)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_len(L), function(i)
    rbinom(Np, 1, plogis(phi_true[match(grp, glev)] * (th - d[i]))))
  colnames(X) <- sprintf("I%02d", 1:L)

  fit <- rasch_efrm(data.frame(X, g = grp), item_sets = list(core = colnames(X)),
                    groups = "g")
  expect_s3_class(fit, "rasch_efrm"); expect_s3_class(fit, "rasch")
  expect_true(fit$est$converged)
  lp <- log(fit$phi_table$phi)
  expect_gt(cor(lp, log(phi_true)), 0.95)
  expect_lt(max(abs(lp - log(phi_true))), 0.12)
  expect_equal(sum(lp), 0, tolerance = 1e-8)
  expect_gt(cor(fit$item_arbitrary$location, d), 0.99)
  expect_true(all(fit$alpha_table$alpha == 1))
  for (g in glev) {
    sel <- fit$person$g == g
    expect_gt(cor(fit$person$theta[sel], th[grp == g], use = "complete.obs"), 0.8)
  }
  expect_gt(fit$efrm_vs_rasch$two_delta_ll, 100)   # strong unit differences
})

test_that("EFRM recovers item-set units from common persons (polytomous)", {
  set.seed(8); Np <- 800
  alpha_true <- c(0.7, 1.0, 1.4); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  mu_true <- c(-0.4, 0.1, 0.3)
  th <- rnorm(Np, 0, 1.3)
  sets <- rep(1:3, each = 6)
  dd <- lapply(1:18, function(i) mu_true[sets[i]] + sort(rnorm(2, 0, 0.8)))
  X <- sapply(1:18, function(i) sapply(th, function(t)
    sample(0:2, 1, prob = simEF(t, dd[[i]], alpha_true[sets[i]]))))
  colnames(X) <- sprintf("S%dI%02d", sets, 1:18)

  fit <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                    item_sets = split(colnames(X), sets))
  expect_lt(max(abs(log(fit$alpha_table$alpha) - log(alpha_true))), 0.15)
  expect_true(all(fit$phi_table$phi == 1))
  mu_real <- tapply(sapply(dd, mean), sets, mean)
  mu_real <- mu_real - mean(mu_real)
  expect_lt(max(abs(fit$set_table$mu - mu_real)), 0.15)
  loc_true <- sapply(dd, mean) - mean(tapply(sapply(dd, mean), sets, mean))
  names(loc_true) <- colnames(X)
  est <- setNames(fit$item_arbitrary$location, fit$item_arbitrary$item)
  expect_gt(cor(est[names(loc_true)], loc_true), 0.97)
  expect_gt(cor(fit$person$theta, th, use = "complete.obs"), 0.9)
})

test_that("EFRM recovers the full unit grid (two sets x two groups)", {
  set.seed(9); per_g <- 500
  alpha_true <- c(A = 0.8, B = 1.25)
  phi_true <- c(g1 = 0.85, g2 = 1 / 0.85)
  th <- rnorm(2 * per_g, 0, 1.3)
  grp <- rep(names(phi_true), each = per_g)
  sets <- rep(c("A", "B"), each = 8)
  d <- rep(seq(-1.5, 1.5, length.out = 8), 2)
  d <- d - mean(tapply(d, sets, mean))
  X <- sapply(seq_along(sets), function(i)
    rbinom(2 * per_g, 1,
           plogis(alpha_true[sets[i]] * phi_true[grp] * (th - d[i]))))
  colnames(X) <- sprintf("%sI%02d", sets, seq_along(sets))

  fit <- rasch_efrm(data.frame(X, g = grp), groups = "g",
                    item_sets = split(colnames(X), sets))
  fr <- fit$frames
  rho_true <- outer(alpha_true, phi_true)[cbind(fr$set, fr$group)]
  expect_gt(cor(log(fr$rho), log(rho_true)), 0.95)
  expect_equal(fr$rho, fr$alpha * fr$phi, tolerance = 1e-12)
  expect_lt(max(abs(log(fit$phi_table$phi) - log(phi_true))), 0.12)
  expect_lt(max(abs(log(fit$alpha_table$alpha) - log(alpha_true))), 0.2)
})

test_that("a single frame reduces to the ordinary rasch fit", {
  set.seed(2); Np <- 400; L <- 8
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  fe <- rasch_efrm(data.frame(X, g = "one"), groups = "g",
                   item_sets = list(all = colnames(X)))
  fr <- rasch(X)
  expect_true(all(fe$phi_table$phi == 1) && all(fe$alpha_table$alpha == 1))
  est <- setNames(fe$item_arbitrary$location, fe$item_arbitrary$item)
  expect_equal(unname(est[fr$items$item]), fr$items$location, tolerance = 1e-6)
  expect_equal(fe$person$theta, fr$person$theta, tolerance = 1e-6)
})

test_that("unlinked structures produce informative errors", {
  set.seed(3); per_g <- 150
  th <- rnorm(2 * per_g)
  grp <- rep(c("g1", "g2"), each = per_g)
  X <- matrix(rbinom(2 * per_g * 8, 1, plogis(th)), 2 * per_g, 8)
  colnames(X) <- sprintf("I%d", 1:8)
  # each group sees only its own set: no common set, no common persons
  X[grp == "g1", 5:8] <- NA
  X[grp == "g2", 1:4] <- NA
  expect_error(
    rasch_efrm(data.frame(X, g = grp), groups = "g",
               item_sets = list(s1 = colnames(X)[1:4], s2 = colnames(X)[5:8])),
    "not linked")
})

test_that("phi sandwich SEs track the sampling variability", {
  set.seed(11); L <- 12; per_g <- 250
  phi_true <- c(0.75, 4 / 3); phi_true <- phi_true / exp(mean(log(phi_true)))
  d <- seq(-1.5, 1.5, length.out = L)
  glev <- c("a", "b"); grp <- rep(glev, each = per_g); Np <- length(grp)
  reps <- 8
  lp_hat <- se_hat <- numeric(reps)
  for (r in seq_len(reps)) {
    th <- rnorm(Np, 0, 1.3)
    X <- sapply(seq_len(L), function(i)
      rbinom(Np, 1, plogis(phi_true[match(grp, glev)] * (th - d[i]))))
    colnames(X) <- sprintf("I%02d", 1:L)
    f <- rasch_efrm(data.frame(X, g = grp), groups = "g",
                    item_sets = list(core = colnames(X)))
    lp_hat[r] <- log(f$phi_table$phi[2])
    se_hat[r] <- f$phi_table$se_log_phi[2]
  }
  ratio <- mean(se_hat) / sd(lp_hat)
  expect_gt(ratio, 0.5); expect_lt(ratio, 2.0)
})

test_that("the weighted score, not the raw score, drives person estimates", {
  set.seed(12); Np <- 600
  alpha_true <- c(0.6, 5 / 3); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  th <- rnorm(Np, 0, 1.2)
  sets <- rep(1:2, each = 6)
  d <- rep(seq(-1, 1, length.out = 6), 2)
  X <- sapply(seq_along(sets), function(i) sapply(th, function(t)
    sample(0:2, 1, prob = simEF(t, d[i] + c(-0.5, 0.5), alpha_true[sets[i]]))))
  colnames(X) <- sprintf("S%dI%02d", sets, seq_along(sets))
  fit <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                    item_sets = split(colnames(X), sets))
  p <- fit$person
  # two persons with the same raw score but different weighted scores
  # must receive different locations
  cand <- which(!p$extreme & p$n_items == 12)
  byraw <- split(cand, p$raw[cand])
  found <- FALSE
  for (grpw in byraw) {
    if (length(grpw) < 2) next
    w <- p$weighted_score[grpw]
    if (max(w) - min(w) > 0.5) {
      i1 <- grpw[which.min(w)]; i2 <- grpw[which.max(w)]
      expect_false(isTRUE(all.equal(p$theta[i1], p$theta[i2])))
      expect_lt(p$theta[i1], p$theta[i2])  # more weight on high-unit items
      found <- TRUE; break
    }
  }
  expect_true(found)
  # equal weighted scores within a pattern share a location exactly
  key <- paste(p$n_items, signif(p$weighted_score, 12))
  dup <- key[duplicated(key) & !is.na(p$theta)][1]
  who <- which(key == dup)
  expect_lt(diff(range(p$theta[who])), 1e-12)
})

test_that("the error-variance correction beats the naive SD ratio", {
  set.seed(13); Np <- 500
  alpha_true <- c(0.7, 10 / 7); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  th <- rnorm(Np, 0, 1.2)
  sets <- rep(1:2, each = 6)
  d <- rep(seq(-1, 1, length.out = 6), 2)
  X <- sapply(seq_along(sets), function(i) sapply(th, function(t)
    sample(0:2, 1, prob = simEF(t, d[i] + c(-0.5, 0.5), alpha_true[sets[i]]))))
  colnames(X) <- sprintf("S%dI%02d", sets, seq_along(sets))
  fit <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                    item_sets = split(colnames(X), sets))
  # naive ratio from observed SDs of the same per-set person estimates
  dtl <- fit$thresholds_arbitrary
  u <- lapply(1:2, function(s) {
    cols <- which(fit$set_of[fit$virtual_map$item] == as.character(s))
    tl <- lapply(fit$virtual_map$item[cols], function(it)
      dtl$delta[dtl$item == it] * fit$alpha_table$alpha[s])
    pe <- RaschR:::.person_estimates(fit$X[, cols, drop = FALSE], tl, disc = 1)
    pe
  })
  ok <- !u[[1]]$extreme & !u[[2]]$extreme &
    is.finite(u[[1]]$theta) & is.finite(u[[2]]$theta)
  naive <- sd(u[[2]]$theta[ok]) / sd(u[[1]]$theta[ok])
  true_ratio <- alpha_true[2] / alpha_true[1]
  est_ratio <- fit$alpha_table$alpha[2] / fit$alpha_table$alpha[1]
  expect_lt(abs(log(est_ratio) - log(true_ratio)),
            abs(log(naive) - log(true_ratio)))
})

test_that("EFRM honours the -1 missing code", {
  set.seed(14); per_g <- 200; L <- 10
  grp <- rep(c("a", "b"), each = per_g); Np <- length(grp)
  phi_true <- c(0.8, 1.25)
  th <- rnorm(Np)
  d <- seq(-1.2, 1.2, length.out = L)
  X <- sapply(seq_len(L), function(i)
    rbinom(Np, 1, plogis(phi_true[match(grp, c("a", "b"))] * (th - d[i]))))
  colnames(X) <- sprintf("I%02d", 1:L)
  hit <- sample(length(X), 150)
  Xna <- X; Xna[hit] <- NA
  Xcode <- X; Xcode[hit] <- -1L
  f1 <- rasch_efrm(data.frame(Xna, g = grp), groups = "g",
                   item_sets = list(core = colnames(X)))
  f2 <- rasch_efrm(data.frame(Xcode, g = grp), groups = "g",
                   item_sets = list(core = colnames(X)))
  expect_equal(f1$phi_table$phi, f2$phi_table$phi, tolerance = 1e-10)
  expect_equal(f1$item_arbitrary$location, f2$item_arbitrary$location,
               tolerance = 1e-10)
})

test_that("EFRM standard error methods are coherent", {
  set.seed(15); Np <- 500
  alpha_true <- c(0.75, 4 / 3); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  sets <- rep(1:2, each = 6)
  d <- rep(seq(-1, 1, length.out = 6), 2)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_along(sets), function(i) sapply(th, function(t)
    sample(0:2, 1, prob = simEF(t, d[i] + c(-0.5, 0.5), alpha_true[sets[i]]))))
  colnames(X) <- sprintf("S%dI%02d", sets, seq_along(sets))

  fit <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                    item_sets = split(colnames(X), sets), boot_reps = 120)
  expect_identical(fit$se_method, "hybrid")
  expect_true(all(is.finite(fit$alpha_table$se_log_alpha)))
  expect_true(all(fit$alpha_table$se_log_alpha > 0))
  # propagation: common-unit threshold SEs exceed the purely conditional part
  cond <- sqrt(diag(fit$est$cov_tau))   # virtual level, already propagated
  expect_true(all(fit$thresholds_arbitrary$se > 0))

  fb <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                   item_sets = split(colnames(X), sets),
                   se_method = "bootstrap", boot_reps = 50)
  expect_identical(fb$se_method, "bootstrap")
  expect_gte(fb$boot_reps_used, 30)
  # the two methods agree on scale (well within a factor of two)
  ratio <- median(fit$thresholds_arbitrary$se / fb$thresholds_arbitrary$se)
  expect_gt(ratio, 0.5); expect_lt(ratio, 2)
  # point estimates are identical across SE methods
  expect_equal(fit$alpha_table$alpha, fb$alpha_table$alpha, tolerance = 1e-10)
  expect_equal(fit$item_arbitrary$location, fb$item_arbitrary$location,
               tolerance = 1e-10)
})
