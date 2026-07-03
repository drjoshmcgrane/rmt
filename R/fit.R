# rmt :: fit statistics
# ===========================================================================
# Test-of-fit statistics for the pairwise analysis: standardised residuals;
# the log-of-mean-square fit residual of Andrich & Marais (2019, ch. 23)
# for items and persons (with its
# untransformed "natural" form); infit and outfit mean squares with
# Wilson-Hilferty standardisations; the class-interval ANOVA item-fit F;
# the item-trait interaction chi-square; the person separation index with
# and without extremes; Cronbach's alpha; targeting; and the test
# information function.
# ===========================================================================

.wh <- function(ms, q) (ms^(1/3) - 1) * (3 / q) + (q / 3)   # mean square -> z

# Because each person's location is estimated from their own responses, the
# expected squared standardised residual in cell (n, i) is close to
# 1 - V_ni / sum_j V_nj rather than 1; mean squares are rescaled by this
# information share so they centre on 1 under fit.
.z2_expectation <- function(mo, Z, disc = NULL) {
  Vobs <- mo$V; Vobs[is.na(Z)] <- NA
  if (!is.null(disc)) Vobs <- sweep(Vobs, 2, disc^2, "*")
  share <- Vobs / rowSums(Vobs, na.rm = TRUE)
  pmax(1 - share, 1e-4)
}

# ---------------------------------------------------------------------------
# The fit residual of Andrich & Marais (2019, ch. 23; see also Andrich 1988).
# The observed cells of non-extreme persons
# carry C - (N + P) model-testing degrees of freedom, where P item parameters
# and N person locations were estimated (with complete dichotomous data this
# is (N-1)(I-1) - (m-1)); apportioned equally, each cell carries
# f_cell = df_total / C. Summing z^2 over an item's (person's) observed
# cells gives Y^2 with E[Y^2] = f (the summed cell df) and model variance
# V[Y^2] = sum(C4/V^2 - 1) (the dichotomous form 2 tanh((b-d)/2) sinh(b-d)
# generalised to ordered categories). The reported fit residual is the
# symmetrising log-of-mean-square transform
#   T2 = f (ln Y^2 - ln f) / sqrt(V[Y^2])    (A&M 2019, eq. 23.14)
# with expectation 0 and variance 1 under fit: negative = too deterministic /
# over-discriminating (Guttman-like), positive = erratic /
# under-discriminating. The conventional flagging value is |T2| > 2.5
# (Andrich & Marais 2019, ch. 15). The untransformed statistic
#   T1 = (Y^2 - f) / sqrt(V[Y^2])
# is the "natural" fit residual and is kept alongside as natural_resid.
# Extreme persons are excluded throughout, exactly as they are set aside
# from the calibrating sample, so their person fit residual is NA and they
# contribute nothing to item fit or to the degrees of freedom.
# ---------------------------------------------------------------------------
.fitres_df <- function(Z, extreme, n_parameters) {
  obs <- !is.na(Z) & !extreme
  C <- sum(obs)
  N_ne <- sum(rowSums(obs) > 0)
  df_total <- C - N_ne - n_parameters
  f_cell <- if (C > 0 && df_total > 0) df_total / C else NA_real_
  list(obs = obs, f_cell = f_cell,
       f_item = f_cell * colSums(obs), f_person = f_cell * rowSums(obs))
}

# One margin's fit residuals from cell sums: y2 = sum z^2, v = sum of
# per-cell variances C4/V^2 - 1, f = summed cell df, n = cells used.
.fitres_transform <- function(y2, v, f, n, min_cells = 2L) {
  ok <- n >= min_cells & !is.na(f) & f > 0 & v > 1e-8 & y2 > 0
  natural <- fit_resid <- rep(NA_real_, length(y2))
  natural[ok] <- (y2[ok] - f[ok]) / sqrt(v[ok])
  fit_resid[ok] <- f[ok] * (log(y2[ok]) - log(f[ok])) / sqrt(v[ok])
  list(fit_resid = fit_resid, natural = natural, df = ifelse(ok, f, NA_real_))
}

.fitres <- function(Z, mo, extreme, n_parameters) {
  dfs <- .fitres_df(Z, extreme, n_parameters)
  z2 <- Z^2; z2[!dfs$obs] <- NA
  vcell <- mo$M4 / mo$V^2 - 1; vcell[!dfs$obs] <- NA
  it <- .fitres_transform(colSums(z2, na.rm = TRUE),
                        colSums(vcell, na.rm = TRUE),
                        dfs$f_item, colSums(dfs$obs))
  pe <- .fitres_transform(rowSums(z2, na.rm = TRUE),
                        rowSums(vcell, na.rm = TRUE),
                        dfs$f_person, rowSums(dfs$obs), min_cells = 3L)
  list(items = it, persons = pe, f_cell = dfs$f_cell)
}

# Class-interval ANOVA item fit (Andrich & Marais 2019, ch. 15): per item,
# a one-way analysis of variance of the standardised residuals
# over the class intervals. Under fit the interval means share a common
# zero mean, so the between-interval F on (G - 1, n - G) degrees of freedom
# tests the same item-trait interaction as the chi-square but through the
# ANOVA calibration. Reported with Benjamini-Hochberg (false discovery
# rate) and Bonferroni (familywise) adjustments across items.
.item_anova <- function(Z, ci, extreme, ci_list = NULL) {
  L <- ncol(Z)
  out <- data.frame(item = colnames(Z), F_anova = NA_real_, df1 = NA_integer_,
                    df2 = NA_integer_, p = NA_real_)
  for (i in seq_len(L)) {
    ci_i <- if (is.null(ci_list)) ci else ci_list[[i]]
    sel <- which(!is.na(Z[, i]) & !is.na(ci_i) & !extreme)
    if (!length(sel)) next
    g <- ci_i[sel]; z <- Z[sel, i]
    keep <- g %in% as.integer(names(which(table(g) >= 2)))
    g <- factor(g[keep]); z <- z[keep]
    G <- nlevels(g); n <- length(z)
    if (G < 2 || n - G < 1) next
    mg <- tapply(z, g, mean); ng <- tabulate(g)
    ssb <- sum(ng * (mg - mean(z))^2)
    ssw <- sum((z - mg[g])^2)
    if (ssw <= 0) next
    out$F_anova[i] <- (ssb / (G - 1)) / (ssw / (n - G))
    out$df1[i] <- G - 1L; out$df2[i] <- n - G
    out$p[i] <- pf(out$F_anova[i], G - 1, n - G, lower.tail = FALSE)
  }
  out$p_adj <- p.adjust(out$p, method = "BH")
  out$p_bonf <- p.adjust(out$p, method = "bonferroni")
  out
}

# Item fit from per-person model moments (observed cells only).
.item_fit <- function(X, Z, mo, disc = NULL) {
  L <- ncol(X)
  E2 <- .z2_expectation(mo, Z, disc)
  out <- data.frame(item = colnames(X), infit_ms = NA_real_, outfit_ms = NA_real_,
                    infit_z = NA_real_, outfit_z = NA_real_, n = NA_integer_)
  for (i in seq_len(L)) {
    ok <- which(!is.na(Z[, i]))
    if (length(ok) < 3) next
    z2 <- Z[ok, i]^2
    V <- mo$V[ok, i]; C4 <- mo$M4[ok, i]; n <- length(ok)
    e2 <- E2[ok, i]
    outfit <- sum(z2) / sum(e2)
    infit  <- sum(z2 * V) / sum(e2 * V)
    qo <- sqrt(max(sum(C4 / V^2) / n^2 - 1 / n, 1e-8))
    qi <- sqrt(max(sum(C4 - V^2) / sum(V)^2, 1e-8))
    out$outfit_ms[i] <- outfit; out$infit_ms[i] <- infit
    out$outfit_z[i] <- .wh(outfit, qo); out$infit_z[i] <- .wh(infit, qi)
    out$n[i] <- n
  }
  out
}

# Person fit residuals: each person's standardised residuals across their
# observed items, summarised exactly as for items.
.person_fit <- function(X, Z, mo, disc = NULL) {
  N <- nrow(X)
  E2 <- .z2_expectation(mo, Z, disc)
  out <- data.frame(infit_ms = rep(NA_real_, N), outfit_ms = NA_real_,
                    outfit_z = NA_real_)
  for (n in seq_len(N)) {
    ok <- which(!is.na(Z[n, ]))
    if (length(ok) < 3) next
    z2 <- Z[n, ok]^2
    V <- mo$V[n, ok]; C4 <- mo$M4[n, ok]; k <- length(ok)
    e2 <- E2[n, ok]
    outfit <- sum(z2) / sum(e2)
    out$outfit_ms[n] <- outfit
    out$infit_ms[n] <- sum(z2 * V) / sum(e2 * V)
    qo <- sqrt(max(sum(C4 / V^2) / k^2 - 1 / k, 1e-8))
    out$outfit_z[n] <- .wh(outfit, qo)
  }
  out
}

# The default number of class intervals (Andrich & Marais 2019, ch. 15):
# as many intervals of at least 50 persons as the non-extreme sample
# allows, at most 10, at least 2.
.default_n_groups <- function(n_ne) max(2L, min(10L, n_ne %/% 50L))

# Allocate locations to n_groups contiguous intervals, as equal-sized as
# possible WITHOUT splitting ties: persons sharing a location are
# indistinguishable (equal raw scores give equal measures), so they belong
# to the same interval and interval sizes are generally unequal, as in the
# worked class-interval tables of Andrich & Marais (2019, ch. 13; sizes
# such as 13/20/16). A boundary falls where the cumulative count comes
# closest to each equal-share target.
.ci_allocate <- function(th, n_groups) {
  ut <- sort(unique(th))
  if (length(ut) <= n_groups) return(match(th, ut))
  cnt <- as.integer(table(factor(th, levels = ut)))
  cum <- cumsum(cnt)
  n <- length(th)
  b <- integer(n_groups - 1L)
  lo <- 1L
  for (gg in seq_len(n_groups - 1L)) {
    target <- n * gg / n_groups
    cand <- seq(lo, length(ut) - (n_groups - gg))
    b[gg] <- cand[which.min(abs(cum[cand] - target))]
    lo <- b[gg] + 1L
  }
  grp_of_ut <- findInterval(seq_along(ut), b + 1L) + 1L
  grp_of_ut[match(th, ut)]
}

# Class intervals over non-extreme person locations. n_groups = NULL
# applies the default rule above.
.class_intervals <- function(theta, extreme, n_groups = NULL) {
  g <- rep(NA_integer_, length(theta))
  use <- which(!is.na(theta) & !extreme)
  if (is.null(n_groups)) n_groups <- .default_n_groups(length(use))
  g[use] <- .ci_allocate(theta[use], n_groups)
  attr(g, "n_groups") <- max(g[use], na.rm = TRUE)
  g
}

# Class intervals compiled per item, the automatic adjustment with missing
# data (Andrich & Marais 2019, ch. 15): each item allocates the persons who
# answered it into intervals of its own, with the group count from the same
# rule applied to that item's responders. Returns a list of allocation vectors, one per item.
.class_intervals_by_item <- function(X, theta, extreme, n_groups = NULL) {
  lapply(seq_len(ncol(X)), function(i) {
    obs <- !is.na(X[, i]) & !is.na(theta) & !extreme
    g <- rep(NA_integer_, length(theta))
    if (!any(obs)) return(g)
    ng <- if (is.null(n_groups)) .default_n_groups(sum(obs)) else n_groups
    g[obs] <- .ci_allocate(theta[obs], ng)
    g
  })
}

# Item-trait interaction chi-square over class intervals, per item. ci is
# the common allocation vector; ci_list, when supplied, gives each item its
# own allocation (per-item basis). The degrees of freedom are per item:
# (number of class intervals contributing at least 2 responders to that
# item) - 1, so items with missing data are tested on the intervals they
# actually reach.
.item_trait <- function(X, Z, mo, ci, adjust_N = NA, ci_list = NULL) {
  L <- ncol(X)
  chi <- setNames(numeric(L), colnames(X))
  used <- integer(L)
  for (i in seq_len(L)) {
    ci_i <- if (is.null(ci_list)) ci else ci_list[[i]]
    G <- suppressWarnings(max(ci_i, na.rm = TRUE))
    if (!is.finite(G)) next
    for (gg in seq_len(G)) {
      sel <- which(ci_i == gg & !is.na(X[, i]))
      if (length(sel) < 2) next
      Obar <- mean(X[sel, i])
      Ebar <- mean(mo$E[sel, i]); Vbar <- mean(mo$V[sel, i])
      chi[i] <- chi[i] + length(sel) * (Obar - Ebar)^2 / Vbar
      used[i] <- used[i] + 1L
    }
  }
  df_i <- pmax(used - 1L, 1L)
  n_used <- sum(!is.na(ci))
  if (!is.na(adjust_N)) chi <- chi * (adjust_N / n_used)
  p <- pchisq(chi, df_i, lower.tail = FALSE)
  p_adj <- p.adjust(p, method = "BH")
  data.frame(item = colnames(X), chisq = chi, df = df_i, p = p,
             p_adj = p_adj, p_bonf = p.adjust(p, method = "bonferroni"),
             misfit = p_adj < 0.05)
}

# Correlation that degrades to NA (rather than erroring) when fewer than 3
# complete pairs are available.
.safe_cor <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

# Distribution summary of a fit-statistic column: mean, SD, skewness, and
# (excess) kurtosis (the summary block of Andrich & Marais 2019, app. C).
.dist_stats <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(list(mean = NA_real_, sd = NA_real_,
                                 skewness = NA_real_, kurtosis = NA_real_))
  m <- mean(x); s <- sd(x); d <- x - m
  list(mean = m, sd = s,
       skewness = mean(d^3) / (mean(d^2))^1.5,
       kurtosis = mean(d^4) / (mean(d^2))^2 - 3)
}

#' Class-interval detail for one item's chi-square test of fit
#'
#' The per-class-interval breakdown behind an item's item-trait chi-square,
#' as dissected in Andrich and Marais (2019, ch. 13): for every class
#' interval the size, the maximum and mean
#' person location, the standardised residual between observed and expected
#' interval means, its squared chi-square component, the observed and
#' expected means (OM, EV), the sample-size-free effect size
#' ES = (OM - EV)/sqrt(mean V), and per response category the observed
#' proportion (OBS.P), the mean model probability (EST.P), and the observed
#' conditional threshold proportion (OBS.T), the proportion scoring k among
#' those scoring k - 1 or k.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or index.
#' @return A list with \code{item}, \code{location}, the \code{intervals}
#'   data frame, the \code{categories} data frame, the whole-sample observed
#'   mean \code{ave}, and the item's total \code{chisq}, \code{df}, and
#'   \code{p}. Intervals with fewer than 2 responders are shown but carry no
#'   chi-square contribution (\code{used = FALSE}), matching the item-trait
#'   computation.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' chisq_detail(rasch(X), "I3")$intervals
#' @export
chisq_detail <- function(fit, item) {
  i <- .item_idx(fit, item)
  # per-item interval allocation when the fit carries one (missing data)
  ci <- if (!is.null(fit$ci_item)) fit$ci_item[[i]] else fit$person$class_interval
  th <- fit$person$theta
  x <- fit$X[, i]; E <- fit$moments$E[, i]; V <- fit$moments$V[, i]
  mi <- length(fit$tau_list[[i]])
  G <- max(ci, na.rm = TRUE)
  iv <- data.frame(interval = seq_len(G), n = 0L, theta_max = NA_real_,
                   theta_mean = NA_real_, obs_mean = NA_real_,
                   exp_value = NA_real_, residual = NA_real_,
                   chisq = NA_real_, es = NA_real_, used = FALSE)
  cats <- expand.grid(interval = seq_len(G), category = 0:mi)
  cats <- cats[order(cats$interval, cats$category), ]
  cats$obs_p <- cats$est_p <- cats$obs_t <- NA_real_
  disc_i <- if (is.null(fit$disc)) 1 else fit$disc[i]
  for (g in seq_len(G)) {
    sel <- which(!is.na(ci) & ci == g & !is.na(x))
    iv$n[g] <- length(sel)
    if (!length(sel)) next
    iv$theta_max[g] <- max(th[sel]); iv$theta_mean[g] <- mean(th[sel])
    OM <- mean(x[sel]); EV <- mean(E[sel]); Vbar <- mean(V[sel])
    iv$obs_mean[g] <- OM; iv$exp_value[g] <- EV
    iv$es[g] <- (OM - EV) / sqrt(Vbar)
    if (length(sel) >= 2) {
      iv$residual[g] <- sqrt(length(sel)) * (OM - EV) / sqrt(Vbar)
      iv$chisq[g] <- iv$residual[g]^2
      iv$used[g] <- TRUE
    }
    P <- vapply(th[sel], function(b)
      item_moments(b, fit$tau_list[[i]], disc = disc_i)$P, numeric(mi + 1))
    est_p <- rowMeans(P)
    obs_n <- as.integer(table(factor(x[sel], levels = 0:mi)))
    rows <- cats$interval == g
    cats$obs_p[rows] <- obs_n / length(sel)
    cats$est_p[rows] <- est_p
    for (k in seq_len(mi)) {
      pair <- obs_n[k] + obs_n[k + 1]
      cats$obs_t[rows][k + 1] <- if (pair > 0) obs_n[k + 1] / pair else NA_real_
    }
  }
  it_row <- fit$item_trait[i, ]
  list(item = fit$items$item[i], location = fit$items$location[i],
       intervals = iv, categories = cats, ave = mean(x, na.rm = TRUE),
       chisq = it_row$chisq, df = it_row$df, p = it_row$p)
}

# Person separation index (separation reliability; Andrich 1982).
.psi <- function(theta, se, keep = TRUE) {
  ok <- !is.na(theta) & !is.na(se) & keep
  if (sum(ok) < 3) return(list(PSI = NA_real_, separation = NA_real_,
                               var_theta = NA_real_, mean_error_var = NA_real_,
                               n = sum(ok)))
  vt <- var(theta[ok]); mse <- mean(se[ok]^2)
  psi <- max((vt - mse) / vt, 0)
  sep <- if (psi < 1) sqrt(psi / (1 - psi)) else Inf
  list(PSI = psi, separation = sep, var_theta = vt, mean_error_var = mse,
       n = sum(ok))
}

# Cronbach's alpha (Cronbach 1951) on complete cases, reported alongside
# the PSI. Alpha has no missing-data form, so the applicable flag carries
# that caveat (the complete-case value is still reported, with its n).
.alpha <- function(X) {
  Xc <- X[stats::complete.cases(X), , drop = FALSE]
  applicable <- nrow(Xc) == nrow(X)
  if (nrow(Xc) < 3 || ncol(Xc) < 2) return(list(alpha = NA_real_, n = nrow(Xc),
                                                applicable = applicable))
  L <- ncol(Xc); vi <- apply(Xc, 2, var); vt <- var(rowSums(Xc))
  list(alpha = L / (L - 1) * (1 - sum(vi) / vt), n = nrow(Xc),
       applicable = applicable)
}

# Qualitative power-of-test-of-fit assessment, driven by the PSI.
.fit_power <- function(psi) {
  if (is.na(psi)) "unknown"
  else if (psi >= 0.9) "excellent"
  else if (psi >= 0.8) "good"
  else if (psi >= 0.7) "reasonable"
  else if (psi >= 0.5) "low"
  else "too low"
}

# Targeting summary: how well item thresholds cover the person distribution.
.targeting <- function(person, thresholds) {
  ok <- !is.na(person$theta)
  th <- person$theta[ok]
  list(person_mean = mean(th), person_sd = sd(th),
       person_mean_noext = mean(person$theta[ok & !person$extreme]),
       item_mean = 0,
       threshold_range = range(thresholds$tau),
       prop_below = mean(th < min(thresholds$tau)),
       prop_above = mean(th > max(thresholds$tau)))
}

#' Test information function
#'
#' Fisher information of the whole test over a grid of person locations, with
#' the corresponding standard error of measurement.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param grid Logit grid over which to evaluate the information.
#' @return A data frame with \code{theta}, \code{info}, and \code{sem}.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' head(test_information(rasch(X)))
#' @export
test_information <- function(fit, grid = seq(-6, 6, by = 0.1)) {
  L <- length(fit$tau_list)
  disc <- if (is.null(fit$disc)) rep(1, L) else fit$disc
  info <- vapply(grid, function(th)
    sum(vapply(seq_len(L), function(i)
      disc[i]^2 * item_moments(th, fit$tau_list[[i]], disc = disc[i])$V, 0)), 0)
  data.frame(theta = grid, info = info, sem = 1 / sqrt(info))
}
