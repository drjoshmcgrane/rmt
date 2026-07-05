# rmt :: plots
# ===========================================================================
# The Rasch diagnostic plot suite in base graphics with a modern flat style:
# item characteristic curves (with group overlay for DIF), category
# probability curves, threshold probability curves, the person-item
# threshold distribution, the threshold map, the test characteristic curve,
# test information and measurement error, item and person fit maps, the
# residual-correlation heatmap, residual principal-component loadings, and
# category frequencies.
# ===========================================================================

.rr <- list(
  ink    = "#1e293b", grid = "#e2e8f0", soft = "#94a3b8",
  blue   = "#2563eb", red  = "#dc2626", amber = "#f59e0b",
  teal   = "#0f766e", purple = "#7c3aed",
  pal    = c("#2563eb", "#dc2626", "#f59e0b", "#0f766e", "#7c3aed",
             "#db2777", "#65a30d", "#475569")
)

# Modern flat canvas: light horizontal grid, open axes. A title is drawn
# only when `main` carries information (item, person, summary figures);
# plot-type names are left to the surrounding context.
.rr_canvas <- function(xlim, ylim, xlab, ylab, main = "", grid_y = TRUE,
                       grid_x = FALSE, yaxis = TRUE, right = 1.5) {
  has_main <- !is.null(main) && nzchar(main)
  op <- par(mar = c(4.2, 4.4, if (has_main) 3.2 else 1.6, right),
            mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15, cex.lab = 1.0)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       main = "", axes = FALSE)
  if (has_main) title(main = main, adj = 0, line = 1.4)
  if (grid_y) abline(h = pretty(ylim), col = .rr$grid, lwd = 0.8)
  if (grid_x) abline(v = pretty(xlim), col = .rr$grid, lwd = 0.8)
  axis(1, col = .rr$grid, col.ticks = .rr$soft)
  if (yaxis) axis(2, col = .rr$grid, col.ticks = .rr$soft)
  invisible(op)
}

.rr_legend <- function(pos, ...) legend(pos, ..., bty = "n", text.col = .rr$ink,
                                        cex = 0.85)

.item_idx <- function(fit, item) if (is.character(item)) match(item, fit$items$item) else item

# Discrimination (frame unit) of column i; 1 unless the fit carries units.
.disc_of <- function(fit, i) if (is.null(fit$disc)) 1 else fit$disc[i]

# ---------------------------------------------------------------------------
# Item characteristic curve, with optional group overlay (the graphical DIF
# display).
# ---------------------------------------------------------------------------
#' Plot an item characteristic curve
#'
#' Draws the model expected-score curve with observed class-interval means
#' overlaid. With \code{group} supplied, observed means are drawn separately
#' per group, the conventional graphical DIF display.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @param group Optional person grouping vector, or one or more names of
#'   factors nominated in the fit, for a DIF overlay; several names give
#'   the factor-combination cells (the factorial display).
#' @param n_groups Number of class intervals for the observed means; by
#'   default the fit's own count, or -- with a \code{group} overlay -- a
#'   count adapted to keep the smallest group's interval cells adequately
#'   filled.
#' @param grid Logit grid over which to draw the model curve.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- sprintf("I%02d", 1:6)
#' plot_icc(rasch(X), "I03")
#' @export
plot_icc <- function(fit, item, group = NULL, n_groups = NULL,
                     grid = seq(-5, 5, 0.05)) {
  i <- .item_idx(fit, item); tau_i <- fit$tau_list[[i]]; mmax <- length(tau_i)
  if (is.character(group) && length(group) < nrow(fit$X) &&
      !is.null(fit$factors) && all(group %in% names(fit$factors)))
    group <- if (length(group) == 1L) fit$factors[[group]] else
      interaction(fit$factors[group], sep = ":", drop = TRUE)
  if (is.null(n_groups))
    n_groups <- if (is.null(group)) fit$n_groups else
      .dif_n_groups(fit, group)
  Ecurve <- vapply(grid, function(th)
    item_moments(th, tau_i, disc = .disc_of(fit, i))$E, 0)
  th <- fit$person$theta; x <- fit$X[, i]; ok <- !is.na(th) & !is.na(x)
  op <- .rr_canvas(range(grid), c(0, mmax), "Person location (logits)",
                   "Expected score",
                   sprintf("%s  (location %.3f)", fit$items$item[i],
                           fit$items$location[i]))
  on.exit(par(op))
  lines(grid, Ecurve, lwd = 3, col = .rr$ink)
  ci <- cut(rank(th[ok], ties.method = "first"), n_groups, labels = FALSE)
  if (is.null(group)) {
    obsTh <- tapply(th[ok], ci, mean); obsX <- tapply(x[ok], ci, mean)
    points(obsTh, obsX, pch = 21, bg = .rr$blue, col = "white", cex = 1.5, lwd = 1.2)
    .rr_legend("topleft", c("Model", "Observed (class intervals)"),
               lwd = c(3, NA), pch = c(NA, 21), pt.bg = c(NA, .rr$blue),
               col = c(.rr$ink, "white"), pt.cex = 1.4)
  } else {
    g <- factor(group)[ok]
    levs <- levels(droplevels(g))
    for (li in seq_along(levs)) {
      sel <- g == levs[li]
      obsTh <- tapply(th[ok][sel], ci[sel], mean)
      obsX <- tapply(x[ok][sel], ci[sel], mean)
      colr <- .rr$pal[(li - 1L) %% length(.rr$pal) + 1L]
      lines(obsTh, obsX, col = colr, lwd = 1.4, lty = 3)
      points(obsTh, obsX, pch = 21, bg = colr, col = "white", cex = 1.5, lwd = 1.2)
    }
    .rr_legend("topleft", levs, lwd = 1.4, lty = 3, pch = 21,
               pt.bg = .rr$pal[seq_along(levs)],
               col = .rr$pal[seq_along(levs)], pt.cex = 1.3)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Category probability curves.
# ---------------------------------------------------------------------------
#' Plot category probability curves
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @param grid Logit grid over which to draw the curves.
#' @param observed Overlay the observed category proportions per class
#'   interval (Andrich and Marais 2019, ch. 20).
#' @param n_groups Class intervals for the observed points.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(1:4, function(i) sapply(th, function(t) sample(0:3, 1, prob = simP(t, c(-1, 0, 1)))))
#' colnames(X) <- sprintf("P%02d", 1:4)
#' plot_ccc(rasch(X), "P01", observed = TRUE)
#' @export
plot_ccc <- function(fit, item, grid = seq(-6, 6, 0.05), observed = FALSE,
                     n_groups = fit$n_groups) {
  i <- .item_idx(fit, item); tau_i <- fit$tau_list[[i]]; mmax <- length(tau_i)
  P <- vapply(grid, function(th)
    item_moments(th, tau_i, disc = .disc_of(fit, i))$P, numeric(mmax + 1))
  ordered <- all(diff(tau_i) > 0) || mmax == 1L
  op <- .rr_canvas(range(grid), c(0, 1), "Person location (logits)",
                   "Category probability",
                   fit$items$item[i])
  on.exit(par(op))
  abline(v = tau_i, lty = 3, col = .rr$soft)
  for (cat in 0:mmax)
    lines(grid, P[cat + 1, ], lwd = 2.6,
          col = .rr$pal[cat %% length(.rr$pal) + 1L])
  if (observed) {
    th <- fit$person$theta; x <- fit$X[, i]; ok <- !is.na(th) & !is.na(x)
    ci <- cut(rank(th[ok], ties.method = "first"), n_groups, labels = FALSE)
    obsTh <- tapply(th[ok], ci, mean)
    for (cat in 0:mmax) {
      obsP <- tapply(x[ok] == cat, ci, mean)
      points(obsTh, obsP, pch = 21, cex = 1.2, lwd = 1.1, col = "white",
             bg = .rr$pal[cat %% length(.rr$pal) + 1L])
    }
  }
  mtext(if (ordered) "thresholds ordered" else "THRESHOLDS DISORDERED",
        side = 3, line = 0.2, adj = 0, cex = 0.8,
        col = if (ordered) .rr$teal else .rr$red)
  .rr_legend("right", paste0("Category ", 0:mmax), lwd = 2.6,
             col = .rr$pal[(0:mmax) %% length(.rr$pal) + 1L])
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Threshold probability curves: conditional adjacent-category ogives, each
# crossing 0.5 at its threshold.
# ---------------------------------------------------------------------------
#' Plot threshold probability curves
#'
#' Conditional probability of success at each threshold,
#' \code{P(X = k | X = k - 1 or k)}, a logistic ogive crossing 0.5 at the
#' threshold location. Disordered thresholds are immediately visible as
#' out-of-sequence ogives. With \code{observed = TRUE} the observed
#' conditional proportions per class interval are overlaid, the direct
#' check on whether each threshold discriminates (and hence whether
#' collapsing categories could ever be justified; Andrich and Marais 2019,
#' ch. 22).
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @param grid Logit grid over which to draw the curves.
#' @param observed Overlay the observed conditional threshold proportions
#'   per class interval.
#' @param n_groups Class intervals for the observed points.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(1:4, function(i) sapply(th, function(t) sample(0:3, 1, prob = simP(t, c(-1, 0, 1)))))
#' colnames(X) <- sprintf("P%02d", 1:4)
#' plot_threshold_prob(rasch(X), "P01")
#' @export
plot_threshold_prob <- function(fit, item, grid = seq(-6, 6, 0.05),
                                observed = FALSE, n_groups = fit$n_groups) {
  i <- .item_idx(fit, item); tau_i <- fit$tau_list[[i]]
  op <- .rr_canvas(range(grid), c(0, 1), "Person location (logits)",
                   "Threshold probability",
                   fit$items$item[i])
  on.exit(par(op))
  abline(h = 0.5, lty = 2, col = .rr$soft)
  for (k in seq_along(tau_i)) {
    colr <- .rr$pal[(k - 1L) %% length(.rr$pal) + 1L]
    lines(grid, plogis(.disc_of(fit, i) * (grid - tau_i[k])), lwd = 2.6, col = colr)
    points(tau_i[k], 0.5, pch = 21, bg = colr, col = "white", cex = 1.4)
  }
  if (observed) {
    # observed conditional threshold proportions per class interval:
    # among persons responding k - 1 or k, the proportion responding k
    # (Andrich & Marais 2019, ch. 22 rescoring check)
    th <- fit$person$theta; x <- fit$X[, i]; ok <- !is.na(th) & !is.na(x)
    ci <- cut(rank(th[ok], ties.method = "first"), n_groups, labels = FALSE)
    for (k in seq_along(tau_i)) {
      colr <- .rr$pal[(k - 1L) %% length(.rr$pal) + 1L]
      inpair <- ok & !is.na(x) & (x == k - 1L | x == k)
      cip <- ci[inpair[ok]]
      if (!sum(inpair)) next
      obsTh <- tapply(th[inpair], cip, mean)
      obsT <- tapply(x[inpair] == k, cip, mean)
      points(obsTh, obsT, pch = 21, cex = 1.2, lwd = 1.1, col = "white",
             bg = colr)
    }
  }
  .rr_legend("topleft", sprintf("threshold %d (%.3f)", seq_along(tau_i), tau_i),
             lwd = 2.6, col = .rr$pal[seq_along(tau_i)])
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Person-item threshold distribution: mirrored histograms on one logit scale.
# ---------------------------------------------------------------------------
#' Plot the person-item threshold distribution
#'
#' The targeting display: the person location distribution above the axis and
#' the item threshold distribution mirrored below it, on a shared logit
#' scale.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param bins Number of histogram bins.
#' @param xlim Optional logit range for the shared scale; persons and
#'   thresholds outside it are omitted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_pimap(rasch(X))
#' @export
plot_pimap <- function(fit, bins = 35, xlim = NULL) {
  th <- fit$person$theta[!is.na(fit$person$theta)]
  tau <- fit$thresholds$tau
  rng <- if (is.null(xlim)) range(c(th, tau)) + c(-0.4, 0.4) else sort(xlim)
  th <- th[th >= rng[1] & th <= rng[2]]
  tau <- tau[tau >= rng[1] & tau <= rng[2]]
  brk <- seq(rng[1], rng[2], length.out = bins + 1)
  hp <- hist(th, breaks = brk, plot = FALSE)
  hi <- hist(tau, breaks = brk, plot = FALSE)
  pp <- hp$counts / sum(hp$counts); pi <- hi$counts / sum(hi$counts)
  ymax <- max(pp) * 1.15; ymin <- -max(pi) * 1.6
  op <- .rr_canvas(rng, c(ymin, ymax), "Location (logits)", "Proportion",
                   "", grid_y = FALSE,
                   yaxis = FALSE)
  on.exit(par(op))
  at <- pretty(c(0, max(c(pp, pi))))
  axis(2, at = c(-rev(at[-1]), at), labels = c(rev(at[-1]), at),
       col = .rr$grid, col.ticks = .rr$soft, cex.axis = 0.85)
  abline(h = 0, col = .rr$ink, lwd = 1)
  rect(brk[-length(brk)], 0, brk[-1], pp, col = .rr$blue, border = "white", lwd = 0.6)
  rect(brk[-length(brk)], -pi, brk[-1], 0, col = .rr$amber, border = "white", lwd = 0.6)
  segments(mean(th), 0, mean(th), ymax * 0.95, col = .rr$ink, lty = 2)
  text(mean(th), ymax * 0.98, sprintf("persons: mean %.2f, SD %.2f", mean(th), sd(th)),
       col = .rr$ink, cex = 0.8, adj = -0.02)
  text(mean(tau), ymin * 0.98, sprintf("thresholds: mean %.2f, SD %.2f", mean(tau), sd(tau)),
       col = .rr$ink, cex = 0.8, adj = -0.02)
  .rr_legend("topleft", c("Persons", "Item thresholds"),
             fill = c(.rr$blue, .rr$amber), border = NA)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Wright map: the conventional vertical person-item map (Wright & Stone
# 1979), person distribution left of a common logit axis, item threshold
# labels stacked to its right.
# ---------------------------------------------------------------------------
#' Plot a Wright map
#'
#' The conventional vertical person-item map (Wright and Stone 1979): the
#' person distribution to the left of a shared logit axis and the item
#' thresholds, labelled by item (and threshold number for polytomous items),
#' stacked to its right.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param bins Number of bins for the person distribution and the threshold
#'   label rows.
#' @param xlim Optional logit range for the shared scale; persons and
#'   thresholds outside it are omitted.
#' @param cex_labels Character expansion for the threshold labels.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @references Wright, B. D., & Stone, M. H. (1979). \emph{Best Test
#'   Design}. Chicago: MESA Press.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_wright(rasch(X))
#' @export
plot_wright <- function(fit, bins = 35, xlim = NULL, cex_labels = 0.8) {
  th <- fit$person$theta[!is.na(fit$person$theta)]
  thr <- fit$thresholds
  rng <- if (is.null(xlim)) range(c(th, thr$tau)) + c(-0.4, 0.4) else sort(xlim)
  th <- th[th >= rng[1] & th <= rng[2]]
  keep <- thr$tau >= rng[1] & thr$tau <= rng[2]
  tv <- thr$tau[keep]
  lab <- if (all(fit$m == 1L)) fit$items$item[thr$item[keep]] else
    paste0(fit$items$item[thr$item[keep]], ".", thr$k[keep])
  brk <- seq(rng[1], rng[2], length.out = bins + 1)
  hp <- hist(th, breaks = brk, plot = FALSE)
  pp <- hp$counts / max(1L, max(hp$counts))
  bin <- pmin(findInterval(tv, brk, rightmost.closed = TRUE), bins)
  split <- 0.42
  op <- par(mar = c(1.2, 4.4, 1.6, 0.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = rng, xlab = "", ylab = "Location (logits)",
       axes = FALSE, main = "")
  abline(h = pretty(rng), col = .rr$grid, lwd = 0.8)
  axis(2, col = .rr$grid, col.ticks = .rr$soft)
  segments(split, rng[1], split, rng[2], col = .rr$ink, lwd = 1)
  segments(split, mean(tv), 0.99, mean(tv), col = .rr$amber, lty = 2)
  rect(split - pp * (split - 0.02), brk[-length(brk)], split, brk[-1],
       col = .rr$blue, border = "white", lwd = 0.6)
  segments(0.02, mean(th), split, mean(th), col = .rr$ink, lty = 2)
  cexl <- cex_labels
  colw <- max(strwidth(lab, cex = cexl)) * 1.2
  x0 <- split + 0.015
  kmax <- max(1L, floor((1 - x0) / colw))
  for (b in unique(bin)) {
    ls <- lab[bin == b][order(tv[bin == b])]
    if (length(ls) > kmax)
      ls <- c(ls[seq_len(kmax - 1L)], paste0("+", length(ls) - kmax + 1L))
    text(x0 + (seq_along(ls) - 1L) * colw, (brk[b] + brk[b + 1L]) / 2, ls,
         cex = cexl, adj = 0, col = .rr$ink)
  }
  text(0.02, rng[2], sprintf("persons: mean %.2f, SD %.2f",
                             mean(th), stats::sd(th)),
       cex = 0.78, adj = c(0, 1), col = .rr$blue, font = 2)
  text(0.99, rng[2], sprintf("thresholds: mean %.2f, SD %.2f",
                             mean(tv), stats::sd(tv)),
       cex = 0.78, adj = c(1, 1), col = .rr$amber, font = 2)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Threshold map: every item's thresholds on the logit scale, by location.
# ---------------------------------------------------------------------------
#' Plot the threshold map
#'
#' Each item's threshold locations on a common logit scale, ordered by item
#' location, with disordered thresholds highlighted.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param order_by_location Order items by their location (the default)
#'   rather than their original sequence.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_threshold_map(rasch(X))
#' @export
plot_threshold_map <- function(fit, order_by_location = TRUE) {
  L <- length(fit$tau_list)
  ord <- if (order_by_location) order(fit$items$location) else seq_len(L)
  rng <- range(fit$thresholds$tau); rng <- rng + c(-0.4, 0.4)
  op <- par(mar = c(4.2, 7.5, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = rng, ylim = c(0.5, L + 0.5), xlab = "Location (logits)",
       ylab = "", axes = FALSE, main = "")
  abline(h = seq_len(L), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  axis(1, col = .rr$grid, col.ticks = .rr$soft)
  axis(2, at = seq_len(L), labels = fit$items$item[ord], cex.axis = 0.75,
       col = .rr$grid, col.ticks = NA)
  for (row in seq_len(L)) {
    i <- ord[row]; tau_i <- fit$tau_list[[i]]
    disord <- !(all(diff(tau_i) > 0) || length(tau_i) == 1L)
    segments(min(tau_i), row, max(tau_i), row, col = .rr$soft, lwd = 1.4)
    points(tau_i, rep(row, length(tau_i)), pch = 21, cex = 1.25,
           bg = if (disord) .rr$red else .rr$amber, col = "white", lwd = 1)
    points(mean(tau_i), row, pch = 23, bg = .rr$blue, col = "white", cex = 1.3)
  }
  .rr_legend("bottomright", c("threshold", "disordered", "item location"),
             pch = c(21, 21, 23), pt.bg = c(.rr$amber, .rr$red, .rr$blue),
             col = "white", pt.cex = 1.2)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Test characteristic curve and test information.
# ---------------------------------------------------------------------------
#' Plot the test characteristic curve
#'
#' Expected total score against person location for the whole instrument.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param grid Logit grid.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_tcc(rasch(X))
#' @export
plot_tcc <- function(fit, grid = seq(-6, 6, 0.05)) {
  Etot <- vapply(grid, function(th)
    sum(vapply(seq_along(fit$tau_list), function(i)
      item_moments(th, fit$tau_list[[i]], disc = .disc_of(fit, i))$E, 0)), 0)
  Smax <- sum(fit$m)
  op <- .rr_canvas(range(grid), c(0, Smax), "Person location (logits)",
                   "Expected total score")
  on.exit(par(op))
  lines(grid, Etot, lwd = 3, col = .rr$blue)
  abline(h = c(0, Smax), lty = 3, col = .rr$soft)
  invisible(NULL)
}

#' Plot the test information function
#'
#' Test information across the logit scale with the standard error of
#' measurement overlaid on a second axis.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param grid Logit grid.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_tif(rasch(X))
#' @export
plot_tif <- function(fit, grid = seq(-6, 6, 0.05)) {
  ti <- test_information(fit, grid)
  op <- .rr_canvas(range(grid), c(0, max(ti$info) * 1.1),
                   "Person location (logits)", "Test information",
                   right = 3.6)
  on.exit(par(op))
  polygon(c(ti$theta, rev(ti$theta)), c(ti$info, rep(0, nrow(ti))),
          col = paste0(.rr$blue, "22"), border = NA)
  lines(ti$theta, ti$info, lwd = 3, col = .rr$blue)
  sem <- ti$sem; sem[!is.finite(sem)] <- NA
  scl <- max(ti$info) * 1.05 / max(sem[ti$theta > -4 & ti$theta < 4], na.rm = TRUE)
  lines(ti$theta, sem * scl, lwd = 2.2, col = .rr$red, lty = 5)
  sem_ticks <- pretty(c(0, max(sem[ti$theta > -4 & ti$theta < 4], na.rm = TRUE)))
  sem_ticks <- sem_ticks[sem_ticks * scl <= max(ti$info) * 1.1]
  axis(4, at = sem_ticks * scl, labels = sem_ticks,
       col = .rr$grid, col.ticks = .rr$soft, col.axis = .rr$red, cex.axis = 0.8)
  mtext("SEM", side = 4, line = 2.3, col = .rr$red, cex = 0.85)
  .rr_legend("topleft", c("Information", "SEM"), lwd = c(3, 2.2),
             lty = c(1, 5), col = c(.rr$blue, .rr$red))
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Item and person fit maps.
# ---------------------------------------------------------------------------
#' Plot the item map (location against fit residual)
#'
#' Items plotted by location and fit residual, with the conventional
#' acceptance band at +/- 2.5 and misfitting items labelled.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param band Fit residual acceptance band.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_item_map(rasch(X))
#' @export
plot_item_map <- function(fit, band = 2.5) {
  d <- fit$items
  ylim <- range(c(d$fit_resid, -band, band), na.rm = TRUE) * 1.2
  op <- .rr_canvas(range(d$location) + c(-0.5, 0.5), ylim,
                   "Item location (logits)", "Fit residual",
                   grid_x = TRUE)
  on.exit(par(op))
  rect(par("usr")[1], -band, par("usr")[2], band,
       col = paste0(.rr$teal, "11"), border = NA)
  abline(h = c(-band, band), lty = 2, col = .rr$soft)
  out <- !is.na(d$fit_resid) & abs(d$fit_resid) > band
  points(d$location, d$fit_resid, pch = 21, cex = 1.7,
         bg = ifelse(out, .rr$red, .rr$blue), col = "white", lwd = 1.2)
  text(d$location, d$fit_resid, d$item, pos = 3, offset = 0.5, cex = 0.7,
       col = ifelse(out, .rr$red, .rr$soft))
  invisible(NULL)
}

#' Plot person fit
#'
#' Person locations against person fit residuals with the +/- 2.5 band; persons
#' beyond the band respond erratically (positive) or too deterministically
#' (negative).
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param band Fit residual acceptance band.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_person_fit(rasch(X))
#' @export
plot_person_fit <- function(fit, band = 2.5) {
  p <- fit$person
  ok <- !is.na(p$theta) & !is.na(p$fit_resid)
  ylim <- range(c(p$fit_resid[ok], -band - 0.5, band + 0.5))
  op <- .rr_canvas(range(p$theta[ok]) + c(-0.3, 0.3), ylim,
                   "Person location (logits)", "Fit residual",
                   grid_x = TRUE)
  on.exit(par(op))
  rect(par("usr")[1], -band, par("usr")[2], band,
       col = paste0(.rr$teal, "11"), border = NA)
  abline(h = c(-band, band), lty = 2, col = .rr$soft)
  out <- abs(p$fit_resid[ok]) > band
  points(p$theta[ok], p$fit_resid[ok], pch = 21, cex = 0.9,
         bg = ifelse(out, .rr$red, paste0(.rr$blue, "99")), col = "white", lwd = 0.5)
  mtext(sprintf("%d of %d persons beyond +/-%.1f", sum(out), sum(ok), band),
        side = 3, line = 0.2, adj = 0, cex = 0.8, col = .rr$ink)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Residual structure.
# ---------------------------------------------------------------------------
#' Plot the residual-correlation heatmap
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_resid_cor(rasch(X))
#' @export
plot_resid_cor <- function(fit) {
  R <- residual_correlations(fit)$matrix; L <- ncol(R)
  diag(R) <- NA   # self-correlations carry no dependence information
  pal <- colorRampPalette(c("#1d4ed8", "#f8fafc", "#dc2626"))(64)
  op <- par(mar = c(6, 6, 3.2, 3), las = 1, col.axis = .rr$ink,
            col.main = .rr$ink, font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  image(1:L, 1:L, R[, L:1, drop = FALSE], col = pal, zlim = c(-1, 1),
        axes = FALSE, xlab = "", ylab = "", main = "")
  axis(1, 1:L, colnames(R), las = 2, cex.axis = 0.65, col = NA, col.ticks = NA)
  axis(2, 1:L, rev(colnames(R)), cex.axis = 0.65, col = NA, col.ticks = NA)
  # compact colour key
  ky <- seq(0.25, 0.75, length.out = 65) * L
  rect(L + 0.62, ky[-65], L + 0.95, ky[-1], col = pal, border = NA, xpd = TRUE)
  text(L + 0.8, c(min(ky), max(ky)) + c(-0.4, 0.4), c("-1", "+1"),
       cex = 0.65, xpd = TRUE, col = .rr$ink)
  invisible(NULL)
}

#' Plot residual principal-component loadings
#'
#' First-contrast loadings against item location; opposing clusters at top and
#' bottom suggest a second dimension.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_pca(rasch(X))
#' @export
plot_pca <- function(fit) {
  pc <- residual_pca(fit)
  ld <- pc$loadings$pc1_loading[match(fit$items$item, pc$loadings$item)]
  loc <- fit$items$location
  op <- .rr_canvas(range(loc) + c(-0.5, 0.5), range(ld) * 1.25,
                   "Item location (logits)", "PC1 loading",
                   sprintf("eigenvalue %.3f  (%.1f%% of residual variance)",
                           pc$first_eigen, 100 * pc$prop[1]), grid_x = TRUE)
  on.exit(par(op))
  abline(h = 0, lty = 2, col = .rr$soft)
  points(loc, ld, pch = 21, cex = 1.6,
         bg = ifelse(ld > 0, .rr$blue, .rr$amber), col = "white", lwd = 1.2)
  text(loc, ld, fit$items$item, pos = 3, offset = 0.45, cex = 0.7, col = .rr$soft)
  invisible(NULL)
}

#' Plot category frequencies
#'
#' Observed response distribution over the categories of one item.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(1:4, function(i) sapply(th, function(t) sample(0:3, 1, prob = simP(t, c(-1, 0, 1)))))
#' colnames(X) <- sprintf("P%02d", 1:4)
#' plot_catfreq(rasch(X), "P01")
#' @export
plot_catfreq <- function(fit, item) {
  i <- .item_idx(fit, item)
  cnt <- fit$thresholds_diag[[i]]$category_counts
  cats <- seq_along(cnt) - 1L
  op <- .rr_canvas(c(-0.6, max(cats) + 0.6), c(0, max(cnt) * 1.12),
                   "Category", "Count",
                   fit$items$item[i],
                   grid_x = FALSE)
  on.exit(par(op))
  rect(cats - 0.38, 0, cats + 0.38, cnt, col = .rr$blue, border = "white")
  text(cats, cnt, cnt, pos = 3, cex = 0.8, col = .rr$ink)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Person characteristic curve: one person against the item difficulties.
# ---------------------------------------------------------------------------
#' Plot a person characteristic curve
#'
#' The person characteristic curve: the probability of success as a
#' function of item location at the person's estimated measure, with the
#' person's observed responses overlaid, grouped into item-difficulty
#' intervals (proportion of maximum score per interval). Erratic responding
#' (for example lucky guessing on hard items by a low-proficiency person)
#' shows as observed points far from the curve, complementing the person
#' fit residual.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param person Row number of the person, or an ID matching
#'   \code{fit$person$id}.
#' @param n_groups Number of item-difficulty intervals for the observed
#'   points (capped by the number of observed items).
#' @param grid Item-location grid over which to draw the curve.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 12)
#' X <- matrix(rbinom(300 * 12, 1, plogis(outer(rnorm(300), d, "-"))), 300, 12)
#' colnames(X) <- paste0("I", 1:12)
#' plot_pcc(rasch(X), person = 1)
#' @export
plot_pcc <- function(fit, person, n_groups = 5, grid = seq(-5, 5, 0.05)) {
  n <- if (is.numeric(person) && length(person) == 1L &&
           person %in% seq_len(nrow(fit$X))) as.integer(person)
       else match(person, fit$person$id)
  if (is.na(n)) stop("person not found")
  th <- fit$person$theta[n]
  if (is.na(th)) stop("no estimate for this person")
  x <- fit$X[n, ]; ok <- !is.na(x)
  if (sum(ok) < 3) stop("fewer than 3 observed responses for this person")
  loc <- fit$items$location; mm <- fit$m
  op <- .rr_canvas(range(grid), c(0, 1), "Item location (logits)",
                   "Probability of success",
                   sprintf("%s  (location %.3f, fit residual %s)",
                           fit$person$id[n], th,
                           ifelse(is.na(fit$person$fit_resid[n]), "NA",
                                  sprintf("%.2f", fit$person$fit_resid[n]))))
  on.exit(par(op))
  lines(grid, plogis(th - grid), lwd = 3, col = .rr$ink)
  abline(v = th, lty = 3, col = .rr$soft)
  g <- cut(rank(loc[ok], ties.method = "first"),
           min(n_groups, max(2, floor(sum(ok) / 2))), labels = FALSE)
  obsL <- tapply(loc[ok], g, mean)
  obsP <- tapply((x[ok] / mm[ok]), g, mean)
  points(obsL, obsP, pch = 21, bg = .rr$blue, col = "white", cex = 1.6, lwd = 1.2)
  .rr_legend("topright", c("Model at person location",
                           "Observed (item intervals)"),
             lwd = c(3, NA), pch = c(NA, 21), pt.bg = c(NA, .rr$blue),
             col = c(.rr$ink, "white"), pt.cex = 1.4)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Kidmap: the person diagnostic map (Wright, Mead & Ludlow 1980).
# ---------------------------------------------------------------------------
#' Plot a kidmap
#'
#' The person diagnostic map (Wright, Mead and Ludlow 1980): item thresholds
#' the person achieved to the right of a vertical logit axis and thresholds
#' not achieved to the left, with the person's location drawn as a dashed
#' line inside its confidence band. Achieved thresholds above the band
#' (unexpected successes) and unachieved thresholds below it (unexpected
#' failures) are highlighted; a clean response pattern shows achieved
#' thresholds below the band and unachieved ones above it.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param person Row number of the person, or an ID matching
#'   \code{fit$person$id}.
#' @param level Confidence level of the band around the person location used
#'   to mark unexpected responses.
#' @param bins Number of vertical bins used to stack the threshold labels.
#' @param xlim Optional logit range; thresholds outside it are omitted.
#' @param cex_labels Character expansion for the threshold labels.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @references Wright, B. D., Mead, R. J., & Ludlow, L. H. (1980).
#'   \emph{KIDMAP: person-by-item interaction mapping} (Research Memorandum
#'   No. 29). Chicago: University of Chicago, MESA Psychometric Laboratory.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 12)
#' X <- matrix(rbinom(300 * 12, 1, plogis(outer(rnorm(300), d, "-"))), 300, 12)
#' colnames(X) <- paste0("I", 1:12)
#' plot_kidmap(rasch(X), person = 1)
#' @export
plot_kidmap <- function(fit, person, level = 0.95, bins = 35, xlim = NULL,
                        cex_labels = 0.8) {
  n <- if (is.numeric(person) && length(person) == 1L &&
           person %in% seq_len(nrow(fit$X))) as.integer(person)
       else match(person, fit$person$id)
  if (is.na(n)) stop("person not found")
  th <- fit$person$theta[n]
  if (is.na(th)) stop("no estimate for this person")
  se <- fit$person$se[n]
  z <- stats::qnorm(1 - (1 - level) / 2)
  band <- if (is.na(se)) 0 else z * se
  x <- fit$X[n, ]
  thr <- fit$thresholds
  obs <- !is.na(x[thr$item])
  ii <- thr$item[obs]; kk <- thr$k[obs]; tv <- thr$tau[obs]
  att <- x[ii] >= kk
  lab <- if (all(fit$m == 1L)) fit$items$item[ii] else
    paste0(fit$items$item[ii], ".", kk)
  rng <- if (is.null(xlim))
    range(c(th - band, th + band, tv)) + c(-0.5, 0.5) else sort(xlim)
  keep <- tv >= rng[1] & tv <= rng[2]
  tv <- tv[keep]; att <- att[keep]; lab <- lab[keep]
  unexp <- (att & tv > th + band) | (!att & tv < th - band)
  brk <- seq(rng[1], rng[2], length.out = bins + 1)
  bin <- pmin(findInterval(tv, brk, rightmost.closed = TRUE), bins)
  split <- 0.5
  op <- par(mar = c(2.6, 4.4, 1.6, 0.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = rng, xlab = "", ylab = "Location (logits)",
       axes = FALSE, main = "")
  abline(h = pretty(rng), col = .rr$grid, lwd = 0.8)
  axis(2, col = .rr$grid, col.ticks = .rr$soft)
  if (band > 0)
    rect(0, th - band, 1, th + band, col = adjustcolor(.rr$blue, 0.10),
         border = NA)
  segments(split, rng[1], split, rng[2], col = .rr$ink, lwd = 1)
  segments(0, th, 1, th, col = .rr$blue, lty = 2, lwd = 1.4)
  cexl <- cex_labels
  colw <- max(strwidth(lab, cex = cexl)) * 1.2
  for (side in c(TRUE, FALSE)) {           # TRUE = achieved (right)
    sel <- att == side
    for (b in unique(bin[sel])) {
      pick <- sel & bin == b
      ls <- lab[pick][order(tv[pick])]
      cols <- ifelse(unexp[pick][order(tv[pick])], .rr$red, .rr$ink)
      kmax <- max(1L, floor((split - 0.015) / colw))
      if (length(ls) > kmax) {
        cols <- c(cols[seq_len(kmax - 1L)], .rr$soft)
        ls <- c(ls[seq_len(kmax - 1L)], paste0("+", length(ls) - kmax + 1L))
      }
      xs <- if (side) split + 0.015 + (seq_along(ls) - 1L) * colw
            else split - 0.015 - (seq_along(ls) - 1L) * colw
      text(xs, (brk[b] + brk[b + 1L]) / 2, ls, cex = cexl,
           adj = if (side) 0 else 1, col = cols, font = 2)
    }
  }
  fr <- fit$person$fit_resid[n]
  text(0.01, rng[2], sprintf("%s: location %.2f (SE %s), fit residual %s",
                             fit$person$id[n], th,
                             ifelse(is.na(se), "NA", sprintf("%.2f", se)),
                             ifelse(is.na(fr), "NA", sprintf("%.2f", fr))),
       cex = 0.78, adj = c(0, 1), col = .rr$ink, font = 2)
  mtext("not achieved", side = 1, line = 0.4, at = split / 2, cex = 0.8,
        col = .rr$soft)
  mtext("achieved", side = 1, line = 0.4, at = (1 + split) / 2, cex = 0.8,
        col = .rr$soft)
  if (any(unexp))
    legend("bottomleft", "unexpected response", bty = "n",
           text.col = .rr$red, cex = 0.85)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Residual statistics distributions of the fit residuals.
# ---------------------------------------------------------------------------
#' Plot the fit residual distribution
#'
#' A histogram of the
#' item or person fit residuals -- the log-transformed statistic or its
#' untransformed natural form -- against the standard normal density they
#' should approximate under fit (Andrich and Marais 2019, ch. 15). The
#' natural residual is visibly skewed
#' (that is why the log transform is reported); both are available.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param what \code{"items"} or \code{"persons"}.
#' @param statistic \code{"fit_resid"} (log-transformed, default) or
#'   \code{"natural"}.
#' @param bins Number of histogram bins.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 10)
#' X <- matrix(rbinom(400 * 10, 1, plogis(outer(rnorm(400), d, "-"))), 400, 10)
#' colnames(X) <- paste0("I", 1:10)
#' plot_resid_dist(rasch(X), what = "persons")
#' @export
plot_resid_dist <- function(fit, what = c("items", "persons"),
                            statistic = c("fit_resid", "natural"), bins = 25) {
  what <- match.arg(what); statistic <- match.arg(statistic)
  v <- if (what == "items") {
    if (statistic == "fit_resid") fit$items$fit_resid else fit$items$natural_resid
  } else {
    if (statistic == "fit_resid") fit$person$fit_resid else fit$person$natural_resid
  }
  v <- v[!is.na(v)]
  if (length(v) < 3) stop("fewer than 3 residuals to display")
  lab <- if (statistic == "fit_resid") "fit residual (log-transformed)"
         else "natural fit residual"
  rng <- range(c(v, -3, 3)); rng <- rng + c(-0.5, 0.5) * diff(rng) * 0.05
  brk <- seq(rng[1], rng[2], length.out = bins + 1)
  h <- hist(v, breaks = brk, plot = FALSE)
  ymax <- max(h$density, dnorm(0)) * 1.15
  op <- .rr_canvas(rng, c(0, ymax), lab, "Density",
                   sprintf("n = %d, mean %.2f, SD %.2f", length(v), mean(v),
                           sd(v)), grid_x = FALSE)
  on.exit(par(op))
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density,
       col = adjustcolor(.rr$blue, 0.45), border = "white")
  xs <- seq(rng[1], rng[2], length.out = 200)
  lines(xs, dnorm(xs), lwd = 2.6, col = .rr$red)
  abline(v = c(-2.5, 2.5), lty = 3, col = .rr$soft)
  .rr_legend("topright", c("Observed", "Standard normal"),
             lwd = c(NA, 2.6), pch = c(22, NA), pt.bg = c(adjustcolor(.rr$blue, 0.45), NA),
             col = c("white", .rr$red), pt.cex = 1.6)
  invisible(NULL)
}
