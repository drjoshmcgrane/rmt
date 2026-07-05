# rmt :: many-facet Rasch model
# ===========================================================================
# The many-facet Rasch model (Linacre 1989) estimated by the same pairwise
# conditional likelihood as the rest of the package. Each combination of an
# item with the levels of the rating facets (for example item x rater) is a
# "virtual item" whose thresholds decompose structurally as
#
#   tau_{(i, r), k} = delta_{ik} + rho_r (+ further facet terms)
#
# so the facet severities enter the design matrix exactly as the rating
# scale structure does, and the person parameter still cancels within every
# pair of virtual items. Identification: item thresholds sum to zero and
# each facet's severities sum to zero. The pairwise conditional likelihood
# is concave in the structural parameters, so Newton-Raphson from zero
# converges globally.
# ===========================================================================

# Pooled fit over a set of columns of the residual matrix (used for facet
# levels and for underlying items across their virtual columns). The fit
# residual is the log-of-mean-square statistic (Andrich & Marais 2019,
# ch. 23) pooled over the group's
# observed cells of non-extreme persons: Y^2 = sum z^2 against the summed
# cell degrees of freedom f = f_cell x cells, symmetrised as in .fitres.
# The published three-facet fit tables report a different margin statistic,
# the mean of the virtual items' fit residuals, which rasch_mfrm() computes
# alongside this pooled form (see its Details); EFRM frame fit keeps the
# pooled form. Infit and outfit mean squares are kept alongside.
.group_col_fit <- function(Z, mo, cols, disc = NULL, extreme = NULL,
                           f_cell = NA_real_) {
  E2 <- .z2_expectation(mo, Z, disc)
  sub <- Z[, cols, drop = FALSE]
  keep <- if (is.null(extreme)) rep(TRUE, nrow(Z)) else !extreme
  ok <- which(!is.na(sub) & keep)
  if (length(ok) < 3)
    return(list(infit_ms = NA_real_, outfit_ms = NA_real_,
                fit_resid = NA_real_, df_fit = NA_real_, n = length(ok)))
  z2 <- sub[ok]^2
  V <- mo$V[, cols, drop = FALSE][ok]
  C4 <- mo$M4[, cols, drop = FALSE][ok]
  e2 <- E2[, cols, drop = FALSE][ok]
  n <- length(ok)
  outfit <- sum(z2) / sum(e2)
  infit <- sum(z2 * V) / sum(e2 * V)
  fr <- df <- NA_real_
  if (is.finite(f_cell) && f_cell > 0) {
    y2 <- sum(z2); v <- sum(C4 / V^2 - 1); f <- f_cell * n
    if (v > 1e-8 && y2 > 0) {
      fr <- f * (log(y2) - log(f)) / sqrt(v)
      df <- f
    }
  }
  if (!is.finite(fr)) {                     # degenerate pooling: fall back
    qi <- sqrt(max(sum(C4 - V^2) / sum(V)^2, 1e-8))
    fr <- .wh(infit, qi)
  }
  list(infit_ms = infit, outfit_ms = outfit, fit_resid = fr, df_fit = df,
       n = n)
}

#' Fit a many-facet Rasch model
#'
#' Estimates the many-facet Rasch model (Linacre 1989) for long-format data
#' in which each row is one scored response carrying a person, an item, a
#' score, and one or more facet levels (for example the rater). Every
#' item-by-facet combination becomes a virtual item whose thresholds are the
#' item's thresholds shifted by the facet severities, and the whole structure
#' is estimated in one pass of the pairwise conditional likelihood, in which
#' the person parameter cancels. Facet severities are reported with standard
#' errors and pooled fit statistics; the returned object is also a full
#' \code{\link{rasch}} fit at the virtual-item level, so every diagnostic
#' table and plot in the package applies to it.
#'
#' @param data Long-format data frame.
#' @param person Name of the person identifier column.
#' @param item Name of the item column.
#' @param score Name of the integer score column (categories from 0; gaps are
#'   collapsed per item with a note).
#' @param facets Character vector naming one or more facet columns (for
#'   example a rater column).
#' @param n_groups Number of class intervals for the item-trait chi-square;
#'   \code{NULL} (the default) applies the class-interval rule of Andrich and
#'   Marais (2019, ch. 15) (at least 50
#'   non-extreme persons per interval, at most 10 intervals, at least 2).
#' @param adjust_N Optional reference sample size for the chi-square.
#' @param na_codes Score values to read as missing (default \code{-1}); any
#'   negative score is also treated as missing.
#' @param items Optional character vector of item score columns for data in
#'   wide format: one row per person-by-facet combination (for example one
#'   row per script per rater) with one column per item or criterion. The
#'   long form (\code{item} + \code{score}) remains available for data
#'   where the facet varies within items.
#' @param interaction Optional name of one facet to interact with the items
#'   (interactive facet mode). Adds item-by-facet terms
#'   \code{gamma[item, level]} with double sum-to-zero constraints on top of
#'   the additive severities, so each level may be more or less severe on
#'   particular items; estimates are returned in \code{interaction_effects}.
#'   The interactive model remains in the Rasch class (all discriminations
#'   equal one and the parameters are additive), but a significant
#'   interaction qualifies specific objectivity in practice: comparisons of
#'   the interacting facet's levels become item-dependent, which is itself
#'   the substantive finding.
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return An object of classes \code{"rasch_mfrm"} and \code{"rasch"}. In
#'   addition to every component of a \code{\link{rasch}} fit (computed over
#'   the virtual items), it carries \code{facet_effects} (per facet: level,
#'   severity, standard error, observation count, pooled fit),
#'   \code{item_effects} (underlying item locations and pooled fit),
#'   \code{item_thresholds} (the structural \code{delta_ik} with standard
#'   errors), and \code{facet_spec}. Two fit residuals are reported per
#'   facet level and per underlying item. \code{fit_resid} is the
#'   facet-margin statistic of the published three-facet fit tables
#'   (Andrich and Marais 2019, ch. 26 and app. C), the mean of the
#'   constituent virtual items'
#'   fit residuals; it weighs each virtual item equally, so an erratic level
#'   shows the average of its per-item misfit. \code{fit_resid_pooled} is
#'   the log-of-mean-square statistic summed over the margin's
#'   observed cells of non-extreme persons, with its degrees of freedom in
#'   \code{df_fit}; it weighs each response equally and is the more
#'   powerful statistic when misfit is spread evenly over the level's
#'   cells.
#' @examples
#' set.seed(1)
#' simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
#' persons <- sprintf("P%03d", 1:120); raters <- paste0("R", 1:4)
#' th <- setNames(rnorm(120, 0, 1.3), persons)
#' rho <- setNames(c(-0.6, -0.2, 0.2, 0.6), raters)
#' tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
#' d <- expand.grid(person = persons, item = names(tau), rater = raters,
#'                  stringsAsFactors = FALSE)
#' d$score <- mapply(function(p, i, r)
#'   sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho[r])), d$person, d$item, d$rater)
#' fit <- rasch_mfrm(d, person = "person", item = "item", score = "score",
#'                   facets = "rater")
#' fit$facet_effects$rater
#' @export
rasch_mfrm <- function(data, person, item = NULL, score = NULL, facets,
                       items = NULL, n_groups = NULL,
                       adjust_N = NA, na_codes = -1, interaction = NULL,
                       maxit = 60, tol = 1e-8) {
  # wide entry: item score columns are melted to the long form internally
  if (!is.null(items)) {
    if (!is.null(item) || !is.null(score))
      stop("give either `items` (wide: one column per item) or `item` + `score` (long)")
    miss <- setdiff(c(person, facets, items), names(data))
    if (length(miss)) stop("column(s) not in data: ", paste(miss, collapse = ", "))
    long <- data.frame(
      ..person = rep(as.character(data[[person]]), length(items)),
      ..item = rep(items, each = nrow(data)),
      ..score = unlist(lapply(items, function(cn)
        suppressWarnings(as.numeric(data[[cn]])))),
      stringsAsFactors = FALSE)
    for (f in facets) long[[f]] <- rep(as.character(data[[f]]), length(items))
    return(rasch_mfrm(long, person = "..person", item = "..item",
                      score = "..score", facets = facets,
                      n_groups = n_groups, adjust_N = adjust_N,
                      na_codes = na_codes, interaction = interaction,
                      maxit = maxit, tol = tol))
  }
  if (is.null(item) || is.null(score))
    stop("give either `items` (wide) or `item` + `score` (long)")
  if (!is.null(interaction)) {
    interaction <- as.character(interaction)[1]
    if (!interaction %in% facets)
      stop("'interaction' must name one of the facets")
  }
  stopifnot(is.data.frame(data))
  need <- c(person, item, score, facets)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("column(s) not in data: ", paste(miss, collapse = ", "))
  notes <- character(0)

  pid <- as.character(data[[person]])
  itm <- as.character(data[[item]])
  sc <- suppressWarnings(as.integer(as.character(data[[score]])))
  n_na <- sum(!is.na(sc) & (sc %in% na_codes | sc < 0))
  sc[sc %in% na_codes | (!is.na(sc) & sc < 0)] <- NA
  if (n_na > 0)
    notes <- c(notes, sprintf("%d response(s) with a missing-data code (%s) set to missing",
                              n_na, paste(unique(c(na_codes, "negative")), collapse = ", ")))
  if (all(is.na(sc))) stop("score column has no usable integer values")
  fac <- lapply(facets, function(f) as.character(data[[f]]))
  names(fac) <- facets

  # rescore each item to consecutive categories from 0
  items_u <- sort(unique(itm))
  item_m <- setNames(integer(length(items_u)), items_u)
  for (it in items_u) {
    sel <- which(itm == it & !is.na(sc))
    obs <- sort(unique(sc[sel]))
    if (length(obs) < 2) stop("item ", it, " has fewer than two observed categories")
    full <- seq(0L, max(obs))
    if (!identical(obs, full)) {
      sc[sel] <- match(sc[sel], obs) - 1L
      notes <- c(notes, sprintf("item %s rescored: observed categories [%s] mapped to 0:%d",
                                it, paste(obs, collapse = ","), length(obs) - 1L))
    }
    item_m[it] <- length(sort(unique(sc[sel]))) - 1L
  }

  # virtual items: item x facet-level combinations present in the data
  fkey <- do.call(paste, c(fac, list(sep = ":")))
  vkey <- paste(itm, fkey, sep = ":")
  vlev <- unique(vkey[order(match(itm, items_u), fkey)])
  vmap <- data.frame(vkey = vlev,
                     item = itm[match(vlev, vkey)],
                     stringsAsFactors = FALSE)
  for (f in facets) vmap[[f]] <- fac[[f]][match(vlev, vkey)]

  persons_u <- unique(pid)
  Xv <- matrix(NA_integer_, length(persons_u), length(vlev),
               dimnames = list(NULL, vlev))
  ri <- match(pid, persons_u); cj <- match(vkey, vlev)
  dup <- duplicated(cbind(ri, cj))
  if (any(dup))
    notes <- c(notes, sprintf("%d duplicate person-by-virtual-item response(s) ignored (first kept)",
                              sum(dup)))
  use <- !dup & !is.na(sc)
  Xv[cbind(ri[use], cj[use])] <- sc[use]

  m_v <- item_m[vmap$item]
  thr_v <- threshold_index(m_v)

  # --- structural design matrix --------------------------------------------
  thr_items <- threshold_index(item_m[items_u])     # delta enumeration
  Md <- nrow(thr_items)
  A_delta <- rbind(diag(Md - 1L), rep(-1, Md - 1L))
  flevs <- lapply(facets, function(f) sort(unique(fac[[f]])))
  names(flevs) <- facets
  A_fac <- lapply(flevs, function(lv) {
    Lf <- length(lv)
    if (Lf < 2) stop("a facet needs at least two levels")
    rbind(diag(Lf - 1L), rep(-1, Lf - 1L))
  })
  Li <- length(items_u)
  A_item <- rbind(diag(Li - 1L), rep(-1, Li - 1L))   # item-location margins
  n_gam <- if (is.null(interaction)) 0L else
    (Li - 1L) * (length(flevs[[interaction]]) - 1L)
  P <- (Md - 1L) + sum(vapply(A_fac, ncol, 1L)) + n_gam
  B <- matrix(0, nrow(thr_v), P)
  for (row in seq_len(nrow(thr_v))) {
    v <- thr_v$item[row]; k <- thr_v$k[row]
    iu <- match(vmap$item[v], items_u)
    drow <- which(thr_items$item == iu & thr_items$k == k)
    B[row, seq_len(Md - 1L)] <- A_delta[drow, ]
    cursor <- Md - 1L
    for (f in facets) {
      lev <- match(vmap[[f]][v], flevs[[f]])
      nc <- ncol(A_fac[[f]])
      B[row, cursor + seq_len(nc)] <- A_fac[[f]][lev, ]
      cursor <- cursor + nc
    }
    if (n_gam > 0L) {
      lev <- match(vmap[[interaction]][v], flevs[[interaction]])
      B[row, cursor + seq_len(n_gam)] <-
        as.vector(outer(A_item[iu, ], A_fac[[interaction]][lev, ]))
    }
  }

  # concave likelihood: Newton-Raphson from zero with step halving
  sol <- .pcml_solve(Xv, thr_v, m_v, B, rep(0, P), maxit = maxit, tol = tol)

  thr_v$tau <- sol$tau; thr_v$se <- sol$se_tau; thr_v$anchored <- FALSE
  est <- list(model = "MFRM", thr = thr_v, cov_tau = sol$cov_tau,
              loglik = sol$loglik, iterations = sol$iterations,
              converged = sol$converged, m = m_v, anchors = NULL,
              n_parameters = P)

  fit <- .assemble_fit("MFRM", Xv, est, persons_u, NULL, n_groups, adjust_N,
                       notes)

  # --- structural effects -----------------------------------------------------
  covb <- sol$cov_beta
  d_idx <- seq_len(Md - 1L)
  delta <- drop(A_delta %*% sol$beta[d_idx])
  cov_d <- A_delta %*% covb[d_idx, d_idx, drop = FALSE] %*% t(A_delta)
  fit$item_thresholds <- data.frame(item = items_u[thr_items$item],
                                    k = thr_items$k, tau = delta,
                                    se = sqrt(pmax(diag(cov_d), 0)))
  item_fit <- lapply(items_u, function(it)
    .group_col_fit(fit$residuals, fit$moments, which(vmap$item == it),
                    extreme = fit$person$extreme,
                    f_cell = fit$summary_stats$df_factor))
  # The facet-margin fit residual of the published three-facet tables
  # (Andrich & Marais 2019, ch. 26 and app. C) is the MEAN of the
  # constituent virtual items' fit residuals; the pooled log-residual over
  # the margin's cells is kept alongside with its degrees of freedom.
  vmean <- function(sel) mean(fit$items$fit_resid[sel], na.rm = TRUE)
  fit$item_effects <- data.frame(
    item = items_u,
    location = vapply(seq_along(items_u), function(i)
      mean(delta[thr_items$item == i]), 0),
    se = vapply(seq_along(items_u), function(i) {
      rows <- which(thr_items$item == i)
      sqrt(mean(cov_d[rows, rows]))
    }, 0),
    n = vapply(item_fit, `[[`, 0, "n"),
    infit_ms = vapply(item_fit, `[[`, 0, "infit_ms"),
    outfit_ms = vapply(item_fit, `[[`, 0, "outfit_ms"),
    fit_resid = vapply(items_u, function(it) vmean(vmap$item == it), 0),
    fit_resid_pooled = vapply(item_fit, `[[`, 0, "fit_resid"),
    df_fit = vapply(item_fit, `[[`, 0, "df_fit"))

  cursor <- Md - 1L
  fit$facet_effects <- list()
  for (f in facets) {
    nc <- ncol(A_fac[[f]]); idx <- cursor + seq_len(nc); cursor <- cursor + nc
    rho <- drop(A_fac[[f]] %*% sol$beta[idx])
    cov_r <- A_fac[[f]] %*% covb[idx, idx, drop = FALSE] %*% t(A_fac[[f]])
    lev_fit <- lapply(flevs[[f]], function(lv)
      .group_col_fit(fit$residuals, fit$moments, which(vmap[[f]] == lv),
                      extreme = fit$person$extreme,
                      f_cell = fit$summary_stats$df_factor))
    fit$facet_effects[[f]] <- data.frame(
      level = flevs[[f]], severity = rho,
      se = sqrt(pmax(diag(cov_r), 0)),
      n = vapply(lev_fit, `[[`, 0, "n"),
      infit_ms = vapply(lev_fit, `[[`, 0, "infit_ms"),
      outfit_ms = vapply(lev_fit, `[[`, 0, "outfit_ms"),
      fit_resid = vapply(flevs[[f]], function(lv) vmean(vmap[[f]] == lv), 0),
      fit_resid_pooled = vapply(lev_fit, `[[`, 0, "fit_resid"),
      df_fit = vapply(lev_fit, `[[`, 0, "df_fit"))
  }
  if (n_gam > 0L) {
    idx <- cursor + seq_len(n_gam)
    R0 <- length(flevs[[interaction]])
    K <- A_fac[[interaction]] %x% A_item        # vec(gamma) = K beta_gamma
    gvec <- drop(K %*% sol$beta[idx])
    cov_g <- K %*% covb[idx, idx, drop = FALSE] %*% t(K)
    fit$interaction_effects <- data.frame(
      item = rep(items_u, R0),
      level = rep(flevs[[interaction]], each = Li),
      gamma = gvec, se = sqrt(pmax(diag(cov_g), 0)))
    fit$interaction <- interaction
  }
  fit$facet_spec <- facets
  fit$virtual_map <- vmap
  class(fit) <- c("rasch_mfrm", "rasch")
  fit
}

#' @export
print.rasch_mfrm <- function(x, ...) {
  cat(sprintf("rmt many-facet analysis: %d items x %s = %d virtual items, %d persons\n",
              nrow(x$item_effects),
              paste(vapply(x$facet_spec, function(f)
                sprintf("%d %s level(s)", nrow(x$facet_effects[[f]]), f), ""),
                collapse = " x "),
              ncol(x$X), nrow(x$X)))
  cat(sprintf("Pairwise conditional ML: %s in %d iterations\n",
              if (x$est$converged) "converged" else "NOT converged",
              x$est$iterations))
  cat(sprintf("PSI %.3f, power of fit: %s\n", x$psi$PSI, x$power_of_fit))
  for (f in x$facet_spec) {
    fe <- x$facet_effects[[f]]
    core <- c("level", "severity", "se", "n", "fit_resid")
    cat(sprintf("\nFacet '%s' severities (logits):\n", f))
    print(.fmt_df(fe[, intersect(core, names(fe))]), row.names = FALSE)
  }
  cat("(pooled fit residuals and their df on fit$facet_effects)\n")
  if (!is.null(x$interaction)) {
    big <- x$interaction_effects
    big <- big[abs(big$gamma) > 1.96 * big$se, , drop = FALSE]
    cat(sprintf("\nItem-by-%s interactions (interactive facet mode): %d significant of %d\n",
                x$interaction, nrow(big), nrow(x$interaction_effects)))
    if (nrow(big)) print(big, digits = 3, row.names = FALSE)
  }
  if (length(x$notes)) cat("\nNotes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

#' Plot facet severities
#'
#' Caterpillar plot of the severity of each level of a facet from a
#' many-facet analysis, with 95 per cent error bars; levels with pooled fit
#' residuals beyond the band are highlighted.
#'
#' @param fit A fitted object from \code{\link{rasch_mfrm}}.
#' @param facet Facet name; defaults to the first facet.
#' @param band Fit residual band beyond which a level is highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' set.seed(1)
#' simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
#' persons <- sprintf("P%03d", 1:120); raters <- paste0("R", 1:4)
#' th <- setNames(rnorm(120, 0, 1.3), persons)
#' rho <- setNames(c(-0.6, -0.2, 0.2, 0.6), raters)
#' tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
#' d <- expand.grid(person = persons, item = names(tau), rater = raters,
#'                  stringsAsFactors = FALSE)
#' d$score <- mapply(function(p, i, r)
#'   sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho[r])), d$person, d$item, d$rater)
#' plot_facets(rasch_mfrm(d, "person", "item", "score", facets = "rater"))
#' }
#' @export
plot_facets <- function(fit, facet = NULL, band = 2.5) {
  if (!inherits(fit, "rasch_mfrm")) stop("plot_facets needs a rasch_mfrm fit")
  if (is.null(facet)) facet <- fit$facet_spec[1]
  fe <- fit$facet_effects[[facet]]
  if (is.null(fe)) stop("no such facet: ", facet)
  fe <- fe[order(fe$severity), ]
  lo <- fe$severity - 1.96 * fe$se; hi <- fe$severity + 1.96 * fe$se
  n <- nrow(fe)
  op <- par(mar = c(4.2, 7.5, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = range(c(lo, hi, 0)) + c(-0.2, 0.2), ylim = c(0.5, n + 0.5),
       xlab = "Severity (logits)", ylab = "", axes = FALSE, main = "")
  title(main = facet, adj = 0, line = 1.4)
  abline(h = seq_len(n), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  axis(1, col = .rr$grid, col.ticks = .rr$soft)
  axis(2, at = seq_len(n), labels = fe$level, cex.axis = 0.8,
       col = .rr$grid, col.ticks = NA)
  misfit <- !is.na(fe$fit_resid) & abs(fe$fit_resid) > band
  segments(lo, seq_len(n), hi, seq_len(n), lwd = 2.2,
           col = ifelse(misfit, .rr$red, .rr$soft))
  points(fe$severity, seq_len(n), pch = 21, cex = 1.5,
         bg = ifelse(misfit, .rr$red, .rr$blue), col = "white", lwd = 1.2)
  if (any(misfit))
    mtext(sprintf("%d level(s) with |fit residual| > %.1f", sum(misfit), band),
          side = 3, line = 0.2, adj = 0, cex = 0.8, col = .rr$red)
  invisible(NULL)
}
