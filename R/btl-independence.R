# rasch :: paired-comparison independence and dimensionality
#
# The Rasch residual-correlation (Q3) and residual-PCA displays operate on the
# persons x items residual matrix, which paired-comparison data do not
# produce (the datum is a comparison, so residuals live on object pairs). The
# analogues that DO fall out of the pair structure are transitivity (are the
# preferences consistent with a single order?) and a decomposition of the
# object-by-object residual preference matrix (is there structured departure
# from the one-dimensional fit, i.e. a second attribute steering contests?).

# weighted "points" object i scored against object j, as a K x K matrix
# (S[i, j]); resp is the per-comparison response (0..m for a, m - x for b).
.btl_scores <- function(ia, ib, resp, w, m, K) {
  S <- matrix(0, K, K)
  fill <- function(rows, cols, val) {
    idx <- (cols - 1L) * K + rows                    # column-major linear index
    ag <- rowsum(val, idx)
    S[as.integer(rownames(ag))] <<- S[as.integer(rownames(ag))] + ag[, 1]
  }
  fill(ia, ib, w * resp)
  fill(ib, ia, w * (m - resp))
  S
}

# count transitive vs circular triads in a tournament: D[i,j] = +1 if i beats
# j, -1 if j beats i, 0 tie; seen[i,j] whether the pair was compared. Returns
# the triple counts and, per object, how many circular triads it sits in.
.btl_triads <- function(D, seen) {
  K <- nrow(D)
  n_tri <- 0L; n_circ <- 0L
  invol <- integer(K)
  if (K >= 3L) for (i in 1:(K - 2L)) for (j in (i + 1L):(K - 1L)) {
    if (!seen[i, j] || D[i, j] == 0) next
    for (k in (j + 1L):K) {
      if (!(seen[i, k] && seen[j, k]) || D[i, k] == 0 || D[j, k] == 0) next
      n_tri <- n_tri + 1L
      # a triple's win-counts are {2,1,0} when transitive, {1,1,1} when circular
      si <- (D[i, j] > 0) + (D[i, k] > 0)
      sj <- (D[i, j] < 0) + (D[j, k] > 0)
      sk <- (D[i, k] < 0) + (D[j, k] < 0)
      if (si == 1L && sj == 1L && sk == 1L) {
        n_circ <- n_circ + 1L
        invol[i] <- invol[i] + 1L; invol[j] <- invol[j] + 1L
        invol[k] <- invol[k] + 1L
      }
    }
  }
  list(n_triples = n_tri, n_circular = n_circ, involvement = invol)
}

#' Transitivity of paired comparisons
#'
#' The single-dimension analogue for paired comparisons of the
#' unidimensionality question. A Bradley-Terry-Luce scale implies that
#' preferences stack into one consistent order: if A beats B and B beats C
#' then A should beat C. A \emph{circular triad} (A beats B, B beats C, C
#' beats A) is a local contradiction, like rock-paper-scissors. A few are
#' sampling noise; many, systematically, mean the comparisons are not being
#' driven by a single attribute. The rate of circular triads is compared with
#' the value expected from pure guessing (one quarter of triples), and, when
#' every pair has been compared, Kendall's coefficient of consistency is
#' reported (Kendall & Babington Smith 1940). With judges, each judge's own
#' consistency is reported too, flagging judges whose choices approach chance.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}}.
#' @param min_triples A judge is reported only if this many complete triples
#'   (all three pairs judged) are available.
#' @return A list of class \code{"rasch_btl_transitivity"}: \code{summary} (one
#'   row: objects, pairs compared, complete triples, circular triads, the
#'   circular rate, the chance rate 0.25, the consistency index
#'   \code{1 - rate/0.25}, and Kendall's \code{zeta} when the design is a
#'   complete round-robin with no exactly-tied pair -- \code{NA} otherwise);
#'   \code{objects} (each object's circular-triad
#'   involvement); \code{judges} (per-judge consistency, when judges exist);
#'   and \code{notes}.
#' @references Kendall, M. G., & Babington Smith, B. (1940). On the method of
#'   paired comparisons. \emph{Biometrika}, 31(3/4), 324-345.
#' @examples
#' set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
#' pr <- t(utils::combn(objs, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 20), b = rep(pr[, 2], each = 20))
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' btl_transitivity(btl(d, "a", "b", "win"))
#' @export
btl_transitivity <- function(fit, min_triples = 5L) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  objs <- fit$objects$object; K <- length(objs); m <- fit$m
  cmp <- fit$comparisons
  ia <- match(cmp$object_a, objs); ib <- match(cmp$object_b, objs)
  notes <- character(0)

  tournament <- function(rows) {
    S <- .btl_scores(ia[rows], ib[rows], cmp$response[rows], cmp$weight[rows],
                     m, K)
    seen <- (S + t(S)) > 0
    D <- sign(S - t(S)); D[!seen] <- 0
    list(D = D, seen = seen)
  }

  tt <- tournament(seq_len(nrow(cmp)))
  tri <- .btl_triads(tt$D, tt$seen)
  n_pairs <- sum(tt$seen[upper.tri(tt$seen)])
  complete <- n_pairs == choose(K, 2)
  ties <- sum(tt$seen & tt$D == 0) / 2
  rate <- if (tri$n_triples) tri$n_circular / tri$n_triples else NA_real_
  # Kendall's coefficient of consistency, defined for a complete round-robin
  zeta <- NA_real_
  if (complete && ties == 0 && K >= 3) {
    dmax <- if (K %% 2L) (K^3 - K) / 24 else (K^3 - 4 * K) / 24
    zeta <- 1 - tri$n_circular / dmax
  }
  if (!complete)
    notes <- c(notes, sprintf(
      "%d of %d object pairs were compared; triads use complete triples only",
      n_pairs, choose(K, 2)))
  if (ties > 0)
    notes <- c(notes, sprintf("%d pair(s) split exactly evenly and are set aside",
                              ties))

  summary <- data.frame(
    n_objects = K, n_pairs = n_pairs, n_triples = tri$n_triples,
    n_circular = tri$n_circular, circular_rate = rate,
    chance_rate = 0.25, consistency = 1 - rate / 0.25, zeta = zeta)

  objects <- data.frame(object = objs, circular_triads = tri$involvement)
  objects <- objects[order(-objects$circular_triads), ]
  rownames(objects) <- NULL

  judges <- NULL
  if (!all(is.na(cmp$judge))) {
    ju <- sort(unique(cmp$judge[!is.na(cmp$judge)]))
    rows <- lapply(ju, function(j) which(cmp$judge == j))
    jt <- lapply(rows, function(rr) {
      tj <- tournament(rr); .btl_triads(tj$D, tj$seen)
    })
    n_tri <- vapply(jt, `[[`, 0L, "n_triples")
    keep <- n_tri >= min_triples
    if (any(keep)) {
      n_c <- vapply(jt, `[[`, 0L, "n_circular")[keep]
      nt <- n_tri[keep]
      jr <- n_c / nt
      judges <- data.frame(
        judge = ju[keep], n_comparisons = vapply(rows[keep], length, 0L),
        n_triples = nt, n_circular = n_c, circular_rate = jr,
        consistency = 1 - jr / 0.25)
      judges <- judges[order(judges$consistency), ]
      rownames(judges) <- NULL
    } else {
      notes <- c(notes, sprintf(
        "no judge reached %d complete triples; per-judge consistency omitted",
        min_triples))
    }
  }

  out <- list(summary = summary, objects = objects, judges = judges,
              notes = notes)
  class(out) <- "rasch_btl_transitivity"
  out
}

#' @export
print.rasch_btl_transitivity <- function(x, ...) {
  s <- x$summary
  cat(sprintf("Paired-comparison transitivity: %d objects, %d complete triples\n",
              s$n_objects, s$n_triples))
  cat(sprintf("Circular triads: %d (%.1f%% of triples; chance %.0f%%) -> consistency %.2f\n",
              s$n_circular, 100 * s$circular_rate, 100 * s$chance_rate,
              s$consistency))
  if (!is.na(s$zeta))
    cat(sprintf("Kendall coefficient of consistency (complete design): %.3f\n",
                s$zeta))
  if (!is.null(x$judges))
    cat(sprintf("Per-judge consistency reported for %d judge(s); least consistent %.2f\n",
                nrow(x$judges), min(x$judges$consistency)))
  for (n in x$notes) cat("Note:", n, "\n")
  invisible(x)
}

# object-by-object residual log-odds matrix (skew-symmetric): observed
# pairwise log-odds minus the model difference beta_i - beta_j, with a
# continuity correction for extreme cells and zero where a pair is unseen.
.btl_resid_matrix <- function(ia, ib, resp, w, m, K, beta) {
  S <- .btl_scores(ia, ib, resp, w, m, K)
  tot <- S + t(S)
  P <- (S + 0.5) / (tot + 1)
  L <- qlogis(P); L[tot == 0] <- 0
  R <- L - outer(beta, beta, "-"); R[tot == 0] <- 0
  (R - t(R)) / 2                                     # enforce skew-symmetry
}

# real skew-symmetric R has eigenvalues in +/- i*lambda pairs; the positive
# lambda are the "bimension" strengths, each a plane of cyclic residual
# structure. Returns strengths (desc) and the leading plane's coordinates.
.btl_bimensions <- function(R) {
  e <- eigen(R)
  lam <- Im(e$values)
  keep <- which(lam > 1e-8)
  if (!length(keep)) return(list(strength = numeric(0), coord = NULL))
  ord <- keep[order(lam[keep], decreasing = TRUE)]
  v <- e$vectors[, ord[1]]
  list(strength = lam[ord], total = sum(R^2),
       coord = cbind(x = Re(v), y = Im(v)))
}

#' Residual dimensionality of paired comparisons
#'
#' The residual-PCA analogue for paired comparisons. The fitted model predicts
#' how often each object should beat each other from their locations; the
#' object-by-object matrix of departures from that prediction (on the
#' log-odds scale) is \emph{skew-symmetric}, so its structure decomposes into
#' rotational planes -- Gower's (1977) bimensions -- rather than the ordinary
#' components of a symmetric residual-correlation matrix. A dominant leading
#' bimension is a coherent \dQuote{swirl} in the residuals (A over-beats B, B
#' over-beats C, C over-beats A): a second attribute steering some contests. A
#' flat spectrum is noise: the single scale suffices. The leading bimension is
#' judged against a reference built by simulating unidimensional data from the
#' fitted model with the observed pair counts (a parametric bootstrap, as in
#' \code{\link{plot_scree}}); an observed strength above the reference is
#' structure the one-dimensional model does not explain. For graded fits the
#' residual log-odds are taken on the points-proportion scale, whose model
#' mean is not exactly \code{plogis(beta_i - beta_j)}; the simulated reference
#' carries the same construction, so the test stays calibrated (verified
#' mildly conservative on model-true graded data) rather than anticonservative.
#' Likewise, when the fit carries within-judge dependence effects
#' (\code{order}), the reference is simulated sequentially through each
#' judge's comparisons WITH the fitted exposure and carry-over coefficients:
#' order effects push the marginal pair rates around in a structured way, and
#' a reference without them would read that structure as a second attribute.
#' The price is power: carry-over and a judge-camp second attribute are
#' partially confounded (both appear as consistent within-judge deviation),
#' so with \code{order} modelled the test is conservative about attributing
#' the ambiguous share to a second dimension.
#' The reference simulates from the point estimates without refitting each
#' replicate, so it carries sampling noise in the responses but not
#' estimation noise in the parameters -- adequate for the screening use
#' here, slightly liberal in tiny designs.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}}.
#' @param reps Model-simulated replicates for the noise reference.
#' @return A list of class \code{"rasch_btl_dim"}: \code{bimensions} (per
#'   bimension: strength and share of residual size; the reference mean, 95th
#'   percentile, and the clears-the-reference flag are reported for the
#'   leading bimension and \code{NA} for the rest);
#'   \code{coords} (each object's position in the leading bimension plane, for
#'   the residual map); \code{leading_structured} (whether bimension 1 clears
#'   its reference); \code{residual_matrix}; and \code{notes}.
#' @references Gower, J. C. (1977). The analysis of asymmetry and orthogonality.
#'   In J. R. Barra et al. (Eds.), \emph{Recent Developments in Statistics}
#'   (pp. 109-123). North-Holland.
#' @examples
#' set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
#' pr <- t(utils::combn(objs, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' btl_dimensionality(btl(d, "a", "b", "win"), reps = 20)
#' @export
btl_dimensionality <- function(fit, reps = 50L) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  objs <- fit$objects$object; K <- length(objs); m <- fit$m
  if (K < 3L) stop("need at least three objects")
  beta <- setNames(fit$objects$location, objs)
  cmp <- fit$comparisons
  ia <- match(cmp$object_a, objs); ib <- match(cmp$object_b, objs)
  w <- cmp$weight
  notes <- character(0)

  R <- .btl_resid_matrix(ia, ib, cmp$response, w, m, K, beta)
  bm <- .btl_bimensions(R)
  if (!length(bm$strength)) stop("no residual structure to decompose")

  # parametric-bootstrap reference: simulate unidimensional data from the
  # fitted model at the observed pair counts and recompute the leading
  # strength. When the fit carries within-judge dependence effects the null
  # must carry them too -- order effects push the marginal pair rates away
  # from plogis(beta_i - beta_j) in a structured way, and a reference without
  # them reads that structure as a second attribute (false positives)
  d_lp <- beta[ia] - beta[ib]
  tau <- if (m > 1L) fit$thresholds$tau else NULL
  dep <- fit$dependence
  seq_sim <- NULL
  if (!is.null(dep)) {
    dd <- fit$dependence_data                       # sorted by judge, order
    sa <- match(dd$object_a, objs); sb <- match(dd$object_b, objs)
    sw <- dd$weight; sjd <- dd$judge
    coef_exp <- dep$estimate[match("exposure", dep$effect)]
    coef_cry <- dep$estimate[match("carry_over", dep$effect)]
    coef_pos <- dep$estimate[match("position", dep$effect)]
    if (is.na(coef_exp)) coef_exp <- 0
    if (is.na(coef_cry)) coef_cry <- 0
    if (is.na(coef_pos)) coef_pos <- 0
    # sequential simulation mirroring .btl_exposure's history rules, with
    # the FITTED coefficients: seen-before indicator and running mean verdict
    seq_sim <- function() {
      cnt <- new.env(hash = TRUE, parent = emptyenv())
      tot <- new.env(hash = TRUE, parent = emptyenv())
      gets <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
      resp <- integer(length(sa))
      for (r in seq_along(sa)) {
        ka <- paste0(sjd[r], "\r", sa[r]); kb <- paste0(sjd[r], "\r", sb[r])
        na_ <- gets(cnt, ka); nb_ <- gets(cnt, kb)
        z_exp <- as.numeric(na_ > 0) - as.numeric(nb_ > 0)
        z_cry <- (if (na_ > 0) gets(tot, ka) / na_ else 0) -
                 (if (nb_ > 0) gets(tot, kb) / nb_ else 0)
        # position is a constant +1 for every (object_a-first) row
        lp <- beta[sa[r]] - beta[sb[r]] + coef_exp * z_exp + coef_cry * z_cry +
              coef_pos
        x <- if (m == 1L) as.integer(stats::runif(1) < stats::plogis(lp))
             else sample.int(m + 1L, 1L, prob = item_moments(lp, tau)$P) - 1L
        resp[r] <- x
        cnt[[ka]] <- na_ + sw[r]; cnt[[kb]] <- nb_ + sw[r]
        tot[[ka]] <- gets(tot, ka) + sw[r] * (2 * x / m - 1)
        tot[[kb]] <- gets(tot, kb) + sw[r] * (2 * (m - x) / m - 1)
      }
      resp
    }
  }
  Pcat <- if (m == 1L || !is.null(seq_sim)) NULL else
    vapply(d_lp, function(dd) item_moments(dd, tau)$P, numeric(m + 1L))
  lead_ref <- vapply(seq_len(reps), function(r) {
    if (is.null(seq_sim)) {
      # a count-weighted row stands for w comparisons: simulate the
      # binomial (multinomial) SUM over those w and pass the mean response
      # with weight w, so the reference carries variance w p(1-p) per row,
      # not the w^2 p(1-p) of one weighted Bernoulli -- the overdispersed
      # reference made the test conservative on aggregated data
      resp <- if (m == 1L)
        stats::rbinom(length(d_lp), as.integer(w), stats::plogis(d_lp)) / w
      else vapply(seq_along(d_lp), function(r) {
        cnt <- stats::rmultinom(1L, as.integer(w[r]), Pcat[, r])
        sum((0:m) * cnt) / w[r]
      }, 0)
      Rr <- .btl_resid_matrix(ia, ib, resp, w, m, K, beta)
    } else {
      resp <- seq_sim()
      Rr <- .btl_resid_matrix(match(fit$dependence_data$object_a, objs),
                              match(fit$dependence_data$object_b, objs),
                              resp, fit$dependence_data$weight, m, K, beta)
    }
    s <- .btl_bimensions(Rr)$strength
    if (length(s)) s[1] else 0
  }, 0)
  ref_mean <- mean(lead_ref); ref_p95 <- stats::quantile(lead_ref, 0.95, names = FALSE)

  nb <- length(bm$strength)
  prop <- 2 * bm$strength^2 / bm$total
  bimensions <- data.frame(
    bimension = seq_len(nb), strength = bm$strength,
    prop_residual = prop,
    ref_mean = c(ref_mean, rep(NA_real_, nb - 1L)),
    ref_p95 = c(ref_p95, rep(NA_real_, nb - 1L)),
    above_reference = c(bm$strength[1] > ref_p95, rep(NA, nb - 1L)))

  coords <- data.frame(object = objs, location = unname(beta),
                       x = bm$coord[, "x"], y = bm$coord[, "y"])
  S <- .btl_scores(ia, ib, cmp$response, w, m, K)
  n_seen <- sum(((S + t(S)) > 0)[upper.tri(diag(K))])
  if (n_seen < choose(K, 2))
    notes <- c(notes, sprintf(
      "%d of %d pairs compared; unseen pairs contribute no residual",
      n_seen, choose(K, 2)))

  out <- list(bimensions = bimensions, coords = coords,
              leading_structured = isTRUE(bm$strength[1] > ref_p95),
              reference = list(mean = ref_mean, p95 = ref_p95, reps = reps),
              residual_matrix = R, notes = notes)
  class(out) <- "rasch_btl_dim"
  out
}

#' @export
print.rasch_btl_dim <- function(x, ...) {
  b <- x$bimensions
  cat(sprintf("Paired-comparison residual dimensionality: %d bimension(s)\n",
              nrow(b)))
  cat(sprintf("Leading bimension strength %.3f (%.0f%% of residual; reference 95%%: %.3f) -> %s\n",
              b$strength[1], 100 * b$prop_residual[1], x$reference$p95,
              if (x$leading_structured) "structured (a second attribute)"
              else "within noise (one scale suffices)"))
  for (n in x$notes) cat("Note:", n, "\n")
  invisible(x)
}

#' Consistency plot for paired-comparison transitivity
#'
#' With \code{by = "judge"} (the default when judges exist), plots each
#' judge's consistency -- one minus the circular-triad rate over the chance
#' rate -- as a dot against the chance line at zero: the individual-judge
#' lens, a judge-fit analogue. With \code{by = "object"}, plots each object's
#' circular-triad involvement instead: the structural lens, showing which
#' objects sit in the most contradictions.
#'
#' @param x A \code{"rasch_btl_transitivity"} object.
#' @param by \code{"auto"} (judges if present, else objects), \code{"judge"},
#'   or \code{"object"}.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @export
plot_btl_transitivity <- function(x, by = c("auto", "judge", "object"), ...) {
  stopifnot(inherits(x, "rasch_btl_transitivity"))
  by <- match.arg(by)
  use_judge <- if (by == "auto") !is.null(x$judges) else by == "judge"
  if (use_judge && is.null(x$judges))
    stop("no per-judge consistency (no judges, or too few compared triples)")
  if (use_judge) {
    j <- x$judges[order(x$judges$consistency), ]
    n <- nrow(j)
    op <- .rr_canvas(c(min(0, min(j$consistency)) - 0.02, 1), c(0.5, n + 0.5),
                     "Consistency  (1 = one clean order, 0 = guessing)", "",
                     yaxis = FALSE, grid_x = TRUE, grid_y = FALSE)
    on.exit(par(op))
    abline(v = 0, col = .rr$red, lty = 2, lwd = 1.5)
    abline(v = 1, col = .rr$soft, lty = 3)
    segments(0, seq_len(n), j$consistency, seq_len(n), col = .rr$grid, lwd = 3)
    points(j$consistency, seq_len(n), pch = 21, bg = .rr$blue, col = "white",
           cex = 1.5)
    text(par("usr")[1], seq_len(n), j$judge, pos = 4, cex = 0.8, col = .rr$ink,
         offset = 0.3)
    .rr_legend("bottomright", "chance", lwd = 1.5, lty = 2, col = .rr$red)
  } else {
    o <- x$objects
    n <- nrow(o)
    op <- .rr_canvas(c(0, max(o$circular_triads, 1)), c(0.5, n + 0.5),
                     "Circular triads the object sits in", "",
                     yaxis = FALSE, grid_x = TRUE, grid_y = FALSE)
    on.exit(par(op))
    segments(0, seq_len(n), rev(o$circular_triads), seq_len(n),
             col = .rr$grid, lwd = 3)
    points(rev(o$circular_triads), seq_len(n), pch = 21, bg = .rr$blue,
           col = "white", cex = 1.5)
    text(0, seq_len(n), rev(o$object), pos = 4, cex = 0.8, col = .rr$ink,
         offset = 0.3)
  }
}

#' Scree of paired-comparison residual bimensions
#'
#' Bimension strengths against the model-simulated noise reference (its mean
#' and 95th percentile band). A leading bar clearing the band is structured
#' residual dependence -- a likely second attribute.
#'
#' @param x A \code{"rasch_btl_dim"} object.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @export
plot_btl_scree <- function(x, ...) {
  stopifnot(inherits(x, "rasch_btl_dim"))
  b <- x$bimensions; k <- nrow(b)
  ref_m <- x$reference$mean; ref_p <- x$reference$p95
  ymax <- max(c(b$strength, ref_p)) * 1.15
  op <- .rr_canvas(c(0.5, k + 0.5), c(0, ymax), "Bimension", "Strength",
                   grid_x = FALSE)
  on.exit(par(op))
  rect(seq_len(k) - 0.32, 0, seq_len(k) + 0.32, b$strength,
       col = ifelse(c(x$leading_structured, rep(FALSE, k - 1)),
                    .rr$blue, .rr$soft), border = NA)
  # noise reference: mean line with a shaded band up to the 95th percentile
  rect(0.5, ref_m, k + 0.5, ref_p, col = "#dc262622", border = NA)
  abline(h = ref_m, col = .rr$red, lty = 5, lwd = 1.6)
  axis(1, at = seq_len(k), col = NA, col.ticks = NA)
  .rr_legend("topright", c("Observed", "Noise reference (mean, 95%)"),
             fill = c(.rr$blue, "#dc262633"), border = NA)
}

#' Residual map of the leading paired-comparison bimension
#'
#' Objects placed in the leading bimension plane. Reading round the swirl, an
#' object sits \dQuote{upstream} of those it over-beats relative to the fitted
#' locations; a clear rotational arrangement is the second attribute, a
#' formless blob near the origin is noise. Point size grows with the object's
#' location on the primary scale.
#'
#' @param x A \code{"rasch_btl_dim"} object.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @export
plot_btl_dim_map <- function(x, ...) {
  stopifnot(inherits(x, "rasch_btl_dim"))
  d <- x$coords
  r <- max(sqrt(d$x^2 + d$y^2), 1e-9)
  lim <- c(-r, r) * 1.25
  op <- .rr_canvas(lim, lim, "Leading bimension", "", grid_y = FALSE)
  on.exit(par(op))
  abline(h = 0, v = 0, col = .rr$grid, lwd = 0.8)
  symbols(0, 0, circles = r, add = TRUE, inches = FALSE, fg = .rr$grid,
          lwd = 0.8)
  loc <- d$location; cex <- 1.2 + 2 * (loc - min(loc)) / (max(loc) - min(loc) + 1e-9)
  points(d$x, d$y, pch = 21, bg = .rr$blue, col = "white", cex = cex)
  text(d$x, d$y, d$object, pos = 3, cex = 0.8, col = .rr$ink, offset = 0.5)
}

#' Unexpected judgements of one judge
#'
#' The paired-comparison counterpart of the kidmap. A judge has no ability to
#' condition on, so the reference is the consensus object scale (the pooled
#' locations). For the nominated judge, each object it met is given a
#' standardised residual oriented to the object -- how much more (\code{z > 0},
#' over-rated) or less (\code{z < 0}, under-rated) that judge favoured it than
#' its consensus location predicts. A surprise is an object the judge treated
#' against its standing: a strong object under-rated, or a weak object
#' over-rated (residual opposite in sign to the location), beyond
#' \code{flag_z} and seen at least \code{min_n} times.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}} with judges.
#' @param judge The judge to profile (a value of the fit's judge column).
#' @param min_n Objects met fewer times are shown but never flagged.
#' @param flag_z Absolute residual at or beyond which a contrary judgement is
#'   flagged unexpected.
#' @return A list of class \code{"rasch_btl_judge"}: \code{objects} (per object
#'   met: location, times met \code{n}, residual \code{z}, \code{surprise} flag
#'   and its \code{type}); \code{all_locations} (every object, for orientation);
#'   the \code{judge} and settings.
#' @examples
#' set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
#' pr <- t(utils::combn(objs, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
#' d$judge <- sample(paste0("J", 1:5), nrow(d), TRUE)
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' judge_surprise(btl(d, "a", "b", "win", judge = "judge"), "J1")
#' @export
judge_surprise <- function(fit, judge, min_n = 2L, flag_z = 1.96) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  cmp <- fit$comparisons
  if (all(is.na(cmp$judge))) stop("no judges in this fit")
  judge <- as.character(judge)
  sel <- which(!is.na(cmp$judge) & as.character(cmp$judge) == judge)
  if (!length(sel)) stop("no comparisons for judge ", judge)
  objs <- fit$objects$object; K <- length(objs); m <- fit$m
  beta <- setNames(fit$objects$location, objs)
  tau <- if (m > 1L) fit$thresholds$tau else NULL
  moments <- function(dd) if (m == 1L) {
    p <- stats::plogis(dd); list(E = p, V = p * (1 - p))
  } else { mo <- item_moments(dd, tau); list(E = mo$E, V = mo$V) }

  d <- cmp[sel, , drop = FALSE]
  ia <- match(d$object_a, objs); ib <- match(d$object_b, objs)
  obs <- exq <- vr <- nn <- numeric(K)
  for (r in seq_len(nrow(d))) {
    a <- ia[r]; b <- ib[r]; w <- d$weight[r]; x <- d$response[r]
    ma <- moments(beta[a] - beta[b])
    obs[a] <- obs[a] + w * x;         exq[a] <- exq[a] + w * ma$E
    vr[a]  <- vr[a]  + w * ma$V;      nn[a]  <- nn[a]  + w
    mb <- moments(beta[b] - beta[a])
    obs[b] <- obs[b] + w * (m - x);   exq[b] <- exq[b] + w * mb$E
    vr[b]  <- vr[b]  + w * mb$V;      nn[b]  <- nn[b]  + w
  }
  keep <- nn > 0
  z <- (obs - exq) / sqrt(pmax(vr, 1e-9))
  o <- data.frame(object = objs, location = unname(beta), n = nn,
                  z = z)[keep, , drop = FALSE]
  # a surprise: residual opposite in sign to the location (the judge pushed
  # the object against its consensus standing), large enough, and seen enough
  o$surprise <- abs(o$z) >= flag_z & o$n >= min_n & o$z * o$location < 0
  o$type <- ifelse(!o$surprise, "",
                   ifelse(o$location > 0, "strong object under-rated",
                          "weak object over-rated"))
  o <- o[order(-abs(o$z)), ]; rownames(o) <- NULL
  structure(list(judge = judge, objects = o, all_locations = beta,
                 n_comparisons = length(sel), flag_z = flag_z, min_n = min_n),
            class = "rasch_btl_judge")
}

#' @export
print.rasch_btl_judge <- function(x, ...) {
  cat(sprintf("Judge %s: %d comparisons over %d objects\n",
              x$judge, x$n_comparisons, nrow(x$objects)))
  s <- x$objects[x$objects$surprise, , drop = FALSE]
  if (nrow(s)) {
    cat("Unexpected judgements:\n")
    for (i in seq_len(nrow(s)))
      cat(sprintf("  %-6s (loc %+.2f): z = %+.2f  [%s]\n",
                  s$object[i], s$location[i], s$z[i], s$type[i]))
  } else cat("No object judged against its consensus standing.\n")
  invisible(x)
}

#' Unexpected judgements of one judge, pair by pair
#'
#' The comparison-level companion of \code{\link{judge_surprise}}. Each pair
#' the judge met is oriented to its stronger object (higher consensus
#' location) and given a standardised residual: \code{z < 0} means the
#' stronger object won less than its lead predicts -- the judge backed the
#' underdog. A matchup is an unexpected judgement when \code{z} falls at or
#' below \code{-flag_z} and the pair was seen at least \code{min_n} times, i.e.
#' the judge favoured the weaker object further than sampling noise explains.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}} with judges.
#' @param judge The judge to profile.
#' @param min_n Pairs met fewer times are shown but never flagged.
#' @param flag_z Absolute residual at or beyond which an upset is flagged.
#' @return A list of class \code{"rasch_btl_judge_pairs"}: \code{pairs} (per
#'   matchup: the stronger and weaker object and their locations, the location
#'   \code{gap}, times met \code{n}, residual \code{z}, the \code{net_winner},
#'   and the \code{surprise} flag); \code{all_locations}; the \code{judge} and
#'   settings.
#' @examples
#' set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
#' pr <- t(utils::combn(objs, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
#' d$judge <- sample(paste0("J", 1:5), nrow(d), TRUE)
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' judge_pair_surprise(btl(d, "a", "b", "win", judge = "judge"), "J1")
#' @export
judge_pair_surprise <- function(fit, judge, min_n = 1L, flag_z = 1.96) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  cmp <- fit$comparisons
  if (all(is.na(cmp$judge))) stop("no judges in this fit")
  judge <- as.character(judge)
  sel <- which(!is.na(cmp$judge) & as.character(cmp$judge) == judge)
  if (!length(sel)) stop("no comparisons for judge ", judge)
  objs <- fit$objects$object; K <- length(objs); m <- fit$m
  beta <- setNames(fit$objects$location, objs)
  tau <- if (m > 1L) fit$thresholds$tau else NULL
  d <- cmp[sel, , drop = FALSE]
  ia <- match(d$object_a, objs); ib <- match(d$object_b, objs)
  S <- .btl_scores(ia, ib, d$response, d$weight, m, K)   # points i scored on j
  N <- (S + t(S)) / m                                     # comparisons per pair
  rows <- list()
  for (i in seq_len(K - 1L)) for (j in (i + 1L):K) {
    if (N[i, j] <= 0) next
    hi <- if (beta[i] >= beta[j]) i else j; lo <- if (hi == i) j else i
    n <- N[hi, lo]; obs <- S[hi, lo]; dd <- beta[hi] - beta[lo]
    mo <- if (m == 1L) { p <- stats::plogis(dd); list(E = p, V = p * (1 - p)) }
          else { z <- item_moments(dd, tau); list(E = z$E, V = z$V) }
    zed <- (obs - n * mo$E) / sqrt(max(n * mo$V, 1e-9))
    rows[[length(rows) + 1L]] <- data.frame(
      object_hi = objs[hi], object_lo = objs[lo],
      loc_hi = unname(beta[hi]), loc_lo = unname(beta[lo]),
      gap = dd, n = n, z = zed,
      net_winner = if (obs >= n * m / 2) objs[hi] else objs[lo],
      surprise = zed <= -flag_z & n >= min_n)
  }
  p <- do.call(rbind, rows)
  p <- p[order(p$z), ]; rownames(p) <- NULL
  structure(list(judge = judge, pairs = p, all_locations = beta,
                 n_comparisons = length(sel), flag_z = flag_z, min_n = min_n),
            class = "rasch_btl_judge_pairs")
}

#' @export
print.rasch_btl_judge_pairs <- function(x, ...) {
  cat(sprintf("Judge %s: %d comparisons over %d matchups\n",
              x$judge, x$n_comparisons, nrow(x$pairs)))
  s <- x$pairs[x$pairs$surprise, , drop = FALSE]
  if (nrow(s)) {
    cat("Unexpected judgements (weaker object favoured beyond its lead):\n")
    for (i in seq_len(nrow(s)))
      cat(sprintf("  %s vs %s  (gap %.2f, z = %+.2f, %s)\n",
                  s$object_hi[i], s$object_lo[i], s$gap[i], s$z[i],
                  if (s$net_winner[i] == s$object_lo[i]) "upset"
                  else "favourite under-performed"))
  } else cat("No matchup went against the consensus beyond noise.\n")
  invisible(x)
}

#' Unexpected-judgement map for one judge (pair level)
#'
#' The judge counterpart of the kidmap, drawn matchup by matchup. Each pair the
#' judge met is a segment on the consensus location axis, spanning its two
#' objects, positioned horizontally by how surprising the verdict was: at zero
#' (the dashed line, inside the shaded band) the stronger object won as its
#' lead predicts; to the left the judge backed the underdog. A filled dot marks
#' the object the judge's verdict favoured, hollow the other -- so an upset is a
#' red segment on the left with its filled dot at the lower end. The rug marks
#' every object's location.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}} with judges.
#' @param judge The judge to map.
#' @param min_n,flag_z Passed to \code{\link{judge_pair_surprise}}.
#' @param ... Unused.
#' @return Called for its plotting side effect; invisibly the
#'   \code{rasch_btl_judge_pairs} object.
#' @export
plot_btl_judge_map <- function(fit, judge, min_n = 1L, flag_z = 1.96, ...) {
  jp <- judge_pair_surprise(fit, judge, min_n = min_n, flag_z = flag_z)
  p <- jp$pairs
  if (!nrow(p)) stop("this judge made no usable comparisons")
  xr <- max(abs(p$z), flag_z * 1.3)
  yr <- range(jp$all_locations)
  op <- .rr_canvas(c(-xr, xr) * 1.08,
                   yr + c(-1, 1) * (0.12 * diff(yr) + 0.2),
                   "Matchup residual   (backed the underdog  <-  0  ->  as expected)",
                   "Object location (logits)",
                   main = sprintf("Judge %s  \u00b7  %d matchups",
                                  jp$judge, nrow(p)))
  on.exit(par(op))
  u <- par("usr")
  rect(-flag_z, u[3], flag_z, u[4], col = "#94a3b81f", border = NA)  # expected
  abline(v = 0, col = .rr$soft, lwd = 1.2, lty = 2)
  rug(jp$all_locations, side = 2, col = .rr$grid, lwd = 1.4)
  # each matchup: a segment spanning its two objects at x = its residual, faint
  # when expected and bold red when the judge backed the underdog
  col <- ifelse(p$surprise, .rr$red, .rr$soft)
  lwd <- ifelse(p$surprise, 2.4, 1)
  segments(p$z, p$loc_lo, p$z, p$loc_hi, col = col, lwd = lwd)
  win_y <- ifelse(p$net_winner == p$object_hi, p$loc_hi, p$loc_lo)
  los_y <- ifelse(p$net_winner == p$object_hi, p$loc_lo, p$loc_hi)
  points(p$z, los_y, pch = 21, bg = "white", col = col, cex = 0.7)  # loser end
  points(p$z, win_y, pch = 21, bg = col, col = "white", cex = 1.2)  # winner end
  if (any(p$surprise)) {
    s <- p[p$surprise, ]
    text(s$z, (s$loc_hi + s$loc_lo) / 2,
         sprintf("%s-%s", s$object_hi, s$object_lo), pos = 2, offset = 0.4,
         cex = 0.7, col = .rr$red)
  }
  .rr_legend("bottomright", c("backed the underdog", "as expected"),
             lwd = c(2.4, 1), col = c(.rr$red, .rr$soft))
  invisible(jp)
}
