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
#'   complete round-robin); \code{objects} (each object's circular-triad
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
#' structure the one-dimensional model does not explain.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}}.
#' @param reps Model-simulated replicates for the noise reference.
#' @return A list of class \code{"rasch_btl_dim"}: \code{bimensions} (per
#'   bimension: strength, share of residual size, the reference mean and 95th
#'   percentile, and whether the observed strength clears the reference);
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
  # fitted model at the observed pair counts and recompute the leading strength
  d_lp <- beta[ia] - beta[ib]
  Pcat <- if (m == 1L) NULL else {
    tau <- fit$thresholds$tau
    vapply(d_lp, function(dd) item_moments(dd, tau)$P, numeric(m + 1L))
  }
  lead_ref <- vapply(seq_len(reps), function(r) {
    resp <- if (m == 1L) stats::rbinom(length(d_lp), 1L, stats::plogis(d_lp))
            else apply(Pcat, 2, function(p) sample.int(m + 1L, 1L, prob = p) - 1L)
    Rr <- .btl_resid_matrix(ia, ib, resp, w, m, K, beta)
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
  n_seen <- sum((.btl_scores(ia, ib, cmp$response, w, m, K) +
                 t(.btl_scores(ia, ib, cmp$response, w, m, K)) > 0)[
                   upper.tri(diag(K))])
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
