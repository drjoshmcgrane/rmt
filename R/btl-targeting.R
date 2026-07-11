# rasch :: paired-comparison targeting and information
# ===========================================================================
# The test-information function of a Rasch calibration says where on the
# scale the instrument measures well; its paired-comparison analogue asks the
# same of a comparison design. The Fisher information one comparison carries
# about the location difference d = beta_a - beta_b is, in this exponential
# family, exactly the variance of its score: P(1 - P) for the dichotomous
# choice, and the graded response variance V for the ordinal extension (the
# score IS the sufficient statistic for d, so its variance IS the
# information). Summing per-comparison information over the comparisons a
# design actually contains gives a DESIGN information for every object -- the
# paired-comparison counterpart of an item's contribution to test
# information, and the quantity a naive (independent-trials) standard error
# would invert. Because information peaks at gap zero, close pairs are the
# informative ones: the same fact that makes an adaptive design spend its
# comparisons on near-neighbours (Pollitt 2012). The recommender here is that
# adaptive step, made honest about its two well-known hazards -- it is a
# greedy one-step heuristic, not an optimal design, and adaptive selection
# inflates a naively computed separation reliability (Bramley 2015).
# ===========================================================================

# per-comparison Fisher information about the location difference d: the
# variance of the (sufficient) score. Dichotomous P(1-P) is vectorised;
# the graded variance comes from item_moments, one difference at a time.
# The variance is symmetric in the orientation of the pair (symmetric
# thresholds make the graded model presentation-order invariant), so a
# comparison carries the same information about each of its two objects.
.btl_info_of_d <- function(d, m, tau) {
  if (m == 1L) {
    p <- stats::plogis(d)
    p * (1 - p)
  } else {
    vapply(d, function(z) item_moments(z, tau)$V, 0)
  }
}

#' Information and targeting of a paired-comparison design
#'
#' The paired-comparison analogue of the test-information function. The
#' Fisher information a single comparison carries about the location
#' difference \code{d = beta_a - beta_b} is, in this exponential family, the
#' variance of its score -- \code{P(1 - P)} for the dichotomous choice and
#' the graded response variance \code{V} for the ordinal extension (the
#' score is the sufficient statistic for \code{d}, so its variance is the
#' information). Weighted by each comparison's replication count and summed
#' over the comparisons the design actually contains, this gives a
#' \emph{design information} for every object: how much the observed
#' comparisons pin its location down, the counterpart of an item's
#' contribution to test information. Because the information peaks at gap
#' zero and falls away with the location gap, near-neighbour contests are
#' the informative ones.
#'
#' The design information inverts to a naive standard error
#' \code{se_naive = 1 / sqrt(information)}, the error one would report if the
#' comparisons were independent trials. The fit's own \code{se} is the
#' Godambe sandwich standard error, clustered by judge when judges are
#' identified; the two agree when comparisons are independent and diverge
#' when judges cluster (repeated verdicts by one judge carry less independent
#' information than their count suggests), so the gap between them is itself a
#' diagnostic of clustering.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}}.
#' @return A list of class \code{"rasch_btl_info"}: \code{objects} (per
#'   object: \code{location}, the fit's \code{se}, \code{n_comparisons},
#'   the design \code{information}, and \code{se_naive}); \code{pairs} (per
#'   observed pair: \code{n}, the mean location \code{gap}, and the pair's
#'   \code{information}); \code{comparisons} (per comparison: the signed
#'   \code{gap}, \code{weight}, and the single-comparison \code{information});
#'   the scalar \code{total} information; \code{m}; the \code{clustered}
#'   flag; and \code{notes}.
#' @references Pollitt, A. (2012). The method of adaptive comparative
#'   judgement. \emph{Assessment in Education}, 19(3), 281-300.
#' @seealso \code{\link{plot_btl_targeting}}, \code{\link{btl_next_pairs}}
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' btl_information(btl(d, "a", "b", "win"))
#' @export
btl_information <- function(fit) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  objs <- fit$objects$object
  K <- length(objs)
  beta <- setNames(fit$objects$location, objs)
  se <- setNames(fit$objects$se, objs)
  m <- fit$m
  tau <- if (m > 1L) fit$thresholds$tau else NULL
  cmp <- fit$comparisons
  ia <- match(cmp$object_a, objs); ib <- match(cmp$object_b, objs)
  w <- cmp$weight
  d <- unname(beta[ia] - beta[ib])            # signed location gap per comparison

  I_row <- .btl_info_of_d(d, m, tau)          # one-comparison information
  Iw <- w * I_row                             # weighted contribution to design

  # per-object design information: the sum of weighted per-comparison
  # information over every comparison the object took part in (as a or as b)
  obj_info <- numeric(K)
  for (idx in list(ia, ib)) {
    ag <- rowsum(Iw, idx)
    obj_info[as.integer(rownames(ag))] <-
      obj_info[as.integer(rownames(ag))] + ag[, 1]
  }

  objects <- data.frame(
    object = objs,
    location = unname(beta),
    se = unname(se),
    n_comparisons = fit$objects$comparisons,
    information = obj_info,
    se_naive = 1 / sqrt(obj_info),
    stringsAsFactors = FALSE)
  rownames(objects) <- NULL

  # per observed (unordered) pair: replications, the location gap (constant
  # within a pair), and the pair's pooled information
  lo <- pmin(ia, ib); hi <- pmax(ia, ib)
  key <- paste(lo, hi)
  n_pair <- tapply(w, key, sum)
  info_pair <- tapply(Iw, key, sum)
  gap_pair <- tapply(abs(d), key, mean)       # mean |gap| = |beta_a - beta_b|
  ix <- do.call(rbind, strsplit(names(n_pair), " "))
  pairs <- data.frame(
    object_a = objs[as.integer(ix[, 1])],
    object_b = objs[as.integer(ix[, 2])],
    n = as.numeric(n_pair),
    gap = as.numeric(gap_pair),
    information = as.numeric(info_pair),
    stringsAsFactors = FALSE)
  pairs <- pairs[order(-pairs$information), ]
  rownames(pairs) <- NULL

  comparisons <- data.frame(
    object_a = cmp$object_a, object_b = cmp$object_b,
    gap = d, weight = w, information = I_row,
    stringsAsFactors = FALSE)

  notes <- paste0(
    "se is the ", if (fit$clustered) "judge-clustered " else "",
    "Godambe sandwich standard error; se_naive = 1/sqrt(information) treats ",
    "comparisons as independent trials",
    if (fit$clustered)
      ", so the two differ when judges cluster their verdicts" else "")

  out <- list(objects = objects, pairs = pairs, comparisons = comparisons,
              total = sum(Iw), m = m, clustered = fit$clustered,
              notes = notes)
  class(out) <- "rasch_btl_info"
  out
}

#' @export
print.rasch_btl_info <- function(x, ...) {
  cat(sprintf(
    "Paired-comparison design information: %d objects, total %.2f\n",
    nrow(x$objects), x$total))
  cat(sprintf("One-comparison Fisher information about the location gap %s\n",
              if (x$m == 1L) "(dichotomous: P(1 - P))"
              else sprintf("(graded, %d categories: response variance)",
                           x$m + 1L)))
  print(.fmt_df(x$objects), row.names = FALSE)
  if (length(x$notes)) cat(sprintf("Note: %s\n", x$notes))
  invisible(x)
}

#' Targeting plot for a paired-comparison design
#'
#' The paired-comparison counterpart of a test-information display. Every
#' object is a dot at its location (x) and its design information (y, the
#' pooled Fisher information of the comparisons it took part in), the dot
#' sized by how many comparisons that is. A reference curve, read on the
#' right axis, traces the information a single \emph{new} comparison would
#' carry against an opponent at each location -- anchored at the centre of
#' the scale, so it peaks at gap zero and falls away with the gap. The curve
#' is the visual explanation of why an adaptive design chases near-neighbour
#' contests: information is bought most cheaply where the two objects are
#' close, and a well-targeted design lifts the low dots by pairing their
#' objects against opponents near them.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}}.
#' @param grid Optional location grid for the reference curve.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @seealso \code{\link{btl_information}}, \code{\link{btl_next_pairs}}
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' plot_btl_targeting(btl(d, "a", "b", "win"))
#' @export
plot_btl_targeting <- function(fit, grid = NULL) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  info <- btl_information(fit)
  o <- info$objects
  m <- fit$m
  tau <- if (m > 1L) fit$thresholds$tau else NULL

  # reference curve: information one new comparison carries against an
  # opponent at location x, anchored at the centre of the scale so gap =
  # anchor - x and the curve peaks at gap 0
  anchor <- stats::median(o$location)
  if (is.null(grid)) {
    rng <- range(o$location) + c(-1, 1)
    grid <- seq(rng[1], rng[2], length.out = 201)
  }
  refI <- .btl_info_of_d(anchor - grid, m, tau)

  ymax <- max(o$information) * 1.15
  op <- .rr_canvas(range(grid), c(0, ymax),
                   "Object location (logits)", "Design information",
                   right = 3.6)
  on.exit(par(op))

  # the reference curve carries its own (single-comparison) scale on axis 4
  scl <- ymax * 0.9 / max(refI)
  abline(v = anchor, col = .rr$soft, lty = 3)
  lines(grid, refI * scl, lwd = 2.2, col = .rr$red, lty = 5)
  ref_ticks <- pretty(c(0, max(refI)))
  ref_ticks <- ref_ticks[ref_ticks * scl <= ymax]
  axis(4, at = ref_ticks * scl, labels = ref_ticks, col = .rr$grid,
       col.ticks = .rr$soft, col.axis = .rr$red, cex.axis = 0.8)
  mtext("Information per comparison", side = 4, line = 2.3, col = .rr$red,
        cex = 0.85)

  # objects: dots sized by comparison count
  nc <- o$n_comparisons
  cex <- 1.1 + 2.2 * (nc - min(nc)) / (max(nc) - min(nc) + 1e-9)
  points(o$location, o$information, pch = 21, bg = .rr$blue, col = "white",
         cex = cex, lwd = 1.2)
  text(o$location, o$information, o$object, pos = 3, offset = 0.6, cex = 0.8,
       col = .rr$ink)

  .rr_legend("topright",
             c("Design information (dot size = comparisons)",
               "Information of one new comparison"),
             pch = c(21, NA), pt.bg = c(.rr$blue, NA), pt.cex = 1.3,
             lwd = c(NA, 2.2), lty = c(NA, 5), col = c("white", .rr$red))
  invisible(NULL)
}

#' Recommend the next informative comparisons (adaptive step)
#'
#' The adaptive comparative judgement step of Pollitt (2012): rank candidate
#' object pairs by the information one additional comparison would carry at
#' the current estimates. That information peaks when the two objects are
#' close in location, so at equal measurement the recommender favours
#' near-neighbour contests. With \code{weight_se = TRUE} (the default) the
#' expected information is multiplied by \code{se_a^2 + se_b^2}, so pairs of
#' poorly measured objects are promoted -- a greedy, one-step
#' D-optimality-flavoured heuristic that spends the next comparison where it
#' most reduces error.
#'
#' Two honest cautions. This is a \emph{greedy} rule that scores each pair on
#' its own immediate information; it is not a full optimal design and can be
#' beaten by a design that plans several comparisons jointly. And adaptive
#' selection is known to inflate a separation (scale) reliability computed
#' naively afterwards, because the design concentrates comparisons where they
#' shrink the errors most: report reliability from an independent or
#' non-adaptive subset, or treat an adaptive reliability as an upper bound
#' (Bramley 2015).
#'
#' @param fit A paired-comparison fit from \code{\link{btl}}.
#' @param n Number of pairs to return.
#' @param weight_se Weight the expected information by \code{se_a^2 + se_b^2}
#'   to prioritise poorly measured objects (default \code{TRUE}). When
#'   \code{FALSE}, pairs are ranked by expected information alone.
#' @return A data frame of the top \code{n} candidate pairs, each oriented to
#'   its stronger object: \code{object_a}, \code{object_b}, the location
#'   \code{gap}, \code{n_existing} (replications already observed for the
#'   pair), \code{expected_information} (of one new comparison), and
#'   \code{priority}. Sorted by \code{priority} (or by
#'   \code{expected_information} when \code{weight_se = FALSE}).
#' @references Pollitt, A. (2012). The method of adaptive comparative
#'   judgement. \emph{Assessment in Education}, 19(3), 281-300. Bramley, T.
#'   (2015). Investigating the reliability of Adaptive Comparative Judgment.
#'   \emph{Cambridge Assessment Research Report}.
#' @seealso \code{\link{btl_information}}, \code{\link{plot_btl_targeting}}
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' btl_next_pairs(btl(d, "a", "b", "win"), n = 5)
#' @export
btl_next_pairs <- function(fit, n = 10, weight_se = TRUE) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  objs <- fit$objects$object
  K <- length(objs)
  if (K < 2L) stop("need at least two objects to recommend a pair")
  beta <- setNames(fit$objects$location, objs)
  se <- setNames(fit$objects$se, objs)
  m <- fit$m
  tau <- if (m > 1L) fit$thresholds$tau else NULL

  # existing replications per unordered pair (keyed on object index)
  cmp <- fit$comparisons
  ia <- match(cmp$object_a, objs); ib <- match(cmp$object_b, objs)
  n_existing <- tapply(cmp$weight, paste(pmin(ia, ib), pmax(ia, ib)), sum)

  # every candidate pair, oriented to its stronger object
  cb <- utils::combn(K, 2)
  i <- cb[1, ]; j <- cb[2, ]
  hi <- ifelse(beta[i] >= beta[j], i, j)
  loo <- ifelse(hi == i, j, i)
  gap <- unname(beta[hi] - beta[loo])
  eI <- .btl_info_of_d(gap, m, tau)
  prio <- if (weight_se) eI * (se[hi]^2 + se[loo]^2) else eI
  ne <- n_existing[paste(pmin(i, j), pmax(i, j))]
  ne[is.na(ne)] <- 0

  df <- data.frame(
    object_a = objs[hi], object_b = objs[loo],
    gap = gap, n_existing = as.numeric(ne),
    expected_information = eI, priority = unname(prio),
    stringsAsFactors = FALSE)
  ord <- if (weight_se) order(-df$priority) else order(-df$expected_information)
  df <- df[ord, ]
  rownames(df) <- NULL
  utils::head(df, n)
}
