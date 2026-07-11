# rasch :: common-object equating for paired comparisons
# ===========================================================================
# The paired-comparison analogue of equate_tests(). Two Bradley-Terry-Luce
# calibrations that share a set of common objects -- the same scripts,
# performances, or products judged by different panels, or by the same panel
# in different years -- each fix their own origin by the sum-zero constraint.
# Because the two constraints are imposed over DIFFERENT object sets, the two
# origins do not coincide even when the objects are unchanged: each scale is
# centred on the mean of a different collection. Equating therefore estimates
# the scale shift between the origins (the precision-weighted mean difference
# over the common objects) and then tests each common object against the
# shifted identity line:
#
#   t_i = (b1_i - b2_i - shift) / sqrt(se1_i^2 + se2_i^2)
#
# Objects that survive define the equating link and carry the second panel's
# whole scale onto the first; objects that fail show drift -- a script the two
# panels valued differently, or a standard that moved between years -- and
# should be reviewed before the link is trusted. This is the standards-
# maintenance use of comparative judgement (Bramley 2007): a common set of
# anchor scripts lets panels judged apart be placed on one scale.
# ===========================================================================

# Coerce the second calibration to (object, location, se). It may be another
# btl fit or a "bank" -- a data frame of previously banked object locations,
# the paired-comparison counterpart of an item bank.
.btl_equate_ref <- function(reference) {
  if (inherits(reference, "rasch_btl"))
    return(data.frame(object = as.character(reference$objects$object),
                      location = reference$objects$location,
                      se = reference$objects$se,
                      stringsAsFactors = FALSE))
  if (is.data.frame(reference) || (is.list(reference) && !is.null(names(reference)))) {
    reference <- as.data.frame(reference, stringsAsFactors = FALSE)
    if (!all(c("object", "location") %in% names(reference)))
      stop("a bank needs columns 'object' and 'location' (and ideally 'se')")
    if (is.null(reference$se)) reference$se <- 0
    return(data.frame(object = as.character(reference$object),
                      location = as.numeric(reference$location),
                      se = as.numeric(reference$se),
                      stringsAsFactors = FALSE))
  }
  stop("`fit2` must be a btl fit or a bank data frame (object, location, se)")
}

#' Equate two paired-comparison calibrations through their common objects
#'
#' Compares the object locations of a Bradley-Terry-Luce fit with those of a
#' second fit (or a banked table of object locations), matched by object name.
#' The use case is standards maintenance and comparative judgement across
#' panels or years: the same scripts, performances, or products are judged by
#' different panels -- or by one panel in successive years -- and a common set
#' of anchor objects is carried through so that the two rounds land on a single
#' scale.
#'
#' Each calibration is identified by the sum-zero constraint, but the two
#' constraints are imposed over \emph{different} object sets, so the origins do
#' not coincide even when the shared objects are unchanged: each scale is
#' centred on the mean of a different collection. A scale shift between the two
#' origins is therefore estimated by the precision-weighted mean difference
#' over the common objects, and each common object is then tested against the
#' shifted identity line. A flagged object shows drift -- a script the two
#' panels valued differently, or a standard that moved between years -- and
#' weakens the equating link; the surviving objects carry the second
#' calibration's whole scale onto the first (\code{loc2 + shift}).
#'
#' The single shift presumes the drifting objects are a \emph{minority}. When
#' most of the common objects have genuinely moved, the precision-weighted
#' shift is pulled toward the movers and the drift tests can invert -- the
#' stable anchors flag as the apparent drifters. Read wholesale flagging
#' (several objects, one direction) as a contaminated link, not as evidence
#' about the individual objects; equate through a vetted anchor subset
#' instead. The \code{shift_se} accounts for the covariance of the location
#' estimates within each calibration (each is sum-zero constrained, so its
#' locations are not independent).
#'
#' @param fit1 A fitted object from \code{\link{btl}}: the calibration whose
#'   scale (origin) the equating targets.
#' @param fit2 A second \code{\link{btl}} fit, or a bank: a data frame with
#'   columns \code{object}, \code{location}, and optionally \code{se}.
#' @param alpha Significance level for the (multiplicity-adjusted) drift tests.
#' @param p_adjust Multiple-comparison adjustment across the common objects,
#'   passed to \code{\link[stats]{p.adjust}} (default \code{"holm"}).
#' @return A list of class \code{"rasch_btl_equate"}: the comparison
#'   \code{table} (per common object: object, both locations and standard
#'   errors, their \code{difference}, the \code{shifted_difference} against the
#'   estimated origin, the pooled \code{se_diff}, \code{t}, raw and adjusted
#'   \code{p}, and the \code{drifting} flag); the estimated \code{shift} and
#'   its \code{shift_se}; \code{equated}, the second calibration's full object
#'   table re-expressed on \code{fit1}'s scale; the number of common objects
#'   \code{n_common}; \code{alpha}; \code{p_adjust}; and \code{notes}.
#' @references Bramley, T. (2007). Paired comparison methods. In P. Newton,
#'   J. Baird, H. Goldstein, H. Patrick, & P. Tymms (Eds.), \emph{Techniques
#'   for monitoring the comparability of examination standards} (pp. 246-294).
#'   London: Qualifications and Curriculum Authority.
#' @examples
#' set.seed(1)
#' beta <- setNames(seq(-2, 2, length.out = 8), paste0("O", 1:8))
#' sim <- function(objs) {
#'   pr <- t(utils::combn(objs, 2))
#'   d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#'   d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#'   btl(d, "a", "b", "win")
#' }
#' eq <- btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)))
#' eq$table
#' @export
btl_equate <- function(fit1, fit2, alpha = 0.05, p_adjust = "holm") {
  if (!inherits(fit1, "rasch_btl"))
    stop("`fit1` must be a paired-comparison (btl) fit")
  cur <- data.frame(object = as.character(fit1$objects$object),
                    location = fit1$objects$location,
                    se = fit1$objects$se, stringsAsFactors = FALSE)
  ref <- .btl_equate_ref(fit2)
  common <- intersect(cur$object, ref$object)
  if (length(common) < 3)
    stop("need at least three common objects to equate paired-comparison scales")
  a <- cur[match(common, cur$object), ]
  b <- ref[match(common, ref$object), ]
  d <- a$location - b$location
  v <- a$se^2 + b$se^2
  w <- 1 / pmax(v, 1e-10)
  # precision-weighted mean difference: the shift between the two sum-zero
  # origins, best estimated where both calibrations are most certain
  c0 <- sum(w * d) / sum(w)
  # the common objects' location estimates are CORRELATED within each
  # sum-zero calibration, so Var(c0) = u' (Sigma1 + Sigma2) u with
  # u = w / sum(w), taken from the stored sandwich covariances (a bank
  # supplies no covariance and contributes diag(se^2) -- conservative)
  covsub <- function(fit, objs_c) {
    if (inherits(fit, "rasch_btl") && !is.null(fit$cov_beta)) {
      i <- match(objs_c, as.character(fit$objects$object))
      fit$cov_beta[i, i, drop = FALSE]
    } else NULL
  }
  u <- w / sum(w)
  S1 <- covsub(fit1, common)
  S2 <- covsub(fit2, common)
  if (is.null(S1)) S1 <- diag(a$se^2, length(common))
  if (is.null(S2)) S2 <- diag(b$se^2, length(common))
  shift_se <- sqrt(max(drop(t(u) %*% (S1 + S2) %*% u), 0))
  se_diff <- sqrt(pmax(v, 1e-10))
  t <- (d - c0) / se_diff
  p <- 2 * pnorm(-abs(t))
  p_adj <- p.adjust(p, method = p_adjust)
  drifting <- p_adj < alpha
  tab <- data.frame(object = common,
                    location_1 = a$location, se_1 = a$se,
                    location_2 = b$location, se_2 = b$se,
                    difference = d, shifted_difference = d - c0,
                    se_diff = se_diff, t = t, p = p, p_adj = p_adj,
                    drifting = drifting, stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  # the second calibration, whole, carried onto fit1's scale
  equated <- data.frame(object = ref$object,
                        location = ref$location + c0,
                        se = ref$se, stringsAsFactors = FALSE)
  rownames(equated) <- NULL
  notes <- sprintf(paste0("Origins differ because each calibration is sum-zero ",
                          "over its own object set; a shift of %.3f logits ",
                          "aligns fit2 to fit1."), c0)
  if (any(drifting))
    notes <- c(notes, sprintf(
      "%d common object(s) drift beyond the shifted link: %s",
      sum(drifting), paste(common[drifting], collapse = ", ")))
  structure(class = "rasch_btl_equate",
            list(table = tab, shift = c0, shift_se = shift_se,
                 equated = equated, n_common = length(common),
                 alpha = alpha, p_adjust = p_adjust, notes = notes))
}

#' Plot a paired-comparison equating comparison
#'
#' Scatter of the two calibrations' common-object locations with the shifted
#' identity line and per-object 95 per cent bands; objects that drift (after
#' the multiplicity adjustment) are highlighted and labelled. The counterpart
#' of \code{\link{plot_equate}} for Bradley-Terry-Luce scales.
#'
#' @param fit1 A fitted object from \code{\link{btl}}.
#' @param fit2 A second \code{\link{btl}} fit, or a bank data frame with columns
#'   \code{object}, \code{location}, and optionally \code{se}.
#' @param ... Passed to \code{\link{btl_equate}} (e.g. \code{alpha},
#'   \code{p_adjust}).
#' @return Called for its plotting side effect; invisibly the
#'   \code{\link{btl_equate}} result.
#' @examples
#' set.seed(1)
#' beta <- setNames(seq(-2, 2, length.out = 8), paste0("O", 1:8))
#' sim <- function(objs) {
#'   pr <- t(utils::combn(objs, 2))
#'   d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#'   d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#'   btl(d, "a", "b", "win")
#' }
#' plot_btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)))
#' @export
plot_btl_equate <- function(fit1, fit2, ...) {
  eq <- btl_equate(fit1, fit2, ...)
  tab <- eq$table
  rng <- range(c(tab$location_1, tab$location_2)) + c(-0.4, 0.4)
  op <- .rr_canvas(rng, rng, "Calibration 2 location (logits)",
                   "Calibration 1 location (logits)",
                   sprintf("%d common objects, shift %.3f, r = %.3f",
                           eq$n_common, eq$shift,
                           stats::cor(tab$location_1, tab$location_2)),
                   grid_x = TRUE)
  on.exit(par(op))
  abline(eq$shift, 1, col = .rr$ink, lwd = 2)
  band <- 1.96 * sqrt(mean(tab$se_1^2 + tab$se_2^2))
  abline(eq$shift + band, 1, lty = 3, col = .rr$soft)
  abline(eq$shift - band, 1, lty = 3, col = .rr$soft)
  segments(tab$location_2, tab$location_1 - 1.96 * tab$se_1,
           tab$location_2, tab$location_1 + 1.96 * tab$se_1,
           col = paste0(.rr$soft, "88"))
  points(tab$location_2, tab$location_1, pch = 21, cex = 1.6,
         bg = ifelse(tab$drifting, .rr$red, .rr$blue), col = "white", lwd = 1.2)
  if (any(tab$drifting))
    text(tab$location_2[tab$drifting], tab$location_1[tab$drifting],
         tab$object[tab$drifting], pos = 3, offset = 0.5, cex = 0.75,
         col = .rr$red)
  invisible(eq)
}

#' @export
print.rasch_btl_equate <- function(x, ...) {
  tab <- x$table
  cat(sprintf(paste0("Common-object equating over %d object(s): shift %.3f ",
                     "(SE %.3f), correlation %.3f, RMSD %.3f\n"),
              x$n_common, x$shift, x$shift_se,
              stats::cor(tab$location_1, tab$location_2),
              sqrt(mean(tab$shifted_difference^2))))
  core <- c("object", "location_1", "location_2", "shifted_difference", "t",
            "p_adj", "drifting")
  print(.fmt_df(tab[, intersect(core, names(tab))]), row.names = FALSE)
  cat(sprintf("%d object(s) drift beyond the %s-adjusted %.0f%% level.\n",
              sum(tab$drifting), x$p_adjust, 100 * (1 - x$alpha)))
  cat("(standard errors and unadjusted columns on $table; fit2 on fit1's scale in $equated)\n")
  invisible(x)
}
