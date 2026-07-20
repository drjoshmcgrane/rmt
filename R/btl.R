# rasch :: Bradley-Terry-Luce paired comparisons
# ===========================================================================
# The Bradley-Terry-Luce model (Bradley & Terry 1952; Luce 1959) for paired
# comparisons: Pr{a beats b} = exp(beta_a - beta_b) / (1 + exp(...)). This
# is the conditional form of the dichotomous Rasch model (Rasch 1961;
# Andrich 1978): given that one of two items is answered correctly, the
# probability the harder one was is exactly of BTL form, and the package's
# pairwise conditional estimation maximises precisely such a likelihood on
# the pair-conditional counts. BTL therefore sits inside the same family
# with the person parameter replaced by exchangeable comparison
# replications, and it is estimated here by the same conventions as
# everything else: Newton-Raphson on the conditional likelihood, sum-zero
# identification, and Godambe sandwich standard errors, clustered by judge
# when judges are identified. Fit follows the same residual logic as the
# rest of the package (Andrich & Marais 2019, ch. 23): per
# comparison z = (y - P)/sqrt(PQ); objects and judges carry the
# log-of-mean-square fit residual over their comparisons with apportioned
# degrees of freedom, and the classical pairwise chi-square compares the
# observed and expected win proportions of every pair.
# ===========================================================================

# Exposure and carry-over covariates from each judge's own history: for
# comparison r, exposure is 1(judge saw object_a before) - 1(saw object_b
# before); carry-over differences the judge's mean prior verdicts on the
# two objects (oriented to each object, scaled to [-1, 1], zero when
# unseen). Both enter the exponent like the location difference does, so
# the dependence is measured in logits (Davidson & Beaver 1977 order-effect
# device; response-dependence logic of Marais & Andrich 2008).
.btl_exposure <- function(a, b, x, m, jd, ord, w = rep(1, length(a))) {
  R <- length(a)
  Fa <- Fb <- Wa <- Wb <- numeric(R)
  cnt <- new.env(hash = TRUE, parent = emptyenv())
  tot <- new.env(hash = TRUE, parent = emptyenv())
  gets <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
  for (r in order(jd, ord)) {
    ka <- paste0(jd[r], "\r", a[r]); kb <- paste0(jd[r], "\r", b[r])
    na_ <- gets(cnt, ka); nb_ <- gets(cnt, kb)
    Fa[r] <- as.numeric(na_ > 0); Fb[r] <- as.numeric(nb_ > 0)
    if (na_ > 0) Wa[r] <- gets(tot, ka) / na_
    if (nb_ > 0) Wb[r] <- gets(tot, kb) / nb_
    # count-weighted rows stand for w identical comparisons, so they enter
    # the history with weight w, consistent with the weighted likelihood
    cnt[[ka]] <- na_ + w[r]; cnt[[kb]] <- nb_ + w[r]
    tot[[ka]] <- gets(tot, ka) + w[r] * (2 * x[r] / m - 1)
    tot[[kb]] <- gets(tot, kb) + w[r] * (2 * (m - x[r]) / m - 1)
  }
  cbind(exposure = Fa - Fb, carry_over = Wa - Wb)
}

# union-find connectivity over the comparison graph
.btl_components <- function(K, ia, ib) {
  parent <- seq_len(K)
  find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
  for (r in seq_along(ia)) {
    ra <- find(ia[r]); rb <- find(ib[r])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(K), find, 1L)
}

#' Fit the Bradley-Terry-Luce model to paired comparisons
#'
#' Estimates object locations from paired-comparison data by conditional
#' maximum likelihood. The Bradley-Terry-Luce model is the conditional form
#' of the dichotomous Rasch model -- within an item pair, given one correct
#' response, the Rasch probability that it was the easier item is exactly
#' of BTL form -- so it belongs to the same measurement family and is
#' estimated by the same conventions as the rest of the package:
#' Newton-Raphson on the person-free likelihood, locations identified by
#' the sum-zero constraint, and Godambe sandwich standard errors, clustered
#' by judge when a judge column is given (so repeated comparisons by the
#' same judge need not be independent). Objects that win or lose every
#' comparison have no finite estimate and are removed with a note, exactly
#' as extreme persons are set aside in a Rasch calibration; the comparison
#' graph must remain connected.
#'
#' Fit is reported at three levels, mirroring the Rasch diagnostics.
#' Per object and (when given) per judge: the log-of-mean-square fit
#' residual of Andrich and Marais (2019, ch. 23) over their comparisons, with apportioned degrees of freedom --
#' an erratic judge or an object of inconsistent quality shows exactly as
#' an erratic person or misfitting item does. Per pair: the classical
#' goodness-of-fit table comparing observed and expected win proportions,
#' with the total chi-square on (pairs used) minus (free location parameters)
#' degrees of freedom. The object separation index is the analogue of the
#' PSI: the proportion of observed location variance not due to error.
#' Anchored objects enter it with their fixed locations and zero error --
#' they are separated with certainty by construction.
#'
#' @param data A data frame with one comparison per row.
#' @param object_a,object_b Names of the columns holding the two objects
#'   compared.
#' @param winner Name of the column holding the winner of each row: its
#'   value must equal the row's \code{object_a} or \code{object_b} entry;
#'   \code{"tie"} or \code{"draw"} marks a tie; anything else (including
#'   blanks) is treated as missing and dropped with a note. Ignored when
#'   \code{response} is given.
#' @param margin Optional name of a column holding the extent of the win
#'   ("a little", "much", ...), as an ordered factor or increasing values;
#'   combined with \code{winner} it assembles the graded response without
#'   any orientation bookkeeping ("B by much" means the same thing
#'   whichever column B sits in). Winner values matching neither object
#'   are ties and form the middle category.
#' @param thresholds \code{"free"} (default) estimates every symmetric
#'   threshold parameter; \code{"pc"} pools them to the spread (linear)
#'   principal component -- the symmetric case of the principal-component
#'   threshold structure, whose even skewness component is structurally
#'   zero here -- so thinly used categories borrow strength from every
#'   response. Both modes report the component decomposition.
#' @param response Optional name of a column holding a graded preference
#'   for \code{object_a} over \code{object_b} -- an ORDERED factor
#'   (\code{factor(..., ordered = TRUE)}, levels worst to best for
#'   \code{object_a}) or integer scores \code{0..m}; a plain factor is
#'   refused, since its alphabetical level order would silently define
#'   (and can reverse) the response scale. Fits the
#'   adjacent-categories ordinal extension of BTL (Tutz 1986; Agresti
#'   1992): a partial-credit structure on the difference of locations with
#'   thresholds constrained symmetric, \code{tau_k = -tau_(m+1-k)}, so the
#'   model is invariant to presentation order. Two categories reproduce
#'   BTL exactly; three give the Davidson (1970) ties model.
#' @param judge Optional name of a judge column; enables the judge fit
#'   table and clusters the sandwich standard errors by judge.
#' @param order Optional name of a column giving each judge's judgment
#'   sequence (timestamps or ranks; requires \code{judge}). Adds the
#'   within-judge dependence analysis: an exposure effect (the advantage,
#'   in logits, of an object the judge has seen before over one they have
#'   not) and a carry-over effect (the pull of the judge's own earlier
#'   verdicts on the same object -- response dependence in the sense of
#'   Marais and Andrich 2008), estimated jointly with the locations and
#'   reported in \code{dependence}. Incompatible with
#'   \code{ties = "half"}.
#' @param position Logical: when \code{TRUE}, \code{object_a} is taken as the
#'   first-presented (left) object of every comparison and a first-position
#'   advantage is estimated -- a single coefficient, in logits, added to
#'   every comparison's location difference, the pure positional form of the
#'   Davidson and Beaver (1977) within-pair order-effect device. It is
#'   reported in \code{dependence} with \code{effect = "position"} (every
#'   comparison is informative, so \code{n_informative} is the total weighted
#'   comparison count) and estimated jointly with the locations, alongside the
#'   exposure and carry-over effects when \code{order} is also given.
#'   Identification comes from triangle closure (K >= 3), so the constant
#'   oriented covariate is estimable even when each pair has a fixed
#'   orientation, though weakly. Note that \code{ties = "half"} duplicates
#'   rows in the same orientation, so the first position stays well defined.
#' @param anchors Optional named numeric vector for equating: names are object
#'   names, values are fixed locations in logits. The named objects are held
#'   exactly at those locations and the remaining objects are estimated freely
#'   with no sum-zero constraint -- the origin and scale come from the anchors,
#'   exactly as an anchored \code{\link{rasch}} calibration works. Anchored
#'   objects report a standard error of zero (their location is a constant, not
#'   an estimate). An anchored object that is undefeated or winless is an error,
#'   not silently removed as a free boundary object would be.
#' @param count Optional name of a column of replication counts (a row
#'   standing for several identical comparisons).
#' @param ties How to treat ties in the dichotomous analysis:
#'   \code{"drop"} (default, removed with a note), \code{"half"} (half a
#'   win each way, a common pragmatic device -- flagged in the notes
#'   because the halves are not independent Bernoulli trials), or
#'   \code{"error"}. With graded responses, code ties as a middle
#'   category instead.
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return A list of class \code{"rasch_btl"}: \code{objects} (location, se,
#'   comparisons, wins -- or the graded \code{score} -- infit and outfit
#'   mean squares, fit residual and its df),
#'   \code{pairs} (per pair: n, observed and expected win proportions --
#'   or mean graded responses --
#'   standardised residual, chi-square component -- the pair chi-squares
#'   treat comparisons as independent and are descriptive under judge
#'   clustering; the object and judge fit residuals and the clustered
#'   standard errors carry the robust inference), \code{judges} (when
#'   given: per judge n, infit, outfit, fit residual, df), \code{total_chisq},
#'   \code{total_df}, \code{total_p}, the object separation index
#'   \code{osi}, \code{loglik}, \code{cl} (the composite-likelihood
#'   information ingredients used by \code{\link{compare_fits}}: the Godambe
#'   effective parameter count and the independent-unit count),
#'   convergence details, and \code{notes}.
#'   Graded fits add \code{thresholds} (the symmetric threshold estimates
#'   with standard errors), \code{m}, and \code{categories}. With an
#'   \code{order} column the within-judge \code{dependence} effects table
#'   carries an \code{n_informative} count, and \code{dependence_data} holds
#'   every comparison with its per-comparison exposure and carry-over
#'   covariates (see \code{\link{plot_btl_dependence}}).
#' @references Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#'   incomplete block designs: I. The method of paired comparisons.
#'   Biometrika, 39, 324-345. Luce, R. D. (1959). Individual Choice
#'   Behavior. Wiley. Andrich, D. (1978). Relationships between the
#'   Thurstone and Rasch approaches to item scaling. Applied Psychological
#'   Measurement, 2, 451-462.
#'
#'   Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response.
#'   Journal of Mathematical Psychology, 30(3), 306-316. Agresti, A.
#'   (1992). Analysis of ordinal paired comparison data. Journal of the
#'   Royal Statistical Society C, 41(2), 287-297. Davidson, R. R. (1970).
#'   On extending the Bradley-Terry model to accommodate ties in paired
#'   comparison experiments. Journal of the American Statistical
#'   Association, 65(329), 317-328.
#'
#'   Davidson, R. R., & Beaver, R. J. (1977). On extending the Bradley-Terry
#'   model to incorporate within-pair order effects. Biometrics, 33(4),
#'   693-702.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pairs <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pairs[, 1], each = 30),
#'                 b = rep(pairs[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' btl(d, object_a = "a", object_b = "b", winner = "win")
#' @export
btl <- function(data, object_a, object_b, winner = NULL, response = NULL,
                margin = NULL, judge = NULL, count = NULL, order = NULL,
                position = FALSE, anchors = NULL,
                ties = c("drop", "half", "error"),
                thresholds = c("free", "pc"), maxit = 60, tol = 1e-8) {
  ties <- match.arg(ties)
  thresholds <- match.arg(thresholds)
  data <- as.data.frame(data)
  if (is.null(winner) && is.null(response))
    stop("give either `winner` (dichotomous) or `response` (graded)")
  if (!is.null(margin) && is.null(winner))
    stop("`margin` requires `winner`")
  if (!is.null(order) && is.null(judge))
    stop("`order` requires `judge`: exposure is a within-judge history")
  for (col in c(object_a, object_b, winner, response, margin, judge, count,
                order))
    if (!col %in% names(data)) stop("column not found: ", col)
  a <- trimws(as.character(data[[object_a]]))
  b <- trimws(as.character(data[[object_b]]))
  jd <- if (is.null(judge)) NULL else as.character(data[[judge]])
  w <- if (is.null(count)) rep(1, nrow(data)) else as.numeric(data[[count]])
  ord <- if (is.null(order)) NULL else as.numeric(data[[order]])
  notes <- character(0)
  if (!is.null(anchors)) {
    if (!is.numeric(anchors) || is.null(names(anchors)) ||
        any(!nzchar(names(anchors))))
      stop("`anchors` must be a named numeric vector (names = object names)")
    if (sum(names(anchors) %in% unique(c(a, b))) < 1L)
      stop("no `anchors` name matches an object in the data")
  }
  # a constant, object_a-oriented covariate for the first-position advantage;
  # appended to the dependence design (alone, or beside exposure/carry-over)
  add_pos <- function(Z, n)
    if (isTRUE(position)) cbind(Z, position = rep(1, n)) else Z

  if (!is.null(response)) {
    xr <- data[[response]]
    if (is.factor(xr)) {
      # a plain factor's alphabetical level order would silently define the
      # response scale (and can reverse it); the order must be explicit
      if (!is.ordered(xr))
        stop("a graded response factor must be ORDERED ",
             "(factor(..., ordered = TRUE) with levels from worst to ",
             "best), or supply integer scores 0..m: an unordered ",
             "factor's alphabetical levels would silently define -- and ",
             "can reverse -- the response scale")
      cats <- levels(xr); x <- as.integer(xr) - 1L
    } else {
      .check_integer_scores(xr, "the graded response")
      xn <- suppressWarnings(as.numeric(as.character(xr)))
      x <- as.integer(xn)
      if (any(x < 0, na.rm = TRUE))
        stop("graded responses must be non-negative integers 0..m")
      cats <- as.character(0:max(x, na.rm = TRUE))
    }
    keep <- !is.na(a) & !is.na(b) & !is.na(x) & a != b & !is.na(w) & w > 0
    if (!is.null(jd)) keep <- keep & !is.na(jd)
    if (!is.null(ord)) keep <- keep & !is.na(ord)
    if (any(!keep)) {
      notes <- c(notes, sprintf(
        "%d row(s) dropped (missing, zero-count, or self-comparison)",
        sum(!keep)))
      a <- a[keep]; b <- b[keep]; x <- x[keep]; w <- w[keep]
      if (!is.null(jd)) jd <- jd[keep]
      if (!is.null(ord)) ord <- ord[keep]
    }
    if (!length(a)) stop("no usable comparisons")
    Z <- if (is.null(ord)) NULL else
      .btl_exposure(a, b, x, length(cats) - 1L, jd, ord, w)
    Z <- add_pos(Z, length(a))
    return(.btl_graded(a, b, x, jd, w, cats, maxit, tol, notes,
                       thr = thresholds, Z = Z, ord = ord, anchors = anchors))
  }

  if (!is.null(margin)) {
    # winner + margin entry: orientation-free by construction. The graded
    # response is assembled from "who won" and "by how much"; a winner value
    # matching neither object is a tie and becomes the middle category.
    mg <- data[[margin]]
    # the margin's level order defines the graded scale, so it must be
    # explicit: an ordered factor (smallest to largest margin) or a numeric
    # magnitude. A plain factor's -- or a character column's -- alphabetical
    # order can silently reverse which margin counts as the big win.
    if (is.factor(mg) && !is.ordered(mg))
      stop("`margin` must be an ORDERED factor (factor(..., ordered = ",
           "TRUE), levels smallest to largest margin) or a numeric ",
           "magnitude: a plain factor's alphabetical level order would ",
           "silently define -- and can reverse -- the margin scale")
    if (is.character(mg))
      stop("`margin` is a character column; supply an ORDERED factor ",
           "(levels smallest to largest margin) or a numeric magnitude, ",
           "so the margin order is explicit rather than alphabetical")
    lv <- if (is.factor(mg)) levels(droplevels(mg)) else
      as.character(sort(unique(mg[!is.na(mg)])))
    q <- length(lv)
    if (q < 1L) stop("`margin` has no usable levels")
    mgi <- match(as.character(mg), lv)
    wn <- trimws(as.character(data[[winner]]))
    is_a <- !is.na(wn) & wn == a
    is_b <- !is.na(wn) & wn == b
    tie <- !is.na(wn) & !is_a & !is_b & tolower(wn) %in% c("tie", "draw")
    miss_wn <- !is.na(wn) & !is_a & !is_b & !tie
    if (any(miss_wn))
      notes <- c(notes, sprintf(
        "%d row(s) with winner matching neither object treated as missing",
        sum(miss_wn)))
    ties_present <- any(tie)
    keep <- !is.na(a) & !is.na(b) & !is.na(wn) & a != b & !is.na(w) & w > 0 &
      (tie | !is.na(mgi)) & !miss_wn
    if (!is.null(jd)) keep <- keep & !is.na(jd)
    if (!is.null(ord)) keep <- keep & !is.na(ord)
    if (any(!keep)) {
      notes <- c(notes, sprintf(
        "%d row(s) dropped (missing winner or margin, zero-count, or self-comparison)",
        sum(!keep)))
      a <- a[keep]; b <- b[keep]; w <- w[keep]
      mgi <- mgi[keep]; is_a <- is_a[keep]; is_b <- is_b[keep]
      tie <- tie[keep]
      if (!is.null(jd)) jd <- jd[keep]
      if (!is.null(ord)) ord <- ord[keep]
    }
    if (!length(a)) stop("no usable comparisons")
    base <- q - 1L + as.integer(ties_present)
    x <- ifelse(is_a, base + mgi, ifelse(is_b, q - mgi, q))
    cats <- c(paste0("worse by ", rev(lv)), if (ties_present) "tie",
              paste0("better by ", lv))
    if (ties_present)
      notes <- c(notes, sprintf("%d tie(s) placed in the middle category",
                                sum(tie)))
    Z <- if (is.null(ord)) NULL else
      .btl_exposure(a, b, as.integer(x), length(cats) - 1L, jd, ord, w)
    Z <- add_pos(Z, length(a))
    return(.btl_graded(a, b, as.integer(x), jd, w, cats, maxit, tol, notes,
                       thr = thresholds, Z = Z, ord = ord, anchors = anchors))
  }

  wn <- trimws(as.character(data[[winner]]))
  keep <- !is.na(a) & !is.na(b) & !is.na(wn) & a != b & !is.na(w) & w > 0
  if (!is.null(jd)) keep <- keep & !is.na(jd)
  if (!is.null(ord)) keep <- keep & !is.na(ord)
  if (any(!keep)) {
    notes <- c(notes, sprintf("%d row(s) dropped (missing, zero-count, or self-comparison)",
                              sum(!keep)))
    a <- a[keep]; b <- b[keep]; wn <- wn[keep]; w <- w[keep]
    if (!is.null(jd)) jd <- jd[keep]
    if (!is.null(ord)) ord <- ord[keep]
  }
  if (!length(a)) stop("no usable comparisons")
  if (!is.null(ord) && ties == "half")
    stop("exposure analysis is incompatible with ties = 'half';",
         " drop ties or code them as a graded middle category")

  # outcome: 1 = a wins, 0 = b wins; an explicit "tie"/"draw" entry is a
  # tie; anything else matching neither object is missing, not a tie
  y <- ifelse(wn == a, 1, ifelse(wn == b, 0, NA))
  is_tie <- is.na(y) & tolower(wn) %in% c("tie", "draw")
  miss <- is.na(y) & !is_tie
  if (any(miss)) {
    notes <- c(notes, sprintf(
      "%d row(s) with winner matching neither object treated as missing",
      sum(miss)))
    sel <- !miss
    a <- a[sel]; b <- b[sel]; y <- y[sel]; w <- w[sel]
    if (!is.null(jd)) jd <- jd[sel]
    if (!is.null(ord)) ord <- ord[sel]
    if (!length(a)) stop("no usable comparisons")
  }
  if (anyNA(y)) {
    n_tie <- sum(is.na(y))
    if (ties == "error") stop(n_tie, " tie(s) present; set ties = 'drop' or 'half'")
    if (ties == "drop") {
      notes <- c(notes, sprintf("%d tie(s) dropped", n_tie))
      sel <- !is.na(y)
      a <- a[sel]; b <- b[sel]; y <- y[sel]; w <- w[sel]
      if (!is.null(jd)) jd <- jd[sel]
      if (!is.null(ord)) ord <- ord[sel]
    } else {
      notes <- c(notes, sprintf("%d tie(s) scored half a win each way (halves are not independent trials)",
                                n_tie))
      t_i <- which(is.na(y))
      a <- c(a, a[t_i]); b <- c(b, b[t_i])
      y[t_i] <- 1; y <- c(y, rep(0, length(t_i)))
      w[t_i] <- w[t_i] / 2; w <- c(w, w[t_i])
      if (!is.null(jd)) jd <- c(jd, jd[t_i])
    }
  }

  if (!is.null(ord) || isTRUE(position)) {
    # exposure/position covariates route through the graded engine, whose two-
    # category case reproduces the dichotomous analysis exactly
    Z <- if (is.null(ord)) NULL else
      .btl_exposure(a, b, as.integer(y), 1L, jd, ord, w)
    Z <- add_pos(Z, length(a))
    return(.btl_graded(a, b, as.integer(y), jd, w, c("0", "1"), maxit, tol,
                       notes, thr = "free", Z = Z, ord = ord,
                       anchors = anchors))
  }

  # the two-category graded engine IS the dichotomous conditional model
  # (their equivalence is tested to machine precision), so one estimator
  # serves both routes; m == 1 results are presented as wins / win
  # proportions inside .btl_graded
  .btl_graded(a, b, as.integer(y), jd, w, c("0", "1"), maxit, tol,
              notes, thr = "free", anchors = anchors)
}

#' @export
print.rasch_btl <- function(x, ...) {
  cat(sprintf("Bradley-Terry-Luce analysis: %d objects, %.0f comparisons%s\n",
              nrow(x$objects), x$n_comparisons,
              if (!is.null(x$judges)) sprintf(", %d judges", nrow(x$judges)) else ""))
  cat(sprintf("Conditional ML: %s in %d iterations; sandwich SEs%s\n",
              if (x$converged) "converged" else "NOT converged", x$iterations,
              if (x$clustered) " clustered by judge" else ""))
  cat(sprintf("Object separation index %.3f; pairwise chi-square %.2f on %d df, p = %s\n",
              x$osi$PSI, x$total_chisq, x$total_df, .fmt_p(x$total_p)))
  if (!is.null(x$anchors))
    cat(sprintf("Anchored at %d object(s) (se = 0): %s\n",
                length(x$anchors), paste(names(x$anchors), collapse = ", ")))
  if (!is.null(x$dependence)) {
    for (r in seq_len(nrow(x$dependence))) {
      # position is a static first-presented advantage, not a within-judge
      # history effect, so it is not labelled "Within-judge"
      lab <- if (x$dependence$effect[r] == "position")
        "First-position advantage" else
        paste("Within-judge", gsub("_", "-", x$dependence$effect[r]))
      # three saved-fit schemas: `t` (current, t reference); `z` with `df`
      # (1.11.4 transitional: t-based inference under the old name); `z`
      # without `df` (older fits whose p-values used the NORMAL reference:
      # keep their own label rather than misrepresent their inference)
      if (!is.null(x$dependence$t)) {
        st_lab <- "t"; st_r <- x$dependence$t[r]
      } else if (!is.null(x$dependence$df)) {
        st_lab <- "t"; st_r <- x$dependence$z[r]
      } else {
        st_lab <- "z"; st_r <- x$dependence$z[r]
      }
      cat(sprintf("%s: %.3f logits (SE %.3f, %s = %.2f, p = %s)\n",
                  lab, x$dependence$estimate[r], x$dependence$se[r],
                  st_lab, st_r, .fmt_p(x$dependence$p[r])))
    }
  }
  if (!is.null(x$thresholds)) {
    cat(sprintf("Graded comparisons in %d categories%s; symmetric thresholds: %s\n",
                x$m + 1L,
                if (!is.null(x$categories) &&
                    !all(x$categories == as.character(0:x$m)))
                  paste0(" (", paste(x$categories, collapse = " < "), ")")
                else "",
                paste(sprintf("%.3f", x$thresholds$tau), collapse = ", ")))
  }
  print(.fmt_df(x$objects[, intersect(c("object", "location", "se",
                                        "comparisons", "wins", "score",
                                        "fit_resid"),
                                      names(x$objects))]), row.names = FALSE)
  if (!is.null(x$judges)) {
    mis <- x$judges[!is.na(x$judges$fit_resid) & abs(x$judges$fit_resid) > 2.5, ]
    cat(sprintf("Judges beyond |fit residual| 2.5: %d%s\n", nrow(mis),
                if (nrow(mis)) paste0(" (", paste(mis$judge, collapse = ", "), ")")
                else ""))
  }
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' Plot Bradley-Terry-Luce object locations
#'
#' Caterpillar plot of the object locations with 95 per cent error bars,
#' misfitting objects highlighted, in the package's house style.
#'
#' @param fit An object from \code{\link{btl}}.
#' @param band Absolute fit-residual value beyond which an object is
#'   highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pairs <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pairs[, 1], each = 30),
#'                 b = rep(pairs[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' plot_btl(btl(d, "a", "b", "win"))
#' @export
plot_btl <- function(fit, band = 2.5) {
  d <- fit$objects[order(fit$objects$location), ]
  k <- nrow(d)
  xlim <- range(c(d$location - 1.96 * d$se, d$location + 1.96 * d$se))
  op <- .rr_canvas(xlim + c(-0.15, 0.15) * diff(xlim), c(0.5, k + 0.5),
                   "Location (logits)", "", grid_y = FALSE, grid_x = TRUE,
                   yaxis = FALSE)
  on.exit(par(op))
  mis <- !is.na(d$fit_resid) & abs(d$fit_resid) > band
  segments(d$location - 1.96 * d$se, seq_len(k),
           d$location + 1.96 * d$se, seq_len(k),
           col = ifelse(mis, .rr$red, .rr$soft), lwd = 2.2)
  points(d$location, seq_len(k), pch = 21, cex = 1.6, lwd = 1.2,
         bg = ifelse(mis, .rr$red, .rr$blue), col = "white")
  text(d$location, seq_len(k), d$object, pos = 3, offset = 0.55, cex = 0.8,
       col = .rr$ink)
  if (any(mis))
    .rr_legend("bottomright", sprintf("|fit residual| > %.1f", band),
               pch = 21, pt.bg = .rr$red, col = "white", pt.cex = 1.4)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Graded paired comparisons: the adjacent-categories (Rasch-type) ordinal
# extension of BTL (Tutz 1986; Agresti 1992). The response is one of m+1
# ordered categories from object_a's perspective ("much worse" ... "much
# better"); category probabilities follow a partial-credit structure on the
# difference beta_a - beta_b with thresholds constrained symmetric,
# tau_k = -tau_{m+1-k}, so the model is invariant to presentation order and
# judge tendencies cancel. m = 1 is exactly BTL; m = 2 is the Davidson
# (1970) ties model. Estimation, identification, sandwich errors, and fit
# follow the package conventions established in btl().
# ---------------------------------------------------------------------------
.btl_graded <- function(a, b, x, jd, w, cats, maxit, tol, notes,
                        thr = "free", Z = NULL, ord = NULL, anchors = NULL) {
  m <- length(cats) - 1L
  if (m < 1L) stop("graded responses need at least two categories")
  # identifiability: empty EXTREME categories leave no finite spread (the
  # data are evidence of infinite spread, as a zero raw score is of an
  # infinite person location); empty interior categories are unidentified
  # under free thresholds but pooled over by the principal-component
  # structure
  check_cats <- function(x, note_interior) {
    xs <- c(x, m - x)
    emp <- which(tabulate(xs + 1L, m + 1L) == 0) - 1L
    if (any(emp %in% c(0L, m)))
      stop("extreme category never used (in either orientation): ",
           paste(cats[intersect(emp, c(0L, m)) + 1L], collapse = ", "),
           "; no finite threshold estimate exists - collapse categories")
    if (length(emp) && thr == "free")
      stop("interior category never used (in either orientation): ",
           paste(cats[emp + 1L], collapse = ", "),
           "; use thresholds = 'pc' (pooled principal-component structure)",
           " or collapse categories")
    if (length(emp) && note_interior)
      notes <<- c(notes, sprintf(
        "interior category unused (%s); thresholds pooled by the principal-component structure",
        paste(cats[emp + 1L], collapse = ", ")))
    invisible(NULL)
  }
  check_cats(x, note_interior = TRUE)

  # objects whose every response sits at the boundary have no finite
  # estimate, as extreme persons are set aside in a Rasch calibration
  removed_any <- FALSE
  repeat {
    objs <- sort(unique(c(a, b)))
    T_of <- setNames(numeric(length(objs)), objs); N_of <- T_of
    for (r in seq_along(a)) {
      T_of[a[r]] <- T_of[a[r]] + w[r] * x[r]
      T_of[b[r]] <- T_of[b[r]] + w[r] * (m - x[r])
      N_of[a[r]] <- N_of[a[r]] + w[r] * m
      N_of[b[r]] <- N_of[b[r]] + w[r] * m
    }
    ext <- names(T_of)[T_of == 0 | T_of == N_of]
    if (!length(ext)) break
    # an anchored object at a boundary is held at a known location, so it is
    # not silently removed the way a free extreme object is; equating cannot
    # proceed on a scale whose anchor has no comparisons to place it against
    if (!is.null(anchors) && any(ext %in% names(anchors)))
      stop("anchored object(s) at a response boundary (undefeated or winless): ",
           paste(intersect(ext, names(anchors)), collapse = ", "),
           "; an anchored object cannot be removed - drop it from `anchors`",
           " or supply comparisons that place it")
    notes <- c(notes, sprintf(
      "object(s) at a response boundary removed (no finite estimate): %s",
      paste(ext, collapse = ", ")))
    sel <- !(a %in% ext) & !(b %in% ext)
    a <- a[sel]; b <- b[sel]; x <- x[sel]; w <- w[sel]
    if (!is.null(jd)) jd <- jd[sel]
    if (!is.null(ord)) ord <- ord[sel]
    if (!is.null(Z)) Z <- Z[sel, , drop = FALSE]
    removed_any <- TRUE
    if (!length(a)) stop("no comparisons remain after removing extreme objects")
  }
  # removing boundary objects can itself empty a category; re-check so the
  # user gets the collapse-categories error, not a singular Newton step
  if (removed_any) check_cats(x, note_interior = FALSE)
  objs <- sort(unique(c(a, b)))
  K <- length(objs)
  if (K < 3) stop("need at least three comparable objects")
  ia <- match(a, objs); ib <- match(b, objs)
  comp <- .btl_components(K, ia, ib)
  if (length(unique(comp)) > 1) {
    parts <- split(objs, comp)
    stop("the comparison graph is disconnected; components: ",
         paste(vapply(parts, paste, "", collapse = ","), collapse = " | "))
  }

  # symmetric-threshold map: tau = Cmat %*% tfree, tau_k = -tau_{m+1-k}.
  # Under thr = "pc" the free symmetric parameters are further pooled to
  # the spread (linear) component alone - Andrich's principal-component
  # structure with the even (skewness) components structurally zero under
  # symmetry - so sparse categories borrow strength from every response.
  if (thr == "pc" && m >= 2L) {
    v1 <- seq_len(m) - (m + 1) / 2
    Cmat <- cbind(v1 / sqrt(sum(v1^2)))
    q <- 1L
  } else {
    q <- m %/% 2L
    Cmat <- matrix(0, m, q)
    for (k in seq_len(q)) { Cmat[k, k] <- 1; Cmat[m + 1L - k, k] <- -1 }
  }
  # location design B (K x nb) and fixed offset beta0 (length K), so that
  # beta = Bmat %*% bfree + beta0. Without anchors this is the sum-zero map
  # (bit-identical to the pre-anchor path); with anchors it is a selection
  # matrix -- identity rows for the free objects, zero rows for the anchored
  # ones -- plus the anchored values as the offset, so the anchored locations
  # are held fixed and the rest float with no sum-zero constraint, the origin
  # and scale coming from the anchors (as in an anchored rasch() calibration).
  anch <- NULL
  if (is.null(anchors)) {
    Bmat <- rbind(diag(K - 1L), rep(-1, K - 1L))
    beta0 <- numeric(K)
  } else {
    anch <- anchors[names(anchors) %in% objs]
    if (!length(anch)) stop("no `anchors` name matches a comparable object")
    apos <- match(names(anch), objs)
    free <- setdiff(seq_len(K), apos)
    if (!length(free)) stop("every object is anchored; nothing to estimate")
    Bmat <- matrix(0, K, length(free))
    Bmat[cbind(free, seq_along(free))] <- 1
    beta0 <- numeric(K); beta0[apos] <- as.numeric(anch)
    # start the free objects AT the anchored scale's centre, as pcml's
    # anchored path shifts its start values: beginning at 0 when the anchors
    # sit several logits away sends the undamped Newton step into the
    # logistic tails and diverges, though the maximum is finite
    beta0[free] <- mean(as.numeric(anch))
    notes <- c(notes, sprintf(
      "%d object(s) anchored; scale origin from anchors", length(anch)))
  }
  nb <- ncol(Bmat)
  pz <- if (is.null(Z)) 0L else ncol(Z)
  Zfull <- Z                         # all effect columns, for the audit table
  if (pz) {
    keepz <- colSums(abs(Z)) > 0
    if (!all(keepz)) {
      notes <- c(notes, sprintf(
        "dependence effect(s) with no informative comparisons dropped: %s",
        paste(colnames(Z)[!keepz], collapse = ", ")))
      Z <- Z[, keepz, drop = FALSE]; pz <- ncol(Z)
    }
  }
  np <- nb + q + pz
  dep <- numeric(pz)
  sc <- 0:m

  # per-row moments for current parameters: probabilities, E, V, mu4,
  # survivor S_k = P(X >= k) and EXc_k = E[X 1(X >= k)], k = 1..m
  moments <- function(beta, tfree, dep) {
    tau <- if (q) drop(Cmat %*% tfree) else numeric(m)
    d <- beta[ia] - beta[ib]
    if (pz) d <- d + drop(Z %*% dep)
    eta <- outer(d, sc) - matrix(rep(c(0, cumsum(tau)), each = length(d)),
                                 length(d), m + 1L)
    eta <- eta - apply(eta, 1, max)
    P <- exp(eta); P <- P / rowSums(P)
    E <- drop(P %*% sc)
    V <- drop(P %*% sc^2) - E^2
    mu4 <- drop(P %*% sc^4) - 4 * E * drop(P %*% sc^3) +
      6 * E^2 * drop(P %*% sc^2) - 3 * E^4
    S <- t(apply(P, 1, function(p) rev(cumsum(rev(p)))))[, -1L, drop = FALSE]
    EXc <- t(apply(P * rep(sc, each = length(d)), 1,
                   function(p) rev(cumsum(rev(p)))))[, -1L, drop = FALSE]
    list(P = P, E = E, V = pmax(V, 1e-12), mu4 = mu4, S = S, EXc = EXc,
         tau = tau)
  }
  cumInd <- outer(x, seq_len(m), ">=") * 1

  # accumulate per-comparison rows (vector or matrix) into indexed slots:
  # rowsum() replaces the interpreted per-row loops that dominated large fits
  acc <- function(v, idx, nrows) {
    out <- matrix(0, nrows, if (is.matrix(v)) ncol(v) else 1L)
    rs <- rowsum(v, idx)
    out[as.integer(rownames(rs)), ] <- rs
    out
  }

  # generic gradient/Hessian over theta = (beta_red, tfree, dep): the
  # covariates enter the exponent multiplied by the score, exactly as the
  # location difference does, so every block follows the same moments
  gH <- function(mo) {
    resE <- w * (x - mo$E)
    g_beta_full <- drop(acc(resE, ia, K) - acc(resE, ib, K))
    g <- drop(crossprod(Bmat, g_beta_full))
    if (q) g <- c(g, drop(crossprod(Cmat, colSums(w * (mo$S - cumInd)))))
    if (pz) g <- c(g, drop(crossprod(Z, resE)))
    H <- matrix(0, np, np)
    hv <- w * mo$V
    Hbb <- matrix(0, K, K)
    diag(Hbb) <- drop(acc(hv, ia, K) + acc(hv, ib, K))
    # off-diagonal cells accumulated over the unordered pair, so both
    # presentation orders of the same pair land in the same cell
    lo <- pmin(ia, ib); hi <- pmax(ia, ib)
    hp <- rowsum(hv, (lo - 1L) * K + hi)
    kk <- as.integer(rownames(hp))
    i0 <- (kk - 1L) %/% K + 1L; j0 <- (kk - 1L) %% K + 1L
    Hbb[cbind(i0, j0)] <- Hbb[cbind(i0, j0)] - hp
    Hbb[cbind(j0, i0)] <- Hbb[cbind(j0, i0)] - hp
    H[1:nb, 1:nb] <- crossprod(Bmat, Hbb %*% Bmat)
    if (q) {
      CovXc <- mo$EXc - mo$E * mo$S              # rows: Cov(X, 1(X>=k))
      wc <- (w * CovXc) %*% Cmat
      Hbt_full <- acc(wc, ib, K) - acc(wc, ia, K)
      ti <- (nb + 1L):(nb + q)
      H[1:nb, ti] <- crossprod(Bmat, Hbt_full)
      H[ti, 1:nb] <- t(H[1:nb, ti])
      # sum_r w_r Cov(1(X>=i), 1(X>=j)) = ws[max(i,j)] - crossprod term,
      # since S_r[max(i, j)] depends only on the larger index
      ws <- colSums(w * mo$S)
      Mcc <- outer(seq_len(m), seq_len(m), function(i, j) ws[pmax(i, j)]) -
        crossprod(mo$S, w * mo$S)
      H[ti, ti] <- crossprod(Cmat, Mcc %*% Cmat)
    }
    if (pz) {
      zi <- (nb + q + 1L):np
      wv <- w * mo$V
      H[zi, zi] <- crossprod(Z, Z * wv)
      Hbz_full <- acc(Z * wv, ia, K) - acc(Z * wv, ib, K)
      H[1:nb, zi] <- crossprod(Bmat, Hbz_full)
      H[zi, 1:nb] <- t(H[1:nb, zi])
      if (q) {
        CovXc <- mo$EXc - mo$E * mo$S
        wc <- (w * CovXc) %*% Cmat
        Htz <- -crossprod(wc, Z)                     # q x pz
        H[(nb + 1L):(nb + q), zi] <- Htz
        H[zi, (nb + 1L):(nb + q)] <- t(Htz)
      }
    }
    list(g = g, H = H, resE = resE)
  }

  beta <- beta0; tfree <- numeric(q)
  repeat {
    for (it in seq_len(maxit)) {
      mo <- moments(beta, tfree, dep)
      gh <- gH(mo)
      step <- solve(gh$H, gh$g)
      # trust-region damp: an undamped Newton step of many logits overshoots
      # into the flat logistic tails (seen with distant anchors); inert for
      # ordinary fits, whose steps are far smaller
      ms <- max(abs(step))
      if (is.finite(ms) && ms > 5) step <- step * (5 / ms)
      beta <- beta + drop(Bmat %*% step[1:nb])
      if (q) tfree <- tfree + step[(nb + 1L):(nb + q)]
      if (pz) dep <- dep + step[(nb + q + 1L):np]
      if (max(abs(step)) < tol) break
    }
    # a dependence effect running to a boundary is separation: its informative
    # comparisons all point one way, so the data are evidence of an infinite
    # effect (as an all-wins object is of an infinite location). The gradient
    # plateaus there, so the usual test would report convergence with an
    # arbitrary, maxit-dependent estimate; instead the column is set aside
    # with a note and the model refitted without it.
    runaway <- if (pz) which(abs(dep) > 10) else integer(0)
    if (!length(runaway)) break
    notes <- c(notes, sprintf(
      "dependence effect(s) separated (one-sided informative comparisons) and dropped: %s",
      paste(colnames(Z)[runaway], collapse = ", ")))
    Z <- Z[, -runaway, drop = FALSE]; pz <- ncol(Z)
    np <- nb + q + pz
    dep <- numeric(pz)
    beta <- beta0; tfree <- numeric(q)
  }
  mo <- moments(beta, tfree, dep)
  tau <- mo$tau
  loglik <- sum(w * log(pmax(mo$P[cbind(seq_along(x), x + 1L)], 1e-300)))
  gh <- gH(mo)
  resE <- gh$resE
  # scale-free convergence: gradient per (count-weighted) comparison, so
  # duplicated or very large data cannot flag a converged fit as unconverged
  converged <- isTRUE(max(abs(gh$g)) < 1e-6 * sum(w))
  if (!converged)
    warning("btl estimation did NOT converge in ", it, " iterations: ",
            "estimates, standard errors, and fit statistics are ",
            "unreliable -- increase maxit or check the comparison design",
            call. = FALSE)

  # Godambe sandwich over the full parameter, clustered by judge (each
  # cluster's score contributions accumulated by rowsum, not per-row loops)
  cl <- if (is.null(jd)) as.character(seq_along(ia)) else jd
  ucl <- unique(cl)
  nc <- length(ucl); cidx <- match(cl, ucl)
  # the clustered sandwich estimates the meat from between-judge variation:
  # with one judge it is identically ~zero (SEs collapse to ~1e-16), and
  # with very few judges it understates -- refuse the former, note the latter
  if (!is.null(jd)) {
    if (nc < 2L)
      stop("judge-clustered standard errors need at least 2 judges (got ",
           nc, "); with a single judge drop judge= so comparisons are ",
           "treated as independent, or supply more judges")
    if (nc < 10L)
      notes <- c(notes, sprintf(
        "only %d judge clusters: judge-clustered standard errors are likely to understate with so few clusters", nc))
  }
  Gm <- matrix(0, nc, np, dimnames = list(ucl, NULL))
  # beta block: per-cluster sums of resE into the winner / loser slots,
  # laid out cluster-major so one rowsum fills the (cluster x object) grid
  gA <- drop(acc(resE, (cidx - 1L) * K + ia, nc * K))
  gB <- drop(acc(resE, (cidx - 1L) * K + ib, nc * K))
  Gm[, 1:nb] <- t(matrix(gA - gB, nrow = K)) %*% Bmat
  if (q) {
    st_tau <- (w * (mo$S - cumInd)) %*% Cmat
    Gm[, (nb + 1L):(nb + q)] <- acc(st_tau, cidx, nc)
  }
  if (pz) Gm[, (nb + q + 1L):np] <- acc(Z * resE, cidx, nc)
  # count-weighted rows with NO judge stand for w INDEPENDENT comparisons,
  # each its own cluster: the meat must accumulate w * (x - E)^2 per row,
  # not (w * (x - E))^2 -- the per-row scores carry w, so divide by sqrt(w)
  # (with a judge the w replicates share the judge's cluster, where the
  # summed score w * (x - E) is exactly right and nothing is rescaled)
  if (is.null(jd)) Gm <- Gm / sqrt(w)
  H <- gh$H
  Hi <- solve(H)
  # CR1 small-sample factor: with G clusters the empirical meat understates
  # by ~G/(G-1); the correction is standard practice and matters exactly
  # where few-cluster inference is already fragile (a note fires below 10)
  cr1 <- if (!is.null(jd) && nc > 1L) nc / (nc - 1) else 1
  # with fewer clusters than parameters the empirical meat is singular by
  # construction: no scalar correction repairs that, so say so
  rank_deficient <- !is.null(jd) && nc <= np
  if (rank_deficient)
    notes <- c(notes, sprintf(
      "%d judge clusters for %d parameters: the clustered covariance is rank-deficient; marginal SEs are reported as consistent but understating estimates, dependence t/p should be read as descriptive, and the OSI is withheld (understated SEs would overstate it) -- reliable clustered inference needs more judges than parameters", nc, np))
  covth <- Hi %*% (crossprod(Gm) * cr1) %*% Hi
  # composite-likelihood information ingredients: tr(H^-1 J) = tr(covth H)
  # is the effective parameter count of the Godambe penalty (Varin & Vidoni
  # 2005); abs() makes it sign-convention free (the eigenvalues of H^-1 J
  # share one sign). Independent units are judges when clustered, else the
  # count-weighted comparisons.
  cl_info <- list(eff_params = abs(sum(diag(covth %*% H))),
                  n_units = if (is.null(jd)) sum(w) else length(ucl),
                  n_parameters = np)
  # anchored objects have a zero row in Bmat, so their location variance is
  # structurally zero (se == 0): the location is a fixed constant, not an estimate
  cov_beta <- Bmat %*% covth[1:nb, 1:nb, drop = FALSE] %*% t(Bmat)
  se <- sqrt(pmax(diag(cov_beta), 0))
  dependence <- NULL
  if (pz) {
    zi <- (nb + q + 1L):np
    dse <- sqrt(pmax(diag(covth)[zi], 0))
    # clustered: the z statistics get a t reference with G - 1 degrees of
    # freedom (the standard few-cluster correction) rather than normal
    # theory, so five judges give honestly wide p-values
    t_df <- if (!is.null(jd)) max(nc - 1L, 1L) else Inf
    dependence <- data.frame(
      effect = colnames(Z), estimate = dep, se = dse,
      t = dep / dse, df = t_df,
      p = 2 * pt(-abs(dep / dse), df = t_df),
      # count-weighted: the number of comparisons (not rows) that carry
      # information about each effect
      n_informative = vapply(seq_len(ncol(Z)), function(j)
        sum(w[Z[, j] != 0]), 0))
    rownames(dependence) <- NULL
  }
  thresholds <- NULL; components <- NULL
  if (q) {
    ti <- (nb + 1L):(nb + q)
    cov_tau <- Cmat %*% covth[ti, ti, drop = FALSE] %*% t(Cmat)
    thresholds <- data.frame(threshold = seq_len(m), tau = tau,
                             se = sqrt(pmax(diag(cov_tau), 0)))
    # principal-component decomposition of the threshold structure: the
    # odd components (spread; kurtosis from five thresholds up) carry the
    # symmetric structure, the even skewness component is structurally
    # zero under presentation-order symmetry
    v1 <- seq_len(m) - (m + 1) / 2
    v1 <- v1 / sqrt(sum(v1^2))
    comp_rows <- list(data.frame(
      component = "spread", estimate = sum(v1 * tau),
      se = sqrt(pmax(drop(t(v1) %*% cov_tau %*% v1), 0))))
    if (m >= 4L) {
      v3 <- (seq_len(m) - (m + 1) / 2)^3
      v3 <- v3 - sum(v3 * v1) * v1
      v3 <- v3 / sqrt(sum(v3^2))
      comp_rows[[2]] <- data.frame(
        component = "kurtosis", estimate = sum(v3 * tau),
        se = if (thr == "pc") 0 else
          sqrt(pmax(drop(t(v3) %*% cov_tau %*% v3), 0)))
    }
    components <- do.call(rbind, comp_rows)
    rownames(components) <- NULL
  } else if (m > 1L) {
    thresholds <- data.frame(threshold = seq_len(m), tau = tau, se = 0)
  }

  # fit: per-comparison z; objects and judges pool their cells
  z <- (x - mo$E) / sqrt(mo$V)
  c4v <- mo$mu4 / mo$V^2 - 1
  n_rows <- sum(w)
  f_cell <- (n_rows - np) / n_rows
  pool <- function(sel) {
    if (sum(w[sel]) < 3)
      return(list(infit_ms = NA_real_, outfit_ms = NA_real_,
                  fit_resid = NA_real_, df = NA_real_, n = sum(w[sel])))
    y2 <- sum(w[sel] * z[sel]^2); f <- f_cell * sum(w[sel])
    # information-weighted infit, as in the dichotomous path but over the
    # graded response variance
    wv <- sum(w[sel] * mo$V[sel])
    infit <- if (wv > 1e-12)
      sum(w[sel] * z[sel]^2 * mo$V[sel]) / (f_cell * wv) else NA_real_
    v <- sum(w[sel] * c4v[sel])
    fr <- if (v > 1e-8 && y2 > 0) f * (log(y2) - log(f)) / sqrt(v) else NA_real_
    list(infit_ms = infit, outfit_ms = y2 / f, fit_resid = fr, df = f,
         n = sum(w[sel]))
  }
  ofit <- lapply(seq_len(K), function(k) pool(ia == k | ib == k))
  score_of <- vapply(seq_len(K), function(k)
    sum(w[ia == k] * x[ia == k]) + sum(w[ib == k] * (m - x[ib == k])), 0)
  objects <- data.frame(object = objs, location = beta, se = se,
                        comparisons = vapply(ofit, `[[`, 0, "n"),
                        score = score_of,
                        infit_ms = vapply(ofit, `[[`, 0, "infit_ms"),
                        outfit_ms = vapply(ofit, `[[`, 0, "outfit_ms"),
                        fit_resid = vapply(ofit, `[[`, 0, "fit_resid"),
                        df_fit = vapply(ofit, `[[`, 0, "df"))
  rownames(objects) <- NULL

  judges <- NULL
  if (!is.null(jd)) {
    ju <- sort(unique(jd))
    jfit <- lapply(ju, function(j) pool(jd == j))
    judges <- data.frame(judge = ju,
                         n = vapply(jfit, `[[`, 0, "n"),
                         infit_ms = vapply(jfit, `[[`, 0, "infit_ms"),
                         outfit_ms = vapply(jfit, `[[`, 0, "outfit_ms"),
                         fit_resid = vapply(jfit, `[[`, 0, "fit_resid"),
                         df_fit = vapply(jfit, `[[`, 0, "df"))
    rownames(judges) <- NULL
  }

  # pairwise goodness of fit on the oriented mean response
  key <- ifelse(ia < ib, paste(ia, ib), paste(ib, ia))
  x_lo <- ifelse(ia < ib, x, m - x)
  E_lo <- ifelse(ia < ib, mo$E, m - mo$E)
  n_pair <- tapply(w, key, sum)
  obs_m <- tapply(w * x_lo, key, sum) / n_pair
  exp_m <- tapply(w * E_lo, key, sum) / n_pair
  v_pair <- tapply(w * mo$V, key, sum)
  zp <- tapply(w * (x_lo - E_lo), key, sum) / sqrt(pmax(v_pair, 1e-12))
  idx <- do.call(rbind, strsplit(names(n_pair), " "))
  pairs <- data.frame(object_a = objs[as.integer(idx[, 1])],
                      object_b = objs[as.integer(idx[, 2])],
                      n = as.numeric(n_pair),
                      obs_mean = as.numeric(obs_m),
                      exp_mean = as.numeric(exp_m),
                      residual = as.numeric(zp),
                      chisq = as.numeric(zp)^2)
  rownames(pairs) <- NULL
  used <- pairs$n >= 2
  # df: informative pairs minus ALL estimated parameters (locations,
  # thresholds, and any position/dependence covariates); when nothing
  # testable remains the total is NA, not a manufactured df = 1
  total_chisq <- sum(pairs$chisq[used])
  total_df <- sum(used) - np
  if (total_df < 1L) { total_chisq <- NA_real_; total_df <- NA_integer_ }
  osi <- if (rank_deficient)
    list(PSI = NA_real_, separation = NA_real_, strata = NA_real_,
         var_theta = NA_real_, mean_error_var = NA_real_, n = 0L)
  else .psi(objects$location, objects$se)

  # two categories ARE the dichotomous conditional model, so an m == 1 fit is
  # presented in dichotomous terms: the score is the win count and the mean
  # graded responses are win proportions (one estimator serves both routes)
  if (m == 1L) {
    names(objects)[names(objects) == "score"] <- "wins"
    names(pairs)[names(pairs) == "obs_mean"] <- "obs_prop"
    names(pairs)[names(pairs) == "exp_mean"] <- "exp_prop"
  }

  # the per-comparison history covariates, so the dependence effects can be
  # interrogated: a comparison is informative for an effect when its covariate
  # is non-zero (the two objects' histories differ)
  dependence_data <- if (is.null(Zfull)) NULL else {
    dd <- data.frame(judge = if (is.null(jd)) NA_character_ else jd,
                     order = if (is.null(ord)) NA_real_ else ord,
                     object_a = a, object_b = b, response = x, weight = w,
                     stringsAsFactors = FALSE)
    # whichever covariate columns the design carried (exposure, carry_over,
    # position), added by name so any subset works
    for (cn in colnames(Zfull)) dd[[cn]] <- Zfull[, cn]
    dd <- dd[order(dd$judge, dd$order), ]; rownames(dd) <- NULL; dd
  }
  out <- list(objects = objects, thresholds = thresholds,
              components = components, thr_structure = thr,
              dependence = dependence, dependence_data = dependence_data,
              pairs = pairs,
              judges = judges, m = m, categories = cats,
              total_chisq = total_chisq, total_df = total_df,
              total_p = pchisq(total_chisq, total_df, lower.tail = FALSE),
              osi = osi, loglik = loglik, iterations = it,
              converged = converged, n_comparisons = n_rows,
              clustered = !is.null(jd), cov_beta = cov_beta, cl = cl_info,
              comparisons = {
                cmp <- data.frame(object_a = a, object_b = b,
                                  response = x, weight = w,
                                  judge = if (is.null(jd))
                                    NA_character_ else jd)
                # row-aligned history covariates, so downstream analyses
                # (btl_dif) can hold the fitted dependence effects fixed
                if (!is.null(Zfull))
                  for (cn in colnames(Zfull)) cmp[[cn]] <- Zfull[, cn]
                cmp
              },
              anchors = anch,
              notes = notes)
  class(out) <- "rasch_btl"
  out
}

#' Plot graded-comparison category curves
#'
#' For a graded paired-comparison fit, the probability of each response
#' category as a function of the location difference
#' \code{beta_a - beta_b}, with the symmetric threshold structure marked.
#' The display is the paired-comparison counterpart of the category
#' probability curves of a polytomous item.
#'
#' @param fit A graded fit from \code{\link{btl}} (with \code{response}).
#' @param grid Difference grid, in logits.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#' P <- vapply(seq_len(nrow(d)), function(r)
#'   item_moments(beta[d$a[r]] - beta[d$b[r]], c(-1, 0, 1))$P, numeric(4))
#' d$grade <- apply(P, 2, function(p) sample(0:3, 1, prob = p))
#' plot_btl_categories(btl(d, "a", "b", response = "grade"))
#' @export
plot_btl_categories <- function(fit, grid = seq(-4, 4, 0.05)) {
  if (is.null(fit$m) || fit$m < 2L)
    stop("category curves need a graded fit (three or more categories)")
  tau <- fit$thresholds$tau
  P <- vapply(grid, function(d) item_moments(d, tau)$P, numeric(fit$m + 1L))
  op <- .rr_canvas(range(grid), c(0, 1),
                   "Location difference (logits)", "Category probability")
  on.exit(par(op))
  abline(v = tau, lty = 3, col = .rr$soft)
  labs <- if (!is.null(fit$categories)) fit$categories else
    as.character(0:fit$m)
  for (cat in 0:fit$m)
    lines(grid, P[cat + 1L, ], lwd = 2.6,
          col = .rr$pal[cat %% length(.rr$pal) + 1L])
  .rr_legend("right", labs, lwd = 2.6,
             col = .rr$pal[(0:fit$m) %% length(.rr$pal) + 1L])
  invisible(NULL)
}

#' Plot an object characteristic curve
#'
#' The paired-comparison counterpart of the item characteristic curve: the
#' model expected response for one object as a function of opponent
#' location (the win probability, or the expected graded response), with
#' the observed mean response against each opponent overlaid at that
#' opponent\'s estimated location. Observed points shrink in toward the
#' curve as the model holds; an object of inconsistent quality shows
#' points straying from it, exactly as a misfitting item does.
#'
#' @param fit An object from \code{\link{btl}}.
#' @param object Object name.
#' @param group Optional judge grouping for a DIF overlay: either one value
#'   per comparison row of \code{fit$comparisons} or a vector named by
#'   judge. Observed means are then drawn separately per group, as
#'   \code{\link{plot_icc}} draws person groups.
#' @param grid Opponent-location grid, in logits.
#' @param min_n An opponent's observed point is drawn only when the object
#'   (or, in the grouped display, that judge group) met it at least this many
#'   times; sparser pairs from incomplete or unbalanced designs are omitted.
#' @return Called for its plotting side effect; invisibly the names of the
#'   opponents drawn (the ungrouped display), or \code{NULL} for the grouped
#'   display.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' plot_btl_icc(btl(d, "a", "b", winner = "win"), "C")
#' @export
plot_btl_icc <- function(fit, object, group = NULL, grid = NULL,
                         min_n = 10) {
  ob <- fit$objects
  if (!object %in% ob$object) stop("no such object: ", object)
  m <- if (is.null(fit$m)) 1L else fit$m
  tau <- if (!is.null(fit$thresholds)) fit$thresholds$tau else numeric(1)
  b_o <- ob$location[ob$object == object]
  if (is.null(grid)) {
    rng <- range(ob$location) + c(-1, 1)
    grid <- seq(rng[1], rng[2], length.out = 201)
  }
  Ecurve <- vapply(grid, function(t) item_moments(b_o - t, tau)$E, 0)
  cm <- fit$comparisons
  gv <- NULL
  if (!is.null(group)) {
    gv <- if (length(group) == nrow(cm)) as.character(group) else {
      if (is.null(names(group)))
        stop("`group` must have one entry per comparison or be named by judge")
      unname(as.character(group)[match(cm$judge, names(group))])
    }
  }
  sel_a <- cm$object_a == object
  sel_b <- cm$object_b == object
  opp <- c(cm$object_b[sel_a], cm$object_a[sel_b])
  resp <- c(cm$response[sel_a], m - cm$response[sel_b])
  wt <- c(cm$weight[sel_a], cm$weight[sel_b])
  gg <- if (is.null(gv)) NULL else c(gv[sel_a], gv[sel_b])
  keep <- opp %in% ob$object
  if (!is.null(gg)) keep <- keep & !is.na(gg)
  opp <- opp[keep]; resp <- resp[keep]; wt <- wt[keep]
  if (!is.null(gg)) gg <- gg[keep]
  obs <- data.frame(
    opponent = tapply(opp, opp, `[`, 1),
    loc = ob$location[match(names(tapply(wt, opp, sum)), ob$object)],
    mean = as.numeric(tapply(wt * resp, opp, sum) / tapply(wt, opp, sum)),
    n = as.numeric(tapply(wt, opp, sum)))
  op <- .rr_canvas(range(grid), c(0, m), "Opponent location (logits)",
                   if (m == 1L) "Probability preferred" else
                     "Expected graded response",
                   sprintf("%s  (location %.3f)", object, b_o))
  on.exit(par(op))
  lines(grid, Ecurve, lwd = 3, col = .rr$ink)
  abline(v = b_o, lty = 3, col = .rr$soft)
  if (is.null(gg)) {
    # a comparator is shown only when the object met it enough times for the
    # observed proportion to be informative; sparser pairs (incomplete or
    # unbalanced designs) are omitted rather than plotted as noise
    shown <- obs[obs$n >= min_n, , drop = FALSE]
    n_omit <- nrow(obs) - nrow(shown)
    points(shown$loc, shown$mean, pch = 21, bg = .rr$blue,
           col = "white", cex = 1.5, lwd = 1.2)
    text(shown$loc, shown$mean, shown$opponent, pos = 3, offset = 0.45,
         cex = 0.72, col = .rr$soft)
    .rr_legend("topright",
               c("Model", "Observed (per opponent)",
                 if (n_omit)
                   sprintf("%d omitted (< %d comparisons)", n_omit, min_n)),
               lwd = c(3, NA, if (n_omit) NA),
               pch = c(NA, 21, if (n_omit) NA),
               pt.bg = c(NA, .rr$blue, if (n_omit) NA),
               col = c(.rr$ink, "white", if (n_omit) .rr$soft),
               pt.cex = 1.3)
  } else {
    # the graphical DIF display: per-opponent means drawn separately for
    # each judge group, as plot_icc draws person groups. A group's point for
    # an opponent is shown only where that group met it enough times
    levs <- sort(unique(gg))
    for (li in seq_along(levs)) {
      sel <- gg == levs[li]
      nn <- tapply(wt[sel], opp[sel], sum)
      om <- tapply(wt[sel] * resp[sel], opp[sel], sum) / nn
      om <- om[nn >= min_n]
      ol <- ob$location[match(names(om), ob$object)]
      colr <- .rr$pal[(li - 1L) %% length(.rr$pal) + 1L]
      oo <- order(ol)
      lines(ol[oo], om[oo], col = colr, lwd = 1.4, lty = 3)
      points(ol, om, pch = 21, bg = colr, col = "white", cex = 1.4,
             lwd = 1.1)
    }
    .rr_legend("topright", c("Model", levs), lwd = c(3, rep(1.4, length(levs))),
               lty = c(1, rep(3, length(levs))),
               pch = c(NA, rep(21, length(levs))),
               pt.bg = c(NA, .rr$pal[seq_along(levs)]),
               col = c(.rr$ink, .rr$pal[seq_along(levs)]), pt.cex = 1.2)
  }
  # the opponents actually drawn (ungrouped display), for inspection and tests
  invisible(if (is.null(gg)) obs[obs$n >= min_n, "opponent"] else NULL)
}

#' Plot a within-judge dependence effect
#'
#' The graphical display of a paired-comparison dependence effect, the
#' counterpart of the DIF characteristic curve. For every comparison the
#' departure of the observed response from what the object locations alone
#' predict is taken, and the contribution of the \emph{other} dependence
#' effect is removed (a partial-residual display); these departures are then
#' averaged in bins of the effect's own history covariate and plotted against
#' it, with the model's fitted contribution overlaid. Observed points that
#' rise with the covariate along the fitted line are the effect the
#' coefficient summarises; a flat, scattered cloud means the estimate rests on
#' little. Only the informative comparisons (a non-zero covariate: the two
#' objects' histories differ) carry the effect, and the count in each bin is
#' printed so a thin exposure tail is visible.
#'
#' @param fit An object from \code{\link{btl}} fitted with an \code{order}
#'   column, so \code{fit$dependence_data} is present.
#' @param effect Which effect to display: \code{"exposure"} (the seen-before
#'   advantage, the default) or \code{"carry_over"} (response dependence).
#' @param bins Number of covariate bins for the continuous carry-over display;
#'   exposure takes its three natural levels (-1, 0, +1).
#' @return Called for its plotting side effect; invisibly a data frame of the
#'   binned covariate value, observed and fitted departure, and bin count.
#' @references Davidson, R. R., & Beaver, R. J. (1977). On extending the
#'   Bradley-Terry model to incorporate within-pair order effects.
#'   \emph{Biometrics}, 33(4), 693-702.
#' @examples
#' set.seed(1)
#' beta <- c(A = -0.8, B = -0.2, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#' d$judge <- sample(sprintf("J%02d", 1:8), nrow(d), TRUE)
#' d <- d[order(d$judge), ]; d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' f <- btl(d, "a", "b", winner = "win", judge = "judge", order = "t")
#' plot_btl_dependence(f, "carry_over")
#' @export
plot_btl_dependence <- function(fit, effect = c("exposure", "carry_over"),
                                bins = 6) {
  effect <- match.arg(effect)
  dd <- fit$dependence_data
  if (is.null(dd))
    stop("no dependence data: fit btl() with an `order` (and `judge`) column")
  if (is.null(fit$dependence))
    stop("no dependence effect was estimable (see fit$notes): ",
         paste(fit$notes, collapse = "; "))
  eff <- fit$dependence[fit$dependence$effect == effect, ]
  if (!nrow(eff))
    stop("effect not estimated (see fit$notes): ", effect)
  bl <- setNames(fit$objects$location, fit$objects$object)
  m <- if (is.null(fit$m)) 1L else fit$m
  tau <- if (!is.null(fit$thresholds)) fit$thresholds$tau else numeric(1)
  dep <- setNames(fit$dependence$estimate, fit$dependence$effect)
  phi <- if ("exposure" %in% names(dep)) dep[["exposure"]] else 0
  psi <- if ("carry_over" %in% names(dep)) dep[["carry_over"]] else 0
  Emom <- function(v) vapply(v, function(t) item_moments(t, tau)$E, 0)

  # partial-residual display: hold everything but this effect at its fitted
  # value, so the plotted departure isolates this covariate's contribution
  base <- unname(bl[dd$object_a] - bl[dd$object_b])
  lin_full <- base + phi * dd$exposure + psi * dd$carry_over
  cov <- dd[[effect]]
  coef_e <- if (effect == "exposure") phi else psi
  E_other <- Emom(lin_full - coef_e * cov)
  obs <- dd$response - E_other                 # observed departure
  fit_lift <- Emom(lin_full) - E_other         # model-fitted departure

  wt <- if (is.null(dd$weight)) rep(1, nrow(dd)) else dd$weight
  if (effect == "exposure") {
    g <- factor(cov, levels = sort(unique(cov))); xb <- as.numeric(levels(g))
  } else {
    bins <- max(2L, as.integer(bins))
    br <- unique(stats::quantile(cov, seq(0, 1, length.out = bins + 1L),
                                 na.rm = TRUE))
    # a heavy mass point (many unseen pairs at 0) can collapse the quantile
    # breaks; fall back to equal-width bins rather than per-value singletons
    if (length(br) <= 2L)
      br <- unique(pretty(range(cov, na.rm = TRUE), bins))
    g <- if (length(br) > 2L) cut(cov, br, include.lowest = TRUE) else
      factor(cov)
    xw <- tapply(wt * cov, g, sum); xb <- as.numeric(xw / tapply(wt, g, sum))
  }
  # count-weighted rows stand for several comparisons: weighted bin means,
  # and the printed n is the number of comparisons, not rows
  ob <- as.numeric(tapply(wt * obs, g, sum) / tapply(wt, g, sum))
  fb <- as.numeric(tapply(wt * fit_lift, g, sum) / tapply(wt, g, sum))
  nb <- as.numeric(tapply(wt, g, sum))
  keep <- !is.na(xb) & !is.na(ob) & nb > 0
  xb <- xb[keep]; ob <- ob[keep]; fb <- fb[keep]; nb <- nb[keep]
  oo <- order(xb); xb <- xb[oo]; ob <- ob[oo]; fb <- fb[oo]; nb <- nb[oo]

  lab <- gsub("_", "-", effect)
  yl <- range(c(ob, fb, 0), na.rm = TRUE); yl <- yl + c(-1, 1) * 0.08 * (diff(yl) + 1e-6)
  xr <- range(xb); xl <- xr + c(-1, 1) * (0.12 * diff(xr) + 0.05)
  op <- .rr_canvas(xl, yl, sprintf("%s covariate", lab),
                   if (m == 1L) "Observed - expected win probability"
                   else "Observed - expected response",
                   sprintf("%s dependence: %.3f logits (SE %.3f, p = %s)",
                           lab, eff$estimate, eff$se, .fmt_p(eff$p)),
                   grid_x = TRUE)
  on.exit(par(op))
  abline(h = 0, lty = 3, col = .rr$soft)
  lines(xb, fb, lwd = 2.6, col = .rr$ink)
  points(xb, ob, pch = 21, bg = .rr$blue, col = "white", cex = 1.7, lwd = 1.2)
  text(xb, ob, nb, pos = 3, offset = 0.6, cex = 0.65, col = .rr$soft)
  .rr_legend("topleft", c("Model", "Observed (n per bin)"),
             lwd = c(2.6, NA), pch = c(NA, 21), pt.bg = c(NA, .rr$blue),
             col = c(.rr$ink, "white"), pt.cex = 1.3)
  invisible(data.frame(covariate = xb, observed = ob, fitted = fb, n = nb))
}

# ---------------------------------------------------------------------------
# DIF for paired comparisons: object-by-judge-group interaction. Judge
# severity cancels within a comparison, so group membership can only reach
# the measurement through object-specific preference - which is DIF, tested
# here by the package's two standard routes: a residual analysis of
# variance (group crossed with opponent-strength bands, the class-interval
# analogue) and resolved locations in logits (the object split into one
# copy per judge group inside a joint fit; Dittrich, Hatzinger &
# Katzenbeisser 1998 model these judge-covariate-by-object terms in the
# log-linear frame).
# ---------------------------------------------------------------------------
#' DIF analysis for paired comparisons
#'
#' Tests whether objects function differently for identifiable groups of
#' judges. One judge factor is analysed on its own; several factors are
#' modelled jointly -- with main effects by default and factor-by-factor
#' interactions optional -- exactly as \code{\link{dif_anova}} treats person
#' factors. For each object the standardised residuals of its comparisons,
#' oriented to the object, are analysed by the judge factor(s) crossed with
#' opponent-strength bands: a term is uniform DIF, its crossing with the band
#' non-uniform DIF, and a significant higher-order group term supersedes the
#' lower-order group terms built from a subset of its factors. Each term
#' flagged for uniform DIF and not superseded is then resolved -- the object
#' split into one copy per cell of the term's factors inside a joint refit --
#' and the differences between the resolved locations reported in logits with
#' judge-clustered Wald tests and the practical-significance flag, mirroring
#' \code{\link{dif_size}}. Fits with within-judge dependence effects
#' (\code{order}) keep those effects in the residual moments and in the
#' refits, so dependence is not mistaken for judge-group DIF; count-weighted
#' comparisons enter all tests with their weights.
#'
#' The screening ANOVA treats JUDGES as the independent units: residuals are
#' aggregated to one weighted mean per judge (per opponent band) and tested
#' in a split-plot design with the judge as the error unit -- group terms
#' between judges, band terms and their interactions within. Testing
#' judge-level factors against comparison-level residuals would
#' pseudo-replicate (a null simulation with judge heterogeneity and
#' arbitrary groups falsely flagged uniform DIF in 6 of 10 datasets); the
#' judge-level design is calibrated, and its power grows with the number of
#' judges per group, not the number of comparisons. Each factor level needs
#' at least two judges, and an object at least four judges overall, to be
#' testable.
#'
#' Each object is resolved against the other objects' common locations. When
#' several objects carry real DIF, resolving them one at a time can spread a
#' large effect onto clean objects as compensating, opposite-signed artificial
#' DIF (Andrich & Hagquist 2012, 2015); read large flags on several objects
#' together with that hazard in mind, and prefer resolving the largest effect
#' first and re-running.
#'
#' @param fit An object from \code{\link{btl}}.
#' @param factors A judge factor, or a named list of them, each either one
#'   value per row of \code{fit$comparisons} or a vector named by judge.
#' @param objects Objects to test; all by default.
#' @param effects \code{"main"} (default) models several factors additively
#'   (each factor's main effect and its band interaction); \code{"factorial"}
#'   also crosses the factors with one another.
#' @param p_adjust Multiplicity adjustment across objects within each term;
#'   the resolved-size probabilities are adjusted in one pool over all
#'   objects, terms, and cell pairs.
#' @param alpha Significance level for adjusted probabilities.
#' @param flag_logits Absolute resolved difference flagged as practically
#'   significant.
#' @param min_n Term cells with fewer comparisons involving the object are
#'   dropped from its resolution, with a note.
#' @param maxit,tol Newton controls for the resolution refits.
#' @return A list of class \code{"rasch_btl_dif"}: \code{summary} (one row per
#'   object and group term with the uniform F, adjusted p and partial
#'   eta-squared -- the term itself -- the non-uniform ones -- the term
#'   crossed with the opponent band -- plus \code{uniform_DIF},
#'   \code{nonuniform_DIF} and \code{superseded} flags); \code{terms} (the
#'   full per-object analysis-of-variance table); \code{levels} (resolved
#'   location and SE per object, term and cell); \code{sizes} (per object,
#'   term and cell pair: difference in logits, SE, z, adjusted p, significance
#'   and practical flags); \code{effects}, \code{factors}, and \code{notes}.
#' @references Andrich, D., & Hagquist, C. (2012). Real and artificial
#'   differential item functioning. \emph{Journal of Educational and
#'   Behavioral Statistics}, 37(3), 387-416.
#'
#'   Dittrich, R., Hatzinger, R., & Katzenbeisser, W. (1998).
#'   Modelling the effect of subject-specific covariates in paired
#'   comparison studies with an application to university rankings.
#'   \emph{Journal of the Royal Statistical Society C}, 47(4), 511-525.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 60), b = rep(pr[, 2], each = 60),
#'                 judge = sample(sprintf("J%02d", 1:12), 360, TRUE))
#' shift <- ifelse(d$judge %in% sprintf("J%02d", 1:6) & d$a == "C", 0.9,
#'          ifelse(d$judge %in% sprintf("J%02d", 1:6) & d$b == "C", -0.9, 0))
#' p <- plogis(beta[d$a] - beta[d$b] + shift)
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' f <- btl(d, "a", "b", winner = "win", judge = "judge")
#' grp <- setNames(rep(c("g1", "g2"), each = 6), sprintf("J%02d", 1:12))
#' btl_dif(f, grp, objects = "C")
#' @export
btl_dif <- function(fit, factors, objects = NULL,
                    effects = c("main", "factorial"),
                    p_adjust = "BH", alpha = 0.05, flag_logits = 0.5,
                    min_n = 20, maxit = 60, tol = 1e-8) {
  effects <- match.arg(effects)
  cm <- fit$comparisons
  if (is.null(cm)) stop("the fit carries no comparisons")
  # a single grouping is promoted to a one-factor list; several judge factors
  # are modelled jointly (main effects by default, interactions if asked)
  if (!is.list(factors)) factors <- list(group = factors)
  if (is.null(names(factors)) || any(!nzchar(names(factors))))
    names(factors) <- paste0("factor", seq_along(factors))
  fnames <- names(factors)
  gvs <- lapply(factors, function(g) {
    if (length(g) == nrow(cm)) as.character(g)
    else {
      if (is.null(names(g)))
        stop("each factor needs one value per comparison or names by judge")
      unname(as.character(g)[match(cm$judge, names(g))])
    }
  })
  # judge-group DIF tests judge ATTRIBUTES: a row-wise factor that varies
  # within a judge has no judge-level value, and the judge-level analysis
  # would silently take whichever row came first
  for (j in seq_along(gvs)) {
    nvar <- tapply(gvs[[j]], cm$judge, function(v)
      length(unique(v[!is.na(v)])))
    if (any(nvar > 1L, na.rm = TRUE))
      stop("factor '", fnames[j], "' varies within judge(s) ",
           paste(names(nvar)[which(nvar > 1L)], collapse = ", "),
           ": judge-group DIF needs judge-constant factors")
  }
  ok <- Reduce(`&`, lapply(gvs, function(g) !is.na(g)))
  safe <- paste0("f", seq_along(fnames))            # syntactic stand-ins
  op <- if (effects == "factorial") " * " else " + "
  tvars <- function(t) strsplit(t, ":", fixed = TRUE)[[1]]

  m <- if (is.null(fit$m)) 1L else fit$m
  cats <- if (!is.null(fit$categories)) fit$categories else c("0", "1")
  thr <- if (!is.null(fit$thr_structure)) fit$thr_structure else "free"
  tau <- if (!is.null(fit$thresholds)) fit$thresholds$tau else numeric(1)
  its <- if (is.null(objects)) fit$objects$object else objects
  bl <- setNames(fit$objects$location, fit$objects$object)
  jd_all <- if (all(is.na(cm$judge))) NULL else cm$judge

  # base-fit moments per comparison, including any fitted within-judge
  # dependence effects: leaving them out would push the dependence structure
  # into the residuals, where a judge-level factor would absorb it as
  # spurious DIF (and the resolved locations would absorb it as spurious
  # magnitude)
  d0 <- bl[cm$object_a] - bl[cm$object_b]
  Zc <- NULL
  if (!is.null(fit$dependence) &&
      all(fit$dependence$effect %in% names(cm))) {
    Zc <- as.matrix(cm[, fit$dependence$effect, drop = FALSE])
    d0 <- d0 + drop(Zc %*% fit$dependence$estimate)
  }
  sc <- 0:m
  eta <- outer(unname(d0), sc) -
    matrix(rep(c(0, cumsum(tau)), each = nrow(cm)), nrow(cm), m + 1L)
  eta <- eta - apply(eta, 1, max)
  P <- exp(eta); P <- P / rowSums(P)
  E <- drop(P %*% sc); V <- pmax(drop(P %*% sc^2) - E^2, 1e-12)

  # per object: the residual ANOVA z ~ (f1 [+/*] fk) * band, one row per term
  notes <- character(0); term_rows <- list()
  for (o in its) {
    sel_a <- cm$object_a == o & ok
    sel_b <- cm$object_b == o & ok
    zo <- c((cm$response[sel_a] - E[sel_a]) / sqrt(V[sel_a]),
            -(cm$response[sel_b] - E[sel_b]) / sqrt(V[sel_b]))
    opp <- c(cm$object_b[sel_a], cm$object_a[sel_b])
    wo <- c(cm$weight[sel_a], cm$weight[sel_b])
    gcols <- lapply(gvs, function(g) c(g[sel_a], g[sel_b]))
    oloc <- bl[opp]
    keep <- !is.na(oloc)
    zo <- zo[keep]; oloc <- oloc[keep]; wo <- wo[keep]
    gcols <- lapply(gcols, function(g) g[keep])
    # a count-weighted row stands for `weight` identical comparisons: the
    # residual z then has variance 1/weight, so weighted least squares is
    # exactly the expanded-rows analysis
    n_o <- sum(wo)
    jo <- c(cm$judge[sel_a], cm$judge[sel_b])[keep]
    d <- data.frame(z = zo, w = wo, judge_unit = factor(jo))
    for (j in seq_along(fnames)) d[[safe[j]]] <- factor(gcols[[j]])
    if (n_o < 10 || any(vapply(safe, function(s)
      nlevels(droplevels(d[[s]])) < 2, TRUE))) next
    nb <- if (n_o >= 90) 3L else if (n_o >= 40) 2L else 1L
    if (nb > 1L) {
      # band breaks at count-weighted quantiles of the opponent location:
      # invariant to whether comparisons arrive expanded or count-weighted
      # (row-rank cuts were not), so aggregated and expanded data land in
      # identical bands
      ord <- order(oloc); cw <- cumsum(wo[ord]) / sum(wo)
      br <- unique(vapply(seq_len(nb - 1L) / nb, function(q)
        oloc[ord][which(cw >= q - 1e-12)[1]], 0))
      d$band <- factor(findInterval(oloc, br, left.open = TRUE) + 1L)
      if (nlevels(droplevels(d$band)) < 2L) d$band <- NULL
    }
    if (is.null(d$band)) nb <- 1L
    # JUDGES, not comparisons, are the independent units here: a judge's
    # comparisons share that judge's idiosyncratic preferences, and testing
    # judge-level factors against comparison-level residuals
    # pseudo-replicates -- a null simulation with judge heterogeneity and
    # arbitrary groups falsely flagged uniform DIF in 6 of 10 datasets.
    # Aggregate to one weighted-mean residual per judge (per opponent
    # band), then a split-plot aov with the judge as the error unit: group
    # terms are tested between judges, band terms and their interactions
    # within judges -- the same design logic as the mixed-design person
    # DIF ANOVA.
    cellkey <- if (nb > 1L) interaction(d$judge_unit, d$band, drop = TRUE)
               else droplevels(d$judge_unit)
    zbar <- tapply(d$z * d$w, cellkey, sum) / tapply(d$w, cellkey, sum)
    firsts <- which(!duplicated(cellkey))
    ag <- d[firsts[match(levels(cellkey), as.character(cellkey[firsts]))],
            c("judge_unit", safe, if (nb > 1L) "band"), drop = FALSE]
    ag$z <- as.numeric(zbar)
    # between-judge tests need at least two judges per level and enough
    # judges overall to leave residual degrees of freedom
    nj_ok <- all(vapply(safe, function(sn)
      min(tapply(as.character(ag$judge_unit), ag[[sn]],
                 function(x) length(unique(x)))) >= 2L, 0) >= 1)
    if (!nj_ok || length(unique(ag$judge_unit)) < 4L) next
    # the same order-invariant machinery as the person DIF ANOVA:
    # between-judge terms by Type II sums of squares on band-centred judge
    # margins (sequential aov let entry order decide which correlated
    # judge factor flagged), band-crossing terms on the judge-by-band mean
    # matrix through orthonormal contrasts with the Greenhouse-Geisser
    # correction
    rhs_terms <- attr(stats::terms(stats::as.formula(
      paste("z ~ (", paste(safe, collapse = op), ")",
            if (nb > 1L) "* band" else ""))), "term.labels")
    bterms_o <- rhs_terms[!vapply(rhs_terms, function(tt)
      "band" %in% .term_vars(tt), TRUE)]
    wterms_o <- setdiff(rhs_terms, bterms_o)
    if (nb > 1L) {
      bmn <- tapply(ag$z, ag$band, mean)
      zc <- ag$z - as.numeric(bmn[as.character(ag$band)])
    } else zc <- ag$z
    jk <- factor(ag$judge_unit)
    pzj <- tapply(zc, jk, mean)
    jfirst <- which(!duplicated(jk))
    pdat_o <- ag[jfirst[match(levels(jk), as.character(jk[jfirst]))],
                 c("judge_unit", safe), drop = FALSE]
    pdat_o$z <- as.numeric(pzj)
    ft_b <- .dif_type2(pdat_o, bterms_o)
    ft_w <- NULL
    if (nb > 1L && length(wterms_o)) {
      Yw <- tapply(ag$z, list(jk, ag$band), mean)
      compl <- rowSums(is.na(Yw)) == 0L
      if (sum(compl) >= 6L) {
        Yb <- Yw[compl, , drop = FALSE]
        pd2 <- pdat_o[match(rownames(Yb),
                            as.character(pdat_o$judge_unit)), ,
                      drop = FALSE]
        ft_w <- .dif_within_tests(Yb, pd2, "band",
                                  list(band = ncol(Yw)), wterms_o,
                                  bterms_o)
      }
    }
    ft <- rbind(ft_b, ft_w)
    if (is.null(ft)) next
    ft <- ft[ft$term != "Residuals" & is.finite(ft$F_value), , drop = FALSE]
    for (k in seq_len(nrow(ft)))
      term_rows[[length(term_rows) + 1L]] <- data.frame(
        object = o, term = ft$term[k], df = ft$df[k], sum_sq = ft$sum_sq[k],
        F_value = ft$F_value[k], p = ft$p[k], resid_ss = ft$resid_ss[k])
  }
  if (!length(term_rows)) stop("no object yielded an estimable DIF ANOVA")
  terms <- do.call(rbind, term_rows); rownames(terms) <- NULL
  terms$eta2_partial <- terms$sum_sq / (terms$sum_sq + terms$resid_ss)
  terms$resid_ss <- NULL
  # adjust across objects within each term
  terms$p_adj <- NA_real_
  for (tt in unique(terms$term)) {
    sel <- terms$term == tt
    terms$p_adj[sel] <- p.adjust(terms$p[sel], method = p_adjust)
  }
  terms$significant <- !is.na(terms$p_adj) & terms$p_adj < alpha
  # a significant higher-order GROUP term supersedes lower-order group terms
  # built from a subset of its factors, within the same object. Band-crossing
  # terms are excluded from the pass: a term's own band interaction is
  # reported WITH it (as non-uniform DIF), so it must not supersede it.
  terms$superseded <- FALSE
  is_group <- !vapply(terms$term, function(t) "band" %in% tvars(t), TRUE)
  for (ob in unique(terms$object)) {
    sel <- which(terms$object == ob & terms$significant & is_group)
    for (i in sel) for (k in sel) if (i != k) {
      vi <- tvars(terms$term[i]); vk <- tvars(terms$term[k])
      if (length(vi) < length(vk) && all(vi %in% vk))
        terms$superseded[i] <- TRUE
    }
  }
  # map a term's syntactic stand-ins (f1..fk) back to the nominated factor
  # names by exact whole-token match, so a factor named like a stand-in ("f1")
  # or like the opponent band cannot be re-substituted or collide. Applied only
  # for display, after all term classification is done on the stand-ins.
  relab <- function(x) vapply(x, function(t) {
    toks <- strsplit(t, ":", fixed = TRUE)[[1]]
    i <- match(toks, safe); toks[!is.na(i)] <- fnames[i[!is.na(i)]]
    # a user factor literally named "band" would otherwise be
    # indistinguishable from the opponent-strength band in the display
    if ("band" %in% fnames)
      toks[is.na(i) & toks == "band"] <- "(opponent band)"
    paste(toks, collapse = ":")
  }, character(1), USE.NAMES = FALSE)

  # compact reading: one row per object and group term (a term not crossing the
  # opponent band), its own effect uniform DIF and its band crossing
  # non-uniform DIF. Classified on the stand-in tokens, so a factor named
  # "band" (held as f_j) is never confused with the band variable.
  gterms <- unique(terms$term)   # Residuals never reaches term_rows
  gterms <- gterms[!vapply(gterms, function(t) "band" %in% tvars(t), TRUE)]
  srows <- list()
  for (ob in unique(terms$object)) for (tt in gterms) {
    u <- terms[terms$object == ob & terms$term == tt, , drop = FALSE]
    if (!nrow(u)) next
    nu <- terms[terms$object == ob & terms$term == paste0(tt, ":band"), ,
                drop = FALSE]
    srows[[length(srows) + 1L]] <- data.frame(
      object = ob, term = tt,
      F_uniform = u$F_value, p_uniform = u$p, p_uniform_adj = u$p_adj,
      eta2_uniform = u$eta2_partial, uniform_DIF = isTRUE(u$significant),
      F_nonuniform = if (nrow(nu)) nu$F_value else NA_real_,
      p_nonuniform = if (nrow(nu)) nu$p else NA_real_,
      p_nonuniform_adj = if (nrow(nu)) nu$p_adj else NA_real_,
      eta2_nonuniform = if (nrow(nu)) nu$eta2_partial else NA_real_,
      nonuniform_DIF = nrow(nu) > 0 && isTRUE(nu$significant),
      superseded = isTRUE(u$superseded))
  }
  summary_tab <- if (length(srows)) do.call(rbind, srows) else NULL

  # resolution: for each flagged, non-superseded group term, resolve the object
  # into one copy per cell of the term's factors and report the location
  # differences in logits (a main-effect term by its levels, an interaction by
  # its factor-combination cells)
  lev_rows <- list(); sz_rows <- list()
  flagged <- if (is.null(summary_tab)) integer(0) else
    which(summary_tab$uniform_DIF & !summary_tab$superseded)
  for (r in flagged) {
    ob <- summary_tab$object[r]; tt <- summary_tab$term[r]; ttd <- relab(tt)
    jf <- match(tvars(tt), safe)
    cell <- do.call(paste, c(lapply(jf, function(j) gvs[[j]]), sep = ":"))
    inv <- ok & (cm$object_a == ob | cm$object_b == ob)
    # cell sizes in comparisons (count-weighted), not rows
    lev_n <- tapply(cm$weight[inv], cell[inv], sum)
    use_lev <- names(lev_n)[lev_n >= min_n]
    if (length(use_lev) < 2) {
      notes <- c(notes, sprintf(
        "%s [%s]: fewer than two cells with %d+ comparisons; not resolved",
        ob, ttd, min_n))
      next
    }
    if (length(use_lev) < length(lev_n))
      notes <- c(notes, sprintf(
        "%s [%s]: cell(s) dropped with fewer than %d comparisons: %s",
        ob, ttd, min_n, paste(setdiff(names(lev_n), use_lev), collapse = ", ")))
    rsel <- ok & (!(cm$object_a == ob | cm$object_b == ob) | cell %in% use_lev)
    a2 <- cm$object_a[rsel]; b2 <- cm$object_b[rsel]; c2 <- cell[rsel]
    a2 <- ifelse(a2 == ob, paste0(ob, " (", c2, ")"), a2)
    b2 <- ifelse(b2 == ob, paste0(ob, " (", c2, ")"), b2)
    # the refit keeps the fitted dependence structure: the history covariates
    # are fixed by the original judgment sequence, so they pass through as-is
    rf <- tryCatch(.btl_graded(
      a2, b2, cm$response[rsel], if (is.null(jd_all)) NULL else jd_all[rsel],
      cm$weight[rsel], cats, maxit, tol, character(0), thr = thr,
      Z = if (is.null(Zc)) NULL else Zc[rsel, , drop = FALSE]),
      error = function(e) NULL)
    if (is.null(rf)) {
      notes <- c(notes, sprintf("%s [%s]: resolution failed", ob, ttd))
      next
    }
    idx <- match(paste0(ob, " (", use_lev, ")"), rf$objects$object)
    if (anyNA(idx)) {
      notes <- c(notes, sprintf("%s [%s]: resolved copies missing", ob, ttd))
      next
    }
    loc <- rf$objects$location[idx]
    vv <- rf$cov_beta[idx, idx, drop = FALSE]
    lev_rows[[length(lev_rows) + 1L]] <- data.frame(
      object = ob, term = tt, level = use_lev, location = loc,
      se = sqrt(pmax(diag(vv), 0)), n = as.numeric(lev_n[use_lev]))
    pr <- t(utils::combn(seq_along(use_lev), 2))
    sz_rows[[length(sz_rows) + 1L]] <- data.frame(
      object = ob, term = tt,
      level_a = use_lev[pr[, 1]], level_b = use_lev[pr[, 2]],
      difference = loc[pr[, 1]] - loc[pr[, 2]],
      se = sqrt(pmax(diag(vv)[pr[, 1]] + diag(vv)[pr[, 2]] -
                     2 * vv[pr], 1e-12)))
  }
  # a summary row can be flagged yet carry no magnitude row; say so rather
  # than leave the omission silent
  if (!is.null(summary_tab)) {
    for (r in which(summary_tab$uniform_DIF & summary_tab$superseded))
      notes <- c(notes, sprintf(
        "%s [%s]: uniform DIF superseded by a higher-order term; not resolved separately",
        summary_tab$object[r], relab(summary_tab$term[r])))
    for (r in which(summary_tab$nonuniform_DIF & !summary_tab$uniform_DIF))
      notes <- c(notes, sprintf(
        "%s [%s]: non-uniform DIF only; no single location difference summarises it, so no magnitude row is reported",
        summary_tab$object[r], relab(summary_tab$term[r])))
  }
  levels_df <- if (length(lev_rows)) do.call(rbind, lev_rows) else NULL
  sizes <- if (length(sz_rows)) do.call(rbind, sz_rows) else NULL
  if (!is.null(sizes)) {
    sizes$z <- sizes$difference / sizes$se
    sizes$p <- 2 * pnorm(-abs(sizes$z))
    sizes$p_adj <- p.adjust(sizes$p, method = p_adjust)
    sizes$significant <- sizes$p_adj < alpha
    sizes$practical <- abs(sizes$difference) >= flag_logits
    rownames(sizes) <- NULL
  }
  if (!is.null(levels_df)) rownames(levels_df) <- NULL
  # relabel stand-ins to the factor names for display, now that all term
  # classification and resolution are done
  terms$term <- relab(terms$term)
  if (!is.null(summary_tab)) {
    summary_tab$term <- relab(summary_tab$term); rownames(summary_tab) <- NULL
  }
  if (!is.null(sizes)) sizes$term <- relab(sizes$term)
  if (!is.null(levels_df)) levels_df$term <- relab(levels_df$term)
  out <- list(summary = summary_tab, terms = terms, levels = levels_df,
              sizes = sizes, effects = effects, factors = fnames,
              alpha = alpha, p_adjust = p_adjust, flag_logits = flag_logits,
              notes = unique(notes))
  class(out) <- "rasch_btl_dif"
  out
}

#' @export
print.rasch_btl_dif <- function(x, ...) {
  nf <- length(x$factors)
  cat(sprintf("DIF for paired comparisons: %d factor(s) [%s], %s effects\n",
              nf, paste(x$factors, collapse = ", "), x$effects))
  cat("Residual ANOVA per object and term (uniform = term; non-uniform = term x opponent band)\n")
  print(.fmt_df(x$summary[, c("object", "term", "F_uniform", "p_uniform_adj",
                              "uniform_DIF", "F_nonuniform",
                              "p_nonuniform_adj", "nonuniform_DIF")]),
        row.names = FALSE)
  if (!is.null(x$sizes)) {
    cat(sprintf("\nResolved locations (logits; %s over %d comparison(s); practical %.2f)\n",
                x$p_adjust, nrow(x$sizes), x$flag_logits))
    print(.fmt_df(x$sizes[, c("object", "term", "level_a", "level_b",
                              "difference", "se", "z", "p_adj", "significant",
                              "practical")]), row.names = FALSE)
  }
  if (length(x$notes)) cat("Notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}
