# rasch :: estimation
# ===========================================================================
# Pairwise conditional maximum likelihood after Andrich & Luo (2003) and
# Zwinderman (1995). For items i, j with maximum scores m_i,
# m_j, the distribution of X_i given the pair total X_i + X_j = r is free of
# the person parameter:
#
#   P(X_i = k | X_i + X_j = r) = exp(-L_i(k) - L_j(r-k)) / sum_k' exp(...)
#
# where L_i(k) = sum_{h<=k} tau_ih is the cumulative threshold sum. The
# pairwise conditional log-likelihood, summed over all item pairs and pair
# totals, is maximised by Newton-Raphson. Standard errors come from the
# observed information of the pseudo-likelihood. The rating scale model is
# the same likelihood under the constraint tau_ik = delta_i + kappa_k,
# imposed through the design matrix. Dichotomous data is the special case
# m_i = 1. Australian English; no em dashes by house style.
# ===========================================================================

#' Enumerate item-category thresholds
#'
#' Builds the index mapping each item-category threshold to a global id, given
#' the maximum score of each item.
#'
#' @param m Integer vector of maximum scores per item (1 for dichotomous items).
#' @return A data frame with columns \code{id}, \code{item}, and \code{k} (the
#'   within-item threshold number).
#' @examples
#' threshold_index(c(1, 3, 2))
#' @export
threshold_index <- function(m) {
  thr <- do.call(rbind, lapply(seq_along(m), function(i)
    if (m[i] >= 1) data.frame(item = i, k = seq_len(m[i])) else NULL))
  thr$id <- seq_len(nrow(thr)); thr[, c("id", "item", "k")]
}

# Cross-tabulated pair counts: for every item pair i < j, the matrix of
# joint category counts over persons observed on both items.
.pair_counts <- function(X, m) {
  L <- ncol(X); out <- vector("list", 0L)
  for (i in seq_len(L - 1)) for (j in (i + 1):L) {
    both <- !is.na(X[, i]) & !is.na(X[, j])
    if (!any(both)) next
    idx <- X[both, i] * (m[j] + 1L) + X[both, j] + 1L
    n <- matrix(tabulate(idx, nbins = (m[i] + 1L) * (m[j] + 1L)),
                nrow = m[i] + 1L, byrow = TRUE)
    out[[length(out) + 1L]] <- list(i = i, j = j, n = n)
  }
  out
}

# ---------------------------------------------------------------------------
# Structural identification checks, run before solving. The pairwise
# conditional likelihood factorises over item pairs, so relative locations
# between two blocks of items are identified only if some person answered
# items in both (the item-pair graph is connected); on a disconnected design
# the likelihood is flat in the between-block shift and Newton lands wherever
# the ridge sends it -- garbage that must be an error, not a result. Anchors
# rescue a block: a block containing an anchored item has its origin fixed.
# ---------------------------------------------------------------------------
.pcml_check_connected <- function(pairs, L, item_names, anchored = integer(0)) {
  edges <- if (length(pairs))
    do.call(rbind, lapply(pairs, function(p) c(p$i, p$j)))
  else matrix(integer(0), 0L, 2L)
  comp <- .btlef_components(L, edges)
  if (length(unique(comp)) == 1L) return(invisible(comp))
  if (length(anchored)) {
    bad <- setdiff(unique(comp), unique(comp[anchored]))
    if (!length(bad)) return(invisible(comp))
    blocks <- vapply(bad, function(cc)
      paste(item_names[comp == cc], collapse = ", "), "")
    stop("the item-pair graph is not connected, and block(s) without an ",
         "anchored item have no identified origin: ",
         paste0("{", blocks, "}", collapse = " "),
         "; link the blocks through common persons or anchor an item in ",
         "every block", call. = FALSE)
  }
  blocks <- tapply(item_names, comp, paste, collapse = ", ")
  stop("the item-pair graph is not connected: no person answered items in ",
       "more than one of the blocks ",
       paste0("{", blocks, "}", collapse = " | "),
       "; relative locations between the blocks are unidentified -- link ",
       "them through common items or persons, or anchor item(s) in every ",
       "block", call. = FALSE)
}

# A threshold k is informed by the persons observed in its two adjacent
# categories; when either count is tiny the conditional estimate can run
# away (a category with one response sends its threshold toward the
# boundary) while the ridged covariance reports a spuriously small standard
# error. Flag such thresholds so the caller can report the estimate with an
# NA standard error and a note naming the cause -- an honest answer, not a
# manufactured one. Categories with zero responses are the caller's problem
# (rasch() rescores them away); the danger zone handled here is 1-2.
.pcml_weak_thresholds <- function(X, m, thr, item_names, min_count = 3L) {
  flag <- logical(nrow(thr)); notes <- character(0)
  for (i in seq_len(ncol(X))) {
    cnt <- tabulate(X[, i] + 1L, nbins = m[i] + 1L)
    weak_k <- which(pmin(cnt[-length(cnt)], cnt[-1]) < min_count)
    if (!length(weak_k)) next
    flag[thr$item == i & thr$k %in% weak_k] <- TRUE
    kc <- which(cnt < min_count) - 1L
    notes <- c(notes, sprintf(
      "item %s: only %s response(s) in category %s; threshold(s) %s and the item location are weakly determined (SE reported as NA) -- consider pc_components or collapsing categories",
      item_names[i], paste(cnt[kc + 1L], collapse = "/"),
      paste(kc, collapse = "/"), paste(weak_k, collapse = "/")))
  }
  list(flag = flag, notes = notes)
}

# ---------------------------------------------------------------------------
# Starting values: weighted least squares on the pairwise log-ratios. Used
# only to seed Newton-Raphson; the returned estimates always come from the
# conditional likelihood itself.
# ---------------------------------------------------------------------------
.start_tau <- function(X, thr, cont = 0.5) {
  M <- nrow(thr)
  D <- matrix(NA_real_, M, M); W <- matrix(0, M, M)
  for (p in seq_len(M)) {
    i <- thr$item[p]; k <- thr$k[p]
    for (q in seq_len(M)) {
      j <- thr$item[q]; l <- thr$k[q]
      if (i == j) next
      both <- !is.na(X[, i]) & !is.na(X[, j])
      cA <- sum(X[both, i] == (k - 1L) & X[both, j] == l)
      cB <- sum(X[both, i] == k        & X[both, j] == (l - 1L))
      if (cA + cB > 0) {
        a <- cA; b <- cB
        if (cA == 0 || cB == 0) { a <- cA + cont; b <- cB + cont }
        D[p, q] <- log(a / b)
        W[p, q] <- (a * b) / (a + b)
      }
    }
  }
  rows <- which(!is.na(D) & upper.tri(D), arr.ind = TRUE)
  if (!nrow(rows)) return(rep(0, M))
  C <- matrix(0, nrow(rows), M); d <- numeric(nrow(rows)); wt <- numeric(nrow(rows))
  for (r in seq_len(nrow(rows))) {
    p <- rows[r, 1]; q <- rows[r, 2]
    C[r, p] <- 1; C[r, q] <- -1; d[r] <- D[p, q]; wt[r] <- W[p, q]
  }
  sw <- sqrt(wt)
  tau <- c(qr.coef(qr((C[, -M, drop = FALSE]) * sw), d * sw), 0)
  tau[is.na(tau)] <- 0
  tau - mean(tau)
}

# Pseudo log-likelihood, gradient, and Hessian over the full threshold vector.
.pcml_glh <- function(tau, thr, pairs, m) {
  M <- nrow(thr)
  g <- numeric(M); H <- matrix(0, M, M); ll <- 0
  cum <- lapply(seq_along(m), function(i) cumsum(tau[thr$item == i]))
  ids <- lapply(seq_along(m), function(i) thr$id[thr$item == i])
  for (pc in pairs) {
    i <- pc$i; j <- pc$j; n <- pc$n
    Li <- c(0, cum[[i]]); Lj <- c(0, cum[[j]])
    idx <- c(ids[[i]], ids[[j]]); mi <- m[i]; mj <- m[j]
    for (r in seq_len(mi + mj - 1L)) {
      ks <- max(0L, r - mj):min(mi, r)
      if (length(ks) < 2L) next
      nk <- n[cbind(ks + 1L, r - ks + 1L)]
      N <- sum(nk)
      if (N == 0) next
      lp <- -(Li[ks + 1L] + Lj[r - ks + 1L])
      lp <- lp - max(lp); elp <- exp(lp); p <- elp / sum(elp)
      ll <- ll + sum(nk * (lp - log(sum(elp))))
      # local coefficients d lp_k / d tau, item i columns then item j columns
      U <- cbind(
        -outer(ks, seq_len(mi), ">="),
        -outer(r - ks, seq_len(mj), ">="))
      storage.mode(U) <- "double"
      g[idx] <- g[idx] + drop(crossprod(U, nk - N * p))
      Ep <- drop(crossprod(U, p))
      S  <- crossprod(U, p * U)
      H[idx, idx] <- H[idx, idx] - N * (S - tcrossprod(Ep))
    }
  }
  list(ll = ll, g = g, H = H)
}

# Godambe sandwich covariance for the pairwise pseudo-likelihood. The naive
# inverse information overstates precision because every response enters
# L - 1 overlapping pairs; the sandwich H^-1 J H^-1 with J the empirical
# covariance of the per-person scores corrects this.
.pcml_sandwich <- function(X, thr, m, tau, pairs) {
  M <- nrow(thr); N <- nrow(X)
  cum <- lapply(seq_along(m), function(i) cumsum(tau[thr$item == i]))
  ids <- lapply(seq_along(m), function(i) thr$id[thr$item == i])
  S <- matrix(0, N, M)
  for (pc in pairs) {
    i <- pc$i; j <- pc$j; mi <- m[i]; mj <- m[j]
    Li <- c(0, cum[[i]]); Lj <- c(0, cum[[j]])
    idx <- c(ids[[i]], ids[[j]])
    # per-cell score vector u_k - ubar_r, indexed by cell (k, l)
    V <- matrix(0, (mi + 1L) * (mj + 1L), mi + mj)
    for (r in seq_len(mi + mj - 1L)) {
      ks <- max(0L, r - mj):min(mi, r)
      if (length(ks) < 2L) next
      lp <- -(Li[ks + 1L] + Lj[r - ks + 1L])
      lp <- lp - max(lp); p <- exp(lp) / sum(exp(lp))
      U <- cbind(-outer(ks, seq_len(mi), ">="),
                 -outer(r - ks, seq_len(mj), ">="))
      storage.mode(U) <- "double"
      ub <- drop(crossprod(U, p))
      V[ks * (mj + 1L) + (r - ks) + 1L, ] <- sweep(U, 2, ub)
    }
    both <- which(!is.na(X[, i]) & !is.na(X[, j]))
    if (!length(both)) next
    cell <- X[both, i] * (mj + 1L) + X[both, j] + 1L
    S[both, idx] <- S[both, idx] + V[cell, , drop = FALSE]
  }
  crossprod(S)
}

# Newton-Raphson on tau = offset + B beta, where B removes the location
# indeterminacy, imposes the rating scale or facet structure, or restricts
# estimation to the unanchored thresholds (offset carrying the anchors).
.pcml_solve <- function(X, thr, m, B, beta0, offset = 0, maxit = 60, tol = 1e-8,
                        pairs = NULL) {
  if (is.null(pairs)) pairs <- .pair_counts(X, m)
  if (!length(pairs)) stop("no informative item pairs: check the data")
  beta <- beta0
  glh <- .pcml_glh(drop(offset + B %*% beta), thr, pairs, m)
  it <- 0L
  for (it in seq_len(maxit)) {
    gb <- drop(crossprod(B, glh$g))
    Hb <- crossprod(B, glh$H %*% B)
    step <- tryCatch(solve(Hb, gb), error = function(e)
      solve(Hb - diag(1e-8, nrow(Hb)), gb))
    # step halving on the pseudo-likelihood
    lam <- 1; ok <- FALSE; g2 <- glh
    for (half in 1:30) {
      cand <- beta - lam * step
      g2 <- .pcml_glh(drop(offset + B %*% cand), thr, pairs, m)
      if (is.finite(g2$ll) && g2$ll >= glh$ll - 1e-12) { ok <- TRUE; break }
      lam <- lam / 2
    }
    if (!ok) break
    done <- max(abs(lam * step)) < tol
    beta <- cand; glh <- g2
    if (done) break
  }
  Hb <- crossprod(B, glh$H %*% B)
  Hinv <- tryCatch(solve(Hb), error = function(e)
    solve(Hb - diag(1e-8, nrow(Hb))))
  J  <- .pcml_sandwich(X, thr, m, drop(offset + B %*% beta), pairs)
  Jb <- crossprod(B, J %*% B)
  covb <- Hinv %*% Jb %*% Hinv
  covt <- B %*% covb %*% t(B)
  list(tau = drop(offset + B %*% beta), beta = beta, cov_beta = covb,
       cov_tau = covt, se_tau = sqrt(pmax(diag(covt), 0)), H_beta = Hb,
       loglik = glh$ll, iterations = it,
       converged = max(abs(drop(crossprod(B, glh$g)))) < 1e-4)
}

#' Estimate Rasch thresholds by pairwise conditional maximum likelihood
#'
#' Maximises the pairwise conditional likelihood, in which the person
#' parameter cancels within every item pair, by Newton-Raphson (Andrich and
#' Luo 2003; Zwinderman 1995). The partial credit model
#' estimates every threshold freely; the rating scale model constrains
#' \code{tau_ik = delta_i + kappa_k} through the design matrix.
#'
#' @param X Persons-by-items integer score matrix (categories from 0). Missing
#'   values are handled by pairwise deletion, so linked booklet designs and
#'   random missingness estimate without imputation; the item-pair graph must
#'   be connected (some person answering items in both of any two blocks),
#'   otherwise relative locations between blocks are unidentified and the fit
#'   stops with an error naming the blocks -- unless \code{anchors} fix an
#'   item in every block, the disjoint-form equating case.
#' @param model \code{"PCM"} or \code{"RSM"}.
#' @param anchors Optional anchor table for equating: a data frame with
#'   columns \code{item} (name or column index), \code{k}, and \code{tau}
#'   (the fixed value). A numeric \code{k} fixes that single threshold
#'   (individual anchoring); \code{k = NA} fixes the item's mean location at
#'   \code{tau} while its thresholds remain free (average anchoring). The
#'   remaining parameters are estimated on the anchored scale and no
#'   recentring is applied. PCM only.
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return A list with the threshold table \code{thr} (columns \code{id},
#'   \code{item}, \code{k}, \code{tau}, \code{se}, \code{anchored}, and
#'   \code{weak} -- \code{TRUE} for a threshold adjacent to a category with
#'   fewer than 3 responses, whose estimate can run toward a boundary while
#'   the ridged covariance understates the error; its \code{se} is reported
#'   as \code{NA} and a note names the item and category), the
#'   threshold covariance matrix \code{cov_tau}, the pairwise conditional
#'   log-likelihood, the iteration count, a convergence flag, \code{notes},
#'   and the max-score vector \code{m}.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' pcml(X)$thr
#' # anchor two items at fixed values (equating)
#' pcml(X, anchors = data.frame(item = c("I1", "I6"), k = 1, tau = c(-1.5, 1.5)))$thr
#' @export
pcml <- function(X, model = c("PCM", "RSM"), anchors = NULL,
                 maxit = 60, tol = 1e-8) {
  model <- match.arg(model)
  X <- as.matrix(X); storage.mode(X) <- "integer"
  m <- apply(X, 2, max, na.rm = TRUE); L <- ncol(X)
  thr <- threshold_index(m); M <- nrow(thr)
  inames <- if (is.null(colnames(X))) paste0("V", seq_len(L)) else colnames(X)
  pairs <- .pair_counts(X, m)
  weak <- .pcml_weak_thresholds(X, m, thr, inames)

  if (!is.null(anchors)) {
    if (model != "PCM") stop("anchoring is supported for the PCM only")
    if (!all(c("item", "k", "tau") %in% names(anchors)))
      stop("anchors needs columns item, k, tau")
    a_item <- if (is.character(anchors$item) || is.factor(anchors$item))
      match(as.character(anchors$item), colnames(X)) else as.integer(anchors$item)
    if (anyNA(a_item)) stop("anchor item(s) not found among the item columns")

    # k = NA anchors the item's mean location (average anchoring) with its
    # thresholds free; a numeric k fixes that single threshold. For a
    # dichotomous item the two coincide.
    is_mean <- is.na(anchors$k)
    conv <- is_mean & m[a_item] == 1L
    anchors$k[conv] <- 1L; is_mean[conv] <- FALSE

    ft_item <- a_item[!is_mean]
    a_id <- thr$id[match(paste(ft_item, anchors$k[!is_mean]),
                         paste(thr$item, thr$k))]
    if (anyNA(a_id)) stop("anchor threshold number(s) out of range for the item")
    if (anyDuplicated(a_id)) stop("duplicate anchor threshold(s)")
    mean_items <- a_item[is_mean]; mean_tau <- anchors$tau[is_mean]
    if (anyDuplicated(mean_items)) stop("duplicate average anchor(s) for an item")
    if (length(intersect(mean_items, ft_item)))
      stop("an item cannot carry both an average anchor and threshold anchors")

    offset <- numeric(M)
    offset[a_id] <- anchors$tau[!is_mean]
    for (j in seq_along(mean_items))
      offset[thr$item == mean_items[j]] <- mean_tau[j]

    # start values shifted onto the anchored scale
    st <- .start_tau(X, thr)
    shifts <- c(anchors$tau[!is_mean] - st[a_id],
                vapply(seq_along(mean_items), function(j)
                  mean_tau[j] - mean(st[thr$item == mean_items[j]]), 0))
    st <- st + mean(shifts)

    # design: identity columns for plain free thresholds; a sum-zero spread
    # block for each average-anchored item; nothing for fixed thresholds
    plain <- which(!(seq_len(M) %in% a_id) & !(thr$item %in% mean_items))
    blocks <- list()
    if (length(plain)) blocks$plain <- list(B = diag(M)[, plain, drop = FALSE],
                                            beta0 = st[plain])
    for (j in seq_along(mean_items)) {
      rows <- which(thr$item == mean_items[j]); mi <- length(rows)
      A <- matrix(0, M, mi - 1L)
      A[rows, ] <- rbind(diag(mi - 1L), rep(-1, mi - 1L))
      s <- st[rows] - mean(st[rows])
      blocks[[paste0("mean", j)]] <- list(B = A, beta0 = s[-mi])
    }
    if (!length(blocks)) stop("at least one parameter must remain free")
    B <- do.call(cbind, lapply(blocks, `[[`, "B"))
    beta0 <- unlist(lapply(blocks, `[[`, "beta0"), use.names = FALSE)

    .pcml_check_connected(pairs, L, inames,
                          anchored = unique(c(ft_item, mean_items)))
    sol <- .pcml_solve(X, thr, m, B, beta0, offset = offset,
                       maxit = maxit, tol = tol, pairs = pairs)
    thr$tau <- sol$tau; thr$se <- sol$se_tau; thr$se[a_id] <- 0
    thr$anchored <- seq_len(M) %in% a_id | thr$item %in% mean_items
    thr$weak <- weak$flag & !thr$anchored
    thr$se[thr$weak] <- NA_real_
    return(list(model = model, thr = thr, cov_tau = sol$cov_tau,
                loglik = sol$loglik, iterations = sol$iterations,
                converged = sol$converged, m = m, anchors = anchors,
                n_parameters = ncol(B), B = B, cov_beta = sol$cov_beta,
                H_beta = sol$H_beta, notes = weak$notes))
  }

  if (model == "RSM") {
    if (length(unique(m)) != 1L) stop("RSM requires equal max score across items")
    mm <- m[1]
    # parameters: delta_1..delta_{L-1}, kappa_1..kappa_{mm-1}; sum-zero each
    P <- (L - 1L) + (mm - 1L)
    B <- matrix(0, M, P)
    for (row in seq_len(M)) {
      i <- thr$item[row]; k <- thr$k[row]
      if (i < L) B[row, i] <- 1 else B[row, seq_len(L - 1L)] <- -1
      if (mm > 1L) {
        if (k < mm) B[row, L - 1L + k] <- B[row, L - 1L + k] + 1
        else B[row, L:(L + mm - 2L)] <- B[row, L:(L + mm - 2L)] - 1
      }
    }
    st <- .start_tau(X, thr)
    del <- vapply(seq_len(L), function(i) mean(st[thr$item == i]), 0)
    kap <- vapply(seq_len(mm), function(k) mean(st[thr$k == k] -
                    del[thr$item[thr$k == k]]), 0)
    del <- del - mean(del); kap <- kap - mean(kap)
    beta0 <- c(del[-L], if (mm > 1L) kap[-mm] else numeric(0))
  } else {
    # sum-zero over all thresholds during estimation; recentred afterwards
    B <- rbind(diag(M - 1L), rep(-1, M - 1L))
    st <- .start_tau(X, thr)
    beta0 <- st[-M] - mean(st)
  }

  .pcml_check_connected(pairs, L, inames)
  sol <- .pcml_solve(X, thr, m, B, beta0, maxit = maxit, tol = tol,
                     pairs = pairs)
  # recentre so the mean item location is zero
  loc <- vapply(seq_len(L), function(i) mean(sol$tau[thr$item == i]), 0)
  sol$tau <- sol$tau - mean(loc)
  thr$tau <- sol$tau; thr$se <- sol$se_tau; thr$anchored <- FALSE
  thr$weak <- weak$flag
  thr$se[thr$weak] <- NA_real_

  list(model = model, thr = thr, cov_tau = sol$cov_tau,
       loglik = sol$loglik, iterations = sol$iterations,
       converged = sol$converged, m = m, anchors = NULL,
       n_parameters = ncol(B), B = B, cov_beta = sol$cov_beta,
       H_beta = sol$H_beta, notes = weak$notes)
}

# ---------------------------------------------------------------------------
# Andrich principal-components reparameterisation (Andrich 1978, 1985; Pedler
# 1987): an optional alternative to the free-threshold pcml() above, useful
# when some categories are sparsely populated. Each item's mi thresholds are
# re-expressed as up to four orthogonal-polynomial components in the
# category score x = 0, ..., mi: location (linear), spread (quadratic),
# skewness (cubic), and kurtosis (quartic),
#
#   L_i(x) = x.omega_1i - x(mi-x).omega_2i - x(mi-x)(2x-mi).omega_3i - ...
#
# with tau_ik = L_i(k) - L_i(k-1) (Andrich 1985, eqs 1.6-1.7). Every
# component's coefficient differences telescope to zero across an item's
# thresholds except location's (which differences to a constant 1), so
# location is exactly the item's mean threshold and carries the same
# across-item additive-shift redundancy as in .start_tau; it is sum-zero
# constrained across items the same way pcml()'s RSM delta is. The
# higher-order components are free per item, capped at an item's own
# threshold count (a dichotomous item has location only). The family stops
# at the quartic (kurtosis) term, so this reproduces pcml()'s free-PCM
# thresholds exactly only while every item has at most 3 thresholds (location
# + spread + skewness then span it exactly); from 4 thresholds on it is
# necessarily a reduced-rank smoothing of the thresholds to a polynomial
# trend across categories, however large n_components is set, because no
# fifth component has been derived (Pedler 1987) and, at exactly 4
# thresholds, the quartic is collinear with the cubic (see .pc_select).
# ---------------------------------------------------------------------------
.pc_gcoefs <- function(mi) {
  x <- 0:mi
  G <- cbind(x,
             -x * (mi - x),
             -x * (mi - x) * (2 * x - mi),
             -x * (mi - x) * (2 * x - mi) * (5 * x^2 - 5 * x * mi + mi^2 + 1))
  G[-1, , drop = FALSE] - G[-nrow(G), , drop = FALSE]
}

# The quartic (kurtosis) column is collinear with the cubic (skewness) column
# at exactly 4 thresholds (it is the only such case up to 14 thresholds,
# checked exhaustively): both are non-zero only at the two threshold rows
# equidistant from the item's centre, where the quartic factor takes the same
# value by symmetry, so kurtosis carries no information beyond skewness
# there. Selecting components by incremental rank, rather than assuming
# n_components - 1 of them are always identified, catches this (and any
# other such coincidence) instead of silently returning an unidentified
# parameter and an unstable standard error.
.pc_select <- function(G, ncomp, tol = 1e-8) {
  if (ncomp < 2L) return(integer(0))
  base <- G[, 1, drop = FALSE]; r0 <- qr(base, tol = tol)$rank
  keep <- integer(0)
  for (l in 2:ncomp) {
    test <- cbind(base, G[, l])
    r1 <- qr(test, tol = tol)$rank
    if (r1 > r0) { keep <- c(keep, l); base <- test; r0 <- r1 }
  }
  keep
}

#' Estimate Rasch thresholds via the Andrich principal-components reparameterisation
#'
#' An optional alternative to \code{\link{pcml}}'s free-threshold estimation,
#' useful when some response categories are sparsely populated. Each item's
#' thresholds are re-expressed as up to four orthogonal polynomial
#' components in the category score: location, spread, skewness, and
#' kurtosis (Andrich 1978, 1985; Pedler 1987). Location is always estimated;
#' spread, skewness, and kurtosis are added in turn as an item's number of
#' thresholds and \code{n_components} allow. Estimation uses the same
#' pairwise conditional likelihood as \code{pcml} (Andrich and Luo 2003), so
#' it inherits the same missing-data handling and sandwich standard errors.
#' The component family stops at the quartic (kurtosis) term, so the
#' reparameterisation is exact, matching \code{pcml}'s free partial credit
#' thresholds and log-likelihood, only while every item has at most 3
#' thresholds (4 categories); from 4 thresholds on \code{pcml_pc} is
#' necessarily a reduced-rank smoothing of the thresholds to a polynomial
#' trend across categories, however large \code{n_components} is set,
#' trading flexibility for the stability that comes from pooling information
#' across all of an item's categories -- useful when a category has low or
#' zero frequency.
#'
#' @param X Persons-by-items integer score matrix (categories from 0).
#'   Missing values are handled by pairwise deletion.
#' @param n_components Maximum number of components per item: 1 (location
#'   only) up to 4 (location, spread, skewness, kurtosis; the highest
#'   derived by Pedler 1987). Capped per item at its own number of
#'   thresholds, and further wherever a component would be collinear with
#'   lower-order ones for that item's threshold count (kurtosis is
#'   unidentified, and dropped, at exactly 4 thresholds).
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return A list with the threshold table \code{thr} (columns \code{id},
#'   \code{item}, \code{k}, \code{tau}, \code{se}), the component table
#'   \code{components} (one row per item, with \code{location},
#'   \code{spread}, \code{skewness}, \code{kurtosis} and their standard
#'   errors, \code{NA} where an item's rank does not support that
#'   component), the threshold covariance matrix \code{cov_tau}, the
#'   pairwise conditional log-likelihood, the iteration count, a convergence
#'   flag, and the max-score vector \code{m}.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' pcml_pc(X)$components
#' @export
pcml_pc <- function(X, n_components = 4, maxit = 60, tol = 1e-8) {
  if (n_components < 1) stop("n_components must be at least 1")
  X <- as.matrix(X); storage.mode(X) <- "integer"
  m <- apply(X, 2, max, na.rm = TRUE); L <- ncol(X)
  thr <- threshold_index(m); M <- nrow(thr)
  ncomp <- pmin(m, n_components, 4L)

  Gs <- lapply(m, .pc_gcoefs)
  keep <- Map(.pc_select, Gs, ncomp)
  extra <- lengths(keep)

  # column layout: (L - 1) sum-zero location columns, then one free column
  # per item for each higher-order component actually identified for it
  P <- (L - 1L) + sum(extra)
  B <- matrix(0, M, P)
  off <- (L - 1L) + cumsum(c(0L, extra))[seq_len(L)]

  for (i in seq_len(L)) {
    rows <- which(thr$item == i)
    G <- Gs[[i]]
    if (i < L) B[rows, i] <- G[, 1] else B[rows, seq_len(L - 1L)] <- -G[, 1]
    if (extra[i] > 0L)
      B[rows, off[i] + seq_len(extra[i])] <- G[, keep[[i]], drop = FALSE]
  }

  st <- .start_tau(X, thr)
  loc <- vapply(seq_len(L), function(i) mean(st[thr$item == i]), 0)
  loc <- loc - mean(loc)
  beta0 <- c(loc[-L], numeric(sum(extra)))

  inames <- if (is.null(colnames(X))) paste0("V", seq_len(L)) else colnames(X)
  pairs <- .pair_counts(X, m)
  .pcml_check_connected(pairs, L, inames)
  sol <- .pcml_solve(X, thr, m, B, beta0, maxit = maxit, tol = tol,
                     pairs = pairs)
  thr$tau <- sol$tau; thr$se <- sol$se_tau

  labs <- c("spread", "skewness", "kurtosis")
  comp <- data.frame(item = colnames(X),
                     location = c(sol$beta[seq_len(L - 1L)],
                                  -sum(sol$beta[seq_len(L - 1L)])),
                     location_se = NA_real_,
                     spread = NA_real_, spread_se = NA_real_,
                     skewness = NA_real_, skewness_se = NA_real_,
                     kurtosis = NA_real_, kurtosis_se = NA_real_)
  Bloc <- rbind(diag(L - 1L), rep(-1, L - 1L))
  cov_loc <- Bloc %*% sol$cov_beta[seq_len(L - 1L), seq_len(L - 1L), drop = FALSE] %*% t(Bloc)
  comp$location_se <- sqrt(pmax(diag(cov_loc), 0))
  for (i in seq_len(L)) if (extra[i] > 0L) {
    cols <- off[i] + seq_len(extra[i])
    for (j in seq_along(keep[[i]])) {
      lab <- labs[keep[[i]][j] - 1L]
      comp[[lab]][i] <- sol$beta[cols[j]]
      comp[[paste0(lab, "_se")]][i] <- sqrt(pmax(sol$cov_beta[cols[j], cols[j]], 0))
    }
  }

  list(model = "PCM", n_components = n_components, thr = thr,
       components = comp, cov_tau = sol$cov_tau, loglik = sol$loglik,
       iterations = sol$iterations, converged = sol$converged, m = m,
       anchors = NULL, n_parameters = ncol(B), B = B,
       cov_beta = sol$cov_beta, H_beta = sol$H_beta)
}
