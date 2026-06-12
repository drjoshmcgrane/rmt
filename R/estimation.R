# RaschR :: estimation
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
.pcml_solve <- function(X, thr, m, B, beta0, offset = 0, maxit = 60, tol = 1e-8) {
  pairs <- .pair_counts(X, m)
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
       cov_tau = covt, se_tau = sqrt(pmax(diag(covt), 0)),
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
#'   values are handled by pairwise deletion.
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
#'   \code{item}, \code{k}, \code{tau}, \code{se}, \code{anchored}), the
#'   threshold covariance matrix \code{cov_tau}, the pairwise conditional
#'   log-likelihood, the iteration count, a convergence flag, and the
#'   max-score vector \code{m}.
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

    sol <- .pcml_solve(X, thr, m, B, beta0, offset = offset,
                       maxit = maxit, tol = tol)
    thr$tau <- sol$tau; thr$se <- sol$se_tau; thr$se[a_id] <- 0
    thr$anchored <- seq_len(M) %in% a_id | thr$item %in% mean_items
    return(list(model = model, thr = thr, cov_tau = sol$cov_tau,
                loglik = sol$loglik, iterations = sol$iterations,
                converged = sol$converged, m = m, anchors = anchors,
                n_parameters = ncol(B)))
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

  sol <- .pcml_solve(X, thr, m, B, beta0, maxit = maxit, tol = tol)
  # recentre so the mean item location is zero
  loc <- vapply(seq_len(L), function(i) mean(sol$tau[thr$item == i]), 0)
  sol$tau <- sol$tau - mean(loc)
  thr$tau <- sol$tau; thr$se <- sol$se_tau; thr$anchored <- FALSE

  list(model = model, thr = thr, cov_tau = sol$cov_tau,
       loglik = sol$loglik, iterations = sol$iterations,
       converged = sol$converged, m = m, anchors = NULL,
       n_parameters = ncol(B))
}
