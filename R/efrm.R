# RaschR :: extended frame of reference model
# ===========================================================================
# Humphry's extended frame of reference model (Humphry 2005; Humphry &
# Andrich 2008). A frame F_sg is one item-set by person-group cell, with
# unit rho_sg = alpha_s * phi_g:
#
#   P(X_ni = x) prop exp( rho_sg * ( x*theta_n - sum_{h<=x} delta_ih ) )
#
# Within a frame all curves are parallel, so the partial credit model holds
# in the frame's natural unit and the pairwise conditional logic of
# Andrich & Luo (2003) applies unchanged within frames.
#
# Humphry (2005) states the model for dichotomous responses; the polytomous
# form above is this package's extension of it. It is characterised by the
# properties the model's logic requires: the unit multiplies the whole
# exponent, so (i) within every frame the partial credit model holds in the
# frame's natural unit (natural thresholds rho*(delta - c), parallel
# curves), which is what makes the pairwise conditional cancellation valid;
# (ii) the weighted score sum_i rho_i x_i remains sufficient for the person
# parameter; and (iii) it reduces exactly to the dichotomous statement when
# every item has two categories and to the ordinary PCM when every unit is
# one. Per-threshold discriminations would destroy (i) and with it Rasch
# measurement; a unit on only part of the exponent would destroy the
# natural-unit reading. A consequence worth noting when interpreting
# results: category widths in natural units scale with the frame unit.
#
# Estimation is in two routes, each assigned to the data structure that
# identifies it:
#
# Route 1 (person-free, within-frame pairwise conditional ML): the centred,
# alpha-absorbed set thresholds dtilde (sum-zero per set) and the person
# group units phi_g, identified because a set taken by two groups shows the
# same threshold pattern at two scales. The structural map
# tau_v = phi_g(v) * dtilde is bilinear, so estimation alternates two exact
# linear Newton steps and finishes with a joint polish for exact sandwich
# standard errors. Frame origins are NOT pairwise-identifiable (an additive
# constant per frame cancels at every pair total) and item-set units
# alpha_s are exactly confounded with the spread of the set's thresholds,
# so neither is a parameter of this stage.
#
# Route 2 (person-side linking): alpha_s and set locations mu_s from
# persons common to two sets. Within set s a person's frame estimate is
# u = alpha_s * (theta - mu_s), so the ratio of error-corrected true-score
# standard deviations over common persons estimates alpha ratios, and the
# offsets give the relative set locations; the log-ratios are reconciled
# over the linking graph by weighted least squares. The error-variance
# correction var_true = var(u_hat) - mean(se^2) removes the attenuation of
# the naive standard-deviation ratio.
#
# Identification: sum dtilde = 0 per set; sum_g log phi_g = 0;
# sum_s log alpha_s = 0; mean_s mu_s = 0. Together these give the
# product-of-units constraint over the frame grid with the arbitrary-unit
# origin at the mean set location. Cross-set item pairs carry no usable
# conditional information when units differ (the person parameter does not
# cancel) and are excluded from the pairwise stage.
# ===========================================================================

# Same-set pair filter: .pair_counts only returns pairs with overlapping
# persons, so cross-group pairs are already absent; cross-set pairs within a
# group must be removed because the person parameter does not cancel there.
.efrm_filter_pairs <- function(pairs, vmap) {
  Filter(function(pc) vmap$set[pc$i] == vmap$set[pc$j], pairs)
}

# Connected components by union-find; x is a list of integer edge pairs.
.efrm_components <- function(n, edges) {
  parent <- seq_len(n)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  for (e in edges) {
    ra <- find(e[1]); rb <- find(e[2])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(n), find, 1L)
}

# Stage 1: block-coordinate Newton on the bilinear map tau = phi_g * dtilde,
# with a joint polish step and Godambe sandwich covariance.
.efrm_solve <- function(Xv, thr_v, m_v, vmap, pairs, drow, A_D,
                        maxit = 50, tol = 1e-7) {
  Mv <- nrow(thr_v)
  glevs <- sort(unique(vmap$group))
  G <- length(glevs)
  gidx <- match(vmap$group[thr_v$item], glevs)   # group of each virtual threshold
  Pd <- ncol(A_D)

  # one inner Newton block on a linear design tau = off + B beta
  newton_block <- function(B, off, beta, ll_prev, positive = FALSE, inner = 4) {
    glh <- .pcml_glh(drop(off + B %*% beta), thr_v, pairs, m_v)
    for (it in seq_len(inner)) {
      gb <- drop(crossprod(B, glh$g))
      Hb <- crossprod(B, glh$H %*% B)
      step <- tryCatch(solve(Hb, gb), error = function(e)
        solve(Hb - diag(1e-8, nrow(Hb)), gb))
      lam <- 1; ok <- FALSE; g2 <- glh
      for (half in 1:30) {
        cand <- beta - lam * step
        if (positive && any(cand <= 0)) { lam <- lam / 2; next }
        g2 <- .pcml_glh(drop(off + B %*% cand), thr_v, pairs, m_v)
        if (is.finite(g2$ll) && g2$ll >= glh$ll - 1e-12) { ok <- TRUE; break }
        lam <- lam / 2
      }
      if (!ok) break
      beta <- cand; glh <- g2
      if (max(abs(lam * step)) < tol) break
    }
    list(beta = beta, ll = glh$ll, glh = glh)
  }

  # start values: log-ratio least squares pooled over groups, phi = 1
  st <- .start_tau(Xv, thr_v)
  dt0 <- vapply(seq_len(max(drow)), function(d) mean(st[drow == d]), 0)
  # recentre per set
  for (s in unique(vmap$set)) {
    rows <- unique(drow[vmap$set[thr_v$item] == s])
    dt0[rows] <- dt0[rows] - mean(dt0[rows])
  }
  beta_d <- qr.coef(qr(A_D), dt0); beta_d[is.na(beta_d)] <- 0
  phi <- rep(1, G)

  ll <- -Inf
  for (outer in seq_len(maxit)) {
    # delta step: tau = B_d beta_d with B_d = phi_g(v) * A_D[drow(v), ]
    B_d <- A_D[drow, , drop = FALSE] * phi[gidx]
    res_d <- newton_block(B_d, 0, beta_d, ll)
    beta_d <- res_d$beta
    dtil <- drop(A_D %*% beta_d)
    # phi step (skip when G = 1): tau = off + B_r phi[-1]
    if (G > 1L) {
      off <- dtil[drow] * (gidx == 1L)
      B_r <- matrix(0, Mv, G - 1L)
      for (g in 2:G) B_r[gidx == g, g - 1L] <- dtil[drow][gidx == g]
      res_r <- newton_block(B_r, off, phi[-1], res_d$ll, positive = TRUE)
      phi <- c(1, res_r$beta)
      ll_new <- res_r$ll
    } else ll_new <- res_d$ll
    if (is.finite(ll) && abs(ll_new - ll) < tol * (abs(ll_new) + 1)) {
      ll <- ll_new; break
    }
    ll <- ll_new
  }

  # joint polish with the exact bilinear Hessian, then sandwich covariance
  tau_hat <- phi[gidx] * dtil[drow]
  glh <- .pcml_glh(tau_hat, thr_v, pairs, m_v)
  build_J <- function(dtil, phi) {
    J <- matrix(0, Mv, Pd + G - 1L)
    J[, seq_len(Pd)] <- A_D[drow, , drop = FALSE] * phi[gidx]
    if (G > 1L) for (g in 2:G)
      J[gidx == g, Pd + g - 1L] <- dtil[drow][gidx == g]
    J
  }
  struct_H <- function(J, glh) {
    H <- crossprod(J, glh$H %*% J)
    if (G > 1L) for (g in 2:G) {
      rows <- which(gidx == g)
      cb <- drop(crossprod(A_D[drow[rows], , drop = FALSE], glh$g[rows]))
      H[seq_len(Pd), Pd + g - 1L] <- H[seq_len(Pd), Pd + g - 1L] + cb
      H[Pd + g - 1L, seq_len(Pd)] <- H[Pd + g - 1L, seq_len(Pd)] + cb
    }
    H
  }
  for (polish in 1:10) {
    J <- build_J(dtil, phi)
    H <- struct_H(J, glh)
    g_full <- drop(crossprod(J, glh$g))
    step <- tryCatch(solve(H, g_full), error = function(e)
      solve(H - diag(1e-8, nrow(H)), g_full))
    cand_d <- beta_d - step[seq_len(Pd)]
    cand_p <- if (G > 1L) phi[-1] - step[Pd + seq_len(G - 1L)] else numeric(0)
    if (G > 1L && any(cand_p <= 0)) break
    dt_c <- drop(A_D %*% cand_d); phi_c <- c(1, cand_p)
    glh_c <- .pcml_glh(phi_c[gidx] * dt_c[drow], thr_v, pairs, m_v)
    if (!is.finite(glh_c$ll) || glh_c$ll < glh$ll - 1e-10) break
    beta_d <- cand_d; dtil <- dt_c; phi <- phi_c; glh <- glh_c
    if (max(abs(step)) < tol) break
  }
  tau_hat <- phi[gidx] * dtil[drow]

  J <- build_J(dtil, phi)
  H <- struct_H(J, glh)
  Hinv <- tryCatch(solve(H), error = function(e) solve(H - diag(1e-8, nrow(H))))
  Jt <- .pcml_sandwich(Xv, thr_v, m_v, tau_hat, pairs)
  covb <- Hinv %*% crossprod(J, Jt %*% J) %*% Hinv
  conv <- max(abs(drop(crossprod(J, glh$g)))) < 1e-3

  # recentre to sum_g log phi = 0; dtilde absorbs the constant
  cc <- mean(log(phi))
  phi_c <- phi * exp(-cc); dtil_c <- dtil * exp(cc)
  # covariance of log phi (ref group has zero variance), centred
  cov_lp <- matrix(0, G, G)
  if (G > 1L) {
    idx <- Pd + seq_len(G - 1L)
    Dl <- diag(1 / phi[-1], G - 1L)
    cov_lp[2:G, 2:G] <- Dl %*% covb[idx, idx, drop = FALSE] %*% Dl
  }
  Ac <- diag(G) - matrix(1 / G, G, G)
  cov_lp <- Ac %*% cov_lp %*% t(Ac)
  cov_dt <- A_D %*% covb[seq_len(Pd), seq_len(Pd), drop = FALSE] %*% t(A_D) *
    exp(2 * cc)

  list(dtilde = dtil_c, phi = setNames(phi_c, glevs),
       se_log_phi = setNames(sqrt(pmax(diag(cov_lp), 0)), glevs),
       cov_dtilde = cov_dt, loglik = glh$ll, iterations = outer,
       converged = conv, gidx = gidx)
}

# Stage 2: alpha_s and set locations mu_s from persons common to set pairs,
# with the error-variance correction and jackknife standard errors.
.efrm_link_sets <- function(u_mat, se_mat, sets_u, min_link_persons) {
  S <- ncol(u_mat)
  edges <- list(); ls_est <- ls_var <- off_est <- off_n <- numeric(0)
  for (a in seq_len(S - 1)) for (b in (a + 1):S) {
    ok <- which(is.finite(u_mat[, a]) & is.finite(u_mat[, b]))
    if (length(ok) < min_link_persons) next
    u1 <- u_mat[ok, a]; u2 <- u_mat[ok, b]
    s1 <- se_mat[ok, a]; s2 <- se_mat[ok, b]; n <- length(ok)
    v1 <- var(u1) - mean(s1^2); v2 <- var(u2) - mean(s2^2)
    if (!is.finite(v1) || !is.finite(v2) || v1 <= 0 || v2 <= 0)
      stop("too little true person variance to link sets '", sets_u[a],
           "' and '", sets_u[b], "'")
    ls <- 0.5 * (log(v2) - log(v1))                  # log(alpha_b / alpha_a)
    # leave-one-out jackknife of the log slope
    S1 <- sum(u1); Q1 <- sum(u1^2); E1 <- sum(s1^2)
    S2 <- sum(u2); Q2 <- sum(u2^2); E2 <- sum(s2^2)
    lsi <- vapply(seq_len(n), function(j) {
      n1 <- n - 1L
      va <- (Q1 - u1[j]^2 - (S1 - u1[j])^2 / n1) / (n1 - 1L) - (E1 - s1[j]^2) / n1
      vb <- (Q2 - u2[j]^2 - (S2 - u2[j])^2 / n1) / (n1 - 1L) - (E2 - s2[j]^2) / n1
      if (va <= 0 || vb <= 0) return(NA_real_)
      0.5 * (log(vb) - log(va))
    }, 0)
    lsi <- lsi[is.finite(lsi)]
    jv <- if (length(lsi) > 2) (length(lsi) - 1) / length(lsi) *
      sum((lsi - mean(lsi))^2) else 1
    edges[[length(edges) + 1L]] <- c(a, b)
    ls_est <- c(ls_est, ls); ls_var <- c(ls_var, max(jv, 1e-8))
    off_est <- c(off_est, mean(u2) - exp(ls) * mean(u1))  # alpha_b (mu_a - mu_b)
    off_n <- c(off_n, n)
  }
  if (!length(edges)) stop("no set pairs share enough persons to link the units")
  comp <- .efrm_components(S, edges)
  if (length(unique(comp)) > 1L)
    stop("item sets are not linked by common persons; relative units (alpha) ",
         "are unidentified between: ",
         paste(tapply(sets_u, comp, paste, collapse = "+"), collapse = " | "))

  # weighted least squares for log alpha (sum-zero) over the edge contrasts
  A <- rbind(diag(S - 1L), rep(-1, S - 1L))
  C <- matrix(0, length(edges), S)
  for (e in seq_along(edges)) { C[e, edges[[e]][2]] <- 1; C[e, edges[[e]][1]] <- -1 }
  w <- 1 / ls_var; M <- C %*% A; sw <- sqrt(w)
  la <- drop(A %*% qr.coef(qr(M * sw), ls_est * sw))
  la[is.na(la)] <- 0
  alpha <- exp(la)
  # per-set se(log alpha) from the WLS covariance
  XtX <- crossprod(M * sw)
  cov_la <- A %*% tryCatch(solve(XtX), error = function(e)
    solve(XtX + diag(1e-8, ncol(XtX)))) %*% t(A)
  # set locations: mu_a - mu_b = off / alpha_b per edge; mean-zero constraint
  dmu <- off_est / alpha[vapply(edges, `[`, 1L, 2)]
  Cm <- matrix(0, length(edges), S)
  for (e in seq_along(edges)) { Cm[e, edges[[e]][1]] <- 1; Cm[e, edges[[e]][2]] <- -1 }
  Mm <- Cm %*% A; swm <- sqrt(off_n)
  mu <- drop(A %*% qr.coef(qr(Mm * swm), dmu * swm))
  mu[is.na(mu)] <- 0

  list(alpha = setNames(alpha, sets_u),
       se_log_alpha = setNames(sqrt(pmax(diag(cov_la), 0)), sets_u),
       mu = setNames(mu, sets_u),
       edges = data.frame(set_a = sets_u[vapply(edges, `[`, 1L, 1)],
                          set_b = sets_u[vapply(edges, `[`, 1L, 2)],
                          n = off_n, log_slope = ls_est,
                          se = sqrt(ls_var)))
}

# Person estimation under unequal frame units: the weighted score
# W = sum_i rho_i x_i is sufficient, so roots are solved once per
# missing-data pattern and unique weighted score.
.efrm_person_estimates <- function(X, tau_list, disc) {
  N <- nrow(X)
  obs <- !is.na(X)
  m <- vapply(tau_list, length, 1L)
  pat <- apply(obs, 1, function(z) paste(which(z), collapse = ","))
  theta <- se <- rep(NA_real_, N)
  raw <- rowSums(X, na.rm = TRUE); raw[rowSums(obs) == 0L] <- NA
  max_raw <- as.numeric(obs %*% m)
  Xw <- sweep(X, 2, disc, "*")
  W <- rowSums(Xw, na.rm = TRUE); W[rowSums(obs) == 0L] <- NA
  Wmax <- as.numeric(obs %*% (disc * m))
  extreme <- !is.na(W) & (W <= 1e-12 | W >= Wmax - 1e-12)

  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    if (!length(cols)) next
    sel <- which(pat == key)
    r <- disc[cols]; tl <- tau_list[cols]
    for (Wu in unique(signif(W[sel], 12))) {
      who <- sel[signif(W[sel], 12) == Wu]
      g <- function(th) {
        mo <- lapply(seq_along(cols), function(j)
          item_moments(th, tl[[j]], disc = r[j]))
        E  <- vapply(mo, `[[`, 0, "E");  V <- vapply(mo, `[[`, 0, "V")
        m3 <- vapply(mo, `[[`, 0, "mu3")
        (Wu - sum(r * E)) + sum(r^3 * m3) / (2 * sum(r^2 * V))
      }
      root <- tryCatch(uniroot(g, c(-30, 30), tol = 1e-9)$root,
                       error = function(e) NA_real_)
      theta[who] <- root
      if (!is.na(root)) {
        V <- vapply(seq_along(cols), function(j)
          item_moments(root, tl[[j]], disc = r[j])$V, 0)
        se[who] <- 1 / sqrt(sum(r^2 * V))
      }
    }
  }
  data.frame(n_items = rowSums(obs), raw = raw, max_raw = max_raw,
             weighted_score = W, theta = theta, se = se, extreme = extreme)
}

#' Fit the extended frame of reference model
#'
#' Estimates Humphry's extended frame of reference model, in which the unit
#' of the latent scale differs across frames (item-set by person-group
#' cells): the response model is
#' \code{P(X = x) prop exp(rho_sg (x theta - sum delta))} with
#' \code{rho_sg = alpha_s phi_g}. Within frames the partial credit model
#' holds in the frame's natural unit, so item thresholds and the person
#' group units \code{phi} are estimated by within-frame pairwise conditional
#' maximum likelihood (the person parameter cancels; Andrich and Luo 2003),
#' jointly across frames through the sets shared by several groups. Item-set
#' units \code{alpha} and set locations are then estimated from persons
#' common to pairs of sets, using error-corrected true-score variances, and
#' reconciled over the linking graph by weighted least squares. Everything
#' is reported in a common arbitrary unit, and the returned object is also a
#' full \code{\link{rasch}} fit at the item-by-group level, so the package's
#' diagnostic tables and plots apply.
#'
#' Humphry (2005) states the model for dichotomous responses. The polytomous
#' form fitted here, with the frame unit multiplying the whole exponent over
#' the item's partial-credit thresholds, is this package's extension of that
#' statement. It is the form characterised by preserving the two properties
#' the model's logic rests on: the partial credit model holds within every
#' frame in the frame's natural unit (so the pairwise conditional
#' cancellation remains valid), and the weighted score remains sufficient
#' for the person parameter. It reduces exactly to the dichotomous model
#' when items are scored 0/1 and to the ordinary partial credit model when
#' all units equal one. One interpretive consequence: category widths in
#' natural units scale with the frame unit, so a high-unit frame makes
#' proportionally sharper category distinctions; frame-level fit and the
#' per-frame category curves are where a violation of this would appear.
#'
#' Estimation order: the within-frame pairwise stage establishes the
#' centred set thresholds and the person-group units \code{phi}; the
#' person-side linking stage then establishes the item-set units
#' \code{alpha} and set locations. The reported item parameters and all
#' person measures are computed only after every unit is established: item
#' thresholds are mapped into the common arbitrary unit using \code{alpha}
#' and the set locations, and person measures are weighted-score weighted
#' likelihood estimates evaluated under the final units
#' \code{rho = alpha * phi}. The per-frame person estimates used inside the
#' linking stage are interim quantities for the unit ratios only and are
#' discarded. The within-frame stage needs no re-estimation once
#' \code{alpha} is known, because the pairwise likelihood is invariant to
#' the within-set rescaling that \code{alpha} represents; threshold
#' standard errors are conditional on the estimated units, whose own
#' uncertainty is reported separately in \code{alpha_table} and
#' \code{phi_table}.
#'
#' @param data Persons-by-items data (matrix or data frame, like
#'   \code{\link{rasch}}), plus a person-group column.
#' @param item_sets A named list mapping set names to item-column names, or
#'   a named character vector mapping item names to set names. Items not
#'   mentioned form their own set \code{"(rest)"} when a list is given.
#' @param groups Name of the person-group column in \code{data}, or a vector
#'   with one entry per person.
#' @param id,factors,items,n_groups,adjust_N,na_codes As in
#'   \code{\link{rasch}}.
#' @param maxit,tol Outer iteration cap and convergence tolerance of the
#'   bilinear pairwise stage.
#' @param min_link_persons Minimum number of common persons required for a
#'   set pair to contribute to the unit linking.
#' @return An object of classes \code{"rasch_efrm"} and \code{"rasch"}. In
#'   addition to the standard components (computed over item-by-group
#'   virtual columns with the frame units carried in \code{disc}), it has
#'   \code{frames} (one row per frame: units, origin, pooled fit),
#'   \code{phi_table}, \code{alpha_table}, \code{set_table},
#'   \code{item_arbitrary} and \code{thresholds_arbitrary} (the structural
#'   parameters in the common unit), \code{score_curves} (per-group
#'   score-to-measure curves, replacing the raw-score table),
#'   \code{efrm_vs_rasch} (fit comparison against the equal-unit model on
#'   the same conditional information), and \code{linking} (the linking
#'   evidence).
#' @examples
#' \donttest{
#' set.seed(1); Np <- 400
#' simP <- function(th, tau, r) { x <- 0:length(tau)
#'   p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
#' grp <- rep(c("A", "B"), each = Np / 2)
#' phi <- c(A = 0.8, B = 1.25)
#' d <- seq(-1.5, 1.5, length.out = 10)
#' X <- sapply(seq_along(d), function(i) sapply(seq_len(Np), function(n)
#'   sample(0:1, 1, prob = simP(rnorm(Np)[n], d[i], phi[grp[n]]))))
#' colnames(X) <- sprintf("I%02d", seq_along(d))
#' fit <- rasch_efrm(data.frame(X, grp = grp), item_sets = list(core = colnames(X)),
#'                   groups = "grp")
#' fit$phi_table
#' }
#' @export
rasch_efrm <- function(data, item_sets, groups, id = NULL, factors = NULL,
                       items = NULL, n_groups = 10, adjust_N = NA,
                       na_codes = -1, maxit = 50, tol = 1e-7,
                       min_link_persons = 30) {
  # --- roles ----------------------------------------------------------------
  id_vec <- NULL; fac_df <- NULL; grp <- NULL; grp_name <- "group"
  if (is.data.frame(data)) {
    nm <- names(data)
    if (is.character(groups) && length(groups) == 1L && groups %in% nm) {
      grp <- data[[groups]]; grp_name <- groups
    }
    if (is.character(id) && length(id) == 1L && id %in% nm) id_vec <- data[[id]]
    else if (!is.null(id) && length(id) == nrow(data)) id_vec <- id
    if (is.character(factors) && all(factors %in% nm))
      fac_df <- data[, factors, drop = FALSE]
    else if (is.data.frame(factors) && nrow(factors) == nrow(data)) fac_df <- factors
    drop_cols <- c(if (is.character(id)) id else NULL,
                   if (is.character(factors)) factors else NULL,
                   if (is.character(groups) && length(groups) == 1L) groups else NULL)
    item_cols <- if (!is.null(items)) intersect(items, nm) else setdiff(nm, drop_cols)
    X <- as.matrix(data[, item_cols, drop = FALSE])
  } else {
    X <- as.matrix(data)
  }
  if (is.null(grp)) {
    if (length(groups) != nrow(X))
      stop("'groups' must name a column of data or give one value per person")
    grp <- groups
  }
  grp <- factor(grp)
  if (nlevels(grp) < 1L) stop("no person groups found")
  if (is.null(id_vec)) id_vec <- seq_len(nrow(X))

  prep <- .prepare_X(X, na_codes = na_codes); X <- prep$X
  notes <- prep$notes
  m_item <- apply(X, 2, max, na.rm = TRUE); L <- ncol(X)

  # --- item sets --------------------------------------------------------------
  if (is.list(item_sets)) {
    set_of <- rep(NA_character_, L); names(set_of) <- colnames(X)
    for (s in names(item_sets)) {
      hit <- intersect(item_sets[[s]], colnames(X))
      set_of[hit] <- s
    }
    if (anyNA(set_of)) set_of[is.na(set_of)] <- "(rest)"
  } else {
    if (is.null(names(item_sets))) stop("item_sets must be a named list or named vector")
    set_of <- as.character(item_sets)[match(colnames(X), names(item_sets))]
    if (anyNA(set_of))
      stop("item(s) missing from the item_sets map: ",
           paste(colnames(X)[is.na(set_of)], collapse = ", "))
    names(set_of) <- colnames(X)
  }
  sets_u <- sort(unique(set_of)); S <- length(sets_u)
  glevs <- levels(grp); G <- length(glevs)

  # --- virtual columns: item x group ------------------------------------------
  vmap <- expand.grid(item = colnames(X), group = glevs,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  vmap$set <- set_of[vmap$item]
  vmap$vkey <- paste(vmap$item, vmap$group, sep = ":")
  Xv <- matrix(NA_integer_, nrow(X), nrow(vmap),
               dimnames = list(NULL, vmap$vkey))
  for (v in seq_len(nrow(vmap))) {
    sel <- which(as.character(grp) == vmap$group[v])
    Xv[sel, v] <- X[sel, vmap$item[v]]
  }
  keep <- colSums(!is.na(Xv)) > 0L
  Xv <- Xv[, keep, drop = FALSE]; vmap <- vmap[keep, , drop = FALSE]
  rownames(vmap) <- NULL
  m_v <- m_item[vmap$item]
  thr_v <- threshold_index(m_v)

  # delta enumeration over underlying items (ordered by set, then item)
  iord <- order(match(set_of, sets_u), colnames(X))
  items_o <- colnames(X)[iord]
  thr_items <- threshold_index(m_item[items_o])
  Md <- nrow(thr_items)
  # map each virtual threshold row to its delta row
  drow <- vapply(seq_len(nrow(thr_v)), function(r) {
    it <- vmap$item[thr_v$item[r]]
    which(thr_items$item == match(it, items_o) & thr_items$k == thr_v$k[r])
  }, 1L)
  # per-set sum-zero blocks for dtilde
  set_of_drow <- set_of[items_o][thr_items$item]
  A_D <- matrix(0, Md, Md - S)
  cursor <- 0L
  for (s in sets_u) {
    rows <- which(set_of_drow == s); Ms <- length(rows)
    if (Ms < 2L) stop("set '", s, "' needs at least two thresholds")
    A_D[rows, cursor + seq_len(Ms - 1L)] <- rbind(diag(Ms - 1L), rep(-1, Ms - 1L))
    cursor <- cursor + Ms - 1L
  }

  # --- pairwise stage ----------------------------------------------------------
  pairs <- .efrm_filter_pairs(.pair_counts(Xv, m_v), vmap)
  if (!length(pairs)) stop("no informative within-frame item pairs")

  # phi-link check: groups joined when a set has pairs in both groups' frames
  frames_p <- unique(data.frame(set = vmap$set[vapply(pairs, `[[`, 1L, "i")],
                                group = vmap$group[vapply(pairs, `[[`, 1L, "i")]))
  edges_g <- list()
  for (s in unique(frames_p$set)) {
    gs <- match(frames_p$group[frames_p$set == s], glevs)
    if (length(gs) > 1L) for (j in 2:length(gs))
      edges_g[[length(edges_g) + 1L]] <- c(gs[1], gs[j])
  }
  comp_g <- .efrm_components(G, edges_g)
  if (length(unique(comp_g)) > 1L)
    stop("person groups are not linked by any common item set; relative units ",
         "(phi) are unidentified between: ",
         paste(tapply(glevs, comp_g, paste, collapse = "+"), collapse = " | "))

  sol <- .efrm_solve(Xv, thr_v, m_v, vmap, pairs, drow, A_D,
                     maxit = maxit, tol = tol)
  phi <- sol$phi; dtil <- sol$dtilde

  # equal-unit comparison on the same conditional information
  B0 <- A_D[drow, , drop = FALSE]
  bd0 <- qr.coef(qr(A_D), dtil); bd0[is.na(bd0)] <- 0
  glh0 <- .pcml_glh(drop(B0 %*% bd0), thr_v, pairs, m_v)
  for (it0 in 1:25) {
    gb <- drop(crossprod(B0, glh0$g)); Hb <- crossprod(B0, glh0$H %*% B0)
    step <- tryCatch(solve(Hb, gb), error = function(e)
      solve(Hb - diag(1e-8, nrow(Hb)), gb))
    lam <- 1; moved <- FALSE
    for (half in 1:30) {
      cand <- bd0 - lam * step
      g2 <- .pcml_glh(drop(B0 %*% cand), thr_v, pairs, m_v)
      if (is.finite(g2$ll) && g2$ll >= glh0$ll - 1e-12) {
        bd0 <- cand; glh0 <- g2; moved <- TRUE; break
      }
      lam <- lam / 2
    }
    if (!moved || max(abs(lam * step)) < 1e-7) break
  }

  # --- person-side linking (alpha, mu) ----------------------------------------
  if (S > 1L) {
    u_mat <- se_mat <- matrix(NA_real_, nrow(X), S,
                              dimnames = list(NULL, sets_u))
    for (si in seq_len(S)) for (g in glevs) {
      cols <- which(vmap$set == sets_u[si] & vmap$group == g)
      if (!length(cols)) next
      # thresholds in dtilde units for these virtual columns
      tl <- lapply(cols, function(v) {
        rows <- which(thr_v$item == v)
        dtil[drow[rows]]
      })
      pe <- .person_estimates(Xv[, cols, drop = FALSE], tl, disc = phi[g])
      sel <- which(pe$n_items > 0 & !pe$extreme)
      u_mat[sel, si] <- pe$theta[sel]; se_mat[sel, si] <- pe$se[sel]
    }
    link <- .efrm_link_sets(u_mat, se_mat, sets_u, min_link_persons)
    alpha <- link$alpha; mu <- link$mu
  } else {
    alpha <- setNames(1, sets_u); mu <- setNames(0, sets_u)
    link <- list(alpha = alpha, se_log_alpha = setNames(0, sets_u), mu = mu,
                 edges = data.frame())
  }

  # --- assembly in arbitrary units ----------------------------------------------
  delta <- dtil / alpha[set_of_drow] + mu[set_of_drow]
  rho_v <- alpha[vmap$set] * phi[vmap$group]
  thr_v$tau <- delta[drow]
  cov_tau <- sol$cov_dtilde[drow, drow, drop = FALSE] /
    tcrossprod(alpha[set_of_drow][drow])
  thr_v$se <- sqrt(pmax(diag(cov_tau), 0))
  thr_v$anchored <- FALSE
  est <- list(model = "EFRM", thr = thr_v, cov_tau = cov_tau,
              loglik = sol$loglik, iterations = sol$iterations,
              converged = sol$converged, m = m_v, anchors = NULL,
              n_parameters = (Md - S) + (G - 1L) + 2L * (S - 1L))

  fac_all <- data.frame(g = as.character(grp), stringsAsFactors = FALSE)
  names(fac_all) <- grp_name
  if (!is.null(fac_df)) fac_all <- cbind(fac_all, fac_df)
  fit <- .assemble_fit("EFRM", Xv, est, id_vec, fac_all, n_groups, adjust_N,
                       notes, disc = rho_v)

  # --- structural tables -----------------------------------------------------------
  fit$phi_table <- data.frame(group = glevs, phi = unname(phi),
                              se_log_phi = unname(sol$se_log_phi))
  fit$alpha_table <- data.frame(set = sets_u, alpha = unname(alpha),
                                se_log_alpha = unname(link$se_log_alpha))
  fit$set_table <- data.frame(set = sets_u, mu = unname(mu),
                              alpha = unname(alpha),
                              n_items = as.integer(table(set_of)[sets_u]))
  fit$thresholds_arbitrary <- data.frame(item = items_o[thr_items$item],
                                         set = set_of_drow,
                                         k = thr_items$k, delta = delta,
                                         se = sqrt(pmax(diag(
                                           sol$cov_dtilde / tcrossprod(alpha[set_of_drow])), 0)))
  fit$item_arbitrary <- do.call(rbind, lapply(seq_along(items_o), function(i) {
    rows <- which(thr_items$item == i)
    data.frame(item = items_o[i], set = set_of[items_o[i]],
               location = mean(delta[rows]),
               se = sqrt(mean((sol$cov_dtilde[rows, rows, drop = FALSE] /
                                 tcrossprod(alpha[set_of_drow][rows])))))
  }))
  rownames(fit$item_arbitrary) <- NULL

  # frame table with pooled fit
  fr <- unique(vmap[, c("set", "group")])
  fit$frames <- do.call(rbind, lapply(seq_len(nrow(fr)), function(j) {
    cols <- which(vmap$set == fr$set[j] & vmap$group == fr$group[j])
    npers <- sum(rowSums(!is.na(Xv[, cols, drop = FALSE])) > 0)
    gf <- .group_col_fit(fit$residuals, fit$moments, cols, disc = rho_v)
    data.frame(set = fr$set[j], group = fr$group[j],
               n_persons = npers, n_items = length(cols),
               alpha = unname(alpha[fr$set[j]]), phi = unname(phi[fr$group[j]]),
               rho = unname(alpha[fr$set[j]] * phi[fr$group[j]]),
               se_log_rho = sqrt(unname(link$se_log_alpha[fr$set[j]])^2 +
                                   unname(sol$se_log_phi[fr$group[j]])^2),
               origin = unname(mu[fr$set[j]]),
               infit_ms = gf$infit_ms, outfit_ms = gf$outfit_ms,
               fit_resid = gf$fit_resid, n_responses = gf$n)
  }))
  rownames(fit$frames) <- NULL

  # equal-unit comparison and score curves
  fit$efrm_vs_rasch <- list(ll_efrm = sol$loglik, ll_equal = glh0$ll,
                            two_delta_ll = 2 * (sol$loglik - glh0$ll),
                            extra_parameters = G - 1L)
  grid <- seq(-6, 6, by = 0.1)
  fit$score_curves <- do.call(rbind, lapply(glevs, function(g) {
    r_i <- alpha[set_of] * phi[g]
    ew <- vapply(grid, function(th) sum(vapply(seq_len(L), function(i)
      r_i[i] * item_moments(th, delta[thr_items$item == match(colnames(X)[i], items_o)],
                            disc = r_i[i])$E, 0)), 0)
    info <- vapply(grid, function(th) sum(vapply(seq_len(L), function(i)
      r_i[i]^2 * item_moments(th, delta[thr_items$item == match(colnames(X)[i], items_o)],
                              disc = r_i[i])$V, 0)), 0)
    data.frame(group = g, theta = grid, expected_score = ew, sem = 1 / sqrt(info))
  }))
  fit$linking <- list(phi_edges = edges_g, alpha_edges = link$edges)
  fit$virtual_map <- vmap
  fit$set_of <- set_of
  class(fit) <- c("rasch_efrm", "rasch")
  fit
}

#' @export
print.rasch_efrm <- function(x, ...) {
  cat(sprintf("RaschR extended frame of reference analysis: %d items in %d set(s) x %d group(s) = %d frames, %d persons\n",
              length(x$set_of), nrow(x$alpha_table), nrow(x$phi_table),
              nrow(x$frames), nrow(x$X)))
  cat(sprintf("Within-frame pairwise conditional ML: %s in %d iterations\n",
              if (x$est$converged) "converged" else "NOT converged",
              x$est$iterations))
  cat(sprintf("PSI %.3f, power of fit: %s\n", x$psi$PSI, x$power_of_fit))
  cat("\nPerson group units (phi):\n")
  print(x$phi_table, digits = 3, row.names = FALSE)
  cat("\nItem set units (alpha) and locations:\n")
  print(merge(x$alpha_table, x$set_table[, c("set", "mu", "n_items")], by = "set"),
        digits = 3, row.names = FALSE)
  cat(sprintf("\nEqual-unit comparison: 2(ll_EFRM - ll_equal) = %.3f with %d extra unit parameter(s)\n",
              x$efrm_vs_rasch$two_delta_ll, x$efrm_vs_rasch$extra_parameters))
  cat("(composite likelihood: descriptive, not a calibrated chi-square)\n")
  if (length(x$notes)) cat("\nNotes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

#' Plot frame units
#'
#' Caterpillar plot of the frame units \code{rho_sg = alpha_s phi_g} on the
#' log scale, grouped by item set and coloured by person group, with 95 per
#' cent error bars; frames with pooled fit residuals beyond the band are
#' highlighted.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param band Pooled fit residual band beyond which a frame is highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?rasch_efrm for a complete simulated example
#' }
#' @export
plot_frames <- function(fit, band = 2.5) {
  if (!inherits(fit, "rasch_efrm")) stop("plot_frames needs a rasch_efrm fit")
  fr <- fit$frames[order(fit$frames$set, fit$frames$rho), ]
  n <- nrow(fr)
  lr <- log(fr$rho)
  lo <- lr - 1.96 * fr$se_log_rho; hi <- lr + 1.96 * fr$se_log_rho
  labs <- paste0(fr$set, " \u00d7 ", fr$group)
  glev <- sort(unique(fr$group))
  colr <- .rr$pal[(match(fr$group, glev) - 1L) %% length(.rr$pal) + 1L]
  op <- par(mar = c(4.2, 9, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = range(c(lo, hi, 0)) + c(-0.1, 0.1), ylim = c(0.5, n + 0.5),
       xlab = "log unit (log rho)", ylab = "", axes = FALSE, main = "")
  title(main = "Frame units", adj = 0, line = 1.4)
  abline(h = seq_len(n), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  axis(1, col = .rr$grid, col.ticks = .rr$soft)
  axis(2, at = seq_len(n), labels = labs, cex.axis = 0.75,
       col = .rr$grid, col.ticks = NA)
  misfit <- !is.na(fr$fit_resid) & abs(fr$fit_resid) > band
  segments(lo, seq_len(n), hi, seq_len(n), lwd = 2.2,
           col = ifelse(misfit, .rr$red, .rr$soft))
  points(lr, seq_len(n), pch = 21, cex = 1.5,
         bg = ifelse(misfit, .rr$red, colr), col = "white", lwd = 1.2)
  .rr_legend("bottomright", glev, pch = 21,
             pt.bg = .rr$pal[seq_along(glev)], col = "white", pt.cex = 1.2)
  if (any(misfit))
    mtext(sprintf("%d frame(s) with |pooled fit residual| > %.1f", sum(misfit), band),
          side = 3, line = 0.2, adj = 0, cex = 0.8, col = .rr$red)
  invisible(NULL)
}

#' Plot an item's characteristic curves across frames
#'
#' The signature display of the extended frame of reference model: the model
#' expected-score curve of one underlying item drawn once per person group
#' (curves fan with the group units), with observed class-interval means per
#' group overlaid.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param item Underlying item name.
#' @param n_groups Number of class intervals for the observed means.
#' @param grid Logit grid.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?rasch_efrm for a complete simulated example
#' }
#' @export
plot_icc_frames <- function(fit, item, n_groups = fit$n_groups,
                            grid = seq(-5, 5, 0.05)) {
  if (!inherits(fit, "rasch_efrm")) stop("plot_icc_frames needs a rasch_efrm fit")
  vm <- fit$virtual_map
  rows <- which(vm$item == item)
  if (!length(rows)) stop("no such item: ", item)
  thr <- fit$thresholds_arbitrary
  tau_i <- thr$delta[thr$item == item]
  mmax <- length(tau_i)
  op <- .rr_canvas(range(grid), c(0, mmax), "Person location (logits, common unit)",
                   "Expected score",
                   paste0("ICC across frames \u2013 ", item))
  on.exit(par(op))
  th_all <- fit$person$theta
  for (j in seq_along(rows)) {
    v <- rows[j]
    rho <- fit$disc[v]
    colr <- .rr$pal[(j - 1L) %% length(.rr$pal) + 1L]
    Ecurve <- vapply(grid, function(t) item_moments(t, tau_i, disc = rho)$E, 0)
    lines(grid, Ecurve, lwd = 2.6, col = colr)
    x <- fit$X[, v]; ok <- !is.na(th_all) & !is.na(x)
    if (sum(ok) >= 2 * n_groups) {
      ci <- cut(rank(th_all[ok], ties.method = "first"),
                min(n_groups, max(2, floor(sum(ok) / 15))), labels = FALSE)
      points(tapply(th_all[ok], ci, mean), tapply(x[ok], ci, mean),
             pch = 21, bg = colr, col = "white", cex = 1.3, lwd = 1)
    }
  }
  .rr_legend("topleft",
             sprintf("%s (rho %.3f)", vm$group[rows], fit$disc[rows]),
             lwd = 2.6, col = .rr$pal[seq_along(rows)])
  invisible(NULL)
}
