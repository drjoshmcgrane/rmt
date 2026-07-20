# rasch :: extended frame of reference for paired comparisons
# ===========================================================================
# The extended frame of reference model of Humphry (2005) and Humphry and
# Andrich (2008) states that the unit of the latent scale is a property of
# the frame of measurement, not a universal constant. rasch_efrm() fits that
# model for the persons-by-items design; btl_efrm() is this package's
# extension of the same idea to the Bradley-Terry-Luce family of paired
# comparisons (Bradley and Terry 1952; Luce 1959), where the frame is a
# judge-panel by object-set cell.
#
# Objects k are partitioned into sets s(k) (S sets), and judges j into panels
# g(j) (G panels). Each object has a true common-scale value
#
#     v_k = alpha_{s(k)} * beta_k + kappa_{s(k)},
#
# with beta_k the within-set (frame-unit) calibration location, alpha_s > 0
# the set unit, and kappa_s the set origin. A comparison judged in panel g
# carries the panel unit phi_g, and:
#
#   * WITHIN a set s (both objects in s):
#         logit P(a beats b) = phi_g * (beta_a - beta_b).
#     The origin kappa_s cancels, and the set unit alpha_s is confounded with
#     the spread of beta -- exactly as the item-set unit is absorbed in the
#     within-frame stage of the Rasch EFRM -- so it is not a parameter here.
#
#   * ACROSS sets (a in A, b in B):
#         logit P(a beats b) = phi_g * (v_a - v_b)
#                            = phi_g * (alpha_A beta_a - alpha_B beta_b
#                                        + kappa_A - kappa_B).
#     The cross-set comparisons place the two sets on one common scale and so
#     identify alpha and kappa.
#
# The lineage of a frame-dependent unit for comparative judgement is
# Thurstone's (1927) varying discriminal dispersion and the varying-precision
# paired-comparison models catalogued by David (1988); the measurement-unit
# reading of it is Humphry's. The paired-comparison form fitted here is this
# package's extension, stated for dichotomous winner data.
#
# Estimation is in two conditional stages, mirroring rasch_efrm():
#
# Stage 1 (within frames): for each set the within-set comparisons, pooled
# over panels, fit the bilinear model logit = rho_{gs} (b_a - b_b) with
# b sum-zero and one reference panel fixed at rho = 1. This is the same
# constrained bilinear maximisation the Rasch EFRM performs in its stage 1;
# the ratios rho_{gs} = phi_g / phi_{ref(s)} estimate the panel units up to
# the set's reference panel, and are reconciled across sets by a
# precision-weighted least squares over the panel-by-set linking graph, then
# normalised to geometric mean one over panels. The reconciled reference-panel
# units put every set's b on the common panel scale, giving beta.
#
# Stage 2 (linking sets): with beta-hat and phi-hat fixed, the cross-set
# comparisons are a low-dimensional maximum likelihood in (log alpha, kappa)
# for the non-reference sets, solved by Newton with the analytic gradient and
# Hessian. The key theoretical point, and the reason this design is worth
# stating, is that the linking uses only comparison OUTCOMES: no distributional
# assumption about the objects is made, so the set units are identified WITHIN
# the conditional (person-free) framework. This is unlike the persons-by-items
# EFRM, whose item-set units are identified only from the person side (their
# distribution), a genuinely distributional step. The paired-comparison design
# supplies its own conditional link and needs no such assumption.
# ===========================================================================

# connected components of an undirected graph on 1..n given a two-column
# integer edge matrix; union-find, as used across the package's linking code
.btlef_components <- function(n, edges) {
  parent <- seq_len(n)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  if (length(edges)) for (r in seq_len(nrow(edges))) {
    ra <- find(edges[r, 1]); rb <- find(edges[r, 2])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(n), find, 1L)
}

# Stage 1 bilinear solve for ONE frame set (or for the pooled single-unit
# model, called with one panel). Objects are indexed 1..K; panel is a
# character vector per comparison; the most-used panel is the reference
# (rho = 1). Parameters are the sum-zero location contrasts and the free
# panels' log discrimination. Fisher scoring with a step-halving line search
# gives the point estimate; a judge-clustered Godambe sandwich gives the
# covariance, exactly as in btl().
.btlef_stage1 <- function(ia, ib, y, panel, jd, K, maxit, tol, rho_fixed = NULL) {
  R <- length(ia)
  pcount <- table(panel)
  present <- names(pcount)
  ref <- present[which.max(as.integer(pcount))]     # reference panel: rho = 1
  free <- if (is.null(rho_fixed)) setdiff(present, ref) else character(0)
  Gf <- length(free)
  pf <- match(panel, free)                           # free-panel index, NA on ref
  rho0 <- if (is.null(rho_fixed)) rep(1, R) else as.numeric(rho_fixed[panel])
  B <- rbind(diag(K - 1L), rep(-1, K - 1L))          # K x (K-1) sum-zero map
  Bd <- B[ia, , drop = FALSE] - B[ib, , drop = FALSE]
  np <- (K - 1L) + Gf

  eval_th <- function(th) {
    bfree <- th[seq_len(K - 1L)]
    beta <- as.numeric(B %*% bfree)
    rho <- rho0
    if (Gf) {
      rr <- exp(th[(K - 1L) + seq_len(Gf)])
      ok <- !is.na(pf); rho[ok] <- rr[pf[ok]]
    }
    d <- beta[ia] - beta[ib]
    p <- plogis(rho * d)
    ll <- sum(ifelse(y == 1, log(pmax(p, 1e-300)), log(pmax(1 - p, 1e-300))))
    list(beta = beta, rho = rho, d = d, p = p, ll = ll)
  }
  design <- function(cur) {
    J <- Bd * cur$rho
    if (Gf) {
      Jr <- matrix(0, R, Gf)
      for (h in seq_len(Gf)) {
        sel <- which(pf == h); Jr[sel, h] <- cur$rho[sel] * cur$d[sel]
      }
      J <- cbind(J, Jr)
    }
    J
  }

  theta <- numeric(np); cur <- eval_th(theta)
  for (it in seq_len(maxit)) {
    J <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
    g <- crossprod(J, u); Fi <- crossprod(J, J * av)
    step <- tryCatch(solve(Fi, g),
                     error = function(e) solve(Fi + diag(1e-8, np), g))
    ms <- max(abs(step))                     # trust-region cap: a near-flat
    if (is.finite(ms) && ms > 5) step <- step * (5 / ms)   # log-rho direction
    lam <- 1; moved <- FALSE                 # cannot run away
    for (half in 1:30) {
      cand <- theta + lam * as.numeric(step); c2 <- eval_th(cand)
      if (is.finite(c2$ll) && c2$ll >= cur$ll - 1e-12) {
        theta <- cand; cur <- c2; moved <- TRUE; break
      }
      lam <- lam / 2
    }
    if (!moved) break
    if (max(abs(lam * step)) < tol) break
  }

  # judge-clustered Godambe sandwich (unclustered when every judge appears once)
  J <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
  Fi <- crossprod(J, J * av)
  bread <- tryCatch(solve(Fi), error = function(e) solve(Fi + diag(1e-8, np)))
  Sr <- J * u
  Sc <- rowsum(Sr, jd)
  cov_theta <- bread %*% crossprod(Sc) %*% bread
  # scale-free convergence: the gradient per comparison, invariant to the
  # number of comparisons -- an absolute threshold flags converged fits as
  # unconverged on large R (and would then misroute them into the screen);
  # a per-observation criterion also stays permissive at a boundary, where
  # the gradient vanishes but a Newton-decrement quadratic would not
  conv <- isTRUE(max(abs(crossprod(J, u))) < 1e-6 * R)

  cov_bb <- B %*% cov_theta[seq_len(K - 1L), seq_len(K - 1L), drop = FALSE] %*% t(B)
  se_beta <- sqrt(pmax(diag(cov_bb), 0))
  rho_p <- if (is.null(rho_fixed)) setNames(rep(1, length(present)), present)
           else rho_fixed[present]
  cov_lrho <- matrix(0, Gf, Gf, dimnames = list(free, free))
  if (Gf) {
    li <- (K - 1L) + seq_len(Gf)
    rho_p[free] <- exp(theta[li])
    cov_lrho <- cov_theta[li, li, drop = FALSE]
    dimnames(cov_lrho) <- list(free, free)
  }
  list(beta = cur$beta, se_beta = se_beta, p = cur$p, ll = cur$ll,
       ref = ref, panels = present, rho = rho_p, free = free,
       cov_lrho = cov_lrho, converged = conv)
}

# Reconcile the per-set panel ratios into one set of panel units phi with
# geometric mean one, by generalised least squares on the panel graph. Each
# set contributes the observations log rho_{gs} = log phi_g - log phi_{ref(s)}
# for its free panels g, with the within-set covariance of those log-ratios
# carried across from stage 1; the observations are independent between sets,
# so the full observation covariance is block-diagonal by set. GLS on this
# gives both the precision-weighted point estimate (the reconciliation over
# the sets where a panel appears) and correctly correlated standard errors.
# Errors informatively when the panels are not connected through shared sets.
.btlef_reconcile_phi <- function(panels_u, blocks) {
  G <- length(panels_u)
  if (G == 1L)
    return(list(phi = setNames(1, panels_u),
                se_log_phi = setNames(NA_real_, panels_u),
                lphi = setNames(0, panels_u)))
  # flatten the per-set blocks into one observation vector, design and
  # block-diagonal covariance
  y <- numeric(0); pan <- ref <- character(0)
  Cov <- matrix(0, 0, 0)
  for (bk in blocks) {
    if (!length(bk$free)) next
    idx <- length(y) + seq_along(bk$free)
    y <- c(y, bk$lrho[bk$free]); pan <- c(pan, bk$free)
    ref <- c(ref, rep(bk$ref, length(bk$free)))
    cb <- bk$cov[bk$free, bk$free, drop = FALSE]
    cb <- cb + diag(1e-10, nrow(cb))                       # numerical floor
    Z <- matrix(0, nrow(Cov) + nrow(cb), ncol(Cov) + ncol(cb))
    if (nrow(Cov)) Z[seq_len(nrow(Cov)), seq_len(ncol(Cov))] <- Cov
    Z[nrow(Cov) + seq_len(nrow(cb)), ncol(Cov) + seq_len(ncol(cb))] <- cb
    Cov <- Z
  }
  if (!length(y))
    stop("the panels cannot be linked: no set contains comparisons from more ",
         "than one panel, so the panel units phi are unidentified")
  ei <- cbind(match(pan, panels_u), match(ref, panels_u))
  comp <- .btlef_components(G, ei)
  if (length(unique(comp)) > 1L)
    stop("the panel-by-set graph is not connected; panel units (phi) are ",
         "unidentified between: ",
         paste(tapply(panels_u, comp, paste, collapse = "+"), collapse = " | "))

  g0 <- panels_u[1]                                        # arbitrary anchor
  cols <- setdiff(panels_u, g0)
  X <- matrix(0, length(y), length(cols))
  for (r in seq_along(y)) {
    cp <- match(pan[r], cols); if (!is.na(cp)) X[r, cp] <- X[r, cp] + 1
    cr <- match(ref[r], cols); if (!is.na(cr)) X[r, cr] <- X[r, cr] - 1
  }
  W <- solve(Cov)                                          # GLS weight
  XtW <- t(X) %*% W
  covred <- solve(XtW %*% X)
  bred <- covred %*% (XtW %*% y)
  lphi <- setNames(numeric(G), panels_u); lphi[cols] <- bred
  cov_full <- matrix(0, G, G, dimnames = list(panels_u, panels_u))
  cov_full[cols, cols] <- covred
  A <- diag(G) - matrix(1 / G, G, G)                       # centre to geo-mean 1
  lphi_c <- as.numeric(A %*% lphi)
  cov_c <- A %*% cov_full %*% t(A)
  list(phi = setNames(exp(lphi_c), panels_u),
       se_log_phi = setNames(sqrt(pmax(diag(cov_c), 0)), panels_u),
       lphi = setNames(lphi_c, panels_u))
}

# Stage 2: cross-set linking. With the frame locations beta and panel units
# phi held fixed, estimate (log alpha, kappa) for the non-reference sets by
# Newton on the cross-set comparison likelihood. Standard errors are the
# inverse observed information, conditional on stage 1 (the stage-1
# uncertainty is not propagated -- see the roxygen note).
.btlef_stage2 <- function(a, b, y, phg, sa, sb, bhat, sets_u, maxit, tol) {
  S <- length(sets_u); free <- sets_u[-1L]; nf <- S - 1L; np <- 2L * nf
  ba <- bhat[a]; bb <- bhat[b]
  fa <- match(sa, free); fb <- match(sb, free)             # NA on the reference set
  R <- length(y)

  eval_th <- function(th) {
    la <- th[seq_len(nf)]; kap <- th[nf + seq_len(nf)]
    alpha <- setNames(rep(1, S), sets_u); kappa <- setNames(rep(0, S), sets_u)
    alpha[free] <- exp(la); kappa[free] <- kap
    va <- alpha[sa] * ba + kappa[sa]; vb <- alpha[sb] * bb + kappa[sb]
    p <- plogis(phg * (va - vb))
    ll <- sum(ifelse(y == 1, log(pmax(p, 1e-300)), log(pmax(1 - p, 1e-300))))
    list(alpha = alpha, kappa = kappa, p = p, ll = ll)
  }
  # derivative design D (R x np): columns log alpha (free sets), then kappa
  design <- function(cur) {
    D <- matrix(0, R, np)
    for (j in seq_len(nf)) {
      s <- free[j]; aj <- cur$alpha[[s]]
      onA <- !is.na(fa) & fa == j; onB <- !is.na(fb) & fb == j
      D[, j] <- phg * (onA * aj * ba - onB * aj * bb)       # d eta / d log alpha_s
      D[, nf + j] <- phg * (onA - onB)                      # d eta / d kappa_s
    }
    D
  }

  theta <- numeric(np); cur <- eval_th(theta)
  for (it in seq_len(maxit)) {
    D <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
    g <- crossprod(D, u); Fi <- crossprod(D, D * av)
    step <- tryCatch(solve(Fi, g),
                     error = function(e) solve(Fi + diag(1e-8, np), g))
    lam <- 1; moved <- FALSE
    for (half in 1:30) {
      cand <- theta + lam * as.numeric(step); c2 <- eval_th(cand)
      if (is.finite(c2$ll) && c2$ll >= cur$ll - 1e-12) {
        theta <- cand; cur <- c2; moved <- TRUE; break
      }
      lam <- lam / 2
    }
    if (!moved) break
    if (max(abs(lam * step)) < tol) break
  }

  # observed information: the la-diagonal carries the curvature d2 eta / d la^2
  D <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
  H <- crossprod(D, D * av)
  for (j in seq_len(nf)) {
    s <- free[j]; aj <- cur$alpha[[s]]
    onA <- !is.na(fa) & fa == j; onB <- !is.na(fb) & fb == j
    curv <- phg * (onA * aj * ba - onB * aj * bb)
    H[j, j] <- H[j, j] - sum(u * curv)
  }
  cov <- tryCatch(solve(H), error = function(e)
    solve(crossprod(D, D * av) + diag(1e-8, np)))
  # scale-free per-comparison gradient criterion (see .btlef_stage1)
  conv <- isTRUE(max(abs(crossprod(D, u))) < 1e-6 * length(y))
  se <- sqrt(pmax(diag(cov), 0))
  list(alpha = cur$alpha, kappa = cur$kappa, p = cur$p, ll = cur$ll,
       cov = cov, se_log_alpha = setNames(se[seq_len(nf)], free),
       se_kappa = setNames(se[nf + seq_len(nf)], free),
       free = free, converged = conv)
}

# pooled log-of-mean-square fit residual over a set of comparisons, using the
# frame model's fitted probabilities (the paired-comparison residual logic of
# btl(): z = (y - p) / sqrt(p(1-p)), Andrich and Marais 2019 ch. 23)
.btlef_frame_fit <- function(y, p) {
  n <- length(y)
  if (n < 3L) return(NA_real_)
  V <- pmax(p * (1 - p), 1e-12)
  z2 <- (y - p)^2 / V
  mu4 <- p * (1 - p)^4 + (1 - p) * p^4
  c4v <- mu4 / V^2 - 1
  y2 <- sum(z2); f <- n; v <- sum(c4v)
  if (v > 1e-8 && y2 > 0) f * (log(y2) - log(f)) / sqrt(v) else NA_real_
}

#' Fit the extended frame of reference model for paired comparisons
#'
#' Estimates object locations from paired comparisons when the unit of the
#' latent scale differs across frames -- judge-panel by object-set cells -- an
#' extension of Humphry's (2005) extended frame of reference model to the
#' Bradley-Terry-Luce family. Objects are partitioned into sets and judges into
#' panels. Each object has a common-scale value \code{v_k = alpha_s beta_k +
#' kappa_s}, with \code{beta_k} the within-set calibration location,
#' \code{alpha_s > 0} the set unit and \code{kappa_s} the set origin; a
#' comparison judged in panel \code{g} carries the panel unit \code{phi_g}.
#' Within a set the comparison logit is \code{phi_g (beta_a - beta_b)} (the set
#' origin cancels and the set unit is confounded with the spread of
#' \code{beta}, so neither is identified within a set); across sets it is
#' \code{phi_g (v_a - v_b)}, which places the sets on one common scale.
#'
#' One modelling convention deserves plain statement. Rewritten with a single
#' latent value per object, the model says within-set-\code{s} contests are
#' judged at discrimination \code{phi_g / alpha_s} while every cross-set
#' contest is judged at exactly \code{phi_g}: the common scale is
#' \emph{defined} as the scale of between-frame judgement. That is an
#' assumption about how discriminal dispersion behaves when unlike objects
#' meet (a Thurstonian question with no assumption-free answer), not a
#' consequence of the theory; the per-frame fit residuals in \code{frames}
#' are where a violation would show.
#'
#' Estimation follows \code{\link{rasch_efrm}} in two conditional stages. In
#' stage one the within-set comparisons of each set, pooled over panels, fit
#' the bilinear model \code{logit = rho_{gs} (b_a - b_b)} with \code{b}
#' sum-zero and the most-used panel fixed at \code{rho = 1}; the ratios
#' \code{rho_{gs} = phi_g / phi_{ref(s)}} estimate the panel units up to each
#' set's reference panel and are reconciled across sets by a precision-weighted
#' least squares over the panel-by-set linking graph, then normalised to
#' geometric mean one. In stage two the cross-set comparisons are a
#' low-dimensional maximum likelihood in \code{(log alpha, kappa)} for the
#' non-reference sets, with \code{alpha_1 = 1} and \code{kappa_1 = 0} fixing
#' the reference set (the first, alphabetically).
#'
#' The set units are identified here WITHIN the conditional framework: the
#' cross-set linking uses only comparison outcomes and makes no distributional
#' assumption about the objects. This is the substantive difference from the
#' persons-by-items EFRM, whose item-set units can only be identified from the
#' person side -- that is, from the distribution of the persons over the linked
#' sets. The paired-comparison design supplies its own conditional link, so no
#' such distributional step is needed. See Humphry (2005) and Humphry and
#' Andrich (2008) for the theory of the unit, and Thurstone (1927) and David
#' (1988) for the varying-discriminal-dispersion lineage from which the
#' frame-dependent unit descends. The paired-comparison form is this package's
#' extension of Humphry's model.
#'
#' The ESTIMATES are always the staged conditional estimator: the within-frame
#' calibrations are invariant to the linking data (the frame-of-reference
#' property), a deliberate trade of some efficiency for invariance, exactly as
#' anchored equating trades. Inference defaults to a parametric bootstrap of
#' the whole pipeline (\code{se_method = "bootstrap"}): winners are resampled
#' from the fitted probabilities and both stages refitted, so the standard
#' errors of \code{phi}, \code{alpha}, \code{kappa}, \code{beta} and the
#' common-scale \code{v} carry every stage's sampling variability -- verified
#' by simulation to restore nominal coverage where the conditional errors
#' cover at a third of their nominal rate on linked designs. The
#' \code{"conditional"} option keeps the fast analytic errors (judge-clustered
#' sandwich for stage one, inverse observed information for stage two,
#' conditional on stage one) for quick inspection; its \code{alpha} and
#' \code{kappa} errors understate, and the fit says so.
#'
#' Two honesty notes on the bootstrap. The default is model-based:
#' replicates are drawn as independent Bernoulli outcomes at the fitted
#' probabilities, which is self-consistent (the model has no judge
#' parameter) but does not carry extra-model dependence within judges. When
#' judges plausibly carry idiosyncratic preferences across their
#' comparisons, use \code{se_method = "judge_bootstrap"}: judges are
#' redrawn with replacement within each panel and the whole pipeline is
#' refitted, so the errors carry both the stage-one uncertainty and the
#' within-judge dependence (it needs enough judges per panel to resample
#' stably; failed replicates are counted and reported). And a parameter that reaches
#' its boundary in some replicates (a set unit driven to zero when a resampled
#' within-set order flips against the cross-set evidence, the signature of a
#' two-object set with a near-even internal pair) has no normal sampling
#' distribution: its standard error is reported as \code{NA} and a note names
#' the parameter and the boundary count rather than manufacturing a number.
#' Relatedly, a set whose within-set contests are all near-even (or all
#' one-sided) carries no stable information about the panel-unit ratios; such
#' sets are screened out of the \code{phi} reconciliation, refit with the
#' panel units held at the reconciled \code{phi} (which the frame model says
#' apply to them regardless), and named in a note.
#'
#' A single set (\code{S = 1}) reduces the model to panel units alone; stage
#' two is skipped and the print states the panel-units model. When
#' additionally \code{G = 1} the fit reduces exactly to \code{\link{btl}} on
#' the same data. The equal-unit (single-unit) comparison refits plain
#' \code{\link{btl}} on all comparisons pooled and reports the descriptive
#' composite log-likelihood difference against the frames model; because that
#' comparison is a composite likelihood, the inference on the units is carried
#' by the Wald tests on \code{log phi_g} and \code{log alpha_s} in
#' \code{phi_table} and \code{alpha_table}.
#'
#' @param data A data frame with one comparison per row.
#' @param object_a,object_b Names of the columns holding the two compared
#'   objects.
#' @param winner Name of the column holding the winner; its value must equal
#'   the row's \code{object_a} or \code{object_b} entry, \code{"tie"} or
#'   \code{"draw"} marks a tie, and anything else is treated as missing.
#' @param judge Name of the judge column (clusters the stage-one standard
#'   errors and defines the panels when \code{panels} is a judge attribute).
#' @param panels Either the name of a judge-attribute column in \code{data} or
#'   a named vector mapping judge to panel.
#' @param object_sets A named list mapping set names to character vectors of
#'   object names; every compared object must belong to exactly one set.
#' @param response Not supported: this first implementation fits dichotomous
#'   winner data only. Supplying it raises an informative error.
#' @param ties \code{"drop"} (default, removed with a note) or \code{"error"}.
#' @param min_link Minimum number of cross-set comparisons a set pair must
#'   supply to be used for linking; sets not reachable from the reference set
#'   through sufficient cross-set pairs raise an error.
#' @param se_method \code{"bootstrap"} (the default): a parametric bootstrap;
#'   \code{"judge_bootstrap"}: judges resampled with replacement within
#'   panels, carrying within-judge dependence;
#'   of the ENTIRE two-stage pipeline -- winners resampled from the fitted
#'   probabilities, both stages refitted \code{boot_reps} times -- so the
#'   reported standard errors carry every source of sampling variability,
#'   including the stage-one uncertainty that flows into the linking.
#'   \code{"conditional"}: the fast analytic errors, exact for \code{beta} and
#'   \code{phi} but conditional on stage one for \code{alpha} and
#'   \code{kappa}, which they therefore UNDERSTATE (verified by simulation);
#'   for quick inspection only.
#' @param boot_reps Bootstrap replicates for \code{se_method = "bootstrap"}.
#' @param maxit,tol Newton iteration cap and convergence tolerance.
#' @return An object of class \code{"rasch_btl_efrm"}: \code{objects} (object,
#'   set, \code{beta_set} and its standard error, common-scale \code{v} and its
#'   standard error), \code{phi_table} (panel units with Wald tests against
#'   \code{log phi = 0}), \code{alpha_table} and \code{kappa_table} (set units
#'   and origins with Wald tests; the reference set carries \code{alpha = 1},
#'   \code{kappa = 0} with no standard error), \code{frames} (panel by set:
#'   unit \code{rho = phi alpha}, comparison count, pooled fit residual),
#'   \code{equal_unit} (the descriptive single-unit comparison), \code{n_cross}
#'   (cross-set comparison counts per set pair), \code{notes} and
#'   \code{converged}.
#' @references Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#'   incomplete block designs: I. The method of paired comparisons.
#'   Biometrika, 39, 324-345.
#'
#'   David, H. A. (1988). The Method of Paired Comparisons (2nd ed.). Griffin.
#'
#'   Humphry, S. M. (2005). Maintaining a common arbitrary unit in social
#'   measurement. PhD thesis, Murdoch University.
#'
#'   Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
#'   Rasch model. Journal of Applied Measurement, 9(3), 249-264.
#'
#'   Luce, R. D. (1959). Individual Choice Behavior. Wiley.
#'
#'   Thurstone, L. L. (1927). A law of comparative judgment. Psychological
#'   Review, 34, 273-286.
#' @examples
#' \donttest{
#' d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
#'                        set_units = c(1, 1.4), set_origins = c(0, 0.8),
#'                        seed = 1)
#' fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
#'                 judge = "judge", panels = "panel",
#'                 object_sets = attr(d, "truth")$object_sets)
#' fit$alpha_table
#' }
#' @export
btl_efrm <- function(data, object_a, object_b, winner, judge, panels,
                     object_sets, response = NULL,
                     ties = c("drop", "error"), min_link = 20,
                     se_method = c("bootstrap", "judge_bootstrap",
                                   "conditional"),
                     boot_reps = 60, maxit = 60, tol = 1e-8) {
  ties <- match.arg(ties)
  se_method <- match.arg(se_method)
  if (!is.null(response))
    stop("btl_efrm fits dichotomous winner data only in this first ",
         "implementation; a graded `response` is not supported. Reduce the ",
         "graded margins to a winner, or use btl() for a single-frame graded ",
         "analysis.")
  data <- as.data.frame(data)
  for (col in c(object_a, object_b, winner, judge))
    if (!col %in% names(data)) stop("column not found: ", col)
  a <- trimws(as.character(data[[object_a]]))
  b <- trimws(as.character(data[[object_b]]))
  wn <- trimws(as.character(data[[winner]]))
  jd <- as.character(data[[judge]])

  # panels: a judge-attribute column, or a named judge -> panel vector
  if (length(panels) == 1L && is.character(panels) && panels %in% names(data)) {
    pan <- as.character(data[[panels]])
  } else if (!is.null(names(panels)) && all(nzchar(names(panels)))) {
    pan <- unname(as.character(panels)[match(jd, names(panels))])
  } else {
    stop("`panels` must name a column of `data` or be a named vector ",
         "mapping judge to panel")
  }
  notes <- character(0)

  keep <- !is.na(a) & !is.na(b) & !is.na(wn) & !is.na(jd) & !is.na(pan) & a != b
  if (any(!keep)) {
    notes <- c(notes, sprintf("%d row(s) dropped (missing or self-comparison)",
                              sum(!keep)))
    a <- a[keep]; b <- b[keep]; wn <- wn[keep]; jd <- jd[keep]; pan <- pan[keep]
  }
  if (!length(a)) stop("no usable comparisons")

  y <- ifelse(wn == a, 1L, ifelse(wn == b, 0L, NA_integer_))
  is_tie <- is.na(y) & tolower(wn) %in% c("tie", "draw")
  miss <- is.na(y) & !is_tie
  if (any(miss)) {
    notes <- c(notes, sprintf(
      "%d row(s) with winner matching neither object treated as missing", sum(miss)))
    sel <- !miss
    a <- a[sel]; b <- b[sel]; y <- y[sel]; jd <- jd[sel]; pan <- pan[sel]
  }
  if (anyNA(y)) {
    nt <- sum(is.na(y))
    if (ties == "error") stop(nt, " tie(s) present; set ties = 'drop'")
    notes <- c(notes, sprintf("%d tie(s) dropped", nt))
    sel <- !is.na(y)
    a <- a[sel]; b <- b[sel]; y <- y[sel]; jd <- jd[sel]; pan <- pan[sel]
  }
  if (!length(a)) stop("no usable comparisons after cleaning")

  # --- object sets ----------------------------------------------------------
  if (!is.list(object_sets) || is.null(names(object_sets)) ||
      any(!nzchar(names(object_sets))))
    stop("`object_sets` must be a named list: set name -> object names")
  objs_all <- sort(unique(c(a, b)))
  set_of <- setNames(rep(NA_character_, length(objs_all)), objs_all)
  multi <- character(0)
  for (s in names(object_sets)) {
    hit <- intersect(as.character(object_sets[[s]]), objs_all)
    for (o in hit) {
      if (!is.na(set_of[o])) multi <- c(multi, o)
      set_of[o] <- s
    }
  }
  if (length(multi))
    stop("object(s) assigned to more than one set: ",
         paste(unique(multi), collapse = ", "))
  if (anyNA(set_of))
    stop("object(s) in the data not found in `object_sets` (every compared ",
         "object must belong to exactly one set): ",
         paste(objs_all[is.na(set_of)], collapse = ", "))
  sets_u <- sort(unique(set_of)); S <- length(sets_u)
  panels_u <- sort(unique(pan)); G <- length(panels_u)
  sa <- set_of[a]; sb <- set_of[b]
  within <- sa == sb

  if (any(table(set_of) < 2L))
    stop("every set needs at least two objects; offending set(s): ",
         paste(names(which(table(set_of) < 2L)), collapse = ", "))

  # --- structural checks that do not depend on the outcomes -----------------
  # panel linkage first (its failure is the more fundamental design flaw):
  # panels are connected when they share within-set comparisons in some set
  if (G > 1L) {
    pe <- unique(do.call(rbind, lapply(sets_u, function(s) {
      pg <- unique(pan[within & sa == s])
      if (length(pg) < 2L) return(NULL)
      t(combn(match(pg, panels_u), 2))
    })))
    pcomp <- .btlef_components(G, if (is.null(pe)) matrix(numeric(0), 0, 2) else pe)
    if (any(pcomp != pcomp[1]))
      stop("the panels cannot be linked: no set contains comparisons from ",
           "more than one panel connecting panel group(s) {",
           paste(panels_u[pcomp != pcomp[1]], collapse = ", "),
           "} to the rest; panel units phi are unidentified across them")
  }
  cross <- which(!within)
  n_cross <- data.frame(set_a = character(0), set_b = character(0),
                        n = integer(0), stringsAsFactors = FALSE)
  if (S > 1L) {
    if (!length(cross))
      stop("no cross-set comparisons: the sets cannot be linked to a common ",
           "scale (set units alpha and origins kappa are unidentified)")
    key <- ifelse(sa[cross] < sb[cross], paste(sa[cross], sb[cross]),
                  paste(sb[cross], sa[cross]))
    tab <- table(key); parts <- do.call(rbind, strsplit(names(tab), " ", fixed = TRUE))
    n_cross <- data.frame(set_a = parts[, 1], set_b = parts[, 2],
                          n = as.integer(tab), stringsAsFactors = FALSE)
    rownames(n_cross) <- NULL
    used <- n_cross$n >= min_link
    edges <- cbind(match(n_cross$set_a[used], sets_u),
                   match(n_cross$set_b[used], sets_u))
    comp <- .btlef_components(S, edges)
    ref_comp <- comp[1]                                 # reference set = sets_u[1]
    if (any(comp != ref_comp))
      stop("set(s) not reachable from the reference set '", sets_u[1],
           "' through cross-set pairs with at least min_link = ", min_link,
           " comparisons: ", paste(sets_u[comp != ref_comp], collapse = ", "),
           " (increase the cross-set data or lower min_link)")
  }

  # --- the staged conditional estimator, callable on any outcome vector -----
  # (one function for the observed data and for every bootstrap replicate, so
  # the resampled pipeline is identical to the reported one)
  fit_once <- function(yy) {
    bhat <- setNames(rep(NA_real_, length(objs_all)), objs_all)
    se_bhat <- bhat
    ref_of_set <- setNames(rep(NA_character_, S), sets_u)
    within_p <- rep(NA_real_, length(a))              # frame-model fitted p
    blocks <- list()                                  # per-set panel-ratio blocks
    ll_within <- 0; s1_conv <- TRUE
    s1 <- list(); dropped <- character(0)
    for (s in sets_u) {
      rows <- which(within & sa == s)
      os <- sort(names(set_of)[set_of == s]); Ks <- length(os)
      ia <- match(a[rows], os); ib <- match(b[rows], os)
      if (anyNA(ia) || anyNA(ib) || any(is.na(match(os, unique(c(a[rows], b[rows]))))))
        stop("set '", s, "' has object(s) with no within-set comparison; ",
             "each object needs at least one comparison inside its own set")
      fit1 <- .btlef_stage1(ia, ib, yy[rows], pan[rows], jd[rows], Ks, maxit, tol)
      s1[[s]] <- list(fit = fit1, rows = rows, os = os, ia = ia, ib = ib)
      # A set carries information about the panel-unit ratios only through its
      # internal separation: when its within-set contests are all near-even
      # (or one-sided), the ratio log rho_gs is the quotient of two near-zero
      # (or unbounded) logits and its estimate diverges with a spurious
      # covariance. Screen such blocks out of the phi reconciliation rather
      # than let them poison it; the set is refit below at the reconciled
      # units, which the frame model says apply to it regardless.
      lr <- log(fit1$rho[fit1$free])
      usable <- isTRUE(fit1$converged) &&
        (!length(fit1$free) ||
           (all(is.finite(lr)) && max(abs(lr)) < 4 &&
            all(is.finite(fit1$cov_lrho)) && all(diag(fit1$cov_lrho) > 0)))
      if (!usable && length(fit1$free)) { dropped <- c(dropped, s); next }
      blocks[[s]] <- list(ref = fit1$ref, free = fit1$free,
                          lrho = setNames(lr, fit1$free),
                          cov = fit1$cov_lrho)
    }
    rec <- tryCatch(.btlef_reconcile_phi(panels_u, blocks), error = function(e) {
      if (length(dropped))
        stop(conditionMessage(e), "\n  (set(s) ", paste(dropped, collapse = ", "),
             " were excluded from the panel-unit reconciliation because their ",
             "within-set comparisons carry no stable panel-ratio information)",
             call. = FALSE)
      stop(e)
    })
    phi <- rec$phi
    for (s in sets_u) {
      rows <- s1[[s]]$rows; os <- s1[[s]]$os
      if (s %in% dropped) {
        # refit the set's locations with the panel units held at the
        # reconciled phi: beta comes out directly on the common scale
        fit1 <- .btlef_stage1(s1[[s]]$ia, s1[[s]]$ib, yy[rows], pan[rows],
                              jd[rows], length(os), maxit, tol,
                              rho_fixed = phi)
        s1[[s]]$fit <- fit1
        bhat[os] <- fit1$beta; se_bhat[os] <- fit1$se_beta
      } else {
        fit1 <- s1[[s]]$fit
        pr <- phi[[fit1$ref]]
        bhat[os] <- fit1$beta / pr; se_bhat[os] <- fit1$se_beta / pr
      }
      within_p[rows] <- fit1$p
      ref_of_set[s] <- fit1$ref
      ll_within <- ll_within + fit1$ll
      s1_conv <- s1_conv && isTRUE(fit1$converged)
    }
    # within-set fitted p on the common scale: logit = phi_g (bhat_a - bhat_b)
    within_p[within] <- plogis(phi[pan[within]] * (bhat[a[within]] - bhat[b[within]]))

    alpha <- setNames(rep(1, S), sets_u); kappa <- setNames(rep(0, S), sets_u)
    se_log_alpha <- setNames(rep(NA_real_, S), sets_u)
    se_kappa <- setNames(rep(NA_real_, S), sets_u)
    cov2 <- NULL; s2_conv <- TRUE; ll_cross <- 0
    p_all <- within_p
    if (S > 1L) {
      st2 <- .btlef_stage2(a[cross], b[cross], yy[cross], phi[pan[cross]],
                           sa[cross], sb[cross], bhat, sets_u, maxit, tol)
      alpha <- st2$alpha; kappa <- st2$kappa; cov2 <- st2$cov
      se_log_alpha[st2$free] <- st2$se_log_alpha
      se_kappa[st2$free] <- st2$se_kappa
      s2_conv <- st2$converged; ll_cross <- st2$ll
      p_all[cross] <- st2$p
    }
    v <- alpha[set_of[objs_all]] * bhat[objs_all] + kappa[set_of[objs_all]]
    list(bhat = bhat, se_bhat = se_bhat, phi = phi,
         se_log_phi = rec$se_log_phi, ref_of_set = ref_of_set,
         alpha = alpha, kappa = kappa, cov2 = cov2,
         se_log_alpha = se_log_alpha, se_kappa = se_kappa, v = v,
         within_p = within_p, p_all = p_all,
         ll_within = ll_within, ll_cross = ll_cross,
         dropped = dropped,
         converged = isTRUE(s1_conv && s2_conv))
  }

  fit0 <- fit_once(y)
  bhat <- fit0$bhat; se_bhat <- fit0$se_bhat
  phi <- fit0$phi; ref_of_set <- fit0$ref_of_set
  alpha <- fit0$alpha; kappa <- fit0$kappa; cov2 <- fit0$cov2
  se_log_phi <- fit0$se_log_phi
  se_log_alpha <- fit0$se_log_alpha; se_kappa <- fit0$se_kappa
  v <- fit0$v; within_p <- fit0$within_p
  ll_within <- fit0$ll_within; ll_cross <- fit0$ll_cross

  # conditional (analytic) delta-method errors for the common-scale values
  se_v <- se_bhat[objs_all]                             # reference set: v = beta
  free <- sets_u[-1L]
  if (S > 1L) for (o in objs_all) {
    s <- set_of[[o]]; if (s == sets_u[1]) next
    j <- match(s, free); idx <- c(j, (S - 1L) + j)
    C2 <- cov2[idx, idx, drop = FALSE]
    gvec <- c(alpha[[s]] * bhat[[o]], 1)                # d v / d(log alpha, kappa)
    var_link <- drop(t(gvec) %*% C2 %*% gvec)
    se_v[[o]] <- sqrt(pmax(alpha[[s]]^2 * se_bhat[[o]]^2 + var_link, 0))
  }

  # --- parametric bootstrap of the whole pipeline (default) -----------------
  # winners resampled from the fitted probabilities, both stages refitted:
  # the SEs then carry stage-one uncertainty into the linking, which the
  # conditional errors omit (and demonstrably understate)
  boot_fail <- 0L
  if (se_method == "judge_bootstrap") {
    # resample JUDGES with replacement within each panel (the panel design
    # is fixed; judges are the sampling units), relabel the copies so
    # clusters stay distinct, and rerun the whole pipeline on each
    # resample: unlike the parametric bootstrap, this carries any
    # extra-model dependence within a judge's comparisons
    jd_rows <- split(seq_along(a), jd)
    pan_of_judge <- vapply(jd_rows, function(r) pan[r[1]], "")
    judges_by_panel <- split(names(jd_rows), pan_of_judge)
    draws <- list()
    for (bb in seq_len(boot_reps)) {
      take <- unlist(lapply(judges_by_panel, function(js)
        sample(js, length(js), replace = TRUE)), use.names = FALSE)
      idx <- unlist(jd_rows[take], use.names = FALSE)
      jd_new <- rep(paste0(take, "#", seq_along(take)),
                    lengths(jd_rows[take]))
      df_b <- data.frame(oa = a[idx], ob = b[idx],
                         win = ifelse(y[idx] == 1L, a[idx], b[idx]),
                         judge = jd_new, stringsAsFactors = FALSE)
      pmap_b <- setNames(pan[idx][!duplicated(jd_new)], unique(jd_new))
      fb <- tryCatch(suppressWarnings(
        btl_efrm(df_b, "oa", "ob", "win", "judge", panels = pmap_b,
                 object_sets = object_sets, ties = "drop",
                 min_link = min_link, se_method = "conditional",
                 maxit = maxit, tol = tol)), error = function(e) NULL)
      if (is.null(fb) || !isTRUE(fb$converged)) { boot_fail <- boot_fail + 1L; next }
      lphi_b <- log(fb$phi_table$phi)[match(panels_u, fb$phi_table$panel)]
      la_b <- log(fb$alpha_table$alpha)[match(sets_u, fb$alpha_table$set)]
      ka_b <- fb$kappa_table$kappa[match(sets_u, fb$kappa_table$set)]
      bh_b <- fb$objects$beta_set[match(objs_all, fb$objects$object)]
      v_b <- fb$objects$v[match(objs_all, fb$objects$object)]
      if (anyNA(c(lphi_b, la_b, ka_b))) { boot_fail <- boot_fail + 1L; next }
      draws[[length(draws) + 1L]] <- c(lphi_b, la_b, ka_b, bh_b, v_b)
    }
    if (length(draws) < max(20L, ceiling(boot_reps / 2)))
      stop("judge bootstrap failed on ", boot_fail, " of ", boot_reps,
           " replicates; too few judges per panel for stable resampling -- ",
           "use se_method = 'bootstrap' or 'conditional'")
    D <- do.call(rbind, draws)
    colnames(D) <- c(paste0("log phi[", panels_u, "]"),
                     paste0("log alpha[", sets_u, "]"),
                     paste0("kappa[", sets_u, "]"),
                     paste0("beta[", objs_all, "]"), paste0("v[", objs_all, "]"))
    n_inf <- colSums(!is.finite(D))
    sds <- rep(NA_real_, ncol(D))
    ok_col <- n_inf == 0L
    sds[ok_col] <- apply(D[, ok_col, drop = FALSE], 2, sd)
    if (any(n_inf > 0L))
      notes <- c(notes, paste0(
        "bootstrap: ", paste(sprintf("%s reached the boundary in %d of %d replicates",
                                     colnames(D)[n_inf > 0L], n_inf[n_inf > 0L],
                                     length(draws)), collapse = "; "),
        "; the parameter is weakly identified and its SE is reported as NA"))
    nO <- length(objs_all)
    se_log_phi <- setNames(sds[seq_len(G)], panels_u)
    se_log_alpha <- setNames(sds[G + seq_len(S)], sets_u)
    se_kappa <- setNames(sds[G + S + seq_len(S)], sets_u)
    se_log_alpha[sets_u[1]] <- NA_real_
    se_kappa[sets_u[1]] <- NA_real_
    se_bhat <- setNames(sds[G + 2L * S + seq_len(nO)], objs_all)
    se_v <- setNames(sds[G + 2L * S + nO + seq_len(nO)], objs_all)
    if (boot_fail > 0)
      notes <- c(notes, sprintf(
        "judge bootstrap: %d of %d replicates failed and were skipped",
        boot_fail, boot_reps))
  }
  if (se_method == "bootstrap") {
    p_hat <- pmin(pmax(fit0$p_all, 1e-8), 1 - 1e-8)
    draws <- list()
    for (bb in seq_len(boot_reps)) {
      yb <- rbinom(length(p_hat), 1L, p_hat)
      fb <- tryCatch(fit_once(yb), error = function(e) NULL)
      if (is.null(fb) || !isTRUE(fb$converged)) { boot_fail <- boot_fail + 1L; next }
      draws[[length(draws) + 1L]] <- c(log(fb$phi), log(fb$alpha), fb$kappa,
                                       fb$bhat, fb$v)
    }
    if (length(draws) < max(20L, ceiling(boot_reps / 2)))
      stop("parametric bootstrap failed on ", boot_fail, " of ", boot_reps,
           " replicates; the design is too sparse for stable resampling -- ",
           "add comparisons or use se_method = 'conditional'")
    D <- do.call(rbind, draws)
    colnames(D) <- c(paste0("log phi[", panels_u, "]"),
                     paste0("log alpha[", sets_u, "]"),
                     paste0("kappa[", sets_u, "]"),
                     paste0("beta[", objs_all, "]"), paste0("v[", objs_all, "]"))
    # a parameter that reaches its boundary in some replicates (a set unit
    # driven to zero when a resampled within-set order flips against the
    # cross-set evidence) has no normal sampling distribution: report NA
    # rather than a standard deviation over infinite draws
    n_inf <- colSums(!is.finite(D))
    sds <- rep(NA_real_, ncol(D))
    ok_col <- n_inf == 0L
    sds[ok_col] <- apply(D[, ok_col, drop = FALSE], 2, sd)
    if (any(n_inf > 0L))
      notes <- c(notes, paste0(
        "bootstrap: ", paste(sprintf("%s reached the boundary in %d of %d replicates",
                                     colnames(D)[n_inf > 0L], n_inf[n_inf > 0L],
                                     length(draws)), collapse = "; "),
        "; the parameter is weakly identified and its SE is reported as NA"))
    nO <- length(objs_all)
    se_log_phi <- setNames(sds[seq_len(G)], panels_u)
    se_log_alpha <- setNames(sds[G + seq_len(S)], sets_u)
    se_kappa <- setNames(sds[G + S + seq_len(S)], sets_u)
    se_log_alpha[sets_u[1]] <- NA_real_                 # reference: fixed at 1 / 0
    se_kappa[sets_u[1]] <- NA_real_
    se_bhat <- setNames(sds[G + 2L * S + seq_len(nO)], objs_all)
    se_v <- setNames(sds[G + 2L * S + nO + seq_len(nO)], objs_all)
    if (boot_fail > 0)
      notes <- c(notes, sprintf(
        "bootstrap: %d of %d replicates failed and were skipped", boot_fail,
        boot_reps))
  }

  # --- equal-unit (single-unit) comparison ----------------------------------
  ll_frames <- ll_within + ll_cross
  single <- tryCatch(
    .btlef_stage1(match(a, objs_all), match(b, objs_all), y,
                  rep("all", length(a)), jd, length(objs_all), maxit, tol),
    error = function(e) NULL)
  ll_single <- if (is.null(single)) NA_real_ else single$ll
  equal_unit <- list(
    loglik_frames = ll_frames, loglik_single = ll_single,
    difference = if (is.na(ll_single)) NA_real_ else ll_frames - ll_single,
    note = paste("descriptive composite-likelihood difference;",
                 "the Wald tests on log phi and log alpha carry the inference"))

  # --- structural tables ----------------------------------------------------
  z_phi <- log(phi) / se_log_phi
  phi_table <- data.frame(panel = panels_u, phi = unname(phi),
                          se_log_phi = unname(se_log_phi),
                          z = unname(z_phi), p = unname(2 * pnorm(-abs(z_phi))),
                          stringsAsFactors = FALSE)
  z_al <- log(alpha) / se_log_alpha
  alpha_table <- data.frame(set = sets_u, alpha = unname(alpha),
                            se_log_alpha = unname(se_log_alpha),
                            z = unname(z_al), p = unname(2 * pnorm(-abs(z_al))),
                            stringsAsFactors = FALSE)
  z_ka <- kappa / se_kappa
  kappa_table <- data.frame(set = sets_u, kappa = unname(kappa),
                            se_kappa = unname(se_kappa),
                            z = unname(z_ka), p = unname(2 * pnorm(-abs(z_ka))),
                            stringsAsFactors = FALSE)

  objects <- data.frame(object = objs_all, set = unname(set_of[objs_all]),
                        beta_set = unname(bhat[objs_all]),
                        se_beta = unname(se_bhat[objs_all]),
                        v = unname(v), se_v = unname(se_v),
                        stringsAsFactors = FALSE)
  rownames(objects) <- NULL

  # frames: one row per panel-by-set cell holding within-set comparisons
  fr <- list()
  for (s in sets_u) for (g in panels_u) {
    rows <- which(within & sa == s & pan == g)
    if (!length(rows)) next
    fr[[length(fr) + 1L]] <- data.frame(
      panel = g, set = s, rho = unname(phi[[g]] * alpha[[s]]),
      n_comparisons = length(rows),
      fit_resid = .btlef_frame_fit(y[rows], within_p[rows]),
      stringsAsFactors = FALSE)
  }
  frames <- if (length(fr)) do.call(rbind, fr) else NULL
  if (!is.null(frames)) rownames(frames) <- NULL

  if (S == 1L)
    notes <- c(notes, "single set: panel-units model (set units alpha not estimated)")
  if (G == 1L && S == 1L)
    notes <- c(notes, "single panel and single set: reduces to btl()")
  if (length(fit0$dropped))
    notes <- c(notes, paste0(
      "set(s) ", paste(fit0$dropped, collapse = ", "), " carry no stable ",
      "panel-ratio information (within-set contests too close to even or ",
      "too one-sided); they were excluded from the phi reconciliation and ",
      "refit with the panel units held at the reconciled phi"))

  out <- list(objects = objects, phi_table = phi_table,
              alpha_table = alpha_table, kappa_table = kappa_table,
              frames = frames, equal_unit = equal_unit, n_cross = n_cross,
              sets = sets_u, panels = panels_u, reference_set = sets_u[1],
              n_comparisons = length(a),
              converged = fit0$converged,
              se_method = se_method,
              boot_reps = if (se_method %in% c("bootstrap", "judge_bootstrap"))
                boot_reps else NA_integer_,
              se_note = if (se_method == "judge_bootstrap")
                paste("standard errors from a judge-resampling bootstrap",
                      "(judges redrawn with replacement within panels, the",
                      "whole pipeline refitted): carries stage-one",
                      "uncertainty AND any extra-model dependence within",
                      "judges")
              else if (se_method == "bootstrap")
                paste("standard errors from a parametric bootstrap of the",
                      "whole two-stage pipeline: the staged conditional",
                      "estimates are unchanged, and their errors carry the",
                      "stage-one uncertainty into the linking")
              else
                paste("CONDITIONAL standard errors: exact for beta and phi,",
                      "but the alpha and kappa errors are conditional on",
                      "stage 1 and UNDERSTATE the total sampling variability;",
                      "use se_method = 'bootstrap' for inference"),
              notes = notes)
  class(out) <- "rasch_btl_efrm"
  out
}

#' @export
print.rasch_btl_efrm <- function(x, ...) {
  cat(sprintf(paste0("Bradley-Terry-Luce extended frame of reference: ",
                     "%d objects in %d set(s) x %d panel(s), %d comparisons\n"),
              nrow(x$objects), nrow(x$alpha_table), nrow(x$phi_table),
              x$n_comparisons))
  cat(sprintf("Two-stage conditional ML: %s; SEs %s\n",
              if (x$converged) "converged" else "NOT converged",
              if (identical(x$se_method, "judge_bootstrap"))
                sprintf("by judge-resampling bootstrap (B = %d)", x$boot_reps)
              else if (identical(x$se_method, "bootstrap"))
                sprintf("by parametric bootstrap (B = %d)", x$boot_reps)
              else "conditional (understate the linking uncertainty)"))
  if (nrow(x$alpha_table) == 1L)
    cat("Model: panel units only (single set; set units not estimated)\n")
  cat("\nPanel units (phi; Wald H0: log phi = 0):\n")
  print(.fmt_df(x$phi_table), row.names = FALSE)
  if (nrow(x$alpha_table) > 1L) {
    cat("\nSet units (alpha) and origins (kappa; reference set = ",
        x$reference_set, "):\n", sep = "")
    at <- merge(x$alpha_table[, c("set", "alpha", "se_log_alpha", "p")],
                x$kappa_table[, c("set", "kappa", "se_kappa")],
                by = "set", sort = FALSE)
    print(.fmt_df(at), row.names = FALSE)
  }
  eu <- x$equal_unit
  if (!is.na(eu$difference))
    cat(sprintf(paste0("\nEqual-unit comparison: ll_frames - ll_single = ",
                       "%.3f (%s)\n"), eu$difference, eu$note))
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' Plot the frame units of a paired-comparison EFRM fit
#'
#' Caterpillar plot of the estimated units on the log scale: one row per panel
#' unit \code{phi_g} and one per set unit \code{alpha_s}, with 95 per cent
#' intervals, the reference (unit one) marked, mirroring
#' \code{\link{plot_frames}} in the package's house style.
#'
#' @param fit A fitted object from \code{\link{btl_efrm}}.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?btl_efrm for a complete simulated example
#' }
#' @export
plot_btl_units <- function(fit) {
  if (!inherits(fit, "rasch_btl_efrm"))
    stop("plot_btl_units needs a rasch_btl_efrm fit")
  ph <- fit$phi_table; al <- fit$alpha_table
  rows <- rbind(
    data.frame(label = paste0("panel: ", ph$panel), kind = "panel",
               est = log(ph$phi), se = ph$se_log_phi, stringsAsFactors = FALSE),
    if (nrow(al) > 1L)
      data.frame(label = paste0("set: ", al$set), kind = "set",
                 est = log(al$alpha), se = al$se_log_alpha,
                 stringsAsFactors = FALSE))
  rows$se[!is.finite(rows$se)] <- 0
  rows <- rows[order(rows$kind, rows$est), ]
  n <- nrow(rows)
  lo <- rows$est - 1.96 * rows$se; hi <- rows$est + 1.96 * rows$se
  colr <- ifelse(rows$kind == "panel", .rr$blue, .rr$purple)
  op <- par(mar = c(4.2, 9, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = range(c(lo, hi, 0)) + c(-0.1, 0.1), ylim = c(0.5, n + 0.5),
       xlab = "log unit", ylab = "", axes = FALSE, main = "")
  abline(h = seq_len(n), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  axis(1, col = .rr$grid, col.ticks = .rr$soft)
  axis(2, at = seq_len(n), labels = rows$label, cex.axis = 0.75,
       col = .rr$grid, col.ticks = NA)
  segments(lo, seq_len(n), hi, seq_len(n), lwd = 2.2, col = .rr$soft)
  points(rows$est, seq_len(n), pch = 21, cex = 1.5, bg = colr,
         col = "white", lwd = 1.2)
  .rr_legend("bottomright", c("panel unit (phi)", "set unit (alpha)"),
             pch = 21, pt.bg = c(.rr$blue, .rr$purple), col = "white",
             pt.cex = 1.2)
  invisible(NULL)
}
