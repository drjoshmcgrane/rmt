# rasch :: top-level analysis
# ===========================================================================
# rasch() ties the engine together: data preparation (with category
# collapsing and constant-item removal, recorded as notes), pairwise
# conditional ML item estimation, Warm WLE person estimation per
# missing-data pattern, the full test-of-fit suite, and the score-to-measure
# table. An ID variable and any number of person factors carry through to
# the person estimates.
# ===========================================================================

# Map missing-data codes (and any negative score) to NA. Valid item scores
# are non-negative integers from zero; by long-standing convention -1 marks a
# missing response, so any value below zero is read as missing.
.apply_na_codes <- function(v, na_codes) {
  v[v %in% na_codes | (!is.na(v) & v < 0)] <- NA
  v
}

# Prepare the item matrix: integer scores from 0, consecutive observed
# categories, no constant items. Returns the matrix plus human-readable notes.
.prepare_X <- function(X, na_codes = -1) {
  notes <- character(0)
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- sprintf("I%02d", seq_len(ncol(X)))
  Xi <- suppressWarnings(apply(X, 2, function(col) as.integer(as.character(col))))
  dim(Xi) <- dim(X); dimnames(Xi) <- dimnames(X)
  bad_num <- colSums(!is.na(X) & is.na(Xi)) > 0
  if (any(bad_num))
    notes <- c(notes, paste0("non-numeric entries set to missing in: ",
                             paste(colnames(X)[bad_num], collapse = ", ")))
  n_na <- sum(!is.na(Xi) & (Xi %in% na_codes | Xi < 0))
  Xi[] <- .apply_na_codes(Xi, na_codes)
  if (n_na > 0) {
    codes <- paste(unique(c(na_codes, "negative")), collapse = ", ")
    notes <- c(notes, sprintf("%d cell(s) with a missing-data code (%s) set to missing",
                              n_na, codes))
  }
  X <- Xi
  const <- apply(X, 2, function(col) length(unique(col[!is.na(col)])) < 2)
  if (any(const)) {
    notes <- c(notes, paste0("dropped constant item(s): ",
                             paste(colnames(X)[const], collapse = ", ")))
    X <- X[, !const, drop = FALSE]
  }
  if (ncol(X) < 2) stop("need at least two non-constant items")
  for (i in seq_len(ncol(X))) {
    v <- X[, i]; obs <- sort(unique(v[!is.na(v)]))
    full <- seq(0L, max(obs))
    if (!identical(obs, full)) {
      X[, i] <- match(v, obs) - 1L
      notes <- c(notes, sprintf("item %s rescored: observed categories [%s] mapped to 0:%d",
                                colnames(X)[i], paste(obs, collapse = ","),
                                length(obs) - 1L))
    }
  }
  list(X = X, notes = notes)
}

#' Fit and diagnose a Rasch model by pairwise conditional estimation
#'
#' Runs a complete Rasch analysis: Andrich and Luo pairwise conditional
#' maximum likelihood item estimation (see \code{\link{pcml}}), Warm weighted
#' likelihood person estimates per missing-data pattern, item and person fit
#' residuals (the log-of-mean-square statistic of Andrich and Marais 2019,
#' ch. 23, with its
#' untransformed natural form and degrees of freedom), infit and outfit,
#' the item-trait interaction chi-square and the class-interval ANOVA
#' item-fit F, the person separation index with and without extremes,
#' Cronbach's alpha, targeting, threshold diagnostics, and the
#' score-to-measure table.
#'
#' The fit residual follows Andrich and Marais (2019, ch. 23) exactly:
#' standardised residuals are
#' squared and summed over each item's (person's) observed cells among
#' non-extreme persons, compared with the summed cell degrees of freedom
#' (the model-testing degrees of freedom, cells minus estimated parameters,
#' apportioned equally over cells), and symmetrised by the
#' log-of-mean-square transform \eqn{f (\ln Y^2 - \ln f)/\sqrt{V[Y^2]}} with
#' model-based variance \eqn{V[Y^2] = \sum (C_4/V^2 - 1)}. Values are
#' approximately N(0,1) under fit; the conventional flagging value is 2.5
#' (Andrich and Marais 2019, ch. 15).
#' Negative values indicate over-discrimination (Guttman-like responses),
#' positive values under-discrimination.
#'
#' @param data Persons-by-items integer score matrix (categories from 0), or a
#'   data frame also containing ID and person-factor columns. Missing values
#'   are allowed.
#' @param model Either \code{"PCM"} (partial credit) or \code{"RSM"} (rating
#'   scale).
#' @param id Optional name of an ID column in \code{data}, or a vector of IDs;
#'   carried through to the person estimates.
#' @param factors Optional character vector of person-factor column names in
#'   \code{data} (for DIF analysis), or a data frame of factors.
#' @param items Optional character vector naming the item columns; by default
#'   every column not named in \code{id} or \code{factors}.
#' @param n_groups Number of class intervals for the item-trait chi-square
#'   and ANOVA item fit. The default \code{NULL} applies the rule of Andrich
#'   and Marais (2019, ch. 15): as
#'   many intervals of at least 50 non-extreme persons as the sample allows,
#'   at most 10, at least 2. The resolved value is stored in
#'   \code{fit$n_groups}.
#' @param adjust_N Optional reference sample size; if supplied, item-trait
#'   chi-squares are rescaled to this size (a sample-size adjustment for the
#'   sensitivity of the chi-square to large samples).
#' @param anchors Optional anchor table for equating: a data frame with
#'   columns \code{item}, \code{k}, and \code{tau} fixing nominated
#'   thresholds at known values; see \code{\link{pcml}}. With anchors in
#'   place the scale origin comes from the anchors, so person measures are
#'   directly comparable across separately analysed datasets.
#' @param na_codes Values to read as missing. Defaults to \code{-1}, the
#'   conventional missing-response code; any negative score is also treated as
#'   missing, since valid category scores start at zero.
#' @param maxit,tol Newton-Raphson iteration cap and convergence
#'   tolerance of the pairwise conditional estimation.
#' @param key Optional multiple-choice scoring key, in any of three forms. (1) A named vector or data frame with columns \code{item}
#'   and \code{key} naming each item's correct option: scored 0/1
#'   (case-insensitive after trimming; blanks become missing). (2) Double
#'   keying: several correct options separated by \code{"/"} (for example
#'   \code{"A/C"}), all scoring 1. (3) Polytomous option scoring (Andrich
#'   and Styles 2011): a data frame with columns \code{item},
#'   \code{option}, and \code{score} assigning an integer score to every
#'   credited option (unlisted options score 0), so informative
#'   distractors receive partial credit and the item is fitted as
#'   polytomous; see \code{\link{distractor_rescore}} for an
#'   evidence-based proposal. Raw responses are retained in \code{fit$mc}
#'   for \code{\link{distractor_analysis}} and
#'   \code{\link{plot_distractors}}.
#' @param pc_components \code{NULL} (default) estimates every PCM threshold
#'   freely. An integer 1 to 4 instead estimates each item's thresholds
#'   through the Andrich principal-components reparameterisation (see
#'   \code{\link{pcml_pc}}): 1 = location only, 2 = + spread (the dispersion
#'   model of Andrich 1982), 3 = + skewness, 4 = + kurtosis (the full
#'   principal-components model; Pedler 1987). Useful when some categories are sparsely
#'   populated; the component estimates are returned in
#'   \code{fit$est$components}. PCM only, and not combinable with anchors.
#' @return An object of class \code{"rasch"}: a list with the item summary
#'   (\code{items}), \code{thresholds} (with standard errors), the person
#'   table (\code{person}, including ID and factors), the score table,
#'   residuals, reliability (\code{psi}, \code{psi_noext}, the item
#'   separation index \code{isi}, \code{alpha}), targeting, item-trait
#'   statistics (\code{item_trait}, \code{item_anova}), the
#'   summary distribution block (\code{summary_stats}: location and fit
#'   residual mean/SD/skewness/kurtosis, fit-location correlations, and the
#'   cell degrees-of-freedom factor), threshold diagnostics, and estimation
#'   details (\code{est}).
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(X, model = "PCM")
#' fit$items
#' fit$psi$PSI
#' @export
rasch <- function(data, model = c("PCM", "RSM"), id = NULL, factors = NULL,
                  items = NULL, n_groups = NULL, adjust_N = NA, anchors = NULL,
                  na_codes = -1, key = NULL, pc_components = NULL,
                  maxit = 60, tol = 1e-8) {
  model <- match.arg(model)
  if (!is.null(pc_components)) {
    if (model != "PCM")
      stop("pc_components applies to the PCM only")
    if (!is.null(anchors))
      stop("pc_components cannot be combined with anchors")
  }

  # --- split data frame into ID, factors, and item columns ---------------
  id_vec <- NULL; fac_df <- NULL
  if (is.data.frame(data)) {
    nm <- names(data)
    if (is.character(id) && length(id) == 1L && id %in% nm) {
      id_vec <- data[[id]]
    } else if (!is.null(id) && length(id) == nrow(data)) id_vec <- id
    if (is.character(factors) && all(factors %in% nm)) {
      fac_df <- data[, factors, drop = FALSE]
    } else if (is.data.frame(factors) && nrow(factors) == nrow(data)) fac_df <- factors
    drop_cols <- c(if (is.character(id)) id else NULL,
                   if (is.character(factors)) factors else NULL)
    item_cols <- if (!is.null(items)) intersect(items, nm) else setdiff(nm, drop_cols)
    X <- as.matrix(data[, item_cols, drop = FALSE])
  } else {
    X <- as.matrix(data)
    if (!is.null(id) && length(id) == nrow(X)) id_vec <- id
    if (is.data.frame(factors) && nrow(factors) == nrow(X)) fac_df <- factors
  }
  if (is.null(id_vec)) id_vec <- seq_len(nrow(X))

  # score multiple-choice items against the key, keeping the raw responses
  mc <- NULL
  if (!is.null(key)) {
    key <- .resolve_key(key)
    sc <- .score_mc(X, key)
    X[, colnames(sc$scored)] <- sc$scored
    mc <- list(key = sc$key, map = sc$map, raw = sc$raw)
  }

  prep <- .prepare_X(X, na_codes = na_codes); X <- prep$X
  if (!is.null(mc)) {
    prep$notes <- c(prep$notes,
                    sprintf("%d item(s) scored 0/1 against the key", ncol(mc$raw)))
    gone <- setdiff(colnames(mc$raw), colnames(X))
    if (length(gone)) mc$raw <- mc$raw[, setdiff(colnames(mc$raw), gone), drop = FALSE]
  }

  if (!is.null(anchors)) {
    a_names <- if (is.character(anchors$item) || is.factor(anchors$item))
      as.character(anchors$item) else colnames(X)[anchors$item]
    gone <- setdiff(a_names, colnames(X))
    if (length(gone))
      stop("anchored item(s) not present after data preparation: ",
           paste(gone, collapse = ", "))
    resc <- grepl("rescored", prep$notes) &
      vapply(prep$notes, function(n) any(vapply(a_names, grepl, TRUE, x = n)), TRUE)
    if (any(resc))
      stop("anchored item(s) were rescored during data preparation; ",
           "anchor values would no longer match the threshold numbering")
    prep$notes <- c(prep$notes,
                    sprintf("%d threshold(s) anchored; scale origin from anchors",
                            nrow(anchors)))
  }

  # --- item estimation ----------------------------------------------------
  est <- if (is.null(pc_components))
    pcml(X, model = model, anchors = anchors, maxit = maxit, tol = tol)
  else pcml_pc(X, n_components = pc_components, maxit = maxit, tol = tol)
  if (!is.null(pc_components))
    prep$notes <- c(prep$notes,
                    sprintf("thresholds estimated through %d principal component(s); see est$components",
                            pc_components))
  fit <- .assemble_fit(model, X, est, id_vec, fac_df, n_groups, adjust_N,
                       c(prep$notes, est$notes))
  fit$mc <- mc
  fit
}

# Post-estimation pipeline shared by rasch(), rasch_mfrm(), and rasch_efrm():
# person estimation, residuals, the full fit suite, and the assembled tables.
# disc is an optional per-column discrimination (frame unit) vector; with
# unequal discriminations the raw score is no longer sufficient, so person
# estimation switches to the weighted-score routine and the score table is
# replaced by per-unit score curves.
.assemble_fit <- function(model, X, est, id_vec, fac_df, n_groups, adjust_N,
                          notes, disc = NULL) {
  m <- est$m; L <- ncol(X)
  thr <- est$thr
  tau_list <- lapply(seq_len(L), function(i) thr$tau[thr$item == i])
  names(tau_list) <- colnames(X)
  equal_disc <- is.null(disc) || length(unique(disc)) == 1L
  disc_v <- if (is.null(disc)) rep(1, L) else disc

  # --- person estimation and residuals -------------------------------------
  person <- if (equal_disc) .person_estimates(X, tau_list, disc = disc_v[1])
            else .efrm_person_estimates(X, tau_list, disc_v)
  mo <- .moment_arrays(person$theta, tau_list, disc = disc_v)
  Z <- (X - mo$E) / sqrt(mo$V)
  colnames(Z) <- colnames(X)

  # --- fit statistics ------------------------------------------------------
  ifit <- .item_fit(X, Z, mo, disc = if (is.null(disc)) NULL else disc_v)
  pfit <- .person_fit(X, Z, mo, disc = if (is.null(disc)) NULL else disc_v)
  n_par <- if (is.null(est$n_parameters)) nrow(est$thr) - 1L else est$n_parameters
  rf <- .fitres(Z, mo, person$extreme, n_par)
  ng_req <- n_groups                       # NULL = the automatic rule
  ci <- .class_intervals(person$theta, person$extreme, n_groups)
  n_groups <- attr(ci, "n_groups")
  # class intervals compiled per item when data are missing (the automatic
  # per-item basis; Andrich & Marais 2019, ch. 15), so every item is tested
  # over intervals of
  # its own responders (with the group-count rule applied per item when
  # automatic)
  ci_list <- if (anyNA(X))
    .class_intervals_by_item(X, person$theta, person$extreme, ng_req)
  else NULL
  it <- .item_trait(X, mo, ci, adjust_N = adjust_N, ci_list = ci_list)
  ia <- .item_anova(Z, ci, person$extreme, ci_list = ci_list)
  psi <- .psi(person$theta, person$se)
  psi_noext <- .psi(person$theta, person$se, keep = !person$extreme)
  alpha <- .alpha(X)

  # --- assembled person table ----------------------------------------------
  parts <- list(data.frame(id = id_vec), fac_df, person,
                data.frame(infit_ms = pfit$infit_ms, outfit_ms = pfit$outfit_ms,
                           outfit_z = pfit$outfit_z,
                           fit_resid = rf$persons$fit_resid,
                           natural_resid = rf$persons$natural,
                           df_fit = rf$persons$df, class_interval = ci))
  person <- do.call(cbind, parts[!vapply(parts, is.null, TRUE)])
  rownames(person) <- NULL

  # --- item table ----------------------------------------------------------
  loc <- vapply(tau_list, mean, 0)
  weak_thr <- if (is.null(thr$weak)) rep(FALSE, nrow(thr)) else thr$weak
  se_loc <- vapply(seq_len(L), function(i) {
    rows <- thr$id[thr$item == i]
    # a weakly determined threshold (sparse adjacent category) makes the
    # ridged covariance block spuriously small: report NA, not a number
    if (any(weak_thr[thr$item == i])) return(NA_real_)
    # anchored items have a structurally zero variance that floating-point
    # noise can render as a tiny negative number on some BLAS builds
    sqrt(max(mean(est$cov_tau[rows, rows]), 0))
  }, 0)
  items_df <- data.frame(item = colnames(X), max = m, location = loc,
                         se = se_loc,
                         fit_resid = rf$items$fit_resid, df_fit = rf$items$df,
                         natural_resid = rf$items$natural,
                         infit_ms = ifit$infit_ms, outfit_ms = ifit$outfit_ms,
                         infit_z = ifit$infit_z, outfit_z = ifit$outfit_z,
                         chisq = it$chisq, df = it$df, p = it$p,
                         p_adj = it$p_adj, p_bonf = it$p_bonf,
                         F_anova = ia$F_anova, p_anova = ia$p)
  rownames(items_df) <- NULL

  # --- score table (complete responders; raw score is only sufficient when
  # --- discriminations are equal) ---------------------------------------------
  sc <- if (equal_disc) {
    pe_full <- person_wle(tau_list, disc = disc_v[1])
    data.frame(score = 0:sum(m), theta = unname(pe_full$theta),
               se = unname(pe_full$se))
  } else NULL
  if (is.null(sc))
    notes <- c(notes, "person measures use the weighted score; per-group score curves replace the raw-score table (see score_curves)")

  # --- threshold diagnostics --------------------------------------------------
  td <- lapply(seq_len(L), function(i) {
    tau_i <- tau_list[[i]]
    grid <- seq(-8, 8, by = 0.05)
    modal <- unique(vapply(grid, function(th)
      which.max(item_moments(th, tau_i, disc = disc_v[i])$P) - 1L, 1L))
    list(item = colnames(X)[i], thresholds = tau_i,
         ordered = all(diff(tau_i) > 0) || length(tau_i) == 1L,
         reversed_at = which(diff(tau_i) <= 0) + 1L,
         never_modal_categories = setdiff(0:length(tau_i), modal),
         category_counts = as.integer(table(factor(X[, i], levels = 0:length(tau_i)))))
  })
  names(td) <- colnames(X)

  out <- list(model = model, X = X, m = m, items = items_df, thresholds = thr,
              tau_list = tau_list, person = person, score_table = sc,
              residuals = Z, moments = mo, n_groups = n_groups,
              ci_item = ci_list,
              item_trait = it, item_anova = ia,
              psi = psi, psi_noext = psi_noext,
              isi = .psi(items_df$location, items_df$se),
              alpha = alpha,
              targeting = .targeting(person, thr),
              power_of_fit = .fit_power(psi$PSI),
              total_chisq = sum(it$chisq), total_df = sum(it$df),
              total_chisq_p = pchisq(sum(it$chisq), sum(it$df), lower.tail = FALSE),
              item_fit_summary = .dist_stats(rf$items$fit_resid),
              person_fit_summary = .dist_stats(rf$persons$fit_resid),
              summary_stats = list(
                item_location = .dist_stats(items_df$location),
                person_location = .dist_stats(person$theta),
                person_location_noext = .dist_stats(person$theta[!person$extreme]),
                cor_item_fit_location = .safe_cor(items_df$location,
                                                  rf$items$fit_resid),
                cor_person_fit_location = .safe_cor(person$theta,
                                                    rf$persons$fit_resid),
                df_factor = rf$f_cell),
              thresholds_diag = td, est = est, notes = notes,
              factors = fac_df, disc = disc)
  class(out) <- "rasch"
  out
}

#' @export
print.rasch <- function(x, ...) {
  cat(sprintf("rasch %s analysis: %d items, %d persons\n",
              x$model, ncol(x$X), nrow(x$X)))
  cat(sprintf("Pairwise conditional ML (Andrich & Luo): %s in %d iterations\n",
              if (x$est$converged) "converged" else "NOT converged",
              x$est$iterations))
  cat(sprintf("PSI %.3f (no extremes %.3f), item SI %.3f, alpha %.3f%s, power of fit: %s\n",
              x$psi$PSI, x$psi_noext$PSI, x$isi$PSI, x$alpha$alpha,
              if (isFALSE(x$alpha$applicable))
                sprintf(" [complete cases only, n = %d]", x$alpha$n) else "",
              x$power_of_fit))
  cat(sprintf("Total item-trait chi-square %.3f on %d df, p = %s\n",
              x$total_chisq, x$total_df, .fmt_p(x$total_chisq_p)))
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' @export
summary.rasch <- function(object, ...) {
  x <- object
  print(x)
  cat(sprintf("\nTargeting: person mean %.3f (SD %.3f); thresholds span %.3f to %.3f\n",
              x$targeting$person_mean, x$targeting$person_sd,
              x$targeting$threshold_range[1], x$targeting$threshold_range[2]))
  cat(sprintf("Item fit residual mean %.3f SD %.3f (skew %.2f, kurt %.2f); person fit residual mean %.3f SD %.3f (skew %.2f, kurt %.2f)\n",
              x$item_fit_summary$mean, x$item_fit_summary$sd,
              x$item_fit_summary$skewness, x$item_fit_summary$kurtosis,
              x$person_fit_summary$mean, x$person_fit_summary$sd,
              x$person_fit_summary$skewness, x$person_fit_summary$kurtosis))
  cat(sprintf("Fit residual-location correlation: items %.3f, persons %.3f; cell df factor %.3f\n",
              x$summary_stats$cor_item_fit_location,
              x$summary_stats$cor_person_fit_location,
              x$summary_stats$df_factor))
  cat(sprintf("Items with adjusted chi-square p < 0.05: %d of %d\n\n",
              sum(x$items$p_adj < 0.05, na.rm = TRUE), nrow(x$items)))
  core <- c("item", "max", "location", "se", "fit_resid", "infit_ms",
            "outfit_ms", "chisq", "df", "p_adj")
  print(.fmt_df(x$items[, intersect(core, names(x$items))]), row.names = FALSE)
  cat("(further columns on fit$items: natural and standardised forms,\n",
      " ANOVA fit, Bonferroni probabilities)\n", sep = "")
  dis <- vapply(x$thresholds_diag, function(d) !d$ordered, TRUE) &
    vapply(x$thresholds_diag, function(d) length(d$thresholds) > 1L, TRUE)
  if (any(dis)) cat(sprintf("\nDisordered thresholds: %s\n", paste(names(dis)[dis], collapse = ", ")))
  invisible(x)
}
