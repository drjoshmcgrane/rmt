# rmt :: differential item functioning
# ===========================================================================
# DIF by analysis of variance of the standardised residuals (Hagquist &
# Andrich 2017). For each item, residuals are analysed by person factor and
# trait class interval: a factor main effect indicates uniform DIF and a
# factor-by-interval interaction indicates non-uniform DIF. With several
# person factors the analysis can be run factor-at-a-time (dif_anova) or as
# one factorial model per item (dif_anova_factorial) with factor-by-factor
# interactions, Tukey HSD comparisons on the significant group terms, and
# the convention that a significant interaction supersedes the main effects
# of the factors involved. Multiplicity across items is handled by
# Benjamini-Hochberg false-discovery-rate adjustment.
# ===========================================================================

.dif_factors <- function(fit, factors) {
  if (is.null(factors)) factors <- fit$factors
  if (is.null(factors)) stop("no person factors supplied or stored in the fit")
  if (is.character(factors) && !is.null(fit$factors) &&
      all(factors %in% names(fit$factors)))
    factors <- fit$factors[, factors, drop = FALSE]
  if (!is.data.frame(factors)) factors <- data.frame(group = factors)
  factors
}

# Class intervals for a DIF analysis are set from the cells the analysis
# actually uses: the residual ANOVA crosses trait intervals with group
# levels (or with the factor-combination cells in the factorial), so the
# interval count is chosen to keep the smallest group's expected cell size
# adequate -- independently of the interval count of the overall fit.
.dif_n_groups <- function(fit, grp, cell_min = 30L) {
  ok <- !is.na(grp) & !is.na(fit$person$theta)
  if (!any(ok)) return(2L)
  n_min <- min(table(droplevels(factor(grp[ok]))))
  max(2L, min(10L, as.integer(n_min) %/% as.integer(cell_min)))
}

.dif_class_intervals <- function(fit, n_groups) {
  ci <- fit$person$class_interval
  if (is.null(ci) || !identical(n_groups, fit$n_groups))
    ci <- .class_intervals(fit$person$theta, fit$person$extreme, n_groups)
  factor(ci)
}

#' Differential item functioning by two-way residual ANOVA
#'
#' For each item and each person factor separately, analyses the
#' standardised residuals by factor group and trait class interval. The
#' group main effect indicates uniform DIF and the group-by-interval
#' interaction indicates non-uniform DIF. Probabilities are adjusted across
#' items within each factor by the Benjamini-Hochberg false-discovery-rate
#' procedure (or any \code{\link[stats]{p.adjust}} method). With several
#' factors, consider \code{\link{dif_anova_factorial}}, which models them
#' jointly.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors A vector (one factor), a data frame of person factors, or a
#'   character vector naming factor columns already nominated in the fit (via
#'   \code{rasch(..., factors = )}). Defaults to every factor stored in the
#'   fit.
#' @param n_groups Number of trait class intervals. By default set per
#'   factor from the smallest group so every interval-by-group cell keeps
#'   about 30 expected responses (between 2 and 10 intervals) --
#'   independently of the interval count of the overall fit, whose rule
#'   guards intervals, not cells. The counts used are returned as the
#'   \code{n_groups} attribute.
#' @param p_adjust Multiplicity adjustment method passed to
#'   \code{\link[stats]{p.adjust}}; \code{"BH"} (default) controls the false
#'   discovery rate, \code{"holm"} or \code{"bonferroni"} the familywise
#'   error rate.
#' @param alpha Significance level applied to the adjusted probabilities.
#' @return A data frame per item and factor with the full two-way table
#'   (Andrich and Marais 2019, ch. 16): the group main effect (uniform
#'   DIF), the
#'   group-by-interval interaction (non-uniform DIF), and the
#'   class-interval main effect, each with its F statistic; raw and
#'   adjusted probabilities and flags for the two DIF terms; and partial
#'   eta-squared effect sizes (\code{eta2_uniform}, \code{eta2_nonuniform}),
#'   the proportion of residual-plus-effect variance the term accounts for.
#'   For DIF magnitude on the logit scale, where practical significance is
#'   judged, see \code{\link{dif_size}}.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' dif_anova(rasch(X), factors = data.frame(group = g))
#' @export
dif_anova <- function(fit, factors = NULL, n_groups = NULL, p_adjust = "BH",
                      alpha = 0.05) {
  Z <- fit$residuals; L <- ncol(Z)
  factors <- .dif_factors(fit, factors)

  res <- list(); ng_used <- integer(0)
  for (fname in names(factors)) {
    grp <- factor(factors[[fname]])
    ng_f <- if (is.null(n_groups)) .dif_n_groups(fit, grp) else n_groups
    ng_used[fname] <- ng_f
    ci <- .dif_class_intervals(fit, ng_f)
    out <- data.frame(factor = fname, item = colnames(Z),
                      F_uniform = NA_real_, p_uniform = NA_real_,
                      eta2_uniform = NA_real_,
                      F_nonuniform = NA_real_, p_nonuniform = NA_real_,
                      eta2_nonuniform = NA_real_,
                      F_class = NA_real_, p_class = NA_real_)
    for (i in seq_len(L)) {
      d <- data.frame(z = Z[, i], g = grp, ci = ci)
      d <- d[stats::complete.cases(d), ]
      if (nrow(d) < 10 || length(unique(d$g)) < 2) next
      a <- tryCatch(stats::anova(stats::lm(z ~ g * ci, data = d)),
                    error = function(e) NULL)
      if (is.null(a)) next
      rn <- rownames(a)
      ss_res <- if ("Residuals" %in% rn) a["Residuals", "Sum Sq"] else NA_real_
      peta <- function(term) a[term, "Sum Sq"] / (a[term, "Sum Sq"] + ss_res)
      if ("g" %in% rn) {
        out$F_uniform[i] <- a["g", "F value"]; out$p_uniform[i] <- a["g", "Pr(>F)"]
        out$eta2_uniform[i] <- peta("g")
      }
      if ("g:ci" %in% rn) {
        out$F_nonuniform[i] <- a["g:ci", "F value"]
        out$p_nonuniform[i] <- a["g:ci", "Pr(>F)"]
        out$eta2_nonuniform[i] <- peta("g:ci")
      }
      if ("ci" %in% rn) {
        out$F_class[i] <- a["ci", "F value"]; out$p_class[i] <- a["ci", "Pr(>F)"]
      }
    }
    out$p_uniform_adj <- p.adjust(out$p_uniform, method = p_adjust)
    out$p_nonuniform_adj <- p.adjust(out$p_nonuniform, method = p_adjust)
    out$uniform_DIF <- !is.na(out$p_uniform_adj) & out$p_uniform_adj < alpha
    out$nonuniform_DIF <- !is.na(out$p_nonuniform_adj) &
      out$p_nonuniform_adj < alpha
    res[[fname]] <- out
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  attr(out, "n_groups") <- ng_used
  out
}

# variables of an ANOVA term label, e.g. "g1:ci" -> c("g1", "ci")
.term_vars <- function(term) strsplit(term, ":", fixed = TRUE)[[1]]

#' Factorial DIF analysis with Tukey comparisons
#'
#' Models all nominated person factors jointly: for each item the
#' standardised residuals are analysed by the full factorial of the person
#' factors crossed with the trait class interval,
#' \code{z ~ (f1 * f2 * ...) * ci}. Terms not involving the class interval
#' are uniform DIF effects (main effects and factor-by-factor interactions);
#' terms involving it are non-uniform. Probabilities are adjusted across
#' items within each term (Benjamini-Hochberg by default). A significant
#' interaction supersedes the main effects (and lower-order interactions) of
#' the factors it involves, which is recorded in the \code{superseded}
#' column; interpret the highest-order significant terms. Tukey HSD
#' comparisons are returned for every significant, non-superseded group term
#' (the cell-mean contrasts for interactions), with Tukey's own familywise
#' adjustment within each term.
#'
#' Sums of squares are sequential (factors in the order given, class
#' interval last), as is conventional for this residual diagnostic; with
#' markedly unbalanced groups the term order matters and the factor-at-a-time
#' \code{\link{dif_anova}} is a useful cross-check.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors As in \code{\link{dif_anova}}; at least one factor, usually
#'   two or more.
#' @param n_groups Number of trait class intervals. By default set from
#'   the smallest factor-combination cell so every interval-by-cell count
#'   keeps about 30 expected responses (between 2 and 10 intervals); the
#'   value used is returned as \code{n_groups}.
#' @param p_adjust Multiplicity adjustment across items within each term;
#'   default \code{"BH"}.
#' @param alpha Significance level applied to the adjusted probabilities.
#' @param effects \code{"factorial"} (default) crosses every person factor
#'   with every other and with the class interval; \code{"main"} fits the
#'   factors additively (each factor's main effect and its interaction with
#'   the class interval, but no factor-by-factor terms).
#' @param sizes Also compute DIF magnitudes in logits (\code{\link{dif_size}})
#'   for every significant, non-superseded group term: the item is resolved
#'   by the term's levels (interaction terms by their cells) and all
#'   pairwise location differences are returned with Holm familywise
#'   adjustment and the practical-significance flag. Each size involves a
#'   re-analysis, so this costs one refit per flagged item-term.
#' @return A list with \code{summary}, the compact reading of the analysis
#'   (one row per item and group term with the uniform F, adjusted p, and
#'   partial eta-squared -- the term itself -- and the non-uniform ones --
#'   the term crossed with class interval -- plus \code{uniform_DIF},
#'   \code{nonuniform_DIF} and \code{superseded} flags); \code{terms},
#'   the complete per-item analysis of
#'   variance table (term, df, sum of squares, mean square, F, partial
#'   eta-squared, raw and adjusted p, significance, supersession, including
#'   the residual row);
#'   and \code{tukey} (per item, term, and level comparison: difference,
#'   95 per cent interval, and Tukey-adjusted p), plus the \code{alpha} and
#'   adjustment used. Tukey comparisons are reported for significant,
#'   non-superseded group terms except two-level main effects, where the
#'   F test is already the only comparison. With \code{sizes = TRUE},
#'   \code{sizes} holds the logit DIF magnitudes per item, term, and level
#'   pair (two-level main effects included, since the single difference is
#'   exactly the DIF size).
#' @examples
#' set.seed(1); n <- 800
#' d <- seq(-1.5, 1.5, length.out = 6)
#' g1 <- rep(c("a", "b"), each = n / 2)
#' g2 <- rep(c("x", "y"), times = n / 2)
#' sh <- matrix(0, n, 6); sh[g1 == "b", 2] <- 0.8
#' X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
#' dif_anova_factorial(fit)$terms
#' @export
dif_anova_factorial <- function(fit, factors = NULL, n_groups = NULL,
                                p_adjust = "BH", alpha = 0.05,
                                effects = c("factorial", "main"),
                                sizes = FALSE) {
  effects <- match.arg(effects)
  Z <- fit$residuals; L <- ncol(Z)
  factors <- .dif_factors(fit, factors)
  if (is.null(n_groups)) {
    cells <- interaction(factors, drop = TRUE)
    n_groups <- .dif_n_groups(fit, cells)
  }
  ci <- .dif_class_intervals(fit, n_groups)

  fnames <- names(factors)
  safe <- paste0("f", seq_along(fnames))           # syntactic stand-ins
  op <- if (effects == "factorial") " * " else " + "
  form <- stats::as.formula(paste("z ~ (", paste(safe, collapse = op), ") * ci"))

  fits <- vector("list", L); rows <- list()
  for (i in seq_len(L)) {
    d <- data.frame(z = Z[, i], ci = ci)
    for (j in seq_along(fnames)) d[[safe[j]]] <- factor(factors[[fnames[j]]])
    d <- d[stats::complete.cases(d), ]
    if (nrow(d) < 10 || any(vapply(safe, function(s)
      length(unique(d[[s]])) < 2, TRUE))) next
    a <- tryCatch(stats::aov(form, data = d), error = function(e) NULL)
    if (is.null(a)) next
    fits[[i]] <- a
    sm <- summary(a)[[1]]
    rows[[length(rows) + 1L]] <- data.frame(
      item = colnames(Z)[i], term = trimws(rownames(sm)), df = sm$Df,
      sum_sq = sm$`Sum Sq`, mean_sq = sm$`Mean Sq`,
      F_value = sm$`F value`, p = sm$`Pr(>F)`)
  }
  if (!length(rows)) stop("no item yielded an estimable factorial ANOVA")
  terms <- do.call(rbind, rows)
  rownames(terms) <- NULL

  # partial eta-squared per term within each item's table
  terms$eta2_partial <- NA_real_
  for (it in unique(terms$item)) {
    sel <- terms$item == it
    ss_res <- terms$sum_sq[sel & terms$term == "Residuals"]
    if (length(ss_res) == 1)
      terms$eta2_partial[sel & terms$term != "Residuals"] <-
        terms$sum_sq[sel & terms$term != "Residuals"] /
        (terms$sum_sq[sel & terms$term != "Residuals"] + ss_res)
  }

  # adjust across items within each term (the residual rows carry no test)
  terms$p_adj <- NA_real_
  for (tt in setdiff(unique(terms$term), "Residuals")) {
    sel <- terms$term == tt
    terms$p_adj[sel] <- p.adjust(terms$p[sel], method = p_adjust)
  }
  terms$significant <- !is.na(terms$p_adj) & terms$p_adj < alpha

  # a significant higher-order interaction supersedes the lower-order terms
  # built from a subset of its variables (within the same item)
  terms$superseded <- FALSE
  for (it in unique(terms$item)) {
    sel <- which(terms$item == it & terms$significant)
    if (length(sel) < 2) next
    vlist <- lapply(terms$term[sel], .term_vars)
    for (a_i in seq_along(sel)) for (b_i in seq_along(sel)) {
      if (a_i == b_i) next
      if (all(vlist[[a_i]] %in% vlist[[b_i]]) &&
          length(vlist[[a_i]]) < length(vlist[[b_i]]))
        terms$superseded[sel[a_i]] <- TRUE
    }
  }

  # Tukey HSD for significant, non-superseded terms that do not involve the
  # class interval (the group structure itself)
  tk <- list()
  for (i in seq_len(L)) {
    a <- fits[[i]]; if (is.null(a)) next
    it <- colnames(Z)[i]
    cand <- terms[terms$item == it & terms$significant & !terms$superseded, ]
    # group terms only; and no comparisons for a two-level main effect,
    # where the F test is already the only contrast
    keep_t <- !vapply(cand$term, function(tt) "ci" %in% .term_vars(tt), TRUE) &
      !(cand$df == 1L & !grepl(":", cand$term, fixed = TRUE))
    cand <- cand$term[keep_t]
    if (!length(cand)) next
    th <- tryCatch(stats::TukeyHSD(a, which = cand), error = function(e) NULL)
    if (is.null(th)) next
    for (tt in names(th)) {
      tb <- as.data.frame(th[[tt]])
      tk[[length(tk) + 1L]] <- data.frame(
        item = it, term = tt, comparison = rownames(tb),
        difference = tb$diff, lower = tb$lwr, upper = tb$upr,
        p_tukey = tb$`p adj`, row.names = NULL)
    }
  }
  tukey <- if (length(tk)) do.call(rbind, tk) else
    data.frame(item = character(), term = character(),
               comparison = character(), difference = numeric(),
               lower = numeric(), upper = numeric(), p_tukey = numeric())

  # map the syntactic stand-ins back to the nominated factor names
  relabel <- function(x) {
    for (j in rev(seq_along(fnames)))
      x <- gsub(paste0("\\bf", j, "\\b"), fnames[j], x)
    x
  }
  terms$term <- relabel(terms$term)
  tukey$term <- relabel(tukey$term)

  # DIF magnitudes in logits for the significant, non-superseded group
  # terms (interaction terms resolved by their cells)
  size_tab <- NULL
  if (isTRUE(sizes)) {
    sz <- list()
    cand <- terms[terms$significant & !terms$superseded &
                  !vapply(terms$term, function(tt)
                    "ci" %in% .term_vars(tt), TRUE), , drop = FALSE]
    for (r in seq_len(nrow(cand))) {
      it <- cand$item[r]; tt <- cand$term[r]
      ds <- tryCatch(dif_size(fit, it, by = .term_vars(tt)),
                     error = function(e) NULL)
      if (is.null(ds)) next
      p <- ds$pairs
      sz[[length(sz) + 1L]] <- data.frame(item = it, term = tt, p,
                                          row.names = NULL)
    }
    size_tab <- if (length(sz)) do.call(rbind, sz) else
      data.frame(item = character(), term = character())
  }

  # compact reading: one row per item and group term, its own effect being
  # uniform DIF and its crossing with the class interval non-uniform DIF
  gterms <- setdiff(unique(terms$term), "Residuals")
  gterms <- gterms[!vapply(gterms, function(tt)
    "ci" %in% .term_vars(tt), TRUE)]
  srows <- list()
  for (it in unique(terms$item)) for (tt in gterms) {
    u <- terms[terms$item == it & terms$term == tt, , drop = FALSE]
    nu <- terms[terms$item == it & terms$term == paste0(tt, ":ci"), ,
                drop = FALSE]
    if (!nrow(u)) next
    srows[[length(srows) + 1L]] <- data.frame(
      item = it, term = tt,
      F_uniform = u$F_value, p_uniform = u$p,
      p_uniform_adj = u$p_adj, eta2_uniform = u$eta2_partial,
      uniform_DIF = isTRUE(u$significant),
      F_nonuniform = if (nrow(nu)) nu$F_value else NA_real_,
      p_nonuniform = if (nrow(nu)) nu$p else NA_real_,
      p_nonuniform_adj = if (nrow(nu)) nu$p_adj else NA_real_,
      eta2_nonuniform = if (nrow(nu)) nu$eta2_partial else NA_real_,
      nonuniform_DIF = nrow(nu) > 0 && isTRUE(nu$significant),
      superseded = isTRUE(u$superseded))
  }
  summary_tab <- do.call(rbind, srows)
  rownames(summary_tab) <- NULL

  out <- list(summary = summary_tab, terms = terms, tukey = tukey,
              n_groups = n_groups, alpha = alpha, p_adjust = p_adjust)
  if (isTRUE(sizes)) out$sizes <- size_tab
  out
}

#' DIF magnitude in logits with pairwise comparisons
#'
#' Quantifies differential item functioning on the measurement scale
#' itself, where practical significance is judged: the item is resolved
#' into one copy per group (or per cell of a factor combination), the
#' model is refitted, and the distance between the resolved locations is
#' the DIF size in logits (Andrich & Marais 2019, ch. 16: a simulated
#' shift of 0.71 was recovered as 0.75 by exactly this method). Every
#' pair of levels is compared with a Wald test using the full sandwich
#' covariance of the resolved locations (the persons behind different
#' levels are disjoint, but the shared calibration of the other items
#' still couples the estimates, so the covariance is used rather than
#' assumed zero), with familywise adjustment over the pairs. Differences
#' at least \code{flag_logits} in absolute size are flagged as practically
#' significant; half a logit is a common working criterion, to be weighed
#' against the test's targeting and purpose.
#'
#' For an interaction, supply several factor names: levels are then the
#' factor-combination cells, which is the post-hoc follow-up to a
#' significant factor-by-factor term in \code{\link{dif_anova_factorial}}.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or index.
#' @param by One or more person-factor names nominated in the fit (several
#'   names give interaction cells), or a grouping vector/data frame with
#'   one entry per person.
#' @param p_adjust Familywise adjustment over the pairwise comparisons;
#'   default \code{"holm"}.
#' @param alpha Significance level for the adjusted probabilities.
#' @param flag_logits Absolute difference flagged as practically
#'   significant.
#' @param min_n Levels with fewer responders to the item are dropped (their
#'   resolved locations would be too unstable to compare), with a note.
#' @return A list of class \code{"rmt_dif_size"}: \code{levels} (resolved
#'   location and SE per level, with its n), \code{pairs} (per comparison:
#'   difference in logits, SE, z, raw and adjusted p, 95 per cent interval,
#'   \code{significant}, \code{practical}), the settings, and any notes.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' dif_size(fit, "I3", by = "grp")
#' @export
dif_size <- function(fit, item, by, p_adjust = "holm", alpha = 0.05,
                     flag_logits = 0.5, min_n = 20) {
  i <- .item_idx(fit, item)
  if (is.na(i)) stop("no such item")
  item <- fit$items$item[i]
  if (is.character(by) && length(by) < nrow(fit$X)) {
    bad <- if (is.null(fit$factors)) by else setdiff(by, names(fit$factors))
    if (length(bad))
      stop("not a person factor nominated in the fit: ",
           paste(bad, collapse = ", "))
  }
  factors <- .dif_factors(fit, by)
  grp <- if (ncol(factors) == 1L) factor(factors[[1]])
         else interaction(factors, sep = ":", drop = TRUE)
  notes <- character(0)

  # drop levels too thin on this item to resolve
  n_lev <- table(grp[!is.na(fit$X[, i]) & !is.na(grp)])
  thin <- names(n_lev)[n_lev < min_n]
  if (length(thin)) {
    notes <- c(notes, sprintf("level(s) dropped with fewer than %d responders: %s",
                              min_n, paste(thin, collapse = ", ")))
    grp <- factor(ifelse(as.character(grp) %in% thin, NA, as.character(grp)))
  }
  if (nlevels(droplevels(grp)) < 2)
    stop("fewer than two usable levels for item ", item)
  grp <- droplevels(grp)

  refit <- split_items(fit, item, by = grp)
  levs <- levels(grp)
  split_names <- paste0(item, " (", levs, ")")
  idx <- match(split_names, refit$items$item)
  if (anyNA(idx))
    stop("resolved item(s) missing after re-analysis (too little data): ",
         paste(split_names[is.na(idx)], collapse = ", "))

  # location covariance from the sandwich: var(mean of a threshold block)
  thr <- refit$thresholds; cv <- refit$est$cov_tau
  block <- lapply(idx, function(k) thr$id[thr$item == k])
  loc <- refit$items$location[idx]
  vloc <- matrix(NA_real_, length(levs), length(levs))
  for (a in seq_along(levs)) for (b in seq_along(levs))
    vloc[a, b] <- mean(cv[block[[a]], block[[b]], drop = FALSE])

  n_item <- as.integer(table(grp[!is.na(fit$X[, i]) & !is.na(grp)])[levs])
  levels_df <- data.frame(level = levs, location = loc,
                          se = sqrt(pmax(diag(vloc), 0)), n = n_item)

  pr <- t(utils::combn(seq_along(levs), 2))
  pairs <- data.frame(
    level_a = levs[pr[, 1]], level_b = levs[pr[, 2]],
    difference = loc[pr[, 1]] - loc[pr[, 2]],
    se = sqrt(pmax(diag(vloc)[pr[, 1]] + diag(vloc)[pr[, 2]] -
                   2 * vloc[cbind(pr[, 1], pr[, 2])], 1e-12)))
  pairs$z <- pairs$difference / pairs$se
  pairs$p <- 2 * pnorm(-abs(pairs$z))
  pairs$p_adj <- p.adjust(pairs$p, method = p_adjust)
  pairs$lower <- pairs$difference - qnorm(0.975) * pairs$se
  pairs$upper <- pairs$difference + qnorm(0.975) * pairs$se
  pairs$significant <- pairs$p_adj < alpha
  pairs$practical <- abs(pairs$difference) >= flag_logits

  out <- list(item = item, by = paste(names(factors), collapse = ":"),
              levels = levels_df, pairs = pairs, alpha = alpha,
              p_adjust = p_adjust, flag_logits = flag_logits, notes = notes)
  class(out) <- "rmt_dif_size"
  out
}

#' @export
print.rmt_dif_size <- function(x, ...) {
  cat(sprintf("DIF size for %s by %s (resolved locations, logits)\n",
              x$item, x$by))
  lv <- x$levels; lv[-1] <- lapply(lv[-1], round, 3)
  print(lv, row.names = FALSE)
  pr <- x$pairs
  num <- vapply(pr, is.numeric, TRUE)
  pr[num] <- lapply(pr[num], round, 3)
  pr$significant <- ifelse(pr$significant, "*", "")
  pr$practical <- ifelse(pr$practical, sprintf(">= %.2f", x$flag_logits), "")
  print(pr, row.names = FALSE)
  cat(sprintf("p adjusted by %s over %d pairwise comparison(s); practical criterion %.2f logits\n",
              x$p_adjust, nrow(pr), x$flag_logits))
  if (length(x$notes)) cat("notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Planned contrasts: the confirmatory alternative to exhaustive post-hoc
# pairwise comparison. The family of questions is derived from the structure
# of the nominated factors (or supplied), estimated on the logit scale from
# resolved item locations, and -- when persons repeat across rows of a
# stacked design -- tested from person-level contrast scores so that
# within-subject dependence is respected.
# ---------------------------------------------------------------------------

# Resolve one item over grouping cells: locations and sandwich covariance.
.dif_resolve <- function(fit, item, grp, min_n) {
  i <- .item_idx(fit, item)
  item <- fit$items$item[i]
  notes <- character(0)
  n_lev <- table(grp[!is.na(fit$X[, i]) & !is.na(grp)])
  thin <- names(n_lev)[n_lev < min_n]
  if (length(thin)) {
    notes <- c(notes, sprintf(
      "%s: level(s) dropped with fewer than %d responders: %s",
      item, min_n, paste(thin, collapse = ", ")))
    grp <- factor(ifelse(as.character(grp) %in% thin, NA, as.character(grp)))
  }
  grp <- droplevels(grp)
  if (nlevels(grp) < 2) return(NULL)
  refit <- split_items(fit, item, by = grp)
  levs <- levels(grp)
  idx <- match(paste0(item, " (", levs, ")"), refit$items$item)
  if (anyNA(idx)) return(NULL)
  thr <- refit$thresholds; cv <- refit$est$cov_tau
  block <- lapply(idx, function(k) thr$id[thr$item == k])
  loc <- refit$items$location[idx]
  vloc <- matrix(NA_real_, length(levs), length(levs))
  for (a in seq_along(levs)) for (b in seq_along(levs))
    vloc[a, b] <- mean(cv[block[[a]], block[[b]], drop = FALSE])
  list(levs = levs, loc = loc, vloc = vloc, notes = notes)
}

# A factor is treated as ordered when declared ordered or when its levels
# parse as numbers (ages, waves, doses).
.dif_is_ordered <- function(f)
  is.ordered(f) || !any(is.na(suppressWarnings(as.numeric(levels(f)))))

# The leading contrast of a factor: the difference for two levels, the
# linear trend for an ordered factor, none for a nominal many-level factor.
.dif_leading <- function(f) {
  K <- nlevels(f)
  if (K == 2L) {
    w <- c(-1, 1); names(w) <- levels(f)
    list(weights = w,
         label = sprintf("%s - %s", levels(f)[2], levels(f)[1]))
  } else if (.dif_is_ordered(f)) {
    sc <- suppressWarnings(as.numeric(levels(f)))
    cp <- if (!any(is.na(sc))) stats::contr.poly(K, scores = sc)
          else stats::contr.poly(K)
    w <- cp[, 1]; names(w) <- levels(f)
    list(weights = w, label = "linear")
  } else NULL
}

# The planned questions a single factor admits.
.dif_factor_contrasts <- function(f, fname) {
  K <- nlevels(f); out <- list()
  if (K == 2L) {
    lead <- .dif_leading(f)
    out[[sprintf("%s: %s", fname, lead$label)]] <- lead$weights
  } else if (.dif_is_ordered(f)) {
    sc <- suppressWarnings(as.numeric(levels(f)))
    cp <- if (!any(is.na(sc))) stats::contr.poly(K, scores = sc)
          else stats::contr.poly(K)
    w1 <- cp[, 1]; names(w1) <- levels(f)
    out[[sprintf("%s: linear", fname)]] <- w1
    w2 <- cp[, 2]; names(w2) <- levels(f)
    out[[sprintf("%s: quadratic", fname)]] <- w2
  } else if (K <= 4L) {
    pr <- utils::combn(levels(f), 2)
    for (j in seq_len(ncol(pr))) {
      w <- stats::setNames(numeric(K), levels(f))
      w[pr[2, j]] <- 1; w[pr[1, j]] <- -1
      out[[sprintf("%s: %s - %s", fname, pr[2, j], pr[1, j])]] <- w
    }
  } else {
    for (l in levels(f)) {
      w <- stats::setNames(rep(-1 / (K - 1), K), levels(f)); w[l] <- 1
      out[[sprintf("%s: %s - others", fname, l)]] <- w
    }
  }
  out
}

# Scale cell weights so the positive and negative parts each sum to one:
# every contrast then reads as a difference between two weighted averages,
# in logits, comparable across contrasts and against the practical flag.
.dif_norm <- function(w) {
  w[is.na(w)] <- 0
  s <- sum(abs(w))
  if (s < 1e-10) return(NULL)
  w * 2 / s
}

# Spread factor-level weights over the design cells (unweighted marginal
# means: each level's weight is shared equally by the cells carrying it).
.dif_cell_weights <- function(cellmap, fw, fname) {
  lev <- as.character(cellmap[[fname]])
  keep <- lev %in% names(fw)
  w <- stats::setNames(numeric(nrow(cellmap)), cellmap$cell)
  w[keep] <- fw[lev[keep]] / as.numeric(table(lev[keep])[lev[keep]])
  w
}

# Derive the planned family from the factor structure.
.dif_contrast_family <- function(factors, cellmap, within_names) {
  fam <- list(); meta <- list()
  for (fname in names(factors)) {
    fc <- .dif_factor_contrasts(factors[[fname]], fname)
    for (nm in names(fc)) {
      w <- .dif_norm(.dif_cell_weights(cellmap, fc[[nm]], fname))
      if (is.null(w)) next
      fam[[nm]] <- w
      meta[[nm]] <- list(factors = fname, fweights = fc[nm],
                         within = fname %in% within_names)
    }
  }
  fns <- names(factors)
  if (length(fns) >= 2) for (a in seq_len(length(fns) - 1))
    for (b in seq(a + 1, length(fns))) {
      la <- .dif_leading(factors[[fns[a]]])
      lb <- .dif_leading(factors[[fns[b]]])
      if (is.null(la) || is.null(lb)) next
      wa <- la$weights[as.character(cellmap[[fns[a]]])]
      wb <- lb$weights[as.character(cellmap[[fns[b]]])]
      key <- paste(cellmap[[fns[a]]], cellmap[[fns[b]]])
      w <- .dif_norm(stats::setNames(
        wa * wb / as.numeric(table(key)[key]), cellmap$cell))
      if (is.null(w)) next
      nm <- sprintf("%s(%s) x %s(%s)", fns[a], la$label, fns[b], lb$label)
      fam[[nm]] <- w
      meta[[nm]] <- list(factors = c(fns[a], fns[b]),
                         fweights = list(la$weights, lb$weights),
                         within = any(c(fns[a], fns[b]) %in% within_names))
    }
  list(family = fam, meta = meta)
}

# Welch test of a linear combination of independent group means.
.welch_contrast <- function(vals, g, w) {
  ok <- !is.na(vals) & !is.na(g)
  vals <- vals[ok]; g <- droplevels(factor(g[ok]))
  w <- w[levels(g)]
  if (any(is.na(w)) || sum(abs(w)) < 1e-10) return(NULL)
  m <- tapply(vals, g, mean); v <- tapply(vals, g, stats::var)
  n <- tapply(vals, g, length)
  if (any(n < 2)) return(NULL)
  vv <- sum(w^2 * v / n)
  if (!is.finite(vv) || vv <= 0) return(NULL)
  df <- vv^2 / sum((w^2 * v / n)^2 / (n - 1))
  t <- sum(w * m) / sqrt(vv)
  list(stat = t, df = df, p = 2 * stats::pt(-abs(t), df))
}

#' Planned DIF contrasts derived from the factor structure
#'
#' The confirmatory alternative to exhaustive post-hoc comparison: instead
#' of every pair of design cells, a small family of one-degree-of-freedom
#' questions is tested, so familywise control costs little power (Maxwell
#' and Delaney 2004, ch. 5). By default the family is derived from the
#' structure of the factors themselves -- a two-level factor contributes its
#' difference; an ordered factor (declared ordered, or with numeric levels
#' such as ages or waves) contributes its linear and quadratic trends; a
#' nominal factor contributes all pairs when it has up to four levels and
#' each-level-against-the-rest otherwise; and every pair of factors with a
#' leading contrast (a difference or a linear trend) contributes the product
#' interaction. Print the returned object to see the family in words before
#' reading the results; a family endorsed in advance of the results is what
#' makes the contrasts planned.
#'
#' Each contrast is estimated in logits from resolved item locations (the
#' item split into one copy per design cell and the model refitted, as in
#' \code{\link{dif_size}}), with cell weights scaled so every estimate is a
#' difference between two weighted averages -- directly comparable to the
#' practical-significance criterion. Because resolution is used, magnitudes
#' are read from a calibration in which compensating artificial DIF has been
#' removed (Andrich and Hagquist 2015).
#'
#' When \code{id} shows that persons repeat across rows (a stacked
#' repeated-measures design), between-row independence fails and the usual
#' tests would be invalid. Significance is then computed from person-level
#' scores of the standardised residuals: a within-subject contrast (for
#' example a trend over time) becomes one contrast score per person, tested
#' against zero; a between-subjects contrast is tested on person-mean
#' residuals; and a between-by-within interaction tests the person contrast
#' scores across the between groups. Logit estimates are still reported from
#' the resolved locations; their standard errors treat rows as independent
#' and are conservative for within-subject differences.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors A data frame of person factors, a character vector naming
#'   factors nominated in the fit, or a single grouping vector. Defaults to
#'   every factor stored in the fit.
#' @param items Item names or indices to test; all items by default.
#' @param within Names of factors that vary within person (for example
#'   time). Detected automatically when \code{id} is supplied and a factor
#'   varies within an id.
#' @param id Person identifier with one entry per row, or the name of a
#'   nominated factor holding it; required for stacked designs where the
#'   same person occupies several rows.
#' @param contrasts \code{"auto"} (derive the family from the factor
#'   structure) or a named list of numeric cell-weight vectors, each named
#'   by the design-cell labels (factor levels joined by \code{":"}).
#'   Weights are rescaled so the positive and negative parts each sum to
#'   one.
#' @param p_adjust Familywise adjustment over the whole family (items by
#'   contrasts); default \code{"holm"}.
#' @param alpha Significance level for the adjusted probabilities.
#' @param flag_logits Absolute estimate flagged as practically significant.
#' @param min_n Cells with fewer responders to an item are dropped from that
#'   item's resolution, with a note.
#' @return A list of class \code{"rmt_dif_contrasts"}: \code{table} (one row
#'   per item and contrast: estimate in logits, SE, statistic, df where a t
#'   test was used, raw and adjusted p, 95 per cent interval,
#'   \code{significant}, \code{practical}, \code{within}), \code{family}
#'   (the derived questions with their cell weights), the settings, and any
#'   \code{notes}.
#' @references Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing
#'   Experiments and Analyzing Data} (2nd ed.). Mahwah, NJ: Erlbaum.
#'
#'   Andrich, D., & Hagquist, C. (2015). Real and artificial differential
#'   item functioning in polytomous items. \emph{Educational and
#'   Psychological Measurement}, 75(2), 185-207.
#'
#'   Hagquist, C., & Andrich, D. (2017). Recent advances in analysis of
#'   differential item functioning in health research using the Rasch
#'   model. \emph{Health and Quality of Life Outcomes}, 15, 181.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' dif_contrasts(fit, items = c("I3", "I5"))
#' @export
dif_contrasts <- function(fit, factors = NULL, items = NULL, within = NULL,
                          id = NULL, contrasts = "auto", p_adjust = "holm",
                          alpha = 0.05, flag_logits = 0.5, min_n = 20) {
  factors <- .dif_factors(fit, factors)
  factors <- as.data.frame(lapply(factors, function(v) {
    f <- droplevels(if (is.ordered(v)) v else factor(v))
    f
  }), check.names = FALSE, stringsAsFactors = FALSE)
  grp <- interaction(factors, sep = ":", drop = TRUE)
  cellmap <- unique(data.frame(cell = as.character(grp), factors,
                               check.names = FALSE))
  cellmap <- cellmap[match(levels(grp), cellmap$cell), , drop = FALSE]

  if (is.character(id) && length(id) == 1L && !is.null(fit$factors) &&
      id %in% names(fit$factors)) id <- fit$factors[[id]]
  if (is.null(within) && !is.null(id) && anyDuplicated(id)) {
    within <- names(factors)[vapply(names(factors), function(fn)
      any(tapply(as.character(factors[[fn]]), id,
                 function(v) length(unique(v)) > 1L)), TRUE)]
  }
  if (is.null(within)) within <- character(0)
  within <- intersect(within, names(factors))
  if (length(within) && is.null(id))
    stop("`within` requires `id` to pair a person's rows")
  paired <- !is.null(id) && anyDuplicated(id) > 0

  if (identical(contrasts, "auto")) {
    fam <- .dif_contrast_family(factors, cellmap, within)
  } else {
    if (!is.list(contrasts) || is.null(names(contrasts)))
      stop("`contrasts` must be \"auto\" or a named list of cell weights")
    fam <- list(family = list(), meta = list())
    for (nm in names(contrasts)) {
      w <- contrasts[[nm]]
      if (is.null(names(w)) || !all(names(w) %in% cellmap$cell))
        stop("weights of contrast '", nm, "' must be named by design cells: ",
             paste(cellmap$cell, collapse = ", "))
      full <- stats::setNames(numeric(nrow(cellmap)), cellmap$cell)
      full[names(w)] <- w
      w <- .dif_norm(full)
      if (is.null(w)) stop("contrast '", nm, "' has no weight")
      fam$family[[nm]] <- w
      fam$meta[[nm]] <- list(factors = names(factors), fweights = NULL,
                             within = FALSE)
    }
  }
  if (!length(fam$family)) stop("no contrasts could be formed")

  its <- if (is.null(items)) fit$items$item else
    fit$items$item[vapply(items, function(x) .item_idx(fit, x), 1L)]
  Z <- fit$residuals
  notes <- character(0)
  rows <- list()

  for (item in its) {
    i <- .item_idx(fit, item)
    rs <- .dif_resolve(fit, item, grp, min_n)
    if (!is.null(rs)) notes <- c(notes, rs$notes)
    for (nm in names(fam$family)) {
      w_full <- fam$family[[nm]]
      mt <- fam$meta[[nm]]
      est <- se <- stat <- df <- p <- NA_real_
      if (!is.null(rs)) {
        w <- w_full[rs$levs]
        w[is.na(w)] <- 0
        if (sum(w > 0) > 0 && sum(w < 0) > 0 &&
            abs(sum(abs(w)) - 2) < 0.5) {   # cells mostly intact
          w <- .dif_norm(w)
          est <- sum(w * rs$loc)
          se <- sqrt(max(drop(t(w) %*% rs$vloc %*% w), 1e-12))
        }
      }
      if (!paired) {
        if (is.finite(est)) {
          stat <- est / se
          p <- 2 * stats::pnorm(-abs(stat))
        }
      } else if (isTRUE(mt$within) && length(mt$factors) == 1L) {
        # within-subject contrast: one score per (complete) person
        fw <- .dif_norm(fam$meta[[nm]]$fweights[[1]])
        lev <- as.character(factors[[mt$factors]])
        zi <- Z[, i]
        ps <- split(seq_along(zi), id)
        psi <- vapply(ps, function(rws) {
          l <- lev[rws]
          if (anyDuplicated(l) || !all(names(fw) %in% l)) return(NA_real_)
          sum(fw * zi[rws][match(names(fw), l)])
        }, 0)
        psi <- psi[!is.na(psi)]
        if (length(psi) >= 10) {
          # sign aligned with the logit estimate: a harder level has the
          # higher resolved location but the lower residuals
          stat <- -mean(psi) / (stats::sd(psi) / sqrt(length(psi)))
          df <- length(psi) - 1
          p <- 2 * stats::pt(-abs(stat), df)
        }
      } else if (isTRUE(mt$within) && length(mt$factors) == 2L) {
        # between-by-within interaction: within contrast scores per person,
        # tested across the between groups
        wf <- intersect(mt$factors, within)[1]
        bf <- setdiff(mt$factors, wf)[1]
        fws <- fam$meta[[nm]]$fweights
        fw <- .dif_norm(fws[[match(wf, mt$factors)]])
        bw <- fws[[match(bf, mt$factors)]]
        lev <- as.character(factors[[wf]])
        blev <- as.character(factors[[bf]])
        zi <- Z[, i]
        ps <- split(seq_along(zi), id)
        psi <- vapply(ps, function(rws) {
          l <- lev[rws]
          if (anyDuplicated(l) || !all(names(fw) %in% l)) return(NA_real_)
          sum(fw * zi[rws][match(names(fw), l)])
        }, 0)
        pb <- vapply(ps, function(rws) blev[rws][1], "")
        ok <- !is.na(psi)
        wc <- .welch_contrast(psi[ok], pb[ok], bw / 2)
        if (!is.null(wc)) { stat <- -wc$stat; df <- wc$df; p <- wc$p }
      } else {
        # between-subjects question in a stacked design: test on person
        # means of the residuals so each person counts once
        zi <- Z[, i]
        pm <- tapply(zi, id, mean, na.rm = TRUE)
        pcell <- tapply(as.character(grp), id, function(v) v[1])
        wc <- .welch_contrast(pm, factor(pcell, levels = cellmap$cell),
                              w_full / 2)
        if (!is.null(wc)) { stat <- -wc$stat; df <- wc$df; p <- wc$p }
      }
      rows[[length(rows) + 1L]] <- data.frame(
        item = item, contrast = nm, within = isTRUE(mt$within),
        estimate = est, se = se, statistic = stat, df = df, p = p)
    }
  }
  tab <- do.call(rbind, rows)
  tab$p_adj <- stats::p.adjust(tab$p, method = p_adjust)
  tab$lower <- tab$estimate - stats::qnorm(0.975) * tab$se
  tab$upper <- tab$estimate + stats::qnorm(0.975) * tab$se
  tab$significant <- !is.na(tab$p_adj) & tab$p_adj < alpha
  tab$practical <- !is.na(tab$estimate) & abs(tab$estimate) >= flag_logits
  rownames(tab) <- NULL

  fam_df <- data.frame(
    contrast = names(fam$family),
    within = vapply(fam$meta, function(m) isTRUE(m$within), TRUE),
    cells = vapply(fam$family, function(w)
      paste(sprintf("%s %+0.2f", names(w)[w != 0], w[w != 0]),
            collapse = ", "), ""))
  rownames(fam_df) <- NULL

  out <- list(table = tab, family = fam_df, within = within,
              paired = paired, alpha = alpha, p_adjust = p_adjust,
              flag_logits = flag_logits, notes = unique(notes))
  class(out) <- "rmt_dif_contrasts"
  out
}

#' @export
print.rmt_dif_contrasts <- function(x, ...) {
  cat("Planned DIF contrasts (", nrow(x$family), " questions x ",
      length(unique(x$table$item)), " items; ", x$p_adjust,
      " over the family)\n", sep = "")
  for (r in seq_len(nrow(x$family)))
    cat(sprintf("  %s%s\n", x$family$contrast[r],
                if (x$family$within[r]) "  [within subjects]" else ""))
  if (x$paired)
    cat("Stacked design: tests use person-level residual scores;",
        "logit SEs are conservative for within contrasts.\n")
  cat("\n")
  tab <- x$table
  show <- tab[, c("item", "contrast", "estimate", "se", "statistic", "p_adj",
                  "significant", "practical")]
  print(.fmt_df(show), row.names = FALSE)
  if (length(x$notes)) cat("\n", paste(x$notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
