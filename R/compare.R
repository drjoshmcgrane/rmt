# rasch :: model comparison
# ===========================================================================
# Side-by-side comparison of fitted models. Three kinds of evidence are
# reported. (1) The pairwise conditional log-likelihood with the number of
# structural parameters and, for fits of the SAME response data (identical
# items, categories, and persons, hence identical conditional information),
# twice the log-likelihood difference from the reference fit. Because the
# likelihood is a composite (pairwise) one, the difference is descriptive
# and is not chi-square calibrated; it is most meaningful for nested
# structures (for example the rating scale model inside the partial credit
# model, or equal units inside the extended frame of reference model).
# (2) Composite-likelihood information criteria that ARE calibrated for the
# pairwise over-counting: CL-AIC (Varin & Vidoni 2005) and CL-BIC (Gao &
# Song 2010) penalise -2 cl with the effective parameter count
# tr(H^-1 J) from the Godambe matrices -- the same quantity whose
# eigenvalues calibrate lr_test() -- instead of the nominal count, with
# the CL-BIC log(n) taken over the independent units (persons; judges for
# paired comparisons). Smaller is better, valid across models of the same
# data whether or not they nest.
# (3) Calibration-free fit descriptors that remain comparable across
# different data preparations: the total item-trait chi-square against its
# degrees of freedom, the spread of the item and person fit residuals
# (ideal SD 1), the person separation index, and Cronbach's alpha.
# ===========================================================================

# Composite-likelihood information criteria for one fit: effective parameter
# count tr(H^-1 J) (sign-convention free via abs: the eigenvalues of H^-1 J
# share one sign), CL-AIC, CL-BIC. NA when the fit does not carry its
# Godambe matrices (MFRM and EFRM assemble their own estimation structures).
.cl_ic <- function(f) {
  if (inherits(f, "rasch_btl")) {
    if (is.null(f$cl)) return(c(eff = NA_real_, aic = NA_real_, bic = NA_real_))
    eff <- f$cl$eff_params; n <- f$cl$n_units; ll <- f$loglik
  } else {
    est <- f$est
    if (is.null(est$cov_beta) || is.null(est$H_beta))
      return(c(eff = NA_real_, aic = NA_real_, bic = NA_real_))
    eff <- abs(sum(diag(est$cov_beta %*% est$H_beta)))
    # independent units: persons contributing at least one informative pair
    # (extreme persons condition every pair onto a boundary total)
    n <- sum(!f$person$extreme & rowSums(!is.na(f$X)) >= 2L)
    ll <- est$loglik
  }
  c(eff = eff, aic = -2 * ll + 2 * eff, bic = -2 * ll + log(n) * eff)
}

#' Compare fitted Rasch models
#'
#' Builds a comparison table for two or more fits from \code{\link{rasch}},
#' \code{\link{rasch_mfrm}}, \code{\link{rasch_efrm}}, or (all together)
#' \code{\link{btl}}. For fits of the
#' same response data (identical item columns, maximum scores, and number of
#' persons) the pairwise conditional log-likelihoods share their conditional
#' information, and twice the difference from the reference fit is reported
#' with the difference in parameter counts; this is descriptive (composite
#' likelihood), and most meaningful for nested structures such as RSM inside
#' PCM.
#'
#' The calibrated comparison is carried by the composite-likelihood
#' information criteria \code{cl_aic} (Varin and Vidoni 2005) and
#' \code{cl_bic} (Gao and Song 2010): \eqn{-2\,cl + c \cdot tr(H^{-1}J)},
#' with \eqn{c = 2} or \eqn{\log n}. Because every response enters every
#' pair its item forms, the pairwise log-likelihood over-counts the data;
#' the effective parameter count \eqn{tr(H^{-1}J)} from the Godambe
#' matrices -- the same quantity whose eigenvalues calibrate
#' \code{\link{lr_test}} -- absorbs exactly that over-counting, where the
#' nominal parameter count would not. \eqn{n} counts independent units:
#' persons contributing at least one informative pair, or judges for
#' paired-comparison fits (count-weighted comparisons when unclustered).
#' Smaller is better; the criteria are valid across models of the same data
#' whether or not they nest, and are \code{NA} (with the reason in the
#' printed note) for MFRM and EFRM fits, which do not carry their Godambe
#' matrices.
#'
#' Across different data preparations (subtests, splits, facet or frame
#' structures) the likelihoods are not comparable and the calibration-free
#' columns carry the comparison: total item-trait chi-square per degree of
#' freedom, item and person fit residual SDs (ideal 1), PSI, and alpha
#' (OSI for paired comparisons).
#'
#' @param ... Two or more fitted objects, ideally named
#'   (\code{compare_fits(PCM = f1, RSM = f2)}). Either all Rasch-family
#'   fits or all \code{btl} fits; for \code{btl}, fits of the same
#'   comparison data (same objects, comparisons, and judges) support the
#'   likelihood columns -- e.g. free versus principal-component thresholds,
#'   with and without a position effect or within-judge dependence.
#' @param reference Index or name of the reference fit for the
#'   log-likelihood difference; defaults to the first.
#' @return A data frame with one row per fit: label, model, persons, items
#'   (judges, objects, comparisons for \code{btl}),
#'   parameters, log-likelihood, \code{eff_params}, \code{cl_aic},
#'   \code{cl_bic}, comparability with the reference,
#'   \code{two_delta_ll} and \code{delta_parameters} (same-data fits only),
#'   chi-square per df, fit residual SDs, PSI, and alpha (OSI for
#'   \code{btl}).
#' @examples
#' set.seed(1)
#' simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(seq(-1, 1, length.out = 6), function(b)
#'   sapply(th, function(t) sample(0:3, 1, prob = simP(t, b + c(-0.8, 0, 0.8)))))
#' colnames(X) <- paste0("R", 1:6)
#' compare_fits(PCM = rasch(X, model = "PCM"), RSM = rasch(X, model = "RSM"))
#' @export
compare_fits <- function(..., reference = 1) {
  fits <- list(...)
  if (length(fits) < 2) stop("supply at least two fits to compare")
  is_btl <- vapply(fits, inherits, TRUE, what = "rasch_btl")
  bad <- !vapply(fits, inherits, TRUE, what = "rasch") & !is_btl
  if (any(bad)) stop("argument(s) ", paste(which(bad), collapse = ", "),
                     " are not rasch or btl fits")
  if (any(is_btl) && !all(is_btl))
    stop("compare either Rasch-family fits or btl fits, not a mixture: ",
         "their likelihoods are over different data")
  labs <- names(fits)
  if (is.null(labs)) labs <- rep("", length(fits))
  labs[labs == ""] <- paste0("fit", seq_along(fits))[labs == ""]
  if (is.character(reference)) reference <- match(reference, labs)
  if (is.na(reference) || reference < 1 || reference > length(fits))
    stop("no such reference fit")

  if (all(is_btl)) {
    sig <- function(f) list(objects = sort(f$objects$object),
                            n = f$n_comparisons,
                            judges = if (is.null(f$judges)) 0L
                                     else nrow(f$judges))
    ref_sig <- sig(fits[[reference]])
    rows <- lapply(seq_along(fits), function(i) {
      f <- fits[[i]]
      ic <- .cl_ic(f)
      dep <- if (is.null(f$dependence)) "" else
        paste0(" + ", paste(f$dependence$effect, collapse = " + "))
      data.frame(
        label = labs[i],
        model = paste0("BTL (", if (max(f$m) > 1L)
          paste0("graded, ", f$thr_structure, " thresholds") else
            "dichotomous", dep, ")"),
        judges = if (is.null(f$judges)) NA_integer_ else nrow(f$judges),
        objects = nrow(f$objects), comparisons = f$n_comparisons,
        parameters = if (is.null(f$cl)) NA_integer_ else f$cl$n_parameters,
        loglik = f$loglik,
        eff_params = unname(ic["eff"]), cl_aic = unname(ic["aic"]),
        cl_bic = unname(ic["bic"]),
        same_data = identical(sig(f), ref_sig),
        two_delta_ll = NA_real_, delta_parameters = NA_integer_,
        chisq_per_df = f$total_chisq / f$total_df,
        OSI = f$osi$PSI)
    })
  } else {
    sig <- function(f) list(items = colnames(f$X), m = unname(f$m),
                            n = nrow(f$X))
    ref_sig <- sig(fits[[reference]])
    rows <- lapply(seq_along(fits), function(i) {
      f <- fits[[i]]
      ic <- .cl_ic(f)
      data.frame(
        label = labs[i], model = f$model,
        persons = nrow(f$X), items = ncol(f$X),
        parameters = if (is.null(f$est$n_parameters)) NA_integer_
                     else f$est$n_parameters,
        loglik = f$est$loglik,
        eff_params = unname(ic["eff"]), cl_aic = unname(ic["aic"]),
        cl_bic = unname(ic["bic"]),
        same_data = identical(sig(f), ref_sig),
        two_delta_ll = NA_real_, delta_parameters = NA_integer_,
        chisq_per_df = f$total_chisq / f$total_df,
        item_fit_sd = f$item_fit_summary$sd,
        person_fit_sd = f$person_fit_summary$sd,
        PSI = f$psi$PSI, alpha = f$alpha$alpha)
    })
  }
  out <- do.call(rbind, rows)
  ref <- out[reference, ]
  cmp <- out$same_data & seq_len(nrow(out)) != reference
  out$two_delta_ll[cmp] <- 2 * (out$loglik[cmp] - ref$loglik)
  out$delta_parameters[cmp] <- out$parameters[cmp] - ref$parameters
  rownames(out) <- NULL
  attr(out, "reference") <- labs[reference]
  attr(out, "note") <- paste0(
    "cl_aic and cl_bic are composite-likelihood information criteria ",
    "(Varin & Vidoni 2005; Gao & Song 2010): -2 cl penalised by the ",
    "effective parameter count tr(H^-1 J), which absorbs the pairwise ",
    "over-counting that the nominal count would not; smaller is better, ",
    "valid across models of the same data",
    if (any(!vapply(fits, function(f)
      is.finite(.cl_ic(f)["eff"]), TRUE)))
      " (NA for MFRM/EFRM fits, which do not carry their Godambe matrices)"
    else "",
    ". two_delta_ll is the raw composite difference against the reference, ",
    "descriptive only. Across different data preparations compare ",
    "chisq_per_df, the fit residual SDs (ideal 1), and the ",
    "separation/reliability columns.")
  class(out) <- c("rasch_compare", "data.frame")
  out
}

#' @export
print.rasch_compare <- function(x, ...) {
  cat(sprintf("Model comparison (reference: %s)\n\n", attr(x, "reference")))
  core <- c("label", "model", "persons", "items", "judges", "objects",
            "eff_params", "cl_aic", "cl_bic", "two_delta_ll",
            "chisq_per_df", "item_fit_sd", "person_fit_sd", "PSI", "alpha",
            "OSI")
  y <- as.data.frame(x)[, intersect(core, names(x)), drop = FALSE]
  print(.fmt_df(y), row.names = FALSE)
  cat("(further columns on the object: loglik, parameters, same_data)\n")
  cat("\n", attr(x, "note"), "\n", sep = "")
  invisible(x)
}

#' Likelihood-ratio test of the partial credit against the rating scale model
#'
#' A likelihood-ratio test in the tradition of Andersen (1973): an
#' unrestricted (partial credit)
#' analysis is compared with the rating re-parameterisation of the same
#' model on the same data. Twice the difference in the pairwise conditional
#' log-likelihoods is referred to a chi-square on the difference in the
#' number of threshold parameters. A non-significant outcome supports
#' adopting the simpler rating parameterisation.
#'
#' The likelihood here is the pairwise composite
#' likelihood, not a full likelihood, and twice its difference is not
#' chi-square distributed: each response enters every pair its item forms,
#' so the raw statistic is inflated. Two statistics are therefore reported.
#' \code{chisq} is the raw composite value with its naive \code{p}, the
#' conventional display. The limiting
#' law of the raw statistic is \eqn{\sum_j \lambda_j \chi^2_1} (Kent 1982;
#' Varin, Reid and Firth 2011) with \eqn{\lambda_j} the eigenvalues of
#' \eqn{(C'H^{-1}C)^{-1}\,C'H^{-1}JH^{-1}C} over the \eqn{r} constrained
#' directions \eqn{C} (the part of the partial-credit threshold space
#' outside the rating subspace), estimated from the same Godambe \eqn{H}
#' and \eqn{J} matrices that supply the sandwich standard errors; matching
#' the mean gives \code{chisq_adj} \eqn{= r W / \sum_j \lambda_j} on
#' \eqn{r} degrees of freedom. Use \code{p_adj} for inference; the naive
#' \code{p} is severely anticonservative and kept only for comparability
#' with conventional software displays.
#'
#' @param fit A \code{"PCM"} fit from \code{\link{rasch}} with equal maximum
#'   scores across items (the rating parameterisation requires them).
#' @param maxit,tol Passed to the rating-scale refit.
#' @return A list of class \code{"rasch_lr"}: raw \code{chisq}, \code{df},
#'   \code{p} (the conventional display); adjusted \code{chisq_adj}, \code{p_adj},
#'   and the eigenvalues \code{lambda}; the two log-likelihoods; and the
#'   rating-scale refit (\code{fit_rsm}).
#' @references Kent, J. T. (1982). Robust properties of likelihood ratio
#'   tests. Biometrika, 69, 19-27. Varin, C., Reid, N. and Firth, D. (2011).
#'   An overview of composite likelihood methods. Statistica Sinica, 21,
#'   5-42.
#' @examples
#' set.seed(1)
#' tau <- c(-0.7, 0.7)
#' X <- sapply(seq(-1, 1, length.out = 6), function(d) vapply(rnorm(300),
#'   function(b) sample(0:2, 1, prob = item_moments(b, tau + d)$P), 0L))
#' colnames(X) <- paste0("Q", 1:6)
#' lr_test(rasch(X, model = "PCM"))
#' @export
lr_test <- function(fit, maxit = 60, tol = 1e-8) {
  if (!identical(fit$model, "PCM"))
    stop("lr_test() compares an unrestricted (PCM) fit with its rating ",
         "re-parameterisation; supply a PCM fit")
  if (length(unique(fit$m)) != 1L)
    stop("the rating parameterisation requires equal maximum scores across items")
  if (max(fit$m) < 2L)
    stop("with dichotomous items the two parameterisations coincide")
  rsm <- rasch(fit$X, model = "RSM", n_groups = fit$n_groups, maxit = maxit,
               tol = tol)
  chisq <- 2 * (fit$est$loglik - rsm$est$loglik)
  df <- fit$est$n_parameters - rsm$est$n_parameters

  # composite-likelihood calibration: eigenvalues of the Godambe ratio over
  # the constrained directions (see Details)
  chisq_adj <- p_adj <- NA_real_; lambda <- NULL
  Bp <- fit$est$B; Br <- rsm$est$B
  if (!is.null(Bp) && !is.null(Br) && !is.null(fit$est$H_beta)) {
    M <- nrow(Bp)
    S <- cbind(Br, rep(1, M))            # rating subspace + the null shift
    Portho <- diag(M) - S %*% solve(crossprod(S), t(S))
    A <- Portho %*% Bp                   # constraint: A beta = 0
    qa <- qr(t(A))
    r <- qa$rank
    if (r > 0) {
      C <- qr.Q(qa)[, seq_len(r), drop = FALSE]
      Hinv <- solve(-fit$est$H_beta)               # Godambe H = -Hessian
      num <- crossprod(C, fit$est$cov_beta %*% C)   # C' H^-1 J H^-1 C
      den <- crossprod(C, Hinv %*% C)               # C' H^-1 C
      lambda <- Re(eigen(solve(den, num), only.values = TRUE)$values)
      chisq_adj <- chisq * r / sum(lambda)
      p_adj <- pchisq(chisq_adj, r, lower.tail = FALSE)
      df <- r
    }
  }
  out <- list(chisq = chisq, df = df,
              p = pchisq(chisq, df, lower.tail = FALSE),
              chisq_adj = chisq_adj, p_adj = p_adj, lambda = lambda,
              loglik_pcm = fit$est$loglik, loglik_rsm = rsm$est$loglik,
              fit_rsm = rsm)
  class(out) <- "rasch_lr"
  out
}

#' @export
print.rasch_lr <- function(x, ...) {
  cat("Likelihood-ratio test: partial credit vs rating parameterisation\n")
  cat(sprintf("  Raw composite chi-square %.3f on %d df, p = %s (conventional display; anticonservative)\n",
              x$chisq, x$df, .fmt_p(x$p)))
  if (is.finite(x$chisq_adj))
    cat(sprintf("  Adjusted chi-square %.3f on %d df, p = %s (Kent 1982 first-order calibration)\n",
                x$chisq_adj, x$df, .fmt_p(x$p_adj)))
  cat(sprintf("  log-likelihood (pairwise composite): PCM %.3f, RSM %.3f\n",
              x$loglik_pcm, x$loglik_rsm))
  invisible(x)
}
