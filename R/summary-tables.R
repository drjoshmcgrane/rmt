# rmt :: summary tables
# ===========================================================================
# The test-of-fit and targeting/reliability summaries as tidy two-column
# tables, so the headline statistics of an analysis can be saved and
# reported rather than read off a text panel.
# ===========================================================================

#' Test-of-fit summary as a table
#'
#' The headline fit statistics of a calibration -- model, estimation,
#' total item-trait chi-square, the item and person fit-residual moments,
#' fit-location correlations, chi-square flag count, and disordered
#' thresholds -- as a two-column table suitable for saving and reporting.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return A data frame with columns \code{statistic} and \code{value}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit_summary_table(rasch(X))
#' @export
fit_summary_table <- function(fit) {
  ss <- fit$summary_stats
  dis <- names(which(vapply(fit$thresholds_diag, function(d)
    !d$ordered && length(d$thresholds) > 1, TRUE)))
  num <- function(x, d = 3) formatC(x, digits = d, format = "f")
  out <- data.frame(statistic = c(
    "Model", "Estimation", "Converged", "Iterations",
    "Total item-trait chi-square", "Degrees of freedom",
    "Item-trait probability", "Class intervals",
    "Item fit residual mean", "Item fit residual SD",
    "Item fit residual skewness", "Item fit residual kurtosis",
    "Person fit residual mean", "Person fit residual SD",
    "Person fit residual skewness", "Person fit residual kurtosis",
    "Fit-location correlation (items)",
    "Fit-location correlation (persons)",
    "Items with adjusted chi-square p < .05",
    "Disordered thresholds"),
    value = c(
      fit$model, "pairwise conditional ML",
      ifelse(isTRUE(fit$est$converged), "yes", "NO"),
      as.character(fit$est$iterations),
      num(fit$total_chisq), as.character(fit$total_df),
      .fmt_p(fit$total_chisq_p), as.character(fit$n_groups),
      num(fit$item_fit_summary$mean, 2), num(fit$item_fit_summary$sd, 2),
      num(fit$item_fit_summary$skewness, 2),
      num(fit$item_fit_summary$kurtosis, 2),
      num(fit$person_fit_summary$mean, 2), num(fit$person_fit_summary$sd, 2),
      num(fit$person_fit_summary$skewness, 2),
      num(fit$person_fit_summary$kurtosis, 2),
      num(ss$cor_item_fit_location), num(ss$cor_person_fit_location),
      sprintf("%d of %d", sum(fit$items$p_adj < 0.05, na.rm = TRUE),
              nrow(fit$items)),
      if (length(dis)) paste(dis, collapse = ", ") else "none"))
  rownames(out) <- NULL
  out
}

#' Targeting and reliability summary as a table
#'
#' The person and item location moments (with and without extreme persons),
#' the threshold range and coverage, and the reliability indices -- PSI,
#' PSI without extremes, item separation, and coefficient alpha -- as a
#' two-column table suitable for saving and reporting.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return A data frame with columns \code{statistic} and \code{value}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' targeting_table(rasch(X))
#' @export
targeting_table <- function(fit) {
  ss <- fit$summary_stats; t <- fit$targeting
  num <- function(x, d = 3) ifelse(is.finite(x),
                                   formatC(x, digits = d, format = "f"), "NA")
  out <- data.frame(statistic = c(
    "Person location mean", "Person location SD",
    "Person location skewness", "Person location kurtosis",
    "Person location mean (no extremes)", "Person location SD (no extremes)",
    "Item location mean (constrained)", "Item location SD",
    "Item location skewness", "Item location kurtosis",
    "Threshold minimum", "Threshold maximum",
    "Persons below threshold range (%)", "Persons above threshold range (%)",
    "PSI", "Separation", "PSI without extremes",
    "n without extremes", "Item separation index", "Coefficient alpha",
    "n complete (alpha)"),
    value = c(
      num(ss$person_location$mean), num(ss$person_location$sd),
      num(ss$person_location$skewness, 2), num(ss$person_location$kurtosis, 2),
      num(ss$person_location_noext$mean), num(ss$person_location_noext$sd),
      num(ss$item_location$mean), num(ss$item_location$sd),
      num(ss$item_location$skewness, 2), num(ss$item_location$kurtosis, 2),
      num(t$threshold_range[1]), num(t$threshold_range[2]),
      num(100 * t$prop_below, 1), num(100 * t$prop_above, 1),
      num(fit$psi$PSI), num(fit$psi$separation),
      num(fit$psi_noext$PSI), as.character(fit$psi_noext$n),
      num(fit$isi$PSI),
      num(fit$alpha$alpha),
      if (is.null(fit$alpha$n)) "NA" else as.character(fit$alpha$n)))
  rownames(out) <- NULL
  out
}
