# rmt :: export
# ===========================================================================
# save_outputs() writes the complete analysis to disk: every table as CSV,
# every plot as PNG (and optionally PDF), and a plain-text summary. The
# Shiny app zips the resulting folder for its "download everything" button.
# ===========================================================================

.rr_device <- function(path, fmt, width, height, dpi) {
  if (fmt == "png") png(path, width = width, height = height, units = "in", res = dpi)
  else pdf(path, width = width, height = height)
}

.rr_save_plot <- function(expr, stem, dir, formats, width, height, dpi) {
  files <- character(0)
  for (fmt in formats) {
    path <- file.path(dir, paste0(stem, ".", fmt))
    ok <- tryCatch({
      .rr_device(path, fmt, width, height, dpi)
      force(expr())
      dev.off()
      TRUE
    }, error = function(e) { try(dev.off(), silent = TRUE); FALSE })
    if (ok) files <- c(files, path)
  }
  files
}

# One plot per element, to a multi-page PDF or a ZIP of PNGs by extension.
.rr_plot_batch <- function(thunks, names, stem, file, width, height, dpi) {
  ext <- tolower(tools::file_ext(file))
  if (!grepl("^(/|[A-Za-z]:)", file)) file <- file.path(getwd(), file)
  if (ext == "pdf") {
    pdf(file, width = width, height = height, onefile = TRUE)
    on.exit(dev.off(), add = TRUE)
    for (f in thunks) tryCatch(f(), error = function(e) invisible())
  } else if (ext == "zip") {
    dir <- tempfile("rmt_plots_"); dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    paths <- character(0)
    for (j in seq_along(thunks)) {
      p <- file.path(dir, paste0(stem, "_",
                                 gsub("[^A-Za-z0-9_.-]", "_", names[j]), ".png"))
      ok <- tryCatch({
        png(p, width = width, height = height, units = "in", res = dpi)
        thunks[[j]](); dev.off(); TRUE
      }, error = function(e) { try(dev.off(), silent = TRUE); FALSE })
      if (ok) paths <- c(paths, p)
    }
    wd <- getwd(); setwd(dir); on.exit(setwd(wd), add = TRUE)
    utils::zip(file, files = basename(paths), flags = "-q9X")
  } else stop("`file` must end in .pdf or .zip")
  invisible(file)
}

#' Save a plot for every item
#'
#' Writes one plot per item -- the item characteristic curve, category
#' probability curves, threshold probability curves, or category
#' frequencies -- to a single multi-page PDF or a ZIP archive of PNGs,
#' chosen by the extension of \code{file}.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param what Which plot: \code{"icc"}, \code{"ccc"}, \code{"tpc"}, or
#'   \code{"cfreq"}.
#' @param file Output path ending in \code{.pdf} (one page per item) or
#'   \code{.zip} (one PNG per item).
#' @param items Item names or indices; all items by default.
#' @param n_groups Class intervals for observed overlays.
#' @param grid Logit grid for the curves.
#' @param observed Overlay observed proportions on the category and
#'   threshold probability curves.
#' @param width,height,dpi Device size in inches and PNG resolution.
#' @return Invisibly, the output path.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' f <- rasch(X)
#' save_item_plots(f, "icc", file.path(tempdir(), "icc_all.pdf"))
#' @export
save_item_plots <- function(fit, what = c("icc", "ccc", "tpc", "cfreq"),
                            file, items = NULL, n_groups = fit$n_groups,
                            grid = seq(-5, 5, 0.05), observed = TRUE,
                            width = 8, height = 5.5, dpi = 300) {
  what <- match.arg(what)
  its <- if (is.null(items)) fit$items$item else items
  draw <- function(it) switch(what,
    icc   = plot_icc(fit, it, n_groups = n_groups, grid = grid),
    ccc   = plot_ccc(fit, it, grid = grid, observed = observed,
                     n_groups = n_groups),
    tpc   = plot_threshold_prob(fit, it, grid = grid, observed = observed,
                                n_groups = n_groups),
    cfreq = plot_catfreq(fit, it))
  thunks <- lapply(its, function(it) function() draw(it))
  .rr_plot_batch(thunks, as.character(its), what, file, width, height, dpi)
}

#' Save a kidmap for every person
#'
#' Writes one kidmap (\code{\link{plot_kidmap}}) per person to a single
#' multi-page PDF or a ZIP archive of PNGs, chosen by the extension of
#' \code{file}. Persons without a location estimate are skipped.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param file Output path ending in \code{.pdf} (one page per person) or
#'   \code{.zip} (one PNG per person).
#' @param persons Row numbers or IDs; all estimated persons by default.
#' @param level Confidence level of the band marking unexpected responses.
#' @param width,height,dpi Device size in inches and PNG resolution.
#' @return Invisibly, the output path.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(60 * 6, 1, plogis(outer(rnorm(60), d, "-"))), 60, 6)
#' colnames(X) <- paste0("I", 1:6)
#' f <- rasch(X)
#' save_person_plots(f, file.path(tempdir(), "kidmaps.pdf"), persons = 1:5)
#' @export
save_person_plots <- function(fit, file, persons = NULL, level = 0.95,
                              width = 8, height = 6, dpi = 300) {
  ps <- if (is.null(persons)) which(!is.na(fit$person$theta)) else persons
  ids <- if (is.numeric(ps)) fit$person$id[ps] else ps
  thunks <- lapply(ps, function(p)
    function() plot_kidmap(fit, p, level = level))
  .rr_plot_batch(thunks, as.character(ids), "kidmap", file,
                 width, height, dpi)
}

#' Save every output of a Rasch analysis to a folder
#'
#' Writes all tables (item statistics, thresholds with standard errors, person
#' estimates including ID and factors, the score-to-measure table, residual
#' correlations, principal-component loadings, category frequencies, and DIF
#' results for every nominated factor) as CSV; every plot, including the
#' per-item characteristic, category, threshold, and frequency plots, as PNG
#' and optionally PDF; and a plain-text analysis summary.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param dir Output directory; created if absent.
#' @param formats Plot formats, any of \code{"png"} and \code{"pdf"}.
#' @param width,height Plot size in inches.
#' @param dpi PNG resolution; the default 300 is publication quality.
#' @param item_plots Also write the per-item plot set (one ICC, category curve,
#'   threshold curve, and frequency chart per item).
#' @return Invisibly, the vector of files written.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' out <- file.path(tempdir(), "rasch-out")
#' save_outputs(rasch(X), out, formats = "png", item_plots = FALSE)
#' @export
save_outputs <- function(fit, dir, formats = c("png", "pdf"), width = 9,
                         height = 6, dpi = 300, item_plots = TRUE) {
  formats <- match.arg(formats, c("png", "pdf"), several.ok = TRUE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tdir <- file.path(dir, "tables"); pdir <- file.path(dir, "plots")
  idir <- file.path(pdir, "items")
  dir.create(tdir, showWarnings = FALSE)
  dir.create(pdir, showWarnings = FALSE)
  if (item_plots) dir.create(idir, showWarnings = FALSE)
  files <- character(0)
  wtab <- function(d, name) {
    path <- file.path(tdir, paste0(name, ".csv"))
    utils::write.csv(d, path, row.names = FALSE)
    files <<- c(files, path)
  }

  # --- tables ---------------------------------------------------------------
  wtab(fit$items, "item_statistics")
  wtab(fit$item_anova, "item_anova_fit")
  thr <- fit$thresholds
  thr$item <- fit$items$item[thr$item]
  wtab(thr[, c("item", "k", "tau", "se")], "thresholds")
  if (!is.null(fit$est$components)) wtab(fit$est$components, "principal_components")
  wtab(fit$person, "person_estimates")
  if (!is.null(fit$score_table)) wtab(score_table(fit), "score_to_measure")
  ctt <- tryCatch(ctt_table(fit), error = function(e) NULL)
  if (!is.null(ctt)) wtab(ctt$table, "traditional_statistics")
  cd_all <- do.call(rbind, lapply(fit$items$item, function(it) {
    cd <- chisq_detail(fit, it)
    cbind(item = cd$item, cd$intervals)
  }))
  wtab(cd_all, "chisq_class_interval_detail")
  rc <- residual_correlations(fit)
  wtab(data.frame(item = rownames(rc$matrix), round(rc$matrix, 4),
                  check.names = FALSE), "residual_correlations")
  wtab(rc$pairs, "q3_statistics")
  if (nrow(rc$flagged)) wtab(rc$flagged, "local_dependence_flagged")
  pc <- residual_pca(fit)
  wtab(pc$loadings_matrix, "pca_loadings")
  wtab(pc$eigen_table, "residual_eigenvalues")
  cf <- do.call(rbind, lapply(fit$thresholds_diag, function(d)
    data.frame(item = d$item, category = seq_along(d$category_counts) - 1L,
               count = d$category_counts)))
  wtab(cf, "category_frequencies")
  gt <- guttman_table(fit)
  wtab(data.frame(id = rownames(gt$matrix), gt$matrix, check.names = FALSE),
       "guttman_ordered_responses")
  if (!is.null(fit$mc)) wtab(distractor_analysis(fit), "distractor_analysis")
  if (inherits(fit, "rasch_mfrm")) {
    wtab(fit$item_effects, "item_effects")
    wtab(fit$item_thresholds, "item_structural_thresholds")
    for (f in fit$facet_spec)
      wtab(fit$facet_effects[[f]], paste0("facet_", gsub("[^A-Za-z0-9_.-]", "_", f)))
    if (!is.null(fit$interaction_effects))
      wtab(fit$interaction_effects, "item_by_facet_interactions")
  }
  if (inherits(fit, "rasch_efrm")) {
    wtab(fit$frames, "frames")
    wtab(fit$phi_table, "group_units_phi")
    wtab(fit$alpha_table, "set_units_alpha")
    wtab(fit$set_table, "set_locations")
    wtab(fit$item_arbitrary, "items_common_unit")
    wtab(fit$thresholds_arbitrary, "thresholds_common_unit")
    wtab(fit$score_curves, "score_curves")
  }
  if (!is.null(fit$factors)) {
    dif <- tryCatch(dif_anova(fit), error = function(e) NULL)
    if (!is.null(dif)) wtab(dif, "dif_anova")
    if (ncol(fit$factors) > 1) {
      fa <- tryCatch(dif_anova_factorial(fit), error = function(e) NULL)
      if (!is.null(fa)) {
        wtab(fa$terms, "dif_factorial_terms")
        if (nrow(fa$tukey)) wtab(fa$tukey, "dif_factorial_tukey")
      }
    }
  }
  if (any(fit$person$extreme)) {
    pe <- tryCatch(person_extrapolated(fit), error = function(e) NULL)
    if (!is.null(pe)) wtab(pe, "person_estimates_extrapolated")
  }

  # --- summary ---------------------------------------------------------------
  spath <- file.path(dir, "summary.txt")
  con <- file(spath, "w")
  sink(con); on.exit({ sink(); close(con) }, add = TRUE)
  summary(fit)
  dt <- dimensionality_test(fit)
  if (is.null(dt$note)) {
    cat(sprintf("\nUnidimensionality t-test: %.1f%% significant (exact 95%% CI %.1f%% to %.1f%%), %s\n",
                100 * dt$prop_significant, 100 * dt$ci[1], 100 * dt$ci[2],
                if (dt$multidimensional) "MULTIDIMENSIONAL" else "consistent with one dimension"))
  }
  cat(sprintf("Average residual correlation: %.3f; %d flagged dependent pair(s)\n",
              rc$average, nrow(rc$flagged)))
  if (!is.null(ctt))
    cat(sprintf("Traditional statistics (complete cases n = %d): raw mean %.2f, SD %.2f, alpha %.3f, SEM %.2f\n",
                ctt$n, ctt$mean, ctt$sd, ctt$alpha, ctt$sem))
  sink(); close(con); on.exit()
  files <- c(files, spath)

  # --- test-level plots --------------------------------------------------------
  sp <- function(f, stem) files <<- c(files,
    .rr_save_plot(f, stem, pdir, formats, width, height, dpi))
  sp(function() plot_pimap(fit), "person_item_distribution")
  sp(function() plot_wright(fit), "wright_map")
  sp(function() plot_threshold_map(fit), "threshold_map")
  sp(function() plot_tcc(fit), "test_characteristic_curve")
  sp(function() plot_tif(fit), "test_information")
  sp(function() plot_item_map(fit), "item_fit_map")
  sp(function() plot_person_fit(fit), "person_fit")
  sp(function() plot_resid_cor(fit), "residual_correlations")
  sp(function() plot_pca(fit), "pca_loadings")
  sp(function() plot_scree(fit), "scree")
  sp(function() plot_guttman(fit), "guttman_scalogram")
  sp(function() plot_resid_dist(fit, "items"), "item_residual_distribution")
  sp(function() plot_resid_dist(fit, "persons"), "person_residual_distribution")
  if (inherits(fit, "rasch_mfrm")) {
    for (f in fit$facet_spec) local({
      f_ <- f
      sp(function() plot_facets(fit, f_),
         paste0("facet_severities_", gsub("[^A-Za-z0-9_.-]", "_", f_)))
    })
  }
  if (inherits(fit, "rasch_efrm")) {
    sp(function() plot_frames(fit), "frame_units")
    if (item_plots) for (it in unique(fit$virtual_map$item)) local({
      it_ <- it
      files <<- c(files, .rr_save_plot(function() plot_icc_frames(fit, it_),
        paste0(gsub("[^A-Za-z0-9_.-]", "_", it_), "_icc_frames"),
        idir, formats, width, height, dpi))
    })
  }

  # --- per-item plots ------------------------------------------------------------
  if (item_plots && !is.null(fit$mc)) {
    for (it in colnames(fit$mc$raw)) local({
      it_ <- it
      files <<- c(files, .rr_save_plot(function() plot_distractors(fit, it_),
        paste0(gsub("[^A-Za-z0-9_.-]", "_", it_), "_options"),
        idir, formats, width, height, dpi))
    })
  }
  if (item_plots) {
    for (it in fit$items$item) {
      safe <- gsub("[^A-Za-z0-9_.-]", "_", it)
      files <- c(files,
        .rr_save_plot(function() plot_icc(fit, it),
                      paste0(safe, "_icc"), idir, formats, width, height, dpi),
        .rr_save_plot(function() plot_ccc(fit, it),
                      paste0(safe, "_categories"), idir, formats, width, height, dpi),
        .rr_save_plot(function() plot_threshold_prob(fit, it),
                      paste0(safe, "_thresholds"), idir, formats, width, height, dpi),
        .rr_save_plot(function() plot_catfreq(fit, it),
                      paste0(safe, "_frequencies"), idir, formats, width, height, dpi))
    }
  }
  invisible(files)
}

# base-R base64 encoder (RFC 4648), used to embed plot images in the
# self-contained HTML report without adding dependencies
.b64 <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  alphabet <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "")[[1]]
  n <- length(raw)
  pad <- (3 - n %% 3) %% 3
  raw <- c(raw, as.raw(rep(0, pad)))
  m <- matrix(as.integer(raw), nrow = 3)
  b1 <- m[1, ] %/% 4
  b2 <- (m[1, ] %% 4) * 16 + m[2, ] %/% 16
  b3 <- (m[2, ] %% 16) * 4 + m[3, ] %/% 64
  b4 <- m[3, ] %% 64
  out <- alphabet[rbind(b1, b2, b3, b4) + 1L]
  if (pad > 0) out[(length(out) - pad + 1):length(out)] <- "="
  paste(out, collapse = "")
}

.report_css <- "
  body { font-family: -apple-system, 'Segoe UI', Roboto, Helvetica, Arial,
         sans-serif; color: #0f172a; margin: 0; background: #f8fafc; }
  .wrap { max-width: 980px; margin: 0 auto; padding: 2.5rem 1.5rem 4rem; }
  h1 { font-size: 1.6rem; margin: 0 0 .25rem; }
  h2 { font-size: 1.15rem; margin: 2.2rem 0 .6rem; padding-bottom: .3rem;
       border-bottom: 2px solid #e2e8f0; }
  .meta { color: #64748b; font-size: .85rem; margin-bottom: 1.5rem; }
  .note { color: #64748b; font-size: .82rem; margin: .3rem 0 .8rem; }
  table { border-collapse: collapse; width: 100%; font-size: .82rem;
          background: #fff; margin: .4rem 0 1rem; }
  th { text-align: left; font-weight: 600; border-bottom: 2px solid #cbd5e1;
       padding: .35rem .55rem; white-space: nowrap; }
  td { border-bottom: 1px solid #eef2f7; padding: .3rem .55rem; }
  td.num { text-align: right; font-variant-numeric: tabular-nums; }
  tr:hover td { background: #f1f5f9; }
  img { max-width: 100%; background: #fff; border: 1px solid #e2e8f0;
        border-radius: 8px; margin: .4rem 0 1rem; }
  .flag { color: #dc2626; font-weight: 600; }
  .chip { display: inline-block; background: #eff6ff; color: #1d4ed8;
          border-radius: 999px; padding: .1rem .6rem; font-size: .78rem;
          margin-right: .35rem; }
"

.html_table <- function(d, digits = 3, max_rows = 500) {
  if (is.null(d) || !nrow(d)) return("")
  trunc_note <- ""
  if (nrow(d) > max_rows) {
    trunc_note <- sprintf("<p class='note'>Showing the first %d of %d rows.</p>",
                          max_rows, nrow(d))
    d <- d[seq_len(max_rows), , drop = FALSE]
  }
  esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", as.character(x)))
  # drop all-FALSE logical flag columns and constant 'max' columns
  drop <- vapply(seq_along(d), function(j)
    (is.logical(d[[j]]) && !any(d[[j]], na.rm = TRUE)) ||
    (names(d)[j] == "max" && length(unique(d[[j]])) == 1L), TRUE)
  d <- d[, !drop, drop = FALSE]
  if (!ncol(d)) return("")
  fd <- .fmt_df(d, digits)
  num <- vapply(d, is.numeric, TRUE)
  cells <- vapply(seq_len(ncol(d)), function(j) {
    v <- fd[[j]]
    if (num[j]) gsub("<", "&lt;", v) else esc(v)
  }, character(nrow(d)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = 1)
  head_html <- paste0("<th>", esc(names(d)), "</th>", collapse = "")
  body_html <- paste(vapply(seq_len(nrow(d)), function(i) {
    paste0("<tr>", paste0("<td", ifelse(num, " class='num'", ""), ">",
                          cells[i, ], "</td>", collapse = ""), "</tr>")
  }, ""), collapse = "\n")
  paste0(trunc_note, "<table><thead><tr>", head_html,
         "</tr></thead><tbody>", body_html, "</tbody></table>")
}

#' Write a self-contained HTML report of a Rasch analysis
#'
#' Builds a single portable HTML file containing the complete analysis:
#' the summary statistics, every diagnostic table, and every test-level
#' plot embedded as an image, styled for reading and sharing. The file has
#' no external dependencies, so it can be e-mailed or archived as the
#' record of an analysis.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param file Path of the HTML file to write.
#' @param title Report title.
#' @param dpi Resolution of the embedded plots.
#' @return Invisibly, \code{file}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' out <- file.path(tempdir(), "report.html")
#' report_html(rasch(X), out)
#' @export
report_html <- function(fit, file, title = "Rasch measurement analysis",
                        dpi = 150) {
  tmp <- tempfile("rmtplots"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  shot <- function(f, name, w = 9, h = 5.4) {
    path <- file.path(tmp, paste0(name, ".png"))
    ok <- tryCatch({
      png(path, width = w, height = h, units = "in", res = dpi)
      f(); dev.off(); TRUE
    }, error = function(e) { try(dev.off(), silent = TRUE); FALSE })
    if (!ok) return("")
    sprintf("<img src='data:image/png;base64,%s' alt='%s'/>", .b64(path), name)
  }
  s <- function(...) paste0(...)
  chips <- s("<span class='chip'>", fit$model, "</span>",
             "<span class='chip'>", nrow(fit$X), " persons</span>",
             "<span class='chip'>", ncol(fit$X), " items</span>",
             sprintf("<span class='chip'>PSI %.3f</span>", fit$psi$PSI),
             sprintf("<span class='chip'>alpha %.3f</span>", fit$alpha$alpha))
  summ <- s(
    sprintf("<p>Pairwise conditional estimation %s in %d iterations. ",
            if (isTRUE(fit$est$converged)) "converged" else "did <b>not</b> converge",
            fit$est$iterations),
    sprintf("Total item-trait chi-square %.2f on %d df (p = %s). ",
            fit$total_chisq, fit$total_df, .fmt_p(fit$total_chisq_p)),
    sprintf("Item fit residual mean %.2f, SD %.2f; person fit residual mean %.2f, SD %.2f. ",
            fit$item_fit_summary$mean, fit$item_fit_summary$sd,
            fit$person_fit_summary$mean, fit$person_fit_summary$sd),
    sprintf("PSI %.3f (%.3f without extremes); item separation %.3f; power of the test of fit: %s.</p>",
            fit$psi$PSI, fit$psi_noext$PSI, fit$isi$PSI, fit$power_of_fit),
    if (length(fit$notes))
      s("<p class='note'>Notes: ", paste(fit$notes, collapse = "; "), "</p>")
    else "")
  rc <- residual_correlations(fit)
  dt <- dimensionality_test(fit)
  dim_html <- if (is.null(dt$note))
    sprintf("<p>%.1f%% of person subset t-tests significant (95%% CI %.1f-%.1f%%): %s.</p>",
            100 * dt$prop_significant, 100 * dt$ci[1], 100 * dt$ci[2],
            if (dt$multidimensional) "<span class='flag'>evidence of multidimensionality</span>"
            else "consistent with one dimension")
  else sprintf("<p class='note'>%s</p>", dt$note)
  ctt <- tryCatch(ctt_table(fit), error = function(e) NULL)

  html <- s(
    "<!DOCTYPE html><html><head><meta charset='utf-8'/>",
    "<meta name='viewport' content='width=device-width, initial-scale=1'/>",
    "<title>", title, "</title><style>", .report_css, "</style></head><body>",
    "<div class='wrap'>",
    "<h1>", title, "</h1>",
    "<p class='meta'>", format(Sys.time(), "%Y-%m-%d %H:%M"),
    " &middot; rmt ", as.character(utils::packageVersion("rmt")), "</p>",
    "<p>", chips, "</p>",
    "<h2>Summary</h2>", summ,
    "<h2>Targeting</h2>",
    shot(function() plot_pimap(fit), "targeting"),
    shot(function() plot_wright(fit), "wright_map"),
    "<h2>Item statistics</h2>",
    .html_table(fit$items[, intersect(c("item", "max", "location", "se",
                                        "fit_resid", "infit_ms", "outfit_ms",
                                        "chisq", "df", "p_adj"),
                                      names(fit$items))]),
    shot(function() plot_item_map(fit), "item_map"),
    "<h2>Thresholds</h2>",
    .html_table({ th <- fit$thresholds
                  th$item <- fit$items$item[th$item]
                  th[, c("item", "k", "tau", "se")] }),
    { dis <- names(which(vapply(fit$thresholds_diag, function(dd)
        !dd$ordered && length(dd$thresholds) > 1L, TRUE)))
      if (length(dis)) sprintf("<p class='flag'>Disordered thresholds: %s.</p>",
                               paste(dis, collapse = ", "))
      else "<p class='note'>All polytomous items have ordered thresholds.</p>" },
    shot(function() plot_threshold_map(fit), "threshold_map"),
    "<h2>Test characteristic and information</h2>",
    shot(function() plot_tcc(fit), "tcc"),
    shot(function() plot_tif(fit), "tif"),
    "<h2>Score to measure</h2>",
    .html_table(score_table(fit)),
    "<h2>Fit residual distributions</h2>",
    shot(function() plot_resid_dist(fit, "items"), "resid_items"),
    shot(function() plot_resid_dist(fit, "persons"), "resid_persons"),
    "<h2>Dimensionality</h2>", dim_html,
    shot(function() plot_scree(fit), "scree"),
    "<h2>Local dependence</h2>",
    sprintf("<p>Average residual correlation %.3f; %d flagged pair(s).</p>",
            rc$average, nrow(rc$flagged)),
    if (nrow(rc$flagged)) .html_table(rc$flagged) else "",
    shot(function() plot_resid_cor(fit), "residcor"),
    if (!is.null(ctt)) s("<h2>Classical companions</h2>",
      sprintf("<p class='note'>Complete cases n = %d; raw mean %.2f, SD %.2f; alpha %.3f; classical SEM %.2f.</p>",
              ctt$n, ctt$mean, ctt$sd, ctt$alpha, ctt$sem),
      .html_table(ctt$table)) else "",
    if (!is.null(fit$factors)) {
      da <- tryCatch(dif_anova(fit), error = function(e) NULL)
      if (!is.null(da)) s("<h2>Differential item functioning</h2>",
        .html_table(da[, intersect(c("factor", "item", "F_uniform",
                                     "p_uniform_adj", "eta2_uniform",
                                     "F_nonuniform", "p_nonuniform_adj",
                                     "eta2_nonuniform", "uniform_DIF",
                                     "nonuniform_DIF"), names(da))])) else ""
    } else "",
    if (!is.null(fit$mc)) s("<h2>Distractor analysis</h2>",
      "<p class='note'>Locations use the rest measure; a distractor whose takers are abler than the keyed option's flags a possible miskey.</p>",
      .html_table(tryCatch(distractor_analysis(fit), error = function(e) NULL))) else "",
    if (inherits(fit, "rasch_mfrm")) s("<h2>Facet severities</h2>",
      paste(vapply(fit$facet_spec, function(f) s("<h3>", f, "</h3>",
        .html_table(fit$facet_effects[[f]][, intersect(c("level", "severity",
          "se", "n", "fit_resid"), names(fit$facet_effects[[f]]))])), ""),
        collapse = "")) else "",
    if (inherits(fit, "rasch_efrm")) s("<h2>Frames and units</h2>",
      .html_table(fit$frames[, intersect(c("set", "group", "rho", "se_log_rho",
        "origin", "fit_resid", "n_responses"), names(fit$frames))])) else "",
    "<h2>Person estimates</h2>",
    .html_table(fit$person[, intersect(c("id", names(fit$factors), "raw",
                                         "max_raw", "theta", "se", "extreme",
                                         "fit_resid"),
                                       names(fit$person))]),
    "</div></body></html>")
  writeLines(html, file, useBytes = TRUE)
  invisible(file)
}
