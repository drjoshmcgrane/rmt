# rasch Shiny GUI
# ---------------------------------------------------------------------------
# A modern bslib interface to the full rasch analysis: data upload with ID,
# person-factor, and item column nomination; pairwise conditional ML
# estimation (Andrich & Luo 2003); the complete test-of-fit suite;
# every diagnostic plot with per-plot PNG and PDF downloads; and one-click
# export of all tables and plots as a ZIP archive.
# Launch with rasch::run_app(), or shiny::runApp() from this folder.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(bsicons)
})

if (requireNamespace("rasch", quietly = TRUE)) {
  library(rasch)
} else {
  rdir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  if (dir.exists(rdir)) {
    for (f in list.files(rdir, "\\.R$", full.names = TRUE)) source(f)
  } else stop("Install rasch, or run the app from inst/shiny in the source tree")
}

# --- demo data: 10 polytomous items, one disordered, DIF on Q05 -------------
.demo_data <- function(seed = 11, Np = 1200) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  mvec <- rep(c(2, 3), length.out = 10)
  tau_true <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.9)))
  tau_true[[2]] <- c(1.2, -1.3, 0.6)                       # disordered item
  th <- rnorm(Np, 0, 1.4)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  sex <- sample(c("female", "male"), Np, replace = TRUE)
  X <- sapply(seq_along(mvec), function(i) {
    sft <- if (i == 5) ifelse(grp == "focal", 0.9, 0) else numeric(Np)  # uniform DIF
    sapply(seq_len(Np), function(n) sample(0:mvec[i], 1, prob = simP(th[n] - sft[n], tau_true[[i]])))
  })
  colnames(X) <- sprintf("Q%02d", seq_along(mvec))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, sex = sex, check.names = FALSE)
}

# dichotomous demo: 15 multiple-choice items (raw A-D responses), DIF planted
# on I05 by group, and I07 deliberately miskeyed (true correct C, key says A)
.demo_dich <- function(seed = 41, Np = 1000) {
  set.seed(seed)
  d <- seq(-2, 2, length.out = 15)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  sex <- sample(c("female", "male"), Np, replace = TRUE)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_along(d), function(i) {
    sft <- if (i == 5) ifelse(grp == "focal", 0.8, 0) else 0
    correct <- if (i == 7) "C" else "A"
    ok <- rbinom(Np, 1, plogis(th - d[i] - sft))
    ifelse(ok == 1, correct,
           sample(setdiff(c("A", "B", "C", "D"), correct), Np, replace = TRUE))
  })
  colnames(X) <- sprintf("I%02d", seq_along(d))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, sex = sex, check.names = FALSE)
}

# the demo key: all "A" (so I07 is the discoverable miskey)
.demo_dich_key <- function()
  setNames(rep("A", 15), sprintf("I%02d", 1:15))

# paired-comparison demo: 8 essays compared pairwise by 10 judges, with
# judge J09 answering at random (discoverable in the judge fit table).
# Besides the winner column it carries a graded `preference` column (four
# ordered categories) simulated from the same object locations, so the
# graded-response role can be pointed at it, and a `margin` column ("a
# little" < "much") derived from the preference for the winner + margin
# entry path.
.demo_btl <- function(seed = 47, reps = 26) {
  set.seed(seed)
  # a moderate object spread: extreme objects have near-zero residual variance
  # in paired comparisons, which would manufacture spurious DIF, so the range
  # is kept modest and the planted DIF is put on the central objects
  beta <- setNames(seq(-1.0, 1.0, length.out = 8), sprintf("E%02d", 1:8))
  pr <- t(utils::combn(names(beta), 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps),
                  object_b = rep(pr[, 2], each = reps),
                  stringsAsFactors = FALSE)
  d$judge <- sprintf("J%02d", sample(1:10, nrow(d), replace = TRUE))
  # two judge factors (each constant within judge) so DIF can be modelled by
  # one factor or several jointly: panel splits judges 1-5 / 6-10, experience
  # splits the odd / even judges independently
  d$panel <- ifelse(d$judge %in% sprintf("J%02d", 1:5), "panel A", "panel B")
  d$experience <- ifelse(d$judge %in% sprintf("J%02d", c(1, 3, 5, 7, 9)),
                         "expert", "novice")
  # judgment order: process each judge's comparisons in sequence so the
  # within-judge history (exposure) is well defined
  d <- d[sample(nrow(d)), ]
  d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  d <- d[order(d$judge, d$t), ]
  rownames(d) <- NULL

  # Several signals are built in so the diagnostics have something to find.
  # (1) Exposure: an object already met by the judge gains `expo` logits (a
  # seen-before advantage). (2) Panel DIF: panel A over-rewards E04, so its
  # location differs by panel. (3) Experience DIF: experts over-reward E05, a
  # second, independent judge factor. The two factors point at different
  # objects so each is cleanly attributable. Judge J09 answers at random (a
  # misfitting judge), as before.
  expo <- 0.7
  dif_panel <- "E04"; dif_exp_obj <- "E05"
  dif <- 1.2         # panel effect
  # larger than the panel effect: expert J09 answers at random (diluting it),
  # and the dependence-adjusted DIF screen absorbs the share of a sequential
  # effect that the carry-over covariate can carry
  dif_exp <- 1.8
  tau <- c(-1.1, 0, 1.1)
  lev <- c("much worse", "a little worse", "a little better", "much better")
  seen <- new.env(parent = emptyenv())
  winner <- character(nrow(d)); pref <- integer(nrow(d))
  for (r in seq_len(nrow(d))) {
    j <- d$judge[r]; a <- d$object_a[r]; b <- d$object_b[r]
    ba <- beta[[a]]; bb <- beta[[b]]
    if (d$panel[r] == "panel A") {
      ba <- ba + dif * (a == dif_panel); bb <- bb + dif * (b == dif_panel)
    }
    if (d$experience[r] == "expert") {
      ba <- ba + dif_exp * (a == dif_exp_obj); bb <- bb + dif_exp * (b == dif_exp_obj)
    }
    ba <- ba + expo * isTRUE(get0(paste(j, a), seen, ifnotfound = FALSE))
    bb <- bb + expo * isTRUE(get0(paste(j, b), seen, ifnotfound = FALSE))
    if (j == "J09") { p <- 0.5; Pp <- rep(0.25, 4) }
    else { p <- plogis(ba - bb); Pp <- item_moments(ba - bb, tau)$P }
    winner[r] <- if (runif(1) < p) a else b
    pref[r] <- sample.int(4, 1, prob = Pp)
    assign(paste(j, a), TRUE, seen); assign(paste(j, b), TRUE, seen)
  }
  d$winner <- winner
  d$preference <- factor(lev[pref], levels = lev, ordered = TRUE)
  # margin of win as an ordered factor (extreme categories are "much" wins)
  d$margin <- factor(ifelse(d$preference %in% c("much worse", "much better"),
                            "much", "a little"),
                     levels = c("a little", "much"), ordered = TRUE)
  d <- d[sample(nrow(d)), ]   # present in random row order (t keeps the order)
  rownames(d) <- NULL
  d
}

# rating scale demo: common step structure, item locations vary
.demo_rsm <- function(seed = 51, Np = 1000) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  loc <- seq(-1.2, 1.2, length.out = 8)
  step <- c(-0.9, 0.0, 0.9)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(loc, function(b) sapply(th, function(t)
    sample(0:3, 1, prob = simP(t, b + step))))
  colnames(X) <- sprintf("R%02d", seq_along(loc))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, check.names = FALSE)
}

# rated (MFRM) demo, wide layout: 5 item columns, 6 raters (one erratic),
# incomplete design — one row per person-by-rater combination. The responses
# are simulated in long form (same structure and seed as always) and
# reshaped, so results are unchanged.
.demo_mfrm <- function(seed = 21, Np = 250) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%04d", seq_len(Np)); raters <- paste0("Rater_", 1:6)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  rho <- setNames(c(-0.9, -0.4, -0.1, 0.1, 0.4, 0.9), raters)
  tau <- list(Essay = c(-1.2, 0.2, 1.1), Argument = c(-0.8, 0.5, 1.3),
              Evidence = c(-1.5, -0.2, 0.9), Style = c(-0.6, 0.4, 1.2),
              Mechanics = c(-1.0, 0.0, 1.0))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  seen <- unlist(lapply(persons, function(p) paste(p, sample(raters, 3))))
  d <- d[paste(d$person, d$rater) %in% seen, ]
  d$score <- mapply(function(p, i, r) {
    if (r == "Rater_6" && runif(1) < 0.2) return(sample(0:3, 1))  # erratic rater
    sample(0:3, 1, prob = simP(th[p], tau[[i]] + rho[r]))
  }, d$person, d$item, d$rater)
  # wide: one row per person-by-rater with one column per item
  w <- reshape(d, idvar = c("person", "rater"), timevar = "item",
               v.names = "score", direction = "wide")
  names(w) <- sub("^score\\.", "", names(w))
  w <- w[order(w$person, w$rater), c("person", "rater", names(tau))]
  rownames(w) <- NULL
  w
}

# frames demo: 2 person groups x 3 item sets with distinct units
.demo_efrm <- function(seed = 31, per_g = 350) {
  set.seed(seed)
  simP <- function(th, tau, r) { x <- 0:length(tau); p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
  glev <- c("year5", "year7"); grp <- rep(glev, each = per_g); Np <- length(grp)
  phi <- c(year5 = 0.8, year7 = 1.25)
  sets <- rep(c("Number", "Algebra", "Space"), each = 6)
  alpha <- c(Number = 0.75, Algebra = 1.0, Space = 4 / 3)
  th <- rnorm(Np, 0, 1.3) + ifelse(grp == "year7", 0.5, 0)
  d <- as.numeric(sapply(c(-0.3, 0.1, 0.2), function(m) m + seq(-1.2, 1.2, length.out = 6)))
  X <- sapply(seq_along(sets), function(i) sapply(seq_len(Np), function(n)
    sample(0:2, 1, prob = simP(th[n], d[i] + c(-0.5, 0.5),
                               alpha[sets[i]] * phi[grp[n]]))))
  colnames(X) <- sprintf("%s_%02d", sets, seq_along(sets))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X, year_group = grp,
             check.names = FALSE)
}

NONE <- "(none)"
# the sentinel VALUE stays "(none)" (the server compares against it), but it
# is always displayed as "None"; selects with no meaningful pre-fit choice
# use empty choices plus a selectize placeholder instead of a sentinel row
NONE_CH <- c(None = "(none)")

# null-coalescing helper: defined here so the app does not depend on the
# base R version that introduced it (R >= 4.4)
`%||%` <- function(a, b) if (is.null(a)) b else a

# sanitise a proportion-type input (alpha level, chance probability): fall
# back to the default outside the open interval (0, 1)
clamp01 <- function(x, default)
  if (is.null(x) || is.na(x) || x <= 0 || x >= 1) default else x

# p-values as text: "%.3f" alone prints a misleading 0.000 for tiny p
fmt_p <- function(p)
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))

# value-box guard: NULL, NA, NaN, and Inf must never reach a conditional or
# a sprintf; such values display as an em dash on a neutral theme
finite1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

# p-values phrased for prose ("p = 0.018" / "p < 0.001"), matching fmt_p's
# thresholds; used by the curated stat rows
p_lab <- function(p) {
  if (!finite1(p)) "p = NA"
  else if (p < 0.001) "p < 0.001"
  else sprintf("p = %.3f", p)
}

# curated stat-box rows (Summary page): one label-value line per statistic
stat_row <- function(label, value)
  div(class = "stat-row",
      span(class = "stat-label", label),
      span(class = "stat-value", value))
stat_rows <- function(...) div(class = "stat-rows", ...)

# measurement-themed value-box glyphs: inline SVG stroked with currentColor,
# so each glyph inherits its box's text colour in light and dark themes
.glyph_body <- list(
  # bell curve: the person distribution
  distribution = '<path d="M2 20 C7 20 8.5 5.5 12 5.5 S17 20 22 20"/>',
  # logit scale with alternating major/minor ticks
  ruler = '<path d="M2 15 H22"/><path d="M5 15 V9"/><path d="M9.5 15 V11.5"/><path d="M14 15 V9"/><path d="M18.5 15 V11.5"/>',
  # two distinct person distributions: separation / reliability
  separation = '<path d="M1 20 C4.5 20 5 10 7.5 10 S10.5 20 14 20"/><path d="M10 20 C13.5 20 14 5.5 16.5 5.5 S19.5 20 23 20"/>',
  alpha = '<text x="12" y="17.5" text-anchor="middle" font-size="17" font-style="italic" font-family="Georgia, serif" fill="currentColor" stroke="none">&#945;</text>',
  chisq = '<text x="11" y="17" text-anchor="middle" font-size="14" font-style="italic" font-family="Georgia, serif" fill="currentColor" stroke="none">&#967;&#178;</text>',
  # magnifier over a curve: power of the test of fit
  power = '<circle cx="10" cy="10" r="6.2"/><path d="M14.6 14.6 L20.5 20.5"/><path d="M6.8 11.5 Q10 6 13.2 11.5"/>',
  # data matrix rows / columns
  grid = '<rect x="3" y="4" width="18" height="16" rx="1.5"/><path d="M3 9.3 H21 M3 14.6 H21"/>',
  columns = '<rect x="3" y="4" width="18" height="16" rx="1.5"/><path d="M9 4 V20 M15 4 V20"/>',
  # 2x2 grid with a crossed-out cell: missing responses
  missing = '<rect x="3" y="4" width="18" height="16" rx="1.5"/><path d="M3 12 H21 M12 4 V20"/><path d="M14.5 14.5 L18.5 18.5 M18.5 14.5 L14.5 18.5" stroke-width="1.3"/>',
  # location span between two ends of the scale
  range = '<path d="M4 6 V18 M20 6 V18 M4 12 H20"/>',
  # a point off the trend: misfit
  outlier = '<path d="M3 17 C9 15.5 15 14.5 21 12.5"/><circle cx="16.5" cy="5.5" r="1.7" fill="currentColor" stroke="none"/>',
  # crossing solid/dashed step lines: reversed thresholds
  disorder = '<path d="M4 7 H10 L16 17 H20"/><path d="M4 17 H10 L16 7 H20" stroke-dasharray="2.6 2.2"/>',
  # two objects with a double-headed arrow: a paired comparison
  pair = '<circle cx="5.2" cy="12" r="2.7"/><circle cx="18.8" cy="12" r="2.7"/><path d="M9.3 12 H14.7 M11.2 10 L9.2 12 L11.2 14 M12.8 10 L14.8 12 L12.8 14"/>',
  # judge's balance
  balance = '<path d="M12 5 V19 M6 5 H18"/><path d="M6 5 L3.5 11 H8.5 Z"/><path d="M18 5 L15.5 11 H20.5 Z"/><path d="M9 19 H15"/>',
  # object locations as a podium
  podium = '<rect x="3.5" y="11" width="5" height="9" rx="1"/><rect x="9.5" y="6" width="5" height="14" rx="1"/><rect x="15.5" y="14" width="5" height="6" rx="1"/>')
glyph <- function(name)
  HTML(sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" style="width:100%%;height:100%%">%s</svg>',
               .glyph_body[[name]]))

# helpers for the "R code for this analysis" disclosure
qstr <- function(x) paste0('"', x, '"')
qvec <- function(x)
  if (length(x) == 1L) qstr(x) else
    paste0("c(", paste(qstr(x), collapse = ", "), ")")

theme <- bs_theme(
  version = 5, preset = "shiny",
  bg = "#f8fafc", fg = "#0f172a",
  primary = "#2563eb", secondary = "#64748b",
  success = "#0f766e", danger = "#dc2626", warning = "#d97706",
  base_font = font_google("Inter"),
  heading_font = font_google("Inter"),
  code_font = font_google("JetBrains Mono"),
  "navbar-bg" = "#0f172a",
  "headings-font-weight" = "600",
  "font-size-base" = "0.925rem"
)

css <- HTML("
  .navbar-brand { font-weight: 700; letter-spacing: .02em; }
  /* two-line wordmark: the package name over the spelled-out tagline */
  .app-brand { display: inline-flex; flex-direction: column;
               justify-content: center; line-height: 1.06; }
  .app-brand-name { font-weight: 700; font-size: 1.05rem; letter-spacing: .02em; }
  .app-brand-sub { font-weight: 400; font-size: .70rem; opacity: .6;
                   letter-spacing: .03em; }
  .card-header { font-weight: 600; }
  .value-box-title { font-size: .72rem; text-transform: uppercase; letter-spacing: .04em; white-space: nowrap; }
  .value-box-value { font-size: 1.45rem; }
  pre, .shiny-text-output { white-space: pre-wrap; font-size: .82rem; }
  .btn-xs { padding: .1rem .5rem; font-size: .75rem; }
  .form-label { font-weight: 600; font-size: .85rem; }
  /* APA-style tables: tabular numerals, no vertical rules, no zebra,
     a strong rule under the header row */
  table.dataTable { font-size: .85rem; width: auto !important; max-width: 100%; }
  table.dataTable td { font-variant-numeric: tabular-nums; }
  table.dataTable thead th {
    font-weight: 600;
    border-bottom: 2px solid var(--bs-emphasis-color) !important;
  }
  table.dataTable th, table.dataTable td {
    border-left: none !important; border-right: none !important;
  }
  table.dataTable.table-striped > tbody > tr:nth-of-type(odd) > * {
    --bs-table-accent-bg: transparent;
  }
  .table-note { color: var(--bs-secondary-color); font-size: .8rem; }
  /* curated stat boxes (Summary page): label-value rows with a hairline
     separator; tabular numerals keep the values aligned in both themes */
  .stat-rows { display: flex; flex-direction: column; }
  .stat-row { display: flex; justify-content: space-between; align-items: baseline; gap: 1rem; padding: .45rem 0; border-bottom: 1px solid var(--bs-border-color-translucent); }
  .stat-row:last-child { border-bottom: 0; }
  .stat-label { color: var(--bs-secondary-color); font-size: .85rem; }
  .stat-value { font-weight: 600; font-variant-numeric: tabular-nums; text-align: right; }
  .stat-head { color: var(--bs-secondary-color); font-size: .85rem; margin-bottom: .35rem; }
  /* card headers: title left, action chips right, with a small gap
     between chips (the chips row right-aligns even when there is no
     title, via margin-left:auto) */
  .rasch-card-header { display: flex; align-items: center; justify-content: space-between; gap: .5rem; }
  .rasch-chips { display: flex; align-items: center; flex-wrap: wrap; gap: .35rem; margin-left: auto; }
  /* inline form controls that sit on a flex row with buttons: strip the
     bottom margin Shiny's containers carry */
  .rasch-inline-check .form-group, .rasch-inline-check .shiny-input-container,
  .rasch-inline-check .checkbox,
  .rasch-inline-select .form-group, .rasch-inline-select .shiny-input-container {
    margin-bottom: 0; width: auto;
  }
  .rasch-inline-select select.form-select { padding: .15rem 1.6rem .15rem .5rem; font-size: .78rem; }
  /* collapsed advanced-settings disclosure inside the sidebar accordion */
  .rasch-advanced { margin-top: .5rem; }
  .rasch-advanced summary { cursor: pointer; font-size: .8rem; font-weight: 600; color: var(--bs-secondary-color); margin-bottom: .35rem; }
  .empty-state { text-align: center; padding: 3rem 1rem; color: var(--bs-secondary-color); }
  .shiny-output-error-validation {
    text-align: center; padding: 2rem 1rem; color: var(--bs-secondary-color);
  }
  .nav-status .badge { font-weight: 500; }
  /* the DT bottom elements (info + pager) float; without clearance they
     collide with the collapsed R-code footer and can render outside the
     card. Clear the footer, give the wrapper self-clearing bottom room. */
  div.dataTables_wrapper { padding-bottom: .25rem; }
  div.dataTables_wrapper::after { content: ''; display: block; clear: both; }
  .rcode { clear: both; margin-top: .5rem; }
  .rcode summary { cursor: pointer; font-size: .78rem; color: var(--bs-secondary-color); }
  .rcode pre { font-size: .78rem; background: var(--bs-tertiary-bg); border: 1px solid var(--bs-border-color); border-radius: 6px; padding: .5rem .75rem; margin: .35rem 0 0; white-space: pre-wrap; }
  .rcode-copy { float: right; margin-top: .35rem; }
")

# collapsed per-output "R code" footer (jamovi-style syntax mode): shows the
# exact rasch call reproducing the output, updating with the current selections;
# the server registers a matching `<id>_code` renderText for every output
rcode_details <- function(id)
  tags$details(class = "rcode",
    tags$summary(bs_icon("code-slash"), " R code"),
    tags$button(class = "btn btn-outline-secondary btn-xs rcode-copy",
                type = "button",
                onclick = sprintf(
                  "navigator.clipboard.writeText(document.getElementById('%s_code').innerText)",
                  id),
                "Copy"),
    verbatimTextOutput(paste0(id, "_code"), placeholder = FALSE))

# header info-circle tooltip used across the card helpers
info_icon <- function(info)
  tooltip(bs_icon("info-circle", class = "ms-1 text-secondary"), info)

# card header as a full-width flex bar: title (when given) on the left,
# action chips right-aligned with a small gap between them. Cards that sit
# inside an accordion panel that already names them pass title = NULL: the
# header then renders buttons-only (the info icon stays in the chips row).
card_header_bar <- function(title = NULL, buttons = NULL, info = NULL)
  card_header(class = "rasch-card-header",
    if (!is.null(title)) span(title, if (!is.null(info)) info_icon(info)),
    div(class = "rasch-chips",
        if (is.null(title) && !is.null(info)) info_icon(info),
        buttons))

# Card with a plot and PNG/PDF download buttons in the header. The body is
# non-fillable so flex sizing can never compress the fixed-height plot
# (the cause of the squashed plots), and a percentage height is avoided
# because it races the layout and renders a zero-height device.
# data-bs-theme is pinned to light because base plots draw on white.
# `info` adds a header tooltip; `controls` takes small inputs rendered
# before the download chips; `extra` takes further header buttons (e.g.
# the batch all-persons downloads on the kidmap card). title = NULL
# renders a buttons-only header (for cards inside named accordion panels).
plotCard <- function(id, title = NULL, height = "560px", info = NULL,
                     controls = NULL, extra = NULL) {
  card(
    full_screen = TRUE,
    `data-bs-theme` = "light",
    card_header_bar(title, info = info, buttons = tagList(
      controls,
      downloadButton(paste0(id, "_png"), "PNG", class = "btn-outline-secondary btn-xs"),
      downloadButton(paste0(id, "_pdf"), "PDF", class = "btn-outline-secondary btn-xs"),
      extra)),
    card_body(plotOutput(id, height = height), rcode_details(id),
              padding = 8, fillable = FALSE)
  )
}

# `info` adds a header tooltip defining the key statistic; `footer` takes a
# small UI slot rendered under the table (dynamic interpretation notes);
# title = NULL renders a buttons-only header
tableCard <- function(id, title = NULL, note = NULL, info = NULL,
                      footer = NULL, controls = NULL) {
  card(
    full_screen = TRUE,
    card_header_bar(title, info = info, buttons = tagList(
      controls,
      downloadButton(paste0(id, "_csv"), "CSV",
                     class = "btn-outline-secondary btn-xs"))),
    # non-fillable body: DT outputs are fill items, and flex sizing inside a
    # natural-height card crops the last rows and the info/pager strip
    card_body(if (!is.null(note)) p(class = "text-muted small mb-2", note),
              DTOutput(id),
              if (!is.null(footer)) div(class = "table-note mt-2", footer),
              rcode_details(id),
              padding = 12, fillable = FALSE)
  )
}


# curated stat-box card: the body is a uiOutput of label-value rows built by
# the server; the CSV chip downloads the COMPLETE summary table (never the
# curated display), and the code footer names the call that builds it
statCard <- function(id, title = NULL, info = NULL, footer = NULL) {
  card(
    full_screen = TRUE,
    card_header_bar(title, info = info,
      buttons = downloadButton(paste0(id, "_csv"), "CSV",
                               class = "btn-outline-secondary btn-xs")),
    card_body(uiOutput(id),
              if (!is.null(footer)) div(class = "table-note mt-2", footer),
              rcode_details(id),
              padding = 12, fillable = FALSE)
  )
}

# compact header switch revealing every column of a curated table
cols_switch <- function(id)
  div(class = "small text-secondary",
      input_switch(id, "All columns", value = FALSE))

# card header with an info-circle tooltip (for non-table cards)
info_header <- function(title, info)
  card_header(span(title,
    tooltip(bs_icon("info-circle", class = "ms-1 text-secondary"), info)))

# ---------------------------------------------------------------------------
# Panels are built as objects and assembled into the workflow-ordered navbar
# (with Independence / Invariance / More menus) at the end of the UI section.
# ----------------------------------------------------------------- DATA --
panel_data <- nav_panel("Data", value = "p_data", icon = bs_icon("database"),
    layout_sidebar(
      sidebar = sidebar(width = 330,
        h6("Data source"),
        fileInput("file", NULL, accept = c(".csv", ".txt", ".tsv"),
                  buttonLabel = "Browse…", placeholder = "CSV / TSV file"),
        selectInput("demo_choice", "Example dataset",
                    c("None" = "none",
                      "Multiple choice, dichotomous" = "dich",
                      "Polytomous (PCM)" = "pcm",
                      "Rating scale (RSM)" = "rsm",
                      "Ratings by raters (MFRM)" = "mfrm",
                      "Item sets x groups (EFRM)" = "efrm",
                      "Paired comparisons (BTL)" = "btl")),
        accordion(
          id = "run_settings", multiple = TRUE,
          open = c("Data roles", "Model"),
          accordion_panel("Model", icon = bs_icon("diagram-2"),
            radioButtons("model_type", NULL,
                         c("Rasch" = "rasch",
                           "Many-facet (MFRM)" = "mfrm",
                           "Extended frames (EFRM)" = "efrm",
                           "Paired comparisons (BTL)" = "btl"))),
          accordion_panel("Data roles", icon = bs_icon("table"),
            conditionalPanel("input.model_type == 'rasch'",
              selectInput("id_col", "ID variable", NONE_CH),
              selectizeInput("factor_cols", "Person factors (DIF groups)", NULL,
                             multiple = TRUE,
                             options = list(placeholder = "none")),
              selectizeInput("item_cols", "Item columns", NULL, multiple = TRUE,
                             options = list(placeholder = "all remaining"))
            ),
            conditionalPanel("input.model_type == 'efrm'",
              selectInput("ef_id", "ID variable", NONE_CH),
              selectInput("ef_group",
                          span("Person group column",
                               info_icon("None treats all persons as one group, so the units differ by item set only.")),
                          NONE_CH),
              selectizeInput("ef_items", "Item columns", NULL, multiple = TRUE,
                             options = list(placeholder = "all remaining")),
              fileInput("ef_sets", "Item-set map (CSV: item,set)",
                        accept = ".csv", placeholder = "optional"),
              checkboxInput("ef_prefix", "Infer sets from item-name prefix", TRUE),
              p(class = "text-muted small",
                "Each item-set by group cell is a frame with its own unit. Group units come from the person-free pairwise comparisons; set units from persons common to the sets.")
            ),
            conditionalPanel("input.model_type == 'mfrm'",
              radioButtons("lp_layout", "Data layout",
                           c("Items in columns (wide)" = "wide",
                             "One score per row (long)" = "long")),
              p(class = "text-muted small",
                "Wide: one row per person-by-facet combination (e.g. one row per script per rater), one column per item or criterion. Long: person, item, and score columns."),
              selectInput("lp_person", "Person column", NONE_CH),
              conditionalPanel("input.lp_layout == 'long'",
                selectInput("lp_item", "Item column", NONE_CH),
                selectInput("lp_score", "Score column", NONE_CH)),
              selectizeInput("lp_facets", "Facet columns (e.g. rater)", NULL,
                             multiple = TRUE,
                             options = list(placeholder = "choose at least one")),
              conditionalPanel("input.lp_layout == 'wide'",
                selectizeInput("lp_items_wide", "Item columns", NULL,
                               multiple = TRUE,
                               options = list(placeholder = "all remaining"))),
              radioButtons("lp_structure", "Facet structure",
                           c("Additive" = "additive",
                             "Interactive (item-by-facet)" = "interactive")),
              p(class = "text-muted small",
                "Additive: one severity per facet level. Interactive: additionally estimates item-by-facet effects (a rater harsh on particular items) — qualifies invariance."),
              conditionalPanel("input.lp_structure == 'interactive'",
                selectInput("lp_interaction", "Interacting facet", NULL)),
              p(class = "text-muted small",
                "Each item x facet combination is calibrated jointly; facet severities are reported with SEs and fit.")
            ),
            conditionalPanel("input.model_type == 'btl'",
              h6("One comparison per row"),
              selectInput("bt_a", "Object A column", NONE_CH),
              selectInput("bt_b", "Object B column", NONE_CH),
              conditionalPanel("!input.bt_response",
                selectInput("bt_win", "Winner column", NONE_CH),
                selectizeInput("bt_margin", "Margin of win (optional)", NULL,
                               options = list(placeholder = "none — dichotomous")),
                p(class = "text-muted small",
                  "The extent of the win (e.g. a little / much) as an ordered factor or increasing values; with the winner column it forms graded categories with no orientation bookkeeping. A winner value of \"tie\" or \"draw\" marks a tie (middle category); any other value matching neither object is treated as missing and the row dropped.")),
              selectizeInput("bt_response", "Graded response (optional)", NULL,
                             options = list(placeholder = "none — use winner")),
              p(class = "text-muted small",
                "A graded preference for the first object (e.g. much worse … much better), as an ordered factor or scores 0..m; overrides the winner column. Ties belong in a middle category."),
              selectInput("bt_judge", "Judge column (optional)", NONE_CH),
              conditionalPanel("input.bt_judge && input.bt_judge != '(none)'",
                selectizeInput("bt_order", "Judgment order (optional)", NULL,
                               options = list(placeholder = "none")),
                p(class = "text-muted small",
                  "Each judge's judgment sequence (timestamps or ranks). Enables the within-judge dependence analysis: exposure (a seen-before advantage) and carry-over (the pull of the judge's own earlier verdicts)."),
                selectizeInput("bt_jfactors", "Judge factors (optional)", NULL,
                               multiple = TRUE,
                               options = list(placeholder = "none")),
                p(class = "text-muted small",
                  "Nominate judge groupings (columns constant within judge, e.g. judge sex or background) to test differential object functioning (DIF) by judge group.")),
              conditionalPanel("!input.bt_response && !input.bt_margin",
                radioButtons("bt_ties", "Ties",
                             c("Drop" = "drop", "Half a win each" = "half"))),
              p(class = "text-muted small",
                "The Bradley-Terry-Luce model: the conditional (person-free) form of the dichotomous Rasch model, estimated by the same conventions. A judge column enables the judge fit table and clusters the standard errors by judge. Results appear on the Summary, Items, and Persons pages.")
            )),
          accordion_panel("Estimation options", icon = bs_icon("gear"),
            conditionalPanel("input.model_type == 'rasch'",
              radioButtons("thr_structure", "Threshold structure",
                           c("Partial credit (item-specific)" = "pcm",
                             "Rating scale (common across items)" = "rsm")),
              p(class = "text-muted small",
                "Dichotomous items are the one-threshold special case and need no setting. The rating parameterisation requires equal maximum scores; lr_test() compares the two."),
              conditionalPanel("input.thr_structure == 'pcm'",
                radioButtons("thr_mode", "Threshold estimation",
                             c("Free thresholds" = "free",
                               "Principal components (Andrich)" = "pc")),
                conditionalPanel("input.thr_mode == 'pc'",
                  selectInput("pc_rank", "Components",
                              c("Location only" = "1",
                                "+ spread (equal spread)" = "2",
                                "+ skewness" = "3",
                                "+ kurtosis (full PC)" = "4"),
                              selected = "4"),
                  p(class = "text-muted small",
                    "Thresholds follow a polynomial trend across categories; useful with sparse categories. Anchors cannot be combined with this option.")))),
            conditionalPanel(
              "input.model_type == 'btl' && (input.bt_response || input.bt_margin)",
              radioButtons("bt_thr", "Threshold structure",
                           c("Free symmetric" = "free",
                             "Principal components (spread)" = "pc")),
              p(class = "text-muted small",
                "PC pools the symmetric thresholds to the spread component so thinly used categories borrow strength; free estimates each threshold pair.")),
            checkboxInput("ng_auto", "Automatic class intervals (at least 50 per interval)", TRUE),
            conditionalPanel("!input.ng_auto",
              sliderInput("ng", "Class intervals", min = 2, max = 16, value = 8)),
            conditionalPanel("input.model_type == 'efrm'",
              selectInput("ef_se", "Standard errors",
                          c("Hybrid (fast)" = "hybrid",
                            "Full person bootstrap (slow, exact)" = "bootstrap")),
              numericInput("ef_reps", "Bootstrap replicates", value = 200,
                           min = 50, step = 50)),
            tags$details(class = "rasch-advanced",
              tags$summary("Advanced"),
              numericInput("run_adjN",
                           span("Adjust chi-square to N",
                                info_icon("Recomputes the item-trait and item chi-squares as if the sample size were N; useful with very large samples, where trivial misfit reaches significance.")),
                           value = NA, min = 50),
              numericInput("maxit", "Maximum iterations", value = 60, min = 5, step = 5),
              numericInput("tol", "Convergence criterion", value = 1e-8,
                           min = 1e-12, step = 1e-8),
              conditionalPanel("input.model_type == 'btl'",
                selectInput("bt_count", "Count column (optional)", NONE_CH)))),
          conditionalPanel("input.model_type == 'rasch'",
            accordion_panel("Scoring & anchors", icon = bs_icon("key"),
              fileInput("key_file",
                        span("Scoring key (CSV)",
                             info_icon("Columns item,key for a multiple-choice key — use \"A/C\" for a double key — or item,option,score for polytomous option scoring.")),
                        accept = ".csv", placeholder = "optional"),
              fileInput("anchor_file", "Anchors for equating (CSV: item,k,tau)",
                        accept = ".csv", placeholder = "optional"),
              radioButtons("anchor_type", "Anchor as",
                           c("Individual thresholds" = "individual",
                             "Average item locations" = "average"),
                           inline = TRUE),
              p(class = "text-muted small mt-1",
                "Anchors match by item name; rows for items not present are ignored. Individual anchoring fixes each listed threshold; average anchoring fixes each item's mean location (thresholds stay free). Save an anchor file from the Items page of a previous analysis.")))
        ),
        input_task_button("run", "Estimate", icon = bs_icon("play-fill"),
                          type = "primary", class = "w-100 btn-lg mt-2"),
        conditionalPanel("output.has_override",
          uiOutput("override_status"),
          actionButton("reset_override", "Reset overrides",
                       class = "btn-outline-warning w-100 mt-1")),
        p(class = "text-muted small mt-3",
          "Estimation: pairwise conditional maximum likelihood (Andrich & Luo 2003).",
          "Person measures: Warm weighted likelihood.")
      ),
      uiOutput("data_main")
    )
  )

# -------------------------------------------------------------- SUMMARY --
# bottom row built by the server: the likelihood-ratio card only applies to
# a PCM fit whose items share a common maximum score, so it hides otherwise
.lr_card <- function()
  card(info_header("Likelihood-ratio test (PCM vs rating)",
         "Compares the partial credit model against the more parsimonious rating parameterisation with common thresholds; a non-significant result supports the rating model."),
    card_body(
      p(class = "text-muted small",
        "Refits the same data with the rating (common threshold structure) parameterisation and compares the pairwise conditional log-likelihoods. A non-significant outcome supports adopting the simpler rating model; use the adjusted statistic for inference."),
      input_task_button("run_lr", "Run likelihood-ratio test",
                        type = "primary"),
      verbatimTextOutput("lr_txt"),
      rcode_details("lr")))

panel_summary <- nav_panel("Summary", value = "p_summary", icon = bs_icon("clipboard-data"),
    # Rasch fits (hidden while a paired-comparison fit is active)
    conditionalPanel("output.is_btl != true",
    uiOutput("vboxes"),
    # stat-box cards sit inside plain divs: the grid row would otherwise
    # stretch them to equal height and pad the shorter card mid-row
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      div(statCard("fitsum_tbl", "Test of fit",
        info = "The total item-trait chi-square tests the invariance of item ordering across the trait; a significant result means at least one item's difficulty is not invariant across class intervals (Andrich & Marais 2019). The CSV download carries the complete test-of-fit summary, including the fit residual moments and fit-location correlations.",
        footer = uiOutput("fitsum_notes"))),
      div(statCard("targeting_tbl", "Targeting & reliability",
        info = "How well the item thresholds cover the person distribution, with the reliability indices; the CSV download carries the complete targeting summary, including the location moments and item separation."))
    ),
    # server-rendered: the likelihood-ratio card only when it applies
    uiOutput("summary_bottom"),
    # test-level displays (the scale-range control drives the characteristic
    # curve and information plots; default -6..6 matches the plot defaults)
    div(class = "mt-3 mb-2", style = "max-width: 300px;",
        sliderInput("ts_rng", "Scale range (logits)", min = -8, max = 8,
                    value = c(-6, 6), step = 0.5, width = "100%")),
    accordion(id = "test_acc", open = "test_tcc", class = "mb-3",
      accordion_panel("Test characteristic curve", value = "test_tcc",
        plotCard("tcc")),
      accordion_panel("Test information", value = "test_tif",
        plotCard("tif",
          info = "Information across the scale, with the standard error of measurement (SEM = 1/sqrt(information)) overlaid.")),
      accordion_panel("Guttman scalogram", value = "test_guttman",
        plotCard("guttman", height = "640px")))),
    # paired-comparison (BTL) fits: the headline value boxes and the
    # test-of-fit summary table
    conditionalPanel("output.is_btl == true",
      uiOutput("btl_boxes"),
      layout_columns(col_widths = breakpoints(sm = 12, xl = 6),
        div(statCard("btl_fitsum_tbl", "Test of fit",
          info = "The pairwise chi-square tests the observed against the expected win proportions over every pair of objects; the object separation index is the paired-comparison counterpart of the PSI. Within-judge dependence effects (exposure and carry-over) appear when a judgment-order column was nominated. The CSV download carries the complete summary.",
          footer = uiOutput("btl_fitsum_notes")))))
  )

# ---------------------------------------------------------------- ITEMS --
panel_items <- nav_panel("Items", value = "p_items", icon = bs_icon("list-check"),
    # Rasch fits (hidden while a paired-comparison fit is active)
    conditionalPanel("output.is_btl != true",
    uiOutput("items_vboxes"),
    div(class = "mb-2 d-flex align-items-center gap-3 flex-wrap",
        div(class = "rasch-inline-check",
            tooltip(checkboxInput("show_obs", "Observed points", TRUE,
                                  width = "auto"),
                    "Show the observed class-interval points on the category and threshold curves.")),
        downloadButton("dl_anchors", "Save anchors (CSV: item,k,tau)",
                       class = "btn-outline-secondary btn-sm")),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      tableCard("items_tbl", "Item statistics",
        controls = cols_switch("items_full"),
                "Click a row to explore that item on the right. Fit residual ~ N(0,1) under fit.",
                info = "Cells are highlighted where a statistic indicates misfit: |fit residual| > 2.5, adjusted chi-square p < 0.05, outfit mean square outside 0.7-1.3 and infit outside the tighter 0.8-1.2 (conventional working bands; infit is information-weighted, so it varies less). No single flag column - read each statistic on its own terms.",
                footer = uiOutput("items_note")),
      navset_card_underline(
        id = "items_nav",
        # the tab strip stays clean: the selected-item title, display
        # settings, and batch downloads all live on the controls row below
        # (no display utility classes on the conditionalPanels themselves:
        # Bootstrap's !important would beat the inline display:none toggle)
        header = div(class = "d-flex align-items-center gap-3 flex-wrap",
          uiOutput("sel_item_title", inline = TRUE),
          # display settings: hidden on the Frequencies and Chi-square tabs,
          # where neither control affects the output
          conditionalPanel(
            "input.items_nav != 'Frequencies' && input.items_nav != 'Chi-square'",
            div(class = "d-flex align-items-center gap-4 flex-wrap",
              div(style = "width: 200px;",
                  sliderInput("ex_ng", "Class intervals", min = 2, max = 16,
                              value = 8, step = 1, width = "100%")),
              div(style = "width: 260px;",
                  sliderInput("ex_rng", "Scale range (logits)", min = -8, max = 8,
                              value = c(-5, 5), step = 0.5, width = "100%")))),
          # batch downloads follow the active tab's plot type; the Chi-square
          # tab has no plot, so the buttons hide there
          conditionalPanel("input.items_nav != 'Chi-square'",
            class = "ms-auto",
            div(class = "rasch-chips",
                downloadButton("items_all_pdf", "PDF (all items)",
                               class = "btn-outline-secondary btn-xs"),
                downloadButton("items_all_zip", "ZIP (all items)",
                               class = "btn-outline-secondary btn-xs")))),
        full_screen = TRUE,
        nav_panel("ICC",
                  plotOutput("icc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("icc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("icc_pdf", "PDF", class = "btn-outline-secondary btn-xs")),
                  rcode_details("icc")),
        nav_panel("Categories",
                  plotOutput("ccc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("ccc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("ccc_pdf", "PDF", class = "btn-outline-secondary btn-xs")),
                  rcode_details("ccc")),
        nav_panel("Thresholds",
                  plotOutput("tpc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("tpc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("tpc_pdf", "PDF", class = "btn-outline-secondary btn-xs")),
                  rcode_details("tpc")),
        nav_panel("Frequencies",
                  plotOutput("cfreq", height = "440px"),
                  div(class = "text-end",
                      downloadButton("cfreq_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("cfreq_pdf", "PDF", class = "btn-outline-secondary btn-xs")),
                  rcode_details("cfreq")),
        nav_panel("Chi-square",
                  uiOutput("chisq_caption"),
                  DTOutput("chisq_int_tbl"),
                  h6("Response categories by class interval", class = "mt-3"),
                  DTOutput("chisq_cat_tbl"),
                  div(class = "text-end mt-2",
                      downloadButton("chisq_int_csv", "Intervals CSV",
                                     class = "btn-outline-secondary btn-xs"),
                      downloadButton("chisq_cat_csv", "Categories CSV",
                                     class = "btn-outline-secondary btn-xs")),
                  rcode_details("chisq")))),
    accordion(id = "items_acc", open = "items_thrmap", class = "mt-3 mb-3",
      accordion_panel("Threshold map", value = "items_thrmap",
        plotCard("thrmap")),
      accordion_panel("Item fit map", value = "items_imap",
        plotCard("imap",
          info = "Item locations plotted against their fit residuals; items beyond |2.5| warrant inspection.")),
      accordion_panel("Fit residual distribution", value = "items_rdist",
        plotCard("rdist_i")),
      accordion_panel("Traditional statistics", value = "items_ctt",
        card(
          full_screen = TRUE,
          card_header_bar(
            buttons = downloadButton("ctt_tbl_csv", "CSV",
                                     class = "btn-outline-secondary btn-xs")),
          card_body(uiOutput("ctt_head"), DTOutput("ctt_tbl"),
                    rcode_details("ctt_tbl"),
                    padding = 12, fillable = FALSE)))),
    uiOutput("pc_comp_ui"),
    conditionalPanel("output.has_mc == true",
    layout_columns(col_widths = 12,
      tableCard("distractor_tbl", "Distractor analysis",
                "Locations use the rest measure; a distractor whose takers are abler than the keyed option's flags a possible miskey."),
      plotCard("distractor_plot", "Option curves"),
      card(card_header_bar("Polytomous option scoring (Andrich & Styles 2011)",
             buttons = conditionalPanel("output.has_rescore == true",
               downloadButton("dl_rescore", "Key CSV",
                              class = "btn-outline-secondary btn-xs"))),
           card_body(fillable = FALSE,
             p(class = "text-muted",
               "Propose partial credit for informative distractors from the rest-measure evidence. Review substantively, download, edit if needed, and upload as the key (item,option,score) to refit."),
             layout_columns(col_widths = c(3, 3, 3, 3),
               numericInput("rescore_min_n", "Min takers", 20, min = 5, step = 5),
               numericInput("rescore_z", "Separation z", 1.96, min = 0.5, step = 0.1),
               div(class = "mt-4",
                   input_task_button("rescore_go", "Propose option scores",
                                     type = "primary")),
               conditionalPanel("output.has_rescore == true", class = "mt-4",
                                cols_switch("rescore_full"))),
             conditionalPanel("output.has_rescore != true",
               p(class = "text-muted small mb-0",
                 "Run to see the rest-measure evidence and the proposed key.")),
             conditionalPanel("output.has_rescore == true",
               DT::DTOutput("rescore_tbl"),
               rcode_details("rescore_tbl"))))))),
    # paired-comparison (BTL) fits: the object side of the analysis
    conditionalPanel("output.is_btl == true",
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
        tableCard("btl_obj_tbl", "Object locations and fit",
          controls = cols_switch("btl_full"),
                  "Click a row to plot that object on the right. Conditional (person-free) estimation with sum-zero identification and sandwich standard errors; infit and outfit are the information-weighted and unweighted mean squares over each object's comparisons, and the fit residual is the log mean square (Andrich & Marais 2019).",
          info = "Cells are flagged where a statistic indicates misfit: outfit mean square outside 0.7-1.3, infit outside the tighter 0.8-1.2, and |fit residual| > 2.5."),
        plotCard("btl_occ", "Object characteristic curve",
          info = "The paired-comparison counterpart of the item characteristic curve: the model expected response for the selected object against opponent location (the win probability, or the expected graded response), with the observed mean response per opponent overlaid at that opponent's location. Opponents met too few times (sparse designs) are omitted. Points straying from the curve flag inconsistent quality, exactly as a misfitting item does.",
          extra = downloadButton("btl_occ_all_pdf", "PDF (all objects)",
                                 class = "btn-outline-secondary btn-xs"))),
      accordion(id = "btl_items_acc", open = "btl_caterpillar",
                class = "mt-3 mb-3",
        accordion_panel("Object caterpillar", value = "btl_caterpillar",
          plotCard("btl_plot")),
        # graded (ordinal) fits only: hidden entirely for dichotomous fits
        conditionalPanel("output.btl_graded == true",
          accordion_panel("Symmetric thresholds", value = "btl_thresholds",
            tableCard("btl_thr_tbl",
                      note = "Adjacent-categories thresholds of the graded structure, constrained symmetric (tau_k = -tau_(m+1-k)) so the model is invariant to presentation order.")),
          accordion_panel("Threshold components", value = "btl_components",
            tableCard("btl_comp_tbl",
                      note = "Spread is the linear component; the skewness component is structurally zero under presentation-order symmetry. Under the PC structure kurtosis is constrained to zero.")),
          accordion_panel("Category probability curves", value = "btl_catcurves",
            plotCard("btl_cats",
              info = "The probability of each graded response category as a function of the location difference between the two objects; the paired-comparison counterpart of a polytomous item's category curves."))),
        accordion_panel("Pairwise fit", value = "btl_pairs",
          tableCard("btl_pairs_tbl",
                    note = "Observed against expected win proportions (mean graded responses for a graded fit) for every pair; the total chi-square tests the BTL structure."))))
  )

# -------------------------------------------------------------- PERSONS --
panel_persons <- nav_panel("Persons", value = "p_persons", icon = bs_icon("people"),
    # Rasch fits (hidden while a paired-comparison fit is active)
    conditionalPanel("output.is_btl != true",
    uiOutput("persons_vboxes"),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      tableCard("person_tbl", "Person estimates",
          controls = cols_switch("persons_full"),
                "Warm WLE location and SE per person, with raw score, fit statistics, and your ID and factor columns. Click a row to draw that person's kidmap on the right."),
      plotCard("kidmap", "Kidmap",
        info = "The person diagnostic map (Wright, Mead & Ludlow 1980): thresholds the person achieved print to the right of the logit axis, thresholds not achieved to the left; the dashed line inside its confidence band is the person location. Achieved thresholds above the band and unachieved thresholds below it are unexpected responses.",
        controls = div(class = "d-flex align-items-center gap-1 me-1",
          span(class = "small text-secondary", "Confidence"),
          div(class = "rasch-inline-select",
              selectInput("kid_level", NULL,
                          c("90%" = "0.9", "95%" = "0.95", "99%" = "0.99"),
                          selected = "0.95", width = "85px"))),
        extra = tagList(
          downloadButton("kidmap_all_pdf", "PDF (all persons)",
                         class = "btn-outline-secondary btn-xs"),
          downloadButton("kidmap_all_zip", "ZIP (all persons)",
                         class = "btn-outline-secondary btn-xs")))),
    accordion(id = "persons_acc", open = "persons_pfit", class = "mt-3",
      accordion_panel("Person fit", value = "persons_pfit",
        plotCard("pfit")),
      accordion_panel("Fit residual distribution", value = "persons_rdist",
        plotCard("rdist_p")))),
    # paired-comparison (BTL) fits: the judges are the persons here, so the
    # page carries their fit and their transitivity consistency -- the two
    # judge-level lenses (offered only when a judge column was nominated)
    conditionalPanel("output.is_btl == true && output.has_judges == true",
      accordion(id = "btl_judge_acc", open = "btl_judge_fit",
        accordion_panel(
          title = span("Judge fit",
            info_icon("An erratic judge carries a large positive fit residual, exactly as an erratic person does; the log-of-mean-square residual and infit/outfit are pooled over the judge's comparisons.")),
          value = "btl_judge_fit",
          tableCard("btl_judges_tbl",
                    controls = cols_switch("btl_judges_full"))),
        accordion_panel(
          title = span("Judge consistency",
            info_icon("The paired-comparison counterpart of person fit. A judge whose choices form many preference loops (prefers A over B, B over C, then C over A) is internally inconsistent - not measuring on a single scale. Consistency is 1 minus the judge's circular-triad rate over the chance rate; 1 is one clean order, 0 is guessing.")),
          value = "btl_judge_consistency",
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(5, 7)),
            tableCard("btl_trans_judges_tbl", title = "Consistency by judge",
                      note = "Judges sorted least consistent first."),
            plotCard("btl_judge_consist", title = "Consistency dotplot",
                     info = "Each judge against the chance line at 0 and the clean-order line at 1. Judges near or below chance contribute little signal.",
                     height = "460px")))))
  )

# ------------------------------------------------------------ TARGETING --
panel_targeting <- nav_panel("Targeting", value = "p_targeting", icon = bs_icon("bullseye"),
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        sliderInput("tg_bins", "Histogram bins", min = 10, max = 60,
                    value = 35, step = 1),
        sliderInput("tg_rng", "Scale range (logits)", min = -10, max = 10,
                    value = c(-5, 5), step = 0.5),
        p(class = "text-muted small",
          "Targeting compares the person distribution with the item threshold distribution on the common logit scale. A well-targeted test places its thresholds where the persons are, so measurement error stays small across the range of the sample; gaps or offsets between the two distributions show where precision is lost.")),
      layout_columns(col_widths = 12,
        plotCard("pim_p", "Person-item threshold distribution"),
        plotCard("wright", "Wright map", height = "640px"))
    )
  )

# ------------------------------------------------------------------ DIF --
panel_dif <- nav_panel("DIF", value = "p_dif", icon = bs_icon("sliders"),
    # Rasch fits: person-factor DIF (hidden while a BTL fit is active)
    conditionalPanel("output.is_btl != true",
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        conditionalPanel("output.dif_multifactor == true",
          radioButtons("dif_effects", "Model",
                       c("Main effects" = "main",
                         "With interactions" = "factorial"))),
        numericInput("dif_alpha", "Significance level (alpha)", value = 0.05,
                     min = 0.001, max = 0.5, step = 0.01),
        selectInput("dif_padj", "Multiplicity adjustment",
                    c("Benjamini-Hochberg" = "BH",
                      "Holm" = "holm",
                      "Bonferroni" = "bonferroni",
                      "None" = "none")),
        p(class = "text-muted small",
          "ANOVA of standardised residuals: a factor effect is uniform DIF, a factor-by-class-interval term is non-uniform DIF. One factor is analysed one-way; several are modelled jointly with main effects by default. Adding interactions lets a significant interaction supersede its main effects. Click a row to see its characteristic curves by group and, below, the pairwise comparisons that resolve that term."),
        hr(),
        input_task_button("make_split", "Resolve the selected item",
                          type = "primary", class = "w-100"),
        p(class = "text-muted small mt-2",
          "Splits the selected analysis-of-variance row's item into independent copies by that row's factor(s) and re-analyses (the override)."),
        input_task_button("resolve_all", "Resolve all DIF automatically",
                          type = "primary", class = "w-100 mt-2"),
        p(class = "text-muted small mt-2",
          "Splits DIF items one at a time, largest effect first, refitting until no item shows significant DIF or the anchor set would fall too low (Andrich & Hagquist 2012)."),
        conditionalPanel("output.has_override_dif",
          actionButton("reset_split", "Reset to original data",
                       class = "btn-outline-warning w-100 mt-2"))),
      accordion(id = "dif_acc", open = "dif_anova",
        accordion_panel("DIF analysis of variance", value = "dif_anova",
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
            tableCard("dif_tbl",
              controls = cols_switch("dif_full"),
                      info = "ANOVA of standardised residuals: a significant factor effect indicates uniform DIF, a significant factor-by-class-interval interaction indicates non-uniform DIF (Andrich & Marais 2019). Click a row to see that item's characteristic curves by group (right) and, below, the pairwise comparisons for the selected item and term.",
                      footer = uiOutput("dif_note")),
            plotCard("dif_icc", "Characteristic curves by group",
              info = "The item characteristic curve drawn separately for each level of the person factor (in factorial mode, each factor-combination cell), with the observed class-interval means overlaid: the graphical display of DIF for the item of the selected table row."))),
        # collapsed by default: the DT output suspends while hidden, so the
        # per-item terms table is only computed when the panel is first opened
        accordion_panel("Full ANOVA table", value = "dif_full_panel",
          tableCard("dif_full_tbl",
                    note = "The complete per-item ANOVA: every model term with its df, sums of squares, mean squares, F, and adjusted probability.")),
        accordion_panel("Pairwise comparisons", value = "dif_pairwise",
          card(card_body(fillable = FALSE,
                 p(class = "text-muted",
                   "Resolves the item of the selected analysis-of-variance row by that row's factor (interaction rows resolve by the interaction cells) and reports the pairwise location differences in logits - the DIF magnitude - with Holm familywise adjustment. A two-level factor gives the single magnitude row."),
                 layout_columns(col_widths = c(4, 4, 4),
                   numericInput("dif_size_flag", "Practical criterion (logits)",
                                0.5, min = 0.1, step = 0.1),
                   numericInput("dif_size_minn", "Min responders", 20,
                                min = 5, step = 5),
                   div(class = "mt-4 d-flex justify-content-end align-items-start",
                       downloadButton("dl_dif_size", "CSV",
                                      class = "btn-outline-secondary btn-xs"))),
                 uiOutput("dif_levels_note"),
                 DT::DTOutput("dif_size_tbl"),
                 rcode_details("dif_size_tbl")))),
        # shown only after an automatic run: the trace of the splits that
        # resolved the DIF (the resolved fit is the active override)
        accordion_panel("Automatic resolution", value = "dif_resolve",
          conditionalPanel("output.has_resolve == true",
            card(
              card_header_bar(
                buttons = downloadButton("resolve_tbl_csv", "CSV",
                                         class = "btn-outline-secondary btn-xs")),
              card_body(fillable = FALSE,
                uiOutput("resolve_summary"),
                DT::DTOutput("resolve_tbl"))))),
        accordion_panel(
          title = span("Planned contrasts",
                       info_icon("Planned one-degree-of-freedom questions derived from the factor structure, tested with familywise control over the small planned family instead of all cell pairs (Maxwell & Delaney 2004). Estimates are DIF magnitudes in logits from resolved item locations. With a person ID and repeated rows, time-like factors are treated within-subjects via person-level residual scores.")),
          value = "dif_contrasts",
          card(card_body(fillable = FALSE,
                 p(class = "text-muted",
                   "Derives the family of questions from the factors themselves - a two-level factor contributes its difference, an ordered factor its linear and quadratic trends, a nominal factor its level comparisons, and factor pairs their product interaction - then tests the whole family at once."),
                 layout_columns(col_widths = c(4, 4, 4),
                   selectizeInput("pc_items", "Items", NULL, multiple = TRUE,
                                  options = list(placeholder = "all items")),
                   selectInput("pc_id", "Person ID (repeated measures)",
                               c("None" = "")),
                   div(class = "mt-4",
                       input_task_button("pc_run", "Derive and test contrasts",
                                         type = "primary"))),
                 conditionalPanel("output.has_contr != true",
                   p(class = "text-muted small mb-0",
                     "Run to see the derived family and its tests.")),
                 conditionalPanel("output.has_contr == true",
                   uiOutput("contr_family"),
                   div(class = "d-flex justify-content-end align-items-center gap-3",
                       cols_switch("contr_full"),
                       downloadButton("contr_tbl_csv", "CSV",
                                      class = "btn-outline-secondary btn-xs")),
                   DT::DTOutput("contr_tbl"),
                   rcode_details("contr_tbl"))))))
    )),
    # paired-comparison (BTL) fits: differential object functioning by
    # judge group (Bradley-Terry counterpart of the person-factor analysis)
    conditionalPanel("output.is_btl == true",
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectizeInput("bdif_factors", "Judge factors", NULL, multiple = TRUE,
                       options = list(placeholder = "nominate judge factors on the Data page")),
        conditionalPanel("output.bdif_multifactor == true",
          radioButtons("bdif_effects", "Model",
                       c("Main effects" = "main",
                         "With interactions" = "factorial"))),
        numericInput("bdif_alpha", "Significance level (alpha)", value = 0.05,
                     min = 0.001, max = 0.5, step = 0.01),
        input_task_button("bdif_run", "Run DIF analysis",
                          type = "primary", class = "w-100"),
        p(class = "text-muted small mt-2",
          "ANOVA of standardised residuals by judge factor: a factor effect indicates uniform DIF, a factor-by-opponent-band interaction non-uniform DIF. One factor is analysed on its own; several are modelled jointly with main effects by default, interactions optional. Each term flagged for uniform DIF (and not superseded by a higher-order term) is then resolved into one copy per cell inside a joint refit and the location differences reported in logits; withheld terms are named in the notes.")),
      accordion(id = "bdif_acc", open = "bdif_anova",
        accordion_panel("DIF analysis of variance", value = "bdif_anova",
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
            tableCard("bdif_anova_tbl",
              info = "ANOVA of the standardised residuals of each object's comparisons, oriented to the object: a significant judge-group effect indicates uniform DIF, a significant group-by-opponent-band interaction non-uniform DIF; probabilities are adjusted across objects by Benjamini-Hochberg. Click a row to see that object's characteristic curves by judge group on the right.",
              footer = uiOutput("bdif_notes")),
            plotCard("bdif_occ", "Characteristic curves by group",
              info = "The object characteristic curve with the observed mean response per opponent overlaid separately for each judge group: the graphical display of DIF for the object of the selected table row."))),
        accordion_panel("DIF magnitude in logits", value = "bdif_size_panel",
          tableCard("bdif_sizes_tbl",
            note = "Pairwise differences between the resolved per-group locations, in logits, with Benjamini-Hochberg adjustment across objects; differences of at least 0.5 logits are flagged as practically significant.")))
    ))
  )

# --------------------------------------------------------------- FACETS --
panel_facets <- nav_panel("Facets", value = "p_facets", icon = bs_icon("person-badge"),
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectizeInput("facet_sel", "Facet", NULL,
                       options = list(placeholder = "run a many-facet analysis")),
        p(class = "text-muted small",
          "Severities from the joint calibration (positive = more severe). Pooled fit residuals beyond +/-2.5 flag inconsistent levels. Many-facet (MFRM) analyses only.")),
      tableCard("facet_tbl", "Facet severities and fit",
        controls = cols_switch("facets_full"),
        footer = uiOutput("facet_structure_note")),
      conditionalPanel("output.has_interaction == true",
        tableCard("facet_int_tbl", "Item-by-facet interactions",
                  "Estimated because the interactive facet structure was chosen in the data roles; gamma is the extra severity of a level on a particular item, and significant terms qualify the invariance of the facet's severities across items.")),
      plotCard("facet_plot", "Severity caterpillar plot")
    )
  )

# -------------------------------------------------------------- EQUATING --
panel_equating <- nav_panel("Equating", value = "p_equating", icon = bs_icon("arrow-left-right"),
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        radioButtons("eq_source", "Reference",
                     c("Uploaded calibration CSV" = "csv",
                       "A kept fit from Compare" = "kept")),
        conditionalPanel("input.eq_source == 'csv'",
          fileInput("eq_file", "Reference calibration (CSV: item,location,se)",
                    accept = ".csv")),
        conditionalPanel("input.eq_source == 'kept'",
          selectizeInput("eq_kept", "Kept fit", NULL,
                         options = list(placeholder = "keep a fit on the Compare page"))),
        radioButtons("eq_shift", "Scale alignment",
                     c("Allow a shift between origins" = "mean",
                       "Compare raw locations (anchored scales)" = "none")),
        downloadButton("dl_calib", "Save current calibration (CSV)",
                       class = "btn-outline-secondary w-100"),
        p(class = "text-muted small mt-2",
          "Common items (matched by name) are tested against the shifted identity line; flagged items show drift and weaken the equating link. Save a calibration now to equate a future analysis against it.")),
      card(
        full_screen = TRUE,
        card_header_bar("Common-item comparison",
          buttons = conditionalPanel("output.has_eq == true",
            div(class = "rasch-chips",
                cols_switch("eq_full"),
                downloadButton("eq_tbl_csv", "CSV",
                               class = "btn-outline-secondary btn-xs")))),
        card_body(
          conditionalPanel("output.has_eq != true",
            p(class = "text-muted small mb-0",
              "Upload a reference calibration (or choose a kept fit) to see the common-item comparison.")),
          conditionalPanel("output.has_eq == true",
            DTOutput("eq_tbl"), rcode_details("eq_tbl")),
          padding = 12, fillable = FALSE)),
      card(
        full_screen = TRUE,
        `data-bs-theme` = "light",
        card_header_bar("Equating plot",
          buttons = conditionalPanel("output.has_eq == true",
            div(class = "rasch-chips",
                downloadButton("eq_plot_png", "PNG", class = "btn-outline-secondary btn-xs"),
                downloadButton("eq_plot_pdf", "PDF", class = "btn-outline-secondary btn-xs")))),
        card_body(
          conditionalPanel("output.has_eq != true",
            p(class = "text-muted small mb-0",
              "Upload a reference calibration (or choose a kept fit) to see the equating plot.")),
          conditionalPanel("output.has_eq == true",
            plotOutput("eq_plot", height = "560px"), rcode_details("eq_plot")),
          padding = 8, fillable = FALSE))
    )
  )

# --------------------------------------------------------------- FRAMES --
panel_frames <- nav_panel("Frames", value = "p_frames", icon = bs_icon("grid-3x3"),
    layout_sidebar(
      sidebar = sidebar(width = 290, open = "always",
        selectizeInput("frame_item", "Item for ICC across frames", NULL,
                       options = list(placeholder = "run a frames analysis")),
        p(class = "text-muted small",
          "Units rho = alpha (set) x phi (group) on a common arbitrary scale. Within a frame all curves are parallel; across frames they fan with the unit. Extended frame of reference analyses only.")),
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(7, 5)),
        tableCard("frame_tbl", "Frames: units, origins, pooled fit",
                  controls = cols_switch("frames_full")),
        div(tableCard("phi_tbl", "Person group units (phi)"),
            tableCard("alpha_tbl", "Item set units (alpha) and locations"))),
      layout_columns(col_widths = 12,
        plotCard("frame_plot", "Frame units"),
        plotCard("frame_icc", "ICC across frames")),
      card(card_header("Equal-unit comparison"),
           card_body(verbatimTextOutput("efrm_cmp")))
    )
  )

# -------------------------------------------------- INDEPENDENCE: TRAIT --
panel_dim <- nav_panel("Trait", value = "p_dim", icon = bs_icon("diagram-3"),
    # paired-comparison (BTL) fits: the Rasch residual-PCA suite needs a
    # persons x items residual matrix that paired comparisons do not produce,
    # so it hides and the pair-structure analogues take its place
    conditionalPanel("output.is_btl == true",
      p(class = "text-muted small",
        "Dimensionality for paired comparisons: is one scale enough to explain the contests? Two reads - a decomposition of the leftover (residual) preferences, and the rate of preference loops (Kendall & Babington Smith 1940)."),
      accordion(id = "btl_dim_acc", open = "btl_dim_swirl",
        accordion_panel(
          title = span("Residual dimensions",
            info_icon("The paired-comparison counterpart of residual PCA. The model predicts how often each object should beat each other; the object-by-object table of departures is skew-symmetric, so it decomposes into rotational planes (bimensions; Gower 1977), not ordinary components. A leading bimension clearing the model-simulated noise band is a coherent swirl in the leftovers - A over-beats B, B over-beats C, C over-beats A - a second attribute steering some contests.")),
          value = "btl_dim_swirl",
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
            plotCard("btl_scree", title = "Bimension strengths",
                     info = "Each bimension's strength against the mean and 95th-percentile band of data simulated from the fitted one-scale model. A bar clearing the band is structure the single scale does not explain.",
                     height = "460px"),
            plotCard("btl_dim_map", title = "Leading residual map",
                     info = "Objects in the leading bimension plane. A rotational arrangement is the second attribute; a formless cloud at the centre is noise. Point size grows with the object's location on the main scale.",
                     height = "460px")),
          tableCard("btl_bimensions_tbl", title = "Bimensions",
                    note = "Strength, share of the total residual, and the noise reference (mean and 95th percentile) for the leading bimension.")),
        accordion_panel(
          title = span("Preference loops",
            info_icon("The single-dimension check. If one attribute drives the contests, preferences stack into one order: A beats B and B beats C implies A beats C. A loop (A beats B, B beats C, C beats A) is a contradiction, like rock-paper-scissors. The loop rate is set against pure guessing (a quarter of triples); consistency is 1 minus loop-rate over chance. Each judge's own consistency is on the Persons tab.")),
          value = "btl_dim_loops",
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(7, 5)),
            tableCard("btl_trans_tbl", title = "Transitivity summary",
                      note = "Circular triads out of the complete triples, the loop rate against the 25% chance rate, and Kendall's coefficient of consistency when every pair was compared."),
            plotCard("btl_involve_plot", title = "Objects in loops",
                     info = "How many circular triads each object sits in - the objects whose order is least stable, and the likeliest seat of a second attribute.",
                     height = "460px"))))),
    conditionalPanel("output.is_btl != true",
    p(class = "text-muted small",
      "Trait dependence (dimensionality) threatens local independence: more than one trait driving the responses (Marais & Andrich 2008)."),
    accordion(id = "dim_acc", open = "dim_components",
      accordion_panel("Residual components", value = "dim_components",
        p(class = "text-muted small",
          "Loadings of each item on the leading residual principal components, with a biplot of the first two - typically the only interpretable contrasts. Items with opposing signs and large magnitude on PC1 mark a possible second dimension; PC2 separates them further."),
        layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
          tableCard("loadings_tbl", title = "Loadings",
                    note = "First 10 components shown."),
          plotCard("pca_biplot", title = "Biplot (PC1 vs PC2)",
                   height = "auto"))),
      accordion_panel("Scree plot", value = "dim_scree",
        plotCard("scree")),
      accordion_panel(
        title = span("Unidimensionality t-test",
                     info_icon("Smith's test: each person is measured separately on the two item subsets and the estimates compared by t-test; unidimensionality is questioned when clearly more than 5% of tests are significant.")),
        value = "dim_ttest",
        layout_columns(col_widths = breakpoints(sm = 12, xl = c(4, 8)),
          div(
            h6("t-test item subsets"),
            div(class = "mb-2 d-flex align-items-center gap-2",
              span(class = "small text-secondary", "Automatic split component"),
              div(class = "rasch-inline-select",
                  selectInput("pca_component", NULL, choices = 1, selected = 1,
                              width = "80px"))),
            selectizeInput("dim_pos", "Subset A", NULL, multiple = TRUE,
                           options = list(placeholder = "positive loadings on the selected component")),
            selectizeInput("dim_neg", "Subset B", NULL, multiple = TRUE,
                           options = list(placeholder = "negative loadings on the selected component")),
            input_task_button("dim_apply", "Run t-test",
                              type = "primary", class = "w-100"),
            p(class = "text-muted small mt-2",
              "Leave both empty (and press the button) to return to the split from the selected component. Persons extreme on either subset are excluded; the proportion of significant tests carries an exact binomial confidence interval.")),
          card(card_body(verbatimTextOutput("dim_txt"), rcode_details("dim"))))),
      accordion_panel("Magnitude of multidimensionality", value = "dim_magnitude",
        card(
          full_screen = TRUE,
          card_body(
            p(class = "text-muted small",
              "Compares reliability with all items treated as independent (run1) against the subtest analysis in which each subset becomes one polytomous super-item (Andrich 2016). c is the unique-variance loading, rho the latent correlation between the subsets, and A the proportion of common variance. Uses the manual subsets above if set, otherwise the selected component's split; every item must belong to a subset."),
            div(input_task_button("dm_run", "Estimate from current subsets",
                                  type = "primary")),
            conditionalPanel("output.has_dm != true",
              p(class = "text-muted small mb-0 mt-2",
                "Run to see the resolved reliability comparison.")),
            conditionalPanel("output.has_dm == true",
              div(class = "d-flex justify-content-end",
                  downloadButton("dm_tbl_csv", "CSV",
                                 class = "btn-outline-secondary btn-xs")),
              DTOutput("dm_tbl"), rcode_details("dm_tbl")),
            padding = 12, fillable = FALSE))),
      accordion_panel("Eigenvalues", value = "dim_eigen",
        tableCard("eigen_tbl", note = "First 10 eigenvalues shown."))))
  )

# -------------------------------------------------- INDEPENDENCE: LOCAL --
panel_ld <- nav_panel("Local", value = "p_ld", icon = bs_icon("link-45deg"),
    # paired-comparison (BTL) fits: within-judge dependence estimated from
    # the judgment order; the Rasch Q3 suite hides while a BTL fit is active
    conditionalPanel("output.is_btl == true",
      p(class = "text-muted small",
        "Local (response) dependence threatens local independence (Marais & Andrich 2008): a judge's own history pulling their later judgments, over and above the object locations."),
      conditionalPanel("output.has_btl_dep != true",
        card(card_body(p(class = "text-muted small mb-0",
          "Nominate a judgment-order column in the Data roles to estimate within-judge dependence.")))),
      conditionalPanel("output.has_btl_dep == true",
        layout_columns(col_widths = breakpoints(sm = 12, lg = c(5, 7)),
          tableCard("btl_dep_tbl", "Within-judge dependence",
            info = "Exposure is the seen-before advantage: the benefit, in logits, an object gains once the judge has already met it. Carry-over is response dependence (Marais & Andrich): the judge's own earlier verdicts on an object pull the later one. Both are estimated jointly with the object locations. Informative comparisons is how many comparisons carry information about each effect - those where the two objects' histories differ (a non-zero covariate).",
            note = "An effect estimated on few informative comparisons carries a wide standard error; the plot and the comparison table below show which comparisons drive it."),
          plotCard("btl_dep_plot", "Dependence effect",
            info = "The counterpart of the DIF characteristic curve: the observed departure from the location-only prediction, binned by the effect's history covariate, with the model's fitted contribution overlaid and the count in each bin printed. Observed points rising with the covariate along the line are the effect; a flat, sparse cloud means the estimate rests on little.",
            controls = div(class = "d-flex align-items-center gap-1 me-1",
              span(class = "small text-secondary", "Effect"),
              div(class = "rasch-inline-select",
                  selectInput("btl_dep_effect", NULL,
                              c("Exposure" = "exposure",
                                "Carry-over" = "carry_over"),
                              width = "130px"))))),
        accordion(class = "mt-3",
          accordion_panel("Comparison covariates", value = "btl_dep_comps_panel",
            tableCard("btl_dep_comps",
              note = "Every comparison in judgment order with its exposure and carry-over covariates; a comparison is informative for an effect when its covariate is non-zero."))))),
    conditionalPanel("output.is_btl != true",
    p(class = "text-muted small",
      "Local (response) dependence threatens local independence (Marais & Andrich 2008): responses depending on one another directly, over and above the trait."),
    accordion(id = "ld_acc", open = "ld_cormat",
      accordion_panel("Residual Correlations (Q3 statistics)", value = "ld_cormat",
        p(class = "text-muted small",
          "Yen's (1984) Q3 is the correlation of the standardised residuals for an item pair; Q3* subtracts the average off-diagonal Q3, so 0 is the local-independence baseline and a pair well above it signals response dependence. Each matrix shows the lower triangle only, beside its heatmap. A pair is shown in red when it clears the flag threshold under its own rule: |Q3| for the raw matrix (Yen 1993) and Q3* for the adjusted matrix (Christensen, Makransky & Horton 2017); 0.2 is conventional for both."),
        numericInput("ld_flag",
                     "Flag threshold (|Q3| or Q3* at or above this value)",
                     value = 0.2, min = 0.05, max = 0.9, step = 0.05,
                     width = "420px"),
        layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
          tableCard("cormat_q3_tbl", title = "Q3 correlations",
                    info = "The residual correlation of every item pair, with 1.00 on the diagonal (Yen 1984). Pairs with |Q3| at or above the threshold are red (Yen 1993)."),
          plotCard("rcor_q3", height = "auto")),
        layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
          tableCard("cormat_q3s_tbl", title = "Adjusted Q3 (Q3*)",
                    info = "Each Q3 less the average off-diagonal Q3: 0 marks local independence. Pairs with Q3* at or above the threshold are red (Christensen, Makransky & Horton 2017)."),
          plotCard("rcor_q3s", height = "auto"))),
      accordion_panel("Response dependence magnitude", value = "ld_dep",
        card(
          full_screen = TRUE,
          card_body(
            p(class = "text-muted small",
              "Resolves the dependent item by the categories of the independent item and re-analyses (Andrich & Kreiner); d is the size of the dependence in logits, half the split of the resolved thresholds. Both items must share the same maximum score."),
            div(class = "d-flex gap-3 flex-wrap align-items-end",
              selectizeInput("dep_item", "Dependent item", NULL, width = "190px",
                             options = list(placeholder = "run an analysis first")),
              selectizeInput("ind_item", "Independent item", NULL, width = "190px",
                             options = list(placeholder = "run an analysis first")),
              div(class = "mb-3",
                  input_task_button("run_dep", "Estimate d", type = "primary"))),
            conditionalPanel("output.has_dep != true",
              p(class = "text-muted small mb-0",
                "Run to see the resolved magnitudes.")),
            conditionalPanel("output.has_dep == true",
              div(class = "d-flex justify-content-end",
                  downloadButton("dep_tbl_csv", "CSV",
                                 class = "btn-outline-secondary btn-xs")),
              verbatimTextOutput("dep_txt"),
              DTOutput("dep_tbl"), rcode_details("dep_tbl")),
            padding = 12, fillable = FALSE))),
      accordion_panel("Subtest (combine dependent items)", value = "ld_subtest",
        card(
          card_body(
            p(class = "text-muted small",
              "Select two or more items to merge into one polytomous super-item and re-analyse; the dependence is absorbed into the subtest."),
            selectizeInput("subtest_items", NULL, NULL, multiple = TRUE,
                           options = list(placeholder = "items to combine")),
            div(input_task_button("make_subtest", "Combine and re-analyse",
                                  type = "primary")),
            conditionalPanel("output.has_override_subtest",
              div(class = "mt-2",
                actionButton("reset_subtest", "Reset to original data",
                             class = "btn-outline-warning w-100"))),
            uiOutput("subtest_status")))),
      accordion_panel("Spread test (LUB)", value = "ld_spread",
        card(
          full_screen = TRUE,
          card_body(
            p(class = "text-muted small",
              "Spread below the least upper bound indicates dependence among subtest members (Andrich 1985). Polytomous items only; typically applied after combining items into a subtest."),
            div(class = "mb-2",
                input_task_button("run_spread", "Run spread test",
                                  type = "primary")),
            conditionalPanel("output.has_spread != true",
              p(class = "text-muted small mb-0",
                "Run to see the spread of each item against its least upper bound.")),
            conditionalPanel("output.has_spread == true",
              div(class = "d-flex justify-content-end",
                  downloadButton("spread_tbl_csv", "CSV",
                                 class = "btn-outline-secondary btn-xs")),
              DTOutput("spread_tbl"),
              rcode_details("spread_tbl")),
            padding = 12, fillable = FALSE)))))
  )

# ------------------------------------------------------------- GUESSING --
panel_guess <- nav_panel("Guessing", value = "p_guess", icon = bs_icon("question-diamond"),
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        numericInput("guess_chance", "Chance success probability",
                     value = 0.25, min = 0.05, max = 0.95, step = 0.05),
        selectizeInput("guess_anchors", "Anchor items (common origin)", NULL,
                       multiple = TRUE,
                       options = list(placeholder = "automatic: least-affected third")),
        input_task_button("run_guess", "Run tailored analysis",
                          type = "primary", class = "w-100"),
        p(class = "text-muted small mt-2",
          "The tailored procedure of Andrich, Marais and Humphry (2012): every response whose modelled success probability falls below the chance level is set to missing and the test is re-calibrated on a common origin. Difficult items becoming harder in the tailored calibration signals guessing. Dichotomous analyses only.")),
      layout_columns(col_widths = 12,
        card(card_header("Tailored analysis"),
             card_body(verbatimTextOutput("guess_txt"))),
        tableCard("guess_tbl", "Initial vs tailored calibration",
                  "shift = tailored minus origin-equated location; z > 1.96 flags items significantly harder after tailoring (a guessing signature).")),
      plotCard("guess_plot", "Tailored vs origin-equated calibrations")
    )
  )

# -------------------------------------------------------------- COMPARE --
panel_compare <- nav_panel("Compare", value = "p_compare", icon = bs_icon("columns-gap"),
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        input_task_button("keep_fit", "Keep current fit for comparison",
                          type = "primary", class = "w-100"),
        actionButton("clear_fits", "Clear kept fits",
                     class = "btn-outline-secondary w-100 mt-2"),
        selectizeInput("cmp_ref", "Reference fit", NULL,
                       options = list(placeholder = "keep at least two fits")),
        p(class = "text-muted small mt-3",
          "Run an analysis, keep it, change the model or settings, run again, and keep that too. For fits of the same data the pairwise conditional log-likelihoods are compared directly (descriptive, composite likelihood; most meaningful for nested structures such as RSM inside PCM). Across different data preparations, compare the calibration-free columns: chi-square per df, fit residual SDs (ideal 1), PSI, and alpha.")),
      card(
        full_screen = TRUE,
        card_header_bar("Model comparison",
          buttons = conditionalPanel("output.has_cmp == true",
            div(class = "rasch-chips",
                cols_switch("cmp_full"),
                downloadButton("cmp_tbl_csv", "CSV",
                               class = "btn-outline-secondary btn-xs")))),
        card_body(
          conditionalPanel("output.has_cmp != true",
            p(class = "text-muted small mb-0",
              "Keep at least two fits (run, keep, change the settings, run and keep again) to see the comparison.")),
          conditionalPanel("output.has_cmp == true",
            p(class = "text-muted small mb-2",
              "Reference for the log-likelihood comparison is the fit chosen in the sidebar."),
            DTOutput("cmp_tbl"), rcode_details("cmp_tbl")),
          padding = 12, fillable = FALSE))
    )
  )

# --------------------------------------------------------------- EXPORT --
panel_export <- nav_panel("Export", value = "p_export", icon = bs_icon("download"),
    conditionalPanel("output.is_btl == true",
      card(card_body(class = "empty-state",
        bs_icon("trophy", size = "2rem",
                class = "text-secondary d-block mx-auto mb-2"),
        p("The HTML report and ZIP archive cover Rasch analyses. For a paired-comparison (BTL) analysis, download each table as CSV from its card on the Summary, Items, and Persons pages.")))),
    conditionalPanel("output.is_btl != true",
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(card_header(span(bs_icon("file-earmark-text"),
                            " Analysis report (single HTML file)")),
        card_body(
          p("A self-contained HTML report of the current analysis: the summary statistics, every diagnostic table, and the test-level plots embedded as images. One portable file, ready to e-mail or archive."),
          p(class = "text-muted small",
            "Available for Rasch analyses (dichotomous, PCM, RSM, MFRM, EFRM); paired-comparison (BTL) fits are not covered."),
          downloadButton("dl_report", "Download report (HTML)",
                         class = "btn-primary btn-lg", icon = icon("file")))),
      card(card_header("Download everything"),
        card_body(
          p("One archive containing every table (CSV), every plot (in the formats chosen), and a plain-text analysis summary."),
          checkboxGroupInput("exp_formats", "Plot formats",
                             c("PNG" = "png", "PDF" = "pdf"), selected = c("png", "pdf"),
                             inline = TRUE),
          checkboxInput("exp_items", "Include the per-item plot set (ICC, categories, thresholds, frequencies for every item)", TRUE),
          downloadButton("dl_zip", "Download all results (ZIP)", class = "btn-primary btn-lg"))),
      card(card_header("What is included"),
        card_body(tags$ul(
          tags$li("Item statistics with the item ANOVA fit table, thresholds with SEs, person estimates (with ID and factors), score-to-measure table with score frequencies"),
          tags$li("Chi-square class-interval detail for every item, traditional (CTT) statistics, and the principal components table when estimated"),
          tags$li("Residual correlations, flagged dependent pairs, PCA loadings, category frequencies, DIF ANOVA for every factor"),
          tags$li("Person-item distribution, threshold map, TCC, TIF, item and person fit maps, item and person fit residual distributions, residual heatmap, PCA plot"),
          tags$li("Per-item ICC, category curves, threshold curves, and frequency charts"),
          tags$li("For many-facet analyses: facet severities with SEs and fit, structural item thresholds, and severity caterpillar plots"),
          tags$li("summary.txt with the full test-of-fit report"))))
    ))
  )

# ------------------------------------------------------------ ASSEMBLY --
# Workflow order: data -> summary -> items -> persons -> test, then the
# independence, invariance, and utility menus (the two requirements of
# measurement); status chips and the dark-mode toggle sit at the right of
# the navbar.
ui <- page_navbar(
  id = "nav",
  title = span(class = "app-brand",
               span("rasch", class = "app-brand-name"),
               span("Rasch Measurement Theory", class = "app-brand-sub")),
  theme = theme,
  # normal scrolling pages: never compress content to fit the viewport
  fillable = FALSE,
  header = tagList(
    tags$head(tags$style(css),
      # nav visibility by data-value: shiny::hideTab (behind bslib's
      # nav_hide) does not reach nav_panels nested inside a nav_menu
      # dropdown, so the server toggles entries itself through this handler;
      # it covers top-level links, dropdown items, and menu toggles alike
      tags$script(HTML("
        Shiny.addCustomMessageHandler('rasch-nav-vis', function(msg) {
          document.querySelectorAll('.navbar a[data-value]').forEach(function(a) {
            if (a.getAttribute('data-value') !== msg.value) return;
            var li = a.closest('li');
            (li || a).style.display = msg.show ? '' : 'none';
          });
        });
      "))),
    busyIndicatorOptions(spinner_type = "ring2")),
  panel_data,
  panel_summary,
  panel_items,
  panel_persons,
  panel_targeting,
  nav_menu("Independence", value = "menu_independence",
    panel_ld,
    panel_dim),
  nav_menu("Invariance", value = "menu_invariance",
    panel_dif,
    panel_equating,
    panel_guess,
    panel_facets,
    panel_frames),
  nav_menu("More", value = "menu_more",
    panel_compare,
    panel_export),
  nav_spacer(),
  nav_item(uiOutput("nav_status")),
  nav_item(downloadLink("dl_report_nav", label = bs_icon("file-earmark-text"),
                        class = "nav-link px-2",
                        title = "Analysis report (HTML)")),
  nav_item(input_dark_mode())
)

server <- function(input, output, session) {

  # ------------------------------------------------------------- data in --
  # picking an example dataset also selects the matching model; uploading a
  # file clears the example selection
  observeEvent(input$demo_choice, {
    dc <- input$demo_choice
    if (!identical(dc, "none")) {
      updateRadioButtons(session, "model_type",
                         selected = if (dc %in% c("dich", "pcm", "rsm"))
                           "rasch" else dc)
      if (dc %in% c("dich", "pcm", "rsm"))
        updateRadioButtons(session, "thr_structure",
                           selected = if (identical(dc, "rsm")) "rsm" else "pcm")
    }
  }, ignoreInit = TRUE)
  observeEvent(input$file,
    updateSelectInput(session, "demo_choice", selected = "none"))

  raw_data <- reactive({
    if (!identical(input$demo_choice %||% "none", "none"))
      return(switch(input$demo_choice,
                    dich = .demo_dich(), rsm = .demo_rsm(),
                    mfrm = .demo_mfrm(), efrm = .demo_efrm(),
                    btl = .demo_btl(), .demo_data()))
    req(input$file)
    ext <- tolower(tools::file_ext(input$file$name))
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    read.csv(input$file$datapath, sep = sep, check.names = FALSE,
             stringsAsFactors = FALSE)
  })

  anchors_in <- reactive({
    if (is.null(input$anchor_file)) return(NULL)
    a <- tryCatch(read.csv(input$anchor_file$datapath, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(a) || !all(c("item", "k", "tau") %in% names(a))) {
      showNotification("Anchor file needs columns item, k, tau - ignored.",
                       type = "warning")
      return(NULL)
    }
    a
  })

  observeEvent(raw_data(), {
    df <- raw_data(); nm <- names(df)
    guess_id <- nm[grepl("^id$|_id$|^person", tolower(nm))][1]
    guess_fac <- intersect(nm, c("group", "sex", "gender", "site", "country", "age_group"))
    updateSelectInput(session, "id_col", choices = c(NONE_CH, nm),
                      selected = if (!is.na(guess_id)) guess_id else NONE)
    updateSelectizeInput(session, "factor_cols", choices = nm, selected = guess_fac)
    updateSelectizeInput(session, "item_cols", choices = nm,
                         selected = setdiff(nm, c(guess_id, guess_fac)))
    # MFRM role guesses
    g_per <- nm[grepl("person|candidate|student|^id$|_id$", tolower(nm))][1]
    g_itm <- nm[grepl("item|task|criterion|question", tolower(nm))][1]
    g_sco <- nm[grepl("score|rating|grade|mark", tolower(nm))][1]
    g_fac <- setdiff(nm[grepl("rater|judge|marker|occasion|time", tolower(nm))],
                     c(g_per, g_itm, g_sco))
    updateSelectInput(session, "lp_person", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_per)) g_per else NONE)
    updateSelectInput(session, "lp_item", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_itm)) g_itm else NONE)
    updateSelectInput(session, "lp_score", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_sco)) g_sco else NONE)
    updateSelectizeInput(session, "lp_facets", choices = nm, selected = g_fac)
    # wide layout: item columns = the remaining columns (like the Rasch
    # item_cols guess), excluding the guessed person and facet columns
    updateSelectizeInput(session, "lp_items_wide", choices = nm,
                         selected = setdiff(nm, c(g_per, g_fac)))
    # frames layout guesses
    g_grp <- nm[grepl("group|year|grade|cohort|class$", tolower(nm))][1]
    updateSelectInput(session, "ef_id", choices = c(NONE_CH, nm),
                      selected = if (!is.na(guess_id)) guess_id else NONE)
    updateSelectInput(session, "ef_group", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_grp)) g_grp else NONE)
    updateSelectizeInput(session, "ef_items", choices = nm,
                         selected = setdiff(nm, c(guess_id, g_grp)))
    # paired-comparison guesses
    g_a <- nm[grepl("^a$|object_a|left|first|option_a", tolower(nm))][1]
    g_b <- nm[grepl("^b$|object_b|right|second|option_b", tolower(nm))][1]
    g_w <- nm[grepl("win|preferred|chosen|better", tolower(nm))][1]
    g_j <- nm[grepl("judge|rater|marker", tolower(nm))][1]
    g_c <- nm[grepl("^count$|^n$|freq", tolower(nm))][1]
    g_o <- nm[grepl("^t$|order|seq|trial|round|^time$", tolower(nm))][1]
    # judge factors are columns constant within each judge (and not another
    # role) -- exactly the shape of a panel or rater-group variable, so they
    # are offered for judge-group DIF by default
    g_f <- character(0)
    if (!is.na(g_j)) {
      cand <- setdiff(nm, stats::na.omit(c(g_a, g_b, g_w, g_j, g_o, g_c)))
      if (length(cand))
        g_f <- cand[vapply(cand, function(cn) all(tapply(as.character(df[[cn]]),
                     df[[g_j]], function(v) length(unique(v)) == 1L)), TRUE)]
    }
    updateSelectInput(session, "bt_a", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_a)) g_a else NONE)
    updateSelectInput(session, "bt_b", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_b)) g_b else NONE)
    updateSelectInput(session, "bt_win", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_w)) g_w else NONE)
    # empty choice = the "none — …" placeholder (clearable)
    updateSelectizeInput(session, "bt_response", choices = c("", nm),
                         selected = "")
    updateSelectizeInput(session, "bt_margin", choices = c("", nm),
                         selected = "")
    updateSelectInput(session, "bt_judge", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_j)) g_j else NONE)
    # empty choice = the "none" placeholder (clearable)
    updateSelectizeInput(session, "bt_order", choices = c("", nm),
                         selected = if (!is.na(g_o)) g_o else "")
    updateSelectizeInput(session, "bt_jfactors", choices = nm,
                         selected = g_f)
    updateSelectInput(session, "bt_count", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_c)) g_c else NONE)
  })

  # item -> set map: uploaded CSV wins; otherwise infer from the item-name
  # prefix (the part before trailing digits/separators)
  ef_setmap <- reactive({
    its <- if (length(input$ef_items)) input$ef_items else
      setdiff(names(raw_data()),
              c(if (!is.null(input$ef_id) && input$ef_id != NONE) input$ef_id,
                if (!is.null(input$ef_group) && input$ef_group != NONE) input$ef_group))
    if (!is.null(input$ef_sets)) {
      mp <- tryCatch(read.csv(input$ef_sets$datapath, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      if (!is.null(mp) && all(c("item", "set") %in% names(mp))) {
        out <- setNames(as.character(mp$set), mp$item)
        miss <- setdiff(its, names(out))
        if (length(miss)) out[miss] <- "(rest)"
        return(out[its])
      }
      showNotification("Item-set CSV needs columns item,set - using prefixes.",
                       type = "warning")
    }
    if (isTRUE(input$ef_prefix)) {
      pref <- sub("[_. -]*[0-9]+$", "", its)
      pref[pref == ""] <- "(rest)"
      return(setNames(pref, its))
    }
    setNames(rep("all", length(its)), its)
  })

  # the interacting facet (shown in interactive facet mode) is chosen from
  # the nominated facet columns; a single facet is preselected automatically
  observeEvent(input$lp_facets, {
    fs <- input$lp_facets
    sel <- if (!is.null(input$lp_interaction) && input$lp_interaction %in% fs)
      input$lp_interaction else if (length(fs)) fs[1] else character(0)
    updateSelectInput(session, "lp_interaction",
                      choices = if (length(fs)) fs else character(0),
                      selected = sel)
  }, ignoreNULL = FALSE)

  # keep the wide-mode item choices free of the chosen person / facet columns
  observeEvent(c(input$lp_person, input$lp_facets), {
    df <- raw_data(); nm <- names(df)
    taken <- c(if (!is.null(input$lp_person) && input$lp_person != NONE)
                 input$lp_person,
               input$lp_facets)
    sel <- setdiff(if (length(input$lp_items_wide)) input$lp_items_wide else nm,
                   taken)
    updateSelectizeInput(session, "lp_items_wide",
                         choices = setdiff(nm, taken), selected = sel)
  }, ignoreInit = TRUE)

  # keep item choices free of the chosen ID / factor columns
  observeEvent(c(input$id_col, input$factor_cols), {
    df <- raw_data(); nm <- names(df)
    taken <- c(if (!is.null(input$id_col) && input$id_col != NONE) input$id_col,
               input$factor_cols)
    sel <- setdiff(if (length(input$item_cols)) input$item_cols else nm, taken)
    updateSelectizeInput(session, "item_cols", choices = setdiff(nm, taken),
                         selected = sel)
  }, ignoreInit = TRUE)

  # Data page main area: an empty-state hero before any data is loaded, the
  # summary strip + preview + R-code disclosure once data is in
  .demo_labels <- c(dich = "Multiple choice, dichotomous",
                    pcm = "Polytomous (PCM)",
                    rsm = "Rating scale (RSM)",
                    mfrm = "Ratings by raters (MFRM)",
                    efrm = "Item sets x groups (EFRM)",
                    btl = "Paired comparisons (BTL)")
  .demo_chip_labels <- c(dich = "Multiple choice", pcm = "Polytomous (PCM)",
                         rsm = "Rating scale", mfrm = "Ratings (MFRM)",
                         efrm = "Frames (EFRM)", btl = "Paired comparisons")
  output$data_main <- renderUI({
    if (identical(input$demo_choice %||% "none", "none") && is.null(input$file)) {
      div(class = "mx-auto", style = "max-width: 760px; margin-top: 8vh;",
        card(
          card_body(class = "empty-state",
            bs_icon("clipboard-data", size = "3rem", class = "text-primary mb-2"),
            h2("Welcome to rasch"),
            p(class = "lead mb-4",
              "Rasch measurement: pairwise conditional estimation, the complete test-of-fit suite, and every diagnostic table and plot — with one-click export."),
            div(class = "d-flex justify-content-center gap-2 flex-wrap mb-4",
              actionButton("hero_demo", "Try an example dataset",
                           icon = icon("table"), class = "btn-primary btn-lg"),
              tags$button(class = "btn btn-outline-primary btn-lg", type = "button",
                          onclick = "document.getElementById('file').click();",
                          bs_icon("upload"), " Upload data")),
            p(class = "text-muted small mb-2", "Or start from a specific example:"),
            div(class = "d-flex flex-wrap gap-1 justify-content-center",
              lapply(names(.demo_chip_labels), function(k)
                actionButton(paste0("demo_chip_", k), .demo_chip_labels[[k]],
                             class = "btn-outline-secondary btn-sm"))))))
    } else {
      tagList(
        uiOutput("data_strip"),
        card(card_header("Data preview"),
             card_body(uiOutput("data_info"), DTOutput("preview"),
                       padding = 12, fillable = FALSE)),
        accordion(id = "rcode_acc", open = FALSE, class = "mt-3",
          accordion_panel("R code for this analysis", icon = bs_icon("code-slash"),
            p(class = "text-muted small mb-2",
              "The exact rasch call reproducing the current run; updates on every estimation."),
            verbatimTextOutput("rcode_fit"))))
    }
  })
  observeEvent(input$hero_demo,
    updateSelectInput(session, "demo_choice", selected = "pcm"))
  lapply(c("dich", "pcm", "rsm", "mfrm", "efrm", "btl"), function(k)
    observeEvent(input[[paste0("demo_chip_", k)]],
      updateSelectInput(session, "demo_choice", selected = k)))
  output$data_strip <- renderUI({
    df <- raw_data()
    vals <- as.matrix(df)
    miss <- 100 * mean(is.na(vals) | trimws(vals) == "", na.rm = FALSE)
    layout_column_wrap(width = "160px", fill = FALSE, class = "mb-3",
      value_box("Rows", format(nrow(df), big.mark = ","),
                showcase = glyph("grid"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Columns", ncol(df), showcase = glyph("columns"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Missing", sprintf("%.1f%%", miss),
                showcase = glyph("missing"),
                showcase_layout = "left center",
                theme = if (miss > 20) "warning" else "secondary"))
  })
  output$data_info <- renderUI({
    df <- raw_data()
    p(class = "text-muted",
      sprintf("%d rows x %d columns.%s Nominate the column roles in the sidebar, then press Estimate. Missing responses may be left blank or coded as -1; any negative score is read as missing.",
              nrow(df), ncol(df),
              if (nrow(df) > 200) " First 200 rows shown in the preview." else ""))
  })
  output$preview <- renderDT({
    datatable(head(raw_data(), 200), rownames = FALSE, style = "bootstrap5",
              class = "table-sm compact hover order-column",
              options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })
  output$rcode_fit <- renderText({
    validate(need(!is.null(rcode_str()),
                  "Run an analysis to see the reproducible R code."))
    rcode_str()
  })

  # ----------------------------------------------------------------- fit --
  override_fit <- reactiveVal(NULL)
  override_desc <- reactiveVal(NULL)
  # the exact rasch call reproducing the current run (built alongside the fit)
  rcode_str <- reactiveVal(NULL)
  # clear any subtest/split override as soon as a fresh run is requested;
  # fit() short-circuits on the override, so analysis() cannot clear it itself
  observeEvent(input$run, { override_fit(NULL); override_desc(NULL) },
               priority = 10)
  output$has_override <- reactive(!is.null(override_fit()))
  outputOptions(output, "has_override", suspendWhenHidden = FALSE)
  # per-kind flags, so each local reset button shows (and undoes) only its own
  # restructure: the subtest button reverts a subtest, the DIF button reverts a
  # split or an automatic resolution. The single override_fit means only one
  # kind is ever active, so exactly one local button is visible at a time.
  output$has_override_subtest <- reactive(
    grepl("^subtest", override_desc() %||% ""))
  output$has_override_dif <- reactive(
    grepl("^(split|auto-resolved)", override_desc() %||% ""))
  outputOptions(output, "has_override_subtest", suspendWhenHidden = FALSE)
  outputOptions(output, "has_override_dif", suspendWhenHidden = FALSE)
  output$override_status <- renderUI({
    if (is.null(override_desc())) return(NULL)
    p(class = "text-warning small mb-1 mt-2",
      paste("Active override -", override_desc()))
  })
  # clearing an override returns every page to the base analysis; the automatic
  # DIF-resolution trace is dropped too, so nothing stale lingers. The same
  # reset is offered locally beside each transforming action (subtest, split,
  # automatic resolution) as well as on the Data page.
  reset_to_original <- function() {
    override_fit(NULL); override_desc(NULL); resolve_res(NULL)
    showNotification("Reset to the original data; showing the base analysis.",
                     type = "message", duration = 5)
  }
  observeEvent(input$reset_override, reset_to_original())
  observeEvent(input$reset_split,   reset_to_original())
  observeEvent(input$reset_subtest, reset_to_original())

  # estimation runs as a side-effecting observer that stores the completed
  # fit in fit_val; analysis() is a pure accessor, so every reader keeps
  # working unchanged
  fit_val <- reactiveVal(NULL)
  analysis <- reactive(fit_val())
  btl_fit <- reactiveVal(NULL)
  # shared estimation-control resolution (used by every model branch and
  # the reproducible-code footers)
  est_opts <- reactive(list(maxit = max(5, input$maxit %||% 60),
                            tol = max(1e-12, input$tol %||% 1e-8)))
  observeEvent(input$run, {
    df <- raw_data()
    # automatic class intervals pass NULL; rasch() resolves the rule and
    # reports the value used in fit$n_groups
    ng <- if (isTRUE(input$ng_auto %||% TRUE)) NULL else input$ng
    # chi-square sample-size adjustment is applied inside the fit, so every
    # tab and export sees the same adjusted statistics
    adjN <- if (!is.null(input$run_adjN) && !is.na(input$run_adjN) &&
                input$run_adjN > 0) input$run_adjN else NA
    # reproducible-code pieces (spliced into the branch-specific call below)
    src_line <- if (!identical(input$demo_choice %||% "none", "none"))
      paste0("# dat: the \"", .demo_labels[[input$demo_choice]],
             "\" example dataset generated by the app")
    else
      sprintf('dat <- read.csv("%s"%s, check.names = FALSE)',
              input$file$name,
              if (tolower(tools::file_ext(input$file$name)) %in%
                  c("tsv", "txt")) ', sep = "\\t"' else "")
    code_args_common <- c(
      if (!is.null(ng)) paste0("n_groups = ", ng),
      if (!is.na(adjN)) paste0("adjust_N = ", adjN))
    eo <- est_opts()
    code_est <- c(paste0("maxit = ", eo$maxit),
                  paste0("tol = ", format(eo$tol)))
    code_call <- NULL
    code_notes <- character(0)
    withProgress(message = "Estimating (pairwise conditional ML)…", value = 0.3, {
      fit <- tryCatch({
        if (identical(input$model_type, "btl")) {
          # a graded response column overrides the winner column (and the
          # ties rule: graded ties belong in a middle category); otherwise a
          # margin column combines with the winner into the graded response
          # (winner values matching neither object = ties = middle category)
          bt_graded <- !is.null(input$bt_response) && nzchar(input$bt_response)
          bt_marg <- !bt_graded && !is.null(input$bt_margin) &&
            nzchar(input$bt_margin)
          bt_thr <- input$bt_thr %||% "free"
          # the judgment-order column enables the within-judge dependence
          # analysis (exposure and carry-over); it only exists with a judge
          bt_ord <- if (!is.null(input$bt_order) && nzchar(input$bt_order) &&
                        !is.null(input$bt_judge) && input$bt_judge != NONE)
            input$bt_order else NULL
          if (any(c(input$bt_a, input$bt_b) == NONE) ||
              (!bt_graded && identical(input$bt_win, NONE)))
            stop("nominate the object A, object B, and winner (or graded response) columns")
          code_call <- paste0("bt <- btl(dat,\n  ", paste(c(
            paste0("object_a = ", qstr(input$bt_a)),
            paste0("object_b = ", qstr(input$bt_b)),
            if (bt_graded) paste0("response = ", qstr(input$bt_response))
            else paste0("winner = ", qstr(input$bt_win)),
            if (bt_marg) paste0("margin = ", qstr(input$bt_margin)),
            if (!is.null(input$bt_judge) && input$bt_judge != NONE)
              paste0("judge = ", qstr(input$bt_judge)),
            if (!is.null(bt_ord)) paste0("order = ", qstr(bt_ord)),
            if (!is.null(input$bt_count) && input$bt_count != NONE)
              paste0("count = ", qstr(input$bt_count)),
            if (!bt_graded && !bt_marg)
              paste0("ties = ", qstr(input$bt_ties %||% "drop")),
            if ((bt_graded || bt_marg) && identical(bt_thr, "pc"))
              'thresholds = "pc"',
            code_est), collapse = ",\n  "), ")")
          # one shared argument list; the entry path (graded response,
          # winner + margin, winner only) contributes its own arguments
          bt_args <- c(
            list(df, object_a = input$bt_a, object_b = input$bt_b,
                 judge = if (!is.null(input$bt_judge) &&
                             input$bt_judge != NONE) input$bt_judge else NULL,
                 order = bt_ord,
                 count = if (!is.null(input$bt_count) &&
                             input$bt_count != NONE) input$bt_count else NULL,
                 maxit = eo$maxit, tol = eo$tol),
            if (bt_graded)
              list(response = input$bt_response, thresholds = bt_thr)
            else if (bt_marg)
              list(winner = input$bt_win, margin = input$bt_margin,
                   thresholds = bt_thr)
            else
              list(winner = input$bt_win, ties = input$bt_ties %||% "drop"))
          do.call(btl, bt_args)
        } else if (identical(input$model_type, "efrm")) {
          sm <- ef_setmap()
          code_call <- paste0("fit <- rasch_efrm(dat,\n  ", paste(c(
            paste0("item_sets = ",
                   paste(deparse(sm), collapse = "\n    ")),
            if (is.null(input$ef_group) || input$ef_group == NONE)
              'groups = rep("(all)", nrow(dat))'
            else paste0("groups = ", qstr(input$ef_group)),
            if (!is.null(input$ef_id) && input$ef_id != NONE)
              paste0("id = ", qstr(input$ef_id)),
            paste0("items = ", qvec(names(sm))),
            code_args_common,
            code_est,
            paste0("se_method = ", qstr(input$ef_se %||% "hybrid")),
            if (!is.null(input$ef_reps) && !is.na(input$ef_reps))
              paste0("boot_reps = ", max(50, input$ef_reps))),
            collapse = ",\n  "), ")")
          rasch_efrm(df,
                     item_sets = sm,
                     groups = if (is.null(input$ef_group) ||
                                  input$ef_group == NONE)
                       rep("(all)", nrow(df)) else input$ef_group,
                     id = if (!is.null(input$ef_id) && input$ef_id != NONE)
                       input$ef_id else NULL,
                     items = names(sm),
                     n_groups = ng, adjust_N = adjN,
                     maxit = eo$maxit, tol = eo$tol,
                     se_method = input$ef_se %||% "hybrid",
                     boot_reps = if (!is.null(input$ef_reps) &&
                                     !is.na(input$ef_reps))
                       max(50, input$ef_reps) else NULL)
        } else if (identical(input$model_type, "mfrm")) {
          # the interaction is passed only in interactive facet mode, and
          # only when the interacting facet is one of the chosen facets
          lp_int <- if (identical(input$lp_structure %||% "additive",
                                  "interactive") &&
                        !is.null(input$lp_interaction) &&
                        input$lp_interaction %in% input$lp_facets)
            input$lp_interaction else NULL
          if (identical(input$lp_layout %||% "wide", "wide")) {
            if (identical(input$lp_person, NONE) || !length(input$lp_facets) ||
                !length(input$lp_items_wide))
              stop("nominate the person column, the item columns, and at least one facet column")
            code_call <- paste0("fit <- rasch_mfrm(dat,\n  ", paste(c(
              paste0("person = ", qstr(input$lp_person)),
              paste0("facets = ", qvec(input$lp_facets)),
              paste0("items = ", qvec(input$lp_items_wide)),
              code_args_common,
              if (!is.null(lp_int)) paste0("interaction = ", qstr(lp_int)),
              code_est), collapse = ",\n  "), ")")
            rasch_mfrm(df, person = input$lp_person,
                       facets = input$lp_facets, items = input$lp_items_wide,
                       n_groups = ng, adjust_N = adjN, interaction = lp_int,
                       maxit = eo$maxit, tol = eo$tol)
          } else {
            if (any(c(input$lp_person, input$lp_item, input$lp_score) == NONE) ||
                !length(input$lp_facets))
              stop("nominate the person, item, score, and at least one facet column")
            code_call <- paste0("fit <- rasch_mfrm(dat,\n  ", paste(c(
              paste0("person = ", qstr(input$lp_person)),
              paste0("item = ", qstr(input$lp_item)),
              paste0("score = ", qstr(input$lp_score)),
              paste0("facets = ", qvec(input$lp_facets)),
              code_args_common,
              if (!is.null(lp_int)) paste0("interaction = ", qstr(lp_int)),
              code_est), collapse = ",\n  "), ")")
            rasch_mfrm(df, person = input$lp_person, item = input$lp_item,
                       score = input$lp_score, facets = input$lp_facets,
                       n_groups = ng, adjust_N = adjN, interaction = lp_int,
                       maxit = eo$maxit, tol = eo$tol)
          }
        } else {
          idc <- if (!is.null(input$id_col) && input$id_col != NONE) input$id_col else NULL
          fac <- if (length(input$factor_cols)) input$factor_cols else NULL
          its <- if (length(input$item_cols)) input$item_cols else NULL
          # multiple-choice key: uploaded CSV, or the demo key for the demo
          mc_key <- NULL
          if (!is.null(input$key_file)) {
            kf <- tryCatch(read.csv(input$key_file$datapath,
                                    stringsAsFactors = FALSE),
                           error = function(e) NULL)
            if (!is.null(kf) && (all(c("item", "key") %in% names(kf)) ||
                                 all(c("item", "option", "score") %in% names(kf))))
              mc_key <- kf
            else showNotification("Key CSV needs columns item,key or item,option,score - ignored.",
                                  type = "warning")
          } else if (identical(input$demo_choice, "dich")) {
            mc_key <- .demo_dich_key()
          }
          # anchors match by item name; rows for absent items are ignored
          anc <- anchors_in()
          if (!is.null(anc)) {
            cand <- if (!is.null(its)) its else setdiff(names(df), c(idc, fac))
            present <- as.character(anc$item) %in% cand
            if (!all(present))
              showNotification(sprintf("%d anchor row(s) ignored (items not in this dataset)",
                                       sum(!present)), type = "warning")
            anc <- anc[present, , drop = FALSE]
            if (!nrow(anc)) anc <- NULL
            # average anchoring: collapse to one mean-location anchor per item
            if (!is.null(anc) && identical(input$anchor_type, "average")) {
              mu <- tapply(anc$tau, as.character(anc$item), mean)
              anc <- data.frame(item = names(mu), k = NA, tau = as.numeric(mu))
            }
          }
          # threshold structure: partial credit (item-specific) or rating
          # (common); dichotomous items are the one-threshold special case
          rsm_on <- identical(input$thr_structure %||% "pcm", "rsm")
          # principal-components (Andrich) threshold estimation; partial
          # credit route only, and not combinable with anchors
          pcc <- if (!rsm_on && identical(input$thr_mode, "pc"))
            as.integer(input$pc_rank %||% "4") else NULL
          if (!is.null(pcc) && !is.null(anc)) {
            showNotification("Anchors are ignored under principal-components threshold estimation.",
                             type = "warning", duration = 8)
            anc <- NULL
          }
          code_notes <- c(
            if (!is.null(anc))
              "# anchor rows for items not in the data are dropped before fitting",
            if (!is.null(anc) && identical(input$anchor_type, "average"))
              "# average anchoring: the anchor file is collapsed to one mean location per item")
          code_call <- paste0("fit <- rasch(dat,\n  ", paste(c(
            paste0("model = ", qstr(if (rsm_on) "RSM" else "PCM")),
            if (!is.null(idc)) paste0("id = ", qstr(idc)),
            if (!is.null(fac)) paste0("factors = ", qvec(fac)),
            if (!is.null(its)) paste0("items = ", qvec(its)),
            code_args_common,
            if (!is.null(anc) && !is.null(input$anchor_file))
              paste0("anchors = read.csv(", qstr(input$anchor_file$name), ")"),
            if (!is.null(mc_key) && !is.null(input$key_file))
              paste0("key = read.csv(", qstr(input$key_file$name), ")")
            else if (!is.null(mc_key))
              'key = setNames(rep("A", 15), sprintf("I%02d", 1:15))',
            if (!is.null(pcc)) paste0("pc_components = ", pcc),
            code_est), collapse = ",\n  "), ")")
          rasch(df, model = if (rsm_on) "RSM" else "PCM",
                id = idc, factors = fac, items = its,
                n_groups = ng, adjust_N = adjN, anchors = anc,
                key = mc_key, pc_components = pcc,
                maxit = eo$maxit, tol = eo$tol)
        }
      }, error = function(e) e)
    })
    if (inherits(fit, "error")) {
      showNotification(paste("Analysis failed:", conditionMessage(fit)),
                       type = "error", duration = 10)
      btl_fit(NULL)
      fit_val(NULL)
      return(invisible(NULL))
    }
    if (!is.null(code_call))
      rcode_str(paste(c("library(rasch)", "", src_line,
                        if (length(code_notes)) c("", code_notes), "",
                        code_call), collapse = "\n"))
    # routine handling notes are informational; only real problems warn
    if (length(fit$notes))
      showNotification(paste(fit$notes, collapse = "\n"), type = "message",
                       duration = 8)
    conv <- if (!is.null(fit$est)) fit$est$converged else fit$converged
    if (!isTRUE(conv))
      showNotification("Estimation did not converge; consider raising the maximum iterations or loosening the convergence criterion.",
                       type = "warning", duration = 10)
    override_fit(NULL); override_desc(NULL)
    # paired-comparison results render on the Summary / Items / Persons
    # pages (each page's Rasch variant hides while a BTL fit is current,
    # and vice versa); the Rasch outputs suspend meanwhile
    if (inherits(fit, "rasch_btl")) {
      btl_fit(fit)
      fit_val(NULL)
      try(nav_select("nav", "p_summary", session = session), silent = TRUE)
      return(invisible(NULL))
    }
    btl_fit(NULL)
    try(nav_select("nav", "p_summary", session = session), silent = TRUE)
    # reactiveVal skips notification when the new value is identical to the
    # old; clearing first guarantees the on-fit reset observer fires even
    # when a re-run reproduces the same fit
    fit_val(NULL)
    fit_val(fit)
  })
  fit <- reactive({
    f <- override_fit()
    if (is.null(f)) f <- analysis()
    req(f); f
  })
  # the same override-first resolution without req(): NULL before any run,
  # for UI that must render quietly in that state (navbar chips, report)
  fit_or_null <- function() {
    f <- override_fit()
    if (is.null(f)) f <- tryCatch(analysis(), error = function(e) NULL)
    f
  }

  output$has_mc <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    !is.null(f) && !is.null(f$mc)
  })
  outputOptions(output, "has_mc", suspendWhenHidden = FALSE)

  output$sel_item_title <- renderUI(span(class = "fw-semibold",
    tryCatch(sel_item(), error = function(e) "Selected item")))

  # only offer the pages that apply to the current analysis: Facets needs a
  # many-facet fit, Frames an extended-frames fit, and Guessing a
  # dichotomous one. Summary and Items serve Rasch and paired-comparison
  # (BTL) fits alike (each page shows the matching variant); Persons needs
  # person estimates (Rasch) or a judge column (BTL). Everything else stays.
  # Every nav_panel and nav_menu carries an explicit value, and visibility is
  # driven by those values through the rasch-nav-vis handler (shiny::hideTab,
  # which backs bslib::nav_hide, cannot reach entries inside a nav_menu).
  observe({
    f <- tryCatch(fit(), error = function(e) NULL)
    bf <- btl_fit()
    show <- function(value, on)
      session$sendCustomMessage("rasch-nav-vis",
                                list(value = value, show = isTRUE(on)))
    show("p_facets", inherits(f, "rasch_mfrm"))
    show("p_frames", inherits(f, "rasch_efrm"))
    show("p_guess", !is.null(f) && !inherits(f, "rasch_mfrm") &&
           !inherits(f, "rasch_efrm") && max(f$m) == 1L)
    rasch_on <- !is.null(f)
    btl_on <- !is.null(bf)
    # judge-factor DIF applies to a paired-comparison fit once judge
    # factors are nominated in the Data roles
    btl_dif_on <- btl_on && length(input$bt_jfactors) > 0
    show("p_summary", rasch_on || btl_on)
    show("p_items", rasch_on || btl_on)
    show("p_persons", rasch_on || (btl_on && !is.null(bf$judges)))
    for (tgt in c("p_targeting", "p_equating"))
      show(tgt, rasch_on)
    # the Trait tab now carries paired-comparison dimensionality too
    # (transitivity loops + the residual bimension swirl)
    show("p_dim", rasch_on || btl_on)
    # Local dependence: the Rasch Q3 suite, or the within-judge dependence
    # analysis of a paired-comparison fit
    show("p_ld", rasch_on || btl_on)
    # DIF needs at least one person factor in the fit (Rasch), or a
    # paired-comparison fit with judge factors nominated
    show("p_dif", (rasch_on && !is.null(f$factors) &&
                     length(names(f$factors)) > 0) || btl_dif_on)
    # Compare applies to Rasch fits only; for a paired-comparison fit the
    # More menu offers Export alone (its BTL empty state points to the
    # per-card CSV downloads)
    show("p_compare", rasch_on)
    # menu headers hide too when everything inside them is hidden
    show("menu_independence", rasch_on || btl_on)
    show("menu_invariance", rasch_on || btl_dif_on)
    show("menu_more", rasch_on || btl_on)
  })

  # ------------------------------------------------ UI visibility flags --
  # "nothing to show" areas hide instead of rendering empty; each flag pairs
  # with a conditionalPanel in the UI (same pattern as has_mc)
  output$has_interaction <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    inherits(f, "rasch_mfrm") && !is.null(f$interaction)
  })
  outputOptions(output, "has_interaction", suspendWhenHidden = FALSE)
  output$is_btl <- reactive(!is.null(btl_fit()))
  outputOptions(output, "is_btl", suspendWhenHidden = FALSE)
  # a paired-comparison fit with a judge column: the Persons page shows the
  # judge fit table (and is offered in the navbar) only then
  output$has_judges <- reactive({
    b <- btl_fit()
    !is.null(b) && !is.null(b$judges)
  })
  outputOptions(output, "has_judges", suspendWhenHidden = FALSE)
  # empty-state flags for the run-on-demand cards: before a run the card
  # shows only its controls and one muted line (no table, plot, or download)
  output$has_dep <- reactive(!is.null(dep_res()))
  outputOptions(output, "has_dep", suspendWhenHidden = FALSE)
  output$has_spread <- reactive(!is.null(spread_res()))
  outputOptions(output, "has_spread", suspendWhenHidden = FALSE)
  output$has_contr <- reactive(!is.null(contr_res()))
  outputOptions(output, "has_contr", suspendWhenHidden = FALSE)
  output$has_dm <- reactive(!is.null(dm_res()))
  outputOptions(output, "has_dm", suspendWhenHidden = FALSE)
  output$has_rescore <- reactive(!is.null(rescore_res()))
  outputOptions(output, "has_rescore", suspendWhenHidden = FALSE)
  output$has_cmp <- reactive(length(kept_fits()) >= 2)
  outputOptions(output, "has_cmp", suspendWhenHidden = FALSE)
  # equating results exist once a reference is available (upload or kept fit)
  output$has_eq <- reactive({
    if (identical(input$eq_source, "kept"))
      !is.null(input$eq_kept) && nzchar(input$eq_kept) &&
        input$eq_kept %in% names(kept_fits())
    else !is.null(input$eq_file)
  })
  outputOptions(output, "has_eq", suspendWhenHidden = FALSE)

  observeEvent(fit(), {
    its <- fit()$items$item
    updateSelectizeInput(session, "subtest_items", choices = its, selected = character(0))
    fac <- names(fit()$factors)
    updateSelectizeInput(session, "pc_items", choices = its,
                         selected = character(0))
    updateSelectInput(session, "pc_id", choices = c("None" = "", fac),
                      selected = "")
    fs <- if (inherits(fit(), "rasch_mfrm")) fit()$facet_spec else character(0)
    updateSelectizeInput(session, "facet_sel", choices = fs,
                         selected = if (length(fs)) fs[1] else character(0))
    fi <- if (inherits(fit(), "rasch_efrm"))
      unique(fit()$virtual_map$item) else character(0)
    updateSelectizeInput(session, "frame_item", choices = fi,
                         selected = if (length(fi)) fi[1] else character(0))
    updateSelectizeInput(session, "dim_pos", choices = its, selected = character(0))
    updateSelectizeInput(session, "dim_neg", choices = its, selected = character(0))
    # residual principal components available for the t-test default split:
    # 1..min(10, n_items - 1)
    updateSelectInput(session, "pca_component",
                      choices = seq_len(max(1L, min(10L, length(its) - 1L))),
                      selected = 1)
    updateSelectInput(session, "dep_item", choices = its,
                      selected = its[min(2L, length(its))])
    updateSelectInput(session, "ind_item", choices = its, selected = its[1])
    updateSelectizeInput(session, "guess_anchors", choices = its,
                         selected = character(0))
    # explorer class intervals start at the fit's own rule
    updateSliderInput(session, "ex_ng", value = fit()$n_groups)
    # targeting range: the padded person + threshold data range
    r <- range(c(fit()$person$theta, fit()$thresholds$tau), na.rm = TRUE)
    updateSliderInput(session, "tg_rng",
                      value = c(floor((r[1] - 0.4) * 2) / 2,
                                ceiling((r[2] + 0.4) * 2) / 2))
    # results computed on request belong to the fit they came from
    lr_res(NULL); dep_res(NULL); spread_res(NULL); dm_res(NULL); guess_res(NULL)
    contr_res(NULL); rescore_res(NULL)
    # an automatic resolution sets the override fit itself, so its trace must
    # survive its own refit; a fresh run or any other override clears it
    if (!grepl("^auto-resolved", override_desc() %||% "")) resolve_res(NULL)
    # manual dimensionality subsets too: they name items of the previous fit
    dim_subsets(NULL)
  })

  observeEvent(input$make_subtest, {
    req(length(input$subtest_items) >= 2)
    res <- tryCatch(combine_items(fit(), list(input$subtest_items)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Subtest failed:", conditionMessage(res)), type = "error")
    } else {
      override_fit(res)
      override_desc(paste("subtest:", paste(input$subtest_items, collapse = " + ")))
      showNotification("Re-analysed with the subtest in place. Use Reset to original data (or run again) to return to the base fit.",
                       type = "message", duration = 8)
    }
  })
  output$subtest_status <- renderUI({
    if (is.null(override_desc())) return(NULL)
    p(class = "text-success small mt-2", paste("Active:", override_desc()))
  })

  observeEvent(input$make_split, {
    f <- fit()
    it <- dif_sel_item()
    vars <- dif_sel_vars()
    req(it %in% f$items$item,
        !is.null(f$factors), length(vars) >= 1, all(vars %in% names(f$factors)))
    # one factor -> split by the factor name; an interaction row -> split by
    # the factor-combination cells, so both uniform and non-uniform DIF resolve
    by <- if (length(vars) == 1L) vars
          else interaction(f$factors[vars], sep = ":", drop = TRUE)
    res <- tryCatch(split_items(f, it, by = by), error = function(e) e)
    lab <- paste(vars, collapse = ":")
    if (inherits(res, "error")) {
      showNotification(paste("Split failed:", conditionMessage(res)), type = "error")
    } else {
      override_fit(res)
      override_desc(sprintf("split: item %s by %s", it, lab))
      showNotification(
        sprintf("Re-analysed with %s split by %s. Use Reset to original data (or run again) to return to the base fit.",
                it, lab),
        type = "message", duration = 8)
    }
  })

  # automatic iterative DIF resolution: split the largest-effect item, refit,
  # repeat until no item flags DIF (or the anchor set would fall too low). The
  # resolved fit becomes the active override; the trace is kept for the panel.
  resolve_res <- reactiveVal(NULL)
  observeEvent(input$resolve_all, {
    rr <- tryCatch(
      resolve_dif(fit(), factors = names(fit()$factors), alpha = dif_alpha(),
                  p_adjust = input$dif_padj %||% "BH"),
      error = function(e) e)
    if (inherits(rr, "error")) {
      showNotification(paste("Automatic resolution failed:", conditionMessage(rr)),
                       type = "error")
    } else {
      resolve_res(rr)
      override_fit(rr$fit)
      override_desc(sprintf("auto-resolved DIF: %d split(s)", rr$n_splits))
      showNotification(
        sprintf("Automatic DIF resolution: %d split(s); %s. Use Reset to original data (or run again) to return to the base fit.",
                rr$n_splits, rr$stopped),
        type = "message", duration = 8)
    }
  })
  output$has_resolve <- reactive(!is.null(resolve_res()))
  outputOptions(output, "has_resolve", suspendWhenHidden = FALSE)
  output$resolve_summary <- renderUI({
    rr <- resolve_res(); req(!is.null(rr))
    p(class = "text-muted small mb-2",
      sprintf("%d split(s); %s; %d item(s) still flag DIF.",
              rr$n_splits, rr$stopped, rr$n_remaining_dif))
  })
  output$resolve_tbl <- DT::renderDT({
    rr <- resolve_res(); req(!is.null(rr))
    d <- rr$splits
    if (nrow(d)) { d$eta2 <- round(d$eta2, 3); d$magnitude <- round(d$magnitude, 3) }
    num_dt(d)
  })
  output$resolve_tbl_csv <- downloadHandler(
    filename = function() "dif_resolution.csv",
    content = function(file) {
      rr <- resolve_res(); req(!is.null(rr))
      write.csv(rr$splits, file, row.names = FALSE)
    })

  sel_item <- reactive({
    f <- fit()
    i <- input$items_tbl_rows_selected
    if (length(i)) f$items$item[i] else f$items$item[1]
  })

  # ------------------------------------------------------- plot plumbing --
  # per-output "R code" disclosure: `code` is a function returning the exact
  # rasch call reproducing the output (it may read reactives, so the snippet
  # follows the current selections). Rendering is never suspended: the text
  # must be ready when the collapsed <details> footer is opened.
  register_code <- function(id, code) {
    cid <- paste0(id, "_code")
    output[[cid]] <- renderText(code())
    outputOptions(output, cid, suspendWhenHidden = FALSE)
  }
  # `px` optionally sets a reactive on-screen height (a function returning
  # pixels), for plots whose natural size grows with the data (e.g. a matrix
  # heatmap that should track its item count and the table beside it)
  register_plot <- function(id, fun, w = 9, h = 6, code = NULL, px = NULL) {
    output[[id]] <- if (is.null(px)) renderPlot(fun(), res = 96)
                    else renderPlot(fun(), res = 96, height = px)
    if (!is.null(code)) register_code(id, code)
    for (fmt in c("png", "pdf")) local({
      fmt_ <- fmt
      output[[paste0(id, "_", fmt_)]] <- downloadHandler(
        filename = function() paste0("rasch_", id, ".", fmt_),
        content = function(file) {
          # 300 dpi PNG (and vector PDF) for publication
          if (fmt_ == "png") png(file, width = w, height = h, units = "in", res = 300)
          else pdf(file, width = w, height = h)
          fun(); dev.off()
        })
    })
  }
  # `csv_name` overrides the conventional rasch_<id>.csv download filename;
  # the CSV content is always the full table from `fun` (never the curated
  # on-screen display)
  register_table <- function(id, fun, dt_fun, code = NULL, csv_name = NULL) {
    output[[id]] <- renderDT(dt_fun())
    if (!is.null(code)) register_code(id, code)
    output[[paste0(id, "_csv")]] <- downloadHandler(
      filename = function() csv_name %||% paste0("rasch_", id, ".csv"),
      content = function(file) write.csv(fun(), file, row.names = FALSE))
  }
  # the curated stat boxes (test of fit, targeting, BTL test of fit) share
  # one registration: `ui_fun` builds the on-screen label-value rows, while
  # the CSV chip always downloads the COMPLETE table from `csv_fun`
  register_stat_box <- function(id, csv_fun, csv_name, ui_fun, code = NULL) {
    output[[id]] <- renderUI(ui_fun())
    if (!is.null(code)) register_code(id, code)
    output[[paste0(id, "_csv")]] <- downloadHandler(
      filename = function() csv_name,
      content = function(file) write.csv(csv_fun(), file, row.names = FALSE))
  }
  # paired PDF/ZIP batch downloads (all-items explorer plots, all-persons
  # kidmaps): one handler per extension, same content function for both
  # (the writer picks its format from the download file's extension)
  register_batch_download <- function(prefix, base, content) {
    for (ext in c("pdf", "zip")) local({
      ext_ <- ext
      output[[paste0(prefix, "_", ext_)]] <- downloadHandler(
        filename = function() paste0(base(), ".", ext_),
        content = content)
    })
  }
  # APA-leaning DT wrapper: Bootstrap 5 skin, right-aligned numerics, paging
  # controls only when the table needs them. Any fit_resid / infit_ms /
  # outfit_ms column is auto-flagged in the theme danger colour; `p_bold`
  # bolds p-values < .05.
  # curated display columns: fit objects carry every statistic, but the
  # tables show a readable core; the per-table "detailed columns" switch
  # reveals the rest (CSV downloads always contain everything)
  CORE <- list(
    items = c("item", "location", "se", "fit_resid", "infit_ms", "outfit_ms",
              "chisq", "df", "p_adj"),
    person = c("id", "raw", "max_raw", "theta", "se", "extreme", "fit_resid"),
    dif_fact = c("item", "term", "F_uniform", "p_uniform_adj", "eta2_uniform",
                 "F_nonuniform", "p_nonuniform_adj", "eta2_nonuniform"),
    facet = c("level", "severity", "se", "n", "fit_resid"),
    btl_obj = c("object", "location", "se", "comparisons", "wins", "score",
                "infit_ms", "outfit_ms", "fit_resid"),
    btl_judge = c("judge", "n", "infit_ms", "outfit_ms", "fit_resid", "misfit"),
    equate = c("item", "location_1", "location_2", "adj_difference", "t",
               "p_adj", "drift"),
    contrast = c("item", "contrast", "estimate", "se", "statistic", "p_adj"),
    frames = c("set", "group", "rho", "se_log_rho", "origin", "fit_resid",
               "n_responses"),
    compare = c("label", "model", "persons", "items", "two_delta_ll",
                "chisq_per_df", "item_fit_sd", "person_fit_sd", "PSI", "alpha"),
    rescore = c("item", "option", "keyed", "n", "prop", "mean_location",
                "z_sep", "proposed"))
  curate <- function(d, which, full = FALSE, extra = NULL) {
    if (isTRUE(full)) return(d)
    keep <- c(CORE[[which]], extra)
    d[, intersect(keep, names(d)), drop = FALSE]
  }
  # display headers for the tables; downloads always keep the raw names
  DISPLAY_NAMES <- c(
    fit_resid = "Fit resid", fit_resid_pooled = "Pooled fit resid",
    natural_resid = "Natural resid", infit_ms = "Infit MS",
    outfit_ms = "Outfit MS", infit_z = "Infit z", outfit_z = "Outfit z",
    se = "SE", theta = "Location", max_raw = "Max score", raw = "Raw score",
    n_items = "Items", chisq = "Chi-sq", df_fit = "Fit df", p = "p",
    p_adj = "Adj. p", p_bonf = "Bonf. p", p_anova = "ANOVA p",
    F_anova = "ANOVA F", F_uniform = "Uniform F",
    F_nonuniform = "Non-uniform F", p_uniform = "Uniform p",
    p_nonuniform = "Non-uniform p", p_uniform_adj = "Uniform adj. p",
    p_nonuniform_adj = "Non-uniform adj. p",
    eta2_uniform = "Uniform η²",
    eta2_nonuniform = "Non-uniform η²",
    eta2_partial = "Partial η²",
    uniform_DIF = "Uniform DIF", nonuniform_DIF = "Non-uniform DIF",
    superseded = "Superseded", sum_sq = "Sum Sq", mean_sq = "Mean Sq",
    F_value = "F",
    mean_location = "Mean location", point_biserial = "Point-biserial",
    se_location = "SE", z_sep = "Separation z",
    alpha_drop = "α if deleted", item_total = "Item-total r",
    item_rest = "Item-rest r", di = "Discrimination", cum_pct = "Cum. %",
    exp_prop = "Expected", obs_prop = "Observed", obs_mean = "Observed mean",
    exp_mean = "Expected mean",
    exp_value = "Expected value", theta_mean = "Mean location",
    theta_max = "Max location", chisq_per_df = "Chi-sq/df",
    two_delta_ll = "2Δ log-lik", se_log_phi = "SE (log φ)",
    se_log_alpha = "SE (log α)", se_log_rho = "SE (log ρ)",
    mu = "Origin", comparisons = "Comparisons",
    obs_p = "Observed", est_p = "Expected", obs_t = "Threshold prop.",
    item_a = "Item A", item_b = "Item B", q3 = "Q3", q3_star = "Q3*",
    flagged = "Flagged", estimate = "Estimate (logits)",
    statistic = "Statistic", significant = "Significant",
    practical = "Practical", effect = "Effect",
    n_informative = "Informative comparisons")
  # p-value columns render as "<0.001" / 3 dp on the client, so sorting
  # still uses the raw value; detection runs on the ORIGINAL column names
  P_COL_RE <- "^p$|^p_|_p$|^prob$|p_tukey|p_anova|p_adj|p_bonf|p_uniform|p_nonuniform"
  P_RENDER <- DT::JS("function(data,type,row){ if(type==='display'){ if(data===null||data==='') return ''; var x=Number(data); return x<0.001 ? '&lt;0.001' : x.toFixed(3);} return data; }")
  # fit flags, consistent across every model table: a fit residual beyond
  # |2.5|, an outfit mean square outside 0.7-1.3, and an infit mean square
  # outside the tighter 0.8-1.2 (infit is information-weighted, so it varies
  # less under fit; conventional working bands). `orig` is the ORIGINAL
  # column-name vector of the displayed frame (positions match the widget).
  flag_fit_cols <- function(dt, orig) {
    redc <- "var(--bs-danger)"
    for (fc in which(orig == "fit_resid"))
      dt <- formatStyle(dt, fc, color = styleInterval(
        c(-2.5, 2.5), c(redc, "inherit", redc)))
    for (oc in which(orig == "outfit_ms"))
      dt <- formatStyle(dt, oc, color = styleInterval(
        c(0.7, 1.3), c(redc, "inherit", redc)))
    for (ic in which(orig == "infit_ms"))
      dt <- formatStyle(dt, ic, color = styleInterval(
        c(0.8, 1.2), c(redc, "inherit", redc)))
    dt
  }
  num_dt <- function(d, digits = 3, p_bold = NULL,
                     page_len = 15, paging = NULL, ...) {
    orig <- names(d)
    # unname: which() on named logicals yields named position vectors, and
    # jsonlite warns whenever one reaches the widget payload
    num <- unname(vapply(d, is.numeric, TRUE))
    # integer-valued columns (counts, whole-number df) show no decimals;
    # fractional df columns fail the test and keep the 3-dp rounding
    intcol <- unname(vapply(d, function(v)
      is.numeric(v) && all(is.na(v) | v == round(v)), TRUE))
    pcol <- num & grepl(P_COL_RE, orig)
    # pager and info line appear only when the table overflows one page;
    # `paging` can force either way per table
    if (is.null(paging)) paging <- nrow(d) > page_len
    opts <- list(pageLength = if (paging) page_len else max(nrow(d), 1L),
                 scrollX = TRUE,
                 dom = if (paging) "tip" else "t")
    cdefs <- list()
    if (any(num))
      cdefs[[length(cdefs) + 1L]] <- list(className = "dt-right",
                                          targets = which(num) - 1L)
    for (j in which(pcol))
      cdefs[[length(cdefs) + 1L]] <- list(targets = j - 1L, render = P_RENDER)
    if (length(cdefs)) opts$columnDefs <- cdefs
    # display-only renaming; formatting targets are column positions, so
    # they stay tied to the original names computed above
    hit <- orig %in% names(DISPLAY_NAMES)
    # unname: a named names-vector rides into the widget payload and makes
    # jsonlite warn on every table render
    names(d)[hit] <- unname(DISPLAY_NAMES[orig[hit]])
    dt <- datatable(d, rownames = FALSE, style = "bootstrap5",
                    class = "table-sm compact hover order-column", ...,
                    options = opts)
    rnd <- which(num & !intcol & !pcol)
    if (length(rnd)) dt <- formatRound(dt, rnd, digits)
    dt <- flag_fit_cols(dt, orig)
    for (pc in which(orig %in% p_bold))
      dt <- formatStyle(dt, pc,
                        fontWeight = styleInterval(0.05, c("bold", "normal")))
    dt
  }
  # Lower-triangular display of a symmetric matrix: the redundant strictly-upper
  # cells are blanked so each item pair is read once, the first column carries
  # the item name, and NA cells (e.g. the Q3* diagonal) show empty. Values are
  # pre-formatted, so the columns are right-aligned by column definition rather
  # than by numeric class.
  tri_dt <- function(M, digits = 2, flagged = NULL) {
    # tiny magnitudes (including negative zero) print as a clean 0.00
    M[!is.na(M) & abs(M) < 0.5 * 10^(-digits)] <- 0
    disp <- formatC(M, format = "f", digits = digits)
    dim(disp) <- dim(M); dimnames(disp) <- dimnames(M)
    # flagged pairs (Q3* above the threshold) are shown red in place of the
    # heatmap's flag mark; wrap before the triangle is blanked so only kept
    # cells carry the span
    if (!is.null(flagged)) {
      hot <- flagged & !is.na(M) & lower.tri(M, diag = FALSE)
      hot[is.na(hot)] <- FALSE
      disp[hot] <- sprintf(
        '<span style="color:var(--bs-danger);font-weight:600">%s</span>',
        disp[hot])
    }
    disp[upper.tri(M)] <- ""                 # redundant upper triangle
    disp[is.na(M)] <- ""                     # empty diagonal / missing pairs
    df <- data.frame(item = rownames(M), disp, check.names = FALSE,
                     stringsAsFactors = FALSE)
    # escape only the item column (user-controlled labels); the numeric cells
    # intentionally carry the red-flag <span>, so they must not be escaped
    dt <- datatable(df, rownames = FALSE, escape = "item", style = "bootstrap5",
              class = "table hover order-column",
              options = list(paging = FALSE, ordering = FALSE,
                             scrollX = TRUE, dom = "t",
                             columnDefs = list(list(className = "dt-right",
                               targets = seq_len(ncol(M))))))
    # a roomier grid than the compact default, so the matrix reads at a weight
    # closer to its heatmap alongside it
    formatStyle(dt, names(df), fontSize = "1rem",
                paddingTop = "7px", paddingBottom = "7px")
  }
  # Red-highlight a triggering value in place of a boolean flag column (the
  # items table pattern): colour the cell with the theme danger colour when
  # it crosses the threshold, leaving in-range values at the inherited text
  # colour. `d` is the DISPLAYED data frame, so column positions match the
  # rendered table; each helper is a no-op when the column is absent.
  DANGER <- "var(--bs-danger)"
  style_lo_red <- function(dt, d, col, cut)     # red when value < cut (e.g. p)
    Reduce(function(x, j) formatStyle(x, j,
      color = styleInterval(cut, c(DANGER, "inherit"))),
      which(names(d) == col), dt)
  style_mag_red <- function(dt, d, col, mag)    # red when |value| >= mag
    Reduce(function(x, j) formatStyle(x, j,
      color = styleInterval(c(-mag, mag), c(DANGER, "inherit", DANGER))),
      which(names(d) == col), dt)

  # ------------------------------------------------- navbar status chips --
  # compact badges once a fit exists: model, N persons, N items, PSI
  # (objects/comparisons/OSI for a paired-comparison fit); nothing before
  # the first run
  output$nav_status <- renderUI({
    chip <- function(txt, kind = "secondary")
      span(class = paste0("badge text-bg-", kind), txt)
    b <- btl_fit()
    if (!is.null(b)) {
      osi <- b$osi$PSI
      return(div(class = "nav-status d-flex align-items-center gap-1 px-2",
        chip("BTL", "primary"),
        chip(paste(nrow(b$objects), "objects")),
        chip(sprintf("%.0f comparisons", b$n_comparisons)),
        chip(if (finite1(osi)) sprintf("OSI %.2f", osi) else "OSI —",
             if (!finite1(osi)) "secondary"
             else if (osi >= 0.7) "success" else "warning")))
    }
    f <- fit_or_null()
    if (is.null(f)) return(NULL)
    psi <- f$psi$PSI
    # the chip states what was actually fitted: PCM/RSM only make sense for
    # polytomous items, so an all-dichotomous fit reads "Dichotomous"
    model_lab <- if (inherits(f, "rasch_mfrm")) "MFRM"
      else if (inherits(f, "rasch_efrm")) "EFRM"
      else if (max(f$m) == 1L) "Dichotomous"
      else f$model
    div(class = "nav-status d-flex align-items-center gap-1 px-2",
      chip(model_lab, "primary"),
      chip(paste(nrow(f$X), "persons")),
      chip(paste(ncol(f$X), "items")),
      chip(if (finite1(psi)) sprintf("PSI %.2f", psi) else "PSI —",
           if (!finite1(psi)) "secondary"
           else if (psi >= 0.7) "success" else "warning"))
  })

  # -------------------------------------------------------------- summary --
  output$vboxes <- renderUI({
    f <- fit()
    layout_column_wrap(width = "185px", fill = FALSE, class = "mb-3",
      value_box("Persons", nrow(f$X), showcase = glyph("distribution"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Items", ncol(f$X), showcase = glyph("ruler"),
                showcase_layout = "left center", theme = "primary"),
      value_box(span("PSI",
                     tooltip(bs_icon("info-circle", class = "ms-1"),
                             "Person Separation Index: the proportion of variance in person estimates not attributable to measurement error; 0.7 is a conventional minimum for distinguishing groups of persons (Andrich & Marais 2019).")),
                if (finite1(f$psi$PSI)) sprintf("%.3f", f$psi$PSI) else "—",
                showcase = glyph("separation"),
                showcase_layout = "left center",
                theme = if (!finite1(f$psi$PSI)) "secondary"
                        else if (f$psi$PSI >= 0.7) "success" else "warning",
                p(class = "small mb-0",
                  if (finite1(f$psi_noext$PSI))
                    sprintf("%.3f no extremes", f$psi_noext$PSI)
                  else "— no extremes")),
      value_box("Alpha",
                if (finite1(f$alpha$alpha)) sprintf("%.3f", f$alpha$alpha)
                else "—",
                showcase = glyph("alpha"),
                showcase_layout = "left center",
                theme = if (!finite1(f$alpha$alpha)) "secondary"
                        else if (f$alpha$alpha >= 0.7) "success" else "warning",
                p(class = "small mb-0",
                  if (isFALSE(f$alpha$applicable))
                    sprintf("complete cases only (n = %d)", f$alpha$n)
                  else sprintf("n = %d complete", f$alpha$n))),
      value_box("Item-trait p",
                if (finite1(f$total_chisq_p)) fmt_p(f$total_chisq_p) else "—",
                showcase = glyph("chisq"),
                showcase_layout = "left center",
                # neutral even when significant: red is reserved for
                # in-table cell highlighting
                theme = if (finite1(f$total_chisq_p) &&
                            f$total_chisq_p >= 0.05) "success"
                        else "secondary"),
      value_box("Power of fit", f$power_of_fit,
                showcase = glyph("power"),
                showcase_layout = "left center", theme = "secondary")
    )
  })

  # test-of-fit and targeting/reliability summaries as curated stat boxes
  # (values read off the fit object directly); the CSV chips download the
  # COMPLETE fit_summary_table / targeting_table with raw column names.
  # The score-to-measure table stays available via score_table(fit) and the
  # everything-ZIP; thresholds live on the Items explorer Thresholds tab.
  register_stat_box("fitsum_tbl",
    csv_fun = function() fit_summary_table(fit()),
    csv_name = "fit_summary.csv",
    ui_fun = function() {
      f <- fit()
      est <- f$est %||% f
      dis <- names(which(vapply(f$thresholds_diag, function(d)
        !d$ordered && length(d$thresholds) > 1, TRUE)))
      conv <- if (isTRUE(est$converged))
        sprintf("converged in %d iterations", est$iterations)
      else span(class = "text-danger",
                sprintf("did not converge in %d iterations", est$iterations))
      tagList(
        div(class = "stat-head",
            f$model, " · pairwise conditional ML · ", conv),
        stat_rows(
          stat_row("Item-trait chi-square",
                   sprintf("%.2f on %d df, %s", f$total_chisq, f$total_df,
                           p_lab(f$total_chisq_p))),
          stat_row("Item fit residual",
                   sprintf("mean %.2f, SD %.2f", f$item_fit_summary$mean,
                           f$item_fit_summary$sd)),
          stat_row("Person fit residual",
                   sprintf("mean %.2f, SD %.2f", f$person_fit_summary$mean,
                           f$person_fit_summary$sd)),
          stat_row("Item location-residual correlation",
                   { v <- f$summary_stats$cor_item_fit_location
                     if (finite1(v)) sprintf("%.3f", v) else "NA" }),
          stat_row("Person location-residual correlation",
                   { v <- f$summary_stats$cor_person_fit_location
                     if (finite1(v)) sprintf("%.3f", v) else "NA" }),
          stat_row("Items with adj. chi-square p < .05",
                   sprintf("%d of %d",
                           sum(f$items$p_adj < 0.05, na.rm = TRUE),
                           nrow(f$items))),
          stat_row("Disordered thresholds",
                   if (length(dis)) paste(dis, collapse = ", ") else "none")))
    },
    code = function() "fit_summary_table(fit)")
  # routine handling notes (the old text panel printed fit$notes)
  output$fitsum_notes <- renderUI({
    f <- fit()
    if (!length(f$notes)) return(NULL)
    sprintf("Note. %s.", paste(f$notes, collapse = "; "))
  })
  register_stat_box("targeting_tbl",
    csv_fun = function() targeting_table(fit()),
    csv_name = "targeting.csv",
    ui_fun = function() {
      f <- fit(); ss <- f$summary_stats; tg <- f$targeting
      psi_txt <- if (!finite1(f$psi$PSI)) "—"
      else if (finite1(f$psi_noext$PSI))
        sprintf("%.3f (%.3f without extremes)", f$psi$PSI, f$psi_noext$PSI)
      else sprintf("%.3f", f$psi$PSI)
      alpha_txt <- if (finite1(f$alpha$alpha)) sprintf("%.3f", f$alpha$alpha)
      else sprintf("not applicable (%d complete cases)", f$alpha$n %||% 0L)
      stat_rows(
        stat_row("Person location",
                 sprintf("mean %.2f, SD %.2f", ss$person_location$mean,
                         ss$person_location$sd)),
        stat_row("Item location",
                 sprintf("SD %.2f (mean constrained to 0)",
                         ss$item_location$sd)),
        stat_row("Threshold range",
                 sprintf("%.2f to %.2f", tg$threshold_range[1],
                         tg$threshold_range[2])),
        stat_row("Persons beyond thresholds",
                 sprintf("%.1f%% below · %.1f%% above",
                         100 * tg$prop_below, 100 * tg$prop_above)),
        stat_row("PSI", psi_txt),
        stat_row("Item reliability",
                 { v <- f$isi$PSI
                   if (finite1(v)) sprintf("%.3f", v) else "NA" }),
        stat_row("Person strata",
                 { v <- f$psi$strata
                   if (finite1(v)) sprintf("%.1f", v) else "NA" }),
        stat_row("Coefficient alpha", alpha_txt))
    },
    code = function() "targeting_table(fit)")

  # traditional (CTT) statistics: shown on the Items page (last accordion
  # panel), with the header line, CSV, and code footer registered here
  ctt_res <- reactive(tryCatch(ctt_table(fit()), error = function(e) e))
  register_code("ctt_tbl", function() "ctt_table(fit)$table")
  output$ctt_head <- renderUI({
    ct <- ctt_res()
    if (inherits(ct, "error"))
      return(p(class = "text-muted",
               paste("Traditional statistics unavailable:", conditionMessage(ct))))
    # total-score summaries need complete responders; with structural
    # missingness (linked forms) there may be none, so the header falls back
    # to the available-case framing instead of printing NA values
    if (!finite1(ct$mean))
      return(p(class = "small mb-2", HTML(sprintf(
        "Per-item statistics use available cases (item n %d&ndash;%d). Too few complete responders (n = %d) for the total-score mean, SD, SEM%s.",
        ct$n_range[1], ct$n_range[2], ct$n,
        if (finite1(ct$alpha)) sprintf("; alpha <b>%.3f</b>", ct$alpha)
        else ", and alpha"))))
    p(class = "small mb-2", HTML(sprintf(
      "Raw score mean <b>%.2f</b>, SD <b>%.2f</b>; alpha <b>%.3f</b>; SEM <b>%.2f</b> (one value for all persons) &mdash; complete cases n = %d.",
      ct$mean, ct$sd, ct$alpha, ct$sem, ct$n)))
  })
  output$ctt_tbl <- renderDT({
    ct <- ctt_res()
    validate(need(!inherits(ct, "error"),
                  "No traditional statistics for this fit."))
    num_dt(ct$table)
  })
  output$ctt_tbl_csv <- downloadHandler(
    filename = function() "rasch_ctt_tbl.csv",
    content = function(file) {
      ct <- ctt_res(); req(!inherits(ct, "error"))
      write.csv(ct$table, file, row.names = FALSE)
    })

  # likelihood-ratio test of PCM against the rating parameterisation; only
  # meaningful for a PCM fit whose items share a common maximum score > 1.
  # The bottom row is server-built so the card hides when it does not apply.
  lr_res <- reactiveVal(NULL)
  output$summary_bottom <- renderUI({
    f <- fit()
    lr_applies <- identical(f$model, "PCM") && length(unique(f$m)) == 1L &&
      max(f$m) >= 2L
    if (!lr_applies) return(NULL)
    layout_columns(col_widths = breakpoints(sm = 12, xl = 6), .lr_card())
  })
  register_code("lr", function() "lr_test(fit)")
  observeEvent(input$run_lr, {
    f <- fit()
    r <- withProgress(message = "Refitting with the rating parameterisation…",
                      value = 0.4,
                      tryCatch(lr_test(f), error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("LR test failed:", conditionMessage(r)),
                       type = "error", duration = 10)
    else lr_res(r)
  })
  output$lr_txt <- renderPrint({
    r <- lr_res()
    validate(need(!is.null(r), "Press the button to run the test."))
    cat(sprintf("Raw composite chi-square %.3f on %d df, p = %s (conventional display; anticonservative)\n",
                r$chisq, r$df, fmt_p(r$p)))
    if (is.finite(r$chisq_adj))
      cat(sprintf("Adjusted chi-square %.3f, p = %s (Kent 1982 calibration)\n",
                  r$chisq_adj, fmt_p(r$p_adj)))
    else cat("Adjusted (Kent 1982) statistic unavailable for this fit.\n")
    cat(sprintf("log-likelihood (pairwise composite): PCM %.3f, rating %.3f\n",
                r$loglik_pcm, r$loglik_rsm))
  })
  # ---------------------------------------------------------------- items --
  output$items_vboxes <- renderUI({
    f <- fit()
    mis <- sum(f$items$p_adj < 0.05, na.rm = TRUE)
    dis <- sum(vapply(f$thresholds_diag, function(d)
      !d$ordered && length(d$thresholds) > 1, TRUE))
    layout_column_wrap(width = "200px", fill = FALSE, class = "mb-3",
      value_box("Items", nrow(f$items), showcase = glyph("ruler"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Adj. chi-square p < .05", mis,
                showcase = glyph("chisq"),
                showcase_layout = "left center",
                theme = if (mis > 0) "secondary" else "success"),
      value_box("Disordered thresholds", dis,
                showcase = glyph("disorder"),
                showcase_layout = "left center",
                theme = if (dis > 0) "warning" else "success"))
  })
  output$items_note <- renderUI({
    f <- fit(); d <- f$items
    dis <- names(which(vapply(f$thresholds_diag, function(x)
      !x$ordered && length(x$thresholds) > 1, TRUE)))
    sprintf("Note. %d of %d items beyond |fit residual| 2.5; %d with adjusted chi-square p < .05; disordered thresholds: %s.",
            sum(abs(d$fit_resid) > 2.5, na.rm = TRUE), nrow(d),
            sum(d$p_adj < 0.05, na.rm = TRUE),
            if (length(dis)) paste(dis, collapse = ", ") else "none")
  })
  # any chi-square sample-size adjustment is applied inside the fit (the
  # adjust_N run control), so the table shows the fit's own statistics
  register_table("items_tbl", function() fit()$items, function() {
    d <- curate(fit()$items, "items", full = isTRUE(input$items_full))
    dt <- num_dt(d, page_len = 25, selection = "single",
                 p_bold = c("p_adj", "p_anova"))
    # fit residual, infit and outfit are flagged by num_dt; the adjusted
    # chi-square probability turns red at .05 as well (no single flag column)
    for (j in which(names(d) == "p_adj"))
      dt <- formatStyle(dt, j, color = styleInterval(
        0.05, c("var(--bs-danger)", "inherit")))
    dt
  }, code = function() "fit$items")

  # per-class-interval breakdown of the selected item's chi-square
  chisq_res <- reactive(chisq_detail(fit(), sel_item()))
  output$chisq_caption <- renderUI({
    cd <- chisq_res()
    p(class = "small mb-2", HTML(sprintf(
      "<b>%s</b> (location %.3f): total chi-square <b>%.3f</b> on %d df, p = %s; whole-sample mean = %.3f. Intervals with fewer than 2 responders carry no chi-square contribution.",
      cd$item, cd$location, cd$chisq, cd$df, fmt_p(cd$p), cd$ave)))
  })
  output$chisq_int_tbl <- renderDT({
    d <- chisq_res()$intervals
    d$Excluded <- ifelse(d$used, "", "*")
    d$used <- NULL
    d$theta_max <- NULL   # in the CSV download, not the default display
    num_dt(d)
  })
  output$chisq_cat_tbl <- renderDT(num_dt(chisq_res()$categories))
  output$chisq_int_csv <- downloadHandler(
    filename = function() paste0("rasch_chisq_intervals_", sel_item(), ".csv"),
    content = function(file)
      write.csv(chisq_res()$intervals, file, row.names = FALSE))
  output$chisq_cat_csv <- downloadHandler(
    filename = function() paste0("rasch_chisq_categories_", sel_item(), ".csv"),
    content = function(file)
      write.csv(chisq_res()$categories, file, row.names = FALSE))
  register_code("chisq", function()
    sprintf('chisq_detail(fit, "%s")', sel_item()))

  # principal-components estimates (only for pc_components fits)
  output$pc_comp_ui <- renderUI({
    if (is.null(fit()$est$components)) return(NULL)
    tableCard("pc_tbl", "Principal components (location/spread/skewness/kurtosis)",
              "Andrich principal-components threshold estimates with standard errors; NA where an item's number of thresholds does not support the component.")
  })
  register_table("pc_tbl", function() fit()$est$components, function() {
    validate(need(!is.null(fit()$est$components),
                  "Run with principal-components threshold estimation to see the components."))
    num_dt(fit()$est$components)
  }, code = function() "fit$est$components")

  # explorer display settings (inline on the tab-strip controls row):
  # class intervals and scale
  # range, resolved with fallbacks; the code footers add n_groups / grid
  # only when they differ from the defaults, keeping default snippets minimal
  ex_ng <- reactive({
    ng <- input$ex_ng
    if (is.null(ng) || is.na(ng)) fit()$n_groups else as.integer(ng)
  })
  ex_rng <- reactive({
    r <- input$ex_rng
    if (is.null(r) || length(r) != 2L || anyNA(r)) c(-5, 5) else as.numeric(r)
  })
  ex_grid <- reactive(seq(ex_rng()[1], ex_rng()[2], 0.05))
  ex_code_args <- reactive(paste0(c(
    if (ex_ng() != fit()$n_groups) sprintf(", n_groups = %d", ex_ng()),
    if (!isTRUE(all.equal(ex_rng(), c(-5, 5))))
      sprintf(", grid = seq(%g, %g, 0.05)", ex_rng()[1], ex_rng()[2])),
    collapse = ""))
  # second code-footer line pointing at the matching all-items batch export
  ex_batch_line <- function(what)
    sprintf('\n# all items: save_item_plots(fit, "%s", "%s_all_items.pdf")',
            what, what)
  register_plot("icc",  function()
    plot_icc(fit(), sel_item(), n_groups = ex_ng(), grid = ex_grid()),
    code = function() paste0(sprintf('plot_icc(fit, "%s"%s)',
                                     sel_item(), ex_code_args()),
                             ex_batch_line("icc")))
  register_plot("ccc",  function()
    plot_ccc(fit(), sel_item(), observed = isTRUE(input$show_obs),
             n_groups = ex_ng(), grid = ex_grid()),
    code = function() paste0(sprintf('plot_ccc(fit, "%s", observed = %s%s)',
                                     sel_item(), isTRUE(input$show_obs),
                                     ex_code_args()),
                             ex_batch_line("ccc")))
  register_plot("tpc",  function()
    plot_threshold_prob(fit(), sel_item(), observed = isTRUE(input$show_obs),
                        n_groups = ex_ng(), grid = ex_grid()),
    code = function() paste0(
      sprintf('plot_threshold_prob(fit, "%s", observed = %s%s)',
              sel_item(), isTRUE(input$show_obs), ex_code_args()),
      ex_batch_line("tpc")))
  register_plot("cfreq", function() plot_catfreq(fit(), sel_item()),
                code = function() paste0(
                  sprintf('plot_catfreq(fit, "%s")', sel_item()),
                  ex_batch_line("cfreq")))
  # all-items batch downloads for the active explorer tab, honouring the
  # class-interval, range, and observed-points controls
  ex_what <- reactive(switch(input$items_nav %||% "ICC",
    ICC = "icc", Categories = "ccc", Thresholds = "tpc",
    Frequencies = "cfreq", "icc"))
  register_batch_download("items_all",
    base = function() paste0(ex_what(), "_all_items"),
    content = function(file)
      withProgress(message = "Drawing every item…", value = 0.4,
        save_item_plots(fit(), ex_what(), file, n_groups = ex_ng(),
                        grid = ex_grid(),
                        observed = isTRUE(input$show_obs))))
  mc_dat <- reactive({
    f <- fit()
    validate(need(!is.null(f$mc),
                  "Provide a multiple-choice key (CSV: item,key) to see distractor analysis."))
    distractor_analysis(f)
  })
  register_table("distractor_tbl", function() mc_dat(), function() {
    d <- mc_dat()
    d$keyed <- ifelse(d$keyed, "*", "")
    d$flag <- ifelse(d$flag, "MISKEY?", "")
    num_dt(d)
  }, code = function() "distractor_analysis(fit)")
  # the plotted item: the selected one if it is multiple-choice, else the
  # first MC item (the code disclosure mirrors this resolution)
  distractor_item <- reactive({
    f <- fit()
    req(!is.null(f$mc))
    if (sel_item() %in% colnames(f$mc$raw)) sel_item() else
      colnames(f$mc$raw)[1]
  })
  register_plot("distractor_plot", function() {
    f <- fit()
    validate(need(!is.null(f$mc),
                  "Provide a multiple-choice key (CSV: item,key) to see option curves."))
    plot_distractors(f, distractor_item())
  }, code = function()
    sprintf('plot_distractors(fit, "%s")', distractor_item()))

  # polytomous option-scoring proposal (Andrich & Styles 2011)
  rescore_res <- reactiveVal(NULL)
  observeEvent(input$rescore_go, {
    f <- fit()
    if (is.null(f$mc)) {
      showNotification("Provide a multiple-choice key first.", type = "warning")
      return()
    }
    res <- tryCatch(distractor_rescore(f,
                                       min_n = max(2, input$rescore_min_n %||% 20),
                                       z = max(0.1, input$rescore_z %||% 1.96)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Rescore proposal failed:", conditionMessage(res)),
                       type = "error")
      return()
    }
    rescore_res(res)
    n_cred <- sum(res$option_scores$score > 0 &
                  !res$evidence$keyed[match(paste(res$option_scores$item,
                                                  res$option_scores$option),
                                            paste(res$evidence$item,
                                                  res$evidence$option))])
    showNotification(sprintf("%d distractor(s) proposed for partial credit.",
                             n_cred), type = "message")
  })
  output$rescore_tbl <- DT::renderDT({
    res <- rescore_res()
    validate(need(!is.null(res), "Run the proposal to see the evidence table."))
    d <- curate(res$evidence, "rescore", full = isTRUE(input$rescore_full))
    if ("keyed" %in% names(d)) d$keyed <- ifelse(d$keyed, "*", "")
    num_dt(d)
  })
  register_code("rescore_tbl", function()
    sprintf('distractor_rescore(fit, min_n = %s, z = %s)$option_scores',
            max(2, input$rescore_min_n %||% 20),
            max(0.1, input$rescore_z %||% 1.96)))
  output$dl_rescore <- downloadHandler(
    filename = function() "option_scores.csv",
    content = function(file) {
      res <- rescore_res()
      if (is.null(res)) stop("run the proposal first")
      write.csv(res$option_scores, file, row.names = FALSE)
    })

  # -------------------------------------------------------------- persons --
  output$persons_vboxes <- renderUI({
    f <- fit(); d <- f$person
    mis <- sum(abs(d$fit_resid) > 2.5, na.rm = TRUE)
    layout_column_wrap(width = "200px", fill = FALSE, class = "mb-3",
      value_box("Persons", nrow(d), showcase = glyph("distribution"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Extreme scores", sum(d$extreme, na.rm = TRUE),
                showcase = glyph("range"),
                showcase_layout = "left center", theme = "secondary"),
      value_box("Misfitting persons", mis,
                showcase = glyph("outlier"),
                showcase_layout = "left center",
                theme = if (mis > 0) "secondary" else "success"))
  })
  register_table("person_tbl", function() fit()$person, function() {
    fac <- names(fit()$factors)
    d <- curate(fit()$person, "person", full = isTRUE(input$persons_full),
                extra = fac)
    dt <- datatable(d, rownames = FALSE, filter = "top", selection = "single",
                    style = "bootstrap5",
                    class = "table-sm compact hover order-column",
                    options = list(pageLength = 15, scrollX = TRUE, dom = "tip")) |>
      formatRound(names(d)[vapply(d, is.numeric, TRUE) &
                           !names(d) %in% c("raw", "max_raw", "n_items", "class_interval")], 3)
    # the same fit flags as every other table (fit residual, and the infit /
    # outfit mean squares revealed by the detailed-columns switch)
    dt <- flag_fit_cols(dt, names(d))
    dt
  }, code = function() "fit$person")
  # the person drawn on the kidmap: the selected table row, defaulting to 1
  sel_person <- reactive({
    n <- input$person_tbl_rows_selected
    if (length(n)) n[1] else 1L
  })
  # confidence level for the kidmap band (header select; default 95%);
  # the code footers mention `level` only when it differs from the default
  kid_level <- reactive({
    lv <- suppressWarnings(as.numeric(input$kid_level %||% "0.95"))
    if (!isTRUE(is.finite(lv)) || lv <= 0 || lv >= 1) 0.95 else lv
  })
  kid_level_arg <- reactive(
    if (isTRUE(all.equal(kid_level(), 0.95))) ""
    else sprintf(", level = %g", kid_level()))
  register_plot("kidmap", function()
    plot_kidmap(fit(), person = sel_person(), level = kid_level()),
    code = function() paste0(
      sprintf('plot_kidmap(fit, person = %d%s)', sel_person(),
              kid_level_arg()),
      sprintf('\n# all persons: save_person_plots(fit, "kidmaps.pdf"%s)',
              kid_level_arg())))
  # batch kidmaps: multi-page PDF or ZIP of PNGs, chosen by the extension
  # (the download temp file carries the filename's extension)
  register_batch_download("kidmap_all",
    base = function() "kidmaps",
    content = function(file)
      withProgress(message = "Drawing a kidmap for every person…",
                   value = 0.4,
                   save_person_plots(fit(), file, level = kid_level())))
  register_plot("rdist_p", function() plot_resid_dist(fit(), "persons"),
                code = function() 'plot_resid_dist(fit, "persons")')
  register_plot("pfit",  function() plot_person_fit(fit()),
                code = function() "plot_person_fit(fit)")

  # ------------------------------------------------------- targeting plots --
  # sidebar-controlled bins and scale range shared by both targeting plots
  tg_bins <- reactive({
    b <- input$tg_bins
    if (is.null(b) || is.na(b)) 35L else as.integer(b)
  })
  tg_rng <- reactive({
    r <- input$tg_rng
    if (is.null(r) || length(r) != 2L || anyNA(r)) NULL else as.numeric(r)
  })
  tg_code <- function(fun) function() {
    r <- tg_rng() %||% c(-5, 5)
    sprintf("%s(fit, bins = %d, xlim = c(%g, %g))", fun, tg_bins(), r[1], r[2])
  }
  register_plot("pim_p", function()
    plot_pimap(fit(), bins = tg_bins(), xlim = tg_rng()),
    code = tg_code("plot_pimap"))
  register_plot("wright", function()
    plot_wright(fit(), bins = tg_bins(), xlim = tg_rng()), h = 7.5,
    code = tg_code("plot_wright"))

  # ---------------------------------------------- test & item map plots --
  # thrmap and imap render on the Items page; tcc/tif/guttman on Test
  register_plot("thrmap", function() plot_threshold_map(fit()), h = 7,
                code = function() "plot_threshold_map(fit)")
  register_plot("imap",   function() plot_item_map(fit()),
                code = function() "plot_item_map(fit)")
  register_plot("rdist_i", function() plot_resid_dist(fit(), "items"),
                code = function() 'plot_resid_dist(fit, "items")')
  # Test-page scale range (default -6..6 matches the functions' own default,
  # so the code footers add `grid` only when the slider has been moved)
  ts_rng <- reactive({
    r <- input$ts_rng
    if (is.null(r) || length(r) != 2L || anyNA(r)) c(-6, 6) else as.numeric(r)
  })
  ts_grid <- reactive(seq(ts_rng()[1], ts_rng()[2], 0.05))
  ts_code_arg <- reactive(
    if (isTRUE(all.equal(ts_rng(), c(-6, 6)))) ""
    else sprintf(", grid = seq(%g, %g, 0.05)", ts_rng()[1], ts_rng()[2]))
  register_plot("tcc",    function() plot_tcc(fit(), grid = ts_grid()),
                code = function() paste0("plot_tcc(fit", ts_code_arg(), ")"))
  register_plot("tif",    function() plot_tif(fit(), grid = ts_grid()),
                code = function() paste0("plot_tif(fit", ts_code_arg(), ")"))
  register_plot("guttman", function() plot_guttman(fit()), h = 7,
                code = function() "plot_guttman(fit)")

  # ------------------------------------------------------------------ DIF --
  dif_alpha <- reactive(clamp01(input$dif_alpha, 0.05))
  # the Model toggle (main effects vs interactions) is immaterial with a
  # single nominated factor, where the model is always the one-way ANOVA
  dif_multi <- reactive(!is.null(fit()$factors) &&
                          length(names(fit()$factors)) > 1)
  output$dif_multifactor <- reactive(dif_multi())
  outputOptions(output, "dif_multifactor", suspendWhenHidden = FALSE)
  # one merged reactive: one factor -> one-way ANOVA (one row per item); several
  # factors -> joint model, main effects by default or factor-by-factor
  # interactions when requested (effects is ignored with a single factor)
  dif_res <- reactive({
    f <- fit(); req(!is.null(f$factors))
    dif_anova(f, effects = input$dif_effects %||% "main",
              p_adjust = input$dif_padj %||% "BH", alpha = dif_alpha())
  })
  # code footer: omit the effects argument when there is only one factor
  dif_effects_arg <- function()
    if (dif_multi()) sprintf('effects = "%s", ', input$dif_effects %||% "main")
    else ""
  register_table("dif_tbl", function() dif_res()$summary, function() {
    d <- curate(dif_res()$summary, "dif_fact", full = isTRUE(input$dif_full))
    if ("superseded" %in% names(d))
      d$superseded <- ifelse(d$superseded, "(superseded)", "")
    dt <- num_dt(d, selection = "single")
    # no boolean DIF flag columns: the adjusted probability turns red when it
    # crosses alpha (uniform = factor effect, non-uniform = factor x interval)
    dt <- style_lo_red(dt, d, "p_uniform_adj", dif_alpha())
    style_lo_red(dt, d, "p_nonuniform_adj", dif_alpha())
  }, code = function()
    sprintf('dif_anova(fit, %sp_adjust = "%s", alpha = %s)$summary',
            dif_effects_arg(), input$dif_padj %||% "BH", dif_alpha()))
  # full per-item ANOVA table: every model term, computed lazily when its
  # disclosure is first switched on (the DT renders only once visible)
  register_table("dif_full_tbl", function() dif_res()$terms, function() {
    d <- dif_res()$terms
    d$significant <- NULL                     # red adjusted p replaces the flag
    d$superseded <- ifelse(d$superseded, "(superseded)", "")
    style_lo_red(num_dt(d), d, "p_adj", dif_alpha())
  }, code = function()
    sprintf('dif_anova(fit, %sp_adjust = "%s", alpha = %s)$terms',
            dif_effects_arg(), input$dif_padj %||% "BH", dif_alpha()))
  output$dif_note <- renderUI({
    r <- dif_res(); d <- r$summary
    sig <- sum(d$uniform_DIF | d$nonuniform_DIF, na.rm = TRUE)
    base <- sprintf("Note. %d of %d terms significant after adjustment. Class intervals: %s (from the smallest cell).",
                    sig, nrow(d),
                    if (is.null(r$n_groups)) "NA" else r$n_groups)
    if (identical(r$effects, "factorial")) {
      sup <- sum(d$superseded, na.rm = TRUE)
      if (sup)
        base <- paste0(base,
          sprintf(" %d main-effect term(s) superseded by an interaction.", sup))
    }
    within <- r$within[!is.na(r$within)]
    if (length(within))
      base <- paste0(base,
        sprintf(" Within-subject factor(s) tested by repeated-measures ANOVA: %s.",
                paste(within, collapse = ", ")))
    base
  })
  # the items of the DIF summary in rendered row order (curate only drops
  # columns, so the order is preserved). The selected row drives the group-ICC
  # and the pairwise-comparisons panel; nothing selected defaults to the top row.
  dif_tbl_items <- reactive(dif_res()$summary$item)
  # the item and term of the currently selected summary row (top row by
  # default); the term's factor variables (":"-separated for interactions)
  # feed dif_size and the group-ICC
  dif_sel_row <- reactive({
    r <- input$dif_tbl_rows_selected
    its <- dif_tbl_items()
    if (length(r) && !is.na(r[1]) && r[1] >= 1 && r[1] <= length(its))
      r[1] else 1L
  })
  dif_sel_item <- reactive({
    its <- dif_tbl_items(); req(length(its) >= 1); its[dif_sel_row()]
  })
  dif_sel_term <- reactive({
    tm <- dif_res()$summary$term; req(length(tm) >= 1); tm[dif_sel_row()]
  })
  dif_sel_vars <- reactive(strsplit(dif_sel_term(), ":", fixed = TRUE)[[1]])
  # the group-ICC uses the selected term's factor(s); plot_icc accepts several
  # factor names, so an interaction row overlays the factor-combination cells
  register_plot("dif_icc", function() {
    f <- fit()
    req(dif_sel_item() %in% f$items$item, !is.null(f$factors))
    plot_icc(f, dif_sel_item(), group = dif_sel_vars())
  }, code = function()
    sprintf('plot_icc(fit, "%s", group = %s)', dif_sel_item() %||% "",
            qvec(dif_sel_vars())))

  # Pairwise comparisons for the selected DIF row: the row's item resolved by
  # the row's term (interaction terms resolve by their cells), giving the
  # pairwise location differences in logits (the DIF magnitude) with Holm
  # familywise adjustment. Recomputes on row selection and the two criteria.
  dif_size_res <- reactive({
    f <- fit(); req(!is.null(f$factors))
    vars <- dif_sel_vars(); req(length(vars) >= 1, all(vars %in% names(f$factors)))
    req(dif_sel_item() %in% f$items$item)
    flg <- max(0.05, input$dif_size_flag %||% 0.5)
    mn <- max(2, input$dif_size_minn %||% 20)
    tryCatch(dif_size(f, dif_sel_item(), by = vars,
                      flag_logits = flg, min_n = mn, alpha = dif_alpha()),
             error = function(e) e)
  })
  # the resolved-location summary above the pairwise table; a muted line
  # appears only when dif_size genuinely errors (e.g. too few responders)
  output$dif_levels_note <- renderUI({
    ds <- dif_size_res()
    if (inherits(ds, "error"))
      return(p(class = "text-muted small mb-0", conditionMessage(ds)))
    p(class = "text-muted small mb-2",
      sprintf("%s resolved by %s: %s.", ds$item, ds$by,
              paste(sprintf("%s %.2f", ds$levels$level, ds$levels$location),
                    collapse = ", ")))
  })
  # always the pairwise magnitude table when dif_size succeeds: one row for a
  # two-level factor (the single effect size), more for > 2 levels
  output$dif_size_tbl <- DT::renderDT({
    ds <- dif_size_res()
    req(!inherits(ds, "error"))
    d <- ds$pairs[, c("level_a", "level_b", "difference", "se", "z",
                      "lower", "upper", "p_adj")]
    dt <- style_mag_red(num_dt(d), d, "difference", ds$flag_logits)
    style_lo_red(dt, d, "p_adj", ds$alpha)
  })
  register_code("dif_size_tbl", function() {
    # carry the panel's settings so the snippet reproduces the shown table
    flg <- max(0.05, input$dif_size_flag %||% 0.5)
    mn <- max(2, input$dif_size_minn %||% 20)
    extra <- paste0(
      if (flg != 0.5) sprintf(", flag_logits = %s", flg) else "",
      if (mn != 20) sprintf(", min_n = %s", mn) else "",
      if (dif_alpha() != 0.05) sprintf(", alpha = %s", dif_alpha()) else "")
    sprintf('dif_size(fit, "%s", by = %s%s)', dif_sel_item() %||% "",
            if (length(dif_sel_vars()) == 1L) qstr(dif_sel_vars())
            else qvec(dif_sel_vars()), extra)
  })
  output$dl_dif_size <- downloadHandler(
    filename = function() "dif_sizes.csv",
    content = function(file) {
      ds <- dif_size_res()
      if (inherits(ds, "error")) stop(conditionMessage(ds))
      write.csv(ds$pairs, file, row.names = FALSE)
    })

  # planned DIF contrasts: the family is derived from the factor structure
  # and every question tested at once; the nominated ID column (if any) is
  # reserved for pairing person rows and excluded from the contrast factors
  contr_res <- reactiveVal(NULL)
  observeEvent(input$pc_run, {
    f <- fit()
    idn <- input$pc_id %||% ""
    res <- tryCatch({
      if (is.null(f$factors) || !length(names(f$factors)))
        stop("nominate at least one person factor on the Data page")
      fac <- if (nzchar(idn)) setdiff(names(f$factors), idn)
             else names(f$factors)
      if (!length(fac))
        stop("no factors left to contrast once the ID column is reserved")
      its <- if (length(input$pc_items)) input$pc_items else NULL
      withProgress(message = "Resolving items…",
                   detail = "one refit per item", value = 0.15,
        dif_contrasts(f, factors = fac, items = its,
                      id = if (nzchar(idn)) f$factors[[idn]] else NULL))
    }, error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Planned contrasts:", conditionMessage(res)),
                       type = "warning", duration = 8)
      contr_res(NULL)
    } else contr_res(res)
  })
  # the derived family in words, shown once a run has completed
  output$contr_family <- renderUI({
    r <- contr_res()
    if (is.null(r)) return(NULL)
    tagList(
      h6("The planned family"),
      div(class = "small font-monospace mb-2",
          lapply(seq_len(nrow(r$family)), function(i)
            div(paste0(r$family$contrast[i],
                       if (r$family$within[i]) "  [within subjects]" else "")))),
      if (isTRUE(r$paired))
        p(class = "text-muted small mb-2",
          "Stacked design: tests use person-level residual scores."),
      if (length(r$notes))
        p(class = "text-muted small mb-2", paste(r$notes, collapse = "; ")))
  })
  register_table("contr_tbl", function() {
    r <- contr_res(); req(!is.null(r)); r$table
  }, function() {
    r <- contr_res()
    validate(need(!is.null(r),
                  "Press the button to derive the planned family from the factor structure and test it."))
    d <- curate(r$table, "contrast", full = isTRUE(input$contr_full))
    if ("within" %in% names(d)) d$within <- ifelse(d$within, "*", "")
    # no significant/practical flags: the estimate turns red at the practical
    # criterion, the adjusted p at alpha
    dt <- style_mag_red(num_dt(d), d, "estimate", r$flag_logits %||% 0.5)
    style_lo_red(dt, d, "p_adj", r$alpha %||% 0.05)
  }, code = function() {
    its <- input$pc_items
    sprintf('dif_contrasts(fit%s%s)',
            if (length(its))
              paste0(', items = c("', paste(its, collapse = '", "'), '")')
            else "",
            if (nzchar(input$pc_id %||% ""))
              paste0(', id = "', input$pc_id, '"') else "")
  }, csv_name = "dif_contrasts.csv")

  # ------------------------------------------------------------------- BTL --
  bfit <- reactive({
    validate(need(!is.null(btl_fit()),
                  "Run a paired-comparisons (BTL) analysis from the Data page to see results here."))
    btl_fit()
  })
  output$btl_boxes <- renderUI({
    f <- bfit()
    layout_column_wrap(width = "185px", fill = FALSE, class = "mb-3",
      value_box("Objects", nrow(f$objects), showcase = glyph("podium"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Comparisons", sprintf("%.0f", f$n_comparisons),
                showcase = glyph("pair"),
                showcase_layout = "left center", theme = "primary"),
      if (!is.null(f$judges))
        value_box("Judges", nrow(f$judges),
                  showcase = glyph("balance"),
                  showcase_layout = "left center", theme = "primary"),
      value_box("Object separation",
                if (finite1(f$osi$PSI)) sprintf("%.3f", f$osi$PSI) else "—",
                showcase = glyph("separation"),
                showcase_layout = "left center",
                theme = if (!finite1(f$osi$PSI)) "secondary"
                        else if (f$osi$PSI >= 0.7) "success" else "warning"),
      value_box("Pairwise fit p",
                if (finite1(f$total_p)) fmt_p(f$total_p) else "—",
                showcase = glyph("chisq"),
                showcase_layout = "left center",
                theme = if (finite1(f$total_p) && f$total_p >= 0.05)
                  "success" else "secondary"))
  })
  # test-of-fit stat box (Summary page): the paired-comparison headline set
  # read off the fit; the CSV chip downloads the COMPLETE table from
  # fit_summary_table()'s rasch_btl method
  register_stat_box("btl_fitsum_tbl",
    csv_fun = function() fit_summary_table(bfit()),
    csv_name = "fit_summary.csv",
    ui_fun = function() {
      f <- bfit()
      graded <- !is.null(f$m) && f$m > 1L
      model_lab <- if (graded)
        sprintf("Graded paired comparisons (%d categories)", f$m + 1L)
      else "Paired comparisons (BTL)"
      conv <- if (isTRUE(f$converged))
        sprintf("converged in %d iterations", f$iterations)
      else span(class = "text-danger",
                sprintf("did not converge in %d iterations", f$iterations))
      design <- paste(c(
        sprintf("%d objects", nrow(f$objects)),
        sprintf("%.0f comparisons", f$n_comparisons),
        if (!is.null(f$judges)) sprintf("%d judges", nrow(f$judges))),
        collapse = " · ")
      dep_rows <- if (!is.null(f$dependence))
        lapply(seq_len(nrow(f$dependence)), function(r)
          stat_row(sprintf("Within-judge %s",
                           gsub("_", "-", f$dependence$effect[r])),
                   sprintf("%.2f logits (%s)", f$dependence$estimate[r],
                           p_lab(f$dependence$p[r]))))
      tagList(
        div(class = "stat-head",
            model_lab, " · conditional ML · ", conv,
            if (isTRUE(f$clustered)) " · SEs clustered by judge"),
        stat_rows(
          stat_row("Pairwise chi-square",
                   sprintf("%.2f on %d df, %s", f$total_chisq, f$total_df,
                           p_lab(f$total_p))),
          stat_row("Design", design),
          stat_row("Object separation index",
                   if (finite1(f$osi$PSI)) sprintf("%.3f", f$osi$PSI)
                   else "—"),
          if (graded && !is.null(f$thr_structure))
            stat_row("Threshold structure",
                     if (identical(f$thr_structure, "pc"))
                       "principal components (spread)"
                     else "free symmetric"),
          dep_rows))
    },
    code = function() "fit_summary_table(bt)")
  # routine handling notes (the old text panel printed bt$notes)
  output$btl_fitsum_notes <- renderUI({
    f <- bfit()
    if (!length(f$notes)) return(NULL)
    sprintf("Note. %s.", paste(f$notes, collapse = "; "))
  })
  register_table("btl_obj_tbl", function() bfit()$objects,
                 function() {
    d <- curate(bfit()$objects, "btl_obj", full = isTRUE(input$btl_full))
    # fit residual, infit and outfit are flagged by num_dt, as on every table
    num_dt(d, selection = "single")
  }, code = function() "# bt from the Data page\nbt$objects")
  # the object selected by clicking a row of the table drives the object
  # characteristic curve on the right (master-detail, as the item table does)
  sel_object <- reactive({
    b <- bfit(); i <- input$btl_obj_tbl_rows_selected
    if (length(i)) b$objects$object[i] else b$objects$object[1]
  })
  register_table("btl_pairs_tbl", function() bfit()$pairs,
                 function() {
    d <- bfit()$pairs
    d$chisq <- NULL   # residual^2; redundant on screen, kept in the CSV
    num_dt(d)
  }, code = function() "bt$pairs")
  register_table("btl_judges_tbl", function() {
    validate(need(!is.null(bfit()$judges), "No judge column was nominated."))
    bfit()$judges
  }, function() {
    validate(need(!is.null(bfit()$judges), "No judge column was nominated."))
    d <- bfit()$judges
    d$misfit <- ifelse(!is.na(d$fit_resid) & abs(d$fit_resid) > 2.5, "*", "")
    num_dt(curate(d, "btl_judge", full = isTRUE(input$btl_judges_full)))
  }, code = function() "bt$judges")
  register_plot("btl_plot", function() plot_btl(bfit()),
                code = function() "# bt from the Data page\nplot_btl(bt)")
  # object characteristic curve: model expected response against opponent
  # location with per-opponent observed means (dichotomous and graded fits)
  observeEvent(btl_fit(), {
    b <- btl_fit()
    if (!is.null(b)) {
      jf <- input$bt_jfactors
      updateSelectizeInput(session, "bdif_factors",
                           choices = if (length(jf)) jf else character(0),
                           selected = jf)
    }
    # judge-group DIF results belong to the fit they came from
    bdif_res(NULL)
  }, ignoreNULL = FALSE)
  register_plot("btl_occ", function() {
    o <- sel_object(); req(o %in% bfit()$objects$object)
    plot_btl_icc(bfit(), o)
  }, w = 8, h = 5.5, code = function()
    paste0(sprintf('plot_btl_icc(bt, "%s")', sel_object() %||% ""),
           "\n# all objects: one page each in the PDF download"))
  output$btl_occ_all_pdf <- downloadHandler(
    filename = function() "occ_all_objects.pdf",
    content = function(file) {
      b <- bfit()
      pdf(file, width = 8, height = 5.5, onefile = TRUE)
      on.exit(dev.off(), add = TRUE)
      # a failed object still gets its page: a placeholder naming the
      # error, so no object silently vanishes from the batch
      for (o in b$objects$object)
        tryCatch(plot_btl_icc(b, o), error = function(e) {
          plot.new()
          text(0.5, 0.5, sprintf("%s: %s", o, conditionMessage(e)))
        })
    })
  # graded (ordinal) fits carry thresholds and category curves; the flag
  # hides both cards entirely for dichotomous fits
  output$btl_graded <- reactive({
    b <- btl_fit()
    !is.null(b) && !is.null(b$m) && b$m >= 2
  })
  outputOptions(output, "btl_graded", suspendWhenHidden = FALSE)
  register_table("btl_thr_tbl", function() bfit()$thresholds, function() {
    validate(need(!is.null(bfit()$thresholds),
                  "Dichotomous fit: no threshold structure to show."))
    num_dt(bfit()$thresholds)
  }, code = function() "bt$thresholds")
  register_table("btl_comp_tbl", function() bfit()$components, function() {
    validate(need(!is.null(bfit()$components),
                  "Dichotomous fit: no threshold components to show."))
    num_dt(bfit()$components)
  }, code = function() "bt$components")
  register_plot("btl_cats", function() plot_btl_categories(bfit()),
                code = function() "plot_btl_categories(bt)")

  # within-judge dependence (Independence > Local): exposure and carry-over
  # effects estimated when a judgment-order column was nominated
  # keyed on the comparison-level data, not the effects table: when every
  # effect is dropped (no informative comparisons, or separation) the panel
  # must still show the audit table that explains why
  output$has_btl_dep <- reactive({
    b <- btl_fit()
    !is.null(b) && !is.null(b$dependence_data)
  })
  outputOptions(output, "has_btl_dep", suspendWhenHidden = FALSE)
  register_table("btl_dep_tbl", function() {
    validate(need(!is.null(bfit()$dependence_data),
                  "Nominate a judgment-order column in the Data roles to estimate within-judge dependence."))
    bfit()$dependence
  }, function() {
    validate(need(!is.null(bfit()$dependence_data),
                  "Nominate a judgment-order column in the Data roles to estimate within-judge dependence."))
    validate(need(!is.null(bfit()$dependence),
                  paste("No dependence effect was estimable:",
                        paste(bfit()$notes, collapse = "; "))))
    num_dt(bfit()$dependence, p_bold = "p")
  }, code = function() "bt$dependence")
  # the graphical display of the selected dependence effect
  register_plot("btl_dep_plot", function() {
    b <- bfit(); req(!is.null(b$dependence_data))
    validate(need(!is.null(b$dependence),
                  paste("No dependence effect was estimable:",
                        paste(b$notes, collapse = "; "))))
    e <- input$btl_dep_effect %||% "exposure"
    validate(need(e %in% b$dependence$effect,
                  "This effect had no informative comparisons (or separated) and was dropped; see the notes."))
    plot_btl_dependence(b, e)
  }, w = 8, h = 5.5, code = function()
    sprintf('plot_btl_dependence(bt, "%s")', input$btl_dep_effect %||% "exposure"))
  # every comparison with its history covariates, for deep interrogation
  register_table("btl_dep_comps", function() bfit()$dependence_data,
                 function() {
    validate(need(!is.null(bfit()$dependence_data),
                  "Nominate a judgment-order column in the Data roles."))
    num_dt(bfit()$dependence_data, page_len = 20)
  }, code = function() "bt$dependence_data")

  # ---------------------- BTL trait dimensionality (loops + swirl) --------
  # the pair-structure analogues of transitivity (one order?) and residual
  # PCA (a second attribute steering contests?); both cached per fit
  btl_trans <- reactive({ b <- bfit(); req(!is.null(b)); btl_transitivity(b) })
  btl_dim   <- reactive({ b <- bfit(); req(!is.null(b)); btl_dimensionality(b) })
  register_plot("btl_scree", function() plot_btl_scree(btl_dim()),
                w = 7, h = 5, code = function()
                  "plot_btl_scree(btl_dimensionality(bt))")
  register_plot("btl_dim_map", function() plot_btl_dim_map(btl_dim()),
                w = 7, h = 5.5, code = function()
                  "plot_btl_dim_map(btl_dimensionality(bt))")
  register_table("btl_bimensions_tbl", function() btl_dim()$bimensions,
                 function() num_dt(btl_dim()$bimensions),
                 code = function() "btl_dimensionality(bt)$bimensions")
  register_table("btl_trans_tbl", function() btl_trans()$summary,
                 function() num_dt(btl_trans()$summary),
                 code = function() "btl_transitivity(bt)$summary")
  # structural lens (Trait tab): which objects sit in the most loops
  register_plot("btl_involve_plot",
                function() plot_btl_transitivity(btl_trans(), by = "object"),
                w = 7, h = 5, code = function()
                  'plot_btl_transitivity(btl_transitivity(bt), by = "object")')
  # per-judge lens (Persons tab): each judge's consistency
  register_plot("btl_judge_consist", function() {
    tr <- btl_trans()
    validate(need(!is.null(tr$judges),
      "No judge reached enough compared triples for a consistency estimate."))
    plot_btl_transitivity(tr, by = "judge")
  }, w = 7, h = 5, code = function()
    'plot_btl_transitivity(btl_transitivity(bt), by = "judge")')
  register_table("btl_trans_judges_tbl",
                 function() btl_trans()$judges,
                 function() {
                   j <- btl_trans()$judges
                   validate(need(!is.null(j),
                     "Per-judge consistency needs judges with enough compared triples."))
                   num_dt(j)
                 }, code = function() "btl_transitivity(bt)$judges")

  # -------------------------------------------- BTL DIF by judge group --
  # the judge grouping handed to btl_dif() / plot_btl_icc(): the nominated
  # factor's value on each judge's first row, named by judge
  # the nominated judge factors as a named list of judge -> level maps (each
  # factor is constant within judge, so the map is well defined); btl_dif
  # models several factors jointly
  bdif_factor_maps <- function() {
    if (is.null(btl_fit()))
      stop("run a paired-comparisons (BTL) analysis first")
    fcs <- input$bdif_factors
    if (is.null(fcs) || !length(fcs))
      stop("choose one or more judge factors in the sidebar")
    df <- raw_data()
    jc <- input$bt_judge
    if (is.null(jc) || identical(jc, NONE) || !jc %in% names(df))
      stop("judge-group DIF needs the judge column nominated on the Data page")
    jd <- as.character(df[[jc]]); first <- !duplicated(jd)
    setNames(lapply(fcs, function(fc) {
      if (!fc %in% names(df)) stop("column not found: ", fc)
      setNames(as.character(df[[fc]])[first], jd[first])
    }), fcs)
  }
  # the factors of one displayed ANOVA term, matched against a run's factor
  # maps: a factor name that itself contains ":" must match the whole term
  # before the term is split into pieces
  bdif_term_vars <- function(term, maps) {
    if (term %in% names(maps)) return(term)
    intersect(strsplit(term, ":", fixed = TRUE)[[1]], names(maps))
  }
  # the judge -> cell map for one ANOVA term of the DISPLAYED run (frozen at
  # run time, so later sidebar edits cannot silently regroup the overlay)
  bdif_term_group <- function(term) {
    r <- bdif_res()
    if (is.null(r) || is.null(r$run_maps))
      stop("run the DIF analysis first")
    maps <- r$run_maps
    vars <- bdif_term_vars(term, maps)
    if (!length(vars))
      stop("the selected term no longer matches the run's factors; re-run")
    js <- names(maps[[1]])
    setNames(do.call(paste, c(lapply(vars, function(v) maps[[v]][js]),
                              sep = ":")), js)
  }
  # backtick nonsyntactic column names in emitted code
  bq <- function(s) ifelse(grepl("^[a-zA-Z.][a-zA-Z0-9._]*$", s),
                           s, sprintf("`%s`", s))
  # the reproducible-code lines building the DISPLAYED run's factor list
  bdif_code_grp <- function() {
    r <- bdif_res()
    fcs <- if (!is.null(r)) r$factors else (input$bdif_factors %||% "factor")
    jc <- (if (!is.null(r)) r$run_judge_col else input$bt_judge) %||% "judge"
    one <- function(fc) sprintf(
      "%s = setNames(as.character(dat$%s), dat$%s)[!duplicated(dat$%s)]",
      bq(fc), bq(fc), bq(jc), bq(jc))
    paste0("factors <- list(\n  ",
           paste(vapply(fcs, one, ""), collapse = ",\n  "), ")")
  }
  bdif_alpha <- reactive(clamp01(input$bdif_alpha, 0.05))
  bdif_effects <- reactive(input$bdif_effects %||% "main")
  # the displayed run's own settings, for styling and snippets (falling back
  # to the live sidebar before any run)
  bdif_shown_alpha <- function() {
    r <- bdif_res(); if (!is.null(r)) r$alpha else bdif_alpha()
  }
  bdif_shown_effects <- function() {
    r <- bdif_res(); if (!is.null(r)) r$effects else bdif_effects()
  }
  output$bdif_multifactor <- reactive(length(input$bdif_factors) > 1L)
  outputOptions(output, "bdif_multifactor", suspendWhenHidden = FALSE)
  # judge factors nominated (or changed) after the run still reach the factor
  # select; the results themselves are computed on request only
  observeEvent(input$bt_jfactors, {
    jf <- input$bt_jfactors
    keep <- intersect(input$bdif_factors, jf)
    updateSelectizeInput(session, "bdif_factors",
                         choices = if (length(jf)) jf else character(0),
                         selected = if (length(keep)) keep else jf)
  }, ignoreNULL = FALSE)
  bdif_res <- reactiveVal(NULL)
  observeEvent(input$bdif_run, {
    maps <- tryCatch(bdif_factor_maps(), error = function(e) e)
    if (inherits(maps, "error")) {
      showNotification(conditionMessage(maps), type = "error", duration = 10)
      bdif_res(NULL); return()
    }
    r <- withProgress(message = "Resolving objects by judge group…",
                      value = 0.4,
                      tryCatch(btl_dif(btl_fit(),
                                       factors = maps,
                                       effects = bdif_effects(),
                                       alpha = bdif_alpha()),
                               error = function(e) e))
    # freeze the run's configuration alongside its results, so the overlay
    # grouping and the reproducible snippets keep describing THIS table even
    # if the sidebar changes before the next run
    if (!inherits(r, "error")) {
      r$run_maps <- maps
      r$run_judge_col <- input$bt_judge
    }
    if (inherits(r, "error")) {
      showNotification(paste("DIF analysis failed:", conditionMessage(r)),
                       type = "error", duration = 10)
      bdif_res(NULL)
    } else bdif_res(r)
  })
  register_table("bdif_anova_tbl", function() {
    r <- bdif_res(); req(!is.null(r)); r$summary
  }, function() {
    r <- bdif_res()
    validate(need(!is.null(r),
                  "Choose one or more judge factors in the sidebar and run the DIF analysis."))
    d <- r$summary[, intersect(c("object", "term", "F_uniform",
                                 "p_uniform_adj", "F_nonuniform",
                                 "p_nonuniform_adj", "superseded"),
                               names(r$summary)), drop = FALSE]
    # a superseded row's flags are read on its higher-order term instead
    d$superseded <- ifelse(d$superseded, "(superseded)", "")
    # no boolean DIF flags: the adjusted probabilities turn red at the RUN's
    # alpha (frozen with the result, not the live sidebar)
    dt <- style_lo_red(num_dt(d, selection = "single"), d,
                       "p_uniform_adj", bdif_shown_alpha())
    style_lo_red(dt, d, "p_nonuniform_adj", bdif_shown_alpha())
  }, code = function()
    paste0("# dat: the comparison data; bt: the fit from the Data page\n",
           bdif_code_grp(), "\n",
           sprintf('btl_dif(bt, factors, effects = "%s", alpha = %s)$summary',
                   bdif_shown_effects(), bdif_shown_alpha())))
  register_table("bdif_sizes_tbl", function() {
    r <- bdif_res(); req(!is.null(r), !is.null(r$sizes)); r$sizes
  }, function() {
    r <- bdif_res()
    validate(need(!is.null(r),
                  "Choose one or more judge factors in the sidebar and run the DIF analysis."))
    validate(need(!is.null(r$sizes),
                  "No object could be resolved (see the notes on the analysis-of-variance panel)."))
    d <- r$sizes[, intersect(c("object", "term", "level_a", "level_b",
                               "difference", "se", "z", "p_adj"),
                             names(r$sizes)), drop = FALSE]
    # no significant/practical flags: the difference turns red at 0.5 logits,
    # the adjusted p at the run's alpha
    dt <- style_mag_red(num_dt(d), d, "difference", 0.5)
    style_lo_red(dt, d, "p_adj", bdif_shown_alpha())
  }, code = function()
    paste0(bdif_code_grp(), "\n",
           sprintf('btl_dif(bt, factors, effects = "%s", alpha = %s)$sizes',
                   bdif_shown_effects(), bdif_shown_alpha())))
  output$bdif_notes <- renderUI({
    r <- bdif_res()
    if (is.null(r) || !length(r$notes)) return(NULL)
    sprintf("Note. %s.", paste(r$notes, collapse = "; "))
  })
  # the (object, term) whose by-group curve is shown: the selected
  # analysis-of-variance row (master-detail, as the Rasch DIF page drives its
  # ICC); the overlay groups judges by that term's factors
  sel_bdif <- reactive({
    r <- bdif_res(); req(!is.null(r))
    i <- input$bdif_anova_tbl_rows_selected
    row <- if (length(i)) r$summary[i, ] else r$summary[1, ]
    list(object = row$object, term = row$term)
  })
  register_plot("bdif_occ", function() {
    b <- bfit(); sb <- sel_bdif(); req(sb$object %in% b$objects$object)
    # surface the builder's own message (missing judge column, unchosen
    # factor, …) instead of one generic hint
    grp <- tryCatch(bdif_term_group(sb$term), error = function(e) e)
    if (inherits(grp, "error"))
      validate(need(FALSE, conditionMessage(grp)))
    plot_btl_icc(b, sb$object, group = grp)
  }, w = 8, h = 5.5, code = function() {
    sb <- sel_bdif(); r <- bdif_res()
    vars <- if (!is.null(r) && !is.null(r$run_maps))
      bdif_term_vars(sb$term, r$run_maps) else
      strsplit(sb$term, ":", fixed = TRUE)[[1]]
    # setNames keeps the judge names, so the emitted grp is a judge -> cell map
    # (as bdif_term_group builds it) rather than an unnamed vector
    paste0(bdif_code_grp(), "\n",
           sprintf('grp <- setNames(do.call(paste, c(factors[c(%s)], sep = ":")), names(factors[[1]]))\n',
                   paste(sprintf('"%s"', vars), collapse = ", ")),
           sprintf('plot_btl_icc(bt, "%s", group = grp)', sb$object %||% ""))
  })

  # ---------------------------------------------------------------- facets --
  facet_dat <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a many-facet (MFRM) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    f$facet_effects[[input$facet_sel]]
  })
  # additive fits carry no interaction card; the footer says so (and how to
  # get one), so toggling the facet structure visibly changes this page
  output$facet_structure_note <- renderUI({
    f <- tryCatch(fit(), error = function(e) NULL)
    if (!inherits(f, "rasch_mfrm") || !is.null(f$interaction)) return(NULL)
    "Additive structure: no item-by-facet terms estimated (choose Interactive in the data roles to test them)."
  })
  register_table("facet_tbl", function() facet_dat(), function()
    # num_dt flags the facet fit residual beyond |2.5|, as on every model table
    num_dt(curate(facet_dat(), "facet", full = isTRUE(input$facets_full))),
    code = function()
    sprintf('fit$facet_effects[["%s"]]', input$facet_sel %||% ""))
  register_plot("facet_plot", function() {
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a many-facet (MFRM) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    plot_facets(f, input$facet_sel)
  }, code = function()
    sprintf('plot_facets(fit, "%s")', input$facet_sel %||% ""))
  facet_int <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm") && !is.null(f$interaction),
                  "Choose the interactive facet structure in the data roles and re-estimate to test item-by-facet interactions."))
    f$interaction_effects
  })
  register_table("facet_int_tbl", function() facet_int(), function() {
    d <- facet_int()
    d$significant <- ifelse(abs(d$gamma) > 1.96 * d$se, "*", "")
    num_dt(d)
  }, code = function() "fit$interaction_effects")

  # ----------------------------------------------------------------- equating --
  eq_ref <- reactive({
    req(input$eq_file)
    a <- tryCatch(read.csv(input$eq_file$datapath, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    validate(need(!is.null(a) && all(c("item", "location") %in% names(a)),
                  "The reference CSV needs columns item, location (and ideally se)."))
    a
  })
  # reference: an uploaded calibration CSV, or a fit kept on the Compare page
  eq_reference <- reactive({
    if (identical(input$eq_source, "kept")) {
      k <- kept_fits()
      validate(need(length(k) >= 1,
                    "Keep a fit on the Compare page to use it as the equating reference."))
      validate(need(!is.null(input$eq_kept) && input$eq_kept %in% names(k),
                    "Choose a kept fit in the sidebar."))
      k[[input$eq_kept]]
    } else {
      validate(need(!is.null(input$eq_file),
                    "Upload a reference calibration (item, location, se) to equate against."))
      eq_ref()
    }
  })
  eq_res <- reactive(equate_tests(fit(), eq_reference(), shift = input$eq_shift))
  register_table("eq_tbl", function() eq_res()$table, function() {
    d <- curate(eq_res()$table, "equate", full = isTRUE(input$eq_full))
    if ("drift" %in% names(d)) d$drift <- ifelse(d$drift, "*", "")
    num_dt(d)
  }, code = function()
    sprintf('eq <- equate_tests(fit, reference, shift = "%s")\neq$table  # reference: data.frame(item, location, se) or another fit',
            input$eq_shift %||% "mean"))
  register_plot("eq_plot", function()
    plot_equate(fit(), eq_reference(), shift = input$eq_shift),
    code = function()
      sprintf('plot_equate(fit, reference, shift = "%s")',
              input$eq_shift %||% "mean"))
  output$dl_anchors <- downloadHandler(
    filename = function() format(Sys.time(), "rasch_anchors_%Y%m%d_%H%M.csv"),
    content = function(file) {
      f <- fit()
      thr <- f$thresholds
      write.csv(data.frame(item = f$items$item[thr$item], k = thr$k,
                           tau = thr$tau), file, row.names = FALSE)
    })

  output$dl_calib <- downloadHandler(
    filename = function() format(Sys.time(), "rasch_calibration_%Y%m%d_%H%M.csv"),
    content = function(file) {
      f <- fit()
      write.csv(data.frame(item = f$items$item, location = f$items$location,
                           se = f$items$se), file, row.names = FALSE)
    })

  # ---------------------------------------------------------------- frames --
  efrm_fit <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_efrm"),
                  "Run a frames (extended frame of reference) analysis to see results here."))
    f
  })
  register_table("frame_tbl", function() efrm_fit()$frames,
                 function() num_dt(curate(efrm_fit()$frames, "frames",
                                          full = isTRUE(input$frames_full))),
                 code = function() "fit$frames")
  register_table("phi_tbl", function() efrm_fit()$phi_table,
                 function() num_dt(efrm_fit()$phi_table),
                 code = function() "fit$phi_table")
  # keep the fit's original set order: merge() sorts by the key, so it is
  # re-matched to fit$set_table$set
  efrm_alpha_tbl <- reactive({
    f <- efrm_fit()
    d <- merge(f$alpha_table, f$set_table[, c("set", "mu", "n_items")],
               by = "set", sort = FALSE)
    d <- d[stats::na.omit(match(f$set_table$set, d$set)), , drop = FALSE]
    rownames(d) <- NULL
    d
  })
  register_table("alpha_tbl", function() efrm_alpha_tbl(),
                 function() num_dt(efrm_alpha_tbl()),
                 code = function() "fit$alpha_table")
  register_plot("frame_plot", function() plot_frames(efrm_fit()),
                code = function() "plot_frames(fit)")
  register_plot("frame_icc", function() {
    f <- efrm_fit()
    req(input$frame_item %in% f$virtual_map$item)
    plot_icc_frames(f, input$frame_item)
  }, code = function()
    sprintf('plot_icc_frames(fit, "%s")', input$frame_item %||% ""))
  output$efrm_cmp <- renderPrint({
    f <- efrm_fit(); cmp <- f$efrm_vs_rasch
    cat(sprintf("Pairwise conditional log-likelihood: frames model %.3f, equal units %.3f\n",
                cmp$ll_efrm, cmp$ll_equal))
    cat(sprintf("2 x improvement: %.3f with %d extra unit parameter(s)\n",
                cmp$two_delta_ll, cmp$extra_parameters))
    cat("(composite likelihood: descriptive; informative for ",
        cmp$informative_for, ")\n", sep = "")
    if (!is.null(cmp$unit_tests)) {
      cat("\nWald tests of the units (H0: unit = 1):\n")
      print(cmp$unit_tests, digits = 3, row.names = FALSE)
    }
    cat(sprintf("\nItem fit residual SD under the frames model: %.3f\n",
                f$item_fit_summary$sd))
  })

  # --------------------------------------------------------- dimensionality --
  dim_subsets <- reactiveVal(NULL)
  observeEvent(input$dim_apply, {
    if (length(input$dim_pos) >= 2 && length(input$dim_neg) >= 2) {
      if (length(intersect(input$dim_pos, input$dim_neg))) {
        showNotification("The two subsets must be disjoint.", type = "error")
      } else dim_subsets(list(pos = input$dim_pos, neg = input$dim_neg))
    } else if (!length(input$dim_pos) && !length(input$dim_neg)) {
      dim_subsets(NULL)
      showNotification(sprintf("Ran the t-test on the automatic split (residual component %d).",
                               pca_k()),
                       type = "message")
    } else {
      showNotification("Nominate at least two items in each subset (or leave both empty).",
                       type = "warning")
    }
    # the magnitude table is computed from the subsets in force at ITS run;
    # a changed split makes it stale
    dm_res(NULL)
  })
  observeEvent(input$pca_component, dm_res(NULL), ignoreInit = TRUE)
  # the residual principal component that, when no manual subsets are named,
  # defines the t-test default split
  pca_k <- reactive({
    k <- suppressWarnings(as.integer(input$pca_component %||% 1))
    if (is.na(k) || k < 1L) 1L else k
  })
  dim_res <- reactive({
    s <- dim_subsets()
    if (is.null(s)) dimensionality_test(fit(), component = pca_k())
    else dimensionality_test(fit(), items_positive = s$pos,
                             items_negative = s$neg)
  })
  output$dim_txt <- renderPrint({
    dt <- dim_res()
    if (!is.null(dt$note)) { cat(dt$note); return(invisible()) }
    cat(sprintf("Item split: %s\n", dt$split))
    cat(sprintf("First residual eigenvalue: %.3f\n", dt$first_eigenvalue))
    cat(sprintf("Significant person t-tests: %.1f%%  (exact 95%% CI %.1f%% to %.1f%%, n = %d)\n",
                100 * dt$prop_significant, 100 * dt$ci[1], 100 * dt$ci[2], dt$n))
    cat(sprintf("Persons excluded (extreme on a subset): %d\n", dt$n_excluded_extreme))
    cat(sprintf("Verdict: %s\n", if (dt$multidimensional)
      "lower CI exceeds 5% - unidimensionality is questionable"
      else "consistent with unidimensionality"))
    if (!is.null(dt$paired_t))
      cat(sprintf("Paired t-test of subset means: mean difference %.3f, t = %.2f (df %.0f), p = %s\n",
                  dt$paired_t$mean_difference, dt$paired_t$t,
                  dt$paired_t$df, fmt_p(dt$paired_t$p)))
    cat("\nSubset A items:\n ", paste(dt$items_positive, collapse = ", "), "\n")
    cat("Subset B items:\n ", paste(dt$items_negative, collapse = ", "), "\n")
  })
  register_code("dim", function() {
    s <- dim_subsets()
    if (is.null(s)) {
      k <- pca_k()
      return(if (k == 1L) "dimensionality_test(fit)"
             else sprintf("dimensionality_test(fit, component = %d)", k))
    }
    sprintf("dimensionality_test(fit, items_positive = %s,\n  items_negative = %s)",
            qvec(s$pos), qvec(s$neg))
  })

  # magnitude of multidimensionality (Andrich 2016): needs every item in a
  # subset; the component split satisfies this by construction, manual subsets may not
  dm_res <- reactiveVal(NULL)
  observeEvent(input$dm_run, {
    f <- fit()
    s <- dim_subsets()
    if (is.null(s)) {
      dr <- dim_res()
      if (!is.null(dr$note)) {
        showNotification(paste("No usable subsets:", dr$note), type = "warning")
        return()
      }
      s <- list(pos = dr$items_positive, neg = dr$items_negative)
    }
    allit <- c(s$pos, s$neg)
    if (!setequal(allit, f$items$item) || anyDuplicated(allit) > 0) {
      left <- setdiff(f$items$item, allit)
      showNotification(paste0(
        "The magnitude estimate needs every item assigned to exactly one subset. ",
        if (length(left)) paste0("Unassigned: ", paste(left, collapse = ", "), ". ") else "",
        "Adjust subsets A and B (or leave both empty for the selected component's split) and re-run the t-test first."),
        type = "warning", duration = 10)
      return()
    }
    r <- withProgress(message = "Subtest re-analysis…", value = 0.4,
                      tryCatch(dimensionality_magnitude(f, list(s$pos, s$neg)),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Magnitude estimate failed:", conditionMessage(r)),
                       type = "error", duration = 10)
    else dm_res(r)
  })
  output$dm_tbl <- renderDT({
    r <- dm_res()
    validate(need(!is.null(r),
                  "Press the button; the current subsets (manual, or the selected component's split) are combined into super-items and the two reliability calculations compared."))
    num_dt(r$table)
  })
  output$dm_tbl_csv <- downloadHandler(
    filename = function() "rasch_dimensionality_magnitude.csv",
    content = function(file) {
      r <- dm_res(); req(!is.null(r))
      write.csv(r$table, file, row.names = FALSE)
    })
  register_code("dm_tbl", function()
    "dimensionality_magnitude(fit, list(subset_a, subset_b))$table")
  register_plot("scree", function() plot_scree(fit()),
                code = function() "plot_scree(fit)")
  register_table("loadings_tbl", function() residual_pca(fit())$loadings_matrix,
                 function() {
    d <- residual_pca(fit())$loadings_matrix   # up to the first 10 components
    # roomier rows than the compact default, so the table reads at a weight
    # closer to the biplot beside it (matching the local dependence page)
    num_dt(d, paging = FALSE) |>
      formatStyle(names(d), fontSize = "1rem",
                  paddingTop = "7px", paddingBottom = "7px")
  }, code = function() "residual_pca(fit)$loadings_matrix")
  register_table("eigen_tbl", function() residual_pca(fit())$eigen_table,
                 function() {
    d <- residual_pca(fit())$eigen_table
    d$proportion <- 100 * d$proportion
    d$cumulative <- 100 * d$cumulative
    names(d)[match(c("proportion", "cumulative"), names(d))] <-
      c("Proportion %", "Cumulative %")
    num_dt(d) |> formatRound(c("Proportion %", "Cumulative %"), 1)
  }, code = function() "residual_pca(fit)$eigen_table")
  # biplot of the first two residual components; its card grows with the item
  # count so it stays level with the loadings table beside it (as on the local
  # dependence page)
  # an item-panel height that grows with the item count so a table and its
  # plot stay level, clamped so it never runs tiny or huge (~520px at 10 items)
  item_panel_px <- function(n) as.integer(max(400L, min(820L, 260L + 26L * n)))
  biplot_px <- function()
    item_panel_px(tryCatch(nrow(fit()$items), error = function(e) 10L))
  register_plot("pca_biplot", function() plot_pca_biplot(fit()),
                w = 7, h = 7, px = biplot_px,
                code = function() "plot_pca_biplot(fit)")

  # -------------------------------------------------------- local dependence --
  # the residual correlations as two paired matrices, each with its heatmap:
  # the Yen Q3 correlation (diagonal 1.00) and the adjusted Q3* (each Q3 less
  # the average off-diagonal, diagonal empty). Both tables show the lower
  # triangle only, so each item pair is read once, matching the heatmaps.
  rc_all <- reactive(residual_correlations(fit()))
  q3_mat <- reactive(rc_all()$matrix)
  q3s_mat <- reactive(rc_all()$star_matrix)
  # each matrix is flagged by its own rule at the shared threshold: the raw
  # matrix by |Q3| (Yen 1993), the adjusted matrix by Q3* (Christensen,
  # Makransky & Horton 2017). Flagging the raw matrix by the Q3* rule would
  # redden pairs whose own Q3 is nowhere near the cut.
  q3_flag <- reactive(abs(q3_mat()) >= ld_flag())
  q3s_flag <- reactive(q3s_mat() >= ld_flag())
  register_table("cormat_q3_tbl",
                 function() data.frame(item = rownames(q3_mat()), q3_mat(),
                                       check.names = FALSE),
                 function() tri_dt(q3_mat(), flagged = q3_flag()),
                 code = function() "residual_correlations(fit)$matrix")
  register_table("cormat_q3s_tbl",
                 function() data.frame(item = rownames(q3s_mat()), q3s_mat(),
                                       check.names = FALSE),
                 function() tri_dt(q3s_mat(), flagged = q3s_flag()),
                 code = function() "residual_correlations(fit)$star_matrix")
  # on-screen height grows with the item count (item_panel_px), so the triangle
  # keeps readable cells and stays level with the table beside it
  cormat_px <- function()
    item_panel_px(tryCatch(ncol(q3_mat()), error = function(e) 10L))
  register_plot("rcor_q3", function() plot_resid_cor(fit(), stat = "q3"),
                w = 8, h = 7, px = cormat_px,
                code = function() 'plot_resid_cor(fit, stat = "q3")')
  register_plot("rcor_q3s", function() plot_resid_cor(fit(), stat = "q3star"),
                w = 8, h = 7, px = cormat_px,
                code = function() 'plot_resid_cor(fit, stat = "q3star")')
  # the flag threshold (shared by both matrices' red highlighting); defaults to
  # the conventional Q3* > 0.2 (Christensen, Makransky & Horton 2017)
  ld_flag <- reactive({
    fl <- input$ld_flag
    if (is.null(fl) || is.na(fl) || fl <= 0) 0.2 else fl
  })

  # response dependence magnitude (Andrich & Kreiner resolved-item refit)
  dep_res <- reactiveVal(NULL)
  observeEvent(input$run_dep, {
    f <- fit()
    req(input$dep_item %in% f$items$item, input$ind_item %in% f$items$item)
    r <- withProgress(message = "Resolving and re-analysing…", value = 0.4,
                      tryCatch(dependence_magnitude(f,
                                                    dependent = input$dep_item,
                                                    independent = input$ind_item),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Dependence estimate failed:", conditionMessage(r)),
                       type = "error", duration = 10)
    else dep_res(r)
  })
  output$dep_txt <- renderPrint({
    r <- dep_res()
    validate(need(!is.null(r),
                  "Choose the dependent and independent items and press the button."))
    cat(sprintf("Dependence of %s on %s: d = %.3f logits (se %.3f), z = %.2f, p = %s\n",
                r$dependent, r$independent, r$d, r$se, r$z, fmt_p(r$p)))
  })
  output$dep_tbl <- renderDT({
    r <- dep_res()
    validate(need(!is.null(r), ""))
    num_dt(r$thresholds)
  })
  register_code("dep_tbl", function()
    sprintf('dependence_magnitude(fit, dependent = "%s", independent = "%s")',
            input$dep_item %||% "", input$ind_item %||% ""))
  output$dep_tbl_csv <- downloadHandler(
    filename = function() "rasch_dependence_thresholds.csv",
    content = function(file) {
      r <- dep_res(); req(!is.null(r))
      write.csv(r$thresholds, file, row.names = FALSE)
    })

  # spread test against Andrich's least upper bounds
  spread_res <- reactiveVal(NULL)
  observeEvent(input$run_spread, {
    r <- withProgress(message = "Principal-components refit…", value = 0.4,
                      tryCatch(spread_test(fit()), error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Spread test failed:", conditionMessage(r)),
                       type = "error", duration = 10)
    else spread_res(r)
  })
  output$spread_tbl <- renderDT({
    r <- spread_res()
    validate(need(!is.null(r),
                  "Press the button to run the spread test (needs at least one polytomous item, e.g. after combining a subtest)."))
    d <- r
    d$dependent <- ifelse(d$dependent, "*", "")
    num_dt(d)
  })
  register_code("spread_tbl", function() "spread_test(fit)")
  output$spread_tbl_csv <- downloadHandler(
    filename = function() "rasch_spread_test.csv",
    content = function(file) {
      r <- spread_res(); req(!is.null(r))
      write.csv(r, file, row.names = FALSE)
    })

  # ---------------------------------------------------------------- guessing --
  guess_res <- reactiveVal(NULL)
  observeEvent(input$run_guess, {
    f <- fit()
    if (max(f$m) != 1L) {
      showNotification("Tailored analysis applies to dichotomous analyses only.",
                       type = "warning")
      return()
    }
    ch <- clamp01(input$guess_chance, 0.25)
    anc <- if (length(input$guess_anchors)) input$guess_anchors else NULL
    r <- withProgress(message = "Tailored analysis (three re-analyses)…",
                      value = 0.3,
                      tryCatch(tailored_analysis(f, chance = ch,
                                                 anchor_items = anc),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Tailored analysis failed:", conditionMessage(r)),
                       type = "error", duration = 10)
    else guess_res(r)
  })
  output$guess_txt <- renderPrint({
    validate(need(max(fit()$m) == 1L,
                  "Run a dichotomous (multiple-choice) analysis to use the tailored guessing procedure."))
    r <- guess_res()
    validate(need(!is.null(r),
                  "Set the chance level in the sidebar and press the button."))
    cat(sprintf("Responses set to missing (P below chance %.2f): %d\n",
                r$chance, r$n_removed))
    cat("Anchor items for the common origin:",
        paste(r$anchor_items, collapse = ", "), "\n")
    up <- sum(r$table$z > 1.96, na.rm = TRUE)
    cat(sprintf("Items significantly harder in the tailored analysis (z > 1.96): %d of %d\n",
                up, nrow(r$table)))
    cat("Verdict:", if (up > 0) "guessing indicated"
        else "no guessing signature", "\n")
  })
  register_table("guess_tbl", function() {
    r <- guess_res(); req(!is.null(r)); r$table
  }, function() {
    validate(need(max(fit()$m) == 1L,
                  "Run a dichotomous (multiple-choice) analysis to use the tailored guessing procedure."))
    r <- guess_res()
    validate(need(!is.null(r), "Run the tailored analysis to see the comparison."))
    num_dt(r$table)
  }, code = function()
    sprintf("ta <- tailored_analysis(fit, chance = %s)\nta$table",
            clamp01(input$guess_chance, 0.25)))
  register_plot("guess_plot", function() {
    validate(need(max(fit()$m) == 1L,
                  "Run a dichotomous (multiple-choice) analysis to use the tailored guessing procedure."))
    r <- guess_res()
    validate(need(!is.null(r), "Run the tailored analysis to see the equating plot."))
    plot_equate(r$tailored, r$origin_equated, shift = "none")
  }, code = function()
    'plot_equate(ta$tailored, ta$origin_equated, shift = "none")')

  # ---------------------------------------------------------------- compare --
  kept_fits <- reactiveVal(list())
  observeEvent(input$keep_fit, {
    f <- fit()
    k <- kept_fits()
    lab <- sprintf("%d_%s", length(k) + 1L, f$model)
    k[[lab]] <- f
    kept_fits(k)
    updateSelectizeInput(session, "cmp_ref", choices = names(k),
                         selected = if (!is.null(input$cmp_ref) &&
                                        input$cmp_ref %in% names(k))
                           input$cmp_ref else names(k)[1])
    updateSelectizeInput(session, "eq_kept", choices = names(k),
                         selected = if (!is.null(input$eq_kept) &&
                                        input$eq_kept %in% names(k))
                           input$eq_kept else names(k)[length(k)])
    showNotification(sprintf("Kept '%s' (%d fit(s) held).", lab, length(k)),
                     type = "message", duration = 5)
  })
  observeEvent(input$clear_fits, {
    kept_fits(list())
    updateSelectizeInput(session, "cmp_ref", choices = character(0),
                         selected = character(0))
    updateSelectizeInput(session, "eq_kept", choices = character(0),
                         selected = character(0))
    showNotification("Cleared kept fits.", type = "message", duration = 4)
  })
  cmp_res <- reactive({
    k <- kept_fits()
    validate(need(length(k) >= 2,
                  "Keep at least two fits (run, keep, change settings, run, keep) to compare."))
    ref <- if (!is.null(input$cmp_ref) && input$cmp_ref %in% names(k))
      input$cmp_ref else 1
    as.data.frame(do.call(compare_fits, c(k, list(reference = ref))))
  })
  register_table("cmp_tbl", function() cmp_res(), function() {
    d <- cmp_res()
    if ("same_data" %in% names(d))
      d$same_data <- ifelse(d$same_data, "yes", "no")
    num_dt(curate(d, "compare", full = isTRUE(input$cmp_full)))
  }, code = function()
    "compare_fits(fit_a = f1, fit_b = f2, reference = 1)  # keep fits, then compare")

  # ----------------------------------------------------------------- export --
  # single-file HTML report; one content function feeds both the Export-tab
  # button and the navbar icon link (Rasch fits only; BTL is notified)
  report_content <- function(file) {
    if (!is.null(btl_fit())) {
      showNotification("The HTML report covers Rasch analyses; paired-comparison (BTL) fits are not yet supported.",
                       type = "warning", duration = 8)
      stop("report unavailable for a BTL fit")
    }
    f <- fit_or_null()
    if (is.null(f)) {
      showNotification("Run an analysis first, then download the report.",
                       type = "warning", duration = 8)
      stop("no fit to report")
    }
    withProgress(message = "Building the HTML report…", value = 0.4,
                 report_html(f, file))
  }
  output$dl_report <- downloadHandler(
    filename = function() "rasch_report.html", content = report_content)
  output$dl_report_nav <- downloadHandler(
    filename = function() "rasch_report.html", content = report_content)

  output$dl_zip <- downloadHandler(
    filename = function() format(Sys.time(), "rasch_results_%Y%m%d_%H%M.zip"),
    content = function(file) {
      f <- fit()
      tmp <- file.path(tempdir(), paste0("rasch_", as.integer(Sys.time())))
      withProgress(message = "Writing all tables and plots…", value = 0.4, {
        save_outputs(f, tmp,
                     formats = if (length(input$exp_formats)) input$exp_formats else "png",
                     item_plots = isTRUE(input$exp_items))
      })
      owd <- setwd(tmp); on.exit(setwd(owd), add = TRUE)
      utils::zip(zipfile = file, files = list.files(".", recursive = TRUE),
                 flags = "-r9Xq")
    })
}

shinyApp(ui, server)
