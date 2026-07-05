# rmt Shiny GUI
# ---------------------------------------------------------------------------
# A modern bslib interface to the full rmt analysis: data upload with ID,
# person-factor, and item column nomination; pairwise conditional ML
# estimation (Andrich & Luo 2003); the complete test-of-fit suite;
# every diagnostic plot with per-plot PNG and PDF downloads; and one-click
# export of all tables and plots as a ZIP archive.
# Launch with rmt::run_app(), or shiny::runApp() from this folder.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(bsicons)
})

if (requireNamespace("rmt", quietly = TRUE)) {
  library(rmt)
} else {
  rdir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  if (dir.exists(rdir)) {
    for (f in list.files(rdir, "\\.R$", full.names = TRUE)) source(f)
  } else stop("Install rmt, or run the app from inst/shiny in the source tree")
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
.demo_btl <- function(seed = 47, reps = 22) {
  set.seed(seed)
  beta <- setNames(seq(-1.4, 1.4, length.out = 8), sprintf("E%02d", 1:8))
  pr <- t(utils::combn(names(beta), 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps),
                  object_b = rep(pr[, 2], each = reps))
  d$judge <- sprintf("J%02d", sample(1:10, nrow(d), replace = TRUE))
  p <- plogis(beta[d$object_a] - beta[d$object_b])
  p[d$judge == "J09"] <- 0.5
  d$winner <- ifelse(runif(nrow(d)) < p, d$object_a, d$object_b)
  # graded preference for object_a: adjacent-categories probabilities on the
  # same locations with symmetric thresholds (J09 stays random)
  lev <- c("much worse", "a little worse", "a little better", "much better")
  P <- vapply(seq_len(nrow(d)), function(r)
    item_moments(beta[d$object_a[r]] - beta[d$object_b[r]],
                 c(-1.1, 0, 1.1))$P, numeric(4))
  P[, d$judge == "J09"] <- 0.25
  d$preference <- factor(lev[apply(P, 2, function(pp) sample(4, 1, prob = pp))],
                         levels = lev, ordered = TRUE)
  # margin of win as an ordered factor, derived from the graded preference
  # (the extreme categories are "much" wins), for the winner + margin path
  d$margin <- factor(ifelse(d$preference %in% c("much worse", "much better"),
                            "much", "a little"),
                     levels = c("a little", "much"), ordered = TRUE)
  d <- d[sample(nrow(d)), ]
  # judgment order: each judge's sequence 1..n_j in presentation order, so
  # the within-judge dependence analysis can be demonstrated
  d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  # a judge factor (constant within judge) for DIF by judge group
  d$panel <- ifelse(d$judge %in% sprintf("J%02d", 1:5),
                    "panel A", "panel B")
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

# long-format rated demo: 5 items, 6 raters (one erratic), incomplete design
.demo_long <- function(seed = 21, Np = 250) {
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
  rownames(d) <- NULL
  d
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
# defined here (not mid-server) because observers created early reference it
FACTORIAL <- "(all factors: factorial)"

# p-values as text: "%.3f" alone prints a misleading 0.000 for tiny p
fmt_p <- function(p)
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))

# value-box guard: NULL, NA, NaN, and Inf must never reach a conditional or
# a sprintf; such values display as an em dash on a neutral theme
finite1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

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
  .card-header { font-weight: 600; }
  .value-box-title { font-size: .72rem; text-transform: uppercase; letter-spacing: .04em; white-space: nowrap; }
  .value-box-value { font-size: 1.45rem; }
  pre, .shiny-text-output { white-space: pre-wrap; font-size: .82rem; }
  .btn-xs { padding: .1rem .5rem; font-size: .75rem; }
  .form-label { font-weight: 600; font-size: .85rem; }
  /* APA-style tables: tabular numerals, no vertical rules, no zebra,
     a strong rule under the header row */
  table.dataTable { font-size: .85rem; }
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
  /* card headers: title left, action chips right, with a small gap
     between chips (the chips row right-aligns even when there is no
     title, via margin-left:auto) */
  .rmt-card-header { display: flex; align-items: center; justify-content: space-between; gap: .5rem; }
  .rmt-chips { display: flex; align-items: center; flex-wrap: wrap; gap: .35rem; margin-left: auto; }
  /* inline form controls that sit on a flex row with buttons: strip the
     bottom margin Shiny's containers carry */
  .rmt-inline-check .form-group, .rmt-inline-check .shiny-input-container,
  .rmt-inline-check .checkbox,
  .rmt-inline-select .form-group, .rmt-inline-select .shiny-input-container {
    margin-bottom: 0; width: auto;
  }
  .rmt-inline-select select.form-select { padding: .15rem 1.6rem .15rem .5rem; font-size: .78rem; }
  /* collapsed advanced-settings disclosure inside the sidebar accordion */
  .rmt-advanced { margin-top: .5rem; }
  .rmt-advanced summary { cursor: pointer; font-size: .8rem; font-weight: 600; color: var(--bs-secondary-color); margin-bottom: .35rem; }
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
# exact rmt call reproducing the output, updating with the current selections;
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
  card_header(class = "rmt-card-header",
    if (!is.null(title)) span(title, if (!is.null(info)) info_icon(info)),
    div(class = "rmt-chips",
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
                      "Ratings, long format (MFRM)" = "mfrm",
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
              h6("One row per response"),
              selectInput("lp_person", "Person column", NONE_CH),
              selectInput("lp_item", "Item column", NONE_CH),
              selectInput("lp_score", "Score column", NONE_CH),
              selectizeInput("lp_facets", "Facet columns (e.g. rater)", NULL,
                             multiple = TRUE,
                             options = list(placeholder = "choose at least one")),
              selectInput("lp_interaction", "Item-by-facet interaction (optional)", NONE_CH),
              p(class = "text-muted small",
                "Each item x facet combination is calibrated jointly; facet severities are reported with SEs and fit. An interaction lets one facet be more or less severe on particular items.")
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
                "The Bradley-Terry-Luce model: the conditional (person-free) form of the dichotomous Rasch model, estimated by the same conventions. A judge column enables the judge fit table and clusters the standard errors by judge. Results appear on the BTL tab.")
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
            tags$details(class = "rmt-advanced",
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
    uiOutput("vboxes"),
    # DT cards sit inside plain divs: the grid row would otherwise stretch
    # them to equal height and crop the taller table mid-row
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      div(tableCard("fitsum_tbl", "Test of fit",
        info = "The total item-trait chi-square tests the invariance of item ordering across the trait; a significant result means at least one item's difficulty is not invariant across class intervals (Andrich & Marais 2019). The cell df factor scales each item's chi-square degrees of freedom by the proportion of class-interval cells with enough responders to contribute.",
        footer = uiOutput("fitsum_notes"))),
      div(tableCard("targeting_tbl", "Targeting & reliability"))
    ),
    # server-rendered: the likelihood-ratio card only when it applies
    uiOutput("summary_bottom")
  )

# ---------------------------------------------------------------- ITEMS --
panel_items <- nav_panel("Items", value = "p_items", icon = bs_icon("list-check"),
    uiOutput("items_vboxes"),
    div(class = "mb-2 d-flex align-items-center gap-3 flex-wrap",
        div(class = "rmt-inline-check",
            tooltip(checkboxInput("show_obs", "Observed points", TRUE,
                                  width = "auto"),
                    "Show the observed class-interval points on the category and threshold curves.")),
        downloadButton("dl_anchors", "Save anchors (CSV: item,k,tau)",
                       class = "btn-outline-secondary btn-sm")),
    layout_columns(col_widths = c(7, 5),
      tableCard("items_tbl", "Item statistics",
        controls = cols_switch("items_full"),
                "Click a row to explore that item on the right. Fit residual ~ N(0,1) under fit.",
                info = "Cells are highlighted where a statistic indicates misfit: |fit residual| > 2.5, adjusted chi-square p < 0.05, mean squares outside 0.7-1.3 (Wright & Linacre 1994). No single flag column - read each statistic on its own terms.",
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
                  sliderInput("ex_ng", "Class intervals", min = 2, max = 10,
                              value = 8, step = 1, width = "100%")),
              div(style = "width: 260px;",
                  sliderInput("ex_rng", "Scale range (logits)", min = -8, max = 8,
                              value = c(-5, 5), step = 0.5, width = "100%")))),
          # batch downloads follow the active tab's plot type; the Chi-square
          # tab has no plot, so the buttons hide there
          conditionalPanel("input.items_nav != 'Chi-square'",
            class = "ms-auto",
            div(class = "rmt-chips",
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
               rcode_details("rescore_tbl"))))))
  )

# -------------------------------------------------------------- PERSONS --
panel_persons <- nav_panel("Persons", value = "p_persons", icon = bs_icon("people"),
    uiOutput("persons_vboxes"),
    tableCard("person_tbl", "Person estimates",
        controls = cols_switch("persons_full"),
              "Warm WLE location and SE per person, with raw score, fit statistics, and your ID and factor columns. Click a row to draw that person's kidmap below."),
    accordion(id = "persons_acc", open = "persons_kidmap", class = "mt-3",
      accordion_panel("Kidmap", value = "persons_kidmap",
        plotCard("kidmap",
          info = "The person diagnostic map (Wright, Mead & Ludlow 1980): thresholds the person achieved print to the right of the logit axis, thresholds not achieved to the left; the dashed line inside its confidence band is the person location. Achieved thresholds above the band and unachieved thresholds below it are unexpected responses.",
          controls = div(class = "d-flex align-items-center gap-1 me-1",
            span(class = "small text-secondary", "Confidence"),
            div(class = "rmt-inline-select",
                selectInput("kid_level", NULL,
                            c("90%" = "0.9", "95%" = "0.95", "99%" = "0.99"),
                            selected = "0.95", width = "85px"))),
          extra = tagList(
            downloadButton("kidmap_all_pdf", "PDF (all persons)",
                           class = "btn-outline-secondary btn-xs"),
            downloadButton("kidmap_all_zip", "ZIP (all persons)",
                           class = "btn-outline-secondary btn-xs")))),
      accordion_panel("Person fit", value = "persons_pfit",
        plotCard("pfit")),
      accordion_panel("Fit residual distribution", value = "persons_rdist",
        plotCard("rdist_p")))
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

# ----------------------------------------------------------------- TEST --
panel_test <- nav_panel("Test", value = "p_test", icon = bs_icon("graph-up"),
    div(class = "mb-2", style = "max-width: 300px;",
        sliderInput("ts_rng", "Scale range (logits)", min = -8, max = 8,
                    value = c(-6, 6), step = 0.5, width = "100%")),
    accordion(id = "test_acc", open = "test_tcc",
      accordion_panel("Test characteristic curve", value = "test_tcc",
        plotCard("tcc")),
      accordion_panel("Test information", value = "test_tif",
        plotCard("tif",
          info = "Information across the scale, with the standard error of measurement (SEM = 1/sqrt(information)) overlaid.")),
      accordion_panel("Guttman scalogram", value = "test_guttman",
        plotCard("guttman", height = "640px")))
  )

# ------------------------------------------------------------------ DIF --
panel_dif <- nav_panel("DIF", value = "p_dif", icon = bs_icon("sliders"),
    # Rasch fits: person-factor DIF (hidden while a BTL fit is active)
    conditionalPanel("output.is_btl != true",
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectizeInput("dif_factor", "Person factor", NULL,
                       options = list(placeholder = "run an analysis first")),
        conditionalPanel("input.dif_factor == '(all factors: factorial)'",
          radioButtons("dif_effects", "Model",
                       c("Full factorial" = "factorial",
                         "Main effects only" = "main"))),
        selectizeInput("dif_item", "Item for ICC by group", NULL,
                       options = list(placeholder = "run an analysis first")),
        numericInput("dif_alpha", "Significance level (alpha)", value = 0.05,
                     min = 0.001, max = 0.5, step = 0.01),
        selectInput("dif_padj", "Multiplicity adjustment",
                    c("Benjamini-Hochberg" = "BH",
                      "Holm" = "holm",
                      "Bonferroni" = "bonferroni",
                      "None" = "none")),
        p(class = "text-muted small",
          "ANOVA of standardised residuals: factor effects = uniform DIF; factor x class-interval terms = non-uniform DIF. Probabilities are adjusted across items by the chosen method. With several factors, choose the factorial option to model them jointly: significant interactions supersede their main effects, and Tukey HSD compares the levels of each significant group term."),
        hr(),
        conditionalPanel("input.dif_factor != '(all factors: factorial)'",
          input_task_button("make_split", "Resolve: split this item by this factor",
                            type = "primary", class = "w-100"),
          p(class = "text-muted small mt-2",
            "Replaces the selected item with one item per group level (each level keeps only its own responses) and re-analyses; the split locations quantify the DIF. Splitting works one factor at a time; choose a single factor above."))),
      accordion(id = "dif_acc", open = "dif_anova",
        accordion_panel("DIF analysis of variance", value = "dif_anova",
          tableCard("dif_tbl",
            controls = cols_switch("dif_full"),
                    info = "ANOVA of standardised residuals: a significant factor effect indicates uniform DIF, a significant factor-by-class-interval interaction indicates non-uniform DIF (Andrich & Marais 2019).",
                    footer = uiOutput("dif_note"))),
        # collapsed by default: the DT output suspends while hidden, so the
        # per-item terms table is only computed when the panel is first opened
        accordion_panel("Full ANOVA table", value = "dif_full_panel",
          tableCard("dif_full_tbl",
                    note = "The complete per-item ANOVA: every model term with its df, sums of squares, mean squares, F, and adjusted probability.")),
        accordion_panel("DIF magnitude in logits", value = "dif_size_panel",
          card(card_body(fillable = FALSE,
                 p(class = "text-muted",
                   "Resolves the selected item by the selected factor (in factorial mode: by every significant, non-superseded group term) and reports pairwise location differences in logits with Holm familywise adjustment. Differences of at least the criterion are flagged as practically significant."),
                 layout_columns(col_widths = c(3, 3, 3, 3),
                   numericInput("dif_size_flag", "Practical criterion (logits)",
                                0.5, min = 0.1, step = 0.1),
                   numericInput("dif_size_minn", "Min responders per level", 20,
                                min = 5, step = 5),
                   div(class = "mt-4",
                       input_task_button("dif_size_go", "Compute DIF size",
                                         type = "primary")),
                   conditionalPanel("output.has_difsize == true", class = "mt-4",
                     div(class = "d-flex align-items-center gap-3",
                         cols_switch("difsize_full"),
                         downloadButton("dl_dif_size", "CSV",
                                        class = "btn-outline-secondary btn-xs")))),
                 conditionalPanel("output.has_difsize != true",
                   p(class = "text-muted small mb-0",
                     "Run to see the resolved magnitudes.")),
                 conditionalPanel("output.has_difsize == true",
                   DT::DTOutput("dif_size_tbl"),
                   rcode_details("dif_size_tbl"))))),
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
                   rcode_details("contr_tbl"))))),
        accordion_panel("Pairwise comparisons (Tukey)", value = "dif_tukey",
          conditionalPanel("output.dif_is_factorial == true",
            tableCard("dif_tukey_tbl",
                      note = "Pairwise level comparisons for significant, non-superseded group terms (factorial mode).")),
          conditionalPanel("output.dif_is_factorial != true",
            p(class = "text-muted",
              "Choose the factorial option in the sidebar to see Tukey HSD comparisons."))),
        accordion_panel("Characteristic curves by group", value = "dif_icc_panel",
          plotCard("dif_icc")))
    )),
    # paired-comparison (BTL) fits: differential object functioning by
    # judge group (Bradley-Terry counterpart of the person-factor analysis)
    conditionalPanel("output.is_btl == true",
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectizeInput("bdif_factor", "Judge factor", NULL,
                       options = list(placeholder = "nominate judge factors on the Data page")),
        numericInput("bdif_alpha", "Significance level (alpha)", value = 0.05,
                     min = 0.001, max = 0.5, step = 0.01),
        input_task_button("bdif_run", "Run DIF analysis",
                          type = "primary", class = "w-100"),
        p(class = "text-muted small mt-2",
          "ANOVA of standardised residuals by judge group: a group effect indicates uniform DIF, a group-by-opponent-band interaction non-uniform DIF. Each object is then resolved into one copy per judge group inside a joint refit and the location differences reported in logits.")),
      accordion(id = "bdif_acc", open = "bdif_anova",
        accordion_panel("DIF analysis of variance", value = "bdif_anova",
          tableCard("bdif_anova_tbl",
            info = "ANOVA of the standardised residuals of each object's comparisons, oriented to the object: a significant judge-group effect indicates uniform DIF, a significant group-by-opponent-band interaction non-uniform DIF; probabilities are adjusted across objects by Benjamini-Hochberg.")),
        accordion_panel("DIF magnitude in logits", value = "bdif_size_panel",
          tableCard("bdif_sizes_tbl",
            note = "Pairwise differences between the resolved per-group locations, in logits, with Holm familywise adjustment; differences of at least 0.5 logits are flagged as practically significant.",
            footer = uiOutput("bdif_notes"))),
        accordion_panel("Characteristic curves by group", value = "bdif_occ_panel",
          plotCard("bdif_occ",
            info = "The object characteristic curve with the observed mean response per opponent overlaid separately for each judge group: the graphical display of DIF by judge group.",
            controls = div(class = "d-flex align-items-center gap-1 me-1",
              span(class = "small text-secondary", "Object"),
              div(class = "rmt-inline-select",
                  selectizeInput("bdif_obj", NULL, NULL, width = "120px"))))))
    ))
  )

# --------------------------------------------------------------- FACETS --
panel_facets <- nav_panel("Facets", value = "p_facets", icon = bs_icon("person-badge"),
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectizeInput("facet_sel", "Facet", NULL,
                       options = list(placeholder = "run a many-facet analysis")),
        p(class = "text-muted small",
          "Severities from the joint calibration (positive = more severe). Pooled fit residuals beyond +/-2.5 flag inconsistent levels. Long-format analyses only.")),
      tableCard("facet_tbl", "Facet severities and fit",
        controls = cols_switch("facets_full")),
      plotCard("facet_plot", "Severity caterpillar plot"),
      conditionalPanel("output.has_interaction == true",
        tableCard("facet_int_tbl", "Item-by-facet interactions",
                  "Shown when the analysis was run in interactive facet mode; gamma is the extra severity of a level on a particular item."))
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
            div(class = "rmt-chips",
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
            div(class = "rmt-chips",
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
    p(class = "text-muted small",
      "Trait dependence (dimensionality) threatens local independence: more than one trait driving the responses (Marais & Andrich 2008)."),
    accordion(id = "dim_acc", open = "dim_ttest",
      accordion_panel(
        title = span("Unidimensionality t-test",
                     info_icon("Smith's test: each person is measured separately on the two item subsets and the estimates compared by t-test; unidimensionality is questioned when clearly more than 5% of tests are significant.")),
        value = "dim_ttest",
        layout_columns(col_widths = breakpoints(sm = 12, xl = c(4, 8)),
          div(
            h6("t-test item subsets"),
            selectizeInput("dim_pos", "Subset A", NULL, multiple = TRUE,
                           options = list(placeholder = "positive PC1 loadings")),
            selectizeInput("dim_neg", "Subset B", NULL, multiple = TRUE,
                           options = list(placeholder = "negative PC1 loadings")),
            input_task_button("dim_apply", "Run t-test with these subsets",
                              type = "primary", class = "w-100"),
            p(class = "text-muted small mt-2",
              "Leave both empty (and press the button) to return to the first-contrast split. Persons extreme on either subset are excluded; the proportion of significant tests carries an exact binomial confidence interval.")),
          card(card_body(verbatimTextOutput("dim_txt"), rcode_details("dim"))))),
      accordion_panel("Scree", value = "dim_scree",
        plotCard("scree")),
      accordion_panel("First contrast", value = "dim_pca",
        plotCard("pca_plot")),
      accordion_panel("Loadings", value = "dim_loadings",
        tableCard("loadings_tbl", note = "First 10 components shown.",
                  controls = cols_switch("load_full"))),
      accordion_panel("Eigenvalues", value = "dim_eigen",
        tableCard("eigen_tbl", note = "First 10 eigenvalues shown.")),
      accordion_panel("Magnitude of multidimensionality", value = "dim_magnitude",
        card(
          full_screen = TRUE,
          card_body(
            p(class = "text-muted small",
              "Compares reliability with all items treated as independent (run1) against the subtest analysis in which each subset becomes one polytomous super-item (Andrich 2016). c is the unique-variance loading, rho the latent correlation between the subsets, and A the proportion of common variance. Uses the manual subsets above if set, otherwise the current PC1 split; every item must belong to a subset."),
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
            padding = 12, fillable = FALSE))))
  )

# -------------------------------------------------- INDEPENDENCE: LOCAL --
panel_ld <- nav_panel("Local", value = "p_ld", icon = bs_icon("link-45deg"),
    # paired-comparison (BTL) fits: within-judge dependence estimated from
    # the judgment order; the Rasch Q3 suite hides while a BTL fit is active
    conditionalPanel("output.is_btl == true",
      p(class = "text-muted small",
        "Local (response) dependence threatens local independence (Marais & Andrich 2008): a judge's own history pulling their later judgments, over and above the object locations."),
      card(
        full_screen = TRUE,
        card_header_bar("Within-judge dependence",
          info = "Exposure is the seen-before advantage: the benefit, in logits, an object gains once the judge has already met it. Carry-over is response dependence in the Marais & Andrich sense: the judge's own earlier verdicts on an object pull the later one. Both effects are estimated jointly with the object locations.",
          buttons = conditionalPanel("output.has_btl_dep == true",
            div(class = "rmt-chips",
                downloadButton("btl_dep_tbl_csv", "CSV",
                               class = "btn-outline-secondary btn-xs")))),
        card_body(
          conditionalPanel("output.has_btl_dep != true",
            p(class = "text-muted small mb-0",
              "Nominate a judgment-order column in the Data roles to estimate within-judge dependence.")),
          conditionalPanel("output.has_btl_dep == true",
            DTOutput("btl_dep_tbl"), rcode_details("btl_dep_tbl")),
          padding = 12, fillable = FALSE))),
    conditionalPanel("output.is_btl != true",
    p(class = "text-muted small",
      "Local (response) dependence threatens local independence (Marais & Andrich 2008): responses depending on one another directly, over and above the trait."),
    accordion(id = "ld_acc", open = "ld_q3",
      accordion_panel("Q3 statistics", value = "ld_q3",
        numericInput("ld_flag",
                     "Flag threshold (Q3* above this value flags a pair)",
                     value = 0.2, min = 0.05, max = 0.9, step = 0.05,
                     width = "420px"),
        tableCard("rpairs_tbl",
                  info = "Yen's Q3: the residual correlation of an item pair; Q3* is its excess over the average off-diagonal Q3, the conventional criterion for flagging response dependence (Yen 1984).",
                  footer = uiOutput("rpairs_note"))),
      accordion_panel("Residual correlations (heatmap)", value = "ld_rcor",
        plotCard("rcor", height = "640px")),
      accordion_panel("Subtest (combine dependent items)", value = "ld_subtest",
        card(
          card_body(
            p(class = "text-muted small",
              "Select two or more items to merge into one polytomous super-item and re-analyse; the dependence is absorbed into the subtest."),
            selectizeInput("subtest_items", NULL, NULL, multiple = TRUE,
                           options = list(placeholder = "items to combine")),
            div(input_task_button("make_subtest", "Combine and re-analyse",
                                  type = "primary")),
            uiOutput("subtest_status")))),
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

# ------------------------------------------------------------------ BTL --
panel_btl <- nav_panel("BTL", value = "p_btl", icon = bs_icon("trophy"),
    layout_columns(col_widths = 12,
      card(card_header("Paired comparisons (Bradley-Terry-Luce)"),
           card_body(uiOutput("btl_boxes"),
                     verbatimTextOutput("btl_summary"))),
      tableCard("btl_obj_tbl", "Object locations and fit",
        controls = cols_switch("btl_full"),
                "Conditional (person-free) estimation with sum-zero identification and sandwich standard errors; the fit residual is the log-of-mean-square statistic over each object's comparisons (Andrich & Marais 2019)."),
      plotCard("btl_plot", "Object caterpillar"),
      plotCard("btl_occ", "Object characteristic curve",
        info = "The paired-comparison counterpart of the item characteristic curve: the model expected response for the object against opponent location (the win probability, or the expected graded response), with the observed mean response per opponent overlaid at that opponent's location. Points straying from the curve flag inconsistent quality, exactly as a misfitting item does.",
        controls = div(class = "d-flex align-items-center gap-1 me-1",
          span(class = "small text-secondary", "Object"),
          div(class = "rmt-inline-select",
              selectizeInput("btl_occ_obj", NULL, NULL, width = "120px"))),
        extra = downloadButton("btl_occ_all_pdf", "PDF (all objects)",
                               class = "btn-outline-secondary btn-xs")),
      # graded (ordinal) fits only: hidden entirely for dichotomous fits
      conditionalPanel("output.btl_graded == true",
        tableCard("btl_thr_tbl", "Symmetric thresholds",
                  "Adjacent-categories thresholds of the graded structure, constrained symmetric (tau_k = -tau_(m+1-k)) so the model is invariant to presentation order."),
        tableCard("btl_comp_tbl", "Threshold components",
                  "Spread is the linear component; the skewness component is structurally zero under presentation-order symmetry. Under the PC structure kurtosis is constrained to zero."),
        plotCard("btl_cats", "Category probability curves",
          info = "The probability of each graded response category as a function of the location difference between the two objects; the paired-comparison counterpart of a polytomous item's category curves.")),
      tableCard("btl_pairs_tbl", "Pairwise goodness of fit",
                "Observed against expected win proportions (mean graded responses for a graded fit) for every pair; the total chi-square tests the BTL structure."),
      tableCard("btl_judges_tbl", "Judge fit",
                "Available when a judge column is nominated; an erratic judge carries a large positive fit residual, exactly as an erratic person does."))
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
            div(class = "rmt-chips",
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
        p("The HTML report and ZIP archive cover Rasch analyses. For a paired-comparison (BTL) analysis, download each table as CSV from the BTL page.")))),
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
  title = span("rmt"),
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
        Shiny.addCustomMessageHandler('rmt-nav-vis', function(msg) {
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
  panel_test,
  nav_menu("Independence", value = "menu_structure",
    panel_ld,
    panel_dim),
  nav_menu("Invariance", value = "menu_invariance",
    panel_dif,
    panel_equating,
    panel_guess,
    panel_facets,
    panel_frames),
  nav_menu("More", value = "menu_more",
    panel_btl,
    panel_compare,
    panel_export),
  nav_spacer(),
  nav_item(uiOutput("nav_status")),
  nav_item(downloadLink("dl_report_nav", label = bs_icon("file-earmark-text"),
                        class = "nav-link px-2",
                        title = "Analysis report (HTML)")),
  nav_item(input_dark_mode(id = "app_mode"))
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
                    mfrm = .demo_long(), efrm = .demo_efrm(),
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
    # long-format guesses
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
                         selected = "")
    updateSelectizeInput(session, "bt_jfactors", choices = nm,
                         selected = character(0))
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

  observeEvent(input$lp_facets, {
    sel <- if (!is.null(input$lp_interaction) &&
               input$lp_interaction %in% input$lp_facets)
      input$lp_interaction else NONE
    updateSelectInput(session, "lp_interaction",
                      choices = c(NONE_CH, input$lp_facets), selected = sel)
  }, ignoreNULL = FALSE)

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
                    mfrm = "Ratings, long format (MFRM)",
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
            h2("Welcome to rmt"),
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
              "The exact rmt call reproducing the current run; updates on every estimation."),
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
  # the exact rmt call reproducing the current run (built alongside the fit)
  rcode_str <- reactiveVal(NULL)
  # clear any subtest/split override as soon as a fresh run is requested;
  # fit() short-circuits on the override, so analysis() cannot clear it itself
  observeEvent(input$run, { override_fit(NULL); override_desc(NULL) },
               priority = 10)
  output$has_override <- reactive(!is.null(override_fit()))
  outputOptions(output, "has_override", suspendWhenHidden = FALSE)
  output$override_status <- renderUI({
    if (is.null(override_desc())) return(NULL)
    p(class = "text-warning small mb-1 mt-2",
      paste("Active override -", override_desc()))
  })
  observeEvent(input$reset_override, {
    override_fit(NULL); override_desc(NULL)
    showNotification("Override cleared; showing the base analysis.",
                     type = "message", duration = 5)
  })

  analysis <- eventReactive(input$run, {
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
    code_est <- c(paste0("maxit = ", max(5, input$maxit %||% 60)),
                  paste0("tol = ", format(max(1e-12, input$tol %||% 1e-8))))
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
          if (bt_graded)
            btl(df, object_a = input$bt_a, object_b = input$bt_b,
                response = input$bt_response,
                judge = if (!is.null(input$bt_judge) && input$bt_judge != NONE)
                  input$bt_judge else NULL,
                order = bt_ord,
                count = if (!is.null(input$bt_count) && input$bt_count != NONE)
                  input$bt_count else NULL,
                thresholds = bt_thr,
                maxit = max(5, input$maxit %||% 60),
                tol = max(1e-12, input$tol %||% 1e-8))
          else if (bt_marg)
            btl(df, object_a = input$bt_a, object_b = input$bt_b,
                winner = input$bt_win, margin = input$bt_margin,
                judge = if (!is.null(input$bt_judge) && input$bt_judge != NONE)
                  input$bt_judge else NULL,
                order = bt_ord,
                count = if (!is.null(input$bt_count) && input$bt_count != NONE)
                  input$bt_count else NULL,
                thresholds = bt_thr,
                maxit = max(5, input$maxit %||% 60),
                tol = max(1e-12, input$tol %||% 1e-8))
          else
            btl(df, object_a = input$bt_a, object_b = input$bt_b,
                winner = input$bt_win,
                judge = if (!is.null(input$bt_judge) && input$bt_judge != NONE)
                  input$bt_judge else NULL,
                order = bt_ord,
                count = if (!is.null(input$bt_count) && input$bt_count != NONE)
                  input$bt_count else NULL,
                ties = input$bt_ties %||% "drop",
                maxit = max(5, input$maxit %||% 60),
                tol = max(1e-12, input$tol %||% 1e-8))
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
                     maxit = max(5, input$maxit %||% 60),
                     tol = max(1e-12, input$tol %||% 1e-8),
                     se_method = input$ef_se %||% "hybrid",
                     boot_reps = if (!is.null(input$ef_reps) &&
                                     !is.na(input$ef_reps))
                       max(50, input$ef_reps) else NULL)
        } else if (identical(input$model_type, "mfrm")) {
          if (any(c(input$lp_person, input$lp_item, input$lp_score) == NONE) ||
              !length(input$lp_facets))
            stop("nominate the person, item, score, and at least one facet column")
          code_call <- paste0("fit <- rasch_mfrm(dat,\n  ", paste(c(
            paste0("person = ", qstr(input$lp_person)),
            paste0("item = ", qstr(input$lp_item)),
            paste0("score = ", qstr(input$lp_score)),
            paste0("facets = ", qvec(input$lp_facets)),
            code_args_common,
            if (!is.null(input$lp_interaction) && input$lp_interaction != NONE)
              paste0("interaction = ", qstr(input$lp_interaction)),
            code_est), collapse = ",\n  "), ")")
          rasch_mfrm(df, person = input$lp_person, item = input$lp_item,
                     score = input$lp_score, facets = input$lp_facets,
                     n_groups = ng, adjust_N = adjN,
                     interaction = if (!is.null(input$lp_interaction) &&
                                       input$lp_interaction != NONE)
                       input$lp_interaction else NULL,
                     maxit = max(5, input$maxit %||% 60),
                     tol = max(1e-12, input$tol %||% 1e-8))
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
                maxit = max(5, input$maxit %||% 60),
                tol = max(1e-12, input$tol %||% 1e-8))
        }
      }, error = function(e) e)
    })
    if (inherits(fit, "error")) {
      showNotification(paste("Analysis failed:", conditionMessage(fit)),
                       type = "error", duration = NULL)
      return(NULL)
    }
    if (!is.null(code_call))
      rcode_str(paste(c("library(rmt)", "", src_line,
                        if (length(code_notes)) c("", code_notes), "",
                        code_call), collapse = "\n"))
    # routine handling notes are informational; only real problems warn
    if (length(fit$notes))
      showNotification(paste(fit$notes, collapse = "\n"), type = "message",
                       duration = 8)
    conv <- if (!is.null(fit$est)) fit$est$converged else fit$converged
    if (!isTRUE(conv))
      showNotification("Estimation did not converge; consider raising the maximum iterations or loosening the convergence criterion.",
                       type = "warning", duration = NULL)
    override_fit(NULL); override_desc(NULL)
    # paired-comparison results live on their own tab; the Rasch tabs
    # suspend while a BTL analysis is current
    if (inherits(fit, "rmt_btl")) {
      btl_fit(fit)
      try(nav_select("nav", "p_btl", session = session), silent = TRUE)
      return(NULL)
    }
    btl_fit(NULL)
    try(nav_select("nav", "p_summary", session = session), silent = TRUE)
    fit
  })
  btl_fit <- reactiveVal(NULL)
  fit <- reactive({
    f <- override_fit()
    if (is.null(f)) f <- analysis()
    req(f); f
  })

  output$has_mc <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    !is.null(f) && !is.null(f$mc)
  })
  outputOptions(output, "has_mc", suspendWhenHidden = FALSE)

  output$sel_item_title <- renderUI(span(class = "fw-semibold",
    tryCatch(sel_item(), error = function(e) "Selected item")))

  # only offer the pages that apply to the current analysis: Facets needs a
  # many-facet fit, Frames an extended-frames fit, BTL a paired-comparison
  # analysis, and Guessing a dichotomous one. Everything else stays.
  # Every nav_panel and nav_menu carries an explicit value, and visibility is
  # driven by those values through the rmt-nav-vis handler (shiny::hideTab,
  # which backs bslib::nav_hide, cannot reach entries inside a nav_menu).
  observe({
    f <- tryCatch(fit(), error = function(e) NULL)
    bf <- btl_fit()
    show <- function(value, on)
      session$sendCustomMessage("rmt-nav-vis",
                                list(value = value, show = isTRUE(on)))
    show("p_facets", inherits(f, "rasch_mfrm"))
    show("p_frames", inherits(f, "rasch_efrm"))
    show("p_btl", !is.null(bf))
    show("p_guess", !is.null(f) && !inherits(f, "rasch_mfrm") &&
           !inherits(f, "rasch_efrm") && max(f$m) == 1L)
    rasch_on <- !is.null(f)
    # judge-factor DIF applies to a paired-comparison fit once judge
    # factors are nominated in the Data roles
    btl_dif_on <- !is.null(bf) && length(input$bt_jfactors) > 0
    for (tgt in c("p_summary", "p_items", "p_persons", "p_targeting",
                  "p_test", "p_dim", "p_equating"))
      show(tgt, rasch_on)
    # Local dependence: the Rasch Q3 suite, or the within-judge dependence
    # analysis of a paired-comparison fit
    show("p_ld", rasch_on || !is.null(bf))
    # DIF needs at least one person factor in the fit (Rasch), or a
    # paired-comparison fit with judge factors nominated
    show("p_dif", (rasch_on && !is.null(f$factors) &&
                     length(names(f$factors)) > 0) || btl_dif_on)
    # menu headers hide too when everything inside them is hidden
    show("menu_structure", rasch_on || !is.null(bf))
    show("menu_invariance", rasch_on || btl_dif_on)
    show("menu_more", rasch_on || !is.null(bf))
  })

  # ------------------------------------------------ UI visibility flags --
  # "nothing to show" areas hide instead of rendering empty; each flag pairs
  # with a conditionalPanel in the UI (same pattern as has_mc)
  output$dif_is_factorial <- reactive(identical(input$dif_factor, FACTORIAL))
  outputOptions(output, "dif_is_factorial", suspendWhenHidden = FALSE)
  output$has_interaction <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    inherits(f, "rasch_mfrm") && !is.null(f$interaction)
  })
  outputOptions(output, "has_interaction", suspendWhenHidden = FALSE)
  output$is_btl <- reactive(!is.null(btl_fit()))
  outputOptions(output, "is_btl", suspendWhenHidden = FALSE)
  # empty-state flags for the run-on-demand cards: before a run the card
  # shows only its controls and one muted line (no table, plot, or download)
  output$has_dep <- reactive(!is.null(dep_res()))
  outputOptions(output, "has_dep", suspendWhenHidden = FALSE)
  output$has_spread <- reactive(!is.null(spread_res()))
  outputOptions(output, "has_spread", suspendWhenHidden = FALSE)
  output$has_difsize <- reactive(!is.null(dif_size_res()))
  outputOptions(output, "has_difsize", suspendWhenHidden = FALSE)
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
    updateSelectInput(session, "dif_item", choices = its, selected = its[1])
    updateSelectizeInput(session, "subtest_items", choices = its, selected = character(0))
    fac <- names(fit()$factors)
    dif_choices <- if (length(fac)) c(fac, FACTORIAL) else character(0)
    updateSelectizeInput(session, "dif_factor", choices = dif_choices,
                         selected = if (length(dif_choices)) dif_choices[1]
                                    else character(0))
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
    contr_res(NULL)
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
      showNotification("Re-analysed with the subtest in place. Reset the override (Data page) or run again to return to the base fit.",
                       type = "message", duration = 8)
    }
  })
  output$subtest_status <- renderUI({
    if (is.null(override_desc())) return(NULL)
    p(class = "text-success small mt-2", paste("Active:", override_desc()))
  })

  observeEvent(input$make_split, {
    f <- fit()
    req(input$dif_item %in% f$items$item,
        !is.null(f$factors), input$dif_factor %in% names(f$factors))
    res <- tryCatch(split_items(f, input$dif_item, by = input$dif_factor),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Split failed:", conditionMessage(res)), type = "error")
    } else {
      override_fit(res)
      override_desc(sprintf("split: item %s by %s", input$dif_item, input$dif_factor))
      showNotification(
        sprintf("Re-analysed with %s split by %s. Reset the override (Data page) or run again to return to the base fit.",
                input$dif_item, input$dif_factor),
        type = "message", duration = 8)
    }
  })

  sel_item <- reactive({
    f <- fit()
    i <- input$items_tbl_rows_selected
    if (length(i)) f$items$item[i] else f$items$item[1]
  })

  # ------------------------------------------------------- plot plumbing --
  # per-output "R code" disclosure: `code` is a function returning the exact
  # rmt call reproducing the output (it may read reactives, so the snippet
  # follows the current selections). Rendering is never suspended: the text
  # must be ready when the collapsed <details> footer is opened.
  register_code <- function(id, code) {
    cid <- paste0(id, "_code")
    output[[cid]] <- renderText(code())
    outputOptions(output, cid, suspendWhenHidden = FALSE)
  }
  register_plot <- function(id, fun, w = 9, h = 6, code = NULL) {
    output[[id]] <- renderPlot(fun(), res = 96)
    if (!is.null(code)) register_code(id, code)
    for (fmt in c("png", "pdf")) local({
      fmt_ <- fmt
      output[[paste0(id, "_", fmt_)]] <- downloadHandler(
        filename = function() paste0("rmt_", id, ".", fmt_),
        content = function(file) {
          # 300 dpi PNG (and vector PDF) for publication
          if (fmt_ == "png") png(file, width = w, height = h, units = "in", res = 300)
          else pdf(file, width = w, height = h)
          fun(); dev.off()
        })
    })
  }
  register_table <- function(id, fun, dt_fun, code = NULL) {
    output[[id]] <- renderDT(dt_fun())
    if (!is.null(code)) register_code(id, code)
    output[[paste0(id, "_csv")]] <- downloadHandler(
      filename = function() paste0("rmt_", id, ".csv"),
      content = function(file) write.csv(fun(), file, row.names = FALSE))
  }
  # APA-leaning DT wrapper: Bootstrap 5 skin, right-aligned numerics, paging
  # controls only when the table needs them. `fit_col` colours fit residuals
  # beyond |2.5| with the theme danger colour; `p_bold` bolds p-values < .05.
  # curated display columns: fit objects carry every statistic, but the
  # tables show a readable core; the per-table "detailed columns" switch
  # reveals the rest (CSV downloads always contain everything)
  CORE <- list(
    items = c("item", "location", "se", "fit_resid", "infit_ms", "outfit_ms",
              "chisq", "df", "p_adj"),
    person = c("id", "raw", "max_raw", "theta", "se", "extreme", "fit_resid"),
    dif = c("factor", "item", "F_uniform", "p_uniform_adj", "eta2_uniform",
            "F_nonuniform", "p_nonuniform_adj", "eta2_nonuniform",
            "uniform_DIF", "nonuniform_DIF"),
    dif_fact = c("item", "term", "F_uniform", "p_uniform_adj", "eta2_uniform",
                 "uniform_DIF", "F_nonuniform", "p_nonuniform_adj",
                 "eta2_nonuniform", "nonuniform_DIF"),
    facet = c("level", "severity", "se", "n", "fit_resid"),
    btl_obj = c("object", "location", "se", "comparisons", "wins", "score",
                "fit_resid"),
    btl_judge = c("judge", "n", "fit_resid", "misfit"),
    equate = c("item", "location_1", "location_2", "adj_difference", "t",
               "p_adj", "drift"),
    dif_size = c("item", "term", "level_a", "level_b", "difference", "lower",
                 "upper", "p_adj", "significant", "practical"),
    contrast = c("item", "contrast", "estimate", "se", "statistic", "p_adj",
                 "significant", "practical"),
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
  num_dt <- function(d, digits = 3, fit_col = NULL, p_bold = NULL,
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
    for (fc in which(orig %in% fit_col))
      dt <- formatStyle(dt, fc, color = styleInterval(
        c(-2.5, 2.5), c("var(--bs-danger)", "inherit", "var(--bs-danger)")))
    for (pc in which(orig %in% p_bold))
      dt <- formatStyle(dt, pc,
                        fontWeight = styleInterval(0.05, c("bold", "normal")))
    dt
  }

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
    f <- override_fit()
    if (is.null(f)) f <- tryCatch(analysis(), error = function(e) NULL)
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

  # test-of-fit and targeting/reliability summaries as statistic/value
  # tables (the package builds them; the CSVs carry the raw column names).
  # The score-to-measure table stays available via score_table(fit) and the
  # everything-ZIP; thresholds live on the Items explorer Thresholds tab.
  register_table("fitsum_tbl", function() fit_summary_table(fit()),
                 function() {
    d <- fit_summary_table(fit())
    names(d) <- c("Statistic", "Value")
    num_dt(d, paging = FALSE)
  }, code = function() "fit_summary_table(fit)")
  output$fitsum_tbl_csv <- downloadHandler(
    filename = function() "fit_summary.csv",
    content = function(file)
      write.csv(fit_summary_table(fit()), file, row.names = FALSE))
  # routine handling notes (the old text panel printed fit$notes)
  output$fitsum_notes <- renderUI({
    f <- fit()
    if (!length(f$notes)) return(NULL)
    sprintf("Note. %s.", paste(f$notes, collapse = "; "))
  })
  register_table("targeting_tbl", function() targeting_table(fit()),
                 function() {
    d <- targeting_table(fit())
    names(d) <- c("Statistic", "Value")
    num_dt(d, paging = FALSE)
  }, code = function() "targeting_table(fit)")
  output$targeting_tbl_csv <- downloadHandler(
    filename = function() "targeting.csv",
    content = function(file)
      write.csv(targeting_table(fit()), file, row.names = FALSE))

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
    filename = function() "rmt_ctt_tbl.csv",
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
                       type = "error", duration = NULL)
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
    d <- curate(fit()$items, "items", full = isTRUE(input$items_full),
                extra = if (length(unique(fit()$m)) > 1) "max")
    dt <- num_dt(d, page_len = 25, selection = "single",
                 fit_col = "fit_resid", p_bold = c("p_adj", "p_anova"))
    # per-statistic misfit highlighting (no single flag column): adjusted
    # chi-square p < .05, and mean squares outside 0.7-1.3 (Wright &
    # Linacre 1994); |fit residual| > 2.5 is handled by fit_col above
    for (j in which(names(d) == "p_adj"))
      dt <- formatStyle(dt, j, color = styleInterval(
        0.05, c("var(--bs-danger)", "inherit")))
    for (j in which(names(d) %in% c("infit_ms", "outfit_ms")))
      dt <- formatStyle(dt, j, color = styleInterval(
        c(0.7, 1.3), c("var(--bs-danger)", "inherit", "var(--bs-danger)")))
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
    filename = function() paste0("rmt_chisq_intervals_", sel_item(), ".csv"),
    content = function(file)
      write.csv(chisq_res()$intervals, file, row.names = FALSE))
  output$chisq_cat_csv <- downloadHandler(
    filename = function() paste0("rmt_chisq_categories_", sel_item(), ".csv"),
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

  # explorer display settings (gear popover): class intervals and scale
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
  for (ext in c("pdf", "zip")) local({
    ext_ <- ext
    output[[paste0("items_all_", ext_)]] <- downloadHandler(
      filename = function() paste0(ex_what(), "_all_items.", ext_),
      content = function(file)
        withProgress(message = "Drawing every item…", value = 0.4,
          save_item_plots(fit(), ex_what(), file, n_groups = ex_ng(),
                          grid = ex_grid(),
                          observed = isTRUE(input$show_obs))))
  })
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
    if ("fit_resid" %in% names(d))
      dt <- formatStyle(dt, "fit_resid", color = styleInterval(
        c(-2.5, 2.5), c("var(--bs-danger)", "inherit", "var(--bs-danger)")))
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
  for (ext in c("pdf", "zip")) local({
    ext_ <- ext
    output[[paste0("kidmap_all_", ext_)]] <- downloadHandler(
      filename = function() paste0("kidmaps.", ext_),
      content = function(file)
        withProgress(message = "Drawing a kidmap for every person…",
                     value = 0.4,
                     save_person_plots(fit(), file, level = kid_level())))
  })
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
  # (the FACTORIAL constant is defined at the top of the file: observers
  # created before this section reference it)
  dif_alpha <- reactive({
    a <- input$dif_alpha
    if (is.null(a) || is.na(a) || a <= 0 || a >= 1) 0.05 else a
  })
  dif_res <- reactive({
    f <- fit(); req(!is.null(f$factors), length(names(f$factors)) > 0)
    dif_anova(f, p_adjust = input$dif_padj %||% "BH", alpha = dif_alpha())
  })
  dif_fact <- reactive({
    f <- fit(); req(!is.null(f$factors), length(names(f$factors)) >= 1)
    dif_anova_factorial(f, p_adjust = input$dif_padj %||% "BH",
                        alpha = dif_alpha(),
                        effects = input$dif_effects %||% "factorial")
  })
  register_table("dif_tbl", function() {
    if (identical(input$dif_factor, FACTORIAL)) dif_fact()$summary else dif_res()
  }, function() {
    if (identical(input$dif_factor, FACTORIAL)) {
      d <- curate(dif_fact()$summary, "dif_fact",
                  full = isTRUE(input$dif_full))
      d$uniform_DIF <- ifelse(d$uniform_DIF, "*", "")
      d$nonuniform_DIF <- ifelse(d$nonuniform_DIF, "*", "")
      if ("superseded" %in% names(d))
        d$superseded <- ifelse(d$superseded, "(superseded)", "")
      num_dt(d)
    } else {
      d <- dif_res()
      if (!is.null(input$dif_factor) && input$dif_factor %in% d$factor)
        d <- d[d$factor == input$dif_factor, ]
      d <- curate(d, "dif", full = isTRUE(input$dif_full))
      d$uniform_DIF <- ifelse(d$uniform_DIF, "*", "")
      d$nonuniform_DIF <- ifelse(d$nonuniform_DIF, "*", "")
      num_dt(d)
    }
  }, code = function() {
    if (identical(input$dif_factor, FACTORIAL))
      sprintf('dif_anova_factorial(fit, effects = "%s", p_adjust = "%s", alpha = %s)$summary',
              input$dif_effects %||% "factorial", input$dif_padj %||% "BH",
              dif_alpha())
    else
      sprintf('dif_anova(fit, p_adjust = "%s", alpha = %s)',
              input$dif_padj %||% "BH", dif_alpha())
  })
  # full per-item ANOVA table, computed lazily when its disclosure is first
  # switched on (the DT renders only once visible): factorial -> the joint
  # model's terms; single-factor mode -> one factorial fit per nominated
  # factor, stacked with a factor column
  dif_full_dat <- reactive({
    f <- fit()
    if (identical(input$dif_factor, FACTORIAL)) return(dif_fact()$terms)
    req(!is.null(f$factors), length(names(f$factors)) > 0)
    do.call(rbind, lapply(names(f$factors), function(fc)
      cbind(factor = fc,
            dif_anova_factorial(f, factors = fc,
                                p_adjust = input$dif_padj %||% "BH",
                                alpha = dif_alpha())$terms)))
  })
  register_table("dif_full_tbl", function() dif_full_dat(), function() {
    d <- dif_full_dat()
    d$significant <- ifelse(d$significant, "*", "")
    d$superseded <- ifelse(d$superseded, "(superseded)", "")
    num_dt(d)
  }, code = function() {
    if (identical(input$dif_factor, FACTORIAL))
      sprintf('dif_anova_factorial(fit, effects = "%s", p_adjust = "%s", alpha = %s)$terms',
              input$dif_effects %||% "factorial", input$dif_padj %||% "BH",
              dif_alpha())
    else
      sprintf('do.call(rbind, lapply(names(fit$factors), function(f)\n  cbind(factor = f,\n        dif_anova_factorial(fit, factors = f, p_adjust = "%s", alpha = %s)$terms)))',
              input$dif_padj %||% "BH", dif_alpha())
  })
  output$dif_note <- renderUI({
    if (identical(input$dif_factor, FACTORIAL)) {
      fr <- dif_fact()
      d <- fr$terms
      sup <- sum(d$superseded, na.rm = TRUE)
      sprintf("Note. %d of %d terms significant after adjustment%s. Class intervals: %d, set from the smallest factor-combination cell (about 30 responses per interval-by-cell count).",
              sum(d$significant, na.rm = TRUE), nrow(d),
              if (sup) sprintf(" (%d superseded by an interaction)", sup)
              else "", fr$n_groups %||% NA_integer_)
    } else {
      d <- dif_res()
      ng <- attr(d, "n_groups")
      parts <- vapply(split(d, d$factor), function(g)
        sprintf("%s: %d uniform, %d non-uniform", g$factor[1],
                sum(g$uniform_DIF, na.rm = TRUE),
                sum(g$nonuniform_DIF, na.rm = TRUE)), "")
      ci_txt <- if (!is.null(ng))
        sprintf(" Class intervals per factor: %s (set from each factor's smallest group).",
                paste(names(ng), ng, sep = " = ", collapse = ", "))
      else ""
      paste0("Note. Items flagged per factor - ",
             paste(parts, collapse = "; "), ".", ci_txt)
    }
  })
  register_table("dif_tukey_tbl", function() dif_fact()$tukey, function() {
    validate(need(identical(input$dif_factor, FACTORIAL),
                  "Choose the factorial option in the sidebar to see Tukey HSD comparisons."))
    tk <- dif_fact()$tukey
    if (!nrow(tk))
      return(datatable(data.frame(note = "no significant group terms to compare"),
                       rownames = FALSE, style = "bootstrap5",
                       class = "table-sm compact",
                       options = list(dom = "t")))
    num_dt(tk)
  }, code = function() "dif_anova_factorial(fit)$tukey")
  # in factorial mode the graphical display uses the factor-combination
  # cells; plot_icc accepts several factor names for exactly this
  dif_icc_group <- function(f) {
    if (identical(input$dif_factor, FACTORIAL)) names(f$factors)
    else { req(input$dif_factor %in% names(f$factors)); input$dif_factor }
  }
  register_plot("dif_icc", function() {
    f <- fit()
    req(input$dif_item %in% f$items$item, !is.null(f$factors))
    plot_icc(f, input$dif_item, group = dif_icc_group(f))
  }, code = function() {
    g <- if (identical(input$dif_factor, FACTORIAL))
      paste0('c("', paste(names(fit()$factors), collapse = '", "'), '")')
    else sprintf('"%s"', input$dif_factor %||% "")
    sprintf('plot_icc(fit, "%s", group = %s)', input$dif_item %||% "", g)
  })

  # DIF magnitude in logits: single factor -> the selected item and factor;
  # factorial -> sizes for every significant, non-superseded group term
  dif_size_res <- reactiveVal(NULL)
  observeEvent(input$dif_size_go, {
    f <- fit()
    flg <- max(0.05, input$dif_size_flag %||% 0.5)
    mn <- max(2, input$dif_size_minn %||% 20)
    res <- tryCatch({
      if (identical(input$dif_factor, FACTORIAL)) {
        fa <- dif_anova_factorial(f, sizes = TRUE,
                                  p_adjust = input$dif_padj %||% "BH",
                                  alpha = dif_alpha(),
                                  effects = input$dif_effects %||% "factorial")
        sz <- fa$sizes
        if (is.null(sz) || !nrow(sz))
          stop("no significant, non-superseded group terms to size")
        sz$practical <- abs(sz$difference) >= flg
        sz
      } else {
        req(input$dif_item %in% f$items$item,
            input$dif_factor %in% names(f$factors))
        ds <- dif_size(f, input$dif_item, by = input$dif_factor,
                       flag_logits = flg, min_n = mn)
        cbind(item = ds$item, term = ds$by, ds$pairs)
      }
    }, error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("DIF size:", conditionMessage(res)),
                       type = "warning")
      dif_size_res(NULL)
    } else dif_size_res(res)
  })
  output$dif_size_tbl <- DT::renderDT({
    d <- dif_size_res()
    validate(need(!is.null(d), "Compute DIF sizes to see the logit-scale comparisons."))
    d <- curate(d, "dif_size", full = isTRUE(input$difsize_full))
    if ("significant" %in% names(d))
      d$significant <- ifelse(d$significant, "*", "")
    if ("practical" %in% names(d))
      d$practical <- ifelse(d$practical, "PRACTICAL", "")
    num_dt(d)
  })
  register_code("dif_size_tbl", function() {
    if (identical(input$dif_factor, FACTORIAL))
      "dif_anova_factorial(fit, sizes = TRUE)$sizes"
    else sprintf('dif_size(fit, "%s", by = "%s")',
                 input$dif_item %||% "", input$dif_factor %||% "")
  })
  output$dl_dif_size <- downloadHandler(
    filename = function() "dif_sizes.csv",
    content = function(file) {
      d <- dif_size_res()
      if (is.null(d)) stop("compute DIF sizes first")
      write.csv(d, file, row.names = FALSE)
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
    if ("significant" %in% names(d))
      d$significant <- ifelse(d$significant, "*", "")
    if ("practical" %in% names(d))
      d$practical <- ifelse(d$practical, "PRACTICAL", "")
    if ("within" %in% names(d)) d$within <- ifelse(d$within, "*", "")
    num_dt(d)
  }, code = function() {
    its <- input$pc_items
    sprintf('dif_contrasts(fit%s%s)',
            if (length(its))
              paste0(', items = c("', paste(its, collapse = '", "'), '")')
            else "",
            if (nzchar(input$pc_id %||% ""))
              paste0(', id = "', input$pc_id, '"') else "")
  })
  # the conventional register_table CSV name is overridden: the download is
  # the full contrast table under the function's own name
  output$contr_tbl_csv <- downloadHandler(
    filename = function() "dif_contrasts.csv",
    content = function(file) {
      r <- contr_res()
      if (is.null(r)) stop("derive and test the contrasts first")
      write.csv(r$table, file, row.names = FALSE)
    })

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
  output$btl_summary <- renderText({
    f <- bfit()
    paste0(sprintf("Conditional ML: %s in %d iterations; sandwich SEs%s.\n",
                   if (f$converged) "converged" else "NOT converged",
                   f$iterations,
                   if (f$clustered) " clustered by judge" else ""),
           sprintf("Pairwise chi-square %.2f on %d df, p = %s.\n",
                   f$total_chisq, f$total_df, fmt_p(f$total_p)),
           sprintf("Log-likelihood %.3f.", f$loglik),
           if (length(f$notes)) paste0("\nNotes: ",
                                       paste(f$notes, collapse = "; ")) else "")
  })
  register_table("btl_obj_tbl", function() bfit()$objects,
                 function() num_dt(curate(bfit()$objects, "btl_obj",
                                          full = isTRUE(input$btl_full)),
                                   fit_col = "fit_resid"),
                 code = function() "# bt from the Data page\nbt$objects")
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
    num_dt(curate(d, "btl_judge", full = isTRUE(input$btl_full)),
           fit_col = "fit_resid")
  }, code = function() "bt$judges")
  register_plot("btl_plot", function() plot_btl(bfit()),
                code = function() "# bt from the Data page\nplot_btl(bt)")
  # object characteristic curve: model expected response against opponent
  # location with per-opponent observed means (dichotomous and graded fits)
  observeEvent(btl_fit(), {
    b <- btl_fit()
    if (!is.null(b)) {
      updateSelectizeInput(session, "btl_occ_obj",
                           choices = b$objects$object,
                           selected = b$objects$object[1])
      updateSelectizeInput(session, "bdif_obj",
                           choices = b$objects$object,
                           selected = b$objects$object[1])
      jf <- input$bt_jfactors
      updateSelectizeInput(session, "bdif_factor",
                           choices = if (length(jf)) jf else character(0),
                           selected = if (length(jf)) jf[1] else character(0))
    }
    # judge-group DIF results belong to the fit they came from
    bdif_res(NULL)
  }, ignoreNULL = FALSE)
  register_plot("btl_occ", function() {
    req(input$btl_occ_obj %in% bfit()$objects$object)
    plot_btl_icc(bfit(), input$btl_occ_obj)
  }, w = 8, h = 5.5, code = function()
    paste0(sprintf('plot_btl_icc(bt, "%s")', input$btl_occ_obj %||% ""),
           "\n# all objects: one page each in the PDF download"))
  output$btl_occ_all_pdf <- downloadHandler(
    filename = function() "occ_all_objects.pdf",
    content = function(file) {
      b <- bfit()
      pdf(file, width = 8, height = 5.5, onefile = TRUE)
      on.exit(dev.off(), add = TRUE)
      for (o in b$objects$object)
        tryCatch(plot_btl_icc(b, o), error = function(e) NULL)
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
  output$has_btl_dep <- reactive({
    b <- btl_fit()
    !is.null(b) && !is.null(b$dependence)
  })
  outputOptions(output, "has_btl_dep", suspendWhenHidden = FALSE)
  register_table("btl_dep_tbl", function() {
    validate(need(!is.null(bfit()$dependence),
                  "Nominate a judgment-order column in the Data roles to estimate within-judge dependence."))
    bfit()$dependence
  }, function() {
    validate(need(!is.null(bfit()$dependence),
                  "Nominate a judgment-order column in the Data roles to estimate within-judge dependence."))
    num_dt(bfit()$dependence, p_bold = "p")
  }, code = function() "bt$dependence")

  # -------------------------------------------- BTL DIF by judge group --
  # the judge grouping handed to btl_dif() / plot_btl_icc(): the nominated
  # factor's value on each judge's first row, named by judge
  bdif_groups_build <- function() {
    if (is.null(btl_fit()))
      stop("run a paired-comparisons (BTL) analysis first")
    fc <- input$bdif_factor
    if (is.null(fc) || !nzchar(fc))
      stop("choose a judge factor in the sidebar")
    df <- raw_data()
    jc <- input$bt_judge
    if (is.null(jc) || identical(jc, NONE) || !jc %in% names(df))
      stop("judge-group DIF needs the judge column nominated on the Data page")
    if (!fc %in% names(df)) stop("column not found: ", fc)
    jd <- as.character(df[[jc]])
    first <- !duplicated(jd)
    setNames(as.character(df[[fc]])[first], jd[first])
  }
  # the reproducible-code line building the same grouping from the data
  bdif_code_grp <- function()
    sprintf("grp <- setNames(as.character(dat$%s), dat$%s)[!duplicated(dat$%s)]",
            input$bdif_factor %||% "factor", input$bt_judge %||% "judge",
            input$bt_judge %||% "judge")
  bdif_alpha <- reactive({
    a <- input$bdif_alpha
    if (is.null(a) || is.na(a) || a <= 0 || a >= 1) 0.05 else a
  })
  # judge factors nominated (or changed) after the run still reach the
  # factor select; the results themselves are computed on request only
  observeEvent(input$bt_jfactors, {
    jf <- input$bt_jfactors
    sel <- if (!is.null(input$bdif_factor) && input$bdif_factor %in% jf)
      input$bdif_factor else if (length(jf)) jf[1] else character(0)
    updateSelectizeInput(session, "bdif_factor",
                         choices = if (length(jf)) jf else character(0),
                         selected = sel)
  }, ignoreNULL = FALSE)
  bdif_res <- reactiveVal(NULL)
  observeEvent(input$bdif_run, {
    r <- withProgress(message = "Resolving objects by judge group…",
                      value = 0.4,
                      tryCatch(btl_dif(btl_fit(),
                                       groups = bdif_groups_build(),
                                       alpha = bdif_alpha()),
                               error = function(e) e))
    if (inherits(r, "error")) {
      showNotification(paste("DIF analysis failed:", conditionMessage(r)),
                       type = "error", duration = NULL)
      bdif_res(NULL)
    } else bdif_res(r)
  })
  register_table("bdif_anova_tbl", function() {
    r <- bdif_res(); req(!is.null(r)); r$anova
  }, function() {
    r <- bdif_res()
    validate(need(!is.null(r),
                  "Choose a judge factor in the sidebar and run the DIF analysis."))
    d <- r$anova[, intersect(c("object", "n", "F_uniform", "p_uniform_adj",
                               "uniform_DIF", "F_nonuniform",
                               "p_nonuniform_adj", "nonuniform_DIF"),
                             names(r$anova)), drop = FALSE]
    d$uniform_DIF <- ifelse(d$uniform_DIF, "*", "")
    d$nonuniform_DIF <- ifelse(d$nonuniform_DIF, "*", "")
    num_dt(d)
  }, code = function()
    paste0("# dat: the comparison data; bt: the fit from the Data page\n",
           bdif_code_grp(), "\n",
           sprintf("btl_dif(bt, groups = grp, alpha = %s)$anova",
                   bdif_alpha())))
  register_table("bdif_sizes_tbl", function() {
    r <- bdif_res(); req(!is.null(r), !is.null(r$sizes)); r$sizes
  }, function() {
    r <- bdif_res()
    validate(need(!is.null(r),
                  "Choose a judge factor in the sidebar and run the DIF analysis."))
    validate(need(!is.null(r$sizes),
                  "No object could be resolved by this grouping (see the notes on the analysis-of-variance panel)."))
    d <- r$sizes[, intersect(c("object", "level_a", "level_b", "difference",
                               "se", "z", "p_adj", "significant",
                               "practical"), names(r$sizes)), drop = FALSE]
    d$significant <- ifelse(d$significant, "*", "")
    d$practical <- ifelse(d$practical, "PRACTICAL", "")
    num_dt(d)
  }, code = function()
    paste0(bdif_code_grp(), "\n",
           sprintf("btl_dif(bt, groups = grp, alpha = %s)$sizes",
                   bdif_alpha())))
  output$bdif_notes <- renderUI({
    r <- bdif_res()
    if (is.null(r) || !length(r$notes)) return(NULL)
    sprintf("Note. %s.", paste(r$notes, collapse = "; "))
  })
  register_plot("bdif_occ", function() {
    b <- bfit()
    req(input$bdif_obj %in% b$objects$object)
    grp <- tryCatch(bdif_groups_build(), error = function(e) NULL)
    validate(need(!is.null(grp),
                  "Choose a judge factor in the sidebar to overlay the per-group observed means."))
    plot_btl_icc(b, input$bdif_obj, group = grp)
  }, w = 8, h = 5.5, code = function()
    paste0(bdif_code_grp(), "\n",
           sprintf('plot_btl_icc(bt, "%s", group = grp)',
                   input$bdif_obj %||% "")))

  # ---------------------------------------------------------------- facets --
  facet_dat <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a long-format (many-facet) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    f$facet_effects[[input$facet_sel]]
  })
  register_table("facet_tbl", function() facet_dat(), function() {
    d <- curate(facet_dat(), "facet", full = isTRUE(input$facets_full))
    datatable(d, rownames = FALSE, style = "bootstrap5",
              class = "table-sm compact hover order-column",
              options = list(pageLength = 15, scrollX = TRUE,
                             dom = if (nrow(d) > 15) "tip" else "t")) |>
      formatRound(setdiff(names(d)[vapply(d, is.numeric, TRUE)], "n"), 3)
  }, code = function()
    sprintf('fit$facet_effects[["%s"]]', input$facet_sel %||% ""))
  register_plot("facet_plot", function() {
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a long-format (many-facet) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    plot_facets(f, input$facet_sel)
  }, code = function()
    sprintf('plot_facets(fit, "%s")', input$facet_sel %||% ""))
  facet_int <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm") && !is.null(f$interaction),
                  "Run a long-format analysis with an item-by-facet interaction selected."))
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
    filename = function() format(Sys.time(), "rmt_anchors_%Y%m%d_%H%M.csv"),
    content = function(file) {
      f <- fit()
      thr <- f$thresholds
      write.csv(data.frame(item = f$items$item[thr$item], k = thr$k,
                           tau = thr$tau), file, row.names = FALSE)
    })

  output$dl_calib <- downloadHandler(
    filename = function() format(Sys.time(), "rmt_calibration_%Y%m%d_%H%M.csv"),
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
  observeEvent(input$run, dim_subsets(NULL), priority = 9)
  observeEvent(input$dim_apply, {
    if (length(input$dim_pos) >= 2 && length(input$dim_neg) >= 2) {
      if (length(intersect(input$dim_pos, input$dim_neg))) {
        showNotification("The two subsets must be disjoint.", type = "error")
      } else dim_subsets(list(pos = input$dim_pos, neg = input$dim_neg))
    } else if (!length(input$dim_pos) && !length(input$dim_neg)) {
      dim_subsets(NULL)
      showNotification("Reset to the first-contrast split.", type = "message")
    } else {
      showNotification("Nominate at least two items in each subset (or leave both empty).",
                       type = "warning")
    }
  })
  dim_res <- reactive({
    s <- dim_subsets()
    if (is.null(s)) dimensionality_test(fit())
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
    if (is.null(s)) return("dimensionality_test(fit)")
    sprintf("dimensionality_test(fit, items_positive = %s,\n  items_negative = %s)",
            qvec(s$pos), qvec(s$neg))
  })

  # magnitude of multidimensionality (Andrich 2016): needs every item in a
  # subset; the PC1 split satisfies this by construction, manual subsets may not
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
        "Adjust subsets A and B (or leave both empty for the PC1 split) and re-run the t-test first."),
        type = "warning", duration = 10)
      return()
    }
    r <- withProgress(message = "Subtest re-analysis…", value = 0.4,
                      tryCatch(dimensionality_magnitude(f, list(s$pos, s$neg)),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Magnitude estimate failed:", conditionMessage(r)),
                       type = "error", duration = NULL)
    else dm_res(r)
  })
  output$dm_tbl <- renderDT({
    r <- dm_res()
    validate(need(!is.null(r),
                  "Press the button; the current subsets (manual, or the PC1 split) are combined into super-items and the two reliability calculations compared."))
    num_dt(r$table)
  })
  output$dm_tbl_csv <- downloadHandler(
    filename = function() "rmt_dimensionality_magnitude.csv",
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
    d <- residual_pca(fit())$loadings_matrix
    if (!isTRUE(input$load_full))
      d <- d[, intersect(c("item", "PC1", "PC2", "PC3"), names(d)),
             drop = FALSE]
    num_dt(d)
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
  register_plot("pca_plot", function() plot_pca(fit()),
                code = function() "plot_pca(fit)")

  # -------------------------------------------------------- local dependence --
  register_plot("rcor", function() plot_resid_cor(fit()), w = 8, h = 8,
                code = function() "plot_resid_cor(fit)")
  ld_flag <- reactive({
    fl <- input$ld_flag
    if (is.null(fl) || is.na(fl) || fl <= 0) 0.2 else fl
  })
  ld_res <- reactive(residual_correlations(fit(), flag = ld_flag()))
  # Yen's Q3 for every item pair, sorted by Q3; Q3* is the excess over the
  # average off-diagonal Q3, and pairs above the flag threshold are starred
  output$rpairs_tbl <- renderDT({
    d <- ld_res()$pairs
    d$flagged <- ifelse(d$flagged, "*", "")
    num_dt(d)
  })
  register_code("rpairs_tbl", function()
    sprintf("residual_correlations(fit, flag = %s)$pairs", ld_flag()))
  output$rpairs_tbl_csv <- downloadHandler(
    filename = function() "q3_statistics.csv",
    content = function(file)
      write.csv(ld_res()$pairs, file, row.names = FALSE))
  output$rpairs_note <- renderUI(
    sprintf("Average off-diagonal Q3 %.3f; pairs with Q3* above %.1f are flagged (Yen 1984; Christensen, Makransky & Horton 2017).",
            ld_res()$average, ld_flag()))

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
                       type = "error", duration = NULL)
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
    filename = function() "rmt_dependence_thresholds.csv",
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
                       type = "error", duration = NULL)
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
    filename = function() "rmt_spread_test.csv",
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
    ch <- input$guess_chance
    if (is.null(ch) || is.na(ch) || ch <= 0 || ch >= 1) ch <- 0.25
    anc <- if (length(input$guess_anchors)) input$guess_anchors else NULL
    r <- withProgress(message = "Tailored analysis (three re-analyses)…",
                      value = 0.3,
                      tryCatch(tailored_analysis(f, chance = ch,
                                                 anchor_items = anc),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Tailored analysis failed:", conditionMessage(r)),
                       type = "error", duration = NULL)
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
  }, code = function() {
    ch <- input$guess_chance
    if (is.null(ch) || is.na(ch) || ch <= 0 || ch >= 1) ch <- 0.25
    sprintf("ta <- tailored_analysis(fit, chance = %s)\nta$table", ch)
  })
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
    updateSelectInput(session, "cmp_ref", choices = names(k),
                      selected = if (!is.null(input$cmp_ref) &&
                                     input$cmp_ref %in% names(k))
                        input$cmp_ref else names(k)[1])
    updateSelectInput(session, "eq_kept", choices = names(k),
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
    f <- override_fit()
    if (is.null(f)) f <- tryCatch(analysis(), error = function(e) NULL)
    if (is.null(f)) {
      showNotification("Run an analysis first, then download the report.",
                       type = "warning", duration = 8)
      stop("no fit to report")
    }
    withProgress(message = "Building the HTML report…", value = 0.4,
                 report_html(f, file))
  }
  output$dl_report <- downloadHandler(
    filename = function() "rmt_report.html", content = report_content)
  output$dl_report_nav <- downloadHandler(
    filename = function() "rmt_report.html", content = report_content)

  output$dl_zip <- downloadHandler(
    filename = function() format(Sys.time(), "rmt_results_%Y%m%d_%H%M.zip"),
    content = function(file) {
      f <- fit()
      tmp <- file.path(tempdir(), paste0("rmt_", as.integer(Sys.time())))
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
