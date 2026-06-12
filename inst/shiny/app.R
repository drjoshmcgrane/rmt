# RaschR Shiny GUI
# ---------------------------------------------------------------------------
# A modern bslib interface to the full RaschR analysis: data upload with ID,
# person-factor, and item column nomination; pairwise conditional ML
# estimation (Andrich & Luo 2003); the complete test-of-fit suite;
# every diagnostic plot with per-plot PNG and PDF downloads; and one-click
# export of all tables and plots as a ZIP archive.
# Launch with RaschR::run_app(), or shiny::runApp() from this folder.
# ---------------------------------------------------------------------------
library(shiny)
library(bslib)
library(DT)

if (requireNamespace("RaschR", quietly = TRUE)) {
  library(RaschR)
} else {
  rdir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  if (dir.exists(rdir)) {
    for (f in list.files(rdir, "\\.R$", full.names = TRUE)) source(f)
  } else stop("Install RaschR, or run the app from inst/shiny in the source tree")
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

theme <- bs_theme(
  version = 5, bg = "#f8fafc", fg = "#0f172a",
  primary = "#2563eb", secondary = "#64748b",
  success = "#0f766e", danger = "#dc2626", warning = "#f59e0b",
  "navbar-bg" = "#0f172a", "card-border-color" = "#e2e8f0",
  "border-radius" = "0.65rem", "font-size-base" = "0.95rem"
)

css <- HTML("
  .card { box-shadow: 0 1px 3px rgba(15,23,42,.08); }
  .card-header { background: #fff; font-weight: 600; border-bottom: 1px solid #e2e8f0; }
  .navbar-brand { font-weight: 700; letter-spacing: .02em; }
  .value-box-title { font-size: .72rem; text-transform: uppercase; letter-spacing: .04em; white-space: nowrap; }
  .value-box-value { font-size: 1.45rem; }
  pre, .shiny-text-output { white-space: pre-wrap; font-size: .82rem; }
  .btn-xs { padding: .1rem .5rem; font-size: .75rem; }
  .form-label { font-weight: 600; font-size: .85rem; }
  table.dataTable { font-size: .85rem; }
")

# Card with a plot and PNG/PDF download buttons in the header.
plotCard <- function(id, title, height = "430px") {
  card(
    full_screen = TRUE,
    card_header(div(class = "d-flex justify-content-between align-items-center",
      span(title),
      div(class = "btn-group",
          downloadButton(paste0(id, "_png"), "PNG", class = "btn-outline-secondary btn-xs"),
          downloadButton(paste0(id, "_pdf"), "PDF", class = "btn-outline-secondary btn-xs")))),
    card_body(plotOutput(id, height = height), padding = 8)
  )
}

tableCard <- function(id, title, note = NULL) {
  card(
    full_screen = TRUE,
    card_header(div(class = "d-flex justify-content-between align-items-center",
      span(title),
      downloadButton(paste0(id, "_csv"), "CSV", class = "btn-outline-secondary btn-xs"))),
    card_body(if (!is.null(note)) p(class = "text-muted small mb-2", note),
              DTOutput(id), padding = 12)
  )
}

ui <- page_navbar(
  title = span("RaschR"),
  theme = theme,
  header = tags$head(tags$style(css)),

  # ----------------------------------------------------------------- DATA --
  nav_panel("Data",
    layout_sidebar(
      sidebar = sidebar(width = 330,
        h6("Data source"),
        fileInput("file", NULL, accept = c(".csv", ".txt", ".tsv"),
                  buttonLabel = "Browse…", placeholder = "CSV / TSV file"),
        selectInput("demo_choice", "Or pick an example dataset",
                    c("(none)" = "none",
                      "Multiple choice, dichotomous" = "dich",
                      "Polytomous (PCM)" = "pcm",
                      "Rating scale (RSM)" = "rsm",
                      "Ratings, long format (MFRM)" = "mfrm",
                      "Item sets x groups (EFRM)" = "efrm")),
        radioButtons("model_type", "Model",
                     c("Dichotomous" = "dich",
                       "Partial credit (PCM)" = "pcm",
                       "Rating scale (RSM)" = "rsm",
                       "Many-facet (MFRM)" = "mfrm",
                       "Extended frames (EFRM)" = "efrm")),
        hr(),
        conditionalPanel("['dich','pcm','rsm'].indexOf(input.model_type) > -1",
          h6("Column roles"),
          selectInput("id_col", "ID variable", NONE),
          selectizeInput("factor_cols", "Person factors (DIF groups)", NULL,
                         multiple = TRUE,
                         options = list(placeholder = "none selected")),
          selectizeInput("item_cols", "Item columns", NULL, multiple = TRUE,
                         options = list(placeholder = "all remaining columns")),
          hr(),
          fileInput("key_file", "Multiple-choice key (CSV: item,key)",
                    accept = ".csv", placeholder = "optional"),
          fileInput("anchor_file", "Anchors for equating (CSV: item,k,tau)",
                    accept = ".csv", placeholder = "optional"),
          radioButtons("anchor_type", "Anchor as",
                       c("Individual thresholds" = "individual",
                         "Average item locations" = "average"),
                       inline = TRUE),
          p(class = "text-muted small mt-1",
            "Anchors match by item name; rows for items not present are ignored. Individual anchoring fixes each listed threshold; average anchoring fixes each item's mean location (thresholds stay free). Save an anchor file from the Items page of a previous analysis.")
        ),
        conditionalPanel("input.model_type == 'efrm'",
          h6("Column roles"),
          selectInput("ef_id", "ID variable", NONE),
          selectInput("ef_group", "Person group column", NONE),
          selectizeInput("ef_items", "Item columns", NULL, multiple = TRUE,
                         options = list(placeholder = "all remaining columns")),
          fileInput("ef_sets", "Item-set map (CSV: item,set)",
                    accept = ".csv", placeholder = "optional"),
          checkboxInput("ef_prefix", "Infer sets from item-name prefix", TRUE),
          p(class = "text-muted small",
            "Each item-set by group cell is a frame with its own unit. Group units come from the person-free pairwise comparisons; set units from persons common to the sets.")
        ),
        conditionalPanel("input.model_type == 'mfrm'",
          h6("Column roles (one row per response)"),
          selectInput("lp_person", "Person column", NONE),
          selectInput("lp_item", "Item column", NONE),
          selectInput("lp_score", "Score column", NONE),
          selectizeInput("lp_facets", "Facet columns (e.g. rater)", NULL,
                         multiple = TRUE,
                         options = list(placeholder = "choose at least one")),
          selectInput("lp_interaction", "Item-by-facet interaction (optional)", NONE),
          p(class = "text-muted small",
            "Each item x facet combination is calibrated jointly; facet severities are reported with SEs and fit. An interaction lets one facet be more or less severe on particular items.")
        ),
        sliderInput("ng", "Class intervals", min = 2, max = 16, value = 8),
        h6("Estimation"),
        numericInput("maxit", "Maximum iterations", value = 60, min = 5, step = 5),
        numericInput("tol", "Convergence criterion", value = 1e-8,
                     min = 1e-12, step = 1e-8),
        actionButton("run", "Run analysis", class = "btn-primary w-100 btn-lg mt-2"),
        p(class = "text-muted small mt-3",
          "Estimation: pairwise conditional maximum likelihood (Andrich & Luo 2003).",
          "Person measures: Warm weighted likelihood.")
      ),
      card(card_header("Data preview"),
           card_body(uiOutput("data_info"), DTOutput("preview"), padding = 12))
    )
  ),

  # -------------------------------------------------------------- SUMMARY --
  nav_panel("Summary",
    uiOutput("vboxes"),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(card_header("Test of fit"), card_body(verbatimTextOutput("fit_summary"))),
      card(card_header("Targeting & reliability"), card_body(verbatimTextOutput("targeting")))
    ),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      tableCard("score_tbl", "Score-to-measure table",
                "WLE location and SE for every raw score (complete responders)."),
      tableCard("thr_tbl", "Thresholds with standard errors")
    )
  ),

  # ---------------------------------------------------------------- ITEMS --
  nav_panel("Items",
    div(class = "mb-2 d-flex align-items-end gap-3 flex-wrap",
        numericInput("adjN",
                     "Adjust the item-trait chi-square to a reference sample size (blank = off)",
                     value = NA, min = 50, width = "420px"),
        downloadButton("dl_anchors", "Save anchors (CSV: item,k,tau)",
                       class = "btn-outline-secondary mb-3")),
    tableCard("items_tbl", "Item statistics",
              "Click a row to inspect that item's curves below. Location and SE from the pairwise conditional likelihood; fit residual ~ N(0,1) under fit; item-trait chi-square over class intervals; misfit flag uses BH-adjusted probabilities."),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      plotCard("icc", "Item characteristic curve"),
      plotCard("ccc", "Category probability curves")),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      plotCard("tpc", "Threshold probability curves"),
      plotCard("cfreq", "Category frequencies")),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      tableCard("distractor_tbl", "Distractor analysis",
                "Multiple-choice analyses only (provide a key). Locations use the rest measure; a distractor whose takers are abler than the keyed option's flags a possible miskey."),
      plotCard("distractor_plot", "Option curves"))
  ),

  # -------------------------------------------------------------- PERSONS --
  nav_panel("Persons",
    tableCard("person_tbl", "Person estimates",
              "Warm WLE location and SE per person, with raw score, fit statistics, and your ID and factor columns."),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      plotCard("pfit", "Person fit"),
      plotCard("pim_p", "Person-item threshold distribution"))
  ),

  # ----------------------------------------------------------------- TEST --
  nav_panel("Test plots",
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      plotCard("thrmap", "Threshold map"),
      plotCard("imap", "Item map: location by fit residual")),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      plotCard("tcc", "Test characteristic curve"),
      plotCard("tif", "Test information & SEM")),
    plotCard("guttman", "Guttman scalogram", height = "560px")
  ),

  # ------------------------------------------------------------------ DIF --
  nav_panel("DIF",
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectInput("dif_factor", "Person factor", NONE),
        selectInput("dif_item", "Item for ICC by group", NONE),
        p(class = "text-muted small",
          "ANOVA of standardised residuals: factor effects = uniform DIF; factor x class-interval terms = non-uniform DIF. Probabilities are BH-adjusted across items. With several factors, choose the factorial option to model them jointly: significant interactions supersede their main effects, and Tukey HSD compares the levels of each significant group term."),
        hr(),
        actionButton("make_split", "Resolve: split this item by this factor",
                     class = "btn-outline-primary w-100"),
        p(class = "text-muted small mt-2",
          "Replaces the selected item with one item per group level (each level keeps only its own responses) and re-analyses; the split locations quantify the DIF.")),
      tableCard("dif_tbl", "DIF analysis of variance"),
      tableCard("dif_tukey_tbl", "Tukey HSD comparisons",
                "Pairwise level comparisons for significant, non-superseded group terms (factorial mode)."),
      plotCard("dif_icc", "ICC by group (DIF plot)")
    )
  ),

  # --------------------------------------------------------------- FACETS --
  nav_panel("Facets",
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectInput("facet_sel", "Facet", NONE),
        p(class = "text-muted small",
          "Severities from the joint calibration (positive = more severe). Pooled fit residuals beyond +/-2.5 flag inconsistent levels. Long-format analyses only.")),
      tableCard("facet_tbl", "Facet severities and fit"),
      plotCard("facet_plot", "Severity caterpillar plot"),
      tableCard("facet_int_tbl", "Item-by-facet interactions",
                "Shown when the analysis was run in interactive facet mode; gamma is the extra severity of a level on a particular item.")
    )
  ),

  # -------------------------------------------------------------- EQUATING --
  nav_panel("Equating",
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        fileInput("eq_file", "Reference calibration (CSV: item,location,se)",
                  accept = ".csv"),
        radioButtons("eq_shift", "Scale alignment",
                     c("Allow a shift between origins" = "mean",
                       "Compare raw locations (anchored scales)" = "none")),
        downloadButton("dl_calib", "Save current calibration (CSV)",
                       class = "btn-outline-secondary w-100"),
        p(class = "text-muted small mt-2",
          "Common items (matched by name) are tested against the shifted identity line; flagged items show drift and weaken the equating link. Save a calibration now to equate a future analysis against it.")),
      tableCard("eq_tbl", "Common-item comparison"),
      plotCard("eq_plot", "Equating plot")
    )
  ),

  # --------------------------------------------------------------- FRAMES --
  nav_panel("Frames",
    layout_sidebar(
      sidebar = sidebar(width = 290, open = "always",
        selectInput("frame_item", "Item for ICC across frames", NONE),
        p(class = "text-muted small",
          "Units rho = alpha (set) x phi (group) on a common arbitrary scale. Within a frame all curves are parallel; across frames they fan with the unit. Extended frame of reference analyses only.")),
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(7, 5)),
        tableCard("frame_tbl", "Frames: units, origins, pooled fit"),
        div(tableCard("phi_tbl", "Person group units (phi)"),
            tableCard("alpha_tbl", "Item set units (alpha) and locations"))),
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
        plotCard("frame_plot", "Frame units"),
        plotCard("frame_icc", "ICC across frames")),
      card(card_header("Equal-unit comparison"),
           card_body(verbatimTextOutput("efrm_cmp")))
    )
  ),

  # ------------------------------------------------------- DIMENSIONALITY --
  nav_panel("Dimensionality",
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        h6("t-test item subsets"),
        selectizeInput("dim_pos", "Subset A", NULL, multiple = TRUE,
                       options = list(placeholder = "default: positive PC1 loadings")),
        selectizeInput("dim_neg", "Subset B", NULL, multiple = TRUE,
                       options = list(placeholder = "default: negative PC1 loadings")),
        actionButton("dim_apply", "Run t-test with these subsets",
                     class = "btn-outline-primary w-100"),
        p(class = "text-muted small mt-2",
          "Leave both empty (and press the button) to return to the first-contrast split. Persons extreme on either subset are excluded; the proportion of significant tests carries an exact binomial confidence interval.")),
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(5, 7)),
        card(card_header("Unidimensionality t-test (Smith)"),
             card_body(verbatimTextOutput("dim_txt"))),
        plotCard("scree", "Scree of the residual components")),
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(5, 7)),
        plotCard("pca_plot", "Residual first contrast"),
        tableCard("loadings_tbl", "Component loadings (first 10)")),
      tableCard("eigen_tbl", "Residual eigenvalues (first 10)")
    )
  ),

  # ------------------------------------------------------ LOCAL DEPENDENCE --
  nav_panel("Local dependence",
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(7, 5)),
      plotCard("rcor", "Residual correlations", height = "520px"),
      div(
        tableCard("rpairs_tbl", "Flagged dependent pairs",
                  "Pairs more than 0.2 above the average off-diagonal residual correlation."),
        card(card_header("Subtest (combine dependent items)"),
          card_body(
            p(class = "text-muted small",
              "Select two or more items to merge into one polytomous super-item and re-analyse; the dependence is absorbed into the subtest."),
            selectizeInput("subtest_items", NULL, NULL, multiple = TRUE,
                           options = list(placeholder = "items to combine")),
            actionButton("make_subtest", "Combine and re-analyse",
                         class = "btn-outline-primary w-100"),
            uiOutput("subtest_status")))))
  ),

  # -------------------------------------------------------------- COMPARE --
  nav_panel("Compare",
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        actionButton("keep_fit", "Keep current fit for comparison",
                     class = "btn-primary w-100"),
        actionButton("clear_fits", "Clear kept fits",
                     class = "btn-outline-secondary w-100 mt-2"),
        p(class = "text-muted small mt-3",
          "Run an analysis, keep it, change the model or settings, run again, and keep that too. For fits of the same data the pairwise conditional log-likelihoods are compared directly (descriptive, composite likelihood; most meaningful for nested structures such as RSM inside PCM). Across different data preparations, compare the calibration-free columns: chi-square per df, fit residual SDs (ideal 1), PSI, and alpha.")),
      tableCard("cmp_tbl", "Model comparison",
                "Reference for two_delta_ll is the first kept fit.")
    )
  ),

  # --------------------------------------------------------------- EXPORT --
  nav_panel("Export",
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
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
          tags$li("Item statistics, thresholds with SEs, person estimates (with ID and factors), score-to-measure table"),
          tags$li("Residual correlations, flagged dependent pairs, PCA loadings, category frequencies, DIF ANOVA for every factor"),
          tags$li("Person-item distribution, threshold map, TCC, TIF, item and person fit maps, residual heatmap, PCA plot"),
          tags$li("Per-item ICC, category curves, threshold curves, and frequency charts"),
          tags$li("For many-facet analyses: facet severities with SEs and fit, structural item thresholds, and severity caterpillar plots"),
          tags$li("summary.txt with the full test-of-fit report"))))
    )
  )
)

server <- function(input, output, session) {

  # ------------------------------------------------------------- data in --
  # picking an example dataset also selects the matching model; uploading a
  # file clears the example selection
  observeEvent(input$demo_choice, {
    if (!identical(input$demo_choice, "none"))
      updateRadioButtons(session, "model_type", selected = input$demo_choice)
  }, ignoreInit = TRUE)
  observeEvent(input$file,
    updateSelectInput(session, "demo_choice", selected = "none"))

  raw_data <- reactive({
    if (!identical(input$demo_choice %||% "none", "none"))
      return(switch(input$demo_choice,
                    dich = .demo_dich(), rsm = .demo_rsm(),
                    mfrm = .demo_long(), efrm = .demo_efrm(), .demo_data()))
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
    updateSelectInput(session, "id_col", choices = c(NONE, nm),
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
    updateSelectInput(session, "lp_person", choices = c(NONE, nm),
                      selected = if (!is.na(g_per)) g_per else NONE)
    updateSelectInput(session, "lp_item", choices = c(NONE, nm),
                      selected = if (!is.na(g_itm)) g_itm else NONE)
    updateSelectInput(session, "lp_score", choices = c(NONE, nm),
                      selected = if (!is.na(g_sco)) g_sco else NONE)
    updateSelectizeInput(session, "lp_facets", choices = nm, selected = g_fac)
    # frames layout guesses
    g_grp <- nm[grepl("group|year|grade|cohort|class$", tolower(nm))][1]
    updateSelectInput(session, "ef_id", choices = c(NONE, nm),
                      selected = if (!is.na(guess_id)) guess_id else NONE)
    updateSelectInput(session, "ef_group", choices = c(NONE, nm),
                      selected = if (!is.na(g_grp)) g_grp else NONE)
    updateSelectizeInput(session, "ef_items", choices = nm,
                         selected = setdiff(nm, c(guess_id, g_grp)))
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
                      choices = c(NONE, input$lp_facets), selected = sel)
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

  output$data_info <- renderUI({
    if (identical(input$demo_choice %||% "none", "none") && is.null(input$file))
      return(p(class = "text-muted",
               "Upload a CSV/TSV file, or pick an example dataset in the sidebar, to begin."))
    df <- raw_data()
    p(class = "text-muted",
      sprintf("%d rows x %d columns. Nominate the column roles in the sidebar, then run the analysis. Missing responses may be left blank or coded as -1; any negative score is read as missing.",
              nrow(df), ncol(df)))
  })
  output$preview <- renderDT({
    datatable(head(raw_data(), 200), rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })

  # ----------------------------------------------------------------- fit --
  override_fit <- reactiveVal(NULL)
  # clear any subtest/split override as soon as a fresh run is requested;
  # fit() short-circuits on the override, so analysis() cannot clear it itself
  observeEvent(input$run, override_fit(NULL), priority = 10)

  analysis <- eventReactive(input$run, {
    df <- raw_data()
    adjN <- NA   # the chi-square sample-size adjustment is a display option
                 # on the Items page, applied without refitting
    withProgress(message = "Estimating (pairwise conditional ML)…", value = 0.3, {
      fit <- tryCatch({
        if (identical(input$model_type, "efrm")) {
          if (is.null(input$ef_group) || input$ef_group == NONE)
            stop("nominate the person group column")
          rasch_efrm(df,
                     item_sets = ef_setmap(),
                     groups = input$ef_group,
                     id = if (!is.null(input$ef_id) && input$ef_id != NONE)
                       input$ef_id else NULL,
                     items = names(ef_setmap()),
                     n_groups = input$ng, adjust_N = adjN,
                     maxit = max(5, input$maxit %||% 60),
                     tol = max(1e-12, input$tol %||% 1e-8))
        } else if (identical(input$model_type, "mfrm")) {
          if (any(c(input$lp_person, input$lp_item, input$lp_score) == NONE) ||
              !length(input$lp_facets))
            stop("nominate the person, item, score, and at least one facet column")
          rasch_mfrm(df, person = input$lp_person, item = input$lp_item,
                     score = input$lp_score, facets = input$lp_facets,
                     n_groups = input$ng, adjust_N = adjN,
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
            if (!is.null(kf) && all(c("item", "key") %in% names(kf)))
              mc_key <- kf
            else showNotification("Key CSV needs columns item,key - ignored.",
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
          f0 <- rasch(df, model = if (identical(input$model_type, "rsm")) "RSM" else "PCM",
                      id = idc, factors = fac, items = its,
                      n_groups = input$ng, adjust_N = adjN, anchors = anc,
                      key = mc_key,
                      maxit = max(5, input$maxit %||% 60),
                      tol = max(1e-12, input$tol %||% 1e-8))
          if (identical(input$model_type, "dich") && any(f0$m > 1L))
            showNotification("Some items have more than two categories; they were fitted with partial credit thresholds.",
                             type = "warning", duration = 10)
          f0
        }
      }, error = function(e) e)
    })
    if (inherits(fit, "error")) {
      showNotification(paste("Analysis failed:", conditionMessage(fit)),
                       type = "error", duration = NULL)
      return(NULL)
    }
    if (length(fit$notes))
      showNotification(paste(fit$notes, collapse = "\n"), type = "warning", duration = 12)
    override_fit(NULL)
    try(nav_select("Summary", session = session), silent = TRUE)
    fit
  })
  fit <- reactive({
    f <- override_fit()
    if (is.null(f)) f <- analysis()
    req(f); f
  })

  observeEvent(fit(), {
    its <- fit()$items$item
    updateSelectInput(session, "dif_item", choices = its, selected = its[1])
    updateSelectizeInput(session, "subtest_items", choices = its, selected = character(0))
    fac <- names(fit()$factors)
    dif_choices <- if (length(fac) > 1) c(fac, FACTORIAL) else
      if (length(fac)) fac else NONE
    updateSelectInput(session, "dif_factor", choices = dif_choices,
                      selected = dif_choices[1])
    fs <- if (inherits(fit(), "rasch_mfrm")) fit()$facet_spec else NONE
    updateSelectInput(session, "facet_sel", choices = fs, selected = fs[1])
    fi <- if (inherits(fit(), "rasch_efrm"))
      unique(fit()$virtual_map$item) else NONE
    updateSelectInput(session, "frame_item", choices = fi, selected = fi[1])
    updateSelectizeInput(session, "dim_pos", choices = its, selected = character(0))
    updateSelectizeInput(session, "dim_neg", choices = its, selected = character(0))
  })

  observeEvent(input$make_subtest, {
    req(length(input$subtest_items) >= 2)
    res <- tryCatch(combine_items(fit(), list(input$subtest_items)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Subtest failed:", conditionMessage(res)), type = "error")
    } else {
      override_fit(res)
      showNotification("Re-analysed with the subtest in place. Run analysis again to reset.",
                       type = "message", duration = 8)
    }
  })
  output$subtest_status <- renderUI({
    if (is.null(override_fit())) return(NULL)
    p(class = "text-success small mt-2",
      paste("Active:", paste(grep("\\+", fit()$items$item, value = TRUE), collapse = "; ")))
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
      showNotification(
        sprintf("Re-analysed with %s split by %s. Run analysis again to reset.",
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
  register_plot <- function(id, fun, w = 9, h = 6) {
    output[[id]] <- renderPlot(fun(), res = 96)
    for (fmt in c("png", "pdf")) local({
      fmt_ <- fmt
      output[[paste0(id, "_", fmt_)]] <- downloadHandler(
        filename = function() paste0("RaschR_", id, ".", fmt_),
        content = function(file) {
          # 300 dpi PNG (and vector PDF) for publication
          if (fmt_ == "png") png(file, width = w, height = h, units = "in", res = 300)
          else pdf(file, width = w, height = h)
          fun(); dev.off()
        })
    })
  }
  register_table <- function(id, fun, dt_fun) {
    output[[id]] <- renderDT(dt_fun())
    output[[paste0(id, "_csv")]] <- downloadHandler(
      filename = function() paste0("RaschR_", id, ".csv"),
      content = function(file) write.csv(fun(), file, row.names = FALSE))
  }
  num_dt <- function(d, digits = 3, ...) {
    num <- vapply(d, is.numeric, TRUE)
    datatable(d, rownames = FALSE, ...,
              options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
      formatRound(names(d)[num], digits)
  }

  # -------------------------------------------------------------- summary --
  output$vboxes <- renderUI({
    f <- fit()
    layout_column_wrap(width = "165px", fill = FALSE, class = "mb-3",
      value_box("Persons", nrow(f$X), theme = "primary"),
      value_box("Items", ncol(f$X), theme = "primary"),
      value_box("PSI", sprintf("%.3f", f$psi$PSI), theme = "success",
                p(class = "small mb-0", sprintf("%.3f no extremes", f$psi_noext$PSI))),
      value_box("Alpha", sprintf("%.3f", f$alpha$alpha), theme = "success",
                p(class = "small mb-0", sprintf("n = %d complete", f$alpha$n))),
      value_box("Item-trait p", sprintf("%.3f", f$total_chisq_p),
                theme = if (f$total_chisq_p < 0.05) "danger" else "secondary"),
      value_box("Power of fit", f$power_of_fit, theme = "secondary")
    )
  })

  output$fit_summary <- renderPrint({
    f <- fit()
    cat(sprintf("Model: %s  |  Estimation: pairwise conditional ML (%s, %d iterations)\n",
                f$model, if (f$est$converged) "converged" else "NOT CONVERGED",
                f$est$iterations))
    cat(sprintf("Total item-trait chi-square: %.3f on %d df, p = %.3f\n",
                f$total_chisq, f$total_df, f$total_chisq_p))
    cat(sprintf("Item fit residual:   mean %6.2f  SD %5.2f  (ideal 0, 1)\n",
                f$item_fit_summary$mean, f$item_fit_summary$sd))
    cat(sprintf("Person fit residual: mean %6.2f  SD %5.2f  (ideal 0, 1)\n",
                f$person_fit_summary$mean, f$person_fit_summary$sd))
    cat(sprintf("Items flagged misfitting (BH-adjusted): %d of %d\n",
                sum(f$items$misfit, na.rm = TRUE), nrow(f$items)))
    dis <- names(which(vapply(f$thresholds_diag, function(d)
      !d$ordered && length(d$thresholds) > 1, TRUE)))
    cat("Disordered thresholds:", if (length(dis)) paste(dis, collapse = ", ") else "none", "\n")
    if (length(f$notes)) cat("Notes:", paste(f$notes, collapse = "; "), "\n")
  })

  output$targeting <- renderPrint({
    f <- fit(); t <- f$targeting
    cat(sprintf("Person mean (SD):     %6.3f (%.3f) logits\n", t$person_mean, t$person_sd))
    cat(sprintf("Person mean, no extremes: %.3f\n", t$person_mean_noext))
    cat(sprintf("Item mean:             0.00 (constrained)\n"))
    cat(sprintf("Threshold range:      %6.3f to %.3f\n",
                t$threshold_range[1], t$threshold_range[2]))
    cat(sprintf("Persons beyond thresholds: %.1f%% below, %.1f%% above\n",
                100 * t$prop_below, 100 * t$prop_above))
    cat(sprintf("\nPSI: %.3f (separation %.3f)\n", f$psi$PSI, f$psi$separation))
    cat(sprintf("PSI without extremes: %.3f (n = %d)\n", f$psi_noext$PSI, f$psi_noext$n))
    cat(sprintf("Cronbach alpha: %.3f (n = %d complete cases)\n", f$alpha$alpha, f$alpha$n))
  })

  register_table("score_tbl", function() {
    f <- fit()
    if (!is.null(f$score_table)) f$score_table else f$score_curves
  }, function() {
    f <- fit()
    if (!is.null(f$score_table)) {
      datatable(f$score_table, rownames = FALSE,
                options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        formatRound(c("theta", "se"), 3)
    } else {
      validate(need(!is.null(f$score_curves),
                    "No score conversion available for this fit."))
      datatable(f$score_curves, rownames = FALSE,
                caption = "Raw scores are not sufficient under unequal frame units; per-group expected-score curves replace the score table.",
                options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        formatRound(c("theta", "expected_score", "sem"), 3)
    }
  })
  register_table("thr_tbl", function() {
    f <- fit(); d <- f$thresholds
    data.frame(item = f$items$item[d$item], threshold = d$k,
               tau = d$tau, se = d$se)
  }, function() {
    f <- fit(); d <- f$thresholds
    num_dt(data.frame(item = f$items$item[d$item], threshold = d$k,
                      tau = d$tau, se = d$se))
  })

  # ---------------------------------------------------------------- items --
  # item table with the optional chi-square sample-size adjustment applied
  # on display (the chi-square scales linearly in N, so no refit is needed)
  items_view <- reactive({
    d <- fit()$items
    adjN <- input$adjN
    if (!is.null(adjN) && !is.na(adjN) && adjN > 0) {
      n_used <- sum(!is.na(fit()$person$class_interval))
      d$chisq <- d$chisq * adjN / n_used
      d$p <- pchisq(d$chisq, d$df, lower.tail = FALSE)
      d$p_adj <- p.adjust(d$p, method = "BH")
      d$misfit <- d$p_adj < 0.05
    }
    d
  })
  register_table("items_tbl", function() items_view(), function() {
    d <- items_view()
    d$misfit <- ifelse(d$misfit, "*", "")
    num_dt(d, selection = "single") |>
      formatRound(c("p", "p_adj"), 3)
  })
  register_plot("icc",  function() plot_icc(fit(), sel_item()))
  register_plot("ccc",  function() plot_ccc(fit(), sel_item()))
  register_plot("tpc",  function() plot_threshold_prob(fit(), sel_item()))
  register_plot("cfreq", function() plot_catfreq(fit(), sel_item()))
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
  })
  register_plot("distractor_plot", function() {
    f <- fit()
    validate(need(!is.null(f$mc),
                  "Provide a multiple-choice key (CSV: item,key) to see option curves."))
    it <- if (sel_item() %in% colnames(f$mc$raw)) sel_item() else
      colnames(f$mc$raw)[1]
    plot_distractors(f, it)
  })

  # -------------------------------------------------------------- persons --
  register_table("person_tbl", function() fit()$person, function() {
    d <- fit()$person
    datatable(d, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
      formatRound(names(d)[vapply(d, is.numeric, TRUE) &
                           !names(d) %in% c("raw", "max_raw", "n_items", "class_interval")], 3)
  })
  register_plot("pfit",  function() plot_person_fit(fit()))
  register_plot("pim_p", function() plot_pimap(fit()))

  # ------------------------------------------------------------ test plots --
  register_plot("thrmap", function() plot_threshold_map(fit()), h = 7)
  register_plot("imap",   function() plot_item_map(fit()))
  register_plot("tcc",    function() plot_tcc(fit()))
  register_plot("tif",    function() plot_tif(fit()))
  register_plot("guttman", function() plot_guttman(fit()), h = 7)

  # ------------------------------------------------------------------ DIF --
  FACTORIAL <- "(all factors: factorial)"
  dif_res <- reactive({
    f <- fit(); req(!is.null(f$factors), length(names(f$factors)) > 0)
    dif_anova(f)
  })
  dif_fact <- reactive({
    f <- fit(); req(!is.null(f$factors), length(names(f$factors)) > 1)
    dif_anova_factorial(f)
  })
  register_table("dif_tbl", function() {
    if (identical(input$dif_factor, FACTORIAL)) dif_fact()$terms else dif_res()
  }, function() {
    if (identical(input$dif_factor, FACTORIAL)) {
      d <- dif_fact()$terms
      d$significant <- ifelse(d$significant, "*", "")
      d$superseded <- ifelse(d$superseded, "(superseded)", "")
      num_dt(d) |> formatRound(c("p", "p_adj"), 3)
    } else {
      d <- dif_res()
      if (!is.null(input$dif_factor) && input$dif_factor %in% d$factor)
        d <- d[d$factor == input$dif_factor, ]
      d$uniform_DIF <- ifelse(d$uniform_DIF, "*", "")
      d$nonuniform_DIF <- ifelse(d$nonuniform_DIF, "*", "")
      num_dt(d) |> formatRound(c("p_uniform", "p_nonuniform",
                                "p_uniform_adj", "p_nonuniform_adj"), 3)
    }
  })
  register_table("dif_tukey_tbl", function() dif_fact()$tukey, function() {
    validate(need(identical(input$dif_factor, FACTORIAL),
                  "Choose the factorial option in the sidebar to see Tukey HSD comparisons."))
    tk <- dif_fact()$tukey
    if (!nrow(tk))
      return(datatable(data.frame(note = "no significant group terms to compare"),
                       rownames = FALSE, options = list(dom = "t")))
    num_dt(tk) |> formatRound("p_tukey", 3)
  })
  register_plot("dif_icc", function() {
    f <- fit()
    req(input$dif_item %in% f$items$item,
        !is.null(f$factors), input$dif_factor %in% names(f$factors))
    plot_icc(f, input$dif_item, group = input$dif_factor)
  })

  # ---------------------------------------------------------------- facets --
  facet_dat <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a long-format (many-facet) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    f$facet_effects[[input$facet_sel]]
  })
  register_table("facet_tbl", function() facet_dat(), function() {
    d <- facet_dat()
    datatable(d, rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
      formatRound(setdiff(names(d)[vapply(d, is.numeric, TRUE)], "n"), 3)
  })
  register_plot("facet_plot", function() {
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a long-format (many-facet) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    plot_facets(f, input$facet_sel)
  })
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
  })

  # ----------------------------------------------------------------- equating --
  eq_ref <- reactive({
    req(input$eq_file)
    a <- tryCatch(read.csv(input$eq_file$datapath, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    validate(need(!is.null(a) && all(c("item", "location") %in% names(a)),
                  "The reference CSV needs columns item, location (and ideally se)."))
    a
  })
  eq_res <- reactive({
    validate(need(!is.null(input$eq_file),
                  "Upload a reference calibration (item, location, se) to equate against."))
    equate_tests(fit(), eq_ref(), shift = input$eq_shift)
  })
  register_table("eq_tbl", function() eq_res()$table, function() {
    eq <- eq_res()
    d <- eq$table
    d$drift <- ifelse(d$drift, "*", "")
    num_dt(d) |> formatRound(c("p", "p_adj"), 3)
  })
  register_plot("eq_plot", function() {
    req(input$eq_file)
    plot_equate(fit(), eq_ref(), shift = input$eq_shift)
  })
  output$dl_anchors <- downloadHandler(
    filename = function() format(Sys.time(), "RaschR_anchors_%Y%m%d_%H%M.csv"),
    content = function(file) {
      f <- fit()
      thr <- f$thresholds
      write.csv(data.frame(item = f$items$item[thr$item], k = thr$k,
                           tau = thr$tau), file, row.names = FALSE)
    })

  output$dl_calib <- downloadHandler(
    filename = function() format(Sys.time(), "RaschR_calibration_%Y%m%d_%H%M.csv"),
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
                 function() num_dt(efrm_fit()$frames))
  register_table("phi_tbl", function() efrm_fit()$phi_table,
                 function() num_dt(efrm_fit()$phi_table))
  register_table("alpha_tbl", function()
    merge(efrm_fit()$alpha_table, efrm_fit()$set_table[, c("set", "mu", "n_items")],
          by = "set"),
    function() num_dt(merge(efrm_fit()$alpha_table,
                            efrm_fit()$set_table[, c("set", "mu", "n_items")],
                            by = "set")))
  register_plot("frame_plot", function() plot_frames(efrm_fit()))
  register_plot("frame_icc", function() {
    f <- efrm_fit()
    req(input$frame_item %in% f$virtual_map$item)
    plot_icc_frames(f, input$frame_item)
  })
  output$efrm_cmp <- renderPrint({
    f <- efrm_fit(); cmp <- f$efrm_vs_rasch
    cat(sprintf("Pairwise conditional log-likelihood: frames model %.3f, equal units %.3f\n",
                cmp$ll_efrm, cmp$ll_equal))
    cat(sprintf("2 x improvement: %.3f with %d extra unit parameter(s)\n",
                cmp$two_delta_ll, cmp$extra_parameters))
    cat("(composite likelihood: descriptive, not a calibrated chi-square)\n")
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
    cat("\nSubset A items:\n ", paste(dt$items_positive, collapse = ", "), "\n")
    cat("Subset B items:\n ", paste(dt$items_negative, collapse = ", "), "\n")
  })
  register_plot("scree", function() plot_scree(fit()))
  register_table("loadings_tbl", function() residual_pca(fit())$loadings_matrix,
                 function() num_dt(residual_pca(fit())$loadings_matrix))
  register_table("eigen_tbl", function() residual_pca(fit())$eigen_table,
                 function() num_dt(residual_pca(fit())$eigen_table))
  register_plot("pca_plot", function() plot_pca(fit()))

  # -------------------------------------------------------- local dependence --
  register_plot("rcor", function() plot_resid_cor(fit()), w = 8, h = 8)
  register_table("rpairs_tbl", function() residual_correlations(fit())$flagged,
                 function() {
    fl <- residual_correlations(fit())$flagged
    if (!nrow(fl)) fl <- data.frame(note = "no item pairs exceed the flag threshold")
    num_dt(fl)
  })

  # ---------------------------------------------------------------- compare --
  kept_fits <- reactiveVal(list())
  observeEvent(input$keep_fit, {
    f <- fit()
    k <- kept_fits()
    lab <- sprintf("%d_%s", length(k) + 1L, f$model)
    k[[lab]] <- f
    kept_fits(k)
    showNotification(sprintf("Kept '%s' (%d fit(s) held).", lab, length(k)),
                     type = "message", duration = 5)
  })
  observeEvent(input$clear_fits, {
    kept_fits(list())
    showNotification("Cleared kept fits.", type = "message", duration = 4)
  })
  cmp_res <- reactive({
    k <- kept_fits()
    validate(need(length(k) >= 2,
                  "Keep at least two fits (run, keep, change settings, run, keep) to compare."))
    as.data.frame(do.call(compare_fits, k))
  })
  register_table("cmp_tbl", function() cmp_res(), function() {
    d <- cmp_res()
    d$same_data <- ifelse(d$same_data, "yes", "no")
    num_dt(d)
  })

  # ----------------------------------------------------------------- export --
  output$dl_zip <- downloadHandler(
    filename = function() format(Sys.time(), "RaschR_results_%Y%m%d_%H%M.zip"),
    content = function(file) {
      f <- fit()
      tmp <- file.path(tempdir(), paste0("raschr_", as.integer(Sys.time())))
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
