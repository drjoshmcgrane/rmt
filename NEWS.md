# rasch 1.8.0

Three further paired-comparison analogues of the Rasch tool set, and CRAN
housekeeping.

* Anchored estimation for paired comparisons: `btl(anchors = c(Object =
  value, ...))` holds the named objects at fixed locations and estimates
  the rest on that scale (no sum-zero constraint; anchored objects report
  se 0) -- the equating workflow's second half, mirroring
  `rasch(anchors=)`. An anchored object at a response boundary is an
  error rather than a silent removal.
* First-position advantage: `btl(position = TRUE)` estimates the pure
  positional (first/left-presented) advantage of the Davidson & Beaver
  (1977) order-effect device as a constant covariate alongside any
  exposure and carry-over effects; it is identified by triangle closure
  even under fixed presentation orders (more weakly, with honest SEs).
  The dimensionality reference simulation carries the fitted coefficient.
* Common-object equating: `btl_equate()` compares two paired-comparison
  calibrations (or a calibration against a bank) through their common
  objects -- precision-weighted origin shift, per-object drift t-tests
  with multiplicity adjustment, equated locations -- with
  `plot_btl_equate()`; the standards-maintenance workflow of comparative
  judgement across panels or years.
* Targeting and adaptive pairing: `btl_information()` (per-comparison
  Fisher information and per-object design information, with the naive
  1/sqrt(information) SE beside the judge-clustered sandwich),
  `plot_btl_targeting()`, and `btl_next_pairs()` -- the greedy
  most-informative-next-pair step of adaptive comparative judgement
  (Pollitt 2012), with the Bramley (2015) reliability-inflation caution
  documented.
* CRAN housekeeping: pure-ASCII sources, quoted 'shiny' in DESCRIPTION,
  heavy Monte-Carlo recovery tests wrapped in `skip_on_cran()` (they run
  locally and on CI), and `cran-comments.md` drafted from a full
  `--as-cran` check (about 3 minutes; slowest example 1.4s).
* An adversarial verification round over the new features also caught and
  fixed a PRE-EXISTING inference bug: `btl(count=)` inflated every sandwich
  standard error by about the square root of the count (the meat weighted
  aggregated rows twice; point estimates were always correct). Count-
  weighted fits now reproduce the expanded-data standard errors exactly;
  OSI and all z/p values follow. Additionally from that round:
  `btl_equate()`'s shift standard error now uses the stored covariance of
  each calibration (it treated the common-object differences as
  independent, overstating the SE several-fold; 95 per cent coverage
  verified) and documents the majority-drift limitation of single-shift
  equating; `btl_next_pairs()` ranks by the exact one-step reduction in
  total location variance (a Sherman-Morrison update on the fit's
  covariance), verified to beat low-priority pairs empirically; and
  `btl_information()`'s `se_naive` is documented honestly as a
  single-parameter lower bound rather than an independence benchmark.
* Shiny app: the Targeting and Equating tabs gain paired-comparison
  variants (design information + targeting plot + adaptive next-pair
  recommendations; reference-calibration upload with drift table and
  plot), and the Data roles gain a first-position-advantage switch and an
  object-anchors upload.

# rasch 1.7.1

A statistical audit of the new simulation and paired-comparison
independence code (23 findings across two adversarial rounds, all fixed),
plus a vignette and a logo.

* Misfit layers now compose correctly in `simulate_rasch()`: the
  response-dependence and response-style regeneration steps previously
  rebuilt their target items from the base structure, silently erasing any
  DIF or second-dimension loading planted on those items (and leaking a
  DIF-like signal from a DIF-carrying source item into its dependent
  partner). Every regeneration now honours the item's own generating
  structure -- trait, group-shifted thresholds and slope, guessing.
* `btl_dimensionality()` now simulates its noise reference sequentially
  WITH the fitted within-judge dependence effects when the fit carries
  them: real order effects previously inflated the false-positive rate to
  ~36 per cent on one-dimensional data (now ~0-5 per cent). The documented
  price is conservatism, since carry-over and a judge-camp second
  attribute are partially confounded.
* Doc examples pass `id = "id"` so they no longer emit spurious
  dropped-column notes.

* `simulate_rasch()`: `model = "PCM"` now draws each item's threshold
  spacings and span afresh, as the partial credit model allows (previously
  PCM and RSM generated identical rating-scale-structured data);
  `dif` without `n_groups >= 2` is an error instead of a silent no-op;
  polytomous `guessing` warns and is ignored rather than falsely recorded
  as planted; the response-dependence term now carries the guessing
  asymptote of the pair's first item (no more spurious easiness on the
  dependent item); disordered thresholds work with 3 categories and no
  longer shift the item's location; careless responders are no longer
  double-counted as response-style persons in the recorded truth.
* `simulate_mfrm()`: halo raters are capped at the eligible pool (no NA
  raters when erratic raters shrink it); the example runs.
* `sim_recovery()`: person abilities are mean-centred as documented (an
  asymmetric difficulty range no longer masquerades as person bias); the
  many-facet branch reads the item margins (`item_effects`), restoring the
  documented item-difficulty recovery for MFRM fits.
* `btl_dimensionality()`: null calibration verified for graded fits
  (mildly conservative, never anticonservative) and documented; docs also
  clarify the leading-bimension reference columns and the no-ties
  precondition of Kendall's zeta.
* Shiny app (Simulate page): the parameter-recovery card only renders when
  the current fit was estimated on the currently loaded simulation (a
  stale fit can no longer be compared against new truth); the reproducible
  call now carries every argument the generator received (MFRM person/item
  SDs, EFRM person SD, graded categories); non-integer seeds are rounded
  safely.
* New vignette: "Planting misfit and watching the diagnostics fire" — the
  plant-to-detect loop end to end, closing with a small Monte Carlo power
  estimate via `sim_replicate()`.
* A hex logo (an item characteristic curve at its location), reproducible
  from `tools/logo.R`.

# rasch 1.7.0

* Simulation gains the standard controls of a proper simulation tool.
  - Population parameters: the person distribution (mean, SD, shape) and the
    item-difficulty range are now set directly (in the package and the app).
  - More misfit: response styles (extreme / midpoint categories) and
    speededness (a not-reached tail) for \code{simulate_rasch()}, and a halo
    effect for \code{simulate_mfrm()}.
  - \code{sim_replicate()} generates a batch of datasets for Monte Carlo
    power, Type-I, or recovery studies.
  - \code{sim_recovery()} and \code{plot_recovery()} compare the parameters a
    fit recovers with the ones that were planted (correlation, RMSE, bias, and
    a true-vs-estimated scatter) for every layout.
* Shiny app: the Simulate page renames the wide layout to "Rasch", exposes
  the population parameters, shows a parameter-recovery report and the
  true-vs-estimated scatter after Run, and prints the reproducible
  \code{simulate_*()} call.

# rasch 1.6.1

* Shiny app: a Simulate page (under More) exposes all four simulators with a
  data-type selector and per-layout misfit controls. Pressing Simulate
  generates the data, loads it as the current dataset with its roles set, and
  lists what was planted; pressing Run then analyses it so the matching
  diagnostic can be watched firing. Verified end to end for every layout
  (e.g. a planted DIF item is flagged by dif_anova, a many-facet fit recovers
  its rater severities).

# rasch 1.6.0

* Data simulation across all four model families, each with dial-in
  departures from the model -- so a known pathology can be planted and the
  matching diagnostic watched as it fires. Every result carries the true
  parameters (\code{attr(x, "truth")}) and prints what was planted.
  - \code{simulate_rasch()}: dichotomous/PCM/RSM with discrimination misfit,
    guessing, a second dimension, local dependence, DIF (uniform and
    non-uniform), careless responders, disordered thresholds, and missing
    data.
  - \code{simulate_btl()}: paired comparisons with erratic judges, a second
    object attribute (two judge camps), and within-judge exposure/carry-over
    dependence; dichotomous or graded.
  - \code{simulate_mfrm()}: rated data with rater-severity spread, erratic
    raters, and a rater-by-item interaction.
  - \code{simulate_efrm()}: frames with differing set and group units.
  Verified end to end: each planted departure is recovered or flagged by the
  corresponding fit / diagnostic (an app panel follows).

# rasch 1.5.1

* The judge unexpected-judgements map (\code{plot_btl_judge_map}) is now
  drawn at the matchup level: each pair the judge met is a segment spanning
  its two objects on the consensus location axis, placed by how surprising
  the verdict was, with the judge's backed object filled. \code{judge_pair_surprise()}
  returns the per-matchup residuals behind it. The object-level
  \code{judge_surprise()} is retained for per-object residuals.

# rasch 1.5.0

* `judge_surprise()` and `plot_btl_judge_map()`: the paired-comparison
  counterpart of the kidmap. A judge has no ability to condition on, so the
  reference is the consensus object scale; each object a judge met gets a
  standardised residual saying how much more or less that judge favoured it
  than its location predicts. A strong object the judge under-rated, or a
  weak object over-rated, is an unexpected judgement. It catches systematic
  bias -- a distinct failure mode from the noise that judge fit and judge
  consistency flag. In the Shiny app the Persons > Judge fit panel is now a
  master-detail: the judges table drives the selected judge's map, mirroring
  person -> kidmap.

# rasch 1.4.0

* Paired comparisons gain the pair-structure analogues of the Rasch
  independence diagnostics -- the tools that a persons-by-items residual
  matrix would give but paired-comparison data cannot.
  - `btl_transitivity()` counts circular triads (preference loops: A beats B,
    B beats C, C beats A) against the chance rate, reports Kendall's
    coefficient of consistency for complete designs, and gives each judge's
    own consistency (Kendall & Babington Smith 1940) -- the single-dimension
    check, and a judge-fit analogue.
  - `btl_dimensionality()` decomposes the skew-symmetric object-by-object
    residual preference matrix into Gower (1977) bimensions -- the
    paired-comparison counterpart of residual PCA -- and judges the leading
    bimension against a model-simulated noise reference, so a coherent
    residual "swirl" (a second attribute steering some contests) is
    distinguished from noise.
  - `plot_btl_scree()`, `plot_btl_dim_map()`, and `plot_btl_transitivity()`
    display them; the Shiny app's Independence > Trait tab now carries these
    for a BTL fit (it previously hid for paired comparisons).

# rasch 1.3.1

* The package is now called **rasch** (it was developed under the working
  title `rmt`). The name `rmt` belongs to an unrelated package on CRAN, and
  `rasch` says plainly what this is -- and mirrors `mirt::mirt()` with
  `rasch::rasch()`. No functions, arguments, or results change; only the
  package and library name. The auxiliary result classes were unified under
  the `rasch_` prefix (`rasch_btl`, `rasch_dif`, ...) to match the fit
  classes.

# rasch 1.3.0

A full statistical audit of the independence, invariance, and paired-
comparison procedures, with fixes throughout.

* Paired comparisons: `btl()` now interrogates within-judge dependence --
  `dependence_data` holds every comparison's exposure and carry-over
  covariates (count-weighted), `plot_btl_dependence()` displays an effect
  as a partial-residual curve, and a separated effect (one-sided
  informative comparisons) is set aside with a note instead of reported as
  a runaway estimate. The graded engine's threshold-by-dependence Hessian
  block is corrected (free-threshold fits with five or more categories and
  an order column now converge), rows with a missing judge are dropped with
  a note, category checks re-run after boundary objects are removed, and
  the dichotomous path now routes through the one (vectorised) estimator.
* `btl_dif()` holds fitted dependence effects fixed in its residual moments
  and resolution refits (dependence is no longer absorbed as spurious
  judge-group DIF), weights count-aggregated comparisons exactly as
  expanded rows, restricts supersession to group terms, and notes every
  significant term it does not resolve.
* `dif_anova()` keeps each ANOVA term in its own error stratum, so a mixed
  (within-subject) design with missing responses no longer un-flags real
  DIF; factors named like internals (`ci`, `f1`) are safe; a term's own
  class-interval crossing no longer supersedes it (so `resolve_dif()` again
  splits items with both uniform and non-uniform DIF); mixed designs say
  plainly that Tukey comparisons are unavailable; and `dif_size()` notes
  when a within-person factor makes its standard errors conservative.
* `plot_scree()`'s parallel-analysis reference is now model-simulated
  (responses drawn from the calibrated model, persons re-estimated), so
  model-true data sits at the reference instead of always "showing
  structure" (Raiche 2005; Chou & Wang 2010).
* Shiny app: judge-group DIF results freeze their run's factors and alpha
  (later sidebar edits cannot silently regroup the overlay or break the
  reproducible snippets), superseded terms are shown, the within-judge
  dependence panel explains dropped effects, and fit flags are consistent
  on every table (fit residual |2.5|; outfit 0.7-1.3; infit 0.8-1.2).

# rasch 1.2.0

* `plot_pca_biplot()` draws the item loadings on the first two residual
  principal components -- the pair that usually carries any interpretable
  contrast -- on equal axes, coloured by the sign of the PC1 loading.
* `residual_correlations()` now also returns the adjusted-Q3 `star_matrix`
  (each Q3 less the average off-diagonal), and `plot_resid_cor()` gains a
  `stat` argument to colour either the raw Q3 or the adjusted Q3*, drawing
  the lower triangle only.
* Shiny app: the trait- and local-dependence pages pair each table with its
  plot (loadings with a PC1-PC2 biplot; Q3 and Q3* matrices with heatmaps),
  the raw and adjusted correlations are each flagged by their own rule
  (|Q3| and Q3*), and every data-restructuring action (subtest, item split,
  automatic DIF resolution) offers an in-place reset to the original data.


* `dif_anova()` is now the single DIF analysis-of-variance function. One
  factor is analysed one-way; several factors are modelled jointly -- the
  statistically correct treatment -- with main effects by default and
  factor-by-factor interactions optional (`effects = "factorial"`). It
  handles within-subject factors (repeated-measures / mixed ANOVA) and
  returns a classed object with a `summary`, the full `terms` table, and
  Tukey comparisons. The separate `dif_anova_factorial()` is removed.
* `resolve_dif()` resolves DIF iteratively by item splitting, largest
  effect first, to clear artificial DIF.

# rasch 1.0.0

First stable release. The package delivers a complete Rasch Measurement
Theory workflow built entirely from published measurement theory, with a
pairwise conditional estimation core and a Shiny interface that exposes
every analysis with reproducing R code attached to every output.

## Models

* `rasch()`: dichotomous and polytomous (partial credit and rating scale)
  Rasch models by pairwise conditional maximum likelihood (Andrich & Luo
  2003; Zwinderman 1995), with Godambe sandwich standard errors, sum-zero
  identification, Warm (1989) weighted likelihood person estimation, and
  the principal-component threshold parameterisation as an estimation
  option.
* `rasch_mfrm()`: many-facet models with additive facet severities or
  item-by-facet interactions, accepting wide (one column per item) or long
  data.
* `rasch_efrm()`: the extended frame-of-reference model (Humphry) with
  frame units estimated by within-frame pairwise conditional likelihood.
* `btl()`: paired comparisons as the conditional form of the dichotomous
  Rasch model (Bradley & Terry 1952; Andrich 1978), including the
  adjacent-categories graded extension (Tutz 1986; Agresti 1992) with
  symmetric thresholds -- the Davidson (1970) ties model is its
  three-category case -- winner-plus-margin data entry, judge-clustered
  errors, and within-judge exposure and carry-over effects estimated
  jointly with the locations.

## Diagnostics

* The test-of-fit suite follows Andrich & Marais (2019): log-of-mean-square
  fit residuals with apportioned degrees of freedom, the item-trait
  chi-square on automatically sized, tie-preserving class intervals with
  per-interval detail, the ANOVA item fit, and person fit.
* Unidimensionality: residual principal components with parallel analysis,
  the Smith (2002) subset t-test, and magnitude estimation.
* Local dependence: Yen's (1984) Q3 with the Christensen, Makransky &
  Horton (2017) flagging convention, response-dependence magnitude
  (Andrich & Kreiner), subtest formation, and the spread test.
* Guessing: tailored analysis with origin-equated comparison calibrations.
* DIF: two-way residual ANOVA per factor, the joint factorial analysis
  with a compact uniform/non-uniform summary, resolved DIF magnitudes in
  logits with familywise control and a practical-significance criterion,
  planned contrasts derived from the factor structure (with
  repeated-measures support via person-level residual scores), and class
  intervals sized to the cells each analysis actually uses. For paired
  comparisons, `btl_dif()` tests object-by-judge-group interaction by the
  same two routes.
* Equating and invariance: common-item comparison with drift flags, and
  calibration anchoring.

## Display and reporting

* A complete base-graphics plot suite: expected value curves with observed
  class-interval means, category and threshold probability curves, the
  person-item threshold distribution, Wright maps, kidmaps (Wright, Mead &
  Ludlow 1980), object characteristic curves for paired comparisons,
  threshold and item maps, test characteristic and information curves,
  residual displays, and Guttman scalograms -- all with adjustable class
  intervals and scale ranges, and batch export to multi-page PDF or ZIP at
  publication resolution.
* `fit_summary_table()` and `targeting_table()` return the headline
  statistics as tidy tables; `save_outputs()` writes every table and plot
  to disk; `report_html()` produces a single-file HTML report.
* `run_app()` launches the Shiny interface: model-adaptive navigation,
  reproducing R code beneath every output, and one-click export.
