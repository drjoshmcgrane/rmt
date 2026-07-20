# Changelog

## rasch 1.11.2

Edge cases from the fourth review round, verified and closed:

- The BTL dimensionality reference now simulates at the PAIR level:
  count weights are summed per unordered pair and rounded, so half-tie
  rows (weight 0.5 in each direction) recombine into whole comparisons.
  Row-level simulation broke on fractional weights (`as.integer(0.5)` is
  a zero binomial size – every reference draw degenerated to 0). The
  per-replicate leading strengths are exposed as `$reference$draws`.
- `btl_efrm(se_method = "judge_bootstrap")` refuses a panel with a
  single judge (resampling one judge returns the same data every time,
  so its SEs would be a spurious zero) and notes panels with fewer than
  five judges as rough.
- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  notes when there are no more judge clusters than parameters: the
  clustered covariance is then rank-deficient by construction and no
  scalar (CR1) correction repairs it – the SEs likely understate.

## rasch 1.11.1

The follow-up review’s residual findings, all verified and closed:

- The BTL dimensionality reference now simulates binomial (multinomial)
  sums for count-weighted rows instead of one weighted Bernoulli, whose
  w^2 p(1-p) variance overdispersed the reference and made the test
  conservative on aggregated data; aggregated and expanded data now give
  the same reference (verified to \< 1%).
- The pairwise chi-square degrees of freedom subtract every estimated
  parameter (position and dependence covariates were omitted), and a
  design with nothing testable reports NA instead of df = 1.
- Judge-clustered covariance carries the standard CR1 small-sample
  factor G/(G-1), on top of the single-judge error and few-cluster note.
- `btl_efrm(se_method = "judge_bootstrap")`: judges resampled with
  replacement within panels and the whole pipeline refitted, so the
  errors carry within-judge dependence as well as stage-one uncertainty
  – closing the documented gap in the parametric bootstrap. On
  model-true data it agrees with the parametric bootstrap; estimates are
  unchanged by construction.
- The next-pair priority documentation is internally consistent (ranking
  heuristic; the Sherman-Morrison update is exact for information
  matrices, a scoring device on the sandwich).

## rasch 1.11.0

Statistical validity, from a second external review. Every finding was
re-derived or reproduced by simulation before fixing; the fixes are
verified by new calibration tests (test-statistical-validity.R).

- Recentred covariance for mixed max scores.
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  recentres thresholds so the mean item location is zero, but the
  covariance stayed under the sum-of-thresholds constraint – identical
  only when every item has the same maximum score. With mixed
  1/2/3-threshold items, per-item empSD/SE ran 0.90-1.21 by item type;
  the linear transform `(I - 1a') cov (I - a1')` restores 0.96-1.09.
  Affects threshold and item SEs, separation, and equating weights on
  mixed designs.
- Warm WLE discrimination cancellation. The common-discrimination
  weighted score is `disc[(R-E) + mu3/(2V)]` – the discrimination
  cancels; a stray `disc` factor in the correction biased WLEs by up to
  0.26 logits at disc = 0.5 (exact at disc = 1, so ordinary analyses
  were unaffected; the vector-disc EFRM path was already correct).
- Honest chi-square degrees of freedom. An item whose responders occupy
  fewer than two class intervals has no estimable item-trait
  interaction: chi-square, df, and p are now NA (df was manufactured as
  1 with p = 1), and the total test sums over testable items only. The
  `adjust_N` proportional-scaling semantics are now documented
  explicitly.
- Covariance-correct equating drift tests, both families. The shift c0
  is estimated from the same common items the tests examine, and
  locations within a calibration are correlated through the
  identification constraint: drift denominators now use
  `[(I - 1u') Sigma (I - u1')]_ii` with Sigma the sum of the two
  calibrations’ location covariances (a bank contributes diag(se^2)).
  Null calibration: 4.6% rejection at nominal 5%, t-SD 0.97.
- Judge-clustered SEs refuse a single judge (they collapse to ~1e-16)
  and note fewer than 10 clusters as likely understating.
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  now tests at the judge level: residuals aggregate to one weighted mean
  per judge (per opponent band) and enter a split-plot ANOVA with the
  judge as the error unit. The comparison-level ANOVA pseudo-replicated
  – 6 of 10 null datasets with judge heterogeneity and arbitrary groups
  falsely flagged uniform DIF; the judge-level design flags 0 of 10
  while still detecting a planted 1-logit effect. Band breaks moved to
  count-weighted quantiles, making aggregated and expanded data exactly
  equivalent. Power now scales with judges per group – the honest unit.
- Non-convergence now warns loudly at fit time in
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  and [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  (it was only visible in print);
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)’s
  convergence flag is scale-free (gradient per comparison).
- Documentation: the test-of-fit suite’s calibration status is stated
  plainly (approximate, convention-faithful diagnostics; simulation
  shows near-but-not-nominal rejection), with parametric-bootstrap
  references via
  [`sim_replicate()`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md)
  recommended where exact calibration matters; the BTL pair chi-squares
  are marked descriptive under judge clustering; the next-pair priority
  is documented as a ranking heuristic; the BTL dimensionality
  reference’s no-refit design is stated.

## rasch 1.10.3

Input honesty, from an external code review (five findings, all
confirmed and fixed):

- Misspelled `id`, `factors`, or `items` column names in
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  are now errors naming the missing columns; they were silently ignored
  (one misspelled factor silently dropped ALL factors), producing
  valid-looking analyses of the wrong data. Numeric item indices still
  work.
- Fractional scores (e.g. 1.9) are now an error naming the columns;
  [`as.integer()`](https://rdrr.io/r/base/integer.html) used to truncate
  them silently, altering response data. Integer-valued doubles (“2.0”)
  still pass.
- MFRM rows with a missing person, item, or facet identifier are dropped
  with a note; [`paste()`](https://rdrr.io/r/base/paste.html) used to
  coerce `NA` to the literal string “NA”, silently pooling unrelated
  rows under a phantom level.
- [`equate_tests()`](https://drjoshmcgrane.github.io/rasch/reference/equate_tests.md)
  excludes common items whose location or SE is unavailable (e.g. weakly
  determined items with honest NA SEs) from the precision-weighted shift
  and drift tests, notes them, and keeps their rows with NA test columns
  – one such item used to turn every statistic NA.
  [`plot_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_equate.md)
  handles the excluded rows.
- [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  escapes all data-derived text (title, notes, item and facet names) –
  crafted column names could inject markup into reports; table cells
  were already escaped.

## rasch 1.10.2

- App: hover identification is now consistent across the app, through
  one shared mechanism. Beyond the person-fit plot and item fit map,
  hovering identifies the point (or cell) on: the Guttman scalogram
  (person, item, score), both residual-correlation heatmaps (item pair
  and its Q3/Q3\*), the common-item equating plot and the
  tailored-analysis comparison (item, reference and current locations),
  and – for paired comparisons – the unexpected-judgements map (matchup,
  surprise, location), the group characteristic curves’ observed means
  (object, panel, mean), and the equating plot (object, both
  calibrations). Every tooltip is resolved server-side against exactly
  the rows its plot draws. Plots whose points are already labelled
  (kidmap, Wright map, caterpillars) and curve or distribution displays
  are deliberately left without hover.

## rasch 1.10.1

- App: hovering a point on the person-fit plot or the item fit map now
  shows who it is – a floating tooltip with the person ID (or item
  name), location, and fit residual, resolved server-side against
  exactly the subset each plot draws. The navbar wordmark is
  left-justified.

## rasch 1.10.0

- Model comparison, calibrated for composite likelihood.
  [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  gains `cl_aic` and `cl_bic` (Varin & Vidoni 2005; Gao & Song 2010):
  `-2 cl` penalised by the effective parameter count `tr(H^-1 J)` from
  the Godambe matrices – the same quantity whose eigenvalues calibrate
  [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
  – which absorbs the pairwise over-counting a nominal AIC/BIC would
  ignore (empirically 3.5-5x the nominal count). Verified by
  plant-and-detect: the criteria select RSM on rating-structured data,
  PCM under varying threshold spreads, and a planted BTL position effect
  – and reject the richer model when the effect is absent. `n` for
  CL-BIC counts independent units: persons contributing an informative
  pair, or judges (count-weighted comparisons when unclustered).
  MFRM/EFRM fits report NA with the reason noted.
- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  now accepts `btl` fits (all-BTL comparisons: free vs
  principal-component thresholds, with/without a position effect or
  within-judge dependence), reporting judges/objects/comparisons and
  OSI. [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  fits carry their composite-likelihood ingredients in `$cl`; anchored
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  and
  [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  now return their Godambe matrices. The app’s Compare page shows the
  criteria and accepts paired-comparison fits.
- Cross-validation against independent implementations, in the test
  suite whenever the packages are installed (new in Suggests: eRm, sirt,
  psychotools): dichotomous, PCM, and RSM parameters against eRm’s full
  CML (agreement to sampling precision; pairwise SEs at or just above
  CML’s, the documented efficiency price), the same pairwise family
  against sirt::rasch.pairwise (near-identity), and btl() against
  psychotools::btmodel to machine precision (same likelihood).

## rasch 1.9.3

- Missing-data identification, verified and enforced. A booklet-design
  probe (zero complete cases, linked forms, MCAR on top) confirmed the
  pairwise conditional estimator recovers planted parameters with
  nominal coverage – and exposed two silent failure modes at the
  design’s edges, both now honest:
- Disconnected designs error instead of fitting. When no person answers
  items in more than one block, the pairwise likelihood is flat in the
  between-block shift; the fit used to return arbitrary cross-block
  spacing (within-block orders fine, origins meaningless) with no
  warning.
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md),
  [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md),
  and everything built on them now stop with an error naming the
  disconnected blocks. Anchors rescue a block (disjoint-form equating):
  with an anchored item in every block the fit proceeds; with any block
  unanchored it stops and names it.
- Thresholds beside near-empty categories report NA, not a lie. A
  category with 1-2 responses can send its threshold toward a boundary
  (observed: a single top-category response put an item location at
  +13.9 logits) while the ridged covariance printed SE = 0. Such
  thresholds are now flagged (`thr$weak`), their SEs and the
  item-location SE reported as NA, and a note names the item, the
  category, and its count, pointing to `pc_components` or category
  collapsing. Estimates are still returned; categories with zero
  responses are rescored away as before.

## rasch 1.9.2

- Case study: `inst/casestudies/party_blocs_crisis.R` applies
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  to the Tuebingen 2009 party-preference data (via `psychotools`), with
  ideological blocs as object sets and economic-crisis concern as judge
  panels. Running the model against real data surfaced three defects,
  all fixed below.
- A set whose within-set contests are all near-even (or all one-sided)
  carries no stable information about the panel-unit ratios: its
  stage-one ratio estimate diverged and poisoned the phi reconciliation
  with a spurious covariance (on the party data, the right bloc’s single
  CDU/CSU-FDP pair at 55:45 drove phi to 3e4 and a fatal error). Such
  sets are now screened out of the reconciliation, refit with the panel
  units held at the reconciled phi – which the frame model says apply to
  them regardless – and named in a note. Stage one also gains the same
  trust-region step cap as
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).
- Bootstrap replicates in which a set unit reaches its boundary (log
  alpha driven to -Inf when a resampled within-set order flips against
  the cross-set evidence, the signature of a two-object set with a
  near-even internal pair) made the standard error silently NaN, printed
  as blank. The SE is now reported as NA with a note naming the
  parameter and the boundary count: no number is manufactured for a
  sampling distribution that is not normal. A replicate whose
  convergence flag was NA could also crash the bootstrap loop; it now
  counts as a failure.
- The stage convergence checks compared the gradient to an absolute
  threshold, which is scale-dependent: on large designs a converged fit
  could be flagged unconverged and – with the new screen – silently
  rerouted, changing the estimates. Both stages now use a per-comparison
  criterion, verified invariant to 50-fold duplication of the data.
- Documentation: the bootstrap’s two honesty notes (model-based
  resampling versus judge-clustered conditional errors, and
  boundary-unstable parameters) are spelled out in
  [`?btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md).

## rasch 1.9.1

- Shiny app: the Frames page gains its paired-comparison variant –
  nominate a judge-panel column and an object-set map (CSV, or inferred
  from object-name prefixes), estimate the frame units on demand
  (bootstrap standard errors by default), and read the panel units, set
  units and origins, the unit caterpillar, and the per-frame fit. On the
  app’s own paired-comparison demo the two judge panels differ
  significantly in discriminating power – the panel containing the
  erratic judge carries the smaller unit – which is precisely the
  phenomenon the model makes visible.

## rasch 1.9.0

- The extended frame of reference model for paired comparisons:
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  estimates object locations when the unit of the scale differs across
  judge-panel by object-set frames – the package’s extension of
  Humphry’s persons-by-items EFRM to the Bradley-Terry-Luce family.
  Panel units (phi, a panel’s discriminating power) are identified
  conditionally wherever two panels judge a common set; set units and
  origins (alpha, kappa) are identified from cross-set comparisons,
  which are themselves conditional outcomes – so, unlike the
  persons-by-items case, the linking makes no distributional assumption
  about the objects. Everything is estimated by the staged conditional
  estimator: the within-frame calibrations are invariant to the linking
  data (the frame-of-reference property), a deliberate trade of
  efficiency for invariance. The cross-frame discrimination convention
  (cross-set contests judged at exactly phi) is stated plainly in the
  documentation as an assumption, with the per-frame fit residuals as
  its check.
- Inference for
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  defaults to a parametric bootstrap of the whole two-stage pipeline:
  simulation showed the fast conditional errors cover at a third of
  nominal on chain-linked designs (they omit the stage-one uncertainty
  that flows into the linking), while the bootstrap restores 12/12
  coverage and nominal false-flag rates at both small and large designs.
  The conditional errors remain available as `se_method = "conditional"`
  for quick inspection, labelled as understating.
  [`plot_btl_units()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_units.md)
  draws the unit caterpillar;
  [`simulate_btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl_efrm.md)
  generates frame-structured data with planted units for the
  plant-to-detect loop.

## rasch 1.8.0

Three further paired-comparison analogues of the Rasch tool set, and
CRAN housekeeping.

- Anchored estimation for paired comparisons:
  `btl(anchors = c(Object = value, ...))` holds the named objects at
  fixed locations and estimates the rest on that scale (no sum-zero
  constraint; anchored objects report se 0) – the equating workflow’s
  second half, mirroring `rasch(anchors=)`. An anchored object at a
  response boundary is an error rather than a silent removal.
- First-position advantage: `btl(position = TRUE)` estimates the pure
  positional (first/left-presented) advantage of the Davidson & Beaver
  1977. order-effect device as a constant covariate alongside any
        exposure and carry-over effects; it is identified by triangle
        closure even under fixed presentation orders (more weakly, with
        honest SEs). The dimensionality reference simulation carries the
        fitted coefficient.
- Common-object equating:
  [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  compares two paired-comparison calibrations (or a calibration against
  a bank) through their common objects – precision-weighted origin
  shift, per-object drift t-tests with multiplicity adjustment, equated
  locations – with
  [`plot_btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_equate.md);
  the standards-maintenance workflow of comparative judgement across
  panels or years.
- Targeting and adaptive pairing:
  [`btl_information()`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md)
  (per-comparison Fisher information and per-object design information,
  with the naive 1/sqrt(information) SE beside the judge-clustered
  sandwich),
  [`plot_btl_targeting()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md),
  and
  [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  – the greedy most-informative-next-pair step of adaptive comparative
  judgement (Pollitt 2012), with the Bramley (2015)
  reliability-inflation caution documented.
- CRAN housekeeping: pure-ASCII sources, quoted ‘shiny’ in DESCRIPTION,
  heavy Monte-Carlo recovery tests wrapped in `skip_on_cran()` (they run
  locally and on CI), and `cran-comments.md` drafted from a full
  `--as-cran` check (about 3 minutes; slowest example 1.4s).
- An adversarial verification round over the new features also caught
  and fixed a PRE-EXISTING inference bug: `btl(count=)` inflated every
  sandwich standard error by about the square root of the count (the
  meat weighted aggregated rows twice; point estimates were always
  correct). Count- weighted fits now reproduce the expanded-data
  standard errors exactly; OSI and all z/p values follow. Additionally
  from that round:
  [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)’s
  shift standard error now uses the stored covariance of each
  calibration (it treated the common-object differences as independent,
  overstating the SE several-fold; 95 per cent coverage verified) and
  documents the majority-drift limitation of single-shift equating;
  [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  ranks by the exact one-step reduction in total location variance (a
  Sherman-Morrison update on the fit’s covariance), verified to beat
  low-priority pairs empirically; and
  [`btl_information()`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md)’s
  `se_naive` is documented honestly as a single-parameter lower bound
  rather than an independence benchmark.
- Shiny app: the Targeting and Equating tabs gain paired-comparison
  variants (design information + targeting plot + adaptive next-pair
  recommendations; reference-calibration upload with drift table and
  plot), and the Data roles gain a first-position-advantage switch and
  an object-anchors upload.

## rasch 1.7.1

A statistical audit of the new simulation and paired-comparison
independence code (23 findings across two adversarial rounds, all
fixed), plus a vignette and a logo.

- Misfit layers now compose correctly in
  [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md):
  the response-dependence and response-style regeneration steps
  previously rebuilt their target items from the base structure,
  silently erasing any DIF or second-dimension loading planted on those
  items (and leaking a DIF-like signal from a DIF-carrying source item
  into its dependent partner). Every regeneration now honours the item’s
  own generating structure – trait, group-shifted thresholds and slope,
  guessing.

- [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  now simulates its noise reference sequentially WITH the fitted
  within-judge dependence effects when the fit carries them: real order
  effects previously inflated the false-positive rate to ~36 per cent on
  one-dimensional data (now ~0-5 per cent). The documented price is
  conservatism, since carry-over and a judge-camp second attribute are
  partially confounded.

- Doc examples pass `id = "id"` so they no longer emit spurious
  dropped-column notes.

- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md):
  `model = "PCM"` now draws each item’s threshold spacings and span
  afresh, as the partial credit model allows (previously PCM and RSM
  generated identical rating-scale-structured data); `dif` without
  `n_groups >= 2` is an error instead of a silent no-op; polytomous
  `guessing` warns and is ignored rather than falsely recorded as
  planted; the response-dependence term now carries the guessing
  asymptote of the pair’s first item (no more spurious easiness on the
  dependent item); disordered thresholds work with 3 categories and no
  longer shift the item’s location; careless responders are no longer
  double-counted as response-style persons in the recorded truth.

- [`simulate_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_mfrm.md):
  halo raters are capped at the eligible pool (no NA raters when erratic
  raters shrink it); the example runs.

- [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md):
  person abilities are mean-centred as documented (an asymmetric
  difficulty range no longer masquerades as person bias); the many-facet
  branch reads the item margins (`item_effects`), restoring the
  documented item-difficulty recovery for MFRM fits.

- [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md):
  null calibration verified for graded fits (mildly conservative, never
  anticonservative) and documented; docs also clarify the
  leading-bimension reference columns and the no-ties precondition of
  Kendall’s zeta.

- Shiny app (Simulate page): the parameter-recovery card only renders
  when the current fit was estimated on the currently loaded simulation
  (a stale fit can no longer be compared against new truth); the
  reproducible call now carries every argument the generator received
  (MFRM person/item SDs, EFRM person SD, graded categories); non-integer
  seeds are rounded safely.

- New vignette: “Planting misfit and watching the diagnostics fire” —
  the plant-to-detect loop end to end, closing with a small Monte Carlo
  power estimate via
  [`sim_replicate()`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md).

- A hex logo (an item characteristic curve at its location),
  reproducible from `tools/logo.R`.

## rasch 1.7.0

- Simulation gains the standard controls of a proper simulation tool.
  - Population parameters: the person distribution (mean, SD, shape) and
    the item-difficulty range are now set directly (in the package and
    the app).
  - More misfit: response styles (extreme / midpoint categories) and
    speededness (a not-reached tail) for , and a halo effect for .
  - generates a batch of datasets for Monte Carlo power, Type-I, or
    recovery studies.
  - and compare the parameters a fit recovers with the ones that were
    planted (correlation, RMSE, bias, and a true-vs-estimated scatter)
    for every layout.
- Shiny app: the Simulate page renames the wide layout to “Rasch”,
  exposes the population parameters, shows a parameter-recovery report
  and the true-vs-estimated scatter after Run, and prints the
  reproducible call.

## rasch 1.6.1

- Shiny app: a Simulate page (under More) exposes all four simulators
  with a data-type selector and per-layout misfit controls. Pressing
  Simulate generates the data, loads it as the current dataset with its
  roles set, and lists what was planted; pressing Run then analyses it
  so the matching diagnostic can be watched firing. Verified end to end
  for every layout (e.g. a planted DIF item is flagged by dif_anova, a
  many-facet fit recovers its rater severities).

## rasch 1.6.0

- Data simulation across all four model families, each with dial-in
  departures from the model – so a known pathology can be planted and
  the matching diagnostic watched as it fires. Every result carries the
  true parameters () and prints what was planted.
  - : dichotomous/PCM/RSM with discrimination misfit, guessing, a second
    dimension, local dependence, DIF (uniform and non-uniform), careless
    responders, disordered thresholds, and missing data.
  - : paired comparisons with erratic judges, a second object attribute
    (two judge camps), and within-judge exposure/carry-over dependence;
    dichotomous or graded.
  - : rated data with rater-severity spread, erratic raters, and a
    rater-by-item interaction.
  - : frames with differing set and group units. Verified end to end:
    each planted departure is recovered or flagged by the corresponding
    fit / diagnostic (an app panel follows).

## rasch 1.5.1

- The judge unexpected-judgements map () is now drawn at the matchup
  level: each pair the judge met is a segment spanning its two objects
  on the consensus location axis, placed by how surprising the verdict
  was, with the judge’s backed object filled. returns the per-matchup
  residuals behind it. The object-level is retained for per-object
  residuals.

## rasch 1.5.0

- [`judge_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md)
  and
  [`plot_btl_judge_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_judge_map.md):
  the paired-comparison counterpart of the kidmap. A judge has no
  ability to condition on, so the reference is the consensus object
  scale; each object a judge met gets a standardised residual saying how
  much more or less that judge favoured it than its location predicts. A
  strong object the judge under-rated, or a weak object over-rated, is
  an unexpected judgement. It catches systematic bias – a distinct
  failure mode from the noise that judge fit and judge consistency flag.
  In the Shiny app the Persons \> Judge fit panel is now a
  master-detail: the judges table drives the selected judge’s map,
  mirroring person -\> kidmap.

## rasch 1.4.0

- Paired comparisons gain the pair-structure analogues of the Rasch
  independence diagnostics – the tools that a persons-by-items residual
  matrix would give but paired-comparison data cannot.
  - [`btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md)
    counts circular triads (preference loops: A beats B, B beats C, C
    beats A) against the chance rate, reports Kendall’s coefficient of
    consistency for complete designs, and gives each judge’s own
    consistency (Kendall & Babington Smith 1940) – the single-dimension
    check, and a judge-fit analogue.
  - [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
    decomposes the skew-symmetric object-by-object residual preference
    matrix into Gower (1977) bimensions – the paired-comparison
    counterpart of residual PCA – and judges the leading bimension
    against a model-simulated noise reference, so a coherent residual
    “swirl” (a second attribute steering some contests) is distinguished
    from noise.
  - [`plot_btl_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_scree.md),
    [`plot_btl_dim_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dim_map.md),
    and
    [`plot_btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_transitivity.md)
    display them; the Shiny app’s Independence \> Trait tab now carries
    these for a BTL fit (it previously hid for paired comparisons).

## rasch 1.3.1

- The package is now called **rasch** (it was developed under the
  working title `rmt`). The name `rmt` belongs to an unrelated package
  on CRAN, and `rasch` says plainly what this is – and mirrors
  `mirt::mirt()` with
  [`rasch::rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).
  No functions, arguments, or results change; only the package and
  library name. The auxiliary result classes were unified under the
  `rasch_` prefix (`rasch_btl`, `rasch_dif`, …) to match the fit
  classes.

## rasch 1.3.0

A full statistical audit of the independence, invariance, and paired-
comparison procedures, with fixes throughout.

- Paired comparisons:
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) now
  interrogates within-judge dependence – `dependence_data` holds every
  comparison’s exposure and carry-over covariates (count-weighted),
  [`plot_btl_dependence()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dependence.md)
  displays an effect as a partial-residual curve, and a separated effect
  (one-sided informative comparisons) is set aside with a note instead
  of reported as a runaway estimate. The graded engine’s
  threshold-by-dependence Hessian block is corrected (free-threshold
  fits with five or more categories and an order column now converge),
  rows with a missing judge are dropped with a note, category checks
  re-run after boundary objects are removed, and the dichotomous path
  now routes through the one (vectorised) estimator.
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  holds fitted dependence effects fixed in its residual moments and
  resolution refits (dependence is no longer absorbed as spurious
  judge-group DIF), weights count-aggregated comparisons exactly as
  expanded rows, restricts supersession to group terms, and notes every
  significant term it does not resolve.
- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  keeps each ANOVA term in its own error stratum, so a mixed
  (within-subject) design with missing responses no longer un-flags real
  DIF; factors named like internals (`ci`, `f1`) are safe; a term’s own
  class-interval crossing no longer supersedes it (so
  [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  again splits items with both uniform and non-uniform DIF); mixed
  designs say plainly that Tukey comparisons are unavailable; and
  [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  notes when a within-person factor makes its standard errors
  conservative.
- [`plot_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)’s
  parallel-analysis reference is now model-simulated (responses drawn
  from the calibrated model, persons re-estimated), so model-true data
  sits at the reference instead of always “showing structure” (Raiche
  2005; Chou & Wang 2010).
- Shiny app: judge-group DIF results freeze their run’s factors and
  alpha (later sidebar edits cannot silently regroup the overlay or
  break the reproducible snippets), superseded terms are shown, the
  within-judge dependence panel explains dropped effects, and fit flags
  are consistent on every table (fit residual \|2.5\|; outfit 0.7-1.3;
  infit 0.8-1.2).

## rasch 1.2.0

- [`plot_pca_biplot()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pca_biplot.md)
  draws the item loadings on the first two residual principal components
  – the pair that usually carries any interpretable contrast – on equal
  axes, coloured by the sign of the PC1 loading.

- [`residual_correlations()`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md)
  now also returns the adjusted-Q3 `star_matrix` (each Q3 less the
  average off-diagonal), and
  [`plot_resid_cor()`](https://drjoshmcgrane.github.io/rasch/reference/plot_resid_cor.md)
  gains a `stat` argument to colour either the raw Q3 or the adjusted
  Q3\*, drawing the lower triangle only.

- Shiny app: the trait- and local-dependence pages pair each table with
  its plot (loadings with a PC1-PC2 biplot; Q3 and Q3\* matrices with
  heatmaps), the raw and adjusted correlations are each flagged by their
  own rule (\|Q3\| and Q3\*), and every data-restructuring action
  (subtest, item split, automatic DIF resolution) offers an in-place
  reset to the original data.

- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  is now the single DIF analysis-of-variance function. One factor is
  analysed one-way; several factors are modelled jointly – the
  statistically correct treatment – with main effects by default and
  factor-by-factor interactions optional (`effects = "factorial"`). It
  handles within-subject factors (repeated-measures / mixed ANOVA) and
  returns a classed object with a `summary`, the full `terms` table, and
  Tukey comparisons. The separate `dif_anova_factorial()` is removed.

- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  resolves DIF iteratively by item splitting, largest effect first, to
  clear artificial DIF.

## rasch 1.0.0

First stable release. The package delivers a complete Rasch Measurement
Theory workflow built entirely from published measurement theory, with a
pairwise conditional estimation core and a Shiny interface that exposes
every analysis with reproducing R code attached to every output.

### Models

- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md):
  dichotomous and polytomous (partial credit and rating scale) Rasch
  models by pairwise conditional maximum likelihood (Andrich & Luo 2003;
  Zwinderman 1995), with Godambe sandwich standard errors, sum-zero
  identification, Warm (1989) weighted likelihood person estimation, and
  the principal-component threshold parameterisation as an estimation
  option.
- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md):
  many-facet models with additive facet severities or item-by-facet
  interactions, accepting wide (one column per item) or long data.
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md):
  the extended frame-of-reference model (Humphry) with frame units
  estimated by within-frame pairwise conditional likelihood.
- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md):
  paired comparisons as the conditional form of the dichotomous Rasch
  model (Bradley & Terry 1952; Andrich 1978), including the
  adjacent-categories graded extension (Tutz 1986; Agresti 1992) with
  symmetric thresholds – the Davidson (1970) ties model is its
  three-category case – winner-plus-margin data entry, judge-clustered
  errors, and within-judge exposure and carry-over effects estimated
  jointly with the locations.

### Diagnostics

- The test-of-fit suite follows Andrich & Marais (2019):
  log-of-mean-square fit residuals with apportioned degrees of freedom,
  the item-trait chi-square on automatically sized, tie-preserving class
  intervals with per-interval detail, the ANOVA item fit, and person
  fit.
- Unidimensionality: residual principal components with parallel
  analysis, the Smith (2002) subset t-test, and magnitude estimation.
- Local dependence: Yen’s (1984) Q3 with the Christensen, Makransky &
  Horton (2017) flagging convention, response-dependence magnitude
  (Andrich & Kreiner), subtest formation, and the spread test.
- Guessing: tailored analysis with origin-equated comparison
  calibrations.
- DIF: two-way residual ANOVA per factor, the joint factorial analysis
  with a compact uniform/non-uniform summary, resolved DIF magnitudes in
  logits with familywise control and a practical-significance criterion,
  planned contrasts derived from the factor structure (with
  repeated-measures support via person-level residual scores), and class
  intervals sized to the cells each analysis actually uses. For paired
  comparisons,
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  tests object-by-judge-group interaction by the same two routes.
- Equating and invariance: common-item comparison with drift flags, and
  calibration anchoring.

### Display and reporting

- A complete base-graphics plot suite: expected value curves with
  observed class-interval means, category and threshold probability
  curves, the person-item threshold distribution, Wright maps, kidmaps
  (Wright, Mead & Ludlow 1980), object characteristic curves for paired
  comparisons, threshold and item maps, test characteristic and
  information curves, residual displays, and Guttman scalograms – all
  with adjustable class intervals and scale ranges, and batch export to
  multi-page PDF or ZIP at publication resolution.
- [`fit_summary_table()`](https://drjoshmcgrane.github.io/rasch/reference/fit_summary_table.md)
  and
  [`targeting_table()`](https://drjoshmcgrane.github.io/rasch/reference/targeting_table.md)
  return the headline statistics as tidy tables;
  [`save_outputs()`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md)
  writes every table and plot to disk;
  [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  produces a single-file HTML report.
- [`run_app()`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md)
  launches the Shiny interface: model-adaptive navigation, reproducing R
  code beneath every output, and one-click export.
