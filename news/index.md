# Changelog

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
