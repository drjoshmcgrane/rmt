# CRAN comments for rasch 1.11.2

## Note on this submission

This submission replaces rasch_1.10.2, currently awaiting review in the
incoming queue: after that upload, an external code review found
statistical errors in the package (a covariance not transformed to the
recentred parameterisation on mixed max-score designs, a spurious
discrimination factor in the Warm weighted-likelihood correction,
degrees of freedom assigned to untestable items, drift tests ignoring
the estimated equating shift, and a pseudo-replicating judge-group DIF
test). Every finding was verified by simulation, fixed, and locked in
by calibration tests; we would rather the reviewed and published first
release be the corrected one. Apologies for the replacement.

## Summary

rasch provides pairwise conditional maximum likelihood estimation of
dichotomous and polytomous Rasch models (partial credit and rating scale),
with a comprehensive diagnostic suite following Andrich & Marais (2019)
(fit residuals, item-trait interaction, dimensionality, DIF, local
dependence), anchored equating, the many-facet Rasch model, the
Bradley-Terry-Luce model for paired comparisons, the extended frame of
reference model (Humphry, 2005) in both its persons-by-items and
paired-comparison forms, model comparison by composite-likelihood
information criteria, and a data simulator for every model family. It is
implemented from published measurement theory in base R, with no
dependence on other estimation engines, and includes a 'shiny' interface
and a one-call exporter for every table and plot.

References in the Description carry DOIs where the venue assigns them;
the remaining citations (Journal of Applied Measurement articles, books,
and theses) have no DOI.

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.1
* GitHub Actions (r-lib actions): macos-latest (release),
  windows-latest (release), ubuntu-latest (devel, release, oldrel-1)
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE
  New submission. This is the first submission of this package to CRAN.

* checking HTML version of manual ... NOTE (local only)
  "Skipping checking HTML validation: 'tidy' doesn't look like recent
  enough HTML Tidy." This reflects the local machine's 2006-era macOS
  HTML Tidy, not the package; it is not expected on CRAN's check machines.

### Possibly mis-spelled words in DESCRIPTION

All flagged words are correctly spelled: surnames of the statisticians
whose methods are implemented (Rasch, Andrich, Luo, Zwinderman, Godambe,
Pedler, Warm, Marais, Kreiner, Humphry, Linacre, Guttman, Varin, Vidoni)
and standard psychometric terminology (reparameterises [British spelling,
used consistently], infit, familywise, DIF, multidimensionality, subscale,
subtest, distractor(s), rescoring, scalogram).

## Timings

The full --as-cran check runs in about 4 minutes locally; the slowest
example is under 2s. Heavy Monte-Carlo recovery and bootstrap-calibration
tests are wrapped in skip_on_cran() (they run locally and on CI, where the
suite is 1060+ assertions).

## External validation

The estimators are cross-validated in the test suite against independent
implementations whenever those packages are installed (all in Suggests):
eRm (full conditional likelihood), sirt (pairwise), and psychotools
(Bradley-Terry). These tests are skipped via skip_if_not_installed() when
the packages are absent.

## Downstream dependencies

This is a new package; there are no reverse dependencies.
