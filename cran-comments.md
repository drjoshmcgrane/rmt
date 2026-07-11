# CRAN comments for rasch

## Summary

rasch provides pairwise conditional maximum likelihood estimation of
dichotomous and polytomous Rasch models (partial credit and rating scale),
with a comprehensive diagnostic suite following Andrich & Marais (2019)
(fit residuals, item-trait interaction, dimensionality, DIF, local
dependence), anchored equating, the many-facet Rasch model, the
Bradley-Terry-Luce model for paired comparisons, and the extended frame of
reference model (Humphry, 2005). It is implemented from published
measurement theory in base R, with no dependence on other estimation
engines, and includes a 'shiny' interface and a one-call exporter for
every table and plot.

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.1
* GitHub Actions (r-lib actions): macos-latest (release),
  windows-latest (release), ubuntu-latest (devel, release, oldrel-1)

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
Pedler, Warm, Marais, Kreiner, Humphry, Linacre, Guttman) and standard
psychometric terminology (reparameterises [British spelling, used
consistently], infit, familywise, DIF, multidimensionality, subscale,
subtest, distractor(s), rescoring, scalogram).

## Timings

The full --as-cran check runs in about 3 minutes locally; the slowest
example is 1.4s; a handful of heavy Monte-Carlo recovery tests are wrapped
in skip_on_cran() (they run locally and on CI).

## Downstream dependencies

This is a new package; there are no reverse dependencies.
