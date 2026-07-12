# Compare fitted Rasch models

Builds a comparison table for two or more fits from
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
[`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
or (all together)
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md). For
fits of the same response data (identical item columns, maximum scores,
and number of persons) the pairwise conditional log-likelihoods share
their conditional information, and twice the difference from the
reference fit is reported with the difference in parameter counts; this
is descriptive (composite likelihood), and most meaningful for nested
structures such as RSM inside PCM.

## Usage

``` r
compare_fits(..., reference = 1)
```

## Arguments

- ...:

  Two or more fitted objects, ideally named
  (`compare_fits(PCM = f1, RSM = f2)`). Either all Rasch-family fits or
  all `btl` fits; for `btl`, fits of the same comparison data (same
  objects, comparisons, and judges) support the likelihood columns –
  e.g. free versus principal-component thresholds, with and without a
  position effect or within-judge dependence.

- reference:

  Index or name of the reference fit for the log-likelihood difference;
  defaults to the first.

## Value

A data frame with one row per fit: label, model, persons, items (judges,
objects, comparisons for `btl`), parameters, log-likelihood,
`eff_params`, `cl_aic`, `cl_bic`, comparability with the reference,
`two_delta_ll` and `delta_parameters` (same-data fits only), chi-square
per df, fit residual SDs, PSI, and alpha (OSI for `btl`).

## Details

The calibrated comparison is carried by the composite-likelihood
information criteria `cl_aic` (Varin and Vidoni 2005) and `cl_bic` (Gao
and Song 2010): \\-2\\cl + c \cdot tr(H^{-1}J)\\, with \\c = 2\\ or
\\\log n\\. Because every response enters every pair its item forms, the
pairwise log-likelihood over-counts the data; the effective parameter
count \\tr(H^{-1}J)\\ from the Godambe matrices – the same quantity
whose eigenvalues calibrate
[`lr_test`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
– absorbs exactly that over-counting, where the nominal parameter count
would not. \\n\\ counts independent units: persons contributing at least
one informative pair, or judges for paired-comparison fits
(count-weighted comparisons when unclustered). Smaller is better; the
criteria are valid across models of the same data whether or not they
nest, and are `NA` (with the reason in the printed note) for MFRM and
EFRM fits, which do not carry their Godambe matrices.

Across different data preparations (subtests, splits, facet or frame
structures) the likelihoods are not comparable and the calibration-free
columns carry the comparison: total item-trait chi-square per degree of
freedom, item and person fit residual SDs (ideal 1), PSI, and alpha (OSI
for paired comparisons).

## Examples

``` r
set.seed(1)
simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
th <- rnorm(400)
X <- sapply(seq(-1, 1, length.out = 6), function(b)
  sapply(th, function(t) sample(0:3, 1, prob = simP(t, b + c(-0.8, 0, 0.8)))))
colnames(X) <- paste0("R", 1:6)
compare_fits(PCM = rasch(X, model = "PCM"), RSM = rasch(X, model = "RSM"))
#> Model comparison (reference: PCM)
#> 
#>  label model persons items eff_params   cl_aic   cl_bic two_delta_ll
#>    PCM   PCM     400     6     61.697 8108.103 8353.432             
#>    RSM   RSM     400     6     21.817 8051.228 8137.980      -22.886
#>  chisq_per_df item_fit_sd person_fit_sd   PSI alpha
#>         0.829       0.556         0.835 0.745 0.784
#>         0.884       0.600         0.837 0.743 0.784
#> (further columns on the object: loglik, parameters, same_data)
#> 
#> cl_aic and cl_bic are composite-likelihood information criteria (Varin & Vidoni 2005; Gao & Song 2010): -2 cl penalised by the effective parameter count tr(H^-1 J), which absorbs the pairwise over-counting that the nominal count would not; smaller is better, valid across models of the same data. two_delta_ll is the raw composite difference against the reference, descriptive only. Across different data preparations compare chisq_per_df, the fit residual SDs (ideal 1), and the separation/reliability columns.
```
