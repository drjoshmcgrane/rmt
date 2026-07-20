# Fit the extended frame of reference model

Estimates Humphry's extended frame of reference model, in which the unit
of the latent scale differs across frames (item-set by person-group
cells): the response model is
`P(X = x) prop exp(rho_sg (x theta - sum delta))` with
`rho_sg = alpha_s phi_g`. Within frames the partial credit model holds
in the frame's natural unit, so item thresholds and the person group
units `phi` are estimated by within-frame pairwise conditional maximum
likelihood (the person parameter cancels; Andrich and Luo 2003), jointly
across frames through the sets shared by several groups. Item-set units
`alpha` and set locations are then estimated from persons common to
pairs of sets, using error-corrected true-score variances, and
reconciled over the linking graph by weighted least squares. Everything
is reported in a common arbitrary unit, and the returned object is also
a full
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) fit
at the item-by-group level, so the package's diagnostic tables and plots
apply.

## Usage

``` r
rasch_efrm(
  data,
  item_sets,
  groups,
  id = NULL,
  factors = NULL,
  items = NULL,
  n_groups = NULL,
  adjust_N = NA,
  na_codes = -1,
  maxit = 50,
  tol = 1e-07,
  min_link_persons = 30,
  se_method = c("hybrid", "bootstrap"),
  boot_reps = NULL
)
```

## Arguments

- data:

  Persons-by-items data (matrix or data frame, like
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)),
  plus a person-group column.

- item_sets:

  A named list mapping set names to item-column names, or a named
  character vector mapping item names to set names. Items not mentioned
  form their own set `"(rest)"` when a list is given.

- groups:

  Name of the person-group column in `data`, or a vector with one entry
  per person. Several column names may be given: the frames are then
  their crossed cells, per-cell units appear in `phi_table`, and a
  factorial decomposition of the cell units (sum-coded main effects, and
  the interaction when every cell is observed) is returned in
  `phi_factorial`. The decomposition is a generalised least-squares fit
  of the cell log-units using their joint covariance (bootstrap
  replicates when available, otherwise the analytic centred covariance,
  inverted spectrally along its identified directions); coefficient rows
  are descriptive, and inference is carried by the
  multi-degree-of-freedom Wald test per term in `phi_factorial_tests`.
  Group units are checked for identification on the joint information: a
  flat direction along a unit (structural non-identification) is refused
  with an error naming the group, since every common-unit quantity would
  silently depend on it. A unit whose analytic standard error exceeds 5
  log-units (uncertain beyond a factor of about 150) is practically
  uninformative but not structurally unidentified: its estimate is kept
  for sensitivity work, with a warning and a note. Weakly identified
  units with real threshold spread are kept, with standard errors that
  say how weak they are.

- id, factors, items, n_groups, adjust_N, na_codes:

  As in
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- maxit, tol:

  Outer iteration cap and convergence tolerance of the bilinear pairwise
  stage.

- min_link_persons:

  Minimum number of common persons required for a set pair to contribute
  to the unit linking.

- se_method:

  `"hybrid"` (sandwich + linking bootstrap + delta propagation; fast,
  default) or `"bootstrap"` (full person bootstrap of all stages).

- boot_reps:

  Bootstrap replicates; defaults to 300 for the linking bootstrap and
  200 for the full bootstrap.

## Value

An object of classes `"rasch_efrm"` and `"rasch"`. In addition to the
standard components (computed over item-by-group virtual columns with
the frame units carried in `disc`), it has `frames` (one row per frame:
units, origin, pooled fit; under `se_method = "bootstrap"` the
frame-level `se_log_rho` comes from the joint replicate draws of
`log(alpha) + log(phi)`, capturing cross-stage dependence, while the
hybrid fallback combines the stagewise errors as if uncorrelated),
`phi_table`, `alpha_table`, `set_table`, `item_arbitrary` and
`thresholds_arbitrary` (the structural parameters in the common unit),
`score_curves` (per-group score-to-measure curves, replacing the
raw-score table), `efrm_vs_rasch` (fit comparison against the equal-unit
model on the same conditional information), and `linking` (the linking
evidence).

## Details

Humphry (2005) states the model for dichotomous responses. The
polytomous form fitted here, with the frame unit multiplying the whole
exponent over the item's partial-credit thresholds, is this package's
extension of that statement. It is the form characterised by preserving
the two properties the model's logic rests on: the partial credit model
holds within every frame in the frame's natural unit (so the pairwise
conditional cancellation remains valid), and the weighted score remains
sufficient for the person parameter. It reduces exactly to the
dichotomous model when items are scored 0/1 and to the ordinary partial
credit model when all units equal one. One interpretive consequence:
category widths in natural units scale with the frame unit, so a
high-unit frame makes proportionally sharper category distinctions;
frame-level fit and the per-frame category curves are where a violation
of this would appear.

Estimation order: the within-frame pairwise stage establishes the
centred set thresholds and the person-group units `phi`; the person-side
linking stage then establishes the item-set units `alpha` and set
locations. The reported item parameters and all person measures are
computed only after every unit is established: item thresholds are
mapped into the common arbitrary unit using `alpha` and the set
locations, and person measures are weighted-score weighted likelihood
estimates evaluated under the final units `rho = alpha * phi`. The
per-frame person estimates used inside the linking stage are interim
quantities for the unit ratios only and are discarded. The within-frame
stage needs no re-estimation once `alpha` is known, because the pairwise
likelihood is invariant to the within-set rescaling that `alpha`
represents; the units' own uncertainty is reported in `alpha_table` and
`phi_table` and folded into the common-unit standard errors as described
below.

Standard errors: under `se_method = "hybrid"` (default) the group units
carry sandwich standard errors from the pairwise stage, the set units
carry person-bootstrap standard errors from the linking stage, and the
unit uncertainty is propagated into the common-unit threshold and item
standard errors by the delta method (treating the item-side and
person-side information as independent). Under `se_method = "bootstrap"`
all stages are re-estimated on `boot_reps` person resamples and every
standard error and the threshold covariance come from the replicate
spread; slower, but captures all cross-dependencies jointly.

Relation to Humphry (2005, sections 5.3 and 5.4): the two-stage
architecture implemented here operationalises the estimation approach
proposed in section 5.3 of the thesis (conditional estimation within
frames, with the item-set units obtained through person estimates from
linked sets), and retains the thesis's error-variance correction (its
equation 2.29, after Andrich 1982) and its transformation of standard
errors into the common unit by the inverse discrimination. The standard
errors themselves go further than the thesis's section 5.4, which
inverts each diagonal element of the joint-likelihood information
separately and therefore conditions on the remaining parameters,
including the person locations, being treated as known. Here full
covariance matrices are used throughout; the item-side covariance
carries the Godambe sandwich correction required for a pairwise
composite likelihood; the unit uncertainty that the thesis's
transformation treats as fixed is propagated into the common-unit
parameters; and resampling replaces analytic plug-in variances for the
person-side linking stage.

Measurement-theoretic status: within every frame the model is strictly
Rasch, with person-free item comparisons by conditioning. Across frames
it is an argued extension of the theory of the unit (Humphry 2005;
Humphry and Andrich 2008): on this account the unit was always a
frame-dependent empirical property that the ordinary model leaves
implicit, and the extension makes it explicit; the orthodox reading of
Rasch measurement contests this, and applied reports should present it
as an extension rather than settled doctrine. Two concessions are
intrinsic to the model rather than to this implementation: the item-set
units are identified only from the person side (their conditional
identification is impossible, as documented above), so that step uses
distributional information; and person measures rest on weighted-score
sufficiency with estimated weights, whose uncertainty is propagated
rather than ignored.

## Examples

``` r
# \donttest{
set.seed(1); Np <- 400
simP <- function(th, tau, r) { x <- 0:length(tau)
  p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
grp <- rep(c("A", "B"), each = Np / 2)
phi <- c(A = 0.8, B = 1.25)
d <- seq(-1.5, 1.5, length.out = 10)
X <- sapply(seq_along(d), function(i) sapply(seq_len(Np), function(n)
  sample(0:1, 1, prob = simP(rnorm(Np)[n], d[i], phi[grp[n]]))))
colnames(X) <- sprintf("I%02d", seq_along(d))
fit <- rasch_efrm(data.frame(X, grp = grp), item_sets = list(core = colnames(X)),
                  groups = "grp")
fit$phi_table
#>   group       phi se_log_phi
#> 1     A 0.8353562 0.04924023
#> 2     B 1.1970941 0.04924023
# }
```
