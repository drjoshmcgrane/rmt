# Residual dimensionality of paired comparisons

The residual-PCA analogue for paired comparisons. The fitted model
predicts how often each object should beat each other from their
locations; the object-by-object matrix of departures from that
prediction (on the log-odds scale) is *skew-symmetric*, so its structure
decomposes into rotational planes – Gower's (1977) bimensions – rather
than the ordinary components of a symmetric residual-correlation matrix.
A dominant leading bimension is a coherent “swirl” in the residuals (A
over-beats B, B over-beats C, C over-beats A): a second attribute
steering some contests. A flat spectrum is noise: the single scale
suffices. The leading bimension is judged against a reference built by
simulating unidimensional data from the fitted model with the observed
pair counts (a parametric bootstrap, as in
[`plot_scree`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md));
an observed strength above the reference is structure the
one-dimensional model does not explain. For graded fits the residual
log-odds are taken on the points-proportion scale, whose model mean is
not exactly `plogis(beta_i - beta_j)`; the simulated reference carries
the same construction, so the test stays calibrated (verified mildly
conservative on model-true graded data) rather than anticonservative.
Likewise, when the fit carries within-judge dependence effects
(`order`), the reference is simulated sequentially through each judge's
comparisons WITH the fitted exposure and carry-over coefficients: order
effects push the marginal pair rates around in a structured way, and a
reference without them would read that structure as a second attribute.
The price is power: carry-over and a judge-camp second attribute are
partially confounded (both appear as consistent within-judge deviation),
so with `order` modelled the test is conservative about attributing the
ambiguous share to a second dimension. The reference simulates from the
point estimates without refitting each replicate, so it carries sampling
noise in the responses but not estimation noise in the parameters –
adequate for the screening use here, slightly liberal in tiny designs.

## Usage

``` r
btl_dimensionality(fit, reps = 50L)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- reps:

  Model-simulated replicates for the noise reference.

## Value

A list of class `"rasch_btl_dim"`: `bimensions` (per bimension: strength
and share of residual size; the reference mean, 95th percentile, and the
clears-the-reference flag are reported for the leading bimension and
`NA` for the rest); `coords` (each object's position in the leading
bimension plane, for the residual map); `leading_structured` (whether
bimension 1 clears its reference); `residual_matrix`; and `notes`.

## References

Gower, J. C. (1977). The analysis of asymmetry and orthogonality. In J.
R. Barra et al. (Eds.), *Recent Developments in Statistics* (pp.
109-123). North-Holland.

## Examples

``` r
set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
pr <- t(utils::combn(objs, 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
btl_dimensionality(btl(d, "a", "b", "win"), reps = 20)
#> Paired-comparison residual dimensionality: 3 bimension(s)
#> Leading bimension strength 1.370 (91% of residual; reference 95%: 2.440) -> within noise (one scale suffices)
```
