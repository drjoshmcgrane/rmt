# Differential item functioning by two-way residual ANOVA

For each item and each person factor separately, analyses the
standardised residuals by factor group and trait class interval. The
group main effect indicates uniform DIF and the group-by-interval
interaction indicates non-uniform DIF. Probabilities are adjusted across
items within each factor by the Benjamini-Hochberg false-discovery-rate
procedure (or any [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)
method). With several factors, consider
[`dif_anova_factorial`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova_factorial.md),
which models them jointly.

## Usage

``` r
dif_anova(
  fit,
  factors = NULL,
  n_groups = NULL,
  p_adjust = "BH",
  alpha = 0.05,
  id = NULL,
  within = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- factors:

  A vector (one factor), a data frame of person factors, or a character
  vector naming factor columns already nominated in the fit (via
  `rasch(..., factors = )`). Defaults to every factor stored in the fit.

- n_groups:

  Number of trait class intervals. By default set per factor from the
  smallest group so every interval-by-group cell keeps about 30 expected
  responses (between 2 and 10 intervals) – independently of the interval
  count of the overall fit, whose rule guards intervals, not cells. The
  counts used are returned as the `n_groups` attribute.

- p_adjust:

  Multiplicity adjustment method passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html); `"BH"` (default)
  controls the false discovery rate, `"holm"` or `"bonferroni"` the
  familywise error rate.

- alpha:

  Significance level applied to the adjusted probabilities.

- id:

  Optional person identifier (a vector, or the name of a nominated
  factor) for stacked repeated-measures designs. A factor whose levels
  vary within a person is treated as within-subject and tested with a
  person-clustered sandwich, so the within-person dependence the
  between-subjects F ignores is respected. Defaults to the fit's own
  person identifier, so a stacked design is handled automatically.

- within:

  Optional names of factors to treat as within-subject; auto-detected
  from `id` when not given.

## Value

A data frame per item and factor with the full two-way table (Andrich
and Marais 2019, ch. 16): the group main effect (uniform DIF), the
group-by-interval interaction (non-uniform DIF), and the class-interval
main effect, each with its F statistic; raw and adjusted probabilities
and flags for the two DIF terms; and partial eta-squared effect sizes
(`eta2_uniform`, `eta2_nonuniform`), the proportion of
residual-plus-effect variance the term accounts for; and a `within`
indicator per factor. Factors treated as within-subject are named in the
`within` attribute. For DIF magnitude on the logit scale, where
practical significance is judged, see
[`dif_size`](https://drjoshmcgrane.github.io/rmt/reference/dif_size.md).

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
dif_anova(rasch(X), factors = data.frame(group = g))
#>   factor item  F_uniform    p_uniform eta2_uniform F_nonuniform p_nonuniform
#> 1  group   I1  0.1570992 6.919910e-01 0.0002779744    1.0353885   0.40113597
#> 2  group   I2  0.4333293 5.106279e-01 0.0007663668    0.1601699   0.98695468
#> 3  group   I3 27.9738443 1.760845e-07 0.0471755113    1.8762745   0.08278624
#> 4  group   I4  0.4414699 5.066849e-01 0.0007807526    0.4303976   0.85876639
#> 5  group   I5  8.1163049 4.546746e-03 0.0141617065    2.4552343   0.02366131
#> 6  group   I6  1.4228990 2.334275e-01 0.0025120789    0.2161082   0.97167194
#> 7  group   I7  2.4010934 1.218104e-01 0.0042317392    1.6026234   0.14403948
#> 8  group   I8  0.1780050 6.732539e-01 0.0003149538    1.2077998   0.30043851
#>   eta2_nonuniform   F_class    p_class within p_uniform_adj p_nonuniform_adj
#> 1     0.010875695 0.5516798 0.76877969  FALSE  6.919910e-01        0.6418175
#> 2     0.001698031 1.5486895 0.16004610  FALSE  6.808372e-01        0.9869547
#> 3     0.019535788 1.3859574 0.21800630  FALSE  1.408676e-06        0.3311450
#> 4     0.004549799 0.8863650 0.50446152  FALSE  6.808372e-01        0.9869547
#> 5     0.025410743 1.1171114 0.35078915  FALSE  1.818699e-02        0.1892905
#> 6     0.002289699 0.8543115 0.52841764  FALSE  4.668550e-01        0.9869547
#> 7     0.016734209 2.8661923 0.00928187  FALSE  3.248277e-01        0.3841053
#> 8     0.012663765 0.5068670 0.80333411  FALSE  6.919910e-01        0.6008770
#>   uniform_DIF nonuniform_DIF
#> 1       FALSE          FALSE
#> 2       FALSE          FALSE
#> 3        TRUE          FALSE
#> 4       FALSE          FALSE
#> 5        TRUE          FALSE
#> 6       FALSE          FALSE
#> 7       FALSE          FALSE
#> 8       FALSE          FALSE
```
