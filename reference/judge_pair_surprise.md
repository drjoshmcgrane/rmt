# Unexpected judgements of one judge, pair by pair

The comparison-level companion of
[`judge_surprise`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md).
Each pair the judge met is oriented to its stronger object (higher
consensus location) and given a standardised residual: `z < 0` means the
stronger object won less than its lead predicts – the judge backed the
underdog. A matchup is an unexpected judgement when `z` falls at or
below `-flag_z` and the pair was seen at least `min_n` times, i.e. the
judge favoured the weaker object further than sampling noise explains.

## Usage

``` r
judge_pair_surprise(fit, judge, min_n = 1L, flag_z = 1.96)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) with
  judges.

- judge:

  The judge to profile.

- min_n:

  Pairs met fewer times are shown but never flagged.

- flag_z:

  Absolute residual at or beyond which an upset is flagged.

## Value

A list of class `"rasch_btl_judge_pairs"`: `pairs` (per matchup: the
stronger and weaker object and their locations, the location `gap`,
times met `n`, residual `z`, the `net_winner`, and the `surprise` flag);
`all_locations`; the `judge` and settings.

## Examples

``` r
set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
pr <- t(utils::combn(objs, 2))
d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
d$judge <- sample(paste0("J", 1:5), nrow(d), TRUE)
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
judge_pair_surprise(btl(d, "a", "b", "win", judge = "judge"), "J1")
#> Judge J1: 40 comparisons over 15 matchups
#> Unexpected judgements (weaker object favoured beyond its lead):
#>   F vs C  (gap 2.12, z = -2.89, upset)
```
