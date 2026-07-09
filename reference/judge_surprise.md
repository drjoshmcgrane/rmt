# Unexpected judgements of one judge

The paired-comparison counterpart of the kidmap. A judge has no ability
to condition on, so the reference is the consensus object scale (the
pooled locations). For the nominated judge, each object it met is given
a standardised residual oriented to the object – how much more (`z > 0`,
over-rated) or less (`z < 0`, under-rated) that judge favoured it than
its consensus location predicts. A surprise is an object the judge
treated against its standing: a strong object under-rated, or a weak
object over-rated (residual opposite in sign to the location), beyond
`flag_z` and seen at least `min_n` times.

## Usage

``` r
judge_surprise(fit, judge, min_n = 2L, flag_z = 1.96)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) with
  judges.

- judge:

  The judge to profile (a value of the fit's judge column).

- min_n:

  Objects met fewer times are shown but never flagged.

- flag_z:

  Absolute residual at or beyond which a contrary judgement is flagged
  unexpected.

## Value

A list of class `"rasch_btl_judge"`: `objects` (per object met:
location, times met `n`, residual `z`, `surprise` flag and its `type`);
`all_locations` (every object, for orientation); the `judge` and
settings.

## Examples

``` r
set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
pr <- t(utils::combn(objs, 2))
d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
d$judge <- sample(paste0("J", 1:5), nrow(d), TRUE)
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
judge_surprise(btl(d, "a", "b", "win", judge = "judge"), "J1")
#> Judge J1: 40 comparisons over 6 objects
#> No object judged against its consensus standing.
```
