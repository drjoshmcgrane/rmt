# Unexpected-judgement map for one judge

The judge counterpart of
[`plot_icc`](https://drjoshmcgrane.github.io/rasch/reference/plot_icc.md)'s
kidmap. Each object the judge met is placed by its consensus location
(vertical, strong at top) and by the judge's residual for it
(horizontal; zero, the dashed line, is judged as the scale predicts).
The shaded strip is the expected zone; the rug on the axis marks every
object's location. Objects the judge treated against their standing – a
strong object under-rated (upper left) or a weak object over-rated
(lower right) – are drawn in red and labelled; dot size grows with how
often the judge met the object.

## Usage

``` r
plot_btl_judge_map(fit, judge, min_n = 2L, flag_z = 1.96, ...)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) with
  judges.

- judge:

  The judge to map.

- min_n, flag_z:

  Passed to
  [`judge_surprise`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md).

- ...:

  Unused.

## Value

Called for its plotting side effect; invisibly the `rasch_btl_judge`
object.
