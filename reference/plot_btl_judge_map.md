# Unexpected-judgement map for one judge (pair level)

The judge counterpart of the kidmap, drawn matchup by matchup. Each pair
the judge met is a segment on the consensus location axis, spanning its
two objects, positioned horizontally by how surprising the verdict was:
at zero (the dashed line, inside the shaded band) the stronger object
won as its lead predicts; to the left the judge backed the underdog. A
filled dot marks the object the judge's verdict favoured, hollow the
other – so an upset is a red segment on the left with its filled dot at
the lower end. The rug marks every object's location.

## Usage

``` r
plot_btl_judge_map(fit, judge, min_n = 1L, flag_z = 1.96, ...)
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
  [`judge_pair_surprise`](https://drjoshmcgrane.github.io/rasch/reference/judge_pair_surprise.md).

- ...:

  Unused.

## Value

Called for its plotting side effect; invisibly the
`rasch_btl_judge_pairs` object.
