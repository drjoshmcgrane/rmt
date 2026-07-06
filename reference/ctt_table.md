# Traditional (classical test theory) statistics

The classical companion table conventionally reported alongside a Rasch
analysis (Andrich and Marais 2019, chs. 3-5), on complete cases only:
per item the facility (mean score over maximum), the item-total and
corrected item-rest correlations, the discrimination index DI = PRU -
PRL (mean proportion-of-maximum in the upper third of total scores minus
the lower third), and alpha if the item is deleted; plus coefficient
alpha, the raw-score mean, SD, and the classical standard error of
measurement \\s\sqrt{1 - \alpha}\\, which unlike the Rasch SE is one
value for all persons.

## Usage

``` r
ctt_table(fit)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

## Value

A list of class `"rmt_ctt"`: the per-item `table` (`item`, `max`, `n`,
`facility`, `item_total`, `item_rest`, `di`, `alpha_drop`), and the
scalars `alpha`, `n` (complete cases), `mean`, `sd`, and `sem`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
colnames(X) <- paste0("I", 1:8)
ctt_table(rasch(X))
#> Traditional statistics (available cases; item n 400-400; 400 complete)
#> Raw score mean 4.08, SD 1.71 (complete responders); alpha 0.529; SEM 1.17
#>  item max   n facility item_total item_rest    di alpha_drop
#>    I1   1 400    0.848      0.370     0.169 0.298      0.519
#>    I2   1 400    0.787      0.443     0.221 0.373      0.504
#>    I3   1 400    0.672      0.479     0.226 0.484      0.504
#>    I4   1 400    0.575      0.526     0.267 0.589      0.488
#>    I5   1 400    0.432      0.579     0.333 0.663      0.460
#>    I6   1 400    0.338      0.495     0.243 0.497      0.497
#>    I7   1 400    0.250      0.533     0.313 0.527      0.472
#>    I8   1 400    0.175      0.407     0.198 0.316      0.511
```
