# Fit the Bradley-Terry-Luce model to paired comparisons

Estimates object locations from paired-comparison data by conditional
maximum likelihood. The Bradley-Terry-Luce model is the conditional form
of the dichotomous Rasch model – within an item pair, given one correct
response, the Rasch probability that it was the easier item is exactly
of BTL form – so it belongs to the same measurement family and is
estimated by the same conventions as the rest of the package:
Newton-Raphson on the person-free likelihood, locations identified by
the sum-zero constraint, and Godambe sandwich standard errors, clustered
by judge when a judge column is given (so repeated comparisons by the
same judge need not be independent). Objects that win or lose every
comparison have no finite estimate and are removed with a note, exactly
as extreme persons are set aside in a Rasch calibration; the comparison
graph must remain connected. Beyond that, finite estimates exist only
when the directed win graph is strongly connected (Ford 1957): if some
subset of objects never concedes a point to the rest, the likelihood
pushes the two clusters infinitely far apart, so the fit stops with an
error naming the separated objects rather than presenting the
optimiser's boundary values as measures – remove the separated cluster,
or collect comparisons that cross the divide. Anchors relax this
condition: an anchored object is pinned, so the fit only requires every
free object to be tied to an anchor in both win directions (otherwise
the constrained likelihood still recedes along the unanchored cluster,
and the same error results).

## Usage

``` r
btl(
  data,
  object_a,
  object_b,
  winner = NULL,
  response = NULL,
  margin = NULL,
  judge = NULL,
  count = NULL,
  order = NULL,
  position = FALSE,
  anchors = NULL,
  ties = c("drop", "half", "error"),
  thresholds = c("free", "pc"),
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data:

  A data frame with one comparison per row.

- object_a, object_b:

  Names of the columns holding the two objects compared.

- winner:

  Name of the column holding the winner of each row: its value must
  equal the row's `object_a` or `object_b` entry; `"tie"` or `"draw"`
  marks a tie; anything else (including blanks) is treated as missing
  and dropped with a note. Ignored when `response` is given.

- response:

  Optional name of a column holding a graded preference for `object_a`
  over `object_b` – an ORDERED factor (`factor(..., ordered = TRUE)`,
  levels worst to best for `object_a`) or integer scores `0..m`; a plain
  factor is refused, since its alphabetical level order would silently
  define (and can reverse) the response scale. Fits the
  adjacent-categories ordinal extension of BTL (Tutz 1986; Agresti
  1992): a partial-credit structure on the difference of locations with
  thresholds constrained symmetric, `tau_k = -tau_(m+1-k)`, so the model
  is invariant to presentation order. Two categories reproduce BTL
  exactly; three give the Davidson (1970) ties model.

- margin:

  Optional name of a column holding the extent of the win ("a little",
  "much", ...), as an ordered factor or increasing values; combined with
  `winner` it assembles the graded response without any orientation
  bookkeeping ("B by much" means the same thing whichever column B sits
  in). Winner values matching neither object are ties and form the
  middle category.

- judge:

  Optional name of a judge column; enables the judge fit table and
  clusters the sandwich standard errors by judge.

- count:

  Optional name of a column of replication counts (a row standing for
  several identical comparisons).

- order:

  Optional name of a column giving each judge's judgment sequence
  (timestamps or ranks; requires `judge`). Adds the within-judge
  dependence analysis: an exposure effect (the advantage, in logits, of
  an object the judge has seen before over one they have not) and a
  carry-over effect (the pull of the judge's own earlier verdicts on the
  same object – response dependence in the sense of Marais and Andrich
  2008), estimated jointly with the locations and reported in
  `dependence`. Incompatible with `ties = "half"`.

- position:

  Logical: when `TRUE`, `object_a` is taken as the first-presented
  (left) object of every comparison and a first-position advantage is
  estimated – a single coefficient, in logits, added to every
  comparison's location difference, the pure positional form of the
  Davidson and Beaver (1977) within-pair order-effect device. It is
  reported in `dependence` with `effect = "position"` (every comparison
  is informative, so `n_informative` is the total weighted comparison
  count) and estimated jointly with the locations, alongside the
  exposure and carry-over effects when `order` is also given.
  Identification comes from triangle closure (K \>= 3), so the constant
  oriented covariate is estimable even when each pair has a fixed
  orientation, though weakly. Note that `ties = "half"` duplicates rows
  in the same orientation, so the first position stays well defined.

- anchors:

  Optional named numeric vector for equating: names are object names,
  values are fixed locations in logits. The named objects are held
  exactly at those locations and the remaining objects are estimated
  freely with no sum-zero constraint – the origin and scale come from
  the anchors, exactly as an anchored
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  calibration works. Anchored objects report a standard error of zero
  (their location is a constant, not an estimate). An anchored object
  that is undefeated or winless is an error, not silently removed as a
  free boundary object would be.

- ties:

  How to treat ties in the dichotomous analysis: `"drop"` (default,
  removed with a note), `"half"` (half a win each way, a common
  pragmatic device – flagged in the notes because the halves are not
  independent Bernoulli trials), or `"error"`. With graded responses,
  code ties as a middle category instead.

- thresholds:

  `"free"` (default) estimates every symmetric threshold parameter;
  `"pc"` pools them to the spread (linear) principal component – the
  symmetric case of the principal-component threshold structure, whose
  even skewness component is structurally zero here – so thinly used
  categories borrow strength from every response. Both modes report the
  component decomposition.

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

A list of class `"rasch_btl"`: `objects` (location, se, comparisons,
wins – or the graded `score` – infit and outfit mean squares, fit
residual and its df), `pairs` (per pair: n, observed and expected win
proportions – or mean graded responses – standardised residual,
chi-square component – the pair chi-squares treat comparisons as
independent and are descriptive under judge clustering; the object and
judge fit residuals and the clustered standard errors carry the robust
inference), `judges` (when given: per judge n, infit, outfit, fit
residual, df), `total_chisq`, `total_df`, `total_p`, the object
separation index `osi`, `loglik`, `cl` (the composite-likelihood
information ingredients used by
[`compare_fits`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md):
the Godambe effective parameter count and the independent-unit count),
convergence details, and `notes`. Graded fits add `thresholds` (the
symmetric threshold estimates with standard errors), `m`, and
`categories`. With an `order` column the within-judge `dependence`
effects table carries an `n_informative` count, and `dependence_data`
holds every comparison with its per-comparison exposure and carry-over
covariates (see
[`plot_btl_dependence`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dependence.md)).

## Details

Fit is reported at three levels, mirroring the Rasch diagnostics. Per
object and (when given) per judge: the log-of-mean-square fit residual
of Andrich and Marais (2019, ch. 23) over their comparisons, with
apportioned degrees of freedom – an erratic judge or an object of
inconsistent quality shows exactly as an erratic person or misfitting
item does. Per pair: the classical goodness-of-fit table comparing
observed and expected win proportions, with the total chi-square on
(pairs used) minus (free location parameters) degrees of freedom. The
object separation index is the analogue of the PSI: the proportion of
observed location variance not due to error. Anchored objects enter it
with their fixed locations and zero error – they are separated with
certainty by construction.

## References

Bradley, R. A. and Terry, M. E. (1952). Rank analysis of incomplete
block designs: I. The method of paired comparisons. Biometrika, 39,
324-345. Luce, R. D. (1959). Individual Choice Behavior. Wiley. Andrich,
D. (1978). Relationships between the Thurstone and Rasch approaches to
item scaling. Applied Psychological Measurement, 2, 451-462.

Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response.
Journal of Mathematical Psychology, 30(3), 306-316. Agresti, A. (1992).
Analysis of ordinal paired comparison data. Journal of the Royal
Statistical Society C, 41(2), 287-297. Davidson, R. R. (1970). On
extending the Bradley-Terry model to accommodate ties in paired
comparison experiments. Journal of the American Statistical Association,
65(329), 317-328. Ford, L. R. (1957). Solution of a ranking problem from
binary comparisons. American Mathematical Monthly, 64(8), 28-33.

Davidson, R. R., & Beaver, R. J. (1977). On extending the Bradley-Terry
model to incorporate within-pair order effects. Biometrics, 33(4),
693-702.

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pairs <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pairs[, 1], each = 30),
                b = rep(pairs[, 2], each = 30))
p <- plogis(beta[d$a] - beta[d$b])
d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
btl(d, object_a = "a", object_b = "b", winner = "win")
#> Bradley-Terry-Luce analysis: 4 objects, 180 comparisons
#> Conditional ML: converged in 6 iterations; sandwich SEs
#> Object separation index 0.963; pairwise chi-square 1.07 on 3 df, p = 0.783
#>  object location    se comparisons wins fit_resid
#>       A   -1.238 0.214          90   16    -0.111
#>       B   -0.354 0.186          90   36     0.526
#>       C    0.448 0.180          90   56    -0.056
#>       D    1.144 0.209          90   72     0.037
```
