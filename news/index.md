# Changelog

## rasch 1.13.3

Sixteenth review round: two identification edge cases beneath the 1.13.2
checks, and one policy refinement. Development branch only; the CRAN
submission of 1.11.7 is untouched.

- Connectivity now counts only INFORMATIVE response pairs. A pair whose
  every observed total is 0 or the maximum has a single feasible
  conditional allocation and carries no information, so one respondent
  scoring (0, 0) across two otherwise disjoint item blocks no longer
  “connects” them (every intermediate total has at least two
  allocations, so any response off the extreme-total corners is a real
  link). Applied to the item-pair graph of
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)/[`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  and to the MFRM co-observation components.
- As the review prescribed, the graph test is now backed by a
  rank/conditioning check of the projected information at the solution
  in the shared solver. This also catches what no graph can see: a block
  bridged by a SINGLE informative response is perfect separation in the
  conditional pair logit – the pair estimate runs to the boundary and
  the information vanishes there – and previously returned converged
  fits with zero standard errors.
- The Ford (1957) check now respects anchors. Strong connectivity is the
  existence condition for the FREE model only; an anchored fit requires
  each free object to be tied to an anchor in both win directions
  (reachability to and from the anchor set over the points digraph – the
  constrained recession-direction condition). Two anchored components
  joined by a one-directional edge between their anchors is now
  correctly accepted, while removing one anchor still errors.
- EFRM: the practical weak-identification cutoff (SE of log phi above
  5.  is now a loud warning plus a note with the estimates RETAINED,
      supporting reproducible sensitivity analyses; the error is
      reserved for structural non-identification (a flat direction of
      the joint information), where no estimate exists to retain.

## rasch 1.13.2

Fifteenth review round: the 1.13.1 identification guards were the right
idea executed too weakly (or, in one case, too strongly), and this
release replaces heuristics with information-based checks throughout.
Development branch only; the CRAN submission of 1.11.7 is untouched.

- EFRM: the threshold-spread heuristic is gone. Group units are now
  checked on the joint information itself – a flat direction loading on
  a unit, or a unit whose analytic SE exceeds 5 log-units (uncertain
  beyond a factor of ~150), is refused with an error naming the group,
  since every common-unit quantity would silently depend on it. Weakly
  identified units with real spread are KEPT with their honest large
  SEs; the previous heuristic wrongly replaced some of them with NA.
- MFRM: full column rank of the structural design is necessary but not
  sufficient. The fit now also checks that no between-block shift of
  person-disjoint response blocks is expressible by the facet map (a
  column-space intersection test on the co-observation components):
  designs such as two person groups answering disjoint item pairs – or
  one rater per person – previously returned converged fits whose
  between-block contrasts were ridge artefacts, and now error with the
  blocks named. A test fixture with exactly this nested-rater flaw was
  corrected accordingly.
- btl_efrm stage 1: the per-set (locations, panel-ratio) information is
  now checked for rank at the solution; panels observing disjoint object
  pairs leave the ratio underdetermined, and such sets are screened out
  of the panel-unit reconciliation (refit at the reconciled units) or,
  when no ratio information remains anywhere, the fit stops.
- btl_efrm stage 2: rank failures are now CLASSIFIED by where the
  information fails. A flat direction confined to a set’s log-unit – its
  within-set locations are indistinguishable, so the unit has nothing to
  scale – refits with that unit fixed at the conventional 1 and reports
  alpha as NA with a note; the set’s origin kappa and its objects’
  placements remain identified and keep valid SEs (this is the
  GermanParties2009 near-even-set pattern, which a blanket refusal would
  have broken). A flat direction on an origin means the set cannot be
  placed at all and is an error. Bootstrap replicates treat an expected
  NA (the declared unit) as normal and any unexpected rank failure as a
  failed replicate.
- BTL: a Ford (1957) violation is now an error before estimation, not a
  post-fit warning – a separated cluster has no finite maximum
  likelihood locations, and the previous warning still presented the
  optimiser’s boundary values with converged = TRUE.
- EFRM input: data-frame `factors` columns are removed from the item
  matrix (a numeric factor column could previously become an implicit
  `(rest)` item).
- MFRM input: rows dropped for missing person/item/facet identifiers now
  drop from `factors` too (both the column-name and data-frame forms),
  instead of failing with an internal length error.
- Test honesty: the wide-MFRM factor-equivalence test compared a
  nonexistent `$table` component (NULL == NULL, vacuously true); it now
  compares the `$summary` F and p columns, on a properly linked design –
  the equivalence itself was verified to hold.

## rasch 1.13.1

Fourteenth review round: three defects in the 1.13.0 additions and the
review’s nine-item identification backlog. Development branch only; the
CRAN submission of 1.11.7 is untouched.

Defects in the new code:

- EFRM crossed-frame labels: the solver orders cells as sorted strings
  while the frame construction used factor-crossing order, so with
  factor labels that sort differently the per-cell units were assigned
  to the wrong cells (a pure-noise 2x2 gave a spurious region effect at
  p = 2e-5). Cells are now matched by label throughout.
- MFRM
  [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  no longer lets facet severity leak into a between-group magnitude: the
  joint split now pools over the COMMON observed cells with common
  precision weights, so severity differences cancel by construction. A
  design that previously manufactured a -1.75-logit phantom effect (z =
  -6.5) now reports -0.11 (ns) while a planted 0.7 effect is still
  recovered.
- `phi_factorial` was a diagonally-weighted fit with coefficient-wise z
  tests; it is now a generalised least-squares decomposition on the
  JOINT covariance of the cell log-units (bootstrap draws when
  available, otherwise the solver’s centred covariance, inverted
  spectrally along its identified directions, intercept excluded as
  inestimable under the centring), with multi-degree-of-freedom Wald
  tests per term in `phi_factorial_tests`. A planted 1.5x region unit
  ratio is recovered at 1.52 with p = 2e-6 and null factors stay null.

Identification and input-handling backlog:

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) now
  checks Ford’s (1957) existence condition: when the directed win graph
  is not strongly connected the maximum-likelihood locations diverge,
  and the fit warns and notes the separated objects instead of
  presenting the optimiser’s boundary values as measures.
- MFRM: wide input now carries `factors=` through the internal melt
  (character columns are replicated like facets, data-frame factors
  row-wise); duplicate person-by-item-by-facet responses are an error
  (keeping the first made results depend on row order); and a structural
  rank check on the facet design errors when facet levels are confounded
  with the items instead of returning a valid-looking but unidentified
  decomposition.
- EFRM: a group whose frames show no threshold spread beyond estimation
  noise has nothing for its unit to scale – previously the optimiser
  returned phi near 1 with a healthy-looking SE regardless of the true
  unit; such units are now reported NA with a note. The guard is
  signal-to-noise based (spread under 1.5x the pooled threshold SE in
  every set the group answers) and leaves even modest real spreads
  untouched.
- EFRM: under `se_method = "bootstrap"` the frame-level `se_log_rho` now
  comes from the joint replicate draws of log(alpha) + log(phi),
  capturing cross-stage dependence; the hybrid fallback still combines
  the stagewise errors as if uncorrelated, as documented.
- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md):
  stage 2 checks the conditioning of its observed information; when the
  cross-set comparisons barely identify the set units (reciprocal
  condition number below 1e-10), the conditional SEs for alpha, kappa
  and the linked values are withheld with a note rather than reported
  from a near-singular inverse.
- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  frame table: `rho` is now phi / alpha as documented (the within-set
  logit is phi (beta_a - beta_b) with beta = (v - kappa) / alpha), not
  phi \* alpha.
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  validates `item_sets` (named list, no overlapping or unknown items);
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  removes data-frame factor columns from the response matrix even when
  `items` is not given, so a numeric factor can no longer be silently
  treated as an item.

## rasch 1.13.0

The three capability extensions from the review wishlist. Development
branch only; the CRAN submission of 1.11.7 is untouched.

- MFRM DIF pools to the underlying items by default:
  [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  on an MFRM fit now answers “does item A show DIF” by pooling residuals
  over each item’s facet cells (`pool_facets = FALSE` keeps the
  per-virtual-cell tests). With few items, the compensating artificial
  DIF the plant creates on the remaining items is visible – as it should
  be – with the planted item carrying the dominant F.
- [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  works on underlying MFRM items: every facet cell of the item is split
  by the groups in one joint unstructured refit of the virtual matrix
  (the facet decomposition is not reimposed, which the note states), and
  per-level locations are precision-weighted means over the cells with
  the full covariance carried. A planted 0.7-logit sex effect is
  recovered at 0.82 (SE 0.17), and `dif_anova(mf, sizes = TRUE)` now
  yields magnitude rows.
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  accepts several frame-defining factors: the frames are their crossed
  cells, per-cell units appear in `phi_table`, and `phi_factorial`
  reports a factorial decomposition of the cell units (sum-coded mains,
  plus the interaction when every cell is observed), weighted by the
  cells’ unit precisions with the independence approximation stated. All
  frame-defining factors (cell and components) are excluded from DIF
  testing; misspelled group and factor column names now error instead of
  silently falling through.

## rasch 1.12.3

Thirteenth review round: the BTL, EFRM, and MFRM analogues of the DIF
and identification work. Development branch only; the CRAN submission of
1.11.7 is untouched.

- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  runs on the same order-invariant machinery as the person DIF ANOVA:
  between-judge terms by Type II sums of squares on band-centred judge
  margins, band-crossing terms through orthonormal contrasts with the
  Greenhouse-Geisser correction. Verified exactly order-invariant across
  correlated judge factors; the tuned planted detections are unchanged.
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  rejects factors that vary within a judge (the judge-level analysis
  would silently take whichever row came first).
- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  validates within-set object-graph connectivity (relative locations
  inside a disconnected set are unidentified, exactly as in
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md))
  and alpha identification (cross-set comparisons must touch at least
  two objects of every non-reference set; one object identifies only the
  origin).
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  group linkage requires at least two SHARED ITEMS of a common set
  between groups: sharing a set label with disjoint item subsets left
  the unit ratio unidentified but was accepted.
- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  works on EFRM fits: the frame group is excluded (it is the frame
  structure – every frame has its own virtual items – so it has no
  within-item contrast) with a note, and other person factors are tested
  per virtual item; nominating only the frame group is an informative
  error.
- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  gains `factors=` (person factors, constant within person, validated –
  a facet passed as a person factor errors with a pointer to
  `interaction=`), so
  [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  works on MFRM fits directly; unavailable DIF magnitudes are reported
  in `notes` instead of silently returning nothing.
- The frames table’s `se_log_rho` independence approximation (log alpha
  and log phi estimated in different stages; cross-stage covariance
  unavailable) is stated in the documentation.

## rasch 1.12.2

Twelfth review round: the incomplete-panel edges of the DIF engine.
Development branch only; the CRAN submission of 1.11.7 is untouched.

- Between-person margins are centred within each within-cell BY CLASS
  INTERVAL (falling back to the cell mean where a combination is empty):
  cell-only centring removed common occasion effects but not
  occasion-by-trait structure, which with differential missingness
  masqueraded as non-uniform group DIF at F = 214.7; it now reads F =
  0.01 with the occasion terms intact.
- A between level that loses every complete within panel (for example a
  group observed only at baseline) no longer crashes lm() with a
  one-level factor: its within-stratum interactions are reported NA with
  an explanatory note, and every other test proceeds.
- Tukey candidates are filtered to between-only terms, so a significant
  multilevel within term can no longer reach TukeyHSD and crash on its
  empty result.
- Transparency: the count of person panels dropped from the
  within-person tests is returned in `notes`, alongside the
  non-estimable-term note.

## rasch 1.12.1

Eleventh review round on the rebuilt DIF engine: five findings, all
verified and closed. Development branch only; the CRAN submission of
1.11.7 is untouched.

- Multi-within contrast alignment.
  [`interaction()`](https://rdrr.io/r/base/interaction.html) orders
  cells with the FIRST factor fastest while the Kronecker contrast
  matrices assume the LAST fastest; with two within factors a pure w1
  effect therefore loaded onto w1:w2 (F 7.4 on the interaction, 3.2 on
  w1). The within cells are now built in Kronecker order; the same
  effect reads F = 16.9 on w1 and 0.35 on the interaction.
- Comparable between-person margins. Person means for the between tests
  are computed on within-cell-centred residuals: with a common occasion
  effect and one group missing 80% of an occasion, raw means reported
  group DIF at F = 37.6 (p = 2e-9); centred margins give F = 0.29 while
  the occasion effect stays on the occasion term.
- Tukey follow-ups in mixed designs run on the person-level
  between-factor aov (one row per person); ordinary Tukey on repeated
  within cells is no longer offered.
- An NA in a between factor no longer flips it to within-subject
  (detection counts distinct non-missing values).
- The Greenhouse-Geisser correction is fully reported: `terms` carries
  `df_denom` and `gg_epsilon`, the p-value is exactly pf(F, eps*df,
  eps*df_denom), and the stale “sums of squares are sequential” sentence
  is replaced by the Type II description.

## rasch 1.12.0

The multi-factor DIF ANOVA engine rebuilt for statistical validity
(tenth review round; three P1 findings, every one reproduced before
fixing). This release is on the development branch only – it does NOT
supersede the CRAN submission of 1.11.7.

- Order-invariant tests. Between-person terms now use Type II sums of
  squares on person-level residual means: every term is adjusted for
  every term not containing it, with the class interval always among
  them. The old sequential (Type I) tests with factors entered first let
  entry order decide which of two correlated factors flagged – in the
  reproduction, reversing the order LOST a planted 0.8-logit DIF
  entirely – and let a group factor absorb pure trait variance (reported
  p = 2e-25 where the adjusted test gives p = 0.83). Verified exactly
  order-invariant.
- Persons are the units whenever ids repeat. Residuals aggregate to one
  mean per person (per within cell) before any test, and the automatic
  class-interval rule counts persons, not rows: exactly duplicating
  every person previously changed F = 4.06 (unflagged) to F = 8.21
  (flagged); it now changes nothing (verified to 2e-13).
- Within-person terms are tested on person-by-cell means through
  orthonormal contrasts with the Greenhouse-Geisser epsilon (Maxwell &
  Delaney 2004), replacing raw multi-stratum aov: a nonspherical
  four-level null rejected ~9% at nominal 5% before (and could crash on
  tied person means – also fixed, in the interval allocator); it now
  runs at 2.8%, with the spherical null at 6.1% and planted
  within-occasion DIF still detected. Persons missing a within cell are
  dropped from the within tests explicitly rather than projected
  murkily. The two-level case still agrees with its paired-t gold
  standard.
- Input honesty: unknown `within` names error; declaring a factor
  within-subject with no repeated ids errors; and a factor that VARIES
  within persons can no longer be forced between-subjects (the old
  row-level treatment pseudo-replicated).
- The BH family choice (per term, across items) is now documented
  explicitly.

## rasch 1.11.7

Ninth review round: one finding.

- [`print()`](https://rdrr.io/r/base/print.html) distinguishes the three
  saved-fit dependence schemas by their `df` field: current fits label
  the statistic `t`; 1.11.4 transitional fits (`z` name, `df` present,
  t-based p) print as `t`; older fits (`z`, no `df`, normal-reference p)
  KEEP the `z` label – relabelling them `t` misrepresented how their
  p-values were computed (a legacy z = 3.17 with normal p = .002 is not
  the same claim as t = 3.17 with p = .016 at five clusters). The
  regression test now reconstructs both legacy schemas faithfully.

## rasch 1.11.6

Eighth review round: label and migration polish.

- [`print()`](https://rdrr.io/r/base/print.html) reads the dependence
  statistic from either the current `t` column or the pre-1.11.5 `z`
  column, so fits saved before the rename print completely; the printed
  label and the rank-deficiency note both say `t`, matching the t(G - 1)
  reference. No duplicate `z` alias column is added: the package has
  never had a released version with the old name, the rename is
  documented here, and old saved fits are read transparently.

## rasch 1.11.5

Seventh review round: three completions of round-six fixes, one label.

- `margin=` gets the same explicit-ordering rule as `response=`: an
  ORDERED factor (smallest to largest margin) or a numeric magnitude;
  plain factors AND character columns are refused, since alphabetical
  order could silently reverse which margin counts as the big win.
- Invalid graded numeric responses (“abc”, Inf) error through the shared
  guard instead of becoming missing and being dropped – the fit no
  longer succeeds on silently reduced data.
- The simulator’s item resolver rejects fractional indices (4.9 no
  longer truncates to item 4), an empty `second_dim$items`, and a
  vector-valued `rho`.
- The clustered dependence/position statistic is labelled `t`, matching
  its t(G - 1) reference (at infinite df it is the familiar z). The
  reviewer’s independent 400-run null simulation of the cluster-t
  correction found 5.25% rejection at nominal 5%.

## rasch 1.11.4

Sixth review round: four findings, all verified and closed.

- The integer-score guard reads factors through their LABELS:
  `as.numeric(factor)` returns level codes, so a factor score of “1.9”
  slipped past as the integer 3 and was then truncated. Non-numeric and
  non-finite scores now also error at every entry point instead of
  silently becoming missing.
- Graded
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  requires an ORDERED response factor: a plain factor’s alphabetical
  level order silently defined – and could reverse – the response scale
  (worse \< same \< better read as better \< same \< worse). Integer
  scores 0..m remain the alternative.
- Clustered dependence and position tests use a t reference with G - 1
  degrees of freedom (the standard few-cluster correction) instead of
  normal theory; the table carries the df. Five judges now give honestly
  wide p-values.
- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md)
  validates the second-dimension specification: `rho` outside \[-1, 1\]
  and unknown item names/indices are errors (rho = 1.2 used to produce
  an all-missing second dimension that looked like a successful
  simulation), and the item-set resolver now range-checks every
  nominated item list centrally.

## rasch 1.11.3

Fifth review round: seven findings, all verified and closed.

- Every public estimator now rejects fractional scores before integer
  coercion could truncate them:
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md),
  [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md),
  and
  [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  truncated 1.9 to 1 silently, and graded
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  rounded it; only
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  rejected. One shared guard serves all entry points.
- [`item_moments()`](https://drjoshmcgrane.github.io/rasch/reference/item_moments.md)
  is overflow-stable (log-sum-exp): a 31-category item, or wide
  categories with a large discrimination, overflowed the uncentred
  exponent and turned every
  [`person_wle()`](https://drjoshmcgrane.github.io/rasch/reference/person_wle.md)
  location NA.
- When the judge-clustered covariance is rank-deficient (clusters \<=
  parameters) the OSI is withheld as NA – understated SEs would
  overstate the separation reliability – and the note now spells out how
  to read the marginal SEs and dependence tests. The SEs themselves stay
  reported: they are consistent estimates whose understatement is
  disclosed, which is standard few-cluster practice.
- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  rejects a judge assigned to more than one panel (a panel is a judge
  attribute; the judge bootstrap would otherwise silently reclassify
  rows).
- Undefined diagnostics are NA, not numbers: Cronbach’s alpha under zero
  total-score variance (was -Inf), and the omnibus item-trait test when
  no item is testable (was chi-square 0 on 0 df with p = 1).
- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md)’s
  secondary trait keeps the requested mean and sd (combining two mean-mu
  components had shifted the mean by the factor rho + sqrt(1 - rho^2))
  and is returned in the truth metadata.
- The MLE score-table root drops the common discrimination, which
  cancels from the ML equation exactly as it does from Warm’s (ordinary
  unit-discrimination fits were unaffected).

## rasch 1.11.2

Edge cases from the fourth review round, verified and closed:

- The BTL dimensionality reference now simulates at the PAIR level:
  count weights are summed per unordered pair and rounded, so half-tie
  rows (weight 0.5 in each direction) recombine into whole comparisons.
  Row-level simulation broke on fractional weights (`as.integer(0.5)` is
  a zero binomial size – every reference draw degenerated to 0). The
  per-replicate leading strengths are exposed as `$reference$draws`.
- `btl_efrm(se_method = "judge_bootstrap")` refuses a panel with a
  single judge (resampling one judge returns the same data every time,
  so its SEs would be a spurious zero) and notes panels with fewer than
  five judges as rough.
- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  notes when there are no more judge clusters than parameters: the
  clustered covariance is then rank-deficient by construction and no
  scalar (CR1) correction repairs it – the SEs likely understate.

## rasch 1.11.1

The follow-up review’s residual findings, all verified and closed:

- The BTL dimensionality reference now simulates binomial (multinomial)
  sums for count-weighted rows instead of one weighted Bernoulli, whose
  w^2 p(1-p) variance overdispersed the reference and made the test
  conservative on aggregated data; aggregated and expanded data now give
  the same reference (verified to \< 1%).
- The pairwise chi-square degrees of freedom subtract every estimated
  parameter (position and dependence covariates were omitted), and a
  design with nothing testable reports NA instead of df = 1.
- Judge-clustered covariance carries the standard CR1 small-sample
  factor G/(G-1), on top of the single-judge error and few-cluster note.
- `btl_efrm(se_method = "judge_bootstrap")`: judges resampled with
  replacement within panels and the whole pipeline refitted, so the
  errors carry within-judge dependence as well as stage-one uncertainty
  – closing the documented gap in the parametric bootstrap. On
  model-true data it agrees with the parametric bootstrap; estimates are
  unchanged by construction.
- The next-pair priority documentation is internally consistent (ranking
  heuristic; the Sherman-Morrison update is exact for information
  matrices, a scoring device on the sandwich).

## rasch 1.11.0

Statistical validity, from a second external review. Every finding was
re-derived or reproduced by simulation before fixing; the fixes are
verified by new calibration tests (test-statistical-validity.R).

- Recentred covariance for mixed max scores.
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  recentres thresholds so the mean item location is zero, but the
  covariance stayed under the sum-of-thresholds constraint – identical
  only when every item has the same maximum score. With mixed
  1/2/3-threshold items, per-item empSD/SE ran 0.90-1.21 by item type;
  the linear transform `(I - 1a') cov (I - a1')` restores 0.96-1.09.
  Affects threshold and item SEs, separation, and equating weights on
  mixed designs.
- Warm WLE discrimination cancellation. The common-discrimination
  weighted score is `disc[(R-E) + mu3/(2V)]` – the discrimination
  cancels; a stray `disc` factor in the correction biased WLEs by up to
  0.26 logits at disc = 0.5 (exact at disc = 1, so ordinary analyses
  were unaffected; the vector-disc EFRM path was already correct).
- Honest chi-square degrees of freedom. An item whose responders occupy
  fewer than two class intervals has no estimable item-trait
  interaction: chi-square, df, and p are now NA (df was manufactured as
  1 with p = 1), and the total test sums over testable items only. The
  `adjust_N` proportional-scaling semantics are now documented
  explicitly.
- Covariance-correct equating drift tests, both families. The shift c0
  is estimated from the same common items the tests examine, and
  locations within a calibration are correlated through the
  identification constraint: drift denominators now use
  `[(I - 1u') Sigma (I - u1')]_ii` with Sigma the sum of the two
  calibrations’ location covariances (a bank contributes diag(se^2)).
  Null calibration: 4.6% rejection at nominal 5%, t-SD 0.97.
- Judge-clustered SEs refuse a single judge (they collapse to ~1e-16)
  and note fewer than 10 clusters as likely understating.
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  now tests at the judge level: residuals aggregate to one weighted mean
  per judge (per opponent band) and enter a split-plot ANOVA with the
  judge as the error unit. The comparison-level ANOVA pseudo-replicated
  – 6 of 10 null datasets with judge heterogeneity and arbitrary groups
  falsely flagged uniform DIF; the judge-level design flags 0 of 10
  while still detecting a planted 1-logit effect. Band breaks moved to
  count-weighted quantiles, making aggregated and expanded data exactly
  equivalent. Power now scales with judges per group – the honest unit.
- Non-convergence now warns loudly at fit time in
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  and [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  (it was only visible in print);
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)’s
  convergence flag is scale-free (gradient per comparison).
- Documentation: the test-of-fit suite’s calibration status is stated
  plainly (approximate, convention-faithful diagnostics; simulation
  shows near-but-not-nominal rejection), with parametric-bootstrap
  references via
  [`sim_replicate()`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md)
  recommended where exact calibration matters; the BTL pair chi-squares
  are marked descriptive under judge clustering; the next-pair priority
  is documented as a ranking heuristic; the BTL dimensionality
  reference’s no-refit design is stated.

## rasch 1.10.3

Input honesty, from an external code review (five findings, all
confirmed and fixed):

- Misspelled `id`, `factors`, or `items` column names in
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  are now errors naming the missing columns; they were silently ignored
  (one misspelled factor silently dropped ALL factors), producing
  valid-looking analyses of the wrong data. Numeric item indices still
  work.
- Fractional scores (e.g. 1.9) are now an error naming the columns;
  [`as.integer()`](https://rdrr.io/r/base/integer.html) used to truncate
  them silently, altering response data. Integer-valued doubles (“2.0”)
  still pass.
- MFRM rows with a missing person, item, or facet identifier are dropped
  with a note; [`paste()`](https://rdrr.io/r/base/paste.html) used to
  coerce `NA` to the literal string “NA”, silently pooling unrelated
  rows under a phantom level.
- [`equate_tests()`](https://drjoshmcgrane.github.io/rasch/reference/equate_tests.md)
  excludes common items whose location or SE is unavailable (e.g. weakly
  determined items with honest NA SEs) from the precision-weighted shift
  and drift tests, notes them, and keeps their rows with NA test columns
  – one such item used to turn every statistic NA.
  [`plot_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_equate.md)
  handles the excluded rows.
- [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  escapes all data-derived text (title, notes, item and facet names) –
  crafted column names could inject markup into reports; table cells
  were already escaped.

## rasch 1.10.2

- App: hover identification is now consistent across the app, through
  one shared mechanism. Beyond the person-fit plot and item fit map,
  hovering identifies the point (or cell) on: the Guttman scalogram
  (person, item, score), both residual-correlation heatmaps (item pair
  and its Q3/Q3\*), the common-item equating plot and the
  tailored-analysis comparison (item, reference and current locations),
  and – for paired comparisons – the unexpected-judgements map (matchup,
  surprise, location), the group characteristic curves’ observed means
  (object, panel, mean), and the equating plot (object, both
  calibrations). Every tooltip is resolved server-side against exactly
  the rows its plot draws. Plots whose points are already labelled
  (kidmap, Wright map, caterpillars) and curve or distribution displays
  are deliberately left without hover.

## rasch 1.10.1

- App: hovering a point on the person-fit plot or the item fit map now
  shows who it is – a floating tooltip with the person ID (or item
  name), location, and fit residual, resolved server-side against
  exactly the subset each plot draws. The navbar wordmark is
  left-justified.

## rasch 1.10.0

- Model comparison, calibrated for composite likelihood.
  [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  gains `cl_aic` and `cl_bic` (Varin & Vidoni 2005; Gao & Song 2010):
  `-2 cl` penalised by the effective parameter count `tr(H^-1 J)` from
  the Godambe matrices – the same quantity whose eigenvalues calibrate
  [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
  – which absorbs the pairwise over-counting a nominal AIC/BIC would
  ignore (empirically 3.5-5x the nominal count). Verified by
  plant-and-detect: the criteria select RSM on rating-structured data,
  PCM under varying threshold spreads, and a planted BTL position effect
  – and reject the richer model when the effect is absent. `n` for
  CL-BIC counts independent units: persons contributing an informative
  pair, or judges (count-weighted comparisons when unclustered).
  MFRM/EFRM fits report NA with the reason noted.
- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  now accepts `btl` fits (all-BTL comparisons: free vs
  principal-component thresholds, with/without a position effect or
  within-judge dependence), reporting judges/objects/comparisons and
  OSI. [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  fits carry their composite-likelihood ingredients in `$cl`; anchored
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  and
  [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  now return their Godambe matrices. The app’s Compare page shows the
  criteria and accepts paired-comparison fits.
- Cross-validation against independent implementations, in the test
  suite whenever the packages are installed (new in Suggests: eRm, sirt,
  psychotools): dichotomous, PCM, and RSM parameters against eRm’s full
  CML (agreement to sampling precision; pairwise SEs at or just above
  CML’s, the documented efficiency price), the same pairwise family
  against sirt::rasch.pairwise (near-identity), and btl() against
  psychotools::btmodel to machine precision (same likelihood).

## rasch 1.9.3

- Missing-data identification, verified and enforced. A booklet-design
  probe (zero complete cases, linked forms, MCAR on top) confirmed the
  pairwise conditional estimator recovers planted parameters with
  nominal coverage – and exposed two silent failure modes at the
  design’s edges, both now honest:
- Disconnected designs error instead of fitting. When no person answers
  items in more than one block, the pairwise likelihood is flat in the
  between-block shift; the fit used to return arbitrary cross-block
  spacing (within-block orders fine, origins meaningless) with no
  warning.
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md),
  [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md),
  and everything built on them now stop with an error naming the
  disconnected blocks. Anchors rescue a block (disjoint-form equating):
  with an anchored item in every block the fit proceeds; with any block
  unanchored it stops and names it.
- Thresholds beside near-empty categories report NA, not a lie. A
  category with 1-2 responses can send its threshold toward a boundary
  (observed: a single top-category response put an item location at
  +13.9 logits) while the ridged covariance printed SE = 0. Such
  thresholds are now flagged (`thr$weak`), their SEs and the
  item-location SE reported as NA, and a note names the item, the
  category, and its count, pointing to `pc_components` or category
  collapsing. Estimates are still returned; categories with zero
  responses are rescored away as before.

## rasch 1.9.2

- Case study: `inst/casestudies/party_blocs_crisis.R` applies
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  to the Tuebingen 2009 party-preference data (via `psychotools`), with
  ideological blocs as object sets and economic-crisis concern as judge
  panels. Running the model against real data surfaced three defects,
  all fixed below.
- A set whose within-set contests are all near-even (or all one-sided)
  carries no stable information about the panel-unit ratios: its
  stage-one ratio estimate diverged and poisoned the phi reconciliation
  with a spurious covariance (on the party data, the right bloc’s single
  CDU/CSU-FDP pair at 55:45 drove phi to 3e4 and a fatal error). Such
  sets are now screened out of the reconciliation, refit with the panel
  units held at the reconciled phi – which the frame model says apply to
  them regardless – and named in a note. Stage one also gains the same
  trust-region step cap as
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).
- Bootstrap replicates in which a set unit reaches its boundary (log
  alpha driven to -Inf when a resampled within-set order flips against
  the cross-set evidence, the signature of a two-object set with a
  near-even internal pair) made the standard error silently NaN, printed
  as blank. The SE is now reported as NA with a note naming the
  parameter and the boundary count: no number is manufactured for a
  sampling distribution that is not normal. A replicate whose
  convergence flag was NA could also crash the bootstrap loop; it now
  counts as a failure.
- The stage convergence checks compared the gradient to an absolute
  threshold, which is scale-dependent: on large designs a converged fit
  could be flagged unconverged and – with the new screen – silently
  rerouted, changing the estimates. Both stages now use a per-comparison
  criterion, verified invariant to 50-fold duplication of the data.
- Documentation: the bootstrap’s two honesty notes (model-based
  resampling versus judge-clustered conditional errors, and
  boundary-unstable parameters) are spelled out in
  [`?btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md).

## rasch 1.9.1

- Shiny app: the Frames page gains its paired-comparison variant –
  nominate a judge-panel column and an object-set map (CSV, or inferred
  from object-name prefixes), estimate the frame units on demand
  (bootstrap standard errors by default), and read the panel units, set
  units and origins, the unit caterpillar, and the per-frame fit. On the
  app’s own paired-comparison demo the two judge panels differ
  significantly in discriminating power – the panel containing the
  erratic judge carries the smaller unit – which is precisely the
  phenomenon the model makes visible.

## rasch 1.9.0

- The extended frame of reference model for paired comparisons:
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  estimates object locations when the unit of the scale differs across
  judge-panel by object-set frames – the package’s extension of
  Humphry’s persons-by-items EFRM to the Bradley-Terry-Luce family.
  Panel units (phi, a panel’s discriminating power) are identified
  conditionally wherever two panels judge a common set; set units and
  origins (alpha, kappa) are identified from cross-set comparisons,
  which are themselves conditional outcomes – so, unlike the
  persons-by-items case, the linking makes no distributional assumption
  about the objects. Everything is estimated by the staged conditional
  estimator: the within-frame calibrations are invariant to the linking
  data (the frame-of-reference property), a deliberate trade of
  efficiency for invariance. The cross-frame discrimination convention
  (cross-set contests judged at exactly phi) is stated plainly in the
  documentation as an assumption, with the per-frame fit residuals as
  its check.
- Inference for
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  defaults to a parametric bootstrap of the whole two-stage pipeline:
  simulation showed the fast conditional errors cover at a third of
  nominal on chain-linked designs (they omit the stage-one uncertainty
  that flows into the linking), while the bootstrap restores 12/12
  coverage and nominal false-flag rates at both small and large designs.
  The conditional errors remain available as `se_method = "conditional"`
  for quick inspection, labelled as understating.
  [`plot_btl_units()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_units.md)
  draws the unit caterpillar;
  [`simulate_btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl_efrm.md)
  generates frame-structured data with planted units for the
  plant-to-detect loop.

## rasch 1.8.0

Three further paired-comparison analogues of the Rasch tool set, and
CRAN housekeeping.

- Anchored estimation for paired comparisons:
  `btl(anchors = c(Object = value, ...))` holds the named objects at
  fixed locations and estimates the rest on that scale (no sum-zero
  constraint; anchored objects report se 0) – the equating workflow’s
  second half, mirroring `rasch(anchors=)`. An anchored object at a
  response boundary is an error rather than a silent removal.
- First-position advantage: `btl(position = TRUE)` estimates the pure
  positional (first/left-presented) advantage of the Davidson & Beaver
  1977. order-effect device as a constant covariate alongside any
        exposure and carry-over effects; it is identified by triangle
        closure even under fixed presentation orders (more weakly, with
        honest SEs). The dimensionality reference simulation carries the
        fitted coefficient.
- Common-object equating:
  [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  compares two paired-comparison calibrations (or a calibration against
  a bank) through their common objects – precision-weighted origin
  shift, per-object drift t-tests with multiplicity adjustment, equated
  locations – with
  [`plot_btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_equate.md);
  the standards-maintenance workflow of comparative judgement across
  panels or years.
- Targeting and adaptive pairing:
  [`btl_information()`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md)
  (per-comparison Fisher information and per-object design information,
  with the naive 1/sqrt(information) SE beside the judge-clustered
  sandwich),
  [`plot_btl_targeting()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md),
  and
  [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  – the greedy most-informative-next-pair step of adaptive comparative
  judgement (Pollitt 2012), with the Bramley (2015)
  reliability-inflation caution documented.
- CRAN housekeeping: pure-ASCII sources, quoted ‘shiny’ in DESCRIPTION,
  heavy Monte-Carlo recovery tests wrapped in `skip_on_cran()` (they run
  locally and on CI), and `cran-comments.md` drafted from a full
  `--as-cran` check (about 3 minutes; slowest example 1.4s).
- An adversarial verification round over the new features also caught
  and fixed a PRE-EXISTING inference bug: `btl(count=)` inflated every
  sandwich standard error by about the square root of the count (the
  meat weighted aggregated rows twice; point estimates were always
  correct). Count- weighted fits now reproduce the expanded-data
  standard errors exactly; OSI and all z/p values follow. Additionally
  from that round:
  [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)’s
  shift standard error now uses the stored covariance of each
  calibration (it treated the common-object differences as independent,
  overstating the SE several-fold; 95 per cent coverage verified) and
  documents the majority-drift limitation of single-shift equating;
  [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  ranks by the exact one-step reduction in total location variance (a
  Sherman-Morrison update on the fit’s covariance), verified to beat
  low-priority pairs empirically; and
  [`btl_information()`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md)’s
  `se_naive` is documented honestly as a single-parameter lower bound
  rather than an independence benchmark.
- Shiny app: the Targeting and Equating tabs gain paired-comparison
  variants (design information + targeting plot + adaptive next-pair
  recommendations; reference-calibration upload with drift table and
  plot), and the Data roles gain a first-position-advantage switch and
  an object-anchors upload.

## rasch 1.7.1

A statistical audit of the new simulation and paired-comparison
independence code (23 findings across two adversarial rounds, all
fixed), plus a vignette and a logo.

- Misfit layers now compose correctly in
  [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md):
  the response-dependence and response-style regeneration steps
  previously rebuilt their target items from the base structure,
  silently erasing any DIF or second-dimension loading planted on those
  items (and leaking a DIF-like signal from a DIF-carrying source item
  into its dependent partner). Every regeneration now honours the item’s
  own generating structure – trait, group-shifted thresholds and slope,
  guessing.

- [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  now simulates its noise reference sequentially WITH the fitted
  within-judge dependence effects when the fit carries them: real order
  effects previously inflated the false-positive rate to ~36 per cent on
  one-dimensional data (now ~0-5 per cent). The documented price is
  conservatism, since carry-over and a judge-camp second attribute are
  partially confounded.

- Doc examples pass `id = "id"` so they no longer emit spurious
  dropped-column notes.

- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md):
  `model = "PCM"` now draws each item’s threshold spacings and span
  afresh, as the partial credit model allows (previously PCM and RSM
  generated identical rating-scale-structured data); `dif` without
  `n_groups >= 2` is an error instead of a silent no-op; polytomous
  `guessing` warns and is ignored rather than falsely recorded as
  planted; the response-dependence term now carries the guessing
  asymptote of the pair’s first item (no more spurious easiness on the
  dependent item); disordered thresholds work with 3 categories and no
  longer shift the item’s location; careless responders are no longer
  double-counted as response-style persons in the recorded truth.

- [`simulate_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_mfrm.md):
  halo raters are capped at the eligible pool (no NA raters when erratic
  raters shrink it); the example runs.

- [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md):
  person abilities are mean-centred as documented (an asymmetric
  difficulty range no longer masquerades as person bias); the many-facet
  branch reads the item margins (`item_effects`), restoring the
  documented item-difficulty recovery for MFRM fits.

- [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md):
  null calibration verified for graded fits (mildly conservative, never
  anticonservative) and documented; docs also clarify the
  leading-bimension reference columns and the no-ties precondition of
  Kendall’s zeta.

- Shiny app (Simulate page): the parameter-recovery card only renders
  when the current fit was estimated on the currently loaded simulation
  (a stale fit can no longer be compared against new truth); the
  reproducible call now carries every argument the generator received
  (MFRM person/item SDs, EFRM person SD, graded categories); non-integer
  seeds are rounded safely.

- New vignette: “Planting misfit and watching the diagnostics fire” —
  the plant-to-detect loop end to end, closing with a small Monte Carlo
  power estimate via
  [`sim_replicate()`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md).

- A hex logo (an item characteristic curve at its location),
  reproducible from `tools/logo.R`.

## rasch 1.7.0

- Simulation gains the standard controls of a proper simulation tool.
  - Population parameters: the person distribution (mean, SD, shape) and
    the item-difficulty range are now set directly (in the package and
    the app).
  - More misfit: response styles (extreme / midpoint categories) and
    speededness (a not-reached tail) for , and a halo effect for .
  - generates a batch of datasets for Monte Carlo power, Type-I, or
    recovery studies.
  - and compare the parameters a fit recovers with the ones that were
    planted (correlation, RMSE, bias, and a true-vs-estimated scatter)
    for every layout.
- Shiny app: the Simulate page renames the wide layout to “Rasch”,
  exposes the population parameters, shows a parameter-recovery report
  and the true-vs-estimated scatter after Run, and prints the
  reproducible call.

## rasch 1.6.1

- Shiny app: a Simulate page (under More) exposes all four simulators
  with a data-type selector and per-layout misfit controls. Pressing
  Simulate generates the data, loads it as the current dataset with its
  roles set, and lists what was planted; pressing Run then analyses it
  so the matching diagnostic can be watched firing. Verified end to end
  for every layout (e.g. a planted DIF item is flagged by dif_anova, a
  many-facet fit recovers its rater severities).

## rasch 1.6.0

- Data simulation across all four model families, each with dial-in
  departures from the model – so a known pathology can be planted and
  the matching diagnostic watched as it fires. Every result carries the
  true parameters () and prints what was planted.
  - : dichotomous/PCM/RSM with discrimination misfit, guessing, a second
    dimension, local dependence, DIF (uniform and non-uniform), careless
    responders, disordered thresholds, and missing data.
  - : paired comparisons with erratic judges, a second object attribute
    (two judge camps), and within-judge exposure/carry-over dependence;
    dichotomous or graded.
  - : rated data with rater-severity spread, erratic raters, and a
    rater-by-item interaction.
  - : frames with differing set and group units. Verified end to end:
    each planted departure is recovered or flagged by the corresponding
    fit / diagnostic (an app panel follows).

## rasch 1.5.1

- The judge unexpected-judgements map () is now drawn at the matchup
  level: each pair the judge met is a segment spanning its two objects
  on the consensus location axis, placed by how surprising the verdict
  was, with the judge’s backed object filled. returns the per-matchup
  residuals behind it. The object-level is retained for per-object
  residuals.

## rasch 1.5.0

- [`judge_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md)
  and
  [`plot_btl_judge_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_judge_map.md):
  the paired-comparison counterpart of the kidmap. A judge has no
  ability to condition on, so the reference is the consensus object
  scale; each object a judge met gets a standardised residual saying how
  much more or less that judge favoured it than its location predicts. A
  strong object the judge under-rated, or a weak object over-rated, is
  an unexpected judgement. It catches systematic bias – a distinct
  failure mode from the noise that judge fit and judge consistency flag.
  In the Shiny app the Persons \> Judge fit panel is now a
  master-detail: the judges table drives the selected judge’s map,
  mirroring person -\> kidmap.

## rasch 1.4.0

- Paired comparisons gain the pair-structure analogues of the Rasch
  independence diagnostics – the tools that a persons-by-items residual
  matrix would give but paired-comparison data cannot.
  - [`btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md)
    counts circular triads (preference loops: A beats B, B beats C, C
    beats A) against the chance rate, reports Kendall’s coefficient of
    consistency for complete designs, and gives each judge’s own
    consistency (Kendall & Babington Smith 1940) – the single-dimension
    check, and a judge-fit analogue.
  - [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
    decomposes the skew-symmetric object-by-object residual preference
    matrix into Gower (1977) bimensions – the paired-comparison
    counterpart of residual PCA – and judges the leading bimension
    against a model-simulated noise reference, so a coherent residual
    “swirl” (a second attribute steering some contests) is distinguished
    from noise.
  - [`plot_btl_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_scree.md),
    [`plot_btl_dim_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dim_map.md),
    and
    [`plot_btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_transitivity.md)
    display them; the Shiny app’s Independence \> Trait tab now carries
    these for a BTL fit (it previously hid for paired comparisons).

## rasch 1.3.1

- The package is now called **rasch** (it was developed under the
  working title `rmt`). The name `rmt` belongs to an unrelated package
  on CRAN, and `rasch` says plainly what this is – and mirrors
  `mirt::mirt()` with
  [`rasch::rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).
  No functions, arguments, or results change; only the package and
  library name. The auxiliary result classes were unified under the
  `rasch_` prefix (`rasch_btl`, `rasch_dif`, …) to match the fit
  classes.

## rasch 1.3.0

A full statistical audit of the independence, invariance, and paired-
comparison procedures, with fixes throughout.

- Paired comparisons:
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) now
  interrogates within-judge dependence – `dependence_data` holds every
  comparison’s exposure and carry-over covariates (count-weighted),
  [`plot_btl_dependence()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dependence.md)
  displays an effect as a partial-residual curve, and a separated effect
  (one-sided informative comparisons) is set aside with a note instead
  of reported as a runaway estimate. The graded engine’s
  threshold-by-dependence Hessian block is corrected (free-threshold
  fits with five or more categories and an order column now converge),
  rows with a missing judge are dropped with a note, category checks
  re-run after boundary objects are removed, and the dichotomous path
  now routes through the one (vectorised) estimator.
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  holds fitted dependence effects fixed in its residual moments and
  resolution refits (dependence is no longer absorbed as spurious
  judge-group DIF), weights count-aggregated comparisons exactly as
  expanded rows, restricts supersession to group terms, and notes every
  significant term it does not resolve.
- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  keeps each ANOVA term in its own error stratum, so a mixed
  (within-subject) design with missing responses no longer un-flags real
  DIF; factors named like internals (`ci`, `f1`) are safe; a term’s own
  class-interval crossing no longer supersedes it (so
  [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  again splits items with both uniform and non-uniform DIF); mixed
  designs say plainly that Tukey comparisons are unavailable; and
  [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  notes when a within-person factor makes its standard errors
  conservative.
- [`plot_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)’s
  parallel-analysis reference is now model-simulated (responses drawn
  from the calibrated model, persons re-estimated), so model-true data
  sits at the reference instead of always “showing structure” (Raiche
  2005; Chou & Wang 2010).
- Shiny app: judge-group DIF results freeze their run’s factors and
  alpha (later sidebar edits cannot silently regroup the overlay or
  break the reproducible snippets), superseded terms are shown, the
  within-judge dependence panel explains dropped effects, and fit flags
  are consistent on every table (fit residual \|2.5\|; outfit 0.7-1.3;
  infit 0.8-1.2).

## rasch 1.2.0

- [`plot_pca_biplot()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pca_biplot.md)
  draws the item loadings on the first two residual principal components
  – the pair that usually carries any interpretable contrast – on equal
  axes, coloured by the sign of the PC1 loading.

- [`residual_correlations()`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md)
  now also returns the adjusted-Q3 `star_matrix` (each Q3 less the
  average off-diagonal), and
  [`plot_resid_cor()`](https://drjoshmcgrane.github.io/rasch/reference/plot_resid_cor.md)
  gains a `stat` argument to colour either the raw Q3 or the adjusted
  Q3\*, drawing the lower triangle only.

- Shiny app: the trait- and local-dependence pages pair each table with
  its plot (loadings with a PC1-PC2 biplot; Q3 and Q3\* matrices with
  heatmaps), the raw and adjusted correlations are each flagged by their
  own rule (\|Q3\| and Q3\*), and every data-restructuring action
  (subtest, item split, automatic DIF resolution) offers an in-place
  reset to the original data.

- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  is now the single DIF analysis-of-variance function. One factor is
  analysed one-way; several factors are modelled jointly – the
  statistically correct treatment – with main effects by default and
  factor-by-factor interactions optional (`effects = "factorial"`). It
  handles within-subject factors (repeated-measures / mixed ANOVA) and
  returns a classed object with a `summary`, the full `terms` table, and
  Tukey comparisons. The separate `dif_anova_factorial()` is removed.

- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  resolves DIF iteratively by item splitting, largest effect first, to
  clear artificial DIF.

## rasch 1.0.0

First stable release. The package delivers a complete Rasch Measurement
Theory workflow built entirely from published measurement theory, with a
pairwise conditional estimation core and a Shiny interface that exposes
every analysis with reproducing R code attached to every output.

### Models

- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md):
  dichotomous and polytomous (partial credit and rating scale) Rasch
  models by pairwise conditional maximum likelihood (Andrich & Luo 2003;
  Zwinderman 1995), with Godambe sandwich standard errors, sum-zero
  identification, Warm (1989) weighted likelihood person estimation, and
  the principal-component threshold parameterisation as an estimation
  option.
- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md):
  many-facet models with additive facet severities or item-by-facet
  interactions, accepting wide (one column per item) or long data.
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md):
  the extended frame-of-reference model (Humphry) with frame units
  estimated by within-frame pairwise conditional likelihood.
- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md):
  paired comparisons as the conditional form of the dichotomous Rasch
  model (Bradley & Terry 1952; Andrich 1978), including the
  adjacent-categories graded extension (Tutz 1986; Agresti 1992) with
  symmetric thresholds – the Davidson (1970) ties model is its
  three-category case – winner-plus-margin data entry, judge-clustered
  errors, and within-judge exposure and carry-over effects estimated
  jointly with the locations.

### Diagnostics

- The test-of-fit suite follows Andrich & Marais (2019):
  log-of-mean-square fit residuals with apportioned degrees of freedom,
  the item-trait chi-square on automatically sized, tie-preserving class
  intervals with per-interval detail, the ANOVA item fit, and person
  fit.
- Unidimensionality: residual principal components with parallel
  analysis, the Smith (2002) subset t-test, and magnitude estimation.
- Local dependence: Yen’s (1984) Q3 with the Christensen, Makransky &
  Horton (2017) flagging convention, response-dependence magnitude
  (Andrich & Kreiner), subtest formation, and the spread test.
- Guessing: tailored analysis with origin-equated comparison
  calibrations.
- DIF: two-way residual ANOVA per factor, the joint factorial analysis
  with a compact uniform/non-uniform summary, resolved DIF magnitudes in
  logits with familywise control and a practical-significance criterion,
  planned contrasts derived from the factor structure (with
  repeated-measures support via person-level residual scores), and class
  intervals sized to the cells each analysis actually uses. For paired
  comparisons,
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  tests object-by-judge-group interaction by the same two routes.
- Equating and invariance: common-item comparison with drift flags, and
  calibration anchoring.

### Display and reporting

- A complete base-graphics plot suite: expected value curves with
  observed class-interval means, category and threshold probability
  curves, the person-item threshold distribution, Wright maps, kidmaps
  (Wright, Mead & Ludlow 1980), object characteristic curves for paired
  comparisons, threshold and item maps, test characteristic and
  information curves, residual displays, and Guttman scalograms – all
  with adjustable class intervals and scale ranges, and batch export to
  multi-page PDF or ZIP at publication resolution.
- [`fit_summary_table()`](https://drjoshmcgrane.github.io/rasch/reference/fit_summary_table.md)
  and
  [`targeting_table()`](https://drjoshmcgrane.github.io/rasch/reference/targeting_table.md)
  return the headline statistics as tidy tables;
  [`save_outputs()`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md)
  writes every table and plot to disk;
  [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  produces a single-file HTML report.
- [`run_app()`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md)
  launches the Shiny interface: model-adaptive navigation, reproducing R
  code beneath every output, and one-click export.
