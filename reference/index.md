# Package index

## Estimation

Pairwise conditional estimation of the Rasch model family, with Godambe
sandwich standard errors and Warm person measures.

- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  : Fit and diagnose a Rasch model by pairwise conditional estimation
- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  : Fit a many-facet Rasch model
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  : Fit the extended frame of reference model
- [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md) :
  Estimate Rasch thresholds by pairwise conditional maximum likelihood
- [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  : Estimate Rasch thresholds via the Andrich principal-components
  reparameterisation
- [`threshold_index()`](https://drjoshmcgrane.github.io/rasch/reference/threshold_index.md)
  : Enumerate item-category thresholds
- [`item_moments()`](https://drjoshmcgrane.github.io/rasch/reference/item_moments.md)
  : Category-score moments for a polytomous item

## Test of fit and comparison

The fit-residual and chi-square apparatus, targeting and reliability,
and model comparison.

- [`fit_summary_table()`](https://drjoshmcgrane.github.io/rasch/reference/fit_summary_table.md)
  : Test-of-fit summary as a table
- [`targeting_table()`](https://drjoshmcgrane.github.io/rasch/reference/targeting_table.md)
  : Targeting and reliability summary as a table
- [`chisq_detail()`](https://drjoshmcgrane.github.io/rasch/reference/chisq_detail.md)
  : Class-interval detail for one item's chi-square test of fit
- [`test_information()`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md)
  : Test information function
- [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
  : Likelihood-ratio test of the partial credit against the rating scale
  model
- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  : Compare fitted Rasch models
- [`guttman_table()`](https://drjoshmcgrane.github.io/rasch/reference/guttman_table.md)
  : Guttman-ordered response matrix and reproducibility

## Persons

- [`person_wle()`](https://drjoshmcgrane.github.io/rasch/reference/person_wle.md)
  : Warm's weighted likelihood estimates by raw score
- [`person_extrapolated()`](https://drjoshmcgrane.github.io/rasch/reference/person_extrapolated.md)
  : Person measures with extrapolated extreme scores
- [`score_table()`](https://drjoshmcgrane.github.io/rasch/reference/score_table.md)
  : Raw score to measure conversion table

## Invariance and DIF

Differential item functioning over any number of person factors, DIF
magnitudes in logits, resolution by item splitting, and equating.

- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  : Differential item functioning by residual analysis of variance
- [`dif_contrasts()`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
  : Planned DIF contrasts derived from the factor structure
- [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  : DIF magnitude in logits with pairwise comparisons
- [`split_items()`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md)
  : Split items by a person factor to resolve DIF
- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  : Resolve differential item functioning by iterative item splitting
- [`equate_tests()`](https://drjoshmcgrane.github.io/rasch/reference/equate_tests.md)
  : Equate two test calibrations through their common items
- [`tailored_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/tailored_analysis.md)
  : Tailored analysis for guessing

## Independence and dimensionality

Residual principal components, the Smith t-test with magnitude
estimation, Q3 local dependence, and the structural remedies.

- [`residual_correlations()`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md)
  : Residual correlations for local dependence (Yen's Q3)
- [`residual_pca()`](https://drjoshmcgrane.github.io/rasch/reference/residual_pca.md)
  : Principal components of the residual correlations
- [`dimensionality_test()`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md)
  : Residual-component test of unidimensionality
- [`dimensionality_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_magnitude.md)
  : Magnitude of multidimensionality from a subtest analysis
- [`dependence_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dependence_magnitude.md)
  : Estimate the magnitude of response dependence between two items
- [`spread_test()`](https://drjoshmcgrane.github.io/rasch/reference/spread_test.md)
  : Spread-parameter test for dependence within subtests
- [`combine_items()`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md)
  : Combine items into subtests and re-analyse
- [`rack_data()`](https://drjoshmcgrane.github.io/rasch/reference/rack_data.md)
  [`stack_data()`](https://drjoshmcgrane.github.io/rasch/reference/rack_data.md)
  : Reshape repeated measurements for racked or stacked analysis

## Multiple choice and traditional statistics

- [`distractor_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/distractor_analysis.md)
  : Distractor analysis for multiple-choice items
- [`distractor_rescore()`](https://drjoshmcgrane.github.io/rasch/reference/distractor_rescore.md)
  : Propose polytomous option scores from the distractor evidence
- [`ctt_table()`](https://drjoshmcgrane.github.io/rasch/reference/ctt_table.md)
  : Traditional (classical test theory) statistics

## Paired comparisons

The Bradley-Terry-Luce model as the conditional form of the dichotomous
Rasch model, with judge diagnostics, within-judge dependence,
judge-group DIF, and the pair-structure analogues of the independence
diagnostics (transitivity and the residual bimension decomposition).

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) :
  Fit the Bradley-Terry-Luce model to paired comparisons
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  : DIF analysis for paired comparisons
- [`btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md)
  : Transitivity of paired comparisons
- [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  : Residual dimensionality of paired comparisons
- [`judge_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md)
  : Unexpected judgements of one judge
- [`judge_pair_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_pair_surprise.md)
  : Unexpected judgements of one judge, pair by pair
- [`plot_btl()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl.md)
  : Plot Bradley-Terry-Luce object locations
- [`plot_btl_categories()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_categories.md)
  : Plot graded-comparison category curves
- [`plot_btl_icc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_icc.md)
  : Plot an object characteristic curve
- [`plot_btl_dependence()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dependence.md)
  : Plot a within-judge dependence effect
- [`plot_btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_transitivity.md)
  : Consistency plot for paired-comparison transitivity
- [`plot_btl_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_scree.md)
  : Scree of paired-comparison residual bimensions
- [`plot_btl_dim_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dim_map.md)
  : Residual map of the leading paired-comparison bimension
- [`plot_btl_judge_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_judge_map.md)
  : Unexpected-judgement map for one judge (pair level)

## Plots

One plotting function per display, all base graphics.

- [`plot_catfreq()`](https://drjoshmcgrane.github.io/rasch/reference/plot_catfreq.md)
  : Plot category frequencies
- [`plot_ccc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_ccc.md)
  : Plot category probability curves
- [`plot_distractors()`](https://drjoshmcgrane.github.io/rasch/reference/plot_distractors.md)
  : Plot multiple-choice option curves
- [`plot_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_equate.md)
  : Plot a test-equating comparison
- [`plot_facets()`](https://drjoshmcgrane.github.io/rasch/reference/plot_facets.md)
  : Plot facet severities
- [`plot_frames()`](https://drjoshmcgrane.github.io/rasch/reference/plot_frames.md)
  : Plot frame units
- [`plot_guttman()`](https://drjoshmcgrane.github.io/rasch/reference/plot_guttman.md)
  : Plot the Guttman scalogram
- [`plot_icc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_icc.md)
  : Plot an item characteristic curve
- [`plot_icc_frames()`](https://drjoshmcgrane.github.io/rasch/reference/plot_icc_frames.md)
  : Plot an item's characteristic curves across frames
- [`plot_item_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_item_map.md)
  : Plot the item map (location against fit residual)
- [`plot_kidmap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_kidmap.md)
  : Plot a kidmap
- [`plot_pca()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pca.md)
  : Plot residual principal-component loadings
- [`plot_pca_biplot()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pca_biplot.md)
  : Biplot of the first two residual components
- [`plot_pcc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pcc.md)
  : Plot a person characteristic curve
- [`plot_person_fit()`](https://drjoshmcgrane.github.io/rasch/reference/plot_person_fit.md)
  : Plot person fit
- [`plot_pimap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pimap.md)
  : Plot the person-item threshold distribution
- [`plot_resid_cor()`](https://drjoshmcgrane.github.io/rasch/reference/plot_resid_cor.md)
  : Plot the residual-correlation heatmap
- [`plot_resid_dist()`](https://drjoshmcgrane.github.io/rasch/reference/plot_resid_dist.md)
  : Plot the fit residual distribution
- [`plot_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)
  : Scree plot of the residual components with parallel analysis
- [`plot_tcc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_tcc.md)
  : Plot the test characteristic curve
- [`plot_threshold_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_threshold_map.md)
  : Plot the threshold map
- [`plot_threshold_prob()`](https://drjoshmcgrane.github.io/rasch/reference/plot_threshold_prob.md)
  : Plot threshold probability curves
- [`plot_tif()`](https://drjoshmcgrane.github.io/rasch/reference/plot_tif.md)
  : Plot the test information function
- [`plot_wright()`](https://drjoshmcgrane.github.io/rasch/reference/plot_wright.md)
  : Plot a Wright map

## Export and interface

- [`save_outputs()`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md)
  : Save every output of a Rasch analysis to a folder
- [`save_item_plots()`](https://drjoshmcgrane.github.io/rasch/reference/save_item_plots.md)
  : Save a plot for every item
- [`save_person_plots()`](https://drjoshmcgrane.github.io/rasch/reference/save_person_plots.md)
  : Save a kidmap for every person
- [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  : Write a self-contained HTML report of a Rasch analysis
- [`run_app()`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md)
  : Launch the rasch graphical interface
