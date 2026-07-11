# Data simulation: each planted departure must trip its matching diagnostic
# (the plant -> detect loop), and the truth is attached for recovery checks.

test_that("simulate_rasch plants misfit the Rasch diagnostics detect", {
  # discrimination: a central over-discriminating item overfits (low outfit)
  d <- simulate_rasch(600, 11, discrimination = c(rep(1, 5), 3, rep(1, 5)),
                      seed = 1)
  f <- rasch(d)
  expect_lt(f$items$outfit_ms[6], 0.7)
  expect_s3_class(d, "rasch_sim")

  # DIF flags the planted item and (essentially) nothing else
  d <- simulate_rasch(800, 10, dif = list(items = "I05", uniform = 1.2),
                      n_groups = 2, seed = 2)
  s <- dif_anova(rasch(d, factors = "group"))$summary
  expect_true(s$uniform_DIF[s$item == "I05"])
  expect_equal(sum(s$uniform_DIF[s$item != "I05"]), 0L)

  # local dependence flags the planted pair on Q3*
  d <- simulate_rasch(1000, 10,
                      dependence = list(pairs = list(c("I03", "I04")),
                                        strength = 2.5), seed = 4)
  h <- residual_correlations(rasch(d))$flagged
  expect_true(any((h$item_a == "I03" & h$item_b == "I04") |
                  (h$item_a == "I04" & h$item_b == "I03")))

  # careless responders inflate person outfit
  d <- simulate_rasch(600, 15, careless = 0.12, seed = 6)
  po <- rasch(d)$person$outfit_ms; ci <- attr(d, "truth")$careless_idx
  expect_gt(mean(po[ci], na.rm = TRUE), mean(po[-ci], na.rm = TRUE) + 0.5)

  expect_output(print(d), "careless")
})

test_that("simulate_btl plants misfit the paired-comparison diagnostics detect", {
  # erratic judges carry large fit residuals and low consistency
  d <- simulate_btl(8, 12, erratic_judges = 0.17, seed = 1)
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  er <- attr(d, "truth")$erratic
  expect_gt(mean(bt$judges$fit_resid[bt$judges$judge %in% er]),
            mean(bt$judges$fit_resid[!bt$judges$judge %in% er]) + 1)

  # graded comparisons recover the object locations
  d <- simulate_btl(8, 12, model = "graded", n_categories = 4, seed = 4)
  bt <- btl(d, "object_a", "object_b", response = "response", judge = "judge")
  expect_gt(cor(bt$objects$location,
                attr(d, "truth")$location[bt$objects$object]), 0.95)

  # a planted carry-over dependence is recovered
  d <- simulate_btl(6, 10, reps_per_pair = 40,
                    dependence = list(exposure = 0, carry_over = 1.2), seed = 3)
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "order")
  co <- bt$dependence$estimate[bt$dependence$effect == "carry_over"]
  expect_gt(co, 0.5)
})

test_that("simulate_mfrm plants rater severity, misfit, and interaction", {
  d <- simulate_mfrm(120, 5, 6, rater_severity_sd = 0.7, erratic_raters = 0.17,
                     interaction = list(rater = "R3", item = "I2", bias = 1.8),
                     seed = 1)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater", interaction = "rater")
  tr <- attr(d, "truth"); fe <- mf$facet_effects$rater
  rec <- fe$severity[match(names(tr$severity), fe$level)]
  # severities are recovered for the well-behaved raters (the erratic one's
  # true severity is meaningless once it rates at random)
  keep <- !(names(tr$severity) %in% tr$erratic)
  expect_gt(abs(cor(rec[keep], tr$severity[keep])), 0.9)
  er <- tr$erratic
  expect_gt(mean(fe$fit_resid[fe$level %in% er]),
            mean(fe$fit_resid[!fe$level %in% er]) + 1)      # erratic rater misfits
  ie <- mf$interaction_effects                             # interaction at R3xI2
  expect_equal(ie[which.max(abs(ie$gamma)), c("item", "level")],
               data.frame(item = "I2", level = "R3"), ignore_attr = TRUE)
})

test_that("simulate_efrm plants a frame-unit ratio rasch_efrm recovers", {
  d <- simulate_efrm(300, 8, set_unit_ratio = 1.35, seed = 2)
  tr <- attr(d, "truth")
  ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group")
  ratio <- max(ef$alpha_table$alpha) / min(ef$alpha_table$alpha)
  expect_gt(ratio, 1.2); expect_lt(ratio, 1.55)          # ~1.35 recovered
  expect_output(print(d), "set-unit ratio")
})

test_that("the extra misfit types plant detectable signals", {
  # extreme response style: style persons over-use the end categories
  d <- simulate_rasch(600, 12, model = "PCM", n_categories = 4,
                      response_style = list(type = "extreme", prop = 0.3), seed = 1)
  si <- attr(d, "truth")$style_idx; cats <- as.matrix(d[, grep("^I", names(d))])
  expect_gt(mean(cats[si, ] %in% c(0, 3)), mean(cats[-si, ] %in% c(0, 3)) + 0.1)

  # speededness: a missing tail growing toward the last item
  d <- simulate_rasch(800, 15, speeded = 0.5, seed = 2)
  miss <- colMeans(is.na(as.matrix(d[, grep("^I", names(d))])))
  expect_lt(miss[8], 0.02)
  expect_gt(miss[15], 0.3)
  expect_true(miss[15] > miss[13] && miss[13] > miss[11])   # monotone gradient

  # MFRM halo: halo raters barely differentiate items -> large interaction
  d <- simulate_mfrm(140, 6, 6, rater_severity_sd = 0.5, item_sd = 1.3,
                     halo = 0.17, seed = 5)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater", interaction = "rater")
  hr <- attr(d, "truth")$halo; ie <- mf$interaction_effects
  expect_gt(mean(abs(ie$gamma[ie$level %in% hr])),
            2 * mean(abs(ie$gamma[!ie$level %in% hr])))
})

test_that("sim_replicate and sim_recovery support Monte Carlo and recovery", {
  b <- sim_replicate(simulate_rasch, 6, n_persons = 300, n_items = 8, seed = 1)
  expect_s3_class(b, "rasch_sim_batch")
  expect_length(b, 6)
  expect_false(identical(b[[1]]$I01, b[[2]]$I01))          # different datasets

  # recovery: a clean fit gets its planted parameters back
  d <- simulate_rasch(600, 12, seed = 1)
  rec <- sim_recovery(rasch(d), d)
  expect_s3_class(rec, "rasch_recovery")
  s <- rec$summary
  expect_gt(s$correlation[s$parameter == "item difficulty"], 0.95)
  # person ability is noisier (WLE precision from only 12 items limits it)
  expect_gt(s$correlation[s$parameter == "person ability"], 0.75)
  expect_lt(abs(s$bias[s$parameter == "item difficulty"]), 0.1)
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot_recovery(rec))

  # recovery across the other layouts
  d <- simulate_btl(8, 12, seed = 2)
  rb <- sim_recovery(btl(d, "object_a", "object_b", winner = "winner",
                         judge = "judge"), d)
  expect_gt(rb$summary$correlation[1], 0.9)
})

test_that("audit fixes hold: PCM structure, truth honesty, recovery centring", {
  # PCM and RSM now genuinely differ: per-item threshold patterns for PCM,
  # one common pattern for RSM; PCM thresholds stay ordered
  d1 <- simulate_rasch(200, 8, model = "PCM", n_categories = 4, seed = 42)
  d2 <- simulate_rasch(200, 8, model = "RSM", n_categories = 4, seed = 42)
  expect_false(identical(as.matrix(d1[, 2:9]), as.matrix(d2[, 2:9])))
  rel1 <- lapply(attr(d1, "truth")$thresholds, function(t) round(t - mean(t), 3))
  expect_gt(length(unique(rel1)), 1)                       # PCM varies
  rel2 <- lapply(attr(d2, "truth")$thresholds, function(t) round(t - mean(t), 3))
  expect_length(unique(rel2), 1)                           # RSM common
  expect_true(all(vapply(attr(d1, "truth")$thresholds,
                         function(t) !is.unsorted(t), TRUE)))

  # dif without groups is an error, not a silent no-op with a false truth
  expect_error(simulate_rasch(50, 6, dif = list(items = "I03", uniform = 1)),
               "n_groups")
  # polytomous guessing warns and is not recorded as planted
  expect_warning(d <- simulate_rasch(50, 4, model = "PCM", n_categories = 4,
                                     guessing = 0.3, seed = 1), "dichotomous")
  expect_false(any(grepl("guessing", attr(d, "truth")$planted)))
  # disordered thresholds work at 3 categories and preserve the location
  t2 <- rasch:::.sim_thresholds(0.5, 2, 1.2, disordered = TRUE)
  expect_true(is.unsorted(t2)); expect_equal(mean(t2), 0.5)
  # careless overwrite is not double-counted in the truth
  d <- simulate_rasch(300, 10, model = "PCM", n_categories = 4, careless = 0.5,
                      response_style = list(type = "extreme", prop = 0.5),
                      seed = 3)
  tr <- attr(d, "truth")
  expect_length(intersect(tr$style_idx, tr$careless_idx), 0)
  # halo raters never overflow into NA when erratic raters shrink the pool
  d <- simulate_mfrm(30, 4, 5, erratic_raters = 0.4, halo = 0.8, seed = 1)
  expect_false(anyNA(attr(d, "truth")$halo))

  # person ability is centred in recovery: an asymmetric difficulty range
  # must not masquerade as person-ability bias
  d <- simulate_rasch(400, 10, difficulty = c(0, 3), seed = 2)
  r <- sim_recovery(rasch(d), d)
  expect_lt(abs(r$summary$bias[r$summary$parameter == "person ability"]), 0.1)
  # MFRM recovery reports item difficulties from the item margins
  d <- simulate_mfrm(60, 5, 5, seed = 1)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater")
  r <- sim_recovery(mf, d)
  expect_true("item difficulty" %in% r$summary$parameter)
  expect_gt(r$summary$correlation[r$summary$parameter == "item difficulty"], 0.9)
})

test_that("misfit layers compose: dependence and style respect DIF / 2nd dim", {
  # dependence regeneration keeps the target item's DIF
  d <- simulate_rasch(2000, 6, dif = list(items = "I02", uniform = 1.5),
                      n_groups = 2,
                      dependence = list(pairs = list(c("I01", "I02")),
                                        strength = 1), seed = 7)
  g <- attr(d, "truth")$groups
  expect_gt(mean(d$I02[g != "g2"]) - mean(d$I02[g == "g2"]), 0.12)
  # ...and the target item's second dimension
  d <- simulate_rasch(2000, 6, second_dim = list(items = "I02", rho = 0.2),
                      dependence = list(pairs = list(c("I01", "I02")),
                                        strength = 1), seed = 8)
  expect_lt(cor(d$I02, attr(d, "truth")$theta), 0.2)
  # response style keeps DIF for style-affected persons
  d <- simulate_rasch(3000, 6, model = "PCM", n_categories = 4,
                      dif = list(items = "I02", uniform = 2.5), n_groups = 2,
                      response_style = list(type = "extreme", prop = 0.5),
                      seed = 10)
  tr <- attr(d, "truth"); g <- tr$groups
  sty <- seq_len(nrow(d)) %in% tr$style_idx
  expect_gt(mean(d$I02[sty & g != "g2"]) - mean(d$I02[sty & g == "g2"]), 0.8)
})

test_that("btl_dimensionality reference honours fitted dependence effects", {
  skip_on_cran()   # heavy simulation; verified locally and on CI
  # one-dimensional data whose only structure is within-judge order effects,
  # fitted WITH order: the dependence-aware reference must not read the
  # order structure as a second attribute
  flags <- vapply(1:4, function(s) {
    d <- simulate_btl(8, 12, reps_per_pair = 30,
                      dependence = list(exposure = 1, carry_over = 0.8),
                      seed = 300 + s)
    bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
              order = "order")
    isTRUE(btl_dimensionality(bt, reps = 50)$leading_structured)
  }, TRUE)
  expect_lte(sum(flags), 1L)   # was ~36% false-positive before the fix
})
