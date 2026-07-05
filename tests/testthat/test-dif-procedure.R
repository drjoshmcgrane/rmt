# The full DIF procedure: single-factor two-way ANOVA with the
# class-interval row and effect sizes; x-way factorial with interactions;
# post-hoc pairwise comparisons for significant interactions and for
# multi-level main effects, familywise-corrected; and DIF magnitudes in
# logits for practical significance.

sim_dif <- function(n = 900, seed = 2, shifts) {
  # shifts: named list group-label -> per-item logit shift vector
  set.seed(seed)
  d <- seq(-1.5, 1.5, length.out = 8)
  g <- factor(rep(names(shifts), each = n / length(shifts)))
  th <- rnorm(n)
  X <- matrix(NA_integer_, n, 8)
  for (lv in names(shifts))
    X[g == lv, ] <- matrix(rbinom(sum(g == lv) * 8, 1,
      plogis(outer(th[g == lv], d + shifts[[lv]], "-"))), sum(g == lv), 8)
  colnames(X) <- paste0("I", 1:8)
  list(X = X, g = g)
}

test_that("dif_anova reports the full two-way table with effect sizes", {
  s <- sim_dif(shifts = list(a = rep(0, 8), b = c(0, 0, 1, rep(0, 5))))
  fit <- rasch(data.frame(s$X, grp = s$g), factors = "grp")
  da <- dif_anova(fit)
  expect_true(all(c("F_class", "p_class", "eta2_uniform",
                    "eta2_nonuniform") %in% names(da)))
  # the planted item carries by far the largest uniform effect size
  expect_equal(which.max(da$eta2_uniform), 3L)
  expect_gt(da$eta2_uniform[3], max(da$eta2_uniform[-3]) * 3)
  expect_true(da$uniform_DIF[3])
  # eta2 within (0, 1); class-interval F present for every item
  expect_true(all(da$eta2_uniform > 0 & da$eta2_uniform < 1, na.rm = TRUE))
  expect_true(all(is.finite(da$F_class)))
  # familywise option flows through
  db <- dif_anova(fit, p_adjust = "bonferroni")
  expect_true(all(db$p_uniform_adj >= da$p_uniform_adj - 1e-12, na.rm = TRUE))
})

test_that("dif_size recovers a planted uniform DIF in logits", {
  s <- sim_dif(shifts = list(a = rep(0, 8), b = c(0, 0, 0.8, rep(0, 5))))
  fit <- rasch(data.frame(s$X, grp = s$g), factors = "grp")
  ds <- dif_size(fit, "I3", by = "grp")
  expect_equal(nrow(ds$pairs), 1)
  expect_lt(abs(abs(ds$pairs$difference) - 0.8), 0.25)
  expect_true(ds$pairs$significant)
  expect_true(ds$pairs$practical)          # 0.8 > 0.5 logits
  expect_true(ds$pairs$lower < ds$pairs$difference &
              ds$pairs$difference < ds$pairs$upper)
  # a clean item shows neither statistical nor practical DIF
  ds0 <- dif_size(fit, "I6", by = "grp")
  expect_false(ds0$pairs$significant)
  expect_false(ds0$pairs$practical)
})

test_that("multi-level factors get familywise pairwise comparisons in logits", {
  s <- sim_dif(n = 1200, seed = 5,
               shifts = list(a = rep(0, 8),
                             b = c(0, 0.5, rep(0, 6)),
                             c = c(0, 1.0, rep(0, 6))))
  fit <- rasch(data.frame(s$X, grp = s$g), factors = "grp")
  ds <- dif_size(fit, "I2", by = "grp")
  expect_equal(nrow(ds$levels), 3)
  expect_equal(nrow(ds$pairs), 3)          # all pairs of three levels
  d_ac <- ds$pairs$difference[ds$pairs$level_a == "a" & ds$pairs$level_b == "c"]
  d_ab <- ds$pairs$difference[ds$pairs$level_a == "a" & ds$pairs$level_b == "b"]
  expect_lt(abs(abs(d_ac) - 1.0), 0.35)
  expect_lt(abs(abs(d_ab) - 0.5), 0.35)
  expect_gt(abs(d_ac), abs(d_ab))          # graded shifts recovered in order
  # Holm adjustment is monotone in the raw p over the family
  expect_true(all(ds$pairs$p_adj >= ds$pairs$p - 1e-12))
  # the extreme pair is flagged practical, and a-c also significant
  sel <- ds$pairs$level_a == "a" & ds$pairs$level_b == "c"
  expect_true(ds$pairs$significant[sel] && ds$pairs$practical[sel])
})

test_that("factorial procedure: interaction post-hocs and sizes for significant terms", {
  set.seed(7); n <- 1600
  d <- seq(-1.5, 1.5, length.out = 6)
  g1 <- factor(rep(c("a", "b"), each = n / 2))
  g2 <- factor(rep(c("x", "y"), times = n / 2))
  th <- rnorm(n)
  sh <- ifelse(g1 == "b" & g2 == "y", 1.1, 0)   # DIF in one cell only
  X <- matrix(rbinom(n * 6, 1,
    plogis(outer(th, d, "-") - outer(sh, c(0, 0, 1, 0, 0, 0)))), n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
  fa <- dif_anova_factorial(fit, sizes = TRUE)

  # the g1:g2 interaction is significant for the planted item and
  # supersedes the main effects it involves
  t3 <- fa$terms[fa$terms$item == "I3", ]
  expect_true(t3$significant[t3$term == "g1:g2"])
  expect_true(all(c("eta2_partial") %in% names(fa$terms)))
  expect_true(t3$eta2_partial[t3$term == "g1:g2"] > 0)
  sup <- t3$term[t3$superseded]
  expect_true(all(sup %in% c("g1", "g2")))

  # Tukey post-hocs exist for the interaction cells of the planted item
  tk3 <- fa$tukey[fa$tukey$item == "I3" & fa$tukey$term == "g1:g2", ]
  expect_equal(nrow(tk3), 6)               # 4 cells -> 6 pairs
  worst <- tk3$comparison[which.min(tk3$p_tukey)]
  expect_true(grepl("b:y", worst))

  # sizes: logit magnitudes for the significant term, the b:y cell apart
  sz <- fa$sizes[fa$sizes$item == "I3" & fa$sizes$term == "g1:g2", ]
  expect_gt(nrow(sz), 0)
  by_pairs <- sz[sz$level_a == "b:y" | sz$level_b == "b:y", ]
  other_pairs <- sz[!(sz$level_a == "b:y" | sz$level_b == "b:y"), ]
  expect_gt(min(abs(by_pairs$difference)), max(abs(other_pairs$difference)))
  expect_lt(abs(max(abs(by_pairs$difference)) - 1.1), 0.4)
  expect_true(any(by_pairs$practical))
  # clean items produce no size rows
  expect_false(any(fa$sizes$item == "I5"))
})

test_that("dif_size guards: thin levels dropped, unknown factor errors", {
  s <- sim_dif(shifts = list(a = rep(0, 8), b = rep(0, 8)))
  g3 <- as.character(s$g); g3[1:5] <- "tiny"
  fit <- rasch(data.frame(s$X, grp = factor(g3)), factors = "grp")
  ds <- dif_size(fit, "I1", by = "grp", min_n = 20)
  expect_true(any(grepl("tiny", ds$notes)))
  expect_equal(nrow(ds$levels), 2)
  expect_error(dif_size(fit, "I1", by = "nofactor"), "factor")
})

test_that("person_extrapolated continues the score table geometrically", {
  set.seed(23)
  d <- seq(-2, 2, length.out = 10)
  X <- matrix(rbinom(400 * 10, 1, plogis(outer(rnorm(400, 0, 2.2), d, "-"))),
              400, 10)
  colnames(X) <- paste0("I", 1:10)
  fit <- rasch(X)
  skip_if(!any(fit$person$extreme), "no extreme persons in this draw")
  pe <- person_extrapolated(fit)
  st <- score_table(fit, extremes = "extrapolated")
  top <- pe$extreme & pe$raw == pe$max_raw & pe$n_items == 10
  bot <- pe$extreme & pe$raw == 0 & pe$n_items == 10
  if (any(top)) {
    expect_equal(unique(pe$theta_extrapolated[top]),
                 st$theta[nrow(st)], tolerance = 1e-8)
    # the geometric continuation differs from Warm's own extreme value
    # (Warm's finite extreme step exceeds the continued interior steps)
    expect_false(isTRUE(all.equal(unique(pe$theta_extrapolated[top]),
                                  unique(pe$theta[top]))))
    # still beyond the last interior score's measure
    expect_true(all(pe$theta_extrapolated[top] >
                    st$theta[nrow(st) - 1]))
  }
  if (any(bot)) {
    expect_equal(unique(pe$theta_extrapolated[bot]), st$theta[1],
                 tolerance = 1e-8)
    expect_true(all(pe$theta_extrapolated[bot] < st$theta[2]))
  }
  # non-extreme persons unchanged
  ne <- !pe$extreme
  expect_equal(pe$theta_extrapolated[ne], pe$theta[ne])
  expect_equal(pe$se_extrapolated[ne], pe$se[ne])
})

test_that("MFRM facet fit reports margin and pooled statistics with df", {
  set.seed(31)
  simP <- function(th, tau) { x <- 0:length(tau)
    p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%03d", 1:150); raters <- paste0("R", 1:3)
  th <- setNames(rnorm(150, 0, 1.2), persons)
  rho <- setNames(c(-0.5, 0, 0.5), raters)
  tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
  dd <- expand.grid(person = persons, item = names(tau), rater = raters,
                    stringsAsFactors = FALSE)
  dd$score <- mapply(function(p, i, r)
    sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho[r])),
    dd$person, dd$item, dd$rater)
  mf <- rasch_mfrm(dd, person = "person", item = "item", score = "score",
                   facets = "rater")
  fe <- mf$facet_effects$rater
  expect_true(all(c("fit_resid", "fit_resid_pooled", "df_fit") %in% names(fe)))
  expect_true(all(is.finite(fe$fit_resid)))
  # the margin statistic: the mean of the level's virtual items' fits
  vsel <- mf$virtual_map$rater == fe$level[1]
  expect_equal(fe$fit_resid[1],
               mean(mf$items$fit_resid[vsel], na.rm = TRUE),
               tolerance = 1e-10)
  # pooled df equals the cell df factor times the level's cells
  expect_equal(fe$df_fit, mf$summary_stats$df_factor * fe$n, tolerance = 1e-8)
  expect_true(all(is.finite(mf$item_effects$df_fit)))
  expect_true(all(is.finite(mf$item_effects$fit_resid_pooled)))
})

test_that("the factorial summary pivots to uniform/non-uniform per group term", {
  set.seed(1); n <- 600
  g1 <- rep(c("m", "f"), n / 2); g2 <- sample(c("young", "old"), n, TRUE)
  d <- seq(-1.5, 1.5, length.out = 6)
  sh <- matrix(0, n, 6); sh[g1 == "f", 3] <- 0.8
  X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, sex = g1, age = g2), factors = c("sex", "age"))
  fa <- dif_anova_factorial(fit)
  s <- fa$summary
  # one row per item and group term (no ci terms, no residual row)
  expect_setequal(unique(s$term), c("sex", "age", "sex:age"))
  expect_equal(nrow(s), 6 * 3)
  # the pivot agrees with the full table
  u <- fa$terms[fa$terms$item == "I3" & fa$terms$term == "sex", ]
  nu <- fa$terms[fa$terms$item == "I3" & fa$terms$term == "sex:ci", ]
  r <- s[s$item == "I3" & s$term == "sex", ]
  expect_equal(r$F_uniform, u$F_value)
  expect_equal(r$p_uniform_adj, u$p_adj)
  expect_equal(r$F_nonuniform, nu$F_value)
  # the planted uniform DIF is flagged as uniform, not non-uniform
  expect_true(r$uniform_DIF)
  # no misfit flag on the items table any more
  expect_false("misfit" %in% names(fit$items))
})

test_that("DIF class intervals adapt to the cells each analysis uses", {
  set.seed(1); n <- 800
  g1 <- rep(c("m", "f"), n / 2)
  g2 <- sample(c("a", "b", "c", "d"), n, TRUE)
  d <- seq(-1.5, 1.5, length.out = 8)
  X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-"))), n, 8)
  colnames(X) <- paste0("I", 1:8)
  f <- rasch(data.frame(X, sex = g1, age = g2), factors = c("sex", "age"))
  da <- dif_anova(f)
  ng <- attr(da, "n_groups")
  # per factor: smallest level / 30, clamped to 2..10, independent of overall
  expect_equal(unname(ng["sex"]),
               max(2L, min(10L, min(table(g1[!is.na(f$person$theta)])) %/% 30L)))
  expect_equal(unname(ng["age"]),
               max(2L, min(10L, min(table(g2[!is.na(f$person$theta)])) %/% 30L)))
  expect_gt(ng["sex"], ng["age"])   # fewer levels leave bigger cells
  # factorial: set from the smallest factor-combination cell
  fa <- dif_anova_factorial(f)
  cells <- interaction(g1, g2, drop = TRUE)
  expect_equal(fa$n_groups,
               max(2L, min(10L, min(table(cells[!is.na(f$person$theta)])) %/% 30L)))
  expect_lt(fa$n_groups, ng["sex"])
  # explicit n_groups still overrides
  da5 <- dif_anova(f, n_groups = 5)
  expect_true(all(attr(da5, "n_groups") == 5))
})
