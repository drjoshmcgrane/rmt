simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("average anchoring fixes item means with free thresholds", {
  set.seed(9); Np <- 1200
  loc_true <- seq(-1.6, 1.6, length.out = 8) + 0.5
  tt <- lapply(seq_len(8), function(i) loc_true[i] + c(-0.9, 0, 0.9))
  th <- rnorm(Np, 0.5, 1.3)
  X <- sapply(seq_len(8), function(i)
    sapply(th, function(t) sample(0:3, 1, prob = simP(t, tt[[i]]))))
  colnames(X) <- sprintf("P%d", 1:8)

  fit <- rasch(X, anchors = data.frame(item = c("P1", "P8"), k = NA,
                                       tau = loc_true[c(1, 8)]))
  expect_equal(fit$items$location[c(1, 8)], loc_true[c(1, 8)],
               tolerance = 1e-9, ignore_attr = TRUE)
  expect_equal(fit$items$se[c(1, 8)], c(0, 0), tolerance = 1e-6)
  expect_true(all(fit$thresholds$se[fit$thresholds$item == 1] > 0))
  expect_lt(sqrt(mean((fit$items$location[2:7] - loc_true[2:7])^2)), 0.15)

  # mixed: mean anchor and threshold anchor together
  fit2 <- rasch(X, anchors = data.frame(item = c("P1", "P8"), k = c(NA, 2),
                                        tau = c(loc_true[1], tt[[8]][2])))
  expect_equal(mean(fit2$tau_list[[1]]), loc_true[1], tolerance = 1e-9)
  expect_identical(fit2$tau_list[[8]][2], tt[[8]][2])
  expect_error(rasch(X, anchors = data.frame(item = "P1", k = c(NA, 1),
                                             tau = c(0, 0))),
               "both an average anchor and threshold anchors")
})

test_that("split_items resolves planted uniform DIF", {
  set.seed(4); Np <- 1500; L <- 10
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  grp <- rep(c("ref", "foc"), each = Np / 2); th <- rnorm(Np, 0, 1.4)
  shift <- matrix(0, Np, L); shift[grp == "foc", 3] <- 1.0
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-") - shift)), Np, L)
  colnames(X) <- sprintf("G%02d", 1:L)
  fit <- rasch(data.frame(X, grp = grp), factors = "grp", n_groups = 6)
  expect_true(dif_anova(fit)$uniform_DIF[3])

  fit2 <- split_items(fit, "G03", by = "grp")
  locs <- setNames(fit2$items$location, fit2$items$item)
  expect_true(all(c("G03 (foc)", "G03 (ref)") %in% names(locs)))
  expect_equal(unname(locs["G03 (foc)"] - locs["G03 (ref)"]), 1.0, tolerance = 0.3)
  expect_false(any(dif_anova(fit2)$uniform_DIF, na.rm = TRUE))
  expect_true(any(grepl("split", fit2$notes)))
  expect_error(split_items(fit, "G03", by = "nope"), "not a person factor")
})

test_that("equate_tests flags drifted common items only", {
  set.seed(5); L <- 12; d <- seq(-2, 2, length.out = L)
  mk <- function(drift = 0) {
    dd <- d; dd[4] <- dd[4] + drift
    X <- matrix(rbinom(900 * L, 1, plogis(outer(rnorm(900), dd, "-"))), 900, L)
    colnames(X) <- sprintf("I%02d", 1:L); rasch(X)
  }
  f1 <- mk()
  eq0 <- equate_tests(f1, mk())
  expect_equal(sum(eq0$table$drift), 0)
  expect_gt(eq0$correlation, 0.99)

  eq1 <- equate_tests(f1, mk(drift = 0.8))
  expect_identical(eq1$table$item[eq1$table$drift], "I04")

  # reference table path (item bank style)
  bank <- data.frame(item = sprintf("I%02d", 1:L), location = d - mean(d), se = 0.05)
  eqb <- equate_tests(f1, bank)
  expect_equal(eqb$n, L)
  expect_error(equate_tests(f1, data.frame(item = "ZZ", location = 0)),
               "at least two common items")
})

test_that("interactive facet mode recovers a planted item-by-rater effect", {
  set.seed(6); Np <- 500
  persons <- sprintf("P%03d", 1:Np); raters <- paste0("R", 1:4)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4), D = c(0, 0.8))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  d$score <- mapply(function(p, i, r) {
    extra <- if (i == "B" && r == "R2") 0.9 else 0
    sample(0:2, 1, prob = simP(th[p], tau[[i]] + extra))
  }, d$person, d$item, d$rater)

  fit <- rasch_mfrm(d, "person", "item", "score", facets = "rater",
                    interaction = "rater")
  expect_true(fit$est$converged)
  ie <- fit$interaction_effects
  expect_equal(nrow(ie), 16)
  # double sum-to-zero margins
  expect_equal(max(abs(tapply(ie$gamma, ie$item, sum))), 0, tolerance = 1e-8)
  expect_equal(max(abs(tapply(ie$gamma, ie$level, sum))), 0, tolerance = 1e-8)
  top <- ie[which.max(abs(ie$gamma)), ]
  expect_identical(top$item, "B"); expect_identical(top$level, "R2")
  expect_equal(top$gamma, 0.9 * (1 - 1/4 - 1/4 + 1/16), tolerance = 0.2)
  expect_error(rasch_mfrm(d, "person", "item", "score", facets = "rater",
                          interaction = "nope"), "must name one of the facets")
})

test_that("factorial DIF detects the planted factor with Tukey output", {
  set.seed(1); n <- 1200
  d <- seq(-1.5, 1.5, length.out = 8)
  g1 <- rep(c("a", "b"), each = n / 2)
  g2 <- sample(c("x", "y"), n, replace = TRUE)
  th <- rnorm(n)
  X <- sapply(seq_along(d), function(i)
    rbinom(n, 1, plogis(th - d[i] - if (i == 3) ifelse(g1 == "b", 1, 0) else 0)))
  colnames(X) <- paste0("I", 1:8)
  fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))

  df <- dif_anova_factorial(fit)
  t3 <- df$terms[df$terms$item == "I3", ]
  expect_true(t3$significant[t3$term == "g1"])
  expect_false(t3$significant[t3$term == "g2"])
  expect_false(t3$superseded[t3$term == "g1"])
  tk3 <- df$tukey[df$tukey$item == "I3" & df$tukey$term == "g1", ]
  expect_equal(nrow(tk3), 1)
  expect_lt(tk3$p_tukey, 0.001)
  # per-factor analysis with BH agrees on the planted item
  da <- dif_anova(fit)
  expect_true(da$uniform_DIF[da$factor == "g1" & da$item == "I3"])
  expect_true("p_uniform_adj" %in% names(da))
})

test_that("a significant interaction supersedes its main effects", {
  set.seed(2); n <- 1600
  d <- seq(-1, 1, length.out = 6)
  g1 <- rep(c("a", "b"), each = n / 2)
  g2 <- rep(c("x", "y"), times = n / 2)
  th <- rnorm(n)
  # DIF only in the (b, y) cell: pure g1:g2 interaction structure
  sh <- ifelse(g1 == "b" & g2 == "y", 1.2, 0)
  X <- sapply(seq_along(d), function(i)
    rbinom(n, 1, plogis(th - d[i] - if (i == 2) sh else 0)))
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
  df <- dif_anova_factorial(fit)
  t2 <- df$terms[df$terms$item == "I2", ]
  expect_true(t2$significant[t2$term == "g1:g2"])
  # the cell-shift induces main effects too; they must be marked superseded
  for (tt in c("g1", "g2"))
    if (t2$significant[t2$term == tt]) expect_true(t2$superseded[t2$term == tt])
  expect_false(t2$superseded[t2$term == "g1:g2"])
  # Tukey on the interaction compares the four cells
  tki <- df$tukey[df$tukey$item == "I2" & df$tukey$term == "g1:g2", ]
  expect_equal(nrow(tki), 6)   # choose(4, 2) cell contrasts
})

test_that("multiple-choice scoring and miskey detection work", {
  set.seed(5); Np <- 800
  th <- rnorm(Np)
  d <- seq(-1.2, 1.2, length.out = 8)
  keyv <- setNames(rep("A", 8), sprintf("M%02d", 1:8))
  keyv["M04"] <- "B"                                  # miskey: true correct is C
  raw <- sapply(seq_along(d), function(i) {
    correct <- if (i == 4) "C" else "A"
    ok <- rbinom(Np, 1, plogis(th - d[i]))
    ifelse(ok == 1, correct,
           sample(setdiff(c("A", "B", "C", "D"), correct), Np, replace = TRUE))
  })
  colnames(raw) <- sprintf("M%02d", 1:8)
  raw[sample(length(raw), 60)] <- ""                  # blanks become missing

  fit <- rasch(raw, key = keyv)
  expect_false(is.null(fit$mc))
  expect_true(all(fit$m == 1))
  expect_true(any(grepl("scored 0/1 against the key", fit$notes)))
  # scoring matches a manual comparison (blanks NA)
  manual <- ifelse(raw[, "M01"] == "", NA_integer_,
                   as.integer(raw[, "M01"] == "A"))
  expect_identical(unname(fit$X[, "M01"]), manual)

  da <- distractor_analysis(fit)
  m4 <- da[da$item == "M04", ]
  expect_true(m4$flag[m4$option == "C"])              # the real correct option
  expect_false(any(m4$flag[m4$option != "C"]))
  expect_equal(sum(da$flag[da$item != "M04"]), 0)     # clean items stay clean
  # the keyed option carries the top point-biserial on clean items
  clean <- da[da$item == "M01", ]
  expect_identical(clean$option[which.max(clean$point_biserial)], "A")

  # key as a data frame, case-insensitive matching
  fit2 <- rasch(raw, key = data.frame(item = names(keyv), key = tolower(keyv)))
  expect_equal(fit2$items$location, fit$items$location)
  expect_error(distractor_analysis(rasch(fit$X)), "no key")
})

test_that("dimensionality: 10-component PCA, scree, manual subsets, exact CI", {
  set.seed(2); Np <- 1200; L <- 16
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  thA <- rnorm(Np, 0, 1.4); thB <- 0.3 * thA + sqrt(1 - 0.3^2) * rnorm(Np, 0, 1.4)
  XA <- matrix(rbinom(Np * 8, 1, plogis(outer(thA, d[1:8], "-"))), Np, 8)
  XB <- matrix(rbinom(Np * 8, 1, plogis(outer(thB, d[9:16], "-"))), Np, 8)
  X <- cbind(XA, XB); colnames(X) <- sprintf("D%02d", 1:16)
  fit <- rasch(X)

  pc <- residual_pca(fit)
  expect_equal(ncol(pc$loadings_matrix), 11)        # item + PC1..PC10
  expect_equal(nrow(pc$eigen_table), 10)
  expect_true(all(diff(pc$eigen_table$eigenvalue) <= 1e-10))
  expect_equal(pc$eigen_table$cumulative[10],
               sum(pc$eigenvalues[1:10]) / sum(pc$eigenvalues))
  et <- plot_scree(fit)
  expect_equal(nrow(et), 10)

  # default split detects the planted second dimension; exact CI fields present
  dt <- dimensionality_test(fit)
  expect_true(dt$multidimensional)
  expect_identical(dt$split, "first residual contrast")
  expect_true(dt$ci[1] >= 0 && dt$ci[2] <= 1 && dt$ci[1] < dt$ci[2])
  expect_true(dt$n + dt$n_excluded_extreme >= dt$n)

  # manual subsets matching the true structure also detect it
  dtm <- dimensionality_test(fit, items_positive = sprintf("D%02d", 1:8),
                             items_negative = sprintf("D%02d", 9:16))
  expect_true(dtm$multidimensional)
  expect_identical(dtm$split, "manual")
  expect_gt(dtm$prop_significant, 0.05)

  expect_error(dimensionality_test(fit, items_positive = sprintf("D%02d", 1:8)),
               "both item subsets")
  expect_error(dimensionality_test(fit, items_positive = c("D01", "D02"),
                                   items_negative = c("D02", "D03")),
               "disjoint")
  expect_error(dimensionality_test(fit, items_positive = c("D01", "ZZ"),
                                   items_negative = c("D03", "D04")),
               "not in the fit")
})

test_that("compare_fits contrasts nested models on the same data", {
  set.seed(3); Np <- 700
  simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  th <- rnorm(Np)
  loc <- seq(-1, 1, length.out = 6); step <- c(-0.8, 0, 0.8)
  X <- sapply(loc, function(b) sapply(th, function(t)
    sample(0:3, 1, prob = simP(t, b + step))))
  colnames(X) <- paste0("R", 1:6)
  pcm <- rasch(X, model = "PCM"); rsm <- rasch(X, model = "RSM")

  cmp <- compare_fits(PCM = pcm, RSM = rsm)
  expect_s3_class(cmp, "rasch_compare")
  expect_identical(attr(cmp, "reference"), "PCM")
  expect_true(all(cmp$same_data))
  # RSM is nested in PCM: fewer parameters, lower (or equal) loglik
  expect_lt(cmp$parameters[2], cmp$parameters[1])
  expect_lte(cmp$two_delta_ll[2], 1e-8)
  expect_equal(cmp$delta_parameters[2],
               cmp$parameters[2] - cmp$parameters[1])
  # different data -> no loglik comparison, descriptive columns still there
  X2 <- X[, 1:5]
  cmp2 <- compare_fits(full = pcm, short = rasch(X2))
  expect_false(cmp2$same_data[2])
  expect_true(is.na(cmp2$two_delta_ll[2]))
  expect_true(all(is.finite(cmp2$chisq_per_df)))
  expect_error(compare_fits(pcm), "at least two")
})

test_that("maxit and tol are honoured by the estimators", {
  set.seed(4); Np <- 400; L <- 8
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  # a very loose tolerance stops after the first step
  f_loose <- rasch(X, maxit = 60, tol = 10)
  expect_lte(f_loose$est$iterations, 2)
  # tight settings converge as usual and agree with defaults
  f_tight <- rasch(X, maxit = 200, tol = 1e-10)
  f_def <- rasch(X)
  expect_equal(f_tight$items$location, f_def$items$location, tolerance = 1e-6)
})
