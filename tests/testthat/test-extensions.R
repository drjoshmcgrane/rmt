simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("anchored estimation holds anchors and recovers the uncentred scale", {
  set.seed(7); Np <- 1200; L <- 10
  dtrue <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1] + 0.7
  th <- rnorm(Np, 0.7, 1.3)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, dtrue, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)

  fit <- rasch(X, anchors = data.frame(item = c("I01", "I10"), k = 1,
                                       tau = dtrue[c(1, 10)]))
  est <- fit$thresholds$tau
  expect_identical(est[c(1, 10)], dtrue[c(1, 10)])
  expect_identical(fit$thresholds$se[c(1, 10)], c(0, 0))
  expect_lt(sqrt(mean((est[-c(1, 10)] - dtrue[-c(1, 10)])^2)), 0.15)
  # person measures land on the anchored (uncentred) metric
  expect_equal(mean(fit$person$theta, na.rm = TRUE), 0.7, tolerance = 0.15)
  expect_error(rasch(X, model = "RSM",
                     anchors = data.frame(item = "I01", k = 1, tau = 0)),
               "PCM only")
  expect_error(rasch(X, anchors = data.frame(item = "NOPE", k = 1, tau = 0)))
})

test_that("MFRM recovers facet severities and item locations", {
  set.seed(1); Np <- 250
  persons <- sprintf("P%03d", seq_len(Np)); raters <- paste0("R", 1:5)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  rho_true <- setNames(c(-0.8, -0.3, 0, 0.3, 0.8), raters)
  tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4), D = c(0, 0.8))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  seen <- unlist(lapply(persons, function(p) paste(p, sample(raters, 3))))
  d <- d[paste(d$person, d$rater) %in% seen, ]   # incomplete rating design
  d$score <- mapply(function(p, i, r)
    sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho_true[r])),
    d$person, d$item, d$rater)

  fit <- rasch_mfrm(d, "person", "item", "score", facets = "rater")
  expect_s3_class(fit, "rasch_mfrm"); expect_s3_class(fit, "rasch")
  expect_true(fit$est$converged)

  fe <- fit$facet_effects$rater
  expect_gt(cor(fe$severity, rho_true[fe$level]), 0.98)
  expect_equal(sum(fe$severity), 0, tolerance = 1e-6)
  expect_true(all(fe$se > 0 & fe$se < 0.3))

  loc_true <- sapply(tau, mean) - mean(sapply(tau, mean))
  expect_gt(cor(fit$item_effects$location, loc_true[fit$item_effects$item]), 0.97)
  expect_gt(cor(fit$person$theta, th[fit$person$id], use = "complete.obs"), 0.85)
  # the full diagnostic object works at the virtual-item level
  expect_equal(ncol(fit$residuals), 20)
  expect_false(is.na(fit$psi$PSI))
})

test_that("MFRM flags an erratic rater through pooled fit", {
  set.seed(3); Np <- 300
  persons <- sprintf("P%03d", seq_len(Np)); raters <- paste0("R", 1:4)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  tau <- list(A = c(-1, 0, 1), B = c(-0.5, 0.3, 1.0), C = c(-1.1, 0, 0.8))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  d$score <- mapply(function(p, i, r) {
    if (r == "R4" && runif(1) < 0.35) return(sample(0:3, 1))   # erratic
    sample(0:3, 1, prob = simP(th[p], tau[[i]]))
  }, d$person, d$item, d$rater)
  fe <- rasch_mfrm(d, "person", "item", "score", facets = "rater")$facet_effects$rater
  expect_gt(fe$fit_resid[fe$level == "R4"], 2.5)          # erratic: strong underfit
  expect_true(all(fe$fit_resid[fe$level != "R4"] < 2.5))  # others show no underfit
})

test_that("the -1 missing-data code matches NA exactly", {
  set.seed(1); Np <- 400; L <- 8
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  miss <- sample(length(X), 0.05 * length(X))
  Xna <- X; Xna[miss] <- NA
  Xcode <- X; Xcode[miss] <- -1L

  f_na <- rasch(Xna); f_code <- rasch(Xcode)
  expect_equal(f_na$items$location, f_code$items$location)
  expect_equal(sum(is.na(f_code$X)), length(miss))
  expect_gte(min(f_code$X, na.rm = TRUE), 0)
  expect_true(any(grepl("missing-data code", f_code$notes)))

  # negative scores are always missing, even with na_codes opted out,
  # because valid category scores start at zero
  f_keep <- rasch(Xcode, na_codes = integer(0))
  expect_equal(f_keep$items$location, f_na$items$location)
  expect_equal(sum(is.na(f_keep$X)), length(miss))
})

test_that("MFRM honours the -1 missing code", {
  set.seed(2)
  simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%03d", 1:150); raters <- paste0("R", 1:4)
  th <- setNames(rnorm(150, 0, 1.2), persons)
  tau <- list(A = c(-0.5, 0.8), B = c(-1, 0.5), C = c(-0.2, 1))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  d$score <- mapply(function(p, i, r) sample(0:2, 1, prob = simP(th[p], tau[[i]])),
                    d$person, d$item, d$rater)
  d_na <- d; d_code <- d
  hit <- sample(nrow(d), 40)
  d_na$score[hit] <- NA; d_code$score[hit] <- -1L
  f_na <- rasch_mfrm(d_na, "person", "item", "score", facets = "rater")
  f_code <- rasch_mfrm(d_code, "person", "item", "score", facets = "rater")
  expect_equal(f_na$facet_effects$rater$severity,
               f_code$facet_effects$rater$severity, tolerance = 1e-8)
})

test_that("Guttman reproducibility is high for near-deterministic data", {
  set.seed(11); Np <- 300; L <- 10
  d <- seq(-3, 3, length.out = L)
  th <- rnorm(Np, 0, 2.5)
  # very spread persons and items -> near-deterministic Guttman pattern
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  g <- guttman_table(rasch(X))
  expect_true(g$CR >= 0 && g$CR <= 1)
  expect_gt(g$CR, 0.8)
  expect_equal(dim(g$matrix), c(Np, L))
  # items ordered easy to hard across columns
  expect_equal(colnames(g$matrix), rasch(X)$items$item[order(rasch(X)$items$location)])
})

test_that("subtests absorb local dependence", {
  set.seed(1); Np <- 1200; L <- 12
  dl <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np, 0, 1.4)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, dl, "-"))), Np, L)
  colnames(X) <- sprintf("U%02d", 1:L)
  X[, 5] <- ifelse(runif(Np) < 0.9, X[, 4], X[, 5])

  fit <- rasch(X)
  fl <- residual_correlations(fit)$flagged
  expect_true(any(fl$item_a == "U04" & fl$item_b == "U05"))

  fit2 <- combine_items(fit, list(c("U04", "U05")))
  expect_true("U04+U05" %in% fit2$items$item)
  expect_equal(ncol(fit2$X), L - 1)
  expect_equal(fit2$items$max[fit2$items$item == "U04+U05"], 2)
  fl2 <- residual_correlations(fit2)$flagged
  expect_false(any(grepl("U04", fl2$item_a) | grepl("U04", fl2$item_b)))
  expect_true(any(grepl("subtest formed", fit2$notes)))
  expect_error(combine_items(fit, list("U01")), "at least two")
  expect_error(combine_items(fit, list(c("U01", "ZZ"))), "not in the fit")
})

test_that("MFRM handles two facets jointly (rater x occasion)", {
  set.seed(21); Np <- 300
  persons <- sprintf("P%03d", seq_len(Np))
  raters <- paste0("R", 1:3); occs <- paste0("T", 1:2)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  rho_r <- setNames(c(-0.5, 0, 0.5), raters)
  rho_o <- setNames(c(-0.3, 0.3), occs)
  tau <- list(A = c(-1, 1), B = c(-0.4, 1.1), C = c(-1.1, 0.5))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   occasion = occs, stringsAsFactors = FALSE)
  d$score <- mapply(function(p, i, r, o)
    sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho_r[r] + rho_o[o])),
    d$person, d$item, d$rater, d$occasion)

  fit <- rasch_mfrm(d, "person", "item", "score",
                    facets = c("rater", "occasion"))
  expect_true(fit$est$converged)
  expect_setequal(names(fit$facet_effects), c("rater", "occasion"))
  fr <- fit$facet_effects$rater; fo <- fit$facet_effects$occasion
  expect_lt(max(abs(fr$severity - (rho_r - mean(rho_r)))), 0.15)
  expect_lt(max(abs(fo$severity - (rho_o - mean(rho_o)))), 0.15)
  expect_equal(sum(fr$severity), 0, tolerance = 1e-8)
  expect_equal(sum(fo$severity), 0, tolerance = 1e-8)
  # virtual items = item x rater x occasion cells
  expect_equal(ncol(fit$X), 3 * 3 * 2)
})

test_that("wide-format MFRM entry matches the long form exactly", {
  set.seed(5)
  persons <- sprintf("P%03d", 1:120); raters <- c("R1", "R2", "R3")
  items <- sprintf("I%d", 1:5)
  th <- setNames(rnorm(120), persons)
  sev <- c(R1 = -0.4, R2 = 0, R3 = 0.4)
  d <- expand.grid(person = persons, rater = raters, item = items,
                   stringsAsFactors = FALSE)
  dl <- setNames(seq(-1, 1, length.out = 5), items)
  d$score <- rbinom(nrow(d), 1, plogis(th[d$person] - dl[d$item] - sev[d$rater]))
  w <- reshape(d, idvar = c("person", "rater"), timevar = "item",
               direction = "wide")
  names(w) <- sub("^score\\.", "", names(w))
  fl <- rasch_mfrm(d, "person", "item", "score", facets = "rater")
  fw <- rasch_mfrm(w, "person", facets = "rater", items = items)
  expect_equal(fw$items$location, fl$items$location, tolerance = 1e-10)
  expect_equal(fw$facet_effects$rater$severity,
               fl$facet_effects$rater$severity, tolerance = 1e-10)
  # the interactive structure works through the wide entry too
  fwi <- rasch_mfrm(w, "person", facets = "rater", items = items,
                    interaction = "rater")
  expect_false(is.null(fwi$interaction_effects))
  # guardrails
  expect_error(rasch_mfrm(w, "person", item = "x", score = "y",
                          facets = "rater", items = items), "either")
  expect_error(rasch_mfrm(w, "person", facets = "rater"), "either")
})
