# Cross-validation against independent implementations. Three external
# referees, three levels of expected agreement:
# - sirt::rasch.pairwise estimates the SAME pairwise conditional family
#   (Choppin's explicit method): near-identity expected.
# - eRm fits the FULL conditional likelihood (Andersen CML): a different
#   consistent estimator of the same parameters, so agreement to O(1/sqrt(N))
#   with a small pairwise efficiency cost in the standard errors.
# - psychotools::btmodel maximises the same Bradley-Terry likelihood as
#   btl(): agreement to machine precision.

test_that("dichotomous locations match eRm's full CML and sirt's pairwise", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("sirt")
  set.seed(42)
  N <- 800; L <- 15; btrue <- seq(-2, 2, length.out = L)
  th <- rnorm(N)
  X <- matrix(rbinom(N * L, 1, plogis(outer(th, btrue, "-"))), N, L,
              dimnames = list(NULL, sprintf("I%02d", 1:L)))
  f <- rasch(as.data.frame(X))
  ours <- f$items$location[match(colnames(X), f$items$item)]
  se_ours <- f$items$se[match(colnames(X), f$items$item)]

  e <- eRm::RM(X)
  erm_loc <- -e$betapar; erm_loc <- erm_loc - mean(erm_loc)
  expect_gt(cor(ours, erm_loc), 0.999)
  expect_lt(max(abs(ours - erm_loc)), 0.08)
  # pairwise trades a little efficiency for its robustness: SEs at or just
  # above CML, never materially below
  ratio <- se_ours / e$se.beta
  expect_true(all(ratio > 0.97 & ratio < 1.15))

  s <- sirt::rasch.pairwise(X)
  sirt_loc <- s$item$b - mean(s$item$b)
  expect_gt(cor(ours, sirt_loc), 0.9999)
  expect_lt(max(abs(ours - sirt_loc)), 0.02)
})

test_that("PCM and RSM thresholds match eRm's CML decomposition", {
  skip_if_not_installed("eRm")
  set.seed(7)
  N <- 800; L <- 10; btrue <- seq(-1.5, 1.5, length.out = L)
  th <- rnorm(N); m <- 3
  tau <- t(apply(matrix(rnorm(L * m, 0, .8), L, m), 1, sort))
  X <- matrix(0L, N, L, dimnames = list(NULL, sprintf("I%02d", 1:L)))
  for (j in 1:L) for (i in 1:N) {
    d <- th[i] - btrue[j] - tau[j, ]
    p <- c(1, exp(cumsum(d))); X[i, j] <- sample(0:m, 1, prob = p / sum(p))
  }
  f <- rasch(as.data.frame(X))
  tp <- eRm::thresholds(eRm::PCM(X))$threshtable[[1]]
  erm_thr <- as.vector(t(tp[, -1])); erm_thr <- erm_thr - mean(erm_thr)
  our_thr <- f$est$thr$tau - mean(f$est$thr$tau)
  expect_gt(cor(our_thr, erm_thr), 0.999)
  expect_lt(max(abs(our_thr - erm_thr)), 0.15)

  fr <- rasch(as.data.frame(X), model = "RSM")
  tr <- eRm::thresholds(eRm::RSM(X))$threshtable[[1]]
  erm_r <- as.vector(t(tr[, -1])); erm_r <- erm_r - mean(erm_r)
  our_r <- fr$est$thr$tau - mean(fr$est$thr$tau)
  expect_gt(cor(our_r, erm_r), 0.999)
  expect_lt(max(abs(our_r - erm_r)), 0.08)
})

test_that("btl() equals psychotools::btmodel to machine precision", {
  skip_if_not_installed("psychotools")
  data("GermanParties2009", package = "psychotools")
  pref <- GermanParties2009$preference
  bm <- psychotools::btmodel(pref)
  w <- psychotools::itempar(bm)
  m2 <- as.matrix(pref)
  prs <- strsplit(colnames(m2), ":", fixed = TRUE)
  rows <- lapply(seq_along(prs), function(j) {
    a <- prs[[j]][1]; b <- prs[[j]][2]
    data.frame(object_a = a, object_b = b,
               winner = ifelse(m2[, j] == 1, a, b))
  })
  d <- do.call(rbind, rows)
  bt <- btl(d, "object_a", "object_b", "winner")
  pt <- log(as.numeric(w))[match(bt$objects$object, names(w))]
  pt <- pt - mean(pt)
  expect_lt(max(abs(bt$objects$location - pt)), 1e-8)
})
