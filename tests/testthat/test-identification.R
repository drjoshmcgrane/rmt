# Structural identification under missing data: the item-pair connectivity
# requirement (booklet designs), and honest reporting for thresholds whose
# adjacent categories are nearly empty. Both paths were found live by
# probing booklet designs: a disconnected design used to fit silently with
# arbitrary between-block spacing, and a category with one response drove
# its threshold to +12 logits with a reported SE of 0.

test_that("a linked booklet design with zero complete cases recovers", {
  set.seed(11)
  L <- 24; btrue <- seq(-2, 2, length.out = L); N <- 600
  books <- list(1:11, 8:17, 14:24)              # overlaps link the chain
  bk <- sample(1:3, N, TRUE); th <- rnorm(N)
  X <- matrix(NA_integer_, N, L, dimnames = list(NULL, sprintf("I%02d", 1:L)))
  for (i in 1:N) { it <- books[[bk[i]]]
    X[i, it] <- rbinom(length(it), 1, plogis(th[i] - btrue[it])) }
  expect_equal(sum(rowSums(is.na(X)) == 0), 0L)  # no person saw every item
  f <- rasch(as.data.frame(X))
  bt <- btrue - mean(btrue)
  est <- f$items$location[match(colnames(X), f$items$item)]
  se <- f$items$se[match(colnames(X), f$items$item)]
  expect_gt(cor(est, bt), 0.95)
  expect_lt(abs(mean(est - bt)), 0.1)
  expect_gte(sum(abs(est - bt) <= 1.96 * se), L - 3L)   # ~nominal coverage
})

test_that("a disconnected design errors informatively", {
  set.seed(2)
  X <- matrix(rbinom(200 * 10, 1, 0.5), 200, 10,
              dimnames = list(NULL, paste0("I", 1:10)))
  X[1:100, 6:10] <- NA; X[101:200, 1:5] <- NA
  expect_error(rasch(as.data.frame(X)), "not connected.*I1, I2.*I6, I7")
  expect_error(pcml(X), "not connected")
  expect_error(pcml_pc(X, n_components = 1), "not connected")
})

test_that("anchors identify disconnected blocks (disjoint-form equating)", {
  set.seed(3)
  X <- matrix(rbinom(200 * 10, 1, 0.5), 200, 10,
              dimnames = list(NULL, paste0("I", 1:10)))
  X[1:100, 6:10] <- NA; X[101:200, 1:5] <- NA
  # an anchor in every block fixes every origin: this is legitimate
  f <- pcml(X, anchors = data.frame(item = c("I1", "I6"), k = 1,
                                    tau = c(-0.5, 0.5)))
  expect_true(f$converged)
  expect_equal(f$thr$tau[f$thr$item == 1], -0.5)
  # an anchor in only one block leaves the other unidentified
  expect_error(pcml(X, anchors = data.frame(item = "I1", k = 1, tau = -0.5)),
               "without an\\s+anchored item.*I6, I7")
})

test_that("a threshold beside a near-empty category reports NA SE + note", {
  set.seed(4)
  N <- 250; L <- 8; th <- rnorm(N)
  X <- matrix(NA_integer_, N, L, dimnames = list(NULL, paste0("I", 1:L)))
  for (j in 1:L) {
    d <- th - (j - 4.5) / 2
    p1 <- plogis(d); p2 <- plogis(d - 0.8)
    X[, j] <- rbinom(N, 1, p1) + rbinom(N, 1, p2 * 0.8)
  }
  X[, "I8"] <- pmin(X[, "I8"], 1L)               # I8 top category empty so far
  X[1, "I8"] <- 2L                               # exactly one response in cat 2
  f <- rasch(as.data.frame(X))
  expect_true(any(grepl("I8.*only 1 response.*category 2.*weakly determined",
                        f$notes)))
  expect_true(is.na(f$items$se[f$items$item == "I8"]))
  expect_true(all(is.finite(f$items$se[f$items$item != "I8"])))
  # the threshold-level flag and NA sit exactly on the weak threshold
  thr <- f$est$thr
  expect_true(thr$weak[thr$item == 8L & thr$k == 2L])
  expect_true(is.na(thr$se[thr$item == 8L & thr$k == 2L]))
  expect_false(any(thr$weak[thr$item != 8L]))
})
