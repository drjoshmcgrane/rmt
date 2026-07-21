# Composite-likelihood information criteria in compare_fits(): the penalty
# tr(H^-1 J) absorbs the pairwise over-counting, so CL-AIC/CL-BIC (Varin &
# Vidoni 2005; Gao & Song 2010) must select the planted structure where the
# nominal-count criteria could not.

gen_pcm <- function(N, btrue, tau_mat, seed) {
  set.seed(seed)
  L <- length(btrue); m <- ncol(tau_mat); th <- rnorm(N)
  X <- matrix(0L, N, L, dimnames = list(NULL, sprintf("I%02d", 1:L)))
  for (j in 1:L) for (i in 1:N) {
    d <- th[i] - btrue[j] - tau_mat[j, ]
    p <- c(1, exp(cumsum(d)))
    X[i, j] <- sample(0:m, 1, prob = p / sum(p))
  }
  X
}

test_that("CL-ICs select RSM when the rating structure is true", {
  L <- 8; btrue <- seq(-1.2, 1.2, length.out = L)
  tau_rsm <- matrix(rep(c(-0.9, 0, 0.9), each = L), L, 3)
  X <- gen_pcm(500, btrue, tau_rsm, 21)
  cmp <- compare_fits(PCM = rasch(X), RSM = rasch(X, model = "RSM"))
  expect_equal(cmp$label[which.min(cmp$cl_bic)], "RSM")
  expect_equal(cmp$label[which.min(cmp$cl_aic)], "RSM")
  # the effective count sits well above the nominal one: each response
  # enters every pair its item forms
  expect_true(all(cmp$eff_params > cmp$parameters))
  expect_true(all(cmp$eff_params < cmp$parameters * (L - 1)))
})

test_that("CL-ICs select PCM when the threshold spreads truly vary", {
  L <- 8; btrue <- seq(-1.2, 1.2, length.out = L)
  set.seed(5)
  tau_pcm <- t(sapply(1:L, function(j)
    sort(rnorm(3, 0, 1)) * runif(1, 0.4, 1.8)))
  X <- gen_pcm(500, btrue, tau_pcm, 22)
  cmp <- compare_fits(PCM = rasch(X), RSM = rasch(X, model = "RSM"))
  expect_equal(cmp$label[which.min(cmp$cl_bic)], "PCM")
  expect_equal(cmp$label[which.min(cmp$cl_aic)], "PCM")
})

sim_btl_pos <- function(pos_effect, seed) {
  set.seed(seed)
  K <- 8; beta <- seq(-1.5, 1.5, length.out = K)
  names(beta) <- paste0("O", K:1)
  n <- 1200
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  jd <- sample(sprintf("J%02d", 1:15), n, TRUE)
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib] + pos_effect))
  data.frame(object_a = names(beta)[ia], object_b = names(beta)[ib],
             winner = ifelse(win == 1, names(beta)[ia], names(beta)[ib]),
             judge = jd)
}

test_that("CL-ICs detect a planted BTL position effect and reject an absent one", {
  d1 <- sim_btl_pos(0.5, 35)
  b0 <- btl(d1, "object_a", "object_b", "winner", judge = "judge")
  b1 <- btl(d1, "object_a", "object_b", "winner", judge = "judge",
            position = TRUE)
  cmp <- compare_fits(plain = b0, position = b1)
  expect_equal(cmp$label[which.min(cmp$cl_bic)], "position")
  expect_true(all(c("judges", "objects", "OSI") %in% names(cmp)))
  expect_true(cmp$same_data[2])
  d0 <- sim_btl_pos(0, 30)
  c0 <- btl(d0, "object_a", "object_b", "winner", judge = "judge")
  c1 <- btl(d0, "object_a", "object_b", "winner", judge = "judge",
            position = TRUE)
  cmp0 <- compare_fits(plain = c0, position = c1)
  expect_equal(cmp0$label[which.min(cmp0$cl_bic)], "plain")
})

test_that("mixtures are refused; MFRM fits get NA ICs with the reason noted", {
  X <- gen_pcm(150, seq(-1, 1, length.out = 5),
               matrix(rep(c(-0.5, 0.5), each = 5), 5, 2), 9)
  f <- rasch(X)
  d <- sim_btl_pos(0, 11)
  b <- btl(d, "object_a", "object_b", "winner")
  expect_error(compare_fits(f, b), "not a mixture")
  # MFRM: no Godambe matrices on the assembled est -> NA ICs, note says why
  long <- data.frame(person = rep(sprintf("P%03d", 1:120), each = 4),
                     rater = rep(c("R1", "R2"), 240),
                     item = rep(rep(c("A", "B"), each = 2), 120))
  set.seed(13)
  long$score <- rbinom(nrow(long), 2, 0.5)
  mf <- rasch_mfrm(long, person = "person", item = "item", score = "score",
                   facets = "rater")
  cmp <- compare_fits(a = mf, b = mf)
  expect_true(all(is.na(cmp$cl_aic)))
  expect_match(attr(cmp, "note"), "MFRM/EFRM")
})

test_that("compare_fits withholds ICs for an unconverged fit", {
  set.seed(50); N <- 300; L <- 5
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), d, "-"))), N, L)
  colnames(X) <- paste0("I", 1:L)
  good <- rasch(as.data.frame(X))
  bad <- rasch(as.data.frame(X), maxit = 1L)     # forced non-convergence
  expect_false(isTRUE(bad$est$converged))
  cmp <- compare_fits(good = good, bad = bad)
  expect_true("converged" %in% names(cmp))
  expect_true(is.na(cmp$cl_aic[cmp$label == "bad"]))
  expect_true(is.na(cmp$loglik[cmp$label == "bad"]))
})

test_that("compare_fits same_data compares actual responses", {
  set.seed(51); N <- 300; L <- 5
  d <- seq(-1.5, 1.5, length.out = L)
  mk <- function(s) { set.seed(s)
    X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), d, "-"))), N, L)
    colnames(X) <- paste0("I", 1:L); rasch(as.data.frame(X)) }
  f1 <- mk(1); f2 <- mk(2)      # same items/max/N, different responses
  cmp <- compare_fits(a = f1, b = f2)
  expect_false(cmp$same_data[cmp$label == "b"])
  expect_true(is.na(cmp$two_delta_ll[cmp$label == "b"]))
})
