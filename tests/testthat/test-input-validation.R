# Input honesty: misspelled column names, fractional scores, missing
# identifiers, unusable equating items, and HTML escaping. Each of these
# used to fail silently (wrong-data analyses, truncated responses, phantom
# "NA" facet levels, all-NA equating, markup injection in reports).

mkX <- function(n = 150, L = 6, seed = 1) {
  set.seed(seed)
  X <- matrix(rbinom(n * L, 1, plogis(outer(rnorm(n), seq(-1, 1, length.out = L), "-"))),
              n, L, dimnames = list(NULL, paste0("I", 1:L)))
  X
}

test_that("misspelled id, factor, and item columns are errors, not fallbacks", {
  d <- as.data.frame(mkX())
  d$sex <- rep(c("m", "f"), length.out = nrow(d))
  d$pid <- sprintf("P%03d", seq_len(nrow(d)))
  expect_error(rasch(d, id = "person_id", items = paste0("I", 1:6)),
               "id column 'person_id' not found")
  expect_error(rasch(d, factors = c("sex", "agee"), items = paste0("I", 1:6)),
               "factor column\\(s\\) not found.*agee")
  expect_error(rasch(d, items = c("I1", "I2", "Item3")),
               "item column\\(s\\) not found.*Item3")
  # correct names still work, including numeric item indices
  f <- rasch(d, id = "pid", factors = "sex", items = paste0("I", 1:6))
  expect_equal(ncol(f$X), 6L)
  f2 <- rasch(as.data.frame(mkX()), items = 1:6)
  expect_equal(ncol(f2$X), 6L)
})

test_that("fractional scores error instead of silently truncating", {
  X <- mkX()
  Xf <- as.data.frame(X)
  Xf$I3[5] <- 1.9
  expect_error(rasch(Xf), "non-integer score\\(s\\) in: I3.*1\\.9")
  # integer-valued doubles ("2.0") are fine
  Xd <- as.data.frame(X * 1.0)
  expect_s3_class(rasch(Xd), "rasch")
})

test_that("MFRM rows with missing identifiers are dropped with a note", {
  set.seed(3)
  long <- data.frame(person = rep(sprintf("P%03d", 1:80), each = 4),
                     rater = rep(c("R1", "R2"), 160),
                     item = rep(rep(c("A", "B"), each = 2), 80),
                     score = rbinom(320, 2, 0.5))
  long$rater[c(5, 9)] <- NA
  long$person[17] <- NA
  f <- rasch_mfrm(long, person = "person", item = "item", score = "score",
                  facets = "rater")
  expect_true(any(grepl("3 row\\(s\\) dropped: missing person, item, or facet",
                        f$notes)))
  # no phantom "NA" rater level
  expect_false("NA" %in% f$facet_effects$level)
})

test_that("equating excludes unusable items instead of returning all NA", {
  f <- rasch(as.data.frame(mkX(400, 8)))
  ref <- data.frame(item = paste0("I", 1:8),
                    location = f$items$location + 0.3,
                    se = f$items$se)
  ref$se[2] <- NA                       # e.g. a weakly determined item
  eq <- equate_tests(f, ref)
  expect_true(is.finite(eq$shift) && is.finite(eq$rmsd))
  expect_equal(eq$n, 7L)
  expect_true(is.na(eq$table$t[eq$table$item == "I2"]))
  expect_equal(sum(is.finite(eq$table$t)), 7L)
  expect_match(eq$note, "I2")
  expect_no_error(plot_equate(f, ref))
  ref$se[1:7] <- NA                     # fewer than two usable -> error
  expect_error(equate_tests(f, ref), "fewer than two common items")
})

test_that("report_html escapes data-derived text", {
  X <- mkX(200, 6)
  colnames(X) <- c(paste0("I", 1:5), "A<b>&x")
  f <- rasch(as.data.frame(X))
  out <- tempfile(fileext = ".html")
  report_html(f, out, title = "T <script>alert(1)</script>")
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_false(grepl("<script>alert", html, fixed = TRUE))
  expect_true(grepl("&lt;script&gt;alert", html, fixed = TRUE))
  expect_false(grepl("A<b>&x", html, fixed = TRUE))
  expect_true(grepl("A&lt;b&gt;&amp;x", html, fixed = TRUE) ||
              grepl("A&lt;b>&amp;x", html, fixed = TRUE))
  unlink(out)
})

test_that("every public estimator rejects fractional scores", {
  X <- matrix(c(0, 1, 1.9, 1, 0, 1, 0, 1, 1, 0, 1, 0), 6, 2,
              dimnames = list(NULL, c("A", "B")))
  expect_error(pcml(X), "non-integer")
  expect_error(pcml_pc(X), "non-integer")
  long <- data.frame(person = rep(sprintf("P%02d", 1:20), each = 2),
                     item = rep(c("A", "B"), 20),
                     score = c(1.9, rep(c(0, 1, 1, 0), 9), 0, 1, 1))
  expect_error(rasch_mfrm(long, "person", "item", "score", facets = NULL),
               "non-integer")
  d <- data.frame(a = rep("X", 30), b = rep("Y", 30),
                  resp = rep(c(0, 1, 1.5), 10))
  expect_error(btl(d, "a", "b", response = "resp"), "non-integer")
})

test_that("item_moments is overflow-stable and person_wle survives wide items", {
  im <- item_moments(8, seq(-3, 3, length.out = 30))
  expect_true(all(is.finite(unlist(im))))
  expect_equal(sum(im$P), 1, tolerance = 1e-12)
  w <- person_wle(list(seq(-3, 3, length.out = 30)))
  expect_true(all(is.finite(w$theta)))
})

test_that("secondary simulated trait keeps the requested mean and sd", {
  d <- simulate_rasch(n_persons = 20000, n_items = 6, theta_mean = 2,
                      theta_sd = 1, second_dim = list(items = 4:6, rho = 0.5),
                      seed = 1)
  t2 <- attr(d, "truth")$theta2
  expect_false(is.null(t2))
  expect_lt(abs(mean(t2) - 2), 0.05)
  expect_lt(abs(sd(t2) - 1), 0.05)
})

test_that("a judge in two panels is rejected by btl_efrm", {
  d <- simulate_btl_efrm(n_objects_per_set = 5, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 6, reps_within = 15,
                         reps_cross = 15, seed = 5)
  jj <- unique(d$judge)[1]
  d$panel[d$judge == jj][1] <- setdiff(unique(d$panel), d$panel[d$judge == jj][1])[1]
  expect_error(
    btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
             attr(d, "truth")$object_sets, se_method = "conditional"),
    "more than one panel")
})

test_that("OSI is withheld when the clustered covariance is rank-deficient", {
  set.seed(6)
  K <- 10; beta <- seq(-1.5, 1.5, length.out = K); n <- 600
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib]))
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(win == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:5), n, TRUE))
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  expect_true(is.na(f$osi$PSI))
  expect_true(any(grepl("OSI is withheld", f$notes)))
})

test_that("alpha is NA, not -Inf, when the total score is constant", {
  X <- cbind(I1 = rep(c(0L, 1L), 40), I2 = rep(c(1L, 0L), 40))
  f <- suppressWarnings(rasch(X, n_groups = 2))
  expect_true(is.na(f$alpha$alpha))
})

test_that("factor scores cannot bypass the integer guard", {
  long <- data.frame(person = rep(sprintf("P%02d", 1:20), each = 2),
                     item = rep(c("A", "B"), 20),
                     score = factor(c("1.9", rep(c("0", "1", "1", "0"), 9),
                                      "0", "1", "1")))
  expect_error(rasch_mfrm(long, "person", "item", "score", facets = NULL),
               "non-integer")
  expect_error(pcml(matrix(c(0, 1, Inf, 1, 0, 1, 0, 1), 4, 2)), "non-finite")
  expect_error(pcml(matrix(c("0", "1", "abc", "1", "0", "1", "0", "1"), 4, 2)),
               "non-numeric")
})

test_that("graded btl requires an ordered factor", {
  d <- data.frame(a = rep(c("X", "Y", "Z"), 40), b = rep(c("Y", "Z", "X"), 40))
  set.seed(1)
  resp <- sample(c("worse", "same", "better"), 120, TRUE)
  d$plain <- factor(resp)
  d$ord <- factor(resp, levels = c("worse", "same", "better"), ordered = TRUE)
  expect_error(btl(d, "a", "b", response = "plain"), "ORDERED")
  f <- btl(d, "a", "b", response = "ord")
  expect_identical(f$categories, c("worse", "same", "better"))
})

test_that("clustered dependence tests use a t reference with G - 1 df", {
  set.seed(3)
  K <- 6; b <- seq(-1, 1, length.out = K); n <- 800
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(
                    rbinom(n, 1, plogis(b[ia] - b[ib] + 0.4)) == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:6), n, TRUE))
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge",
           position = TRUE)
  dp <- f$dependence
  expect_equal(unique(dp$df), 5L)
  expect_true("t" %in% names(dp))                   # labelled for its reference
  expect_equal(dp$p, 2 * pt(-abs(dp$t), df = 5), tolerance = 1e-12)
  expect_true(all(dp$p >= 2 * pnorm(-abs(dp$t))))   # wider than normal theory
})

test_that("simulate_rasch validates the second-dimension specification", {
  expect_error(simulate_rasch(50, 6, second_dim = list(items = 4:6, rho = 1.2)),
               "correlation in")
  expect_error(simulate_rasch(50, 6, second_dim = list(items = "I99", rho = .5)),
               "unknown item")
})

test_that("margin ordering must be explicit (ordered factor or numeric)", {
  d <- data.frame(a = rep(c("X", "Y", "Z"), 40), b = rep(c("Y", "Z", "X"), 40))
  set.seed(1)
  d$win <- ifelse(runif(120) < .5, d$a, d$b)
  mg <- sample(c("small", "large"), 120, TRUE)
  d$m_plain <- factor(mg); d$m_chr <- mg
  d$m_ord <- factor(mg, levels = c("small", "large"), ordered = TRUE)
  expect_error(btl(d, "a", "b", winner = "win", margin = "m_plain"), "ORDERED")
  expect_error(btl(d, "a", "b", winner = "win", margin = "m_chr"), "character")
  expect_s3_class(btl(d, "a", "b", winner = "win", margin = "m_ord"),
                  "rasch_btl")
})

test_that("invalid graded numeric responses error instead of being dropped", {
  d <- data.frame(a = rep(c("X", "Y", "Z"), 40), b = rep(c("Y", "Z", "X"), 40))
  d$r_chr <- rep(c("0", "1", "abc"), 40)
  d$r_inf <- rep(c(0, 1, Inf), 40)
  expect_error(btl(d, "a", "b", response = "r_chr"), "non-numeric")
  expect_error(btl(d, "a", "b", response = "r_inf"), "non-finite")
})

test_that("simulator rejects malformed second-dimension specifications", {
  expect_error(simulate_rasch(50, 6, second_dim = list(items = 4.9, rho = .5)),
               "whole numbers")
  expect_error(simulate_rasch(50, 6, second_dim = list(items = integer(0),
                                                       rho = .5)),
               "at least one item")
  expect_error(simulate_rasch(50, 6, second_dim = list(items = 4:6,
                                                       rho = c(.5, .6))),
               "single correlation")
})
