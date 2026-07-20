# rasch :: data simulation
#
# Generate data from the model family with dial-in departures from it, so a
# known pathology can be planted and the matching diagnostic watched as it
# fires. Every simulator returns data ready for its fit function, with the
# true generating parameters attached as attr(x, "truth") and a class that
# prints a summary of what was planted.
#
# Shared truth schema (all simulators populate what applies):
#   list(layout, n_*, theta/locations, difficulty/thresholds, discrimination,
#        guessing, groups, planted = <character, human-readable pathologies>)

# null-coalescing helper (package-internal; base R gained %||% only in 4.4)
`%||%` <- function(a, b) if (is.null(a)) b else a

# person locations from one of a few distributions
.sim_theta <- function(n, mean, sd, dist = "normal") {
  z <- switch(dist,
    normal  = stats::rnorm(n),
    uniform = stats::runif(n, -sqrt(3), sqrt(3)),
    skew    = { u <- stats::rgamma(n, 2, 1); (u - 2) / sqrt(2) },
    bimodal = { s <- sample(c(-1, 1), n, TRUE); s * 1.1 + stats::rnorm(n, 0, 0.5) },
    stats::rnorm(n))
  mean + sd * as.numeric(scale(z))
}

# draw one item's responses for a vector of person locations. tau are the
# item's thresholds (length m; dichotomous m = 1). disc scales the whole
# exponent (a departure when != 1); guess is a lower asymptote (dichotomous).
.sim_item <- function(theta, tau, disc = 1, guess = 0) {
  m <- length(tau); xs <- 0:m
  cum <- c(0, cumsum(tau))
  eta <- disc * (outer(theta, xs) - matrix(cum, length(theta), m + 1L, byrow = TRUE))
  eta <- eta - apply(eta, 1, max)
  P <- exp(eta); P <- P / rowSums(P)
  if (guess > 0 && m == 1L) {
    P[, 2] <- guess + (1 - guess) * P[, 2]; P[, 1] <- 1 - P[, 2]
  }
  cs <- t(apply(P, 1, cumsum))
  as.integer(rowSums(stats::runif(length(theta)) > cs))     # category 0..m
}

# thresholds around an item location: evenly spread (rating-scale pattern),
# or with a supplied relative pattern (partial-credit: pattern varies by
# item); optionally deliberately disordered (one threshold dropped below its
# predecessor, re-centred so the item's mean location is preserved)
.sim_thresholds <- function(delta, m, spread, disordered = FALSE,
                            pattern = NULL) {
  if (m == 1L) return(delta)
  step <- if (is.null(pattern)) seq(-1, 1, length.out = m) * spread
          else pattern
  tau <- delta + step - mean(step)
  if (disordered && m >= 2L) {
    i <- max(2L, ceiling(m / 2))           # always has a predecessor to undercut
    tau[i] <- tau[i] - 2.2 * spread
    tau <- tau - mean(tau) + delta         # keep the item's location honest
  }
  tau
}

#' Simulate person-by-item Rasch data with dial-in misfit
#'
#' Generates dichotomous or polytomous (partial credit / rating scale) data
#' from the Rasch model, with optional, individually controllable departures
#' from it -- each of which the package's matching diagnostic is built to
#' detect. The result is a data frame ready for \code{\link{rasch}}, with the
#' true parameters attached as \code{attr(x, "truth")}.
#'
#' @param n_persons,n_items Sample size and test length.
#' @param model \code{"dichotomous"}, \code{"PCM"}, or \code{"RSM"}. Under
#'   \code{"RSM"} every item shares one category-threshold pattern (items
#'   differ by location only); under \code{"PCM"} each item's threshold
#'   spacings and span are drawn afresh, as the partial credit model allows.
#' @param n_categories Response categories for polytomous models (>= 3).
#' @param theta_mean,theta_sd,theta_dist Person distribution: mean, SD, and
#'   shape (\code{"normal"}, \code{"uniform"}, \code{"skew"}, \code{"bimodal"}).
#' @param difficulty Two numbers giving the item-location range (evenly
#'   spaced), or a length-\code{n_items} vector of locations.
#' @param threshold_spread Half-range of the category thresholds about each
#'   item location (polytomous).
#' @param discrimination Scalar or length-\code{n_items}: the slope of each
#'   item. Values above 1 over-discriminate (Guttman-like, negative fit
#'   residual); below 1 under-discriminate (noisy, positive residual). Feeds
#'   infit/outfit and the item-fit F.
#' @param guessing Scalar or length-\code{n_items} lower asymptote
#'   (dichotomous): low-ability persons answer correctly by chance. Feeds
#'   \code{\link{tailored_analysis}}.
#' @param second_dim \code{NULL}, or \code{list(items=, rho=)}: the named items
#'   load on a second trait correlated \code{rho} with the first. Feeds
#'   \code{\link{dimensionality_test}}.
#' @param dependence \code{NULL}, or \code{list(pairs=, strength=)}: each pair's
#'   second item responds partly to the first (response dependence). Feeds
#'   \code{\link{residual_correlations}} / \code{\link{dependence_magnitude}}.
#' @param dif \code{NULL}, or \code{list(items=, uniform=, nonuniform=)}: the
#'   named items function differently for the last person group -- a location
#'   shift (\code{uniform}) and/or a slope change (\code{nonuniform}). Needs
#'   \code{n_groups >= 2}. Feeds \code{\link{dif_anova}} / \code{\link{dif_size}}.
#' @param careless Proportion of persons who answer at random (person misfit;
#'   feeds person infit/outfit).
#' @param response_style \code{NULL}, or \code{list(type=, prop=, strength=)}
#'   with \code{type} \code{"extreme"} or \code{"middle"}: a proportion
#'   \code{prop} of persons favour the end (or middle) categories regardless
#'   of the trait, with distortion \code{strength} (default 1.6) on the
#'   log-probability scale (polytomous; feeds the category diagnostics and
#'   person fit).
#' @param speeded Proportion not-reached at the last item: a growing tail of
#'   missing responses over the final items, as under time pressure (feeds the
#'   item statistics and the missingness pattern).
#' @param disordered \code{NULL} or item names/indices given disordered
#'   thresholds (polytomous; feeds the threshold diagnostics).
#' @param n_groups Number of equal person groups (a \code{group} factor column
#'   is added when > 1, for DIF).
#' @param missing Proportion of responses set missing (completely at random).
#' @param seed Optional RNG seed.
#' @return A data frame of class \code{"rasch_sim"} (item columns
#'   \code{I01}..., an \code{id} column, and a \code{group} column when
#'   grouped), with \code{attr(x, "truth")} holding the generating parameters
#'   and the planted departures.
#' @examples
#' # a clean scale with one over-discriminating item and one DIF item
#' d <- simulate_rasch(400, 12, discrimination = c(3, rep(1, 11)),
#'                     dif = list(items = "I06", uniform = 1), n_groups = 2,
#'                     seed = 1)
#' fit <- rasch(d, id = "id", factors = "group")
#' fit$items[c("item", "infit_ms", "outfit_ms")]   # item 1 misfits
#' dif_anova(fit)$summary                           # item 6 flags
#' @export
simulate_rasch <- function(n_persons = 500, n_items = 20,
                           model = c("dichotomous", "PCM", "RSM"),
                           n_categories = 3, theta_mean = 0, theta_sd = 1,
                           theta_dist = "normal", difficulty = c(-2.5, 2.5),
                           threshold_spread = 1.2, discrimination = 1,
                           guessing = 0, second_dim = NULL, dependence = NULL,
                           dif = NULL, careless = 0, response_style = NULL,
                           speeded = 0, disordered = NULL,
                           n_groups = 1, missing = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  model <- match.arg(model)
  m <- if (model == "dichotomous") 1L else as.integer(n_categories) - 1L
  I <- as.integer(n_items); N <- as.integer(n_persons)
  inm <- sprintf("I%02d", seq_len(I))
  as_idx <- function(x) {
    if (is.null(x) || !length(x)) return(integer(0))
    if (is.numeric(x) && any(!is.na(x) & x != round(x)))
      stop("item index(es) must be whole numbers, got: ",
           paste(x[!is.na(x) & x != round(x)], collapse = ", "))
    i <- if (is.character(x)) match(x, inm) else as.integer(x)
    if (anyNA(i) || any(i < 1L | i > I))
      stop("unknown item name(s)/index(es): ",
           paste(x[is.na(i) | i < 1L | i > I], collapse = ", "))
    i
  }

  # item locations, thresholds, slopes, guessing (with per-item overrides)
  delta <- setNames(if (length(difficulty) == I) difficulty
                    else seq(difficulty[1], difficulty[2], length.out = I), inm)
  disc <- if (length(discrimination) == I) discrimination else rep(discrimination[1], I)
  guess <- if (length(guessing) == I) guessing else rep(guessing[1], I)
  if (m > 1L && any(guess > 0)) {
    warning("guessing applies to dichotomous items only; ignored for ", model)
    guess[] <- 0
  }
  dis_items <- as_idx(disordered)
  # RSM: one common threshold pattern (per-item location only). PCM: the
  # pattern itself varies across items, as the model allows -- each item's
  # spacings are drawn afresh (gated so the RNG stream of the other models
  # is untouched)
  patterns <- vector("list", I)
  if (model == "PCM" && m > 1L) patterns <- lapply(seq_len(I), function(i) {
    # ordered spacings with varied gaps and a varied overall span per item
    p <- cumsum(stats::runif(m, 0.5, 1.5))
    p <- (p - mean(p)) / (max(p) - min(p)) * 2 * threshold_spread *
      stats::runif(1, 0.85, 1.15)
    p
  })
  tau <- lapply(seq_len(I), function(i)
    .sim_thresholds(delta[i], m, threshold_spread, i %in% dis_items,
                    pattern = patterns[[i]]))

  # person locations (primary) and groups
  theta <- .sim_theta(N, theta_mean, theta_sd, theta_dist)
  group <- if (n_groups > 1L)
    factor(sprintf("g%d", (seq_len(N) - 1L) %% n_groups + 1L)) else NULL

  # a second dimension for the nominated items: a correlated latent trait
  theta2 <- NULL; dim_items <- integer(0)
  if (!is.null(second_dim)) {
    dim_items <- as_idx(second_dim$items)
    if (!length(dim_items))
      stop("second_dim$items must name at least one item")
    rho <- second_dim$rho %||% 0.5
    if (length(rho) != 1L || !is.finite(rho) || abs(rho) > 1)
      stop("second_dim$rho must be a single correlation in [-1, 1]")
    # centre both components so the secondary trait keeps the REQUESTED
    # mean and sd: combining two mean-mu variables shifted the mean to
    # mu(rho + sqrt(1 - rho^2))
    z2 <- .sim_theta(N, theta_mean, theta_sd, theta_dist)
    theta2 <- theta_mean + rho * (theta - theta_mean) +
      sqrt(1 - rho^2) * (z2 - theta_mean)
  }

  X <- matrix(NA_integer_, N, I, dimnames = list(NULL, inm))
  dif_items <- as_idx(if (is.null(dif)) NULL else dif$items)
  if (length(dif_items) && n_groups < 2L)
    stop("dif needs n_groups >= 2 (the last group carries the DIF)")
  dif_grp <- if (n_groups > 1L) levels(group)[n_groups] else NA

  # every regeneration of an item must honour that item's OWN generating
  # structure -- its trait (second dimension), its group-shifted thresholds
  # and slope (DIF), its guessing -- or a later misfit layer would silently
  # erase an earlier planted one
  item_pars <- function(i, p) {
    dif_here <- i %in% dif_items && !is.na(dif_grp) && group[p] == dif_grp
    list(tau = tau[[i]] + if (dif_here) dif$uniform %||% 0 else 0,
         disc = disc[i] + if (dif_here) dif$nonuniform %||% 0 else 0)
  }
  gen_item <- function(i, shift = rep(0, N)) {
    th <- (if (i %in% dim_items) theta2 else theta) + shift
    if (i %in% dif_items && !is.na(dif_grp)) {
      g2 <- group == dif_grp
      out <- integer(N)
      out[!g2] <- .sim_item(th[!g2], tau[[i]], disc[i], guess[i])
      out[g2]  <- .sim_item(th[g2], tau[[i]] + (dif$uniform %||% 0),
                            disc[i] + (dif$nonuniform %||% 0), guess[i])
      out
    } else .sim_item(th, tau[[i]], disc[i], guess[i])
  }
  # the model expectation of item i for every person, under the same
  # generating structure (trait, DIF shift, guessing) the responses used
  exp_item <- function(i) {
    th <- if (i %in% dim_items) theta2 else theta
    E <- vapply(seq_len(N), function(p) {
      pp <- item_pars(i, p)
      sum((0:m) * .p_item(th[p], pp$tau, pp$disc))
    }, 0)
    if (m == 1L && guess[i] > 0) E <- guess[i] + (1 - guess[i]) * E
    E
  }

  for (i in seq_len(I)) X[, i] <- gen_item(i)

  # response dependence: the second item of each pair partly follows the
  # first (adds d*(x1 - E1) to its exponent, inducing residual correlation);
  # the regeneration keeps i2's own DIF / second-dimension structure
  dep_pairs <- list()
  if (!is.null(dependence)) {
    d_str <- dependence$strength %||% 1
    for (pp in dependence$pairs) {
      ij <- as_idx(pp); i1 <- ij[1]; i2 <- ij[2]
      # the expectation must match X1's actual generating structure
      # (guessing, DIF, second dimension), or the "residual" x1 - E1 has a
      # systematic mean and leaks an unplanted shift into the second item
      shift <- d_str * (X[, i1] - exp_item(i1)) / m   # per-person carry-over
      X[, i2] <- gen_item(i2, shift = shift)
      dep_pairs[[length(dep_pairs) + 1L]] <- inm[ij]
    }
  }

  # response styles (polytomous): a proportion of persons distort toward the
  # end or the middle categories regardless of the trait; the base
  # probabilities keep each item's own structure (trait, DIF) per person
  style_idx <- integer(0)
  if (!is.null(response_style) && m >= 2L) {
    style_idx <- sample(N, round((response_style$prop %||% 0.15) * N))
    ss <- response_style$strength %||% 1.6; mid <- m / 2
    dev2 <- ((0:m - mid) / mid)^2
    w <- if ((response_style$type %||% "extreme") == "extreme") exp(ss * dev2)
         else exp(-ss * dev2)
    for (p in style_idx) for (i in seq_len(I)) {
      if (is.na(X[p, i])) next
      th_p <- if (i %in% dim_items) theta2[p] else theta[p]
      pp <- item_pars(i, p)
      pr <- .p_item(th_p, pp$tau, pp$disc) * w
      X[p, i] <- sample.int(m + 1L, 1L, prob = pr) - 1L
    }
  }

  # careless responders: answer uniformly at random
  careless_idx <- integer(0)
  if (careless > 0) {
    careless_idx <- sample(N, round(careless * N))
    X[careless_idx, ] <- matrix(sample(0:m, length(careless_idx) * I, TRUE),
                                length(careless_idx), I)
    # careless overwrites any response style; the truth must not double-count
    style_idx <- setdiff(style_idx, careless_idx)
  }

  # speededness: a contiguous not-reached tail over the last items. `speeded`
  # persons drop out somewhere in the final zone, so the missing rate grows
  # linearly to `speeded` at the last item
  if (speeded > 0 && I >= 3L) {
    k <- max(1L, round(0.4 * I)); z0 <- I - k
    for (p in which(stats::runif(N) < speeded)) {
      stop_at <- z0 + sample.int(k, 1L) - 1L              # drop point in zone
      if (stop_at < I) X[p, (stop_at + 1L):I] <- NA
    }
  }
  if (missing > 0) X[sample(length(X), round(missing * length(X)))] <- NA

  out <- data.frame(id = sprintf("P%04d", seq_len(N)), X,
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(group)) out$group <- group

  planted <- character(0)
  if (any(disc != 1)) planted <- c(planted, sprintf("discrimination != 1 on %s",
    paste(inm[disc != 1], collapse = ", ")))
  if (any(guess > 0)) planted <- c(planted, sprintf("guessing on %s",
    paste(inm[guess > 0], collapse = ", ")))
  if (length(dim_items)) planted <- c(planted, sprintf(
    "second dimension (rho %.2f) on %s", second_dim$rho %||% 0.5,
    paste(inm[dim_items], collapse = ", ")))
  if (length(dep_pairs)) planted <- c(planted, sprintf(
    "response dependence: %s", paste(vapply(dep_pairs, paste, "",
                                            collapse = "-"), collapse = "; ")))
  if (length(dif_items)) planted <- c(planted, sprintf(
    "DIF (group %s) on %s: uniform %.2f, non-uniform %.2f", dif_grp,
    paste(inm[dif_items], collapse = ", "), dif$uniform %||% 0,
    dif$nonuniform %||% 0))
  if (length(careless_idx)) planted <- c(planted, sprintf(
    "%d careless responder(s)", length(careless_idx)))
  if (length(style_idx)) planted <- c(planted, sprintf(
    "%s response style on %d person(s)",
    response_style$type %||% "extreme", length(style_idx)))
  if (speeded > 0) planted <- c(planted, sprintf(
    "speededness (%.0f%% not-reached at the last item)", 100 * speeded))
  if (length(dis_items)) planted <- c(planted, sprintf(
    "disordered thresholds on %s", paste(inm[dis_items], collapse = ", ")))
  if (missing > 0) planted <- c(planted, sprintf("%.0f%% missing", 100 * missing))

  attr(out, "truth") <- list(
    layout = "rasch",
    description = sprintf("%s, %d persons x %d items%s", model, N, I,
      if (!is.null(group)) sprintf(", %d groups", nlevels(group)) else ""),
    model = model, n_persons = N, n_items = I,
    theta = theta, theta2 = theta2, difficulty = delta, thresholds = tau,
    discrimination = disc, guessing = guess,
    groups = group, dim_items = inm[dim_items], dif_items = inm[dif_items],
    careless_idx = careless_idx, style_idx = style_idx, planted = planted)
  class(out) <- c("rasch_sim", "data.frame")
  out
}

# category probabilities for one location (used by the dependence term)
.p_item <- function(theta, tau, disc = 1) {
  m <- length(tau); cum <- c(0, cumsum(tau))
  e <- disc * ((0:m) * theta - cum); e <- e - max(e)
  p <- exp(e); p / sum(p)
}

#' @export
print.rasch_sim <- function(x, ...) {
  tr <- attr(x, "truth")
  cat(sprintf("Simulated %s data: %s\n", tr$layout,
              tr$description %||% sprintf("%d rows", nrow(x))))
  if (length(tr$planted)) {
    cat("Planted departures:\n")
    for (p in tr$planted) cat(paste0("  - ", p, "\n"))
  } else cat("Model-conforming (no departures planted).\n")
  invisible(x)
}

#' Simulate paired-comparison (BTL) data with dial-in misfit
#'
#' Generates dichotomous or graded paired comparisons from the
#' Bradley-Terry-Luce model, with optional departures each of which a
#' paired-comparison diagnostic is built to detect. The result is a data frame
#' ready for \code{\link{btl}}, with the truth attached.
#'
#' @param n_objects,n_judges Objects to scale and judges comparing them.
#' @param reps_per_pair Comparisons made of each object pair.
#' @param model \code{"dichotomous"} (a winner) or \code{"graded"} (a rated
#'   margin in \code{n_categories} categories).
#' @param n_categories Categories for the graded model.
#' @param object_sd Spread of the object locations (evenly spaced, sum-zero).
#' @param second_attribute \code{NULL}, or \code{list(rho=)}: half the judges
#'   rank by a second object attribute correlated \code{rho} with the first --
#'   genuine multidimensionality. Feeds \code{\link{btl_dimensionality}} and
#'   \code{\link{btl_transitivity}}.
#' @param erratic_judges Proportion of judges who choose at random. Feeds the
#'   judge fit residual, \code{\link{btl_transitivity}} consistency, and
#'   \code{\link{judge_surprise}}.
#' @param dependence \code{NULL}, or \code{list(exposure=, carry_over=)}:
#'   within-judge order effects (a seen-before advantage and a pull from the
#'   judge's own earlier verdicts). Adds an \code{order} column. Feeds the
#'   dependence effects of \code{\link{btl}}.
#' @param seed Optional RNG seed.
#' @return A data frame of class \code{"rasch_sim"}: \code{object_a},
#'   \code{object_b}, \code{winner} (or \code{response} when graded),
#'   \code{judge}, and \code{order} when dependence is planted; with
#'   \code{attr(x, "truth")}.
#' @examples
#' d <- simulate_btl(8, 12, erratic_judges = 0.15, seed = 1)
#' bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
#' bt$judges          # the erratic judges carry large fit residuals
#' @export
simulate_btl <- function(n_objects = 8, n_judges = 12, reps_per_pair = 25,
                         model = c("dichotomous", "graded"), n_categories = 4,
                         object_sd = 1, second_attribute = NULL,
                         erratic_judges = 0, dependence = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  model <- match.arg(model)
  m <- if (model == "graded") as.integer(n_categories) - 1L else 1L
  K <- as.integer(n_objects); J <- as.integer(n_judges)
  objs <- sprintf("O%d", seq_len(K)); jids <- sprintf("J%d", seq_len(J))
  beta <- setNames(as.numeric(scale(seq_len(K))) * object_sd, objs)
  tau <- if (m > 1L) .sim_thresholds(0, m, 1.2) else NULL

  # a second object attribute (orthogonal part) for the two-camp design
  beta2 <- NULL; camp <- NULL
  if (!is.null(second_attribute)) {
    rho <- second_attribute$rho %||% 0.3
    beta2 <- setNames(rho * beta + sqrt(1 - rho^2) *
      as.numeric(scale(stats::rnorm(K))) * object_sd, objs)
    camp <- setNames(rep(c("a", "b"), length.out = J), jids)
  }
  erratic <- if (erratic_judges > 0)
    jids[seq_len(round(erratic_judges * J))] else character(0)

  pr <- t(utils::combn(objs, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps_per_pair),
                  object_b = rep(pr[, 2], each = reps_per_pair),
                  stringsAsFactors = FALSE)
  d$judge <- sample(jids, nrow(d), TRUE)

  win_prob <- function(a, b, jd) {
    ba <- if (!is.null(camp) && camp[jd] == "b") beta2[a] else beta[a]
    bb <- if (!is.null(camp) && camp[jd] == "b") beta2[b] else beta[b]
    ba - bb
  }
  # dependence needs a per-judge judgment order and running history
  ord <- NULL
  if (!is.null(dependence)) {
    d <- d[order(d$judge), ]
    d$order <- stats::ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
    seen <- new.env(parent = emptyenv()); hs <- new.env(parent = emptyenv())
    hc <- new.env(parent = emptyenv())
    g0 <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
    exq <- dependence$exposure %||% 0; cry <- dependence$carry_over %||% 0
    resp <- integer(nrow(d))
    for (r in seq_len(nrow(d))) {
      j <- d$judge[r]; a <- d$object_a[r]; b <- d$object_b[r]
      ka <- paste(j, a); kb <- paste(j, b)
      lp <- win_prob(a, b, j) +
        exq * (as.numeric(g0(seen, ka) > 0) - as.numeric(g0(seen, kb) > 0)) +
        cry * ((if (g0(hc, ka) > 0) g0(hs, ka) / g0(hc, ka) else 0) -
               (if (g0(hc, kb) > 0) g0(hs, kb) / g0(hc, kb) else 0))
      x <- if (j %in% erratic) sample(0:m, 1)
           else if (m == 1L) as.integer(stats::runif(1) < stats::plogis(lp))
           else sample(0:m, 1, prob = .p_item(lp, tau))
      resp[r] <- x
      assign(ka, g0(seen, ka) + 1, seen); assign(kb, g0(seen, kb) + 1, seen)
      assign(ka, g0(hc, ka) + 1, hc);     assign(kb, g0(hc, kb) + 1, hc)
      assign(ka, g0(hs, ka) + (2 * x / m - 1), hs)
      assign(kb, g0(hs, kb) + (2 * (m - x) / m - 1), hs)
    }
  } else {
    lp <- vapply(seq_len(nrow(d)), function(r)
      win_prob(d$object_a[r], d$object_b[r], d$judge[r]), 0)
    resp <- integer(nrow(d))
    reg <- !(d$judge %in% erratic)
    resp[reg] <- if (m == 1L) as.integer(stats::runif(sum(reg)) < stats::plogis(lp[reg]))
                 else vapply(which(reg), function(r) sample(0:m, 1, prob = .p_item(lp[r], tau)), 0L)
    if (any(!reg)) resp[!reg] <- sample(0:m, sum(!reg), TRUE)
  }

  if (m == 1L) d$winner <- ifelse(resp == 1L, d$object_a, d$object_b)
  else d$response <- resp
  rownames(d) <- NULL

  planted <- character(0)
  if (length(erratic)) planted <- c(planted,
    sprintf("%d erratic judge(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (!is.null(second_attribute)) planted <- c(planted,
    sprintf("second object attribute (rho %.2f), two judge camps",
            second_attribute$rho %||% 0.3))
  if (!is.null(dependence)) planted <- c(planted, sprintf(
    "within-judge dependence: exposure %.2f, carry-over %.2f",
    dependence$exposure %||% 0, dependence$carry_over %||% 0))

  attr(d, "truth") <- list(
    layout = "btl",
    description = sprintf("%s, %d objects, %d judges, %d comparisons",
                          model, K, J, nrow(d)),
    model = model, location = beta, location2 = beta2, camp = camp,
    erratic = erratic, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}

#' Simulate many-facet (rated) data with dial-in misfit
#'
#' Generates ratings from the many-facet Rasch model (Linacre 1989): every
#' rater rates every person on every item, from person ability, item
#' difficulty, and rater severity. Departures each feed an MFRM diagnostic.
#'
#' @param n_persons,n_items,n_raters Facet sizes (fully crossed).
#' @param n_categories Rating categories.
#' @param theta_sd,item_sd Spread of person ability and item difficulty.
#' @param rater_severity_sd Spread of rater severities (the core facet;
#'   recovered in \code{facet_effects}).
#' @param erratic_raters Proportion of raters who rate at random (feeds the
#'   rater fit residual).
#' @param interaction \code{NULL}, or \code{list(rater=, item=, bias=)}: one
#'   rater is unusually harsh (positive) or lenient (negative) on one item.
#'   Feeds the item-by-rater interaction (fit with \code{interaction = }).
#' @param halo Proportion of raters showing a halo effect: they rate by the
#'   person's overall level and barely differentiate items (feeds the rater
#'   fit residual and the item-by-rater interaction).
#' @param seed Optional RNG seed.
#' @return A long data frame of class \code{"rasch_sim"} (\code{person},
#'   \code{item}, \code{rater}, \code{score}) ready for
#'   \code{\link{rasch_mfrm}}, with the truth attached.
#' @examples
#' d <- simulate_mfrm(60, 5, 6, rater_severity_sd = 0.8, seed = 1)
#' mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
#'                  facets = "rater")
#' cor(mf$facet_effects$rater$severity, attr(d, "truth")$severity)  # recovered
#' @export
simulate_mfrm <- function(n_persons = 80, n_items = 5, n_raters = 6,
                          n_categories = 4, theta_sd = 1.2, item_sd = 1,
                          rater_severity_sd = 0.6, erratic_raters = 0,
                          interaction = NULL, halo = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- as.integer(n_categories) - 1L
  N <- as.integer(n_persons); I <- as.integer(n_items); R <- as.integer(n_raters)
  pids <- sprintf("P%03d", seq_len(N)); iids <- sprintf("I%d", seq_len(I))
  rids <- sprintf("R%d", seq_len(R))
  theta <- .sim_theta(N, 0, theta_sd)
  delta <- setNames(seq(-item_sd, item_sd, length.out = I), iids)
  lambda <- setNames(as.numeric(scale(stats::rnorm(R))) * rater_severity_sd, rids)
  base_tau <- .sim_thresholds(0, m, 1.2)
  erratic <- if (erratic_raters > 0) rids[seq_len(round(erratic_raters * R))] else character(0)
  # halo raters (drawn from the end, disjoint from the erratic ones): they
  # rate by the person's overall level, barely differentiating the items;
  # capped at the eligible pool so the truth never records NA raters
  halo_r <- if (halo > 0) {
    pool <- setdiff(rev(rids), erratic)
    pool[seq_len(min(length(pool), round(halo * R)))]
  } else character(0)
  int_bias <- matrix(0, I, R, dimnames = list(iids, rids))
  if (!is.null(interaction))
    int_bias[interaction$item, interaction$rater] <- interaction$bias

  grid <- expand.grid(p = seq_len(N), i = seq_len(I), r = seq_len(R))
  score <- integer(nrow(grid))
  for (i in seq_len(I)) for (r in seq_len(R)) {
    rows <- grid$i == i & grid$r == r
    # item difficulty and rater severity shift the person's thresholds; a halo
    # rater ignores the item's own difficulty (uses the mean instead)
    di <- if (rids[r] %in% halo_r) mean(delta) else delta[i]
    tau_ir <- base_tau + di + lambda[r] + int_bias[i, r]
    score[rows] <- if (rids[r] %in% erratic) sample(0:m, sum(rows), TRUE)
                   else .sim_item(theta[grid$p[rows]], tau_ir)
  }
  d <- data.frame(person = pids[grid$p], item = iids[grid$i],
                  rater = rids[grid$r], score = score, stringsAsFactors = FALSE)

  planted <- character(0)
  if (length(erratic)) planted <- c(planted,
    sprintf("%d erratic rater(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (length(halo_r)) planted <- c(planted,
    sprintf("%d halo rater(s): %s", length(halo_r), paste(halo_r, collapse = ", ")))
  if (!is.null(interaction)) planted <- c(planted, sprintf(
    "rater-by-item bias: %s on %s (%.2f)", interaction$rater,
    interaction$item, interaction$bias))
  planted <- c(planted, sprintf("rater severities SD %.2f", stats::sd(lambda)))

  attr(d, "truth") <- list(
    layout = "mfrm",
    description = sprintf("%d persons x %d items x %d raters (%d ratings)",
                          N, I, R, nrow(d)),
    theta = theta, difficulty = delta, severity = lambda,
    erratic = erratic, halo = halo_r, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}

#' Simulate extended frame-of-reference data with differing units
#'
#' Generates data whose latent unit differs across item-set by person-group
#' frames (Humphry 2005): a person in group g responding to an item in set s
#' does so at the frame unit rho = alpha_set * phi_group scaling the whole
#' exponent. The planted set- and group-unit ratios are recovered by
#' \code{\link{rasch_efrm}}.
#'
#' @param n_per_group Persons in each group.
#' @param items_per_set Items in each set.
#' @param n_sets,n_groups Numbers of item sets and person groups.
#' @param set_unit_ratio,group_unit_ratio Geometric span of the set and group
#'   units across their levels (1 = equal units, i.e. an ordinary Rasch fit).
#' @param theta_sd Spread of person ability.
#' @param seed Optional RNG seed.
#' @return A wide data frame of class \code{"rasch_sim"} (\code{id}, item
#'   columns, \code{group}) with \code{attr(x, "truth")$item_sets} the set map
#'   to pass to \code{\link{rasch_efrm}}.
#' @examples
#' d <- simulate_efrm(300, 8, set_unit_ratio = 1.3, seed = 1)
#' tr <- attr(d, "truth")
#' ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group")
#' ef$alpha_table   # recovers the ~1.3 set-unit ratio
#' @export
simulate_efrm <- function(n_per_group = 300, items_per_set = 8, n_sets = 2,
                          n_groups = 2, set_unit_ratio = 1.3,
                          group_unit_ratio = 1, theta_sd = 1.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  S <- as.integer(n_sets); G <- as.integer(n_groups); K <- as.integer(items_per_set)
  npg <- as.integer(n_per_group)
  # set and group units span the ratio geometrically, normalised to mean 1
  gspan <- function(ratio, n) { u <- exp(seq(0, log(ratio), length.out = n)); u / exp(mean(log(u))) }
  alpha <- gspan(set_unit_ratio, S)
  phi <- gspan(group_unit_ratio, G)
  set_items <- lapply(seq_len(S), function(s) sprintf("S%dI%02d", s, seq_len(K)))
  inm <- unlist(set_items)
  delta <- setNames(rep(seq(-1.5, 1.5, length.out = K), S), inm)
  set_of <- setNames(rep(seq_len(S), each = K), inm)

  grp <- factor(rep(sprintf("g%d", seq_len(G)), each = npg))
  N <- length(grp); theta <- .sim_theta(N, 0, theta_sd)
  X <- matrix(NA_integer_, N, length(inm), dimnames = list(NULL, inm))
  for (col in seq_along(inm)) {
    s <- set_of[inm[col]]; rho <- alpha[s] * phi[as.integer(grp)]  # per-person unit
    X[, col] <- as.integer(stats::runif(N) < stats::plogis(rho * (theta - delta[inm[col]])))
  }
  out <- data.frame(id = sprintf("P%04d", seq_len(N)), X, group = grp,
                    check.names = FALSE, stringsAsFactors = FALSE)

  planted <- sprintf("set-unit ratio %.2f across %d sets", set_unit_ratio, S)
  if (group_unit_ratio != 1)
    planted <- c(planted, sprintf("group-unit ratio %.2f across %d groups",
                                  group_unit_ratio, G))
  attr(out, "truth") <- list(
    layout = "efrm",
    description = sprintf("%d persons, %d sets x %d groups, %d items",
                          N, S, G, length(inm)),
    theta = theta, difficulty = delta, alpha = alpha, phi = phi,
    item_sets = setNames(set_items, sprintf("set%d", seq_len(S))),
    groups = grp, planted = planted)
  class(out) <- c("rasch_sim", "data.frame")
  out
}

#' Replicate a simulation for Monte Carlo studies
#'
#' Calls one of the \code{simulate_*} functions \code{n} times with successive
#' seeds, returning the datasets as a list -- for power, Type-I, or
#' parameter-recovery studies.
#'
#' @param FUN A simulator, e.g. \code{\link{simulate_rasch}}.
#' @param n Number of datasets.
#' @param ... Arguments passed to \code{FUN} (the same each replicate).
#' @param seed Seed of the first replicate (each subsequent one increments it).
#' @return A list of class \code{"rasch_sim_batch"}, one simulated dataset per
#'   element.
#' @examples
#' # 20 datasets with a planted DIF item; how often is it flagged?
#' batch <- sim_replicate(simulate_rasch, 20, n_persons = 400, n_items = 10,
#'                        dif = list(items = "I05", uniform = 0.8), n_groups = 2,
#'                        seed = 1)
#' mean(vapply(batch, function(d)
#'   dif_anova(rasch(d, id = "id", factors = "group"))$summary$uniform_DIF[5], TRUE))
#' @export
sim_replicate <- function(FUN, n, ..., seed = NULL) {
  base <- if (is.null(seed)) sample.int(1e6, 1L) else as.integer(seed)
  reps <- lapply(seq_len(n), function(k) FUN(..., seed = base + k - 1L))
  structure(reps, class = "rasch_sim_batch", n = as.integer(n),
            layout = attr(reps[[1]], "truth")$layout)
}

#' @export
print.rasch_sim_batch <- function(x, ...) {
  cat(sprintf("%d simulated %s datasets (Monte Carlo batch)\n",
              attr(x, "n"), attr(x, "layout")))
  cat("Each element is a simulated dataset; fit and summarise across them, e.g.\n")
  cat("  vapply(batch, function(d) <statistic of fit>, 0)\n")
  invisible(x)
}

#' Parameter recovery of a fit against the simulation truth
#'
#' Compares the parameters recovered by a fit with the ones a
#' \code{simulate_*} function planted (carried on the data as
#' \code{attr(sim, "truth")}): item difficulties and person abilities for a
#' Rasch fit, object locations for a paired-comparison fit, rater severities
#' (with item and person measures) for a many-facet fit, and the set units for
#' a frames fit. Locations are mean-centred before comparison, since the model
#' identifies them only up to an origin.
#'
#' @param fit A fit of the simulated data (\code{\link{rasch}},
#'   \code{\link{btl}}, \code{\link{rasch_mfrm}}, or \code{\link{rasch_efrm}}).
#' @param sim The simulated data (from a \code{simulate_*} function).
#' @return A list of class \code{"rasch_recovery"}: \code{summary} (per
#'   parameter type: n, correlation, RMSE, bias) and \code{pieces} (the true
#'   and estimated values behind each).
#' @examples
#' d <- simulate_rasch(500, 12, seed = 1)
#' sim_recovery(rasch(d, id = "id"), d)$summary
#' @export
sim_recovery <- function(fit, sim) {
  tr <- attr(sim, "truth")
  if (is.null(tr)) stop("`sim` carries no simulation truth")
  pieces <- list()
  add <- function(name, true, est, label = NULL, centre = FALSE) {
    true <- as.numeric(true); est <- as.numeric(est)
    keep <- is.finite(true) & is.finite(est)
    if (!any(keep)) return(invisible())
    if (centre) {                       # location-type: identified up to origin
      true <- true - mean(true[keep]); est <- est - mean(est[keep])
    }
    pieces[[name]] <<- data.frame(
      parameter = name,
      label = if (is.null(label)) NA_character_ else as.character(label)[keep],
      true = true[keep], estimated = est[keep], stringsAsFactors = FALSE)
  }
  lay <- tr$layout
  if (lay == "rasch") {
    ei <- setNames(fit$items$location, fit$items$item)
    cm <- intersect(names(tr$difficulty), names(ei))
    add("item difficulty", tr$difficulty[cm], ei[cm], cm, centre = TRUE)
    add("person ability", tr$theta, fit$person$theta, centre = TRUE)
  } else if (lay == "btl") {
    eo <- setNames(fit$objects$location, fit$objects$object)
    cm <- intersect(names(tr$location), names(eo))
    add("object location", tr$location[cm], eo[cm], cm, centre = TRUE)
  } else if (lay == "mfrm") {
    fe <- fit$facet_effects[[1]]
    es <- setNames(fe$severity, fe$level)
    cm <- intersect(names(tr$severity), names(es))
    add("rater severity", tr$severity[cm], es[cm], cm, centre = TRUE)
    # per-item margins live in item_effects (fit$items holds the virtual
    # item-by-facet combinations, whose names never match)
    ie <- fit$item_effects
    if (!is.null(ie)) {
      ei <- setNames(ie$location, ie$item)
      ci <- intersect(names(tr$difficulty), names(ei))
      add("item difficulty", tr$difficulty[ci], ei[ci], ci, centre = TRUE)
    }
    add("person ability", tr$theta, fit$person$theta, centre = TRUE)
  } else if (lay == "efrm") {
    at <- fit$alpha_table
    # units are identified up to a common scale, so compare on the centred
    # log scale (a ratio); the planted alpha is normalised the same way
    add("set unit (log)", log(tr$alpha), log(at$alpha[match(
      sprintf("set%d", seq_along(tr$alpha)), at$set)]),
      sprintf("set%d", seq_along(tr$alpha)))
  } else stop("unsupported layout: ", lay)

  summ <- do.call(rbind, lapply(pieces, function(d) data.frame(
    parameter = d$parameter[1], n = nrow(d),
    correlation = if (nrow(d) > 2) stats::cor(d$true, d$estimated) else NA_real_,
    rmse = sqrt(mean((d$estimated - d$true)^2)),
    bias = mean(d$estimated - d$true), stringsAsFactors = FALSE)))
  rownames(summ) <- NULL
  structure(list(summary = summ, pieces = pieces, layout = lay),
            class = "rasch_recovery")
}

#' @export
print.rasch_recovery <- function(x, ...) {
  cat("Parameter recovery (planted vs recovered):\n")
  s <- x$summary
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-16s n=%-4d r=%.3f  RMSE=%.3f  bias=%+.3f\n",
                s$parameter[i], s$n[i], s$correlation[i], s$rmse[i], s$bias[i]))
  invisible(x)
}

#' Recovery scatter of planted against recovered parameters
#'
#' One true-versus-estimated panel per parameter type, with the identity line
#' and the correlation and RMSE.
#'
#' @param x A \code{"rasch_recovery"} object.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @export
plot_recovery <- function(x, ...) {
  stopifnot(inherits(x, "rasch_recovery"))
  np <- length(x$pieces)
  op <- par(mfrow = c(1, np), mar = c(4.2, 4.2, 2.4, 1), mgp = c(2.4, 0.7, 0),
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            cex.main = 1.0, font.main = 2)
  on.exit(par(op))
  for (i in seq_len(np)) {
    d <- x$pieces[[i]]; s <- x$summary[i, ]
    rng <- range(c(d$true, d$estimated))
    plot(d$true, d$estimated, xlim = rng, ylim = rng, xlab = "planted",
         ylab = "recovered", main = d$parameter[1], axes = FALSE,
         pch = 21, bg = .rr$blue, col = "white", cex = 1.1)
    abline(0, 1, col = .rr$red, lty = 2, lwd = 1.5)
    axis(1, col = .rr$grid, col.ticks = .rr$soft)
    axis(2, col = .rr$grid, col.ticks = .rr$soft)
    mtext(sprintf("r = %.3f   RMSE = %.2f", s$correlation, s$rmse), 3,
          line = 0.2, cex = 0.8, col = .rr$soft)
  }
}

#' Simulate paired-comparison EFRM data with differing frame units
#'
#' Generates dichotomous paired comparisons whose latent unit differs across
#' judge-panel by object-set frames -- the paired-comparison extension of the
#' extended frame of reference model (Humphry 2005) fitted by
#' \code{\link{btl_efrm}}. Objects in set \code{s} have a within-set
#' calibration location \code{beta}; their common-scale value is
#' \code{v = alpha_s beta + kappa_s}. A comparison judged in panel \code{g}
#' carries the panel unit \code{phi_g}: within a set the comparison logit is
#' \code{phi_g (beta_a - beta_b)}, across sets it is \code{phi_g (v_a - v_b)}.
#' The planted panel units, set units and origins are recovered by
#' \code{\link{btl_efrm}}.
#'
#' @param n_objects_per_set,n_sets Objects in each set and number of sets.
#' @param n_judges_per_panel,n_panels Judges in each panel and number of panels.
#' @param reps_within Replications of each within-set object pair.
#' @param reps_cross Replications of each cross-set object pair.
#' @param panel_units Panel units \code{phi} (length \code{n_panels}); the
#'   default is all one, and any supplied vector is rescaled to geometric
#'   mean one.
#' @param set_units Set units \code{alpha} (length \code{n_sets}); the default
#'   is all one, and \code{alpha_1} is forced to one (the reference set).
#' @param set_origins Set origins \code{kappa} (length \code{n_sets}); the
#'   default is all zero, and \code{kappa_1} is forced to zero.
#' @param object_sd Spread of the within-set calibration locations.
#' @param seed Optional RNG seed.
#' @return A data frame of class \code{"rasch_sim"} with columns
#'   \code{object_a}, \code{object_b}, \code{winner}, \code{judge} and
#'   \code{panel}, and \code{attr(x, "truth")} holding the common-scale values
#'   \code{v}, the per-set \code{beta}, the units \code{phi}, \code{alpha},
#'   \code{kappa}, and the \code{object_sets} map to pass to
#'   \code{\link{btl_efrm}}.
#' @examples
#' d <- simulate_btl_efrm(6, 2, set_units = c(1, 1.4), seed = 1)
#' bt <- btl_efrm(d, "object_a", "object_b", winner = "winner",
#'                judge = "judge", panels = "panel",
#'                object_sets = attr(d, "truth")$object_sets)
#' bt$alpha_table   # recovers the ~1.4 set unit
#' @export
simulate_btl_efrm <- function(n_objects_per_set = 8, n_sets = 2,
                              n_judges_per_panel = 6, n_panels = 2,
                              reps_within = 20, reps_cross = 20,
                              panel_units = NULL, set_units = NULL,
                              set_origins = NULL, object_sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  S <- as.integer(n_sets); G <- as.integer(n_panels)
  Kp <- as.integer(n_objects_per_set); Jp <- as.integer(n_judges_per_panel)

  phi <- if (is.null(panel_units)) rep(1, G) else as.numeric(panel_units)
  if (length(phi) != G) stop("panel_units must have length n_panels")
  phi <- phi / exp(mean(log(phi)))                    # geometric mean one
  alpha <- if (is.null(set_units)) rep(1, S) else as.numeric(set_units)
  if (length(alpha) != S) stop("set_units must have length n_sets")
  alpha <- alpha / alpha[1]                            # alpha_1 = 1
  kappa <- if (is.null(set_origins)) rep(0, S) else as.numeric(set_origins)
  if (length(kappa) != S) stop("set_origins must have length n_sets")
  kappa <- kappa - kappa[1]                            # kappa_1 = 0

  set_nm <- sprintf("set%d", seq_len(S))
  panel_nm <- sprintf("panel%d", seq_len(G))
  objs_by_set <- lapply(seq_len(S), function(s) sprintf("S%dO%02d", s, seq_len(Kp)))
  names(objs_by_set) <- set_nm
  beta <- numeric(0)
  for (s in seq_len(S)) {
    bs <- as.numeric(scale(seq_len(Kp))) * object_sd   # sum-zero, spread object_sd
    beta <- c(beta, setNames(bs, objs_by_set[[s]]))
  }
  set_of <- setNames(rep(seq_len(S), each = Kp), unlist(objs_by_set))
  v <- alpha[set_of] * beta + kappa[set_of]

  judges <- sprintf("J%03d", seq_len(G * Jp))
  panel_of <- setNames(panel_nm[rep(seq_len(G), each = Jp)], judges)

  # assemble the object pairs (within each set, then across every set pair)
  aa <- bb <- character(0)
  for (s in seq_len(S)) {
    pr <- t(utils::combn(objs_by_set[[s]], 2))
    aa <- c(aa, rep(pr[, 1], reps_within)); bb <- c(bb, rep(pr[, 2], reps_within))
  }
  if (S > 1L) for (i in seq_len(S - 1L)) for (j in (i + 1L):S) {
    grid <- expand.grid(oa = objs_by_set[[i]], ob = objs_by_set[[j]],
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    aa <- c(aa, rep(grid$oa, reps_cross)); bb <- c(bb, rep(grid$ob, reps_cross))
  }
  R <- length(aa)
  jd <- sample(judges, R, replace = TRUE)
  pan <- panel_of[jd]
  # frame-dependent logit: within-set uses beta, cross-set uses the common v
  same <- set_of[aa] == set_of[bb]
  lp <- ifelse(same, phi[match(pan, panel_nm)] * (beta[aa] - beta[bb]),
               phi[match(pan, panel_nm)] * (v[aa] - v[bb]))
  win_a <- stats::runif(R) < stats::plogis(lp)
  d <- data.frame(object_a = aa, object_b = bb,
                  winner = ifelse(win_a, aa, bb),
                  judge = jd, panel = unname(pan),
                  stringsAsFactors = FALSE)
  rownames(d) <- NULL

  planted <- character(0)
  if (any(phi != 1)) planted <- c(planted, sprintf(
    "panel units phi = (%s)", paste(sprintf("%.2f", phi), collapse = ", ")))
  if (any(alpha != 1)) planted <- c(planted, sprintf(
    "set units alpha = (%s)", paste(sprintf("%.2f", alpha), collapse = ", ")))
  if (any(kappa != 0)) planted <- c(planted, sprintf(
    "set origins kappa = (%s)", paste(sprintf("%.2f", kappa), collapse = ", ")))
  if (!length(planted)) planted <- "equal units (phi = alpha = 1)"

  attr(d, "truth") <- list(
    layout = "btl_efrm",
    description = sprintf("%d objects (%d sets) x %d panels, %d comparisons",
                          S * Kp, S, G, R),
    v = v, beta = beta, phi = setNames(phi, panel_nm),
    alpha = setNames(alpha, set_nm), kappa = setNames(kappa, set_nm),
    set_of = set_of, object_sets = objs_by_set, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}
