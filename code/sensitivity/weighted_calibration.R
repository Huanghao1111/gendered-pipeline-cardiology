# Weighted isotonic and Platt calibration for the core probability sample.
# Analyst: ZZ H

weighted_brier <- function(y, p, w) {
  stopifnot(length(y) == length(p), length(y) == length(w), all(is.finite(w)), all(w > 0))
  sum(w * (p - y)^2) / sum(w)
}

weighted_isoreg <- function(x, y, w) {
  stopifnot(length(x) == length(y), length(x) == length(w))
  stopifnot(all(is.finite(x)), all(is.finite(y)), all(is.finite(w)), all(w > 0))
  ord <- order(x, -y)
  x <- x[ord]; y <- y[ord]; w <- w[ord]
  n <- length(x)
  means <- numeric(n); weights <- numeric(n); counts <- integer(n)
  starts <- integer(n); ends <- integer(n); nb <- 0L
  for (i in seq_len(n)) {
    nb <- nb + 1L
    means[nb] <- y[i]; weights[nb] <- w[i]; counts[nb] <- 1L
    starts[nb] <- i; ends[nb] <- i
    while (nb > 1L && means[nb - 1L] > means[nb]) {
      new_w <- weights[nb - 1L] + weights[nb]
      new_mean <- (weights[nb - 1L] * means[nb - 1L] + weights[nb] * means[nb]) / new_w
      means[nb - 1L] <- new_mean; weights[nb - 1L] <- new_w
      counts[nb - 1L] <- counts[nb - 1L] + counts[nb]
      ends[nb - 1L] <- ends[nb]
      nb <- nb - 1L
    }
  }
  means <- means[seq_len(nb)]; counts <- counts[seq_len(nb)]
  ends <- ends[seq_len(nb)]
  yf <- rep(means, counts)
  list(x = x, y = y, w = w, yf = yf, knots = ends)
}

isoreg_predictor <- function(fit) {
  n <- length(fit$yf)
  stats::stepfun(x = fit$x[fit$knots], y = c(fit$yf[fit$knots], fit$yf[n]), right = TRUE)
}

platt_fit <- function(x, y, w) {
  suppressWarnings(stats::glm(y ~ qlogis(pmin(pmax(x, 1e-6), 1 - 1e-6)),
                              family = binomial(), weights = w))
}

platt_predictor <- function(fit) {
  function(p) {
    pp <- pmin(pmax(p, 0), 1)
    as.numeric(stats::predict(fit, newdata = data.frame(x = pp), type = "response"))
  }
}

fit_calibration <- function(dt) {
  stopifnot(all(c("p_female", "female_gold", "design_weight") %in% names(dt)))
  x <- dt$p_female; y <- dt$female_gold; w <- dt$design_weight
  iso <- weighted_isoreg(x, y, w)
  iso_raw <- isoreg_predictor(iso)
  iso_f <- function(p) {
    pp <- pmin(pmax(p, min(x, na.rm = TRUE)), max(x, na.rm = TRUE))
    iso_raw(pp)
  }
  platt <- platt_fit(x, y, w)
  platt_f <- platt_predictor(platt)
  b_iso <- weighted_brier(y, iso_f(x), w)
  b_platt <- weighted_brier(y, platt_f(x), w)
  use_iso <- b_iso <= b_platt
  list(selected = if (use_iso) "isotonic" else "platt",
       predict = if (use_iso) iso_f else platt_f,
       brier_isotonic = b_iso, brier_platt = b_platt,
       prevalence = sum(w * y) / sum(w), iso = iso, platt = platt)
}
