# Equivalent fixed-design fallback for sparse bootstrap cells.
# Requires objects defined by 02_ipw_completecase.R.
make_fixed_design <- function(dd) {
  mm <- model.matrix(~ female + ey + index_journal + inst_tier + ctry +
                       academic_age + lp + lc + lh + lco + li + no_preentry_work,
                     data = dd)
  if (!exists("FIXED_KEEP", envir = globalenv())) {
    candidate <- colnames(mm)[-1]
    FIXED_KEEP <<- candidate[colSums(abs(mm[, candidate, drop = FALSE])) > 0]
  }
  original_names <- FIXED_KEEP
  xcols <- paste0("X", seq_along(original_names))
  for (j in seq_along(xcols)) set(dd, j = xcols[j], value = mm[, original_names[j]])
  attr(dd, "female_x") <- xcols[match("female", original_names)]
  attr(dd, "xcols") <- xcols
  dd
}

fit_gcomp_fixed <- function(dd, init1 = NULL, init2 = NULL) {
  dd <- make_fixed_design(dd)
  xcols <- attr(dd, "xcols")
  female_x <- attr(dd, "female_x")
  use_weights <- exists("USE_WEIGHTS", envir = globalenv()) && isTRUE(get("USE_WEIGHTS", envir = globalenv()))
  rhs <- paste(xcols, collapse = " + ")
  f1x <- as.formula(paste("Surv(time_years, ev1) ~", rhs))
  f2x <- as.formula(paste("Surv(time_years, ev2) ~", rhs))
  tvc_x <- setNames(list(DFT), female_x)
  f1 <- if (use_weights) {
    stpm2(f1x, data = dd, df = DFB, tvc = tvc_x, init = init1, weights = w)
  } else {
    stpm2(f1x, data = dd, df = DFB, tvc = tvc_x, init = init1)
  }
  f2 <- if (use_weights) {
    stpm2(f2x, data = dd, df = DFB, init = init2, weights = w)
  } else {
    stpm2(f2x, data = dd, df = DFB, init = init2)
  }
  n <- nrow(dd)
  G <- length(TS)
  nd0 <- dd[, c(xcols, "time_years"), with = FALSE]
  set(nd0, j = female_x, value = 0)
  long <- nd0[rep(1:n, times = G)]
  long[, time_years := rep(TS, each = n)]
  H10 <- matrix(predict(f1, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  H20 <- matrix(predict(f2, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  one <- nd0[rep(1, 2 * G)]
  one[, time_years := rep(TS, times = 2)]
  set(one, j = female_x, value = rep(c(0, 1), each = G))
  d1 <- predict(f1, newdata = one, type = "link")[(G + 1):(2 * G)] -
        predict(f1, newdata = one, type = "link")[1:G]
  d2 <- predict(f2, newdata = one, type = "link")[(G + 1):(2 * G)] -
        predict(f2, newdata = one, type = "link")[1:G]
  H11 <- sweep(H10, 2, exp(d1), `*`)
  H21 <- sweep(H20, 2, exp(d2), `*`)
  rF <- cif_from_H(H11, H21)
  rM <- cif_from_H(H10, H20)
  list(c(riskF = rF, riskM = rM, RD = rF - rM, RR = rF / rM),
       init1 = coef(f1), init2 = coef(f2))
}
