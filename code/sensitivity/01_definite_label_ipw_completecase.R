#!/usr/bin/env Rscript
# Sensitivity analysis (post-lock).
# Preregistration: https://doi.org/10.17605/OSF.IO/CNXS3
#
# Estimand: covariate-standardized 15-year graduation risk difference
# (women minus men), 2005-2009 primary subcohort (>=15y follow-up).
# Model: cause-specific Royston-Parmar (rstpm2), baseline df = 5,
# gender x log(time) interaction df = 3 (graduation cause); cessation
# (>=5y publication gap, absorbing) as competing event, own RP model
# (df = 5, no tvc). G-computation of CIF_grad(15y) under female = 1 vs 0
# over the observed covariate distribution; institution-cluster bootstrap
# (B = 500, outer) around calibrated MI (inner); percentile CI.
# Calibration (isotonic vs Platt, locked Brier rule) fitted ONCE on the
# locked reference standard, held fixed across bootstrap reps.
#
# COMPUTE DEVIATION (deviation_log #17): inner imputations per bootstrap
# rep M_b = 2 (prereg inner M = 20-50 retained for the observed point
# estimate, reported with Rubin pooling). Gender uncertainty propagates
# across the 500 bootstrap reps; the observed-data estimator is unchanged.
# Warm-started fits (init from observed fit) + one-shot grid prediction.
#
# Resumable: per-rep results append to data/bootstrap07_results.csv;
# deterministic per-rep seeds (902607 + b*1000 + m) make reruns identical.

suppressMessages({library(rstpm2); library(data.table); library(yaml); library(parallel)})

cfg <- yaml.load_file("config/config.yaml")
stopifnot(isTRUE(cfg$annotation_lock), file.exists("data/annotation_lock.timestamp"))

MASK_LABELS <- Sys.getenv("MASK_LABELS", "0") == "1"
B      <- as.integer(Sys.getenv("BOOT_B", "500"))  # 500
M_OBS  <- as.integer(Sys.getenv("M_OBS", if (MASK_LABELS) "20" else "1"))
M_BOOT <- as.integer(Sys.getenv("BOOT_M", 2))  # inner M per bootstrap rep (deviation #17)
DFB    <- cfg$model$baseline_df          # 5
DFT    <- cfg$model$gender_time_interaction_df  # 3
FAIL_T <- cfg$model$convergence_fail_threshold  # 0.05
TMAX   <- 15; GRID <- 50
RESF   <- "data/bootstrap07_results.csv"
CORES  <- 1L  # runtime-only serial adaptation for Windows
TS     <- seq(0.05, TMAX, length.out = GRID)

## Data
d0 <- fread("data/analysis_cohort.csv")
d0[, female_gold := fifelse(gender_gold == "F", 1,
                    fifelse(gender_gold == "M", 0, NA_real_))]
d0[, ev1 := as.integer(event_type == 1)]   # graduation first
d0[, ev2 := as.integer(event_type == 2)]   # cessation first (competing)
stopifnot(nrow(d0) == 1783L,
          d0[primary_cohort == TRUE, sum(ev1)] == 866L,
          d0[primary_cohort == TRUE, sum(ev2)] == 534L,
          d0[primary_cohort == TRUE, sum(ev1 == 0L & ev2 == 0L)] == 136L)

## fixed transforms & factor levels (computed once from full primary cohort)
d0[, `:=`(lp = log1p(preentry_pubs), lc = log1p(preentry_citations),
          lh = preentry_hindex,     lco = log1p(preentry_coauthors),
          li = intl_collab_share)]
d0[, no_preentry_work := as.integer(preentry_pubs == 0L)]
d0[, country_incomplete := as.integer(preentry_pubs > 0L & is.na(intl_collab_share))]
ctr <- data.frame(mu = sapply(d0[, .(lp, lc, lh, lco, li)], mean, na.rm = TRUE),
                  sd = sapply(d0[, .(lp, lc, lh, lco, li)], sd, na.rm = TRUE))
stopifnot(all(is.finite(ctr$mu)), all(is.finite(ctr$sd)), all(ctr$sd > 0))
scale_cols <- rownames(ctr)
apply_scale <- function(dd) {
  for (cc in scale_cols) set(dd, j = cc, value = (dd[[cc]] - ctr[cc, "mu"]) / ctr[cc, "sd"])
  set(dd, j = "li", value = ifelse(is.na(dd$li), 0, dd$li))
  dd
}
apply_scale(d0)
COUNTRY_TOPN <- as.integer(Sys.getenv("COUNTRY_TOPN", "3"))
stopifnot(length(COUNTRY_TOPN) == 1L, is.finite(COUNTRY_TOPN), COUNTRY_TOPN >= 1L)
TOPC_ALL <- setdiff(names(sort(-table(d0$inst_country))), "未知")
TOPC  <- TOPC_ALL[seq_len(min(COUNTRY_TOPN, length(TOPC_ALL)))]
CLEV  <- c(TOPC, "Other"); TLEV <- c("高", "中", "低", "缺失")
EYLEV <- sort(unique(d0[primary_cohort == TRUE, entry_year])); JLEV <- sort(unique(d0$index_journal))
fix_levels <- function(dd) {
  dd[, ctry := factor(ifelse(inst_country %in% TOPC, inst_country, "Other"), levels = CLEV)]
  dd[, inst_tier := factor(inst_tier, levels = TLEV)]
  dd[, ey := factor(entry_year, levels = EYLEV)]
  dd[, index_journal := factor(index_journal, levels = JLEV)]
  dd
}
fix_levels(d0)
d0[, cluster := fifelse(is.na(inst_id) | inst_id == "", paste0("solo_", fa_id), inst_id)]
MODC <- c("female", "ey", "index_journal", "inst_tier", "ctry", "academic_age",
          "lp", "lc", "lh", "lco", "li", "no_preentry_work", "country_incomplete", "time_years")

## Calibration and exact core design weights
source("code/sensitivity/weighted_calibration.R")
rs <- d0[core_selected == 1 & !is.na(p_female) & female_gold %in% c(0, 1)]
stopifnot(nrow(rs) == 575L, !anyNA(rs$design_weight), all(rs$design_weight > 0))
cal <- fit_calibration(rs)
brier_iso <- cal$brier_isotonic
brier_platt <- cal$brier_platt
use_iso <- cal$selected == "isotonic"
p_calibrate <- cal$predict
prev_rs <- cal$prevalence
cat(sprintf("calibration: weighted isotonic Brier=%.4f  weighted Platt Brier=%.4f  -> %s | core RS n=%d, weighted prevalence=%.3f
",
            brier_iso, brier_platt, ifelse(use_iso, "ISOTONIC", "PLATT"), nrow(rs), prev_rs))
d0[, p_cal := fifelse(!is.na(p_female), p_calibrate(p_female), prev_rs)]
if (MASK_LABELS) {
  d0[, gender_gold := NA_character_]
  d0[, female_gold := NA_real_]
}

## Core probability-sample IPW and unweighted complete-case analysis.
d0[, w := fifelse(core_selected == 1, design_weight, NA_real_)]
cc_summary <- list(N_full = nrow(d0),
                   N_core_full = sum(d0$core_selected == 1),
                   N_core_primary = sum(d0$core_selected == 1 & d0$primary_cohort == TRUE),
                   weighting = "Exact inclusion-probability design weights; no fitted selection model")
COMPLETE_GENDER_ONLY <- Sys.getenv("COMPLETE_GENDER_ONLY", "0") == "1"
USE_WEIGHTS <- Sys.getenv("USE_WEIGHTS", "1") == "1"
if (COMPLETE_GENDER_ONLY) {
  d0 <- d0[core_selected == 1 & primary_cohort == TRUE & female_gold %in% c(0, 1)]
  if (!USE_WEIGHTS) d0[, w := 1]
} else {
  d0 <- d0[core_selected == 1 & primary_cohort == TRUE]
}
stopifnot(nrow(d0) == if (COMPLETE_GENDER_ONLY) 539L else 689L)
RESULT_TAG <- if (MASK_LABELS) "core_masked" else if (COMPLETE_GENDER_ONLY && USE_WEIGHTS) "core_ipw" else if (COMPLETE_GENDER_ONLY) "core_completecase" else if (USE_WEIGHTS) "core_ipw_full" else "core_unweighted_full"
RESULT_TAG <- paste0(RESULT_TAG, Sys.getenv("RESULT_SUFFIX", ""))
cc_summary$country_top_n <- COUNTRY_TOPN
cc_summary$country_levels <- CLEV
cc_summary$complete_gender_only <- COMPLETE_GENDER_ONLY
cc_summary$N <- nrow(d0)
cc_summary$weight_mean <- mean(d0$w)
cc_summary$weight_min <- min(d0$w)
cc_summary$weight_max <- max(d0$w)
cc_summary$weight_ESS <- sum(d0$w)^2 / sum(d0$w^2)
cc_summary$event1 <- sum(d0$ev1 == 1)
cc_summary$event2 <- sum(d0$ev2 == 1)
cc_summary$censored <- sum(d0$ev1 == 0 & d0$ev2 == 0)

## Model and g-computation
RHS <- "female + ey + index_journal + inst_tier + ctry + academic_age + lp + lc + lh + lco + li + no_preentry_work"
F1 <- as.formula(paste("Surv(time_years, ev1) ~", RHS))
F2 <- as.formula(paste("Surv(time_years, ev2) ~", RHS))

## H(t|x,female=g) for g in {0,1} via one big predict (female=0) + delta(t) trick
## (female enters only main effect + tvc basis, so the F-vs-M difference in the
##  linear predictor is the same for every subject at a given t)
haz_grids <- function(f1, f2, dd) {
  n <- nrow(dd); G <- length(TS)
  nd0 <- dd[, ..MODC]; nd0[, female := 0]
  long <- nd0[rep(1:n, times = G)]; long[, time_years := rep(TS, each = n)]
  H10 <- matrix(predict(f1, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  H20 <- matrix(predict(f2, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  one <- nd0[rep(1, 2 * G)]; one[, time_years := rep(TS, times = 2)]
  one[, female := rep(c(0, 1), each = G)]
  d1 <- predict(f1, newdata = one, type = "link")[(G + 1):(2 * G)] -
        predict(f1, newdata = one, type = "link")[1:G]
  d2 <- predict(f2, newdata = one, type = "link")[(G + 1):(2 * G)] -
        predict(f2, newdata = one, type = "link")[1:G]
  H11 <- sweep(H10, 2, exp(d1), `*`); H21 <- sweep(H20, 2, exp(d2), `*`)
  list(H1F = H11, H1M = H10, H2F = H21, H2M = H20)
}
cif_from_H <- function(H1, H2) {           # midpoint-S (2nd-order; audit 2026-09-07)
  G <- ncol(H1)
  S  <- exp(-(H1 + H2))
  Smid <- (cbind(1, S[, -G]) + S) / 2
  h1 <- t(apply(H1, 1, function(x) diff(c(0, x))))
  mean(rowSums(h1 * Smid))
}
fit_gcomp <- function(dd, init1 = NULL, init2 = NULL) {
  f1 <- if (USE_WEIGHTS) {
    stpm2(F1, data = dd, df = DFB, tvc = list(female = DFT), init = init1, weights = w)
  } else {
    stpm2(F1, data = dd, df = DFB, tvc = list(female = DFT), init = init1)
  }
  f2 <- if (USE_WEIGHTS) {
    stpm2(F2, data = dd, df = DFB, init = init2, weights = w)
  } else {
    stpm2(F2, data = dd, df = DFB, init = init2)
  }
  hg <- haz_grids(f1, f2, dd)
  rF <- cif_from_H(hg$H1F, hg$H2F); rM <- cif_from_H(hg$H1M, hg$H2M)
  list(c(riskF = rF, riskM = rM, RD = rF - rM, RR = rF / rM),
       init1 = coef(f1), init2 = coef(f2))
}
impute_one <- function(dd, m, seed) {
  set.seed(seed)
  dd[, female := fifelse(!is.na(female_gold), female_gold, rbinom(.N, 1, p_cal))]
  dd
}
pool_est <- function(ests) {
  em <- do.call(rbind, ests)
  c(riskF = mean(em[, "riskF"]), riskM = mean(em[, "riskM"]),
    RD = mean(em[, "RD"]), RR = exp(mean(log(em[, "RR"]))))
}
source("code/sensitivity/fixed_design_fallback.R")

## observed-data fits (also provide warm-start inits for bootstrap)
RESULT_BASE <- file.path("data", paste0("results_", RESULT_TAG))
SKIP_OBS <- Sys.getenv("SKIP_OBS", "0") == "1" && file.exists(paste0(RESULT_BASE, ".rds"))
obs0 <- fit_gcomp(impute_one(copy(d0), 0, seed = 902607))
INIT1 <- obs0$init1; INIT2 <- obs0$init2
obs_fixed0 <- fit_gcomp_fixed(impute_one(copy(d0), 0, seed = 902607),
                              init1 = INIT1, init2 = INIT2)
stopifnot(max(abs(obs0[[1]] - obs_fixed0[[1]])) < 1e-4)
FIXED_INIT1 <- obs_fixed0$init1; FIXED_INIT2 <- obs_fixed0$init2
cat(sprintf("fixed-design equivalence check: max absolute difference=%.3g\n",
            max(abs(obs0[[1]] - obs_fixed0[[1]]))))
if (SKIP_OBS) {
  obs_est <- readRDS(paste0(RESULT_BASE, ".rds"))$observed
} else {
  obs_cells <- lapply(1:M_OBS, function(m)
    tryCatch(fit_gcomp(impute_one(copy(d0), m, seed = 902607 + m))[[1]],
             error = function(e) NULL))
  obs_ok <- obs_cells[!vapply(obs_cells, is.null, logical(1))]
  obs_est <- pool_est(obs_ok)
  cat(sprintf("observed (M=%d, %d ok): riskF %.4f | riskM %.4f | RD %+.2f pp | RR %.3f\n",
              M_OBS, length(obs_ok), obs_est["riskF"], obs_est["riskM"],
              100 * obs_est["RD"], obs_est["RR"]))
}

## Institution-cluster bootstrap
clusters <- split(seq_len(nrow(d0)), d0$cluster); CL <- names(clusters)
boot_rep <- function(b) {
  set.seed(902607 + b * 1000)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- copy(d0[idx]); apply_scale(db); fix_levels(db)
  ests <- lapply(1:M_BOOT, function(m) {
    fallback <- 0L
    est <- tryCatch(
      fit_gcomp(impute_one(copy(db), m, seed = 902607 + b * 1000 + m),
                init1 = INIT1, init2 = INIT2)[[1]],
      error = function(e) {
        fallback <<- 1L
        tryCatch(
          fit_gcomp_fixed(impute_one(copy(db), m, seed = 902607 + b * 1000 + m),
                          init1 = FIXED_INIT1, init2 = FIXED_INIT2)[[1]],
          error = function(e2) NULL
        )
      }
    )
    if (is.null(est) || any(!is.finite(est))) return(NULL)
    list(est = est, fallback = fallback)
  })
  ok <- !vapply(ests, is.null, logical(1))
  if (!any(ok)) return(c(b = b, nfail = M_BOOT, riskF = NA, riskM = NA, RD = NA, RR = NA, nfallback = M_BOOT))
  nfallback <- sum(vapply(ests[ok], function(z) z$fallback, numeric(1)))
  c(b = b, nfail = M_BOOT - sum(ok),
    as.list(pool_est(lapply(ests[ok], function(z) z$est))), nfallback = nfallback)
}

PART_DIR <- file.path("work_tmp", paste0(RESULT_TAG, "_parts"))
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)
read_parts <- function() {
  pf <- list.files(PART_DIR, pattern = "\\.csv$", full.names = TRUE)
  if (!length(pf)) return(NULL)
  parts <- rbindlist(lapply(pf, function(f) tryCatch(fread(f), error = function(e) NULL)), fill = TRUE)
  stopifnot(!anyDuplicated(parts$b))
  parts
}
old_parts <- read_parts()
done_b <- if (is.null(old_parts)) integer(0) else old_parts$b
todo <- setdiff(1:B, done_b)
cat(sprintf("bootstrap: %d/%d done, %d to go (%d cores)\n", length(done_b), B, length(todo), CORES))
flush.console()
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  t0 <- Sys.time()
  res <- as.data.table(do.call(rbind, mclapply(part, boot_rep, mc.cores = CORES)))
  pf <- file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1]))
  fwrite(res, pf)
  cat(sprintf("  reps %d-%d done in %.0fs | total %d/%d\n", part[1], part[length(part)],
              difftime(Sys.time(), t0, units = "secs"),
              length(done_b) + max(0, s - 1) + length(part), B))
  flush.console()
}

## ---------------- pooling & output ----------------
br_all <- read_parts()
stopifnot(!is.null(br_all), nrow(br_all) > 0)
br <- br_all[!is.na(RD)]
stopifnot(nrow(br) > 0)
fail_rate <- sum(br_all$nfail) / (nrow(br_all) * M_BOOT)
df_note <- if (fail_rate > FAIL_T) "GUARD TRIGGERED: >5% cells failed; refit df=3 as primary" else ""
ci <- function(x) unname(quantile(x, c(.025, .975), na.rm = TRUE))
rdci <- ci(br$RD); rrci <- exp(ci(log(br$RR)))
gate_open <- (rdci[1] > 0) || (rdci[2] < 0)
out <- list(
  estimand = paste0(RESULT_TAG, ": standardized 15y graduation RD (F - M), core probability sample"),
  analyst = "ZZ H",
  n_primary = nrow(d0),
  calibration = ifelse(use_iso, "isotonic", "platt"),
  brier_isotonic = brier_iso,
  brier_platt = brier_platt,
  observed = list(riskF = unname(obs_est["riskF"]), riskM = unname(obs_est["riskM"]),
                  RD = unname(obs_est["RD"]), RR = unname(obs_est["RR"])),
  RD_ci = rdci,
  RR_ci = rrci,
  riskF_ci = ci(br$riskF),
  riskM_ci = ci(br$riskM),
  bootstrap_requested = B,
  bootstrap_reps = nrow(br),
  bootstrap_failed_reps = nrow(br_all) - nrow(br),
  fixed_design_fallback_reps = if ("nfallback" %in% names(br_all)) sum(br_all$nfallback > 0, na.rm = TRUE) else 0,
  fixed_design_max_abs_difference = max(abs(obs0[[1]] - obs_fixed0[[1]])),
  M_obs = M_OBS,
  M_boot = M_BOOT,
  cell_failure_rate = fail_rate,
  df_note = df_note,
  gate_open = gate_open,
  diagnostics = c(cc_summary, list(label_masking = MASK_LABELS,
                                   calibration_fit_N = nrow(rs),
                                   calibration_fit_before_masking = MASK_LABELS))
)
saveRDS(out, paste0(RESULT_BASE, ".rds"))
fwrite(br_all, file.path("data", paste0("bootstrap_", RESULT_TAG, ".csv")))
jsonlite::write_json(out, paste0(RESULT_BASE, ".json"), auto_unbox = TRUE,
                     pretty = TRUE, digits = NA, na = "null")
cat(readLines(paste0(RESULT_BASE, ".json")), sep = "\n")
