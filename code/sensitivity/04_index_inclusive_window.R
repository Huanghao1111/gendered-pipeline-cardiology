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
# (B = 1000, outer) around calibrated MI (inner); percentile CI.
# Calibration (isotonic vs Platt, locked Brier rule) fitted ONCE on the
# locked reference standard, held fixed across bootstrap reps.
#
# The accepted covariate-window analyses used one inner imputation per
# bootstrap replicate (BOOT_M=1); observed estimates retain 20 imputations.
# Warm-started fits (init from observed fit) + one-shot grid prediction.
#
# Resumable: analysis-specific per-rep chunks use deterministic seeds
# (902607 + b*1000 + m), so interrupted runs can be resumed safely.

suppressMessages({library(rstpm2); library(data.table); library(yaml); library(parallel)})

cfg <- yaml.load_file("config/config.yaml")
stopifnot(isTRUE(cfg$annotation_lock), file.exists("data/annotation_lock.timestamp"))

B      <- as.integer(Sys.getenv("BOOT_B", cfg$model$bootstrap_B))  # 1000
M_OBS  <- min(cfg$model$imputations_M)   # 20 (observed estimate)
M_BOOT <- as.integer(Sys.getenv("BOOT_M", 1))  # accepted window-analysis setting
DFB    <- cfg$model$baseline_df          # 5
DFT    <- cfg$model$gender_time_interaction_df  # 3
FAIL_T <- cfg$model$convergence_fail_threshold  # 0.05
TMAX   <- 15; GRID <- 50
RESF   <- "data/bootstrap07_results.csv"
CORES  <- 1L  # runtime-only serial adaptation for Windows
TS     <- seq(0.05, TMAX, length.out = GRID)

## Data
d0 <- fread("data/analysis_cohort_index_inclusive.csv")[primary_cohort == TRUE]
d0[, female_gold := fifelse(gender_gold == "F", 1,
                    fifelse(gender_gold == "M", 0, NA_real_))]
d0[, ev1 := as.integer(event_type == 1)]   # graduation first
d0[, ev2 := as.integer(event_type == 2)]   # cessation first (competing)
stopifnot(nrow(d0) == 1536L, sum(d0$ev1) == 866L, sum(d0$ev2) == 534L,
          sum(d0$ev1 == 0L & d0$ev2 == 0L) == 136L)

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
TOPC  <- setdiff(names(sort(-table(d0$inst_country))), "未知")[1:8]
CLEV  <- c(TOPC, "Other"); TLEV <- c("高", "中", "低", "缺失")
EYLEV <- sort(unique(d0$entry_year)); JLEV <- sort(unique(d0$index_journal))
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

## Calibration of the core probability sample
source("code/sensitivity/weighted_calibration.R")
rs <- d0[core_selected == 1 & !is.na(p_female) & female_gold %in% c(0, 1)]
stopifnot(nrow(rs) == 495L, !anyNA(rs$design_weight), all(rs$design_weight > 0))
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

## Model and g-computation
## no_preentry_work is structurally zero for everyone after including the entry work.
## It is therefore non-estimable and omitted only in this window variant.
stopifnot(all(d0$no_preentry_work == 0L))
RHS <- "female + ey + index_journal + inst_tier + ctry + academic_age + lp + lc + lh + lco + li"
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
  f1 <- stpm2(F1, data = dd, df = DFB, tvc = list(female = DFT), init = init1)
  f2 <- stpm2(F2, data = dd, df = DFB, init = init2)
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

## observed-data fits (also provide warm-start inits for bootstrap)
SKIP_OBS <- Sys.getenv("SKIP_OBS", "0") == "1" && file.exists("data/results_index_inclusive.rds")
obs0 <- fit_gcomp(impute_one(copy(d0), 0, seed = 902607))
INIT1 <- obs0$init1; INIT2 <- obs0$init2
if (SKIP_OBS) {
  obs_est <- readRDS("data/results_index_inclusive.rds")$observed
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
    tryCatch(fit_gcomp(impute_one(db, m, seed = 902607 + b * 1000 + m),
                       init1 = INIT1, init2 = INIT2)[[1]],
             error = function(e) NULL)
  })
  ok <- !vapply(ests, is.null, logical(1))
  if (!any(ok)) return(c(b = b, nfail = M_BOOT, riskF = NA, riskM = NA, RD = NA, RR = NA))
  c(b = b, nfail = M_BOOT - sum(ok), as.list(pool_est(ests[ok])))
}

dir.create("data/bootstrap_index_inclusive_parts", showWarnings = FALSE)
dir.create("work_tmp/bootstrap_index_inclusive_parts", showWarnings = FALSE)
# seed /tmp parts from the mount (survives only if copied back)
for (f in list.files("data/bootstrap_index_inclusive_parts", pattern = "\\.csv$", full.names = TRUE))
  file.copy(f, file.path("work_tmp/bootstrap_index_inclusive_parts", basename(f)), overwrite = TRUE)
read_parts <- function() {
  pf <- list.files("work_tmp/bootstrap_index_inclusive_parts", pattern = "\\.csv$", full.names = TRUE)
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
  pf <- sprintf("work_tmp/bootstrap_index_inclusive_parts/chunk_%04d.csv", part[1])
  fwrite(res, pf)
  mf <- sprintf("data/bootstrap_index_inclusive_parts/chunk_%04d.csv", part[1])
  copied <- FALSE
  for (att in 1:6) {
    tryCatch({ file.copy(pf, mf, overwrite = TRUE)
               copied <- file.exists(mf) && file.size(mf) == file.size(pf) }, error = function(e) NULL)
    if (copied) break; Sys.sleep(2 * att)
  }
  if (!copied) cat(sprintf("  WARN: mount copy of chunk %04d failed (kept in /tmp)\n", part[1]))
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
ci <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
rdci <- ci(br$RD); rrci <- exp(ci(log(br$RR)))
gate_open <- (rdci[1] > 0) || (rdci[2] < 0)
out <- list(estimand = "standardized 15y graduation RD (F - M), primary cohort 2005-2009",
            n_primary = nrow(d0), calibration = ifelse(use_iso, "isotonic", "platt"),
            brier_isotonic = brier_iso, brier_platt = brier_platt,
            observed = obs_est, RD_ci = rdci, RR_ci = rrci,
            riskF_ci = ci(br$riskF), riskM_ci = ci(br$riskM),
            bootstrap_reps = nrow(br), bootstrap_failed_reps = nrow(br_all) - nrow(br),
            M_obs = M_OBS, M_boot = M_BOOT,
            cell_failure_rate = fail_rate, df_note = df_note, gate_open = gate_open)
saveRDS(out, "data/results_index_inclusive.rds")
json <- sprintf(paste0(
 '{\n  "estimand": "%s",\n  "n_primary": %d,\n  "calibration": "%s",\n  "calibration_brier": {"isotonic": %.4f, "platt": %.4f},\n',
 '  "observed": {"riskF_15y": %.4f, "riskM_15y": %.4f, "RD_pp": %.2f, "RR": %.3f},\n',
 '  "RD_pp_ci95": [%.2f, %.2f],\n  "riskF_ci95": [%.4f, %.4f],\n  "riskM_ci95": [%.4f, %.4f],\n',
 '  "RR_ci95": [%.3f, %.3f],\n  "bootstrap_reps": %d, "bootstrap_failed_reps": %d, "M_obs": %d, "M_boot": %d, "cell_failure_rate": %.4f,\n',
 '  "df_note": "%s",\n  "secondary_gate": "%s"\n}'),
 out$estimand, out$n_primary, out$calibration, brier_iso, brier_platt,
 out$observed["riskF"], out$observed["riskM"], 100 * out$observed["RD"], out$observed["RR"],
 100 * rdci[1], 100 * rdci[2], out$riskF_ci[1], out$riskF_ci[2],
 out$riskM_ci[1], out$riskM_ci[2], rrci[1], rrci[2],
 out$bootstrap_reps, out$bootstrap_failed_reps, M_OBS, M_BOOT, fail_rate, df_note,
 ifelse(gate_open, "OPEN", "CLOSED"))
writeLines(json, "data/results_index_inclusive.json")
cat(json, "\n")
