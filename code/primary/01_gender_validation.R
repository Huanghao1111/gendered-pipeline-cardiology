#!/usr/bin/env Rscript
# Dual-layer gender validation: weighted core calibration and special strata.
# Analyst: ZZ H
suppressMessages({library(data.table); library(jsonlite)})
source("code/primary/weighted_calibration.R")
SEED <- 20261112L; K <- 10L; REPEATS <- 10L; THRESHOLD <- 0.90
OUTDIR <- file.path("data", "gender_validation")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
d0 <- fread("data/analysis_cohort.csv")
stopifnot(!anyDuplicated(d0$fa_id), nrow(d0) == 1783L,
          sum(d0$core_selected == 1) == 800L, sum(d0$special_selected == 1) == 239L,
          sum(d0$core_selected == 1 & d0$special_selected == 1) == 104L)
origin_fun <- function(x) fifelse(is.na(x), "U",
                          fifelse(x %in% c("CN", "KR", "JP"), x,
                          fifelse(x == "other", "other", "U")))
d0[, origin_group := origin_fun(final_name_origin)]
europe <- c("AT", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GB", "GR", "HR", "HU", "IE", "IS", "IT", "LT", "LU", "LV", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK")
east_asia <- c("CN", "KR", "JP", "TW", "HK")
region_fun <- function(country) {
  country <- toupper(trimws(country))
  fifelse(is.na(country) | country == "" | country == "未知", "Unknown",
          fifelse(country %in% c("US", "CA"), "North America",
                  fifelse(country %in% europe, "Europe",
                          fifelse(country %in% east_asia, "East Asia", "Other"))))
}
d0[, inst_region := region_fun(inst_country)]
d0[, y := fifelse(gender_gold == "F", 1, fifelse(gender_gold == "M", 0, NA_real_))]
core <- d0[core_selected == 1]
special <- d0[special_selected == 1]
cal <- core[!is.na(p_female) & !is.na(y) & !is.na(design_weight)]
stopifnot(nrow(core) == 800L, nrow(cal) == 575L,
          sum(core$gender_gold %in% c("F", "M")) == 629L,
          all(cal$design_weight > 0))

weighted_metrics <- function(dt, pvar) {
  p <- dt[[pvar]]; y <- dt$y; w <- dt$design_weight
  keep <- !is.na(p) & !is.na(y) & !is.na(w)
  p <- p[keep]; y <- y[keep]; w <- w[keep]
  high <- p >= THRESHOLD | p <= 1 - THRESHOLD
  correct <- (p >= THRESHOLD & y == 1) | (p <= 1 - THRESHOLD & y == 0)
  wm <- function(x, ww = w) if (sum(ww) == 0) NA_real_ else sum(ww * x) / sum(ww)
  list(N_gold = nrow(dt), N_score = length(p), weight_sum = sum(w),
       score_coverage = length(p) / nrow(dt),
       N_high_call = sum(high), high_call_rate = wm(high),
       high_accuracy = if (any(high)) wm(correct[high], w[high]) else NA_real_,
       sensitivity = if (any(y == 1)) wm((p >= THRESHOLD)[y == 1], w[y == 1]) else NA_real_,
       specificity = if (any(y == 0)) wm((p < THRESHOLD)[y == 0], w[y == 0]) else NA_real_,
       Brier = weighted_brier(y, p, w))
}

apparent_fit <- fit_calibration(cal)
cal[, p_app := apparent_fit$predict(p_female)]
oof <- matrix(NA_real_, nrow = nrow(cal), ncol = REPEATS)
selection_log <- vector("list", REPEATS * K); log_i <- 0L
for (r in seq_len(REPEATS)) {
  set.seed(SEED + r - 1L)
  fold <- integer(nrow(cal))
  strata <- interaction(cal$y, cal$origin_group, drop = TRUE)
  for (s in levels(strata)) {
    idx <- which(strata == s)
    idx <- idx[sample.int(length(idx), replace = FALSE)]
    fold[idx] <- rep(seq_len(K), length.out = length(idx))
  }
  for (k in seq_len(K)) {
    test <- which(fold == k); train <- which(fold != k)
    fit <- fit_calibration(cal[train])
    oof[test, r] <- fit$predict(cal$p_female[test])
    log_i <- log_i + 1L
    selection_log[[log_i]] <- data.table(repeat_id = r, fold = k,
      selected = fit$selected, train_N = length(train), test_N = length(test),
      train_Brier_isotonic = fit$brier_isotonic, train_Brier_platt = fit$brier_platt)
  }
}
stopifnot(!anyNA(oof), all(oof >= 0 & oof <= 1))
selection_log <- rbindlist(selection_log)
fwrite(selection_log, file.path(OUTDIR, "calibration_cv_model_selection.csv"))

make_groups <- function(dd, label) {
  out <- list(list(grouping_type = "Overall", grouping = label, dt = dd))
  for (g in names(split(dd, dd$origin_group)))
    out[[length(out) + 1L]] <- list(grouping_type = "name_origin", grouping = g, dt = dd[origin_group == g])
  for (g in names(split(dd, dd$inst_region)))
    out[[length(out) + 1L]] <- list(grouping_type = "inst_region", grouping = g, dt = dd[inst_region == g])
  out
}
cv_metrics <- function(dt, idx) {
  g <- dt[idx]
  reps <- lapply(seq_len(REPEATS), function(r) {
    tmp <- copy(g); tmp[, p_cv := oof[idx, r]]
    weighted_metrics(tmp, "p_cv")
  })
  getv <- function(name) vapply(reps, function(z) z[[name]], numeric(1))
  list(cv_Brier_mean = mean(getv("Brier")), cv_Brier_sd = sd(getv("Brier")),
       cv_Brier_min = min(getv("Brier")), cv_Brier_max = max(getv("Brier")),
       cv_high_call_rate_mean = mean(getv("high_call_rate")),
       cv_high_accuracy_mean = mean(getv("high_accuracy"), na.rm = TRUE),
       cv_sensitivity_mean = mean(getv("sensitivity"), na.rm = TRUE),
       cv_specificity_mean = mean(getv("specificity"), na.rm = TRUE))
}
performance <- rbindlist(lapply(make_groups(cal, "Core calibration records"), function(spec) {
  idx <- which(cal$fa_id %in% spec$dt$fa_id)
  app <- weighted_metrics(spec$dt, "p_app")
  as.data.table(c(list(grouping_type = spec$grouping_type, grouping = spec$grouping,
                       N_female = sum(spec$dt$y == 1), N_male = sum(spec$dt$y == 0),
                       Brier_raw = weighted_brier(spec$dt$y, spec$dt$p_female, spec$dt$design_weight)),
                  app, cv_metrics(cal, idx)))
}))
fwrite(performance, file.path(OUTDIR, "calibration_stratified_performance.csv"))
raw_groups <- make_groups(core[gender_gold %in% c("F", "M")], "Core reference members")
raw_performance <- rbindlist(lapply(raw_groups, function(spec)
  as.data.table(c(list(grouping_type = spec$grouping_type, grouping = spec$grouping),
                  weighted_metrics(spec$dt, "p_female")))))
fwrite(raw_performance, file.path(OUTDIR, "core_raw_api_performance.csv"))

special_eval <- copy(special[gender_gold %in% c("F", "M")])
special_eval[, p_app := apparent_fit$predict(p_female)]
special_eval[is.na(p_female), p_app := apparent_fit$prevalence]
special_eval[, design_weight := 1]
special_groups <- make_groups(special_eval, "Special sample, descriptive only")
special_performance <- rbindlist(lapply(special_groups, function(spec) {
  raw <- weighted_metrics(spec$dt, "p_female")
  app <- weighted_metrics(spec$dt, "p_app")
  as.data.table(c(list(grouping_type = spec$grouping_type, grouping = spec$grouping,
                       N_female = sum(spec$dt$y == 1), N_male = sum(spec$dt$y == 0)),
                  list(raw_Brier = raw$Brier, calibrated_Brier = app$Brier,
                       raw_high_call_rate = raw$high_call_rate,
                       calibrated_high_call_rate = app$high_call_rate,
                       calibrated_high_accuracy = app$high_accuracy,
                       calibrated_sensitivity = app$sensitivity,
                       calibrated_specificity = app$specificity)))
}))
fwrite(special_performance, file.path(OUTDIR, "special_stratified_performance.csv"))
flow <- rbindlist(list(
  core[, .(sample = "core", N = .N, FM = sum(gender_gold %in% c("F", "M")), U = sum(final_status == "U"), NA_status = sum(final_status == "NA"), with_score_and_FM = sum(!is.na(p_female) & gender_gold %in% c("F", "M")))],
  special[, .(sample = "special", N = .N, FM = sum(gender_gold %in% c("F", "M")), U = sum(final_status == "U"), NA_status = sum(final_status == "NA"), with_score_and_FM = sum(!is.na(p_female) & gender_gold %in% c("F", "M")))],
  d0[core_selected == 1 & special_selected == 1, .(sample = "overlap", N = .N, FM = sum(gender_gold %in% c("F", "M")), U = sum(final_status == "U"), NA_status = sum(final_status == "NA"), with_score_and_FM = sum(!is.na(p_female) & gender_gold %in% c("F", "M")))]
))
fwrite(flow, file.path(OUTDIR, "dual_layer_sample_flow.csv"))
overall <- performance[grouping_type == "Overall"]
summary <- list(analyst = "ZZ H", analysis_date = "2026-11-12",
  design = list(core = "Weighted overall calibration and performance; exact design weights",
                special = "CN/KR/JP stratified description only; not pooled for unweighted overall calibration",
                cross_validation = "Repeated 10-fold CV, 10 repeats; isotonic versus Platt selected within training folds"),
  counts = list(full_N = 1783L, primary_N = 1536L, core_N = 800L, core_FM_N = 629L,
                core_score_and_FM_N = 575L, special_N = 239L, special_FM_N = 132L,
                overlap_N = 104L, primary_core_N = 689L, primary_core_score_and_FM_N = 495L),
  calibration = list(selected = apparent_fit$selected,
                     apparent_Brier_isotonic = apparent_fit$brier_isotonic,
                     apparent_Brier_platt = apparent_fit$brier_platt,
                     weighted_prevalence = apparent_fit$prevalence,
                     cv_Brier_mean = overall$cv_Brier_mean,
                     cv_Brier_sd = overall$cv_Brier_sd,
                     cv_Brier_min = overall$cv_Brier_min,
                     cv_Brier_max = overall$cv_Brier_max,
                     isotonic_selections = sum(selection_log$selected == "isotonic"),
                     platt_selections = sum(selection_log$selected == "platt")),
  caveats = c("Special-sample performance is descriptive and includes 104 people also present in the core; it is not independent validation.",
              "U labels remain unresolved and are not treated as F/M reference labels.",
              "Original in_rs is a legacy marker and is not used for new calibration."))
write_json(summary, file.path(OUTDIR, "validation_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
saveRDS(list(calibration_ids = cal$fa_id, oof_predictions = oof,
             performance = performance, raw_performance = raw_performance,
             special_performance = special_performance, summary = summary),
        file.path(OUTDIR, "validation_results.rds"))
cat(sprintf("Validation complete: %s; apparent Brier %.5f; repeated-CV Brier %.5f (range %.5f to %.5f)\n",
            apparent_fit$selected, overall$Brier, overall$cv_Brier_mean,
            overall$cv_Brier_min, overall$cv_Brier_max))
