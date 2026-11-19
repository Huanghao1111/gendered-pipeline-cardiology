#!/usr/bin/env Rscript
# Repeated internal cross-validation for the 495-person primary calibration set.
# Analyst: ZZ H

suppressMessages({library(data.table); library(jsonlite)})
source("code/primary/weighted_calibration.R")

SEED <- 20261112L
K <- 10L
REPEATS <- 10L
OUTDIR <- file.path("data", "gender_validation_primary495")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

d0 <- fread("data/analysis_cohort.csv")
primary <- tolower(trimws(as.character(d0$primary_cohort))) == "true"
stopifnot(nrow(d0) == 1783L, sum(primary) == 1536L,
          sum(d0$event_type[primary] == 1) == 866L,
          sum(d0$event_type[primary] == 2) == 534L,
          sum(d0$event_type[primary] == 0) == 136L)
d0[, y := fifelse(gender_gold == "F", 1, fifelse(gender_gold == "M", 0, NA_real_))]
d0[, female_gold := y]
cal <- d0[primary & core_selected == 1 & !is.na(p_female) & !is.na(y) & !is.na(design_weight)]
stopifnot(nrow(cal) == 495L, all(cal$design_weight > 0), !anyDuplicated(cal$fa_id))
cal[, origin_group := fifelse(is.na(final_name_origin), "U",
                       fifelse(final_name_origin %in% c("CN", "KR", "JP"), final_name_origin,
                       fifelse(final_name_origin == "other", "other", "U")))]

apparent <- fit_calibration(cal)
stopifnot(apparent$selected == "isotonic")
oof <- matrix(NA_real_, nrow = nrow(cal), ncol = REPEATS)
selection_log <- vector("list", REPEATS * K)
log_i <- 0L
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
    test <- which(fold == k)
    train <- which(fold != k)
    fit <- fit_calibration(cal[train])
    oof[test, r] <- fit$predict(cal$p_female[test])
    log_i <- log_i + 1L
    selection_log[[log_i]] <- data.table(
      repeat_id = r, fold = k, selected = fit$selected,
      train_N = length(train), test_N = length(test),
      train_Brier_isotonic = fit$brier_isotonic,
      train_Brier_platt = fit$brier_platt)
  }
}
stopifnot(!anyNA(oof), all(oof >= 0 & oof <= 1))
selection_log <- rbindlist(selection_log)
stopifnot(nrow(selection_log) == REPEATS * K)

repeat_metrics <- rbindlist(lapply(seq_len(REPEATS), function(r) {
  data.table(repeat_id = r,
             Brier = weighted_brier(cal$y, oof[, r], cal$design_weight))
}))
stopifnot(nrow(repeat_metrics) == REPEATS)
cv_brier <- repeat_metrics$Brier
summary <- list(
  analyst = "ZZ H",
  analysis_date = "2026-11-16",
  calibration_population = "Primary 2005-2009 core members with a name score and definite F/M label",
  n = nrow(cal),
  apparent_brier_isotonic = apparent$brier_isotonic,
  apparent_brier_platt = apparent$brier_platt,
  selected = apparent$selected,
  weighted_prevalence = apparent$prevalence,
  cross_validation = "Repeated stratified 10-fold internal cross-validation with 10 repeats; selection repeated within each training fold",
  seed_start = SEED,
  folds = K,
  repeats = REPEATS,
  isotonic_selections = sum(selection_log$selected == "isotonic"),
  platt_selections = sum(selection_log$selected == "platt"),
  cv_Brier_mean = mean(cv_brier),
  cv_Brier_sd = sd(cv_brier),
  cv_Brier_min = min(cv_brier),
  cv_Brier_max = max(cv_brier)
)
write_json(summary, file.path(OUTDIR, "primary495_cv_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(selection_log, file.path(OUTDIR, "primary495_cv_model_selection.csv"))
fwrite(repeat_metrics, file.path(OUTDIR, "primary495_cv_repeat_metrics.csv"))
oof_long <- rbindlist(lapply(seq_len(REPEATS), function(r) {
  data.table(fa_id = cal$fa_id, repeat_id = r, y = cal$y,
             design_weight = cal$design_weight, p_cv = oof[, r])
}))
fwrite(oof_long, file.path(OUTDIR, "primary495_cv_oof_predictions.csv"))
writeLines(c(
  "PASS primary calibration n = 495",
  "PASS all 100 training-fold selections recorded",
  "PASS all 4950 out-of-fold predictions finite and within [0,1]",
  "PASS apparent calibration fitted on the same 495-person set",
  "PASS all 10 repeat-specific held-out weighted Brier scores calculated"
), file.path(OUTDIR, "primary495_cv_QA.txt"))
print(toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = NA))
