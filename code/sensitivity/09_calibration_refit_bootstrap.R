#!/usr/bin/env Rscript
# Propagate calibration-function estimation uncertainty.
# Each replicate combines a stratified calibration-sample bootstrap with the
# institution-cluster bootstrap and calibrated Bernoulli gender classification.

suppressMessages({
  library(rstpm2)
  library(data.table)
  library(yaml)
  library(parallel)
  library(jsonlite)
})

cfg <- yaml.load_file("config/config.yaml")
stopifnot(isTRUE(cfg$annotation_lock), file.exists("data/annotation_lock.timestamp"))
B <- as.integer(Sys.getenv("BOOT_B", "1000"))
M_BOOT <- 1L
DFB <- cfg$model$baseline_df
DFT <- cfg$model$gender_time_interaction_df
FAIL_T <- cfg$model$convergence_fail_threshold
TMAX <- 15
GRID <- 50
TS <- seq(0.05, TMAX, length.out = GRID)
CORES <- 1L
PART_DIR <- "results/calibration_refit_parts"
dir.create(PART_DIR, recursive = TRUE, showWarnings = FALSE)

d0 <- fread("data/analysis_cohort.csv")[primary_cohort == TRUE]
d0[, female_gold := fifelse(gender_gold == "F", 1,
                            fifelse(gender_gold == "M", 0, NA_real_))]
d0[, ev1 := as.integer(event_type == 1)]
d0[, ev2 := as.integer(event_type == 2)]
stopifnot(nrow(d0) == 1536L, sum(d0$ev1) == 866L, sum(d0$ev2) == 534L,
          sum(d0$ev1 == 0L & d0$ev2 == 0L) == 136L)
d0[, `:=`(lp = log1p(preentry_pubs), lc = log1p(preentry_citations),
          lh = preentry_hindex, lco = log1p(preentry_coauthors),
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
invisible(apply_scale(d0))
TOPC <- setdiff(names(sort(-table(d0$inst_country))), "未知")[1:8]
CLEV <- c(TOPC, "Other")
TLEV <- c("高", "中", "低", "缺失")
EYLEV <- sort(unique(d0$entry_year))
JLEV <- sort(unique(d0$index_journal))
fix_levels <- function(dd) {
  dd[, ctry := factor(ifelse(inst_country %in% TOPC, inst_country, "Other"), levels = CLEV)]
  dd[, inst_tier := factor(inst_tier, levels = TLEV)]
  dd[, ey := factor(entry_year, levels = EYLEV)]
  dd[, index_journal := factor(index_journal, levels = JLEV)]
  dd
}
invisible(fix_levels(d0))
d0[, cluster := fifelse(is.na(inst_id) | inst_id == "", paste0("solo_", fa_id), inst_id)]
MODC <- c("female", "ey", "index_journal", "inst_tier", "ctry", "academic_age",
          "lp", "lc", "lh", "lco", "li", "no_preentry_work", "country_incomplete", "time_years")
RHS <- "female + ey + index_journal + inst_tier + ctry + academic_age + lp + lc + lh + lco + li + no_preentry_work"
F1 <- as.formula(paste("Surv(time_years, ev1) ~", RHS))
F2 <- as.formula(paste("Surv(time_years, ev2) ~", RHS))

region <- function(x) {
  europe <- c("AT", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES",
              "FI", "FR", "GB", "GR", "HR", "HU", "IE", "IS", "IT", "LT",
              "LU", "LV", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK")
  east_asia <- c("CN", "KR", "JP", "TW", "HK")
  x <- toupper(trimws(x))
  ifelse(is.na(x) | x == "" | x == "未知", "Unknown",
         ifelse(x %in% c("US", "CA"), "North America",
                ifelse(x %in% europe, "Europe",
                       ifelse(x %in% east_asia, "East Asia", "Other"))))
}
d0[, inst_region := region(inst_country)]
source("code/primary/weighted_calibration.R")
rs <- d0[core_selected == 1 & !is.na(p_female) & female_gold %in% c(0, 1)]
stopifnot(nrow(rs) == 495L, !anyNA(rs$design_weight), all(rs$design_weight > 0))
stopifnot(length(unique(rs$inst_region)) == 5L)
cal0 <- fit_calibration(rs)
stopifnot(cal0$selected == "isotonic")
d0[, p_cal := fifelse(!is.na(p_female), cal0$predict(p_female), cal0$prevalence)]
cal_strata <- split(seq_len(nrow(rs)), rs$inst_region)
stopifnot(identical(sort(names(cal_strata)),
                    sort(c("North America", "Europe", "East Asia", "Other", "Unknown"))))

haz_grids <- function(f1, f2, dd) {
  n <- nrow(dd)
  G <- length(TS)
  nd0 <- dd[, ..MODC]
  nd0[, female := 0]
  long <- nd0[rep(1:n, times = G)]
  long[, time_years := rep(TS, each = n)]
  H10 <- matrix(predict(f1, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  H20 <- matrix(predict(f2, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  one <- nd0[rep(1, 2 * G)]
  one[, time_years := rep(TS, times = 2)]
  one[, female := rep(c(0, 1), each = G)]
  d1 <- predict(f1, newdata = one, type = "link")[(G + 1):(2 * G)] -
        predict(f1, newdata = one, type = "link")[1:G]
  d2 <- predict(f2, newdata = one, type = "link")[(G + 1):(2 * G)] -
        predict(f2, newdata = one, type = "link")[1:G]
  list(H1F = sweep(H10, 2, exp(d1), `*`), H1M = H10,
       H2F = sweep(H20, 2, exp(d2), `*`), H2M = H20)
}
cif_from_H <- function(H1, H2) {
  G <- ncol(H1)
  S <- exp(-(H1 + H2))
  Smid <- (cbind(1, S[, -G]) + S) / 2
  h1 <- t(apply(H1, 1, function(x) diff(c(0, x))))
  mean(rowSums(h1 * Smid))
}
fit_gcomp <- function(dd, init1 = NULL, init2 = NULL) {
  f1 <- stpm2(F1, data = dd, df = DFB, tvc = list(female = DFT), init = init1)
  f2 <- stpm2(F2, data = dd, df = DFB, init = init2)
  hg <- haz_grids(f1, f2, dd)
  rF <- cif_from_H(hg$H1F, hg$H2F)
  rM <- cif_from_H(hg$H1M, hg$H2M)
  list(c(riskF = rF, riskM = rM, RD = rF - rM, RR = rF / rM),
       init1 = coef(f1), init2 = coef(f2))
}
impute_one <- function(dd, seed) {
  set.seed(seed)
  dd[, female := fifelse(!is.na(female_gold), female_gold, rbinom(.N, 1, p_cal))]
  dd
}

obs0 <- fit_gcomp(impute_one(copy(d0), seed = 902607))
INIT1 <- obs0$init1
INIT2 <- obs0$init2
clusters <- split(seq_len(nrow(d0)), d0$cluster)
CL <- names(clusters)

boot_rep <- function(b) {
  empty <- c(b = b, nfail = 1L, riskF = NA_real_, riskM = NA_real_, RD = NA_real_,
             RR = NA_real_, calibration_selected = NA_character_, calibration_prevalence = NA_real_)
  tryCatch({
    set.seed(902607 + b * 1000)
    idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
    db <- copy(d0[idx])
    apply_scale(db)
    fix_levels(db)

    set.seed(902607 + b * 1000 + 777)
    cal_idx <- unlist(lapply(cal_strata, function(ii) sample(ii, length(ii), replace = TRUE)),
                      use.names = FALSE)
    cal_b <- fit_calibration(rs[cal_idx])
    p_new <- rep(cal_b$prevalence, nrow(db))
    has_score <- !is.na(db$p_female)
    p_new[has_score] <- cal_b$predict(db$p_female[has_score])
    stopifnot(all(is.finite(p_new)), all(p_new >= 0), all(p_new <= 1))
    db[, p_cal := p_new]

    est <- fit_gcomp(impute_one(db, seed = 902607 + b * 1000 + 1L),
                     init1 = INIT1, init2 = INIT2)[[1]]
    c(b = b, nfail = 0L, as.list(est), calibration_selected = cal_b$selected,
      calibration_prevalence = cal_b$prevalence)
  }, error = function(e) empty)
}

read_parts <- function() {
  pf <- list.files(PART_DIR, pattern = "\\.csv$", full.names = TRUE)
  if (!length(pf)) return(NULL)
  parts <- rbindlist(lapply(pf, function(f) tryCatch(fread(f), error = function(e) NULL)), fill = TRUE)
  stopifnot(!anyDuplicated(parts$b))
  parts
}
old_parts <- read_parts()
done_b <- if (is.null(old_parts)) integer(0) else old_parts$b
todo <- setdiff(seq_len(B), done_b)
cat(sprintf("calibration-refit bootstrap: %d/%d done, %d to go\n", length(done_b), B, length(todo)))
flush.console()
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  t0 <- Sys.time()
  res <- as.data.table(do.call(rbind, mclapply(part, boot_rep, mc.cores = CORES)))
  fwrite(res, file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1])))
  cat(sprintf("  reps %d-%d done in %.0fs\n", part[1], part[length(part)],
              difftime(Sys.time(), t0, units = "secs")))
  flush.console()
}

br_all <- read_parts()
stopifnot(!is.null(br_all), nrow(br_all) == B)
br <- br_all[!is.na(RD)]
stopifnot(nrow(br) > 0)
ci <- function(x) unname(quantile(x, c(.025, .975), na.rm = TRUE))
rdci <- ci(br$RD)
rrci <- exp(ci(log(br$RR)))
fail_rate <- sum(br_all$nfail) / nrow(br_all)
summary <- list(
  analysis = "Stratified calibration-refit and institution-cluster bootstrap sensitivity",
  primary_cohort_N = nrow(d0),
  calibration_N = nrow(rs),
  calibration_strata = sort(names(cal_strata)),
  bootstrap_requested = B,
  bootstrap_valid = nrow(br),
  bootstrap_failed = nrow(br_all) - nrow(br),
  failure_rate = fail_rate,
  failure_guard_triggered = fail_rate > FAIL_T,
  original_RD_pp = -7.35,
  original_RD_CI95_pp = c(-12.44, -1.79),
  refit_bootstrap_RD_mean_pp = 100 * mean(br$RD),
  refit_bootstrap_RD_median_pp = 100 * median(br$RD),
  refit_bootstrap_RD_CI95_pp = 100 * rdci,
  refit_bootstrap_riskF_CI95_percent = 100 * ci(br$riskF),
  refit_bootstrap_riskM_CI95_percent = 100 * ci(br$riskM),
  refit_bootstrap_RR_CI95 = rrci,
  calibration_selection_counts = as.list(table(br_all$calibration_selected, useNA = "ifany")),
  calibration_prevalence_range = range(br_all$calibration_prevalence, na.rm = TRUE),
  calibration_resampling = "Independent resampling with replacement within five institutional-region strata; original design weights retained",
  outer_resampling = "Institution-cluster bootstrap using the original seeds and sample size",
  imputation_per_replicate = M_BOOT,
  outcome_model_changed = FALSE,
  accepted_primary_analysis_replaced = FALSE
)
fwrite(br_all, "results/calibration_refit_bootstrap.csv")
write_json(summary, "results/calibration_refit_summary.json", auto_unbox = TRUE,
           pretty = TRUE, digits = 12, na = "string")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
print(summary)
