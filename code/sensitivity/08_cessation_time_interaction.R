#!/usr/bin/env Rscript
# Analyst: ZZ H
# Sensitivity analysis adding the same female time-varying coefficient to
# cessation as used for graduation. The cohort, A0 events, calibration,
# covariates, and seeds are aligned with the primary analysis.

suppressMessages({library(rstpm2); library(data.table); library(yaml); library(jsonlite)})
options(warn = 1)

cfg <- yaml.load_file("config/config.yaml")
stopifnot(isTRUE(cfg$annotation_lock), file.exists("data/annotation_lock.timestamp"))
B <- as.integer(Sys.getenv("BOOT_B", "1000"))
M_OBS <- 20L
M_BOOT <- as.integer(Sys.getenv("BOOT_M", "1"))
DFB <- cfg$model$baseline_df
DFT <- cfg$model$gender_time_interaction_df
TMAX <- 15
GRID <- 50L
TS <- seq(0.05, TMAX, length.out = GRID)
HR_TIMES <- c(1, 5, 10, 15)
HR_NAMES <- paste0("cess_logHR_t", HR_TIMES)
RISK_NAMES <- c("riskF", "riskM", "RD", "RR")
METRICS <- c(RISK_NAMES, HR_NAMES)
PART_DIR <- "data/bootstrap_cessation_tvc_parts"
dir.create(PART_DIR, recursive = TRUE, showWarnings = FALSE)

d0 <- fread("data/analysis_cohort.csv")[primary_cohort == TRUE]
stopifnot(nrow(d0) == 1536L)
stopifnot(d0[, as.integer(table(factor(event_type, levels = c(0, 1, 2))))] == c(136L, 866L, 534L))
d0[, female_gold := fifelse(gender_gold == "F", 1, fifelse(gender_gold == "M", 0, NA_real_))]
d0[, ev1 := as.integer(event_type == 1)]
d0[, ev2 := as.integer(event_type == 2)]
d0[, `:=`(lp = log1p(preentry_pubs), lc = log1p(preentry_citations),
          lh = preentry_hindex, lco = log1p(preentry_coauthors), li = intl_collab_share)]
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
fix_levels(d0)
d0[, cluster := fifelse(is.na(inst_id) | inst_id == "", paste0("solo_", fa_id), inst_id)]
MODC <- c("female", "ey", "index_journal", "inst_tier", "ctry", "academic_age",
          "lp", "lc", "lh", "lco", "li", "no_preentry_work", "country_incomplete", "time_years")
RHS <- "female + ey + index_journal + inst_tier + ctry + academic_age + lp + lc + lh + lco + li + no_preentry_work"
F1 <- as.formula(paste("Surv(time_years, ev1) ~", RHS))
F2 <- as.formula(paste("Surv(time_years, ev2) ~", RHS))

source("code/sensitivity/weighted_calibration.R")
rs <- d0[core_selected == 1 & !is.na(p_female) & female_gold %in% c(0, 1)]
stopifnot(nrow(rs) == 495L, !anyNA(rs$design_weight), all(rs$design_weight > 0))
cal <- fit_calibration(rs)
p_calibrate <- cal$predict
prev_rs <- cal$prevalence
d0[, p_cal := fifelse(!is.na(p_female), p_calibrate(p_female), prev_rs)]
cat(sprintf("Calibration: %s, N=%d, prevalence=%.12f\n", cal$selected, nrow(rs), prev_rs))

impute_one <- function(dd, seed) {
  set.seed(seed)
  dd[, female := fifelse(!is.na(female_gold), female_gold, rbinom(.N, 1, p_cal))]
  dd
}

fit_pair <- function(dd, cessation_tvc = TRUE, init1 = NULL, init2 = NULL) {
  f1 <- stpm2(F1, data = dd, df = DFB, tvc = list(female = DFT), init = init1)
  if (cessation_tvc) {
    f2 <- stpm2(F2, data = dd, df = DFB, tvc = list(female = DFT), init = init2)
  } else {
    f2 <- stpm2(F2, data = dd, df = DFB, init = init2)
  }
  list(f1 = f1, f2 = f2, init1 = coef(f1), init2 = coef(f2))
}

haz_grids <- function(f1, f2, dd) {
  n <- nrow(dd); G <- length(TS)
  nd0 <- dd[, ..MODC]; nd0[, female := 0]
  long <- nd0[rep(1:n, times = G)]; long[, time_years := rep(TS, each = n)]
  H10 <- matrix(predict(f1, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  H20 <- matrix(predict(f2, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  one <- nd0[rep(1, 2 * G)]; one[, time_years := rep(TS, times = 2)]
  one[, female := rep(c(0, 1), each = G)]
  lp1 <- predict(f1, newdata = one, type = "link")
  lp2 <- predict(f2, newdata = one, type = "link")
  d1 <- lp1[(G + 1):(2 * G)] - lp1[1:G]
  d2 <- lp2[(G + 1):(2 * G)] - lp2[1:G]
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

risk_estimates <- function(f1, f2, dd) {
  hg <- haz_grids(f1, f2, dd)
  rF <- cif_from_H(hg$H1F, hg$H2F)
  rM <- cif_from_H(hg$H1M, hg$H2M)
  c(riskF = unname(rF), riskM = unname(rM), RD = unname(rF - rM), RR = unname(rF / rM))
}

cessation_log_hr <- function(fit, dd) {
  nd <- dd[rep(1, 2 * length(HR_TIMES)), ..MODC]
  nd[, time_years := rep(HR_TIMES, times = 2)]
  nd[, female := rep(c(0, 1), each = length(HR_TIMES))]
  lp <- predict(fit, newdata = nd, type = "link")
  values <- lp[(length(HR_TIMES) + 1):(2 * length(HR_TIMES))] - lp[seq_along(HR_TIMES)]
  setNames(values, HR_NAMES)
}

metrics_tvc <- function(fp, dd) {
  c(risk_estimates(fp$f1, fp$f2, dd), cessation_log_hr(fp$f2, dd))
}

pool_metrics <- function(em) {
  em <- as.matrix(em)
  stopifnot(all(METRICS %in% colnames(em)))
  c(riskF = mean(em[, "riskF"]), riskM = mean(em[, "riskM"]),
    RD = mean(em[, "RD"]), RR = exp(mean(log(em[, "RR"]))),
    setNames(colMeans(em[, HR_NAMES, drop = FALSE]), HR_NAMES))
}

warm <- fit_pair(impute_one(copy(d0), 902607), cessation_tvc = TRUE)
INIT1 <- warm$init1
INIT2 <- warm$init2

obs_cells <- list()
obs_errors <- character()
for (m in seq_len(M_OBS)) {
  result <- tryCatch({
    dd <- impute_one(copy(d0), 902607 + m)
    tvc <- fit_pair(dd, cessation_tvc = TRUE)
    ph <- fit_pair(dd, cessation_tvc = FALSE)
    tvc_metrics <- metrics_tvc(tvc, dd)
    ph_risk <- risk_estimates(ph$f1, ph$f2, dd)
    names(ph_risk) <- paste0("PH_", names(ph_risk))
    ll_tvc <- stats4::logLik(tvc$f2)
    ll_ph <- stats4::logLik(ph$f2)
    lr_df <- attr(ll_tvc, "df") - attr(ll_ph, "df")
    lrt <- 2 * (as.numeric(ll_tvc) - as.numeric(ll_ph))
    diagnostics <- c(logLik_tvc = as.numeric(ll_tvc), logLik_ph = as.numeric(ll_ph),
                     AIC_tvc = stats4::AIC(tvc$f2), AIC_ph = stats4::AIC(ph$f2),
                     delta_AIC = stats4::AIC(tvc$f2) - stats4::AIC(ph$f2),
                     LRT = lrt, LRT_df = lr_df,
                     LRT_p = pchisq(lrt, df = lr_df, lower.tail = FALSE),
                     cessation_logHR_range = max(tvc_metrics[HR_NAMES]) - min(tvc_metrics[HR_NAMES]))
    c(list(m = m), as.list(tvc_metrics), as.list(ph_risk), as.list(diagnostics))
  }, error = function(e) e)
  if (inherits(result, "error")) {
    obs_errors <- c(obs_errors, paste0("m=", m, ": ", conditionMessage(result)))
  } else {
    obs_cells[[length(obs_cells) + 1L]] <- as.data.table(result)
  }
}
if (length(obs_errors)) stop(paste(obs_errors, collapse = " | "))
obs_dt <- rbindlist(obs_cells)
stopifnot(nrow(obs_dt) == M_OBS, all(obs_dt$LRT_df == DFT), all(obs_dt$LRT >= -1e-7))
fwrite(obs_dt, "data/observed_cessation_tvc_cells.csv")
obs_est <- pool_metrics(obs_dt[, ..METRICS])
obs_ph <- colMeans(as.matrix(obs_dt[, .(PH_riskF, PH_riskM, PH_RD, PH_RR)]))
names(obs_ph) <- sub("^PH_", "", names(obs_ph))

clusters <- split(seq_len(nrow(d0)), d0$cluster)
CL <- names(clusters)
first_error <- function(x) {
  x <- x[nzchar(x)]
  if (!length(x)) "" else paste(unique(x), collapse = " | ")
}
boot_rep <- function(b) {
  set.seed(902607 + b * 1000)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- copy(d0[idx]); apply_scale(db); fix_levels(db)
  ests <- list(); errors <- character()
  for (m in seq_len(M_BOOT)) {
    result <- tryCatch({
      dd <- impute_one(copy(db), 902607 + b * 1000 + m)
      fp <- fit_pair(dd, cessation_tvc = TRUE, init1 = INIT1, init2 = INIT2)
      metrics_tvc(fp, dd)
    }, error = function(e) e)
    if (inherits(result, "error")) errors <- c(errors, conditionMessage(result))
    else ests[[length(ests) + 1L]] <- result
  }
  if (!length(ests)) return(data.table(b = b, nfail = M_BOOT, error = first_error(errors)))
  pooled <- pool_metrics(do.call(rbind, ests))
  as.data.table(c(list(b = b, nfail = M_BOOT - length(ests)), as.list(pooled),
                  list(error = first_error(errors))))
}

part_reader <- function() {
  files <- list.files(PART_DIR, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  parts <- rbindlist(lapply(files, function(f) tryCatch(fread(f), error = function(e) NULL)), fill = TRUE)
  stopifnot(!anyDuplicated(parts$b))
  parts
}
old <- part_reader()
done_b <- if (is.null(old)) integer(0) else old$b
todo <- setdiff(seq_len(B), done_b)
cat(sprintf("Cessation-TVC bootstrap: %d/%d complete; M_BOOT=%d\n", length(done_b), B, M_BOOT))
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  started <- Sys.time()
  res <- rbindlist(lapply(part, boot_rep), fill = TRUE)
  fwrite(res, file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1])))
  cat(sprintf("Replicates %d-%d written in %.1f seconds\n", part[1], part[length(part)],
              as.numeric(difftime(Sys.time(), started, units = "secs"))))
  flush.console()
}
br_all <- part_reader()
stopifnot(!is.null(br_all), nrow(br_all) == B, setequal(br_all$b, seq_len(B)))
br_all <- br_all[order(b)]
fwrite(br_all, "data/bootstrap_cessation_tvc.csv")
valid <- br_all[is.finite(RD)]
ci <- function(x) unname(quantile(x, c(0.025, 0.975), na.rm = TRUE))
rd_ci <- ci(valid$RD)
rr_ci <- exp(ci(log(valid$RR)))
hr_rows <- lapply(seq_along(HR_TIMES), function(j) {
  metric <- HR_NAMES[j]
  data.table(cause = "cessation", time_years = HR_TIMES[j],
             logHR = unname(obs_est[metric]), HR = exp(unname(obs_est[metric])),
             HR_lo95 = exp(ci(valid[[metric]])[1]), HR_hi95 = exp(ci(valid[[metric]])[2]))
})
hr <- rbindlist(hr_rows)
fwrite(hr, "data/results_cessation_tvc_HR.csv")

observed_comparison <- list(
  logLik_tvc_mean = mean(obs_dt$logLik_tvc), logLik_ph_mean = mean(obs_dt$logLik_ph),
  delta_AIC_mean = mean(obs_dt$delta_AIC), delta_AIC_median = median(obs_dt$delta_AIC),
  delta_AIC_range = range(obs_dt$delta_AIC),
  AIC_tvc_lower_in_imputations = sum(obs_dt$delta_AIC < 0),
  LRT_df = unique(obs_dt$LRT_df),
  LRT_p_median = median(obs_dt$LRT_p), LRT_p_range = range(obs_dt$LRT_p),
  LRT_p_below_0_05_in_imputations = sum(obs_dt$LRT_p < 0.05),
  cessation_logHR_range_mean = mean(obs_dt$cessation_logHR_range),
  cessation_logHR_range_max = max(obs_dt$cessation_logHR_range),
  inference_note = "Imputation-specific likelihood-ratio statistics are reported separately; no pooled P value is claimed."
)

out <- list(
  analyst = "ZZ H", date = "2026-11-16",
  analysis = "Post hoc sensitivity: cessation female TVC natural spline on log time, df=3",
  alignment_with_primary = "The graduation model, cohort, A0 events, covariates, calibration, and seeds are unchanged",
  N = nrow(d0), events = list(graduation = 866L, cessation = 534L, administrative_censoring = 136L),
  calibration = cal$selected, calibration_N = nrow(rs), M_obs = M_OBS, M_boot = M_BOOT,
  baseline_df = DFB, cessation_gender_time_df = DFT,
  observed_tvc = as.list(obs_est), observed_ph_reference = as.list(obs_ph),
  observed_cessation_HR = hr,
  RD_ci = rd_ci, RR_ci = rr_ci,
  riskF_ci = ci(valid$riskF), riskM_ci = ci(valid$riskM),
  bootstrap_requested = B, bootstrap_valid = nrow(valid), failed_cells = sum(br_all$nfail),
  observed_model_comparison = observed_comparison,
  interpretation = "Cessation HRs are time-specific under this sensitivity model, not risk ratios or mechanistic effects. Compare the standardized RD with the locked primary model; do not infer proportional hazards from a non-significant test alone.",
  session = list(R = R.version.string, rstpm2 = as.character(packageVersion("rstpm2")),
                 data.table = as.character(packageVersion("data.table")))
)
saveRDS(out, "data/results_cessation_tvc_sensitivity.rds")
write_json(out, "data/results_cessation_tvc_sensitivity.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("Observed TVC RD %.12f | PH reference RD %.12f | delta %.12f\n",
            obs_est["RD"], obs_ph["RD"], obs_est["RD"] - obs_ph["RD"]))
print(hr)
print(observed_comparison)
cat(sprintf("Bootstrap valid %d/%d; failed cells %d\n", nrow(valid), B, sum(br_all$nfail)))
