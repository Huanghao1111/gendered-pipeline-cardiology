#!/usr/bin/env Rscript
source("code/secondary/common.R")
B <- as.integer(Sys.getenv("BOOT_B3", "500"))
PART_DIR <- "work_tmp/reentry_parts"
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)
EVER_METRICS <- c("gradF", "gradM", "RDgrad")
RE_METRICS <- c("reentryF", "reentryM", "RDreentry",
                "postReentryGradF", "postReentryGradM", "RDpostReentryGrad",
                "gradAfterCessF", "gradAfterCessM", "RDgradAfterCess")
METRIC_MAP <- list(A3 = EVER_METRICS, A0 = RE_METRICS)

set_ever <- function(dd) {
  z <- copy(dd)
  z[, time_years := ever_time_A3]
  z[, ev1 := as.integer(ever_event_A3)]
  z
}
fit_ever <- function(dd, init = NULL) {
  z <- set_ever(dd)
  f <- stpm2(as.formula(paste("Surv(time_years, ev1) ~", RHS)), data = z,
             df = DFB, tvc = list(female = DFT), init = init)
  list(fit = f, init = coef(f))
}
estimate_ever <- function(dd, init = NULL) {
  z <- set_ever(dd)
  f <- stpm2(as.formula(paste("Surv(time_years, ev1) ~", RHS)), data = z,
             df = DFB, tvc = list(female = DFT), init = init)
  nd <- z[, ..MODC]
  nd[, time_years := TMAX]
  nd0 <- copy(nd)[, female := 0]
  nd1 <- copy(nd)[, female := 1]
  H0 <- predict(f, newdata = nd0, type = "cumhaz")
  H1 <- predict(f, newdata = nd1, type = "cumhaz")
  gradM <- mean(1 - exp(-H0))
  gradF <- mean(1 - exp(-H1))
  c(gradF = gradF, gradM = gradM, RDgrad = gradF - gradM)
}
safe_mean <- function(x) if (length(x)) mean(x) else NA_real_
reentry_values <- function(dd, variant) {
  z <- dd[event_type == 2]
  reentry <- z[[paste0("reentry_", variant)]]
  regrad <- z[[paste0("reentry_graduation_", variant)]]
  female <- z$female
  out <- numeric()
  for (sex in c(1, 0)) {
    sel <- female == sex
    label <- if (sex == 1) "F" else "M"
    out[paste0("reentry", label)] <- safe_mean(reentry[sel])
    out[paste0("postReentryGrad", label)] <- safe_mean(regrad[sel & reentry == 1])
    out[paste0("gradAfterCess", label)] <- safe_mean(regrad[sel])
  }
  out["RDreentry"] <- out["reentryF"] - out["reentryM"]
  out["RDpostReentryGrad"] <- out["postReentryGradF"] - out["postReentryGradM"]
  out["RDgradAfterCess"] <- out["gradAfterCessF"] - out["gradAfterCessM"]
  out[RE_METRICS]
}

warm <- fit_ever(impute_one(copy(d0), 932607))
INIT <- warm$init
clusters <- split(seq_len(nrow(d0)), d0$cluster); CL <- names(clusters)

## Observed-data estimates.
ever_cells <- list(); reentry_cells <- list(); errors <- list()
for (m in seq_len(M_OBS)) {
  dd <- impute_one(copy(d0), 932607 + m)
  fit <- tryCatch(estimate_ever(dd), error = function(e) e)
  if (inherits(fit, "error")) {
    errors[[length(errors) + 1L]] <- data.table(m = m, scenario = "A3", message = conditionMessage(fit))
  } else {
    ever_cells[[length(ever_cells) + 1L]] <- as.data.table(c(list(m = m), as.list(fit)))
  }
  values <- reentry_values(dd, "A0")
  reentry_cells[[length(reentry_cells) + 1L]] <- as.data.table(c(list(m = m, scenario = "A0"), as.list(values)))
}
if (length(errors)) stop(paste(vapply(errors, function(x) x$message, character(1)), collapse = " | "))
ever_dt <- rbindlist(ever_cells)
reentry_dt <- rbindlist(reentry_cells)
ever_summary <- ever_dt[, lapply(.SD, mean), .SDcols = EVER_METRICS]
reentry_summary <- reentry_dt[, lapply(.SD, mean), by = scenario, .SDcols = RE_METRICS]
fwrite(ever_dt, "data/observed_ever_cells_A3.csv")
fwrite(ever_summary, "data/observed_ever_summary_A3.csv")
fwrite(reentry_dt, "data/observed_reentry_cells.csv")
fwrite(reentry_summary, "data/observed_reentry_summary.csv")

## Institution-cluster bootstrap; A3 model failures are explicit.
boot_rep <- function(b) {
  set.seed(932607 + b * 1000)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- copy(d0[idx]); apply_scale(db); fix_levels(db)
  ever <- list(); reentry <- list(A0 = list()); errors <- character()
  for (m in seq_len(M_BOOT)) {
    di <- impute_one(copy(db), 932607 + b * 1000 + m)
    fit <- tryCatch(estimate_ever(di, INIT), error = function(e) e)
    if (inherits(fit, "error")) errors <- c(errors, conditionMessage(fit)) else ever[[length(ever) + 1L]] <- fit
    reentry$A0[[length(reentry$A0) + 1L]] <- reentry_values(di, "A0")
  }
  answer <- list(b = b, nfail = M_BOOT - length(ever))
  for (metric in EVER_METRICS) {
    answer[[paste("A3", metric, sep = "_")]] <- if (!length(ever)) NA_real_ else mean(vapply(ever, function(x) unname(x[metric]), numeric(1)))
  }
  for (metric in RE_METRICS) {
    answer[[paste("A0", metric, sep = "_")]] <- mean(vapply(reentry$A0, function(x) unname(x[metric]), numeric(1)), na.rm = TRUE)
  }
  answer$error <- first_error(errors)
  as.data.table(answer)
}

old <- part_reader(PART_DIR)
done_b <- if (is.null(old)) integer(0) else old$b
todo <- setdiff(1:B, done_b)
cat(sprintf("A3/re-entry bootstrap: %d/%d done, M_BOOT=%d\n", length(done_b), B, M_BOOT))
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  t0 <- Sys.time()
  res <- rbindlist(lapply(part, boot_rep), fill = TRUE)
  fwrite(res, file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1])))
  cat(sprintf("A3/re-entry reps %d-%d written in %.1fs; failed cells=%d\n",
              part[1], part[length(part)], as.numeric(difftime(Sys.time(), t0, units = "secs")),
              sum(res$nfail)))
  flush.console()
}
pooled <- part_reader(PART_DIR)
stopifnot(!is.null(pooled), nrow(pooled) == B, setequal(pooled$b, 1:B))
fwrite(pooled, "data/bootstrap_pooled_A3_reentry.csv")
summaries <- list()
for (v in names(METRIC_MAP)) for (metric in METRIC_MAP[[v]]) {
  cc <- paste(v, metric, sep = "_")
  x <- pooled[[cc]]
  summaries[[length(summaries) + 1L]] <- data.table(
    scenario = v, metric = metric, bootstrap_reps = sum(is.finite(x)),
    lower95 = ci(x)[1], upper95 = ci(x)[2])
}
intervals <- rbindlist(summaries)
fwrite(intervals, "data/bootstrap_intervals_A3_reentry.csv")
out <- list(analyst = "ZZ H", date = "2026-11-16", N = nrow(d0),
            M_obs = M_OBS, M_boot = M_BOOT, bootstrap_requested = B,
            bootstrap_valid = sum(is.finite(pooled$A3_RDgrad)),
            failed_cells = sum(pooled$nfail), failure_rate = sum(pooled$nfail) / (B * M_BOOT),
            calibration = cal$selected, brier_isotonic = brier_iso, brier_platt = brier_platt,
            observed_ever = ever_summary, observed_reentry = reentry_summary,
            bootstrap_intervals = intervals)
saveRDS(out, "data/results_A3_reentry.rds")
write_json(out, "data/results_A3_reentry.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(ever_summary)
print(reentry_summary)
print(intervals)
