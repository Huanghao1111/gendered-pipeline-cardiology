#!/usr/bin/env Rscript
source("code/secondary/common.R")
B <- as.integer(Sys.getenv("BOOT_B", "500"))
PART_DIR <- "work_tmp/cause_hr_parts"
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)
HR_TIMES <- c(1, 5, 10, 15)
METRICS <- c(paste0("grad_logHR_t", HR_TIMES), "cess_logHR_constant")

set_a0 <- function(dd) {
  z <- copy(dd)
  z[, ev1 := as.integer(event_type == 1)]
  z[, ev2 := as.integer(event_type == 2)]
  z
}
log_hr_at_times <- function(fit, z) {
  nd <- z[rep(1, 2 * length(HR_TIMES)), ..MODC]
  nd[, time_years := rep(HR_TIMES, times = 2)]
  nd[, female := rep(c(0, 1), each = length(HR_TIMES))]
  eta <- predict(fit, newdata = nd, type = "link")
  eta[(length(HR_TIMES) + 1):(2 * length(HR_TIMES))] - eta[seq_along(HR_TIMES)]
}
estimate_hr <- function(dd, init1 = NULL, init2 = NULL) {
  z <- set_a0(dd)
  fp <- fit_pair(z, init1, init2)
  graduation <- log_hr_at_times(fp$f1, z)
  cessation_all_times <- log_hr_at_times(fp$f2, z)
  if (max(abs(cessation_all_times - cessation_all_times[1])) > 1e-8) {
    stop("Cessation female log-HR unexpectedly varied over time")
  }
  values <- c(graduation, cessation_all_times[1])
  names(values) <- METRICS
  values
}
warm <- fit_pair(set_a0(impute_one(copy(d0), 902607)))
INIT1 <- warm$init1; INIT2 <- warm$init2
clusters <- split(seq_len(nrow(d0)), d0$cluster); CL <- names(clusters)

obs_cells <- list(); obs_errors <- list()
for (m in seq_len(M_OBS)) {
  fit <- tryCatch(estimate_hr(impute_one(copy(d0), 902607 + m)), error = function(e) e)
  if (inherits(fit, "error")) obs_errors[[length(obs_errors) + 1L]] <- conditionMessage(fit)
  else obs_cells[[length(obs_cells) + 1L]] <- as.data.table(c(list(m = m), as.list(fit)))
}
if (length(obs_errors)) stop(paste(unlist(obs_errors), collapse = " | "))
obs_dt <- rbindlist(obs_cells)
obs_log <- unlist(obs_dt[, lapply(.SD, mean), .SDcols = METRICS])
fwrite(obs_dt, "data/observed_cause_specific_HR_cells.csv")

boot_rep <- function(b) {
  set.seed(902607 + b * 1000)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- copy(d0[idx]); apply_scale(db); fix_levels(db)
  ests <- list(); errors <- character()
  for (m in seq_len(M_BOOT)) {
    fit <- tryCatch(estimate_hr(impute_one(copy(db), 902607 + b * 1000 + m), INIT1, INIT2), error = function(e) e)
    if (inherits(fit, "error")) errors <- c(errors, conditionMessage(fit)) else ests[[length(ests) + 1L]] <- fit
  }
  if (!length(ests)) return(data.table(b = b, nfail = M_BOOT, error = first_error(errors)))
  pooled <- colMeans(do.call(rbind, ests))
  as.data.table(c(list(b = b, nfail = M_BOOT - length(ests)), as.list(pooled),
                  list(error = first_error(errors))))
}
old <- part_reader(PART_DIR)
done_b <- if (is.null(old)) integer(0) else old$b
todo <- setdiff(1:B, done_b)
cat(sprintf("Cause-specific HR bootstrap: %d/%d done, M_BOOT=%d\n", length(done_b), B, M_BOOT))
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  res <- rbindlist(lapply(part, boot_rep), fill = TRUE)
  fwrite(res, file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1])))
  cat(sprintf("Cause-specific HR reps %d-%d written\n", part[1], part[length(part)])); flush.console()
}
br_all <- part_reader(PART_DIR)
stopifnot(!is.null(br_all), nrow(br_all) == B, setequal(br_all$b, 1:B))
fwrite(br_all, "data/bootstrap_cause_specific_HR.csv")
rows <- list()
for (metric in METRICS) {
  cause <- if (startsWith(metric, "grad")) "graduation" else "cessation"
  time <- if (startsWith(metric, "grad")) as.integer(sub("grad_logHR_t", "", metric)) else NA_integer_
  rows[[length(rows) + 1L]] <- data.table(
    cause = cause, time_years = time, logHR = obs_log[[metric]], HR = exp(obs_log[[metric]]),
    HR_lo95 = exp(ci(br_all[[metric]])[1]), HR_hi95 = exp(ci(br_all[[metric]])[2]))
}
hr <- rbindlist(rows)
fwrite(hr, "data/results_cause_specific_HR.csv")
out <- list(analyst = "ZZ H", date = "2026-11-16", N = nrow(d0),
            calibration = cal$selected, calibration_N = nrow(rs),
            M_obs = M_OBS, M_boot = M_BOOT, bootstrap_requested = B,
            bootstrap_valid = sum(is.finite(br_all$cess_logHR_constant)),
            failed_cells = sum(br_all$nfail), results = hr,
            interpretation = "Cause-specific female-versus-male HRs. Graduation HRs are time-specific; cessation HR is constant. These are not risk ratios or mechanistic effects.")
saveRDS(out, "data/results_cause_specific_HR.rds")
write_json(out, "data/results_cause_specific_HR.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(hr)
