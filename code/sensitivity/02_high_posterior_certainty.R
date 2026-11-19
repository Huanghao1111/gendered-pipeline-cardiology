#!/usr/bin/env Rscript
source("code/sensitivity/common.R")
M_BOOT <- as.integer(Sys.getenv("BOOT_M", "1"))
B <- as.integer(Sys.getenv("BOOT_B", "1000"))
PART_DIR <- "work_tmp/high_posterior_parts"
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)

d0[, high_posterior := female_gold %in% c(0, 1) | (p_cal >= 0.90 | p_cal <= 0.10)]
high_summary <- list(
  N_full = nrow(d0),
  N_adjudicated = sum(d0$female_gold %in% c(0, 1)),
  N_high_pcal = sum(!(d0$female_gold %in% c(0, 1)) & (d0$p_cal >= 0.90 | d0$p_cal <= 0.10)),
  N_no_score = sum(is.na(d0$p_female)),
  N_high = sum(d0$high_posterior),
  N_event1 = sum(d0$high_posterior & d0$ev1 == 1),
  N_event2 = sum(d0$high_posterior & d0$ev2 == 1),
  N_censored = sum(d0$high_posterior & d0$ev1 == 0 & d0$ev2 == 0)
)
d0 <- d0[high_posterior == TRUE]
high_summary$N <- nrow(d0)
stopifnot(high_summary$N > 0)

obs0 <- est_s1s3(impute_one(copy(d0), 902607))
INIT1 <- obs0$init1; INIT2 <- obs0$init2
clusters <- split(seq_len(nrow(d0)), d0$cluster); CL <- names(clusters)
boot_rep <- function(b) {
  set.seed(902607 + b * 1000)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- copy(d0[idx]); apply_scale(db); fix_levels(db)
  ests <- list(); errors <- character()
  for (m in seq_len(M_BOOT)) {
    fit <- tryCatch(est_s1s3(impute_one(copy(db), 902607 + b * 1000 + m), INIT1, INIT2)[[1]],
                    error = function(e) e)
    if (inherits(fit, "error")) errors <- c(errors, conditionMessage(fit)) else ests[[length(ests) + 1L]] <- fit
  }
  if (!length(ests)) return(data.table(b = b, nfail = M_BOOT, error = first_error(errors)))
  pooled <- colMeans(do.call(rbind, ests))
  as.data.table(c(list(b = b, nfail = M_BOOT - length(ests),
                       riskF = pooled["gradF"], riskM = pooled["gradM"],
                       RD = pooled["RDgrad"], RR = pooled["gradF"] / pooled["gradM"],
                       error = first_error(errors))))
}
old <- part_reader(PART_DIR)
done_b <- if (is.null(old)) integer(0) else old$b
todo <- setdiff(1:B, done_b)
cat(sprintf("High-posterior bootstrap: %d/%d done, N=%d, M_BOOT=%d\n", length(done_b), B, nrow(d0), M_BOOT))
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  res <- rbindlist(lapply(part, boot_rep), fill = TRUE)
  fwrite(res, file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1])))
  cat(sprintf("High-posterior reps %d-%d written\n", part[1], part[length(part)])); flush.console()
}
br_all <- part_reader(PART_DIR)
stopifnot(!is.null(br_all), nrow(br_all) == B, setequal(br_all$b, 1:B))
br <- br_all[is.finite(RD)]
obs_cells <- lapply(1:M_OBS, function(m)
  tryCatch(est_s1s3(impute_one(copy(d0), 902607 + m))[[1]], error = function(e) e))
obs_errors <- vapply(obs_cells, inherits, logical(1), what = "error")
if (any(obs_errors)) stop(paste(vapply(obs_cells[obs_errors], conditionMessage, character(1)), collapse = " | "))
em <- do.call(rbind, obs_cells)
obs_est <- c(riskF = mean(em[, "gradF"]), riskM = mean(em[, "gradM"]),
             RD = mean(em[, "RDgrad"]), RR = exp(mean(log(em[, "gradF"] / em[, "gradM"]))))
out <- list(analyst = "ZZ H", date = "2026-11-16",
            estimand = "High-posterior-certainty standardized 15-year graduation RD (F-M)",
            N = nrow(d0), high_summary = high_summary, calibration = cal$selected,
            calibration_N = nrow(rs), brier_isotonic = brier_iso, brier_platt = brier_platt,
            observed = obs_est, RD_ci = ci(br$RD), RR_ci = exp(ci(log(br$RR))),
            riskF_ci = ci(br$riskF), riskM_ci = ci(br$riskM),
            bootstrap_requested = B, bootstrap_valid = nrow(br),
            failed_cells = sum(br_all$nfail), M_obs = M_OBS, M_boot = M_BOOT)
saveRDS(out, "data/results_high_posterior.rds")
write_json(out, "data/results_high_posterior.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(br_all, "data/bootstrap_high_posterior.csv")
print(out)
