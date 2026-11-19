#!/usr/bin/env Rscript
source("code/sensitivity/common.R")
M_BOOT <- as.integer(Sys.getenv("BOOT_M", "1"))
B <- as.integer(Sys.getenv("BOOT_B", "1000"))
PART_DIR <- "work_tmp/unadjusted_parts"
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)

fit_unadjusted <- function(dd, init1 = NULL, init2 = NULL) {
  f1 <- stpm2(Surv(time_years, ev1) ~ female, data = dd, df = DFB,
              tvc = list(female = DFT), init = init1)
  f2 <- stpm2(Surv(time_years, ev2) ~ female, data = dd, df = DFB, init = init2)
  list(f1 = f1, f2 = f2, init1 = coef(f1), init2 = coef(f2))
}
est_unadjusted <- function(dd, init1 = NULL, init2 = NULL) {
  fp <- fit_unadjusted(dd, init1, init2)
  hg <- haz_grids(fp$f1, fp$f2, dd)
  oF <- occ_from_H(hg$H1F, hg$H2F); oM <- occ_from_H(hg$H1M, hg$H2M)
  list(c(riskF = unname(oF["grad"]), riskM = unname(oM["grad"]),
         RD = unname(oF["grad"] - oM["grad"]), RR = unname(oF["grad"] / oM["grad"])),
       init1 = fp$init1, init2 = fp$init2)
}
obs0 <- est_unadjusted(impute_one(copy(d0), 902607))
INIT1 <- obs0$init1; INIT2 <- obs0$init2
clusters <- split(seq_len(nrow(d0)), d0$cluster); CL <- names(clusters)
boot_rep <- function(b) {
  set.seed(902607 + b * 1000)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- copy(d0[idx]); apply_scale(db); fix_levels(db)
  ests <- list(); errors <- character()
  for (m in seq_len(M_BOOT)) {
    fit <- tryCatch(est_unadjusted(impute_one(copy(db), 902607 + b * 1000 + m), INIT1, INIT2)[[1]],
                    error = function(e) e)
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
cat(sprintf("Unadjusted bootstrap: %d/%d done, M_BOOT=%d\n", length(done_b), B, M_BOOT))
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  res <- rbindlist(lapply(part, boot_rep), fill = TRUE)
  fwrite(res, file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1])))
  cat(sprintf("Unadjusted reps %d-%d written\n", part[1], part[length(part)])); flush.console()
}
br_all <- part_reader(PART_DIR)
stopifnot(!is.null(br_all), nrow(br_all) == B, setequal(br_all$b, 1:B))
br <- br_all[is.finite(RD)]
obs_cells <- lapply(1:M_OBS, function(m)
  tryCatch(est_unadjusted(impute_one(copy(d0), 902607 + m))[[1]], error = function(e) e))
obs_errors <- vapply(obs_cells, inherits, logical(1), what = "error")
if (any(obs_errors)) stop(paste(vapply(obs_cells[obs_errors], conditionMessage, character(1)), collapse = " | "))
obs_est <- colMeans(do.call(rbind, obs_cells)
)
out <- list(analyst = "ZZ H", date = "2026-11-16",
            estimand = "Unadjusted model-based 15-year graduation RD (F-M); no covariates except gender",
            N = nrow(d0), calibration = cal$selected, calibration_N = nrow(rs),
            observed = obs_est, RD_ci = ci(br$RD), RR_ci = exp(ci(log(br$RR))),
            riskF_ci = ci(br$riskF), riskM_ci = ci(br$riskM),
            bootstrap_requested = B, bootstrap_valid = nrow(br),
            failed_cells = sum(br_all$nfail), M_obs = M_OBS, M_boot = M_BOOT)
saveRDS(out, "data/results_unadjusted_RD.rds")
write_json(out, "data/results_unadjusted_RD.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(br_all, "data/bootstrap_unadjusted_RD.csv")
print(out)
