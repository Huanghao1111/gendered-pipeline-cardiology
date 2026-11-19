#!/usr/bin/env Rscript
source("code/secondary/common.R")
B <- as.integer(Sys.getenv("BOOT_B1", "500"))
PART_DIR <- "work_tmp/s1s3_parts"
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)
prim <- readRDS("data/results_newref_primary.rds")
GATE_OPEN <- isTRUE(prim$gate_open)
stopifnot(GATE_OPEN)

obs0 <- est_s1s3(impute_one(copy(d0), 902607))
INIT1 <- obs0$init1; INIT2 <- obs0$init2
clusters <- split(seq_len(nrow(d0)), d0$cluster); CL <- names(clusters)
boot_rep <- function(b) {
  set.seed(912607 + b * 1000)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- copy(d0[idx]); apply_scale(db); fix_levels(db)
  ests <- list(); errors <- character()
  for (m in seq_len(M_BOOT)) {
    fit <- tryCatch(est_s1s3(impute_one(copy(db), 912607 + b * 1000 + m), INIT1, INIT2)[[1]],
                    error = function(e) e)
    if (inherits(fit, "error")) errors <- c(errors, conditionMessage(fit)) else ests[[length(ests) + 1L]] <- fit
  }
  if (!length(ests)) return(data.table(b = b, nfail = M_BOOT, error = first_error(errors)))
  as.data.table(c(list(b = b, nfail = M_BOOT - length(ests)),
                  as.list(colMeans(do.call(rbind, ests))),
                  list(error = first_error(errors))))
}

old <- part_reader(PART_DIR)
done_b <- if (is.null(old)) integer(0) else old$b
todo <- setdiff(1:B, done_b)
cat(sprintf("S1/S3 bootstrap: %d/%d done, M_BOOT=%d\n", length(done_b), B, M_BOOT))
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  t0 <- Sys.time()
  res <- rbindlist(lapply(part, boot_rep), fill = TRUE)
  fwrite(res, file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1])))
  cat(sprintf("S1/S3 reps %d-%d written in %.1fs\n", part[1], part[length(part)],
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  flush.console()
}

br_all <- part_reader(PART_DIR)
stopifnot(!is.null(br_all), nrow(br_all) == B, setequal(br_all$b, 1:B))
br <- br_all[is.finite(RDcess)]
obs_cells <- lapply(1:M_OBS, function(m)
  tryCatch(est_s1s3(impute_one(copy(d0), 902607 + m))[[1]], error = function(e) e))
obs_errors <- vapply(obs_cells, inherits, logical(1), what = "error")
if (any(obs_errors)) stop(paste(vapply(obs_cells[obs_errors], conditionMessage, character(1)), collapse = " | "))
oe <- colMeans(do.call(rbind, obs_cells))

s1 <- list(
  RDcess_pp = 100 * unname(oe["RDcess"]),
  RDcess_ci95_pp = 100 * ci(br$RDcess),
  cessF = unname(oe["cessF"]), cessM = unname(oe["cessM"]),
  primary_gate_open = GATE_OPEN)
s1$exit_channel_significant <- GATE_OPEN &&
  ((s1$RDcess_ci95_pp[1] > 0) || (s1$RDcess_ci95_pp[2] < 0))
s3 <- list(
  occ15_F = c(graduated = unname(oe["gradF"]), active_ungraduated = unname(oe["actF"]),
              inactive = unname(oe["cessF"])),
  occ15_M = c(graduated = unname(oe["gradM"]), active_ungraduated = unname(oe["actM"]),
              inactive = unname(oe["cessM"])),
  RDactive_ci95 = ci(br$RDact),
  occupancy_sum_F = unname(oe["sumF"]), occupancy_sum_M = unname(oe["sumM"]))
out <- list(analyst = "ZZ H", date = "2026-11-16", N = nrow(d0),
            M_obs = M_OBS, M_boot = M_BOOT, bootstrap_requested = B,
            bootstrap_valid = nrow(br), failed_cells = sum(br_all$nfail),
            failure_rate = sum(br_all$nfail) / (B * M_BOOT),
            calibration = cal$selected, brier_isotonic = brier_iso, brier_platt = brier_platt,
            observed = oe, S1 = s1, S3 = s3)
saveRDS(out, "data/results_S1_S3.rds")
write_json(out, "data/results_S1_S3.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(br_all, "data/bootstrap_S1_S3.csv")
cat(sprintf("S1 RDcess %+.3f pp [%.3f, %.3f], significant=%s; S2 formal gate=%s\n",
            s1$RDcess_pp, s1$RDcess_ci95_pp[1], s1$RDcess_ci95_pp[2],
            s1$exit_channel_significant, s1$exit_channel_significant))
print(s3)
