#!/usr/bin/env Rscript
source("code/secondary/common.R")
B <- as.integer(Sys.getenv("BOOT_B2", "200"))
PART_DIR <- "work_tmp/s2_parts"
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)
RHS_TIER <- "female + ey + index_journal + academic_age + lp + lc + lh + lco + li"

tier_est <- function(dd) {
  fits <- list()
  for (tk in c("tier1", "tier2", "tier3")) {
    dd[, evk := as.integer(event_type == 1 & grad_tier == tk)]
    fits[[tk]] <- stpm2(as.formula(paste("Surv(time_years, evk) ~", RHS_TIER)), data = dd, df = DFB)
  }
  dd[, evc := as.integer(event_type == 2)]
  fc <- stpm2(as.formula(paste("Surv(time_years, evc) ~", RHS_TIER)), data = dd, df = DFB)
  n <- nrow(dd); G <- length(TS); cifs <- list()
  for (g in c(1, 0)) {
    nd0 <- dd[, ..MODC]; nd0[, female := g]
    long <- nd0[rep(1:n, times = G)]; long[, time_years := rep(TS, each = n)]
    Hs <- lapply(c(fits, list(cess = fc)), function(f)
      matrix(predict(f, newdata = long, type = "cumhaz"), nrow = n, ncol = G))
    Htot <- Reduce(`+`, Hs)
    S <- exp(-Htot)
    Smid <- (cbind(1, S[, -G]) + S) / 2
    for (nm in names(Hs)[1:3]) {
      dk <- t(apply(Hs[[nm]], 1, function(x) diff(c(0, x))))
      cifs[[paste0(nm, "_", g)]] <- mean(rowSums(dk * Smid))
    }
  }
  c(t1F = cifs[["tier1_1"]], t1M = cifs[["tier1_0"]],
    t2F = cifs[["tier2_1"]], t2M = cifs[["tier2_0"]],
    t3F = cifs[["tier3_1"]], t3M = cifs[["tier3_0"]])
}
tier_est_chk <- function(dd) {
  v <- tier_est(dd)
  if (is.null(v) || length(v) != 6L || any(!is.finite(v))) stop("tier_est invalid return")
  v
}
clusters <- split(seq_len(nrow(d0)), d0$cluster); CL <- names(clusters)
boot_tier <- function(b) {
  set.seed(922607 + b * 1000)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- copy(d0[idx]); apply_scale(db); fix_levels(db)
  ests <- list(); errors <- character()
  for (m in seq_len(M_BOOT)) {
    fit <- tryCatch(tier_est_chk(impute_one(copy(db), 922607 + b * 1000 + m)), error = function(e) e)
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
cat(sprintf("S2 bootstrap: %d/%d done, M_BOOT=%d\n", length(done_b), B, M_BOOT))
for (s in if (length(todo)) seq(1, length(todo), by = 10) else integer(0)) {
  part <- todo[s:min(s + 9, length(todo))]
  t0 <- Sys.time()
  res <- rbindlist(lapply(part, boot_tier), fill = TRUE)
  fwrite(res, file.path(PART_DIR, sprintf("chunk_%04d.csv", part[1])))
  cat(sprintf("S2 reps %d-%d written in %.1fs\n", part[1], part[length(part)],
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  flush.console()
}
br_all <- part_reader(PART_DIR)
stopifnot(!is.null(br_all), nrow(br_all) == B, setequal(br_all$b, 1:B))
br <- br_all[is.finite(t1F) & is.finite(t2F) & is.finite(t3F)]
obs_cells <- lapply(1:M_OBS, function(m)
  tryCatch(tier_est_chk(impute_one(copy(d0), 902607 + m)), error = function(e) e))
obs_errors <- vapply(obs_cells, inherits, logical(1), what = "error")
if (any(obs_errors)) stop(paste(vapply(obs_cells[obs_errors], conditionMessage, character(1)), collapse = " | "))
obs_t <- colMeans(do.call(rbind, obs_cells))
s2 <- data.table(tier = c("tier1", "tier2", "tier3"),
                 riskF = c(obs_t["t1F"], obs_t["t2F"], obs_t["t3F"]),
                 riskM = c(obs_t["t1M"], obs_t["t2M"], obs_t["t3M"]))
s2[, RD_pp := 100 * (riskF - riskM)]
for (prefix in c("t1", "t2", "t3")) {
  delta <- br[[paste0(prefix, "F")]] - br[[paste0(prefix, "M")]]
  s2[tier == sub("t", "tier", prefix), RD_lo95 := 100 * ci(delta)[1]]
  s2[tier == sub("t", "tier", prefix), RD_hi95 := 100 * ci(delta)[2]]
  s2[tier == sub("t", "tier", prefix), RD_lo_bonf := 100 * ci(delta, 1 - 0.05 / 3)[1]]
  s2[tier == sub("t", "tier", prefix), RD_hi_bonf := 100 * ci(delta, 1 - 0.05 / 3)[2]]
}
out <- list(analyst = "ZZ H", date = "2026-11-16", N = nrow(d0),
            M_obs = M_OBS, M_boot = M_BOOT, bootstrap_requested = B,
            bootstrap_valid = nrow(br), failed_cells = sum(br_all$nfail),
            failure_rate = sum(br_all$nfail) / (B * M_BOOT),
            adjustment = "Reduced tier adjustment: female, entry year, index journal, academic age and pre-entry metrics; no inst_tier, country or no_preentry_work",
            intervals = list(descriptive = 0.95, bonferroni_three_tests = 1 - 0.05 / 3),
            tiers = s2)
saveRDS(out, "data/results_S2_tiers.rds")
write_json(out, "data/results_S2_tiers.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(s2, "data/results_S2_tiers.csv")
fwrite(br_all, "data/bootstrap_S2_tiers.csv")
print(s2)
