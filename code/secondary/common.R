suppressMessages({library(rstpm2); library(data.table); library(yaml); library(parallel); library(jsonlite)})

cfg <- yaml.load_file("config/config.yaml")
stopifnot(isTRUE(cfg$annotation_lock), file.exists("data/annotation_lock.timestamp"))
M_OBS <- 20L
M_BOOT <- as.integer(Sys.getenv("BOOT_M", "2"))
DFB <- cfg$model$baseline_df
DFT <- cfg$model$gender_time_interaction_df
CORES <- 1L
TMAX <- 15
GRID <- 50L
TS <- seq(0.05, TMAX, length.out = GRID)

d0 <- fread("data/analysis_cohort.csv")[primary_cohort == TRUE]
stopifnot(nrow(d0) == 1536L)
d0[, female_gold := fifelse(gender_gold == "F", 1, fifelse(gender_gold == "M", 0, NA_real_))]
d0[, ev1 := as.integer(event_type == 1)]
d0[, ev2 := as.integer(event_type == 2)]
stopifnot(sum(d0$ev1) == 866L, sum(d0$ev2) == 534L,
          sum(d0$ev1 == 0L & d0$ev2 == 0L) == 136L)
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

source("code/secondary/weighted_calibration.R")
rs <- d0[core_selected == 1 & !is.na(p_female) & female_gold %in% c(0, 1)]
stopifnot(nrow(rs) == 495L, !anyNA(rs$design_weight), all(rs$design_weight > 0))
cal <- fit_calibration(rs)
brier_iso <- cal$brier_isotonic
brier_platt <- cal$brier_platt
use_iso <- cal$selected == "isotonic"
p_calibrate <- cal$predict
prev_rs <- cal$prevalence
d0[, p_cal := fifelse(!is.na(p_female), p_calibrate(p_female), prev_rs)]
cat(sprintf("calibration: weighted isotonic=%.4f, Platt=%.4f, selected=%s, N=%d\n",
            brier_iso, brier_platt, cal$selected, nrow(rs)))

impute_one <- function(dd, seed) {
  set.seed(seed)
  dd[, female := fifelse(!is.na(female_gold), female_gold, rbinom(.N, 1, p_cal))]
  dd
}

haz_grids <- function(f1, f2, dd) {
  n <- nrow(dd); G <- length(TS)
  nd0 <- dd[, ..MODC]; nd0[, female := 0]
  long <- nd0[rep(1:n, times = G)]; long[, time_years := rep(TS, each = n)]
  H10 <- matrix(predict(f1, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  H20 <- matrix(predict(f2, newdata = long, type = "cumhaz"), nrow = n, ncol = G)
  one <- nd0[rep(1, 2 * G)]; one[, time_years := rep(TS, times = 2)]
  one[, female := rep(c(0, 1), each = G)]
  d1 <- predict(f1, newdata = one, type = "link")[(G + 1):(2 * G)] -
        predict(f1, newdata = one, type = "link")[1:G]
  d2 <- predict(f2, newdata = one, type = "link")[(G + 1):(2 * G)] -
        predict(f2, newdata = one, type = "link")[1:G]
  list(H1F = sweep(H10, 2, exp(d1), `*`), H1M = H10,
       H2F = sweep(H20, 2, exp(d2), `*`), H2M = H20)
}

occ_from_H <- function(H1, H2) {
  G <- ncol(H1)
  S <- exp(-(H1 + H2))
  Smid <- (cbind(1, S[, -G]) + S) / 2
  d1 <- t(apply(H1, 1, function(x) diff(c(0, x))))
  d2 <- t(apply(H2, 1, function(x) diff(c(0, x))))
  c(grad = mean(rowSums(d1 * Smid)), cess = mean(rowSums(d2 * Smid)),
    active = mean(S[, G]))
}

fit_pair <- function(dd, init1 = NULL, init2 = NULL) {
  f1 <- stpm2(as.formula(paste("Surv(time_years, ev1) ~", RHS)), data = dd,
              df = DFB, tvc = list(female = DFT), init = init1)
  f2 <- stpm2(as.formula(paste("Surv(time_years, ev2) ~", RHS)), data = dd, df = DFB, init = init2)
  list(f1 = f1, f2 = f2, init1 = coef(f1), init2 = coef(f2))
}

est_s1s3 <- function(dd, init1 = NULL, init2 = NULL) {
  fp <- fit_pair(dd, init1, init2)
  hg <- haz_grids(fp$f1, fp$f2, dd)
  oF <- occ_from_H(hg$H1F, hg$H2F)
  oM <- occ_from_H(hg$H1M, hg$H2M)
  list(c(gradF = unname(oF["grad"]), gradM = unname(oM["grad"]),
         RDgrad = unname(oF["grad"] - oM["grad"]),
         cessF = unname(oF["cess"]), cessM = unname(oM["cess"]),
         RDcess = unname(oF["cess"] - oM["cess"]),
         actF = unname(oF["active"]), actM = unname(oM["active"]),
         RDact = unname(oF["active"] - oM["active"]),
         sumF = sum(oF), sumM = sum(oM)),
       init1 = fp$init1, init2 = fp$init2)
}

ci <- function(x, conf = 0.95) {
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  unname(quantile(x, probs = probs, na.rm = TRUE))
}

part_reader <- function(part_dir) {
  pf <- list.files(part_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(pf)) return(NULL)
  parts <- rbindlist(lapply(pf, function(f) tryCatch(fread(f), error = function(e) NULL)), fill = TRUE)
  stopifnot(!anyDuplicated(parts$b))
  parts
}

first_error <- function(x) {
  x <- x[nzchar(x)]
  if (!length(x)) "" else paste(unique(x), collapse = " | ")
}
