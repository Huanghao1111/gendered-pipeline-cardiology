#!/usr/bin/env Rscript
suppressMessages(library(data.table))
suppressMessages(library(jsonlite))
source("code/descriptive/common.R")

OUT <- "data"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B <- 1000L
M_OBS_LOCAL <- 20L
d0[, administrative_potential_followup := as.numeric((as.Date("2024-12-31") - as.Date(entry_date)) / 365.25)]
clusters <- split(seq_len(nrow(d0)), d0$cluster)
CL <- names(clusters)

set.seed(20261114)
reps <- vector("list", B)
for (b in seq_len(B)) {
  m <- (b - 1L) %% M_OBS_LOCAL + 1L
  di <- impute_one(copy(d0), 902607 + m)
  idx <- unlist(clusters[sample(CL, length(CL), replace = TRUE)], use.names = FALSE)
  db <- di[idx]
  obsF <- median(db[female == 1, time_years])
  obsM <- median(db[female == 0, time_years])
  potF <- median(db[female == 1, administrative_potential_followup])
  potM <- median(db[female == 0, administrative_potential_followup])
  reps[[b]] <- data.table(
    b = b, m = m,
    observed_female = obsF, observed_male = obsM,
    observed_diff = obsF - obsM,
    potential_female = potF, potential_male = potM,
    potential_diff = potF - potM,
    N_female = db[female == 1, .N], N_male = db[female == 0, .N],
    N_total = nrow(db)
  )
}
reps <- rbindlist(reps)
stopifnot(nrow(reps) == B, all(is.finite(reps$observed_diff)), all(is.finite(reps$potential_diff)))
fwrite(reps, file.path(OUT, "followup_bootstrap_replicates.csv"))

observed_estimates <- rbindlist(lapply(seq_len(M_OBS_LOCAL), function(m) {
  di <- impute_one(copy(d0), 902607 + m)
  data.table(
    m = m,
    observed_female = median(di[female == 1, time_years]),
    observed_male = median(di[female == 0, time_years]),
    potential_female = median(di[female == 1, administrative_potential_followup]),
    potential_male = median(di[female == 0, administrative_potential_followup])
  )
}))
observed_estimates[, `:=`(
  observed_diff = observed_female - observed_male,
  potential_diff = potential_female - potential_male
)]
fwrite(observed_estimates, file.path(OUT, "followup_estimates_by_imputation.csv"))
stopifnot(all(is.finite(observed_estimates$observed_diff)),
          all(is.finite(observed_estimates$potential_diff)),
          nrow(observed_estimates) == M_OBS_LOCAL)

ci95 <- function(x) unname(quantile(x, c(.025, .975), na.rm = TRUE))
p_two <- function(x) {
  n <- length(x)
  min(1, 2 * min((1 + sum(x <= 0)) / (n + 1), (1 + sum(x >= 0)) / (n + 1)))
}
make_row <- function(label, estimate, boot_values) {
  data.table(
    contrast = label,
    point_estimate = estimate,
    ci95_low = ci95(boot_values)[1], ci95_high = ci95(boot_values)[2],
    bootstrap_p_two_sided = p_two(boot_values),
    bootstrap_median = median(boot_values),
    bootstrap_q1 = quantile(boot_values, .25), bootstrap_q3 = quantile(boot_values, .75),
    valid_bootstrap = length(boot_values)
  )
}
summary <- rbind(
  make_row("Observed time to first event or administrative censoring, female minus male median, years",
           mean(observed_estimates$observed_diff), reps$observed_diff),
  make_row("Administrative potential follow-up, female minus male median, years",
           mean(observed_estimates$potential_diff), reps$potential_diff)
)
fwrite(summary, file.path(OUT, "followup_bootstrap_summary.csv"))

qa <- list(
  analyst = "ZZ H",
  date = "2026-11-16",
  status = "EXPLORATORY_FOLLOWUP_MEDIAN_CONTRAST_COMPLETE",
  N = nrow(d0),
  institution_clusters = length(CL),
  bootstrap_replicates = B,
  imputations = M_OBS_LOCAL,
  imputation_rotation = "m = (b - 1) %% 20 + 1",
  seed = 20261114,
  method = "Institution-cluster bootstrap with calibrated Bernoulli gender imputation; unadjusted median differences, percentile 95% CI and sign-based two-sided bootstrap p value",
  observed_point_estimate = mean(observed_estimates$observed_diff),
  potential_point_estimate = mean(observed_estimates$potential_diff),
  counts_by_imputation = reps[, .N, by = m][order(m)][, setNames(as.list(N), paste0("m_", m))],
  checks = list(
    all_replicates_finite = TRUE,
    imputation_counts_equal = all(reps[, .N, by = m]$N == B / M_OBS_LOCAL),
    event_source = "Derived A0 terminal events"
  ),
  interpretation_boundary = "Observed time is jointly determined by graduation, cessation and administrative censoring. This is an exploratory distributional contrast, not a test of follow-up adequacy or of a sex-by-time interaction in the outcome models. Administrative potential follow-up reflects entry dates and the common administrative end."
)
write_json(qa, file.path(OUT, "followup_bootstrap_QA.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(summary)
