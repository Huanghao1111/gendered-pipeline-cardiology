#!/usr/bin/env Rscript
source("code/descriptive/common.R")

region <- function(x) {
  europe <- c("AT", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES",
              "FI", "FR", "GB", "GR", "HR", "HU", "IE", "IS", "IT", "LT",
              "LU", "LV", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK")
  east_asia <- c("CN", "KR", "JP", "TW", "HK")
  x <- toupper(trimws(x))
  ifelse(is.na(x) | x == "" | x == "未知", "Unknown",
         ifelse(x %in% c("US", "CA"), "North America",
                ifelse(x %in% europe, "Europe", ifelse(x %in% east_asia, "East Asia", "Other"))))
}
d0[, inst_region := region(inst_country)]
origin_clean <- toupper(trimws(d0$final_name_origin))
d0[, name_origin_report := fcase(
  origin_clean %in% c("CN", "KR", "JP"), origin_clean,
  origin_clean == "U", "U",
  is.na(origin_clean) | origin_clean == "" | origin_clean == "UNSPECIFIED", "Unspecified",
  default = "other"
)]
d0[, expected_female_contribution := fifelse(!is.na(female_gold), female_gold, p_cal)]
d0[, definite_label := female_gold %in% c(0, 1)]
d0[, high_posterior := definite_label | (p_cal >= 0.90 | p_cal <= 0.10)]

## Overall gender composition and posterior-count convention.
gender_composition <- data.table(
  metric = c("N_total", "definite_F", "definite_M", "definite_F_or_M", "no_definite_F_or_M",
             "no_name_score", "no_name_score_and_no_definite_label", "high_posterior_total",
             "high_pcal_non_adjudicated", "posterior_expected_women", "posterior_expected_women_percent"),
  value = c(nrow(d0), sum(d0$gender_gold == "F"), sum(d0$gender_gold == "M"),
            sum(d0$definite_label), sum(!d0$definite_label), sum(is.na(d0$p_female)),
            sum(is.na(d0$p_female) & !d0$definite_label), sum(d0$high_posterior),
            sum(!d0$definite_label & (d0$p_cal >= 0.90 | d0$p_cal <= 0.10)),
            sum(d0$expected_female_contribution),
            100 * mean(d0$expected_female_contribution))
)
fwrite(gender_composition, "data/gender_composition.csv")

missingness <- rbindlist(list(
  d0[, .(grouping_type = "Overall", grouping = "Primary cohort",
         N = .N, no_score = sum(is.na(p_female)), definite_FM = sum(definite_label),
         score_and_definite = sum(!is.na(p_female) & definite_label),
         high_posterior = sum(high_posterior), expected_women = sum(expected_female_contribution))],
  d0[, .(grouping_type = "inst_region", N = .N, no_score = sum(is.na(p_female)),
         definite_FM = sum(definite_label), score_and_definite = sum(!is.na(p_female) & definite_label),
         high_posterior = sum(high_posterior), expected_women = sum(expected_female_contribution)), by = inst_region][
           , .(grouping_type = "inst_region", grouping = inst_region, N, no_score, definite_FM, score_and_definite, high_posterior, expected_women)],
  d0[, .(grouping_type = "name_origin", N = .N, no_score = sum(is.na(p_female)),
         definite_FM = sum(definite_label), score_and_definite = sum(!is.na(p_female) & definite_label),
         high_posterior = sum(high_posterior), expected_women = sum(expected_female_contribution)), by = name_origin_report][
           , .(grouping_type = "name_origin", grouping = name_origin_report, N, no_score, definite_FM, score_and_definite, high_posterior, expected_women)]
))
missingness[, `:=`(score_missing_percent = 100 * no_score / N,
                   score_coverage_percent = 100 * (N - no_score) / N,
                   definite_FM_percent = 100 * definite_FM / N,
                   expected_women_percent = 100 * expected_women / N,
                   high_posterior_percent = 100 * high_posterior / N)]
fwrite(missingness, "data/regional_gender_missingness.csv")

## Censoring reasons and terminal publication-year distribution.
cens <- d0[event_type == 0]
censoring_summary <- data.table(
  metric = c("N_censored", "last_publication_min_year", "last_publication_max_year",
             "cessation_not_confirmable_by_censor", "graduation_after_censor_excluded", "ever_graduation_through_censor"),
  value = c(nrow(cens), min(cens$last_publication_year_through_2024),
            max(cens$last_publication_year_through_2024),
            sum(cens$last_publication_year_through_2024 >= 2021),
            sum(cens$grad_after_censor), sum(cens$ever_event_A3))
)
fwrite(censoring_summary, "data/censoring_summary.csv")
censoring_by_year <- cens[, .(N = .N, post_censor_graduation_excluded = sum(grad_after_censor),
                              ever_graduation_through_censor = sum(ever_event_A3)),
                          by = last_publication_year_through_2024][order(last_publication_year_through_2024)]
fwrite(censoring_by_year, "data/censoring_by_last_publication_year.csv")

## Follow-up distributions. Sex-specific summaries average imputation-specific quantiles.
d0[, administrative_potential_followup := as.numeric((as.Date("2024-12-31") - as.Date(entry_date)) / 365.25)]
quantiles <- function(x) unname(quantile(x, c(.25, .5, .75), na.rm = TRUE))
overall <- rbind(
  data.table(followup_type = "observed_to_first_event_or_censor", sex = "Overall",
             Q1 = quantiles(d0$time_years)[1], median = quantiles(d0$time_years)[2], Q3 = quantiles(d0$time_years)[3],
             summary_method = "Observed data"),
  data.table(followup_type = "administrative_potential", sex = "Overall",
             Q1 = quantiles(d0$administrative_potential_followup)[1],
             median = quantiles(d0$administrative_potential_followup)[2],
             Q3 = quantiles(d0$administrative_potential_followup)[3], summary_method = "Observed data")
)
cells <- list()
for (m in seq_len(M_OBS)) {
  di <- impute_one(copy(d0), 902607 + m)
  for (sex_value in c(1, 0)) {
    z <- di[female == sex_value]
    cells[[length(cells) + 1L]] <- data.table(
      m = m, sex = ifelse(sex_value == 1, "Female", "Male"), N = nrow(z),
      observed_Q1 = quantiles(z$time_years)[1], observed_median = quantiles(z$time_years)[2], observed_Q3 = quantiles(z$time_years)[3],
      potential_Q1 = quantiles(z$administrative_potential_followup)[1],
      potential_median = quantiles(z$administrative_potential_followup)[2],
      potential_Q3 = quantiles(z$administrative_potential_followup)[3])
  }
}
cells <- rbindlist(cells)
fwrite(cells, "data/followup_by_imputation.csv")
sex_summary <- cells[, .(
  Q1 = mean(observed_Q1), median = mean(observed_median), Q3 = mean(observed_Q3),
  median_min = min(observed_median), median_max = max(observed_median),
  summary_method = "Mean of 20 imputation-specific quantiles"), by = sex]
sex_summary[, followup_type := "observed_to_first_event_or_censor"]
potential_summary <- cells[, .(
  Q1 = mean(potential_Q1), median = mean(potential_median), Q3 = mean(potential_Q3),
  median_min = min(potential_median), median_max = max(potential_median),
  summary_method = "Mean of 20 imputation-specific quantiles"), by = sex]
potential_summary[, followup_type := "administrative_potential"]
followup_summary <- rbind(overall, rbind(sex_summary, potential_summary, fill = TRUE), fill = TRUE)
setcolorder(followup_summary, c("followup_type", "sex", "Q1", "median", "Q3", "median_min", "median_max", "summary_method"))
fwrite(followup_summary, "data/followup_summary.csv")

median_contrast <- cells[, .(female_median = observed_median[sex == "Female"], male_median = observed_median[sex == "Male"]), by = m]
median_contrast[, female_minus_male := female_median - male_median]
fwrite(median_contrast, "data/followup_median_contrast_by_imputation.csv")

qa <- list(
  analyst = "ZZ H", date = "2026-11-16", status = "COMPLETE",
  N = nrow(d0), calibration_N = nrow(rs),
  definite_F = sum(d0$gender_gold == "F"), definite_M = sum(d0$gender_gold == "M"),
  no_definite_label = sum(!d0$definite_label), no_name_score = sum(is.na(d0$p_female)),
  no_score_and_no_definite = sum(is.na(d0$p_female) & !d0$definite_label),
  posterior_expected_women = sum(d0$expected_female_contribution),
  censored = nrow(cens), censoring_last_year_range = range(cens$last_publication_year_through_2024),
  post_censor_graduation_excluded = sum(cens$grad_after_censor),
  high_posterior_N = sum(d0$high_posterior), M_obs = M_OBS)
write_json(qa, "data/descriptive_package_qa.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(gender_composition)
print(missingness)
print(censoring_summary)
print(followup_summary)
