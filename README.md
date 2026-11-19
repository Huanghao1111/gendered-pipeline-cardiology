# From first author to senior author: the gendered pipeline in three top-tier cardiology journals

Analysis code for a longitudinal cohort study of first-author entry and subsequent senior authorship in the *European Heart Journal*, *Circulation*, and *JACC*.

## Study overview

The cohort includes researchers whose first qualifying first-author original article appeared in one of the three journals during 2005–2010. The primary cohort comprises entrants during 2005–2009, with follow-up through 31 December 2024.

Graduation is the first subsequent PubMed-indexed last-author original article with 2–30 authors and no consortium byline. Publication cessation is defined as an adjacent publication-year difference of at least 5 years, equivalent to four complete empty calendar years.

The primary estimand is the covariate-standardized 15-year graduation risk difference, women minus men, estimated using cause-specific Royston–Parmar models and g-computation. Secondary and sensitivity analyses assess three-state outcomes, journal tiers, ever-graduation, re-entry, cause-specific hazards, alternative gender classifications, covariate windows, and calibration procedures.

## Repository structure

- `code/01_clean_articles.py` to `code/06_preentry_pull.py`: article cleaning, cohort construction, gender inference, reference sampling, and metadata extraction.
- `code/cohort_derivation/`: author linkage, covariate derivation, and event construction.
- `code/descriptive/`: cohort characteristics and follow-up summaries.
- `code/primary/`: gender calibration and the primary competing-event analysis.
- `code/secondary/`: secondary outcome analyses.
- `code/sensitivity/`: sensitivity analyses and synthetic-data tests.
- `config/`: analysis settings and blank adjudication templates.

## Data availability

This repository contains analysis code and blank templates only. Frozen raw OpenAlex extracts, the study protocol, statistical analysis plan, and deviation log are available through the [OSF project](https://osf.io/52wcp/).

Cleaned cohorts, author identifiers, manual labels, adjudication records, individual-level analysis datasets, model outputs, and figure-production files are not distributed. OpenAlex metadata are subject to the provider's terms. API keys and contact details are not stored in this repository.

## Requirements

```bash
python -m pip install -r requirements.txt
```

```r
install.packages(c("rstpm2", "data.table", "yaml", "jsonlite"))
```

The R scripts also use the base `parallel` package.

## Analysis workflow

Run commands from the repository root. Local input files should be placed in `data/`; the required locked inputs and environment variables are documented in the scripts.

### 1. Construct the cohort

```bash
python code/01_clean_articles.py
python code/02_build_cohort.py
python code/03_gender_inference.py
python code/04_reference_standard.py
python code/05_outcome_pull.py
python code/06_preentry_pull.py
```

### 2. Derive outcomes and analysis variables

```bash
python code/cohort_derivation/01_link_entry_anchored_records.py
python code/cohort_derivation/02_build_analysis_table.py
python code/cohort_derivation/03_derive_events.py
```

OpenAlex `article` records are screened programmatically. Title patterns flag records for manual review rather than automatically excluding them. Final article-eligibility decisions are read from a local adjudication table based on the blank template in `config/`.

### 3. Run the main analyses

```bash
Rscript code/primary/01_gender_validation.R
Rscript code/primary/01b_primary_calibration_cv.R
Rscript code/descriptive/01_descriptive_characteristics.R
Rscript code/descriptive/02_followup_distribution.R
Rscript code/primary/02_primary_analysis.R
Rscript code/secondary/01_three_state_outcomes.R
Rscript code/secondary/02_journal_tiers.R
Rscript code/secondary/03_ever_graduation_and_reentry.R
Rscript code/secondary/04_cause_specific_hazards.R
```

### 4. Run sensitivity analyses

```bash
Rscript code/sensitivity/01_definite_label_ipw_completecase.R
Rscript code/sensitivity/02_high_posterior_certainty.R
Rscript code/sensitivity/03_unadjusted_model.R
Rscript code/sensitivity/04_index_inclusive_window.R
Rscript code/sensitivity/05_fixed_covariate_window.R
Rscript code/sensitivity/06_full_cohort.R
Rscript code/sensitivity/07_full_cohort_full_calibration.R
Rscript code/sensitivity/08_cessation_time_interaction.R
Rscript code/sensitivity/09_calibration_refit_bootstrap.R
```

Covariate-window inputs can be reconstructed with `rebuild_covariate_window_metrics.py` and `prepare_covariate_windows.py` in `code/sensitivity/`. These scripts require the local source-run files described in their headers.

## Project links

- [OSF project](https://osf.io/52wcp/)
- [Preregistration](https://osf.io/cnxs3/)

## License

Code is released under the MIT License. Source metadata remain subject to the terms of their respective providers.
