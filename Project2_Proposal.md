---
editor_options: 
  markdown: 
    wrap: 72
---

# Project 2 Proposal: Longitudinal CEA Trajectories and Disease Progression Prediction in Colorectal Cancer

**Yanan Fang** **Biostat 629 — Case Studies in Health Big Data** **April
2026**

------------------------------------------------------------------------

## 1. Motivation and Differentiation from Project 1

| Dimension | Project 1 | Project 2 (Proposed) |
|------------------------|------------------------|------------------------|
| Cancer type | NSCLC | **Colorectal Cancer (CRC)** |
| Outcome | Overall survival (time-to-event) | **Disease progression (binary, landmark-based)** |
| Data structure | Cross-sectional (one row per patient) | **Longitudinal** (repeated CEA measurements + progression assessments) |
| Core methods | Cox / LASSO Cox / Survival tree | Mixed-effects logistic, Penalized GEE, Brier Score, Conformal prediction sets |
| Validation | C-index only | **AUC + Brier Score + Calibration** |
| Missing data | Complete-case deletion | **Multiple imputation where appropriate** |

### Why CRC + CEA?

-   CRC is the second largest cancer type in MSK-CHORD (**5,543
    patients**).
-   CEA (carcinoembryonic antigen) is the standard monitoring biomarker
    for CRC, with rich longitudinal data (**\~5,471 CRC patients with
    CEA; median \~12 measurements per patient**).
-   MSI status (stable vs. instable) is a key molecular marker in CRC
    and can serve as an effect modifier.
-   Progression data is abundant, enabling a landmark prediction design.

------------------------------------------------------------------------

## 2. Research Question

> In the MSK-CHORD colorectal cancer cohort, can historical CEA
> trajectory features — measured up to a given clinical time point —
> predict disease progression within the subsequent 6 months?

**Hypothesis:** Patients with rapidly rising or highly variable CEA
trajectories have a higher probability of near-term disease progression,
independent of standard clinical and genomic features.

------------------------------------------------------------------------

## 3. Course Methods Applied

| Course Chapter | Method | Role in This Project |
|------------------------|------------------------|------------------------|
| **Ch1** Mixed Models | Mixed-effects logistic regression | Primary prediction model with random patient intercept to handle within-patient correlation across landmarks |
| **Ch2** Missing Data | Explicit handling of `Unknown` + targeted multiple imputation (MICE) | Treat `PRIOR_MED_TO_MSK` as a 3-level clinical variable (`No`, `Yes`, `Unknown`); reserve MICE for true NA values in other covariates if needed |
| **Ch4** Sparse Regression | Penalized GEE (SCAD) | Extended model: variable selection when adding genomic covariates (gene mutation indicators) |
| **Ch5-6** Model Validation | ROC/AUC + Brier Score + Calibration | Comprehensive model evaluation on held-out test set |
| **Ch9** Conformal Prediction | Split conformal prediction sets | Distribution-free uncertainty quantification for predicted progression status |

Supplementary references: - Raschka (2018) — model evaluation framework,
train/test discipline - ISLR Ch6 — sparse regression / regularization
background

------------------------------------------------------------------------

## 4. Data Sources

All data from **MSK-CHORD 2024** release:

| File | Content | Key Variables |
|------------------------|------------------------|------------------------|
| `data_clinical_patient.txt` | Patient-level demographics and outcomes | Age, sex, stage, smoking, prior treatment, OS |
| `data_clinical_sample.txt` | Sample-level metadata | Cancer type, TMB, MSI type, tumor purity |
| `data_timeline_cea_labs.txt` | Longitudinal CEA measurements | Patient ID, date (`START_DATE`), CEA result |
| `data_timeline_progression.txt` | Longitudinal progression assessments | Patient ID, date, progression (Y/N/Indeterminate) |
| `data_timeline_diagnosis.txt` | Diagnosis dates | Patient ID, `START_DATE` (used to define time origin) |
| `data_timeline_treatment.txt` | Treatment history | Agent, subtype (chemo/targeted/immuno), timing |
| `data_mutations.txt` | Somatic mutation calls | Gene (`Hugo_Symbol`), VAF, read counts |

------------------------------------------------------------------------

## 5. Analysis Plan

### Step 1: Cohort Assembly and Landmark Definition

**Cohort construction:**

1.  **Filter to CRC**: restrict `data_clinical_sample.txt` to
    `CANCER_TYPE == "Colorectal Cancer"`.

2.  **Merge datasets**: link patient clinical data, CEA labs,
    progression records, diagnosis timeline, treatment records, and
    mutation data by `PATIENT_ID`.

3.  **Recover diagnosis-relative time**: in MSK-CHORD, all timeline
    `START_DATE` values are recorded relative to **time zero = first
    tumor sequencing at MSK**. Let `START_DATE_dx` denote the diagnosis
    timeline value from `data_timeline_diagnosis.txt`. For any event
    recorded at raw time `START_DATE_event`, define:

    ```         
    time_from_dx = START_DATE_event - START_DATE_dx
    ```

    so that landmark times, treatment status, and future progression
    windows are all measured consistently from diagnosis.

**Landmark prediction framework:**

Each analysis unit is a **landmark time point** defined by a CEA
measurement occasion. A landmark at time *t* (measured in days from
diagnosis) is eligible for analysis if it satisfies both conditions:

-   **Sufficient history**: the patient has at least 2 prior CEA
    measurements before time *t* (needed to compute slope and
    variability features).
-   **Assessable outcome**: at least one progression assessment exists
    in the window (*t*, *t* + 180 days].

The primary binary outcome at each landmark is:

> Y(t) = 1 if any progression assessment with `PROGRESSION == "Y"`
> occurs in (*t*, *t* + 180]; Y(t) = 0 if at least one assessment with
> `PROGRESSION == "N"` occurs in (*t*, *t* + 180] and no `Y` occurs.

Landmarks are **excluded** from the primary analysis if either:

-   they fail the history requirement,
-   they have no assessment in the 180-day window, or
-   the only assessments in the 180-day window are `Indeterminate`.

This avoids treating unassessable or ambiguous windows as definite
non-events.

**Feature construction — strictly using information up to time *t*:**

| Feature | Definition | Type |
|------------------------|------------------------|------------------------|
| `log_cea_current` | log(CEA) at the landmark time *t* | Time-varying |
| `cea_slope_90d` | Slope of log(CEA) over the 90 days preceding *t* (linear regression on available measurements) | Time-varying |
| `cea_cv` | Coefficient of variation of CEA over the 3 most recent measurements before *t* | Time-varying |
| `time_since_dx` | Days from diagnosis to landmark *t* (`time_from_dx`) | Time-varying |
| `tx_type` | Type of current/most recent systemic treatment at time *t*: chemo, targeted, immuno, none | Time-varying |
| `age` | Age at diagnosis | Time-fixed |
| `sex` | Male vs. Female | Time-fixed |
| `stage` | Highest recorded stage (IV vs. 1-3) | Time-fixed |
| `msi_status` | MSI-Stable vs. MSI-Instable | Time-fixed |
| `tmb` | Tumor mutational burden (nonsynonymous) | Time-fixed |
| `prior_med` | Prior treatment before MSK | Time-fixed |
| Gene indicators | Binary: KRAS, APC, TP53, BRAF, PIK3CA, SMAD4, etc. | Time-fixed |

### Step 2: Missing Data Handling (Ch2)

Missing-data strategy in CRC: - `PRIOR_MED_TO_MSK` is not treated as
ordinary missingness. In the raw data it is recorded as
`No prior medications`, `Prior medications to MSK`, or `Unknown`. The
primary analysis encodes this as a **3-level categorical variable**. -
`TMB_NONSYNONYMOUS`: negligible missingness (\<1%) — complete-case
handling is acceptable. - `MSI_TYPE`: negligible missingness (\<1%) —
complete-case handling is acceptable. - If additional covariates with
**true NA values** remain after cohort construction and have non-trivial
missingness, apply **targeted MICE** (5 imputations; predictive mean
matching for continuous variables, logistic or multinomial models for
categorical variables).

Any MICE procedure is applied **within the training set only** (see Step
4) to prevent information leakage from the test set. The learned
imputation models are then carried forward to the test set.

If MICE is invoked, coefficient estimates from the prediction model are
**pooled across imputations via Rubin's rules**. For performance metrics
(AUC, Brier Score) that do not have simple pooling formulas, we report
the **median and range** across imputations. If no materially missing
covariates remain, the primary analysis proceeds without imputation and
MICE is reported as a sensitivity-analysis framework.

### Step 3: Prediction Models

**Model A — Mixed-effects logistic regression (primary model):**

```         
logit P(Y(t) = 1) = beta_0 + beta_1 * log_cea_current + beta_2 * cea_slope_90d
                   + beta_3 * cea_cv + beta_4 * time_since_dx + beta_5 * tx_type
                   + beta_6 * age + beta_7 * sex + beta_8 * stage + beta_9 * msi_status
                   + beta_10 * prior_med
                   + (1 | patient_id)
```

-   Random intercept accounts for within-patient correlation across
    multiple landmarks.
-   Fitted via `lme4::glmer()` with binomial family.
-   This is the **core model** using clinical + CEA trajectory features
    only.

**Model B — Penalized GEE with genomic covariates (extended model):**

-   Expands Model A's fixed-effect predictor set by adding TMB + binary
    indicators for the top 10 recurrently mutated CRC genes.
-   Uses **Penalized GEE with SCAD penalty** (`PGEE` package) to perform
    variable selection while accounting for within-patient correlation.
-   Working correlation: exchangeable.
-   Tuning parameter selected via BIC or patient-level CV within the
    training set.

The comparison of Model A vs. Model B directly answers: **do genomic
features add predictive value beyond clinical + CEA trajectory
features?**

### Step 4: Validation Framework (Ch5-6)

A single, clearly defined validation pipeline:

```         
Full CRC cohort
  |
  |--- Patient-level random split (70/30)
  |
  v
Training set (70%)                          Test set (30%)
  |                                           |
  |--- If needed: MICE imputation (5 datasets)|--- If used: apply training imputation model
  |--- Fit Models A & B                       |--- Predict on test set (once)
  |--- Internal 5-fold CV within training     |--- Final AUC, Brier Score, Calibration
  |    (patient-level folds, for              |
  |     tuning & internal assessment)         |--- Conformal evaluation (see Step 5)
  |
  |--- For conformal: further split
       training into proper-train (50%)
       + calibration (20%)
```

**Key design decisions:** - **Patient-level splitting**: all landmarks
from the same patient go to the same fold/set. - **Test set used exactly
once** for final performance reporting. - **Internal CV** (within
training set) is used for hyperparameter tuning (PGEE lambda) and
internal model comparison. This does NOT replace the test-set
evaluation. - **Any imputation runs inside the training set** before
model fitting.

**Evaluation metrics on the test set:**

1.  **Discrimination**: ROC curve and AUC on a **patient-independent
    test sample**, defined as one eligible landmark per patient (primary
    analysis; e.g., the first eligible landmark in the test set). As a
    supplementary analysis, also report landmark-level AUC using all
    eligible landmarks with patient-level bootstrap CIs.
2.  **Calibration**: calibration plot — predicted probability vs.
    observed proportion by decile.
3.  **Brier Score**: overall prediction accuracy = (1/n) \* sum(Y_i -
    p_hat_i)\^2.
4.  **Model comparison**: AUC(Model B) vs. AUC(Model A) via DeLong test
    on the patient-independent test sample. As a robustness check,
    compare AUCs via patient-level bootstrap.

### Step 5: Conformal Prediction Sets (Ch9 — optional extension)

Apply split conformal inference to construct **prediction sets** for the
binary progression outcome:

1.  Within the training set, designate a **calibration subset** (the 20%
    carved out in Step 4).
2.  For each calibration observation, compute nonconformity scores. For
    binary classification, a natural score is
    `s_i = 1 - p_hat(Y_i | X_i)` — the model's estimated probability of
    the *wrong* class.
3.  Compute the conformal quantile `q_hat` at level (1 - alpha) from the
    calibration scores.
4.  For each test observation, the **conformal prediction set** at level
    (1 - alpha) includes all labels *y* in {0, 1} for which
    `1 - p_hat(y | x_new) <= q_hat`.

**Interpretation:** - If the model is confident, the prediction set
contains only one label: {0} or {1}. - If the model is uncertain, the
prediction set is {0, 1} — effectively "I don't know." - The coverage
guarantee holds: at least (1 - alpha) of test observations have their
true label in the prediction set.

**What we report:** - Empirical coverage at alpha = 0.10 (target: \>=
90%). - Fraction of prediction sets that are singletons (measure of
informativeness). - Whether singleton fraction differs for progressors
vs. non-progressors.

This is framed as an **extension**, not the core contribution.

------------------------------------------------------------------------

## 6. Expected Deliverables

1.  **R script** (`project2_analysis.R`): fully reproducible analysis
    pipeline
2.  **Report** (\~8 pages): introduction, methods, results with
    tables/figures, discussion
3.  **Presentation** (\~12 slides): for in-class presentation

### Key Figures and Tables (planned)

| \# | Content |
|------------------------------------|------------------------------------|
| Table 1 | CRC cohort characteristics (overall and by progression status at first eligible landmark) |
| Table 2 | Mixed-effects logistic regression: odds ratios, 95% CIs, p-values |
| Table 3 | Penalized GEE: selected variables and coefficients |
| Table 4 | Model comparison on test set: AUC, Brier Score, calibration slope |
| Figure 1 | Spaghetti plot of individual log(CEA) trajectories (sample of \~50 patients), colored by whether they experienced progression |
| Figure 2 | Distribution of CEA 90-day slope for landmarks followed by progression vs. no progression |
| Figure 3 | ROC curves for Model A vs. Model B on test set |
| Figure 4 | Calibration plot for the primary model (Model A) |
| Figure 5 | Conformal prediction sets: singleton fraction by predicted risk decile |

------------------------------------------------------------------------

## 7. Strengths of This Design

1.  **Landmark prediction framework** — avoids the information leakage
    of using future CEA data to predict future events. Each prediction
    uses only information available at the time of the clinical
    decision.
2.  **Principled longitudinal structure** — mixed-effects logistic with
    random intercept handles within-patient correlation, and CEA
    trajectory features (slope, variability) are constructed from
    strictly historical data.
3.  **Clean validation pipeline** — single patient-level train/test
    split, MICE within training only, test set used once. Conformal
    calibration carved from training set.
4.  **Missing-data handling is better aligned with the raw data** —
    `PRIOR_MED_TO_MSK` is treated as a meaningful 3-level clinical
    variable (`No`, `Yes`, `Unknown`), while MICE is reserved for true
    NA values if they remain in other covariates.
5.  **Conformal as principled extension** — correctly framed as
    prediction sets for binary outcomes, not probability intervals.
6.  **Clinical relevance** — CEA is a biomarker oncologists already
    monitor; showing that its *trajectory* (not just current value)
    predicts progression has direct actionable value.

------------------------------------------------------------------------

## 8. Potential Challenges and Mitigations

| Challenge | Mitigation |
|------------------------------------|------------------------------------|
| Informative observation times (sicker patients may have more frequent CEA tests) | Sensitivity analysis: restrict to one landmark per patient (the first eligible), or one per 90-day window |
| Progression labels derived from NLP (noisy) | Sensitivity analysis: restrict to assessments with `NLP_PROGRESSION_PROBABILITY > 0.5`; exclude windows containing only `Indeterminate` labels in the primary analysis |
| Within-patient correlation inflating apparent AUC | Primary AUC uses one eligible landmark per patient; supplementary landmark-level AUC uses patient-level bootstrap |
| PGEE computational cost on large landmark dataset | Use BIC-based lambda selection; if needed, fit on a random subsample of patients |
| Conformal exchangeability assumption | Discuss limitation — landmarks are not fully exchangeable due to temporal ordering within patients; coverage guarantee is approximate |
| Excluding non-assessable landmarks introduces selection | Report the proportion excluded and compare characteristics of included vs. excluded landmarks |
