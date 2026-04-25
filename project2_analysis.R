###############################################################################
## Project 2: Longitudinal CEA Trajectories and Disease Progression
##            Prediction in Colorectal Cancer
## Yanan Fang — Biostat 629, April 2026
###############################################################################

# ── 0. Packages ──────────────────────────────────────────────────────────────
library(data.table)
library(lme4)
library(pROC)
library(ggplot2)
library(glmnet)

data_dir <- "msk_chord_2024"

###############################################################################
## STEP 1: DATA LOADING & COHORT ASSEMBLY
###############################################################################

# ── 1a. Load raw data ────────────────────────────────────────────────────────
read_msk <- function(file) {
  fread(file.path(data_dir, file), sep = "\t", header = TRUE,
        comment.char = "#", na.strings = c("", "NA"))
}

pat   <- read_msk("data_clinical_patient.txt")
samp  <- read_msk("data_clinical_sample.txt")
cea   <- fread(file.path(data_dir, "data_timeline_cea_labs.txt"), sep = "\t")
prog  <- fread(file.path(data_dir, "data_timeline_progression.txt"), sep = "\t")
dx    <- fread(file.path(data_dir, "data_timeline_diagnosis.txt"), sep = "\t")
tx    <- fread(file.path(data_dir, "data_timeline_treatment.txt"), sep = "\t")
mut   <- fread(file.path(data_dir, "data_mutations.txt"), sep = "\t",
               select = c("Hugo_Symbol", "Tumor_Sample_Barcode"))

# ── 1b. Identify CRC patients ───────────────────────────────────────────────
crc_samples <- samp[CANCER_TYPE == "Colorectal Cancer"]
crc_pids    <- unique(crc_samples$PATIENT_ID)
cat("CRC patients:", length(crc_pids), "\n")

# For patients with multiple samples, keep the one with highest coverage
crc_samples <- crc_samples[order(-SAMPLE_COVERAGE)]
crc_samples <- crc_samples[!duplicated(PATIENT_ID)]

# ── 1c. Build patient-level clinical table ───────────────────────────────────
crc_pat <- pat[PATIENT_ID %in% crc_pids,
               .(PATIENT_ID, CURRENT_AGE_DEID, GENDER,
                 STAGE_HIGHEST_RECORDED, SMOKING_PREDICTIONS_3_CLASSES,
                 PRIOR_MED_TO_MSK, OS_MONTHS, OS_STATUS)]

# Merge sample-level features (TMB, MSI)
crc_pat <- merge(crc_pat,
                 crc_samples[, .(PATIENT_ID, SAMPLE_ID, TMB_NONSYNONYMOUS, MSI_TYPE)],
                 by = "PATIENT_ID", all.x = TRUE)

# Clean variable names and coding
setnames(crc_pat, c("CURRENT_AGE_DEID", "GENDER", "STAGE_HIGHEST_RECORDED",
                     "SMOKING_PREDICTIONS_3_CLASSES", "PRIOR_MED_TO_MSK",
                     "TMB_NONSYNONYMOUS", "MSI_TYPE"),
         c("age", "sex", "stage", "smoking", "prior_med", "tmb", "msi_status"))

crc_pat[, sex := factor(sex, levels = c("Female", "Male"))]
crc_pat[, stage := factor(ifelse(stage == "Stage 4", "IV", "I-III"),
                          levels = c("I-III", "IV"))]
crc_pat[, smoking := factor(smoking,
                            levels = c("Former/Current Smoker", "Never", "Unknown"))]
crc_pat[, prior_med := factor(prior_med,
                              levels = c("No prior medications",
                                         "Prior medications to MSK", "Unknown"))]
crc_pat[, msi_status := factor(ifelse(msi_status == "Instable", "MSI-H", "MSS"),
                               levels = c("MSS", "MSI-H"))]

# ── 1d. CRC-specific diagnosis anchoring ──────────────────────────────────────
# FIX: Filter diagnosis records for colorectal-related sites only.
# DX_DESCRIPTION contains ICD topography codes; CRC sites are:
#   C18 (colon), C19 (rectosigmoid), C20 (rectum), C21 (anus/anal canal)
# Also match keywords: COLON, RECTUM, RECTAL, CECUM, SIGMOID, COLORECTAL
crc_dx <- dx[PATIENT_ID %in% crc_pids & SUBTYPE == "Primary"]
crc_dx[, is_crc_dx := grepl("C18|C19|C20|C21|COLON|RECTUM|RECTAL|CECUM|SIGMOID|COLORECT",
                             DX_DESCRIPTION, ignore.case = TRUE)]
crc_dx_filtered <- crc_dx[is_crc_dx == TRUE]

# Patients with CRC-specific dx
cat("Patients with CRC-specific dx:", uniqueN(crc_dx_filtered$PATIENT_ID), "\n")
cat("Patients with only generic dx:", uniqueN(crc_dx[is_crc_dx == FALSE]$PATIENT_ID), "\n")

# Fallback: if no CRC-specific dx found for a patient, use earliest Primary record
# (some descriptions may not contain standard terms)
pids_with_crc_dx <- unique(crc_dx_filtered$PATIENT_ID)
pids_no_crc_dx   <- setdiff(crc_pids, pids_with_crc_dx)
crc_dx_fallback  <- crc_dx[PATIENT_ID %in% pids_no_crc_dx]

crc_dx_all <- rbind(crc_dx_filtered, crc_dx_fallback)
crc_dx_all <- crc_dx_all[order(START_DATE)]
crc_dx_all <- crc_dx_all[!duplicated(PATIENT_ID), .(PATIENT_ID, dx_date = START_DATE)]
cat("Patients with dx anchor after CRC-specific filtering:", nrow(crc_dx_all), "\n")

# ── 1e. Gene mutation indicators ────────────────────────────────────────────
# Map samples to patients
mut_crc <- merge(mut, crc_samples[, .(SAMPLE_ID, PATIENT_ID)],
                 by.x = "Tumor_Sample_Barcode", by.y = "SAMPLE_ID")

# Top 10 genes by patient count
gene_freq <- mut_crc[, .(n_patients = uniqueN(PATIENT_ID)), by = Hugo_Symbol]
gene_freq <- gene_freq[order(-n_patients)]
top_genes <- head(gene_freq$Hugo_Symbol, 10)
cat("Top 10 CRC genes:", paste(top_genes, collapse = ", "), "\n")

# Create binary indicators
gene_mat <- dcast(mut_crc[Hugo_Symbol %in% top_genes],
                  PATIENT_ID ~ Hugo_Symbol,
                  fun.aggregate = function(x) as.integer(length(x) > 0),
                  value.var = "Hugo_Symbol")

# ── 1f. Merge all patient-level data ────────────────────────────────────────
crc_pat <- merge(crc_pat, crc_dx_all, by = "PATIENT_ID", all.x = TRUE)
crc_pat <- merge(crc_pat, gene_mat, by = "PATIENT_ID", all.x = TRUE)

# Fill NA gene indicators with 0 (patient had no mutation in that gene)
for (g in top_genes) {
  set(crc_pat, which(is.na(crc_pat[[g]])), g, 0L)
}

# Remove patients without diagnosis date
cat("Patients without dx date:", sum(is.na(crc_pat$dx_date)), "\n")
crc_pat <- crc_pat[!is.na(dx_date)]
cat("Patients after dx filter:", nrow(crc_pat), "\n")

###############################################################################
## STEP 2: LANDMARK DATASET CONSTRUCTION (vectorized for speed)
###############################################################################

# ── 2a. CEA data for CRC patients ───────────────────────────────────────────
crc_cea <- cea[PATIENT_ID %in% crc_pat$PATIENT_ID]
crc_cea <- merge(crc_cea, crc_pat[, .(PATIENT_ID, dx_date)], by = "PATIENT_ID")
crc_cea[, time_from_dx := START_DATE - dx_date]
crc_cea[, log_cea := log(RESULT + 0.1)]  # small offset to avoid log(0)

# Remove CEA <= 0 or NA
crc_cea <- crc_cea[!is.na(RESULT) & RESULT >= 0]
crc_cea <- crc_cea[order(PATIENT_ID, time_from_dx)]

# Number each CEA within patient
crc_cea[, obs_idx := seq_len(.N), by = PATIENT_ID]
# Patients with < 3 CEA measurements excluded
n_cea_per_pat <- crc_cea[, .N, by = PATIENT_ID]
eligible_pids <- n_cea_per_pat[N >= 3]$PATIENT_ID
crc_cea <- crc_cea[PATIENT_ID %in% eligible_pids]

cat("CEA records for CRC:", nrow(crc_cea), "\n")
cat("Patients with >= 3 CEA:", length(eligible_pids), "\n")

# ── 2b. Progression data for CRC patients ────────────────────────────────────
crc_prog <- prog[PATIENT_ID %in% crc_pat$PATIENT_ID]
crc_prog <- merge(crc_prog, crc_pat[, .(PATIENT_ID, dx_date)], by = "PATIENT_ID")
crc_prog[, time_from_dx := START_DATE - dx_date]

# NOTE on the landmark time t:
#   In the landmark dataset below, `landmark_time` == `time_from_dx` for the
#   CEA draw that defines the landmark. It is the number of days from the CRC
#   primary diagnosis (t = 0 at diagnosis), measured on the same scale as
#   every other timeline variable (progression, treatment). t is therefore a
#   positive integer that varies across landmarks within the same patient;
#   it is NOT zero at the first CEA measurement.

# ── 2b'. Death dates for competing-risk bookkeeping ──────────────────────────
# OS_MONTHS / OS_STATUS in data_clinical_patient give time-from-sequencing, not
# time-from-diagnosis, so we recover a diagnosis-anchored death day only where
# possible. We do NOT change the primary outcome definition; we only *count*
# how many landmarks have a death recorded inside (t, t+180] to quantify the
# competing-risk burden referenced in the report.
#
# Sequencing date (timeline time zero) is day 0 for each patient in the
# timeline files, so a patient's death-from-diagnosis day is:
#   death_from_dx = OS_MONTHS * 30.4375 - dx_date
# where dx_date is measured in days-from-sequencing (negative for diagnoses
# before sequencing, positive after).
crc_death <- crc_pat[!is.na(OS_MONTHS) & OS_STATUS == "1:DECEASED",
                     .(PATIENT_ID,
                       death_from_dx = OS_MONTHS * 30.4375 - dx_date)]

# ── 2c. Treatment data for CRC patients ──────────────────────────────────────
crc_tx <- tx[PATIENT_ID %in% crc_pat$PATIENT_ID]
crc_tx <- merge(crc_tx, crc_pat[, .(PATIENT_ID, dx_date)], by = "PATIENT_ID")
crc_tx[, start_from_dx := START_DATE - dx_date]
crc_tx[, stop_from_dx  := STOP_DATE - dx_date]
# Simplify treatment type
crc_tx[, tx_class := fcase(
  SUBTYPE == "Immuno",       "immuno",
  SUBTYPE == "Targeted",     "targeted",
  SUBTYPE == "Biologic",     "targeted",
  SUBTYPE == "Chemo",        "chemo",
  default = "other"
)]

# ── 2d. Build landmark dataset (vectorized) ─────────────────────────────────
# Landmarks: each CEA measurement from index >= 3 is a potential landmark.
# We only keep it if there is a progression assessment within 180 days.

cat("Building landmark dataset (vectorized)...\n")
t0 <- proc.time()

# Only use CEA observations with index >= 3
lm_base <- crc_cea[obs_idx >= 3, .(PATIENT_ID, landmark_time = time_from_dx,
                                     log_cea_current = log_cea, obs_idx)]

# ── Outcome: non-equi join on progression ──
# For each landmark, check progression in (t, t+180]
lm_base[, t_upper := landmark_time + 180]
prog_joined <- crc_prog[lm_base,
  on = .(PATIENT_ID, time_from_dx > landmark_time, time_from_dx <= t_upper),
  .(PATIENT_ID, landmark_time = i.landmark_time, PROGRESSION, obs_idx = i.obs_idx),
  allow.cartesian = TRUE, nomatch = 0L
]

# Determine outcome per landmark
outcome_dt <- prog_joined[, .(
  has_Y = any(PROGRESSION == "Y"),
  has_N = any(PROGRESSION == "N")
), by = .(PATIENT_ID, landmark_time, obs_idx)]

# Keep landmarks where outcome is determinable
outcome_dt <- outcome_dt[has_Y | has_N]
outcome_dt[, outcome := fifelse(has_Y, 1L, 0L)]
outcome_dt[, c("has_Y", "has_N") := NULL]

# Merge back to get landmark-level identifiers and current CEA
landmarks <- merge(outcome_dt,
                   lm_base[, .(PATIENT_ID, landmark_time, obs_idx, log_cea_current)],
                   by = c("PATIENT_ID", "landmark_time", "obs_idx"))

cat("  Landmarks with assessable outcome:", nrow(landmarks), "\n")

# ── Competing-risk bookkeeping: deaths within the 180-day window ─────────────
# For every landmark with an assessable outcome, flag whether the patient died
# inside (t, t+180]. This is a *diagnostic*, not part of the primary outcome.
landmarks_with_death <- merge(landmarks, crc_death, by = "PATIENT_ID", all.x = TRUE)
landmarks_with_death[, death_in_window := !is.na(death_from_dx) &
                       death_from_dx >  landmark_time &
                       death_from_dx <= landmark_time + 180]
n_death_in_win <- sum(landmarks_with_death$death_in_window)
cat(sprintf(
  "  Deaths recorded inside (t, t+180]: %d / %d landmarks (%.2f%%); among the 0-labelled subset: %d (%.2f%%)\n",
  n_death_in_win, nrow(landmarks_with_death),
  100 * n_death_in_win / nrow(landmarks_with_death),
  sum(landmarks_with_death$death_in_window & landmarks_with_death$outcome == 0),
  100 * sum(landmarks_with_death$death_in_window & landmarks_with_death$outcome == 0) /
    max(sum(landmarks_with_death$outcome == 0), 1)
))
# The second figure answers the reviewer's question: of the landmarks we label
# as "no progression", what fraction actually had a death in the same window?
# If that number is small, miscoding of death-without-assessment as non-event
# has a small effect on the headline AUC.

# ── CEA slope (90d): vectorized per-landmark ──
# For each landmark, fit lm(log_cea ~ time_from_dx) using CEA in [t-90, t]
# Note: both crc_cea and landmarks have 'obs_idx', so rename to avoid collision

lm_for_slope <- landmarks[, .(PATIENT_ID, landmark_time, lm_obs_idx = obs_idx)]

slope_data <- crc_cea[lm_for_slope,
  on = .(PATIENT_ID),
  allow.cartesian = TRUE, nomatch = 0L
][time_from_dx >= landmark_time - 90 & time_from_dx <= landmark_time]

# Need at least 2 observations for slope
slope_counts <- slope_data[, .N, by = .(PATIENT_ID, landmark_time, lm_obs_idx)]
slope_eligible <- slope_counts[N >= 2]

slope_vals <- slope_data[slope_eligible, on = .(PATIENT_ID, landmark_time, lm_obs_idx)][,
  {
    fit <- .lm.fit(cbind(1, time_from_dx), log_cea)
    .(cea_slope_90d = fit$coefficients[2])
  },
  by = .(PATIENT_ID, landmark_time, lm_obs_idx)
]
setnames(slope_vals, "lm_obs_idx", "obs_idx")

# ── CEA CV (last 3 values) ──
# For each landmark at obs_idx k, use observations k-2, k-1, k
lm_with_idx <- landmarks[, .(PATIENT_ID, landmark_time, obs_idx)]
lm_with_idx[, idx_lo := obs_idx - 2L]

cv_data <- crc_cea[lm_with_idx,
  on = .(PATIENT_ID, obs_idx >= idx_lo, obs_idx <= obs_idx),
  .(PATIENT_ID, landmark_time = i.landmark_time, obs_idx = i.obs_idx, RESULT),
  allow.cartesian = TRUE, nomatch = 0L
]

cv_vals <- cv_data[, .(cea_cv = sd(RESULT) / (mean(RESULT) + 0.01)),
                   by = .(PATIENT_ID, landmark_time, obs_idx)]

# ── Treatment status at landmark time ──
# Join treatment data: active if start <= t and (stop >= t or stop is NA)
lm_for_tx <- landmarks[, .(PATIENT_ID, landmark_time, lm_obs_idx = obs_idx)]
tx_at_lm <- crc_tx[lm_for_tx,
  on = .(PATIENT_ID),
  allow.cartesian = TRUE, nomatch = 0L
][start_from_dx <= landmark_time & (is.na(stop_from_dx) | stop_from_dx >= landmark_time)]

# Priority: immuno > targeted > chemo > other
tx_at_lm[, tx_priority := fcase(
  tx_class == "immuno", 1L,
  tx_class == "targeted", 2L,
  tx_class == "chemo", 3L,
  default = 4L
)]

tx_type_dt <- tx_at_lm[, .(tx_type = tx_class[which.min(tx_priority)]),
                        by = .(PATIENT_ID, landmark_time, lm_obs_idx)]
setnames(tx_type_dt, "lm_obs_idx", "obs_idx")

# ── Assemble landmark dataset ──
landmarks <- merge(landmarks, slope_vals,
                   by = c("PATIENT_ID", "landmark_time", "obs_idx"), all.x = TRUE)
landmarks <- merge(landmarks, cv_vals,
                   by = c("PATIENT_ID", "landmark_time", "obs_idx"), all.x = TRUE)
landmarks <- merge(landmarks, tx_type_dt,
                   by = c("PATIENT_ID", "landmark_time", "obs_idx"), all.x = TRUE)

# Fill missing tx_type with "none"
landmarks[is.na(tx_type), tx_type := "none"]

# time_since_dx
landmarks[, time_since_dx := landmark_time]
landmarks[, landmark_id := paste(PATIENT_ID, obs_idx, sep = "__")]

elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Landmark construction took %.1f seconds\n", elapsed))

n_before_slope <- uniqueN(landmarks$PATIENT_ID)
cat("Landmark observations:", nrow(landmarks), "\n")
cat("Patients with landmarks:", n_before_slope, "\n")
cat("Progression rate:", round(mean(landmarks$outcome), 3), "\n")

# ── 2e. Merge patient-level covariates ───────────────────────────────────────
landmarks <- merge(landmarks, crc_pat[, !c("dx_date", "SAMPLE_ID", "OS_MONTHS",
                                           "OS_STATUS")],
                   by = "PATIENT_ID", all.x = TRUE)

landmarks[, tx_type := factor(tx_type,
                              levels = c("none", "chemo", "targeted", "immuno", "other"))]

# Scale time_since_dx to years for interpretability
landmarks[, time_since_dx_yr := time_since_dx / 365.25]

# ── 2f. Remove rows with NA in key features & document attrition ─────────────
n_all    <- nrow(landmarks)
pid_all  <- uniqueN(landmarks$PATIENT_ID)
landmarks <- landmarks[!is.na(cea_slope_90d)]
n_after  <- nrow(landmarks)
pid_after <- uniqueN(landmarks$PATIENT_ID)

cat("\n── Sample attrition due to 90-day slope requirement ──\n")
cat(sprintf("  Before: %d landmarks, %d patients\n", n_all, pid_all))
cat(sprintf("  After:  %d landmarks, %d patients\n", n_after, pid_after))
cat(sprintf("  Dropped: %d landmarks (%.1f%%), %d patients (%.1f%%)\n",
            n_all - n_after, 100 * (n_all - n_after) / n_all,
            pid_all - pid_after, 100 * (pid_all - pid_after) / pid_all))
cat("  NOTE: Remaining cohort is biased toward patients with denser CEA monitoring.\n")

###############################################################################
## STEP 3: TRAIN / TEST SPLIT (patient-level)
###############################################################################

set.seed(629)
all_pids  <- unique(landmarks$PATIENT_ID)
n_train   <- round(0.7 * length(all_pids))
train_pids <- sample(all_pids, n_train)
test_pids  <- setdiff(all_pids, train_pids)

train <- landmarks[PATIENT_ID %in% train_pids]
test  <- landmarks[PATIENT_ID %in% test_pids]

cat("\nTrain:", nrow(train), "landmarks from", length(train_pids), "patients\n")
cat("Test: ", nrow(test),  "landmarks from", length(test_pids),  "patients\n")
cat("Train progression rate:", round(mean(train$outcome), 3), "\n")
cat("Test  progression rate:", round(mean(test$outcome), 3), "\n")

###############################################################################
## STEP 4: DESCRIPTIVE STATISTICS (Table 1 + Figures 1-2)
###############################################################################

# ── Table 1: Cohort characteristics ──────────────────────────────────────────
cat("\n========== TABLE 1: Cohort Characteristics ==========\n")
# Use first eligible landmark per patient for patient-level summary
first_lm <- landmarks[, .SD[which.min(landmark_time)], by = PATIENT_ID]

tab1_vars <- function(dt, label) {
  cat("\n---", label, "(n =", nrow(dt), ") ---\n")
  cat("Age: mean =", round(mean(dt$age, na.rm = TRUE), 1),
      ", SD =", round(sd(dt$age, na.rm = TRUE), 1), "\n")
  cat("Male:", sum(dt$sex == "Male"), "(", round(100 * mean(dt$sex == "Male"), 1), "%)\n")
  cat("Stage IV:", sum(dt$stage == "IV"), "(", round(100 * mean(dt$stage == "IV"), 1), "%)\n")
  cat("MSI-H:", sum(dt$msi_status == "MSI-H", na.rm = TRUE), "(",
      round(100 * mean(dt$msi_status == "MSI-H", na.rm = TRUE), 1), "%)\n")
  cat("Median log(CEA):", round(median(dt$log_cea_current), 2), "\n")
  cat("Median CEA slope (90d):", round(median(dt$cea_slope_90d, na.rm = TRUE), 4), "\n")
  cat("Prior med - No:", sum(dt$prior_med == "No prior medications"), "\n")
  cat("Prior med - Yes:", sum(dt$prior_med == "Prior medications to MSK"), "\n")
  cat("Prior med - Unknown:", sum(dt$prior_med == "Unknown"), "\n")
}

tab1_vars(first_lm, "Overall")
tab1_vars(first_lm[outcome == 1], "Progressors")
tab1_vars(first_lm[outcome == 0], "Non-progressors")

# ── Figure 1: Spaghetti plot of CEA trajectories ────────────────────────────
# FIX: Restrict x-axis to [-2, 10] years for readability
set.seed(42)
fig1_pids_prog <- sample(first_lm[outcome == 1]$PATIENT_ID, min(25, sum(first_lm$outcome == 1)))
fig1_pids_noprog <- sample(first_lm[outcome == 0]$PATIENT_ID, min(25, sum(first_lm$outcome == 0)))
fig1_pids <- c(fig1_pids_prog, fig1_pids_noprog)
fig1_data <- crc_cea[PATIENT_ID %in% fig1_pids]
fig1_data <- merge(fig1_data, first_lm[, .(PATIENT_ID, outcome)], by = "PATIENT_ID")
fig1_data[, prog_label := ifelse(outcome == 1, "Progressor", "Non-progressor")]

p1 <- ggplot(fig1_data, aes(x = time_from_dx / 365.25, y = log_cea,
                             group = PATIENT_ID, color = prog_label)) +
  geom_line(alpha = 0.5) +
  scale_color_manual(values = c("Progressor" = "#D62728", "Non-progressor" = "#2CA02C")) +
  coord_cartesian(xlim = c(-2, 10)) +
  labs(x = "Time from diagnosis (years)", y = "log(CEA)",
       color = "", title = "Figure 1: Individual CEA Trajectories") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
ggsave("fig1_cea_trajectories.pdf", p1, width = 8, height = 5)

# ── Figure 2: CEA slope distribution by outcome ─────────────────────────────
p2 <- ggplot(landmarks, aes(x = factor(outcome, labels = c("No progression", "Progression")),
                            y = cea_slope_90d)) +
  geom_boxplot(fill = c("#2CA02C", "#D62728"), alpha = 0.6, outlier.size = 0.3) +
  coord_cartesian(ylim = quantile(landmarks$cea_slope_90d, c(0.01, 0.99))) +
  labs(x = "", y = "CEA 90-day slope (log scale)",
       title = "Figure 2: CEA Slope by Progression Status") +
  theme_minimal(base_size = 12)
ggsave("fig2_cea_slope_boxplot.pdf", p2, width = 6, height = 5)

###############################################################################
## STEP 5: MODEL 0 / MODEL A — Primary Scientific Question
##         Does adding CEA trajectory improve prediction beyond clinical data?
###############################################################################

cat("\n========== MODEL 0 / MODEL A: Primary Comparison ==========\n")

# Model 0: clinical-only mixed-effects logistic regression
fml_0 <- outcome ~ time_since_dx_yr + tx_type + age + sex + stage + msi_status +
  prior_med + (1 | PATIENT_ID)
fml_0_fe <- ~ time_since_dx_yr + tx_type + age + sex + stage + msi_status + prior_med

# Model A: clinical + CEA mixed-effects logistic regression
fml_a <- outcome ~ log_cea_current + cea_slope_90d + cea_cv +
  time_since_dx_yr + tx_type + age + sex + stage + msi_status + prior_med +
  (1 | PATIENT_ID)
fml_a_fe <- ~ log_cea_current + cea_slope_90d + cea_cv +
  time_since_dx_yr + tx_type + age + sex + stage + msi_status + prior_med

vars_0 <- c("outcome", "time_since_dx_yr", "tx_type", "age", "sex", "stage",
            "msi_status", "prior_med")
vars_a <- c("outcome", "log_cea_current", "cea_slope_90d", "cea_cv",
            "time_since_dx_yr", "tx_type", "age", "sex", "stage",
            "msi_status", "prior_med")

# Fit Model 0 and Model A on the same clinical+CEA complete-case subset so the
# incremental value of CEA is compared fairly.
train_ca <- train[complete.cases(train[, ..vars_a])]
test_ca  <- test[complete.cases(test[, ..vars_a])]
cat("Train shared subset for Model 0 / A:", nrow(train_ca), "landmarks,",
    uniqueN(train_ca$PATIENT_ID), "patients\n")
cat("Test shared subset for Model 0 / A:", nrow(test_ca), "landmarks,",
    uniqueN(test_ca$PATIENT_ID), "patients\n")

model_0 <- glmer(
  fml_0,
  data    = train_ca,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)),
  nAGQ    = 1
)

model_a <- glmer(
  fml_a,
  data    = train_ca,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)),
  nAGQ    = 1
)

cat("\n── Model 0 Summary (clinical-only) ──\n")
print(summary(model_0))
cat("\n── Model A Summary (clinical + CEA) ──\n")
print(summary(model_a))

cat("\n── Model A: Odds Ratios ──\n")
model_a_or <- exp(fixef(model_a))
model_a_ci <- exp(confint(model_a, method = "Wald", parm = "beta_"))
print(cbind(OR = model_a_or, model_a_ci))

###############################################################################
## STEP 6: MODEL B — Genomic Extension
##         Clinical + CEA + genomic LASSO logistic regression
###############################################################################

cat("\n========== MODEL B: Genomic Extension ==========\n")

gene_terms <- paste(top_genes, collapse = " + ")
fml_b <- as.formula(paste0(
  "outcome ~ log_cea_current + cea_slope_90d + cea_cv + ",
  "time_since_dx_yr + tx_type + age + sex + stage + msi_status + prior_med + tmb + ",
  gene_terms
))

vars_b <- c("outcome", "log_cea_current", "cea_slope_90d", "cea_cv",
            "time_since_dx_yr", "tx_type", "age", "sex", "stage",
            "msi_status", "prior_med", "tmb", top_genes, "PATIENT_ID")
train_b <- train[complete.cases(train[, ..vars_b])]
train_b <- train_b[order(PATIENT_ID, landmark_time)]

cat("Train Model B (genomic complete cases):", nrow(train_b), "landmarks,",
    uniqueN(train_b$PATIENT_ID), "patients\n")

X_train_b <- model.matrix(fml_b, data = train_b)[, -1]
y_train_b <- train_b$outcome

fold_ids <- unique(train_b[, .(PATIENT_ID)])
set.seed(629)
fold_ids[, fold := sample(rep(1:5, length.out = .N))]
foldid_vec <- fold_ids$fold[match(train_b$PATIENT_ID, fold_ids$PATIENT_ID)]

cv_lasso <- cv.glmnet(X_train_b, y_train_b, family = "binomial",
                      alpha = 1, foldid = foldid_vec)

cat("\nLASSO lambda.min:", cv_lasso$lambda.min, "\n")
cat("LASSO lambda.1se:", cv_lasso$lambda.1se, "\n")

lasso_coefs <- coef(cv_lasso, s = "lambda.min")
nonzero <- lasso_coefs[lasso_coefs[, 1] != 0, , drop = FALSE]
cat("\n── Model B: Non-zero LASSO coefficients ──\n")
print(nonzero)

###############################################################################
## STEP 7: MODEL EVALUATION ON TEST SET
###############################################################################

cat("\n========== MODEL EVALUATION ON TEST SET ==========\n")

clamp_prob <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

# ── Model 0 / Model A predictions on shared clinical+CEA subset ─────────────
# Use only coefficients that were actually estimated (handle rank-deficiency)
beta_0 <- fixef(model_0)
beta_a <- fixef(model_a)
X_test_0 <- model.matrix(fml_0_fe, data = test_ca)
X_test_a <- model.matrix(fml_a_fe, data = test_ca)

# Align columns: keep only columns that match estimated coefficients
X_test_0 <- X_test_0[, names(beta_0), drop = FALSE]
X_test_a <- X_test_a[, names(beta_a), drop = FALSE]

test_ca[, pred_0 := plogis(as.numeric(X_test_0 %*% beta_0))]
test_ca[, pred_a := plogis(as.numeric(X_test_a %*% beta_a))]

test_1pp_ca <- test_ca[, .SD[which.min(landmark_time)], by = PATIENT_ID]

# ── Model B predictions on genomic complete-case subset ─────────────────────
test_b <- test[complete.cases(test[, ..vars_b])]
cat("Test Model B (genomic complete cases):", nrow(test_b), "landmarks,",
    uniqueN(test_b$PATIENT_ID), "patients\n")
X_test_b <- model.matrix(fml_b, data = test_b)[, -1]
test_b[, pred_b := as.numeric(predict(cv_lasso, newx = X_test_b, s = "lambda.min",
                                      type = "response"))]
test_1pp_b <- test_b[, .SD[which.min(landmark_time)], by = PATIENT_ID]

# Shared subset for Model A vs Model B paired comparison
shared_keys_ab <- intersect(test_ca$landmark_id, test_b$landmark_id)
test_shared_ab <- test_b[landmark_id %in% shared_keys_ab]
X_shared_a <- model.matrix(fml_a_fe, data = test_shared_ab)
X_shared_a <- X_shared_a[, names(beta_a), drop = FALSE]
test_shared_ab[, pred_a := plogis(as.numeric(X_shared_a %*% beta_a))]
test_1pp_shared_ab <- test_shared_ab[, .SD[which.min(landmark_time)], by = PATIENT_ID]

cat("\nTest (1 per patient): Model 0 / A =", nrow(test_1pp_ca),
    ", Model B =", nrow(test_1pp_b),
    ", Shared A/B =", nrow(test_1pp_shared_ab), "\n")

# ── AUC ──────────────────────────────────────────────────────────────────────
roc_0 <- roc(test_1pp_ca$outcome, test_1pp_ca$pred_0, quiet = TRUE)
roc_a <- roc(test_1pp_ca$outcome, test_1pp_ca$pred_a, quiet = TRUE)
roc_b <- roc(test_1pp_b$outcome, test_1pp_b$pred_b, quiet = TRUE)
roc_a_s <- roc(test_1pp_shared_ab$outcome, test_1pp_shared_ab$pred_a, quiet = TRUE)
roc_b_s <- roc(test_1pp_shared_ab$outcome, test_1pp_shared_ab$pred_b, quiet = TRUE)

cat("\n── Discrimination (AUC) ──\n")
cat("Model 0 (clinical-only GLMM): AUC =", round(auc(roc_0), 3),
    sprintf("(95%% CI: %.3f-%.3f)", ci.auc(roc_0)[1], ci.auc(roc_0)[3]), "\n")
cat("Model A (clinical + CEA GLMM): AUC =", round(auc(roc_a), 3),
    sprintf("(95%% CI: %.3f-%.3f)", ci.auc(roc_a)[1], ci.auc(roc_a)[3]), "\n")
cat("Model B (genomic LASSO): AUC =", round(auc(roc_b), 3),
    sprintf("(95%% CI: %.3f-%.3f)", ci.auc(roc_b)[1], ci.auc(roc_b)[3]), "\n")

delong_0a <- roc.test(roc_0, roc_a)
delong_ab <- roc.test(roc_a_s, roc_b_s)
cat("Primary DeLong p-value (Model 0 vs A):", format.pval(delong_0a$p.value, digits = 3), "\n")
cat("Secondary DeLong p-value (Model A vs B, shared subset):",
    format.pval(delong_ab$p.value, digits = 3), "\n")

# ── Brier Score ──────────────────────────────────────────────────────────────
brier_0 <- mean((test_1pp_ca$outcome - test_1pp_ca$pred_0)^2)
brier_a <- mean((test_1pp_ca$outcome - test_1pp_ca$pred_a)^2)
brier_b <- mean((test_1pp_b$outcome - test_1pp_b$pred_b)^2)
cat("\n── Brier Score ──\n")
cat("Model 0:", round(brier_0, 4), "\n")
cat("Model A:", round(brier_a, 4), "\n")
cat("Model B:", round(brier_b, 4), "\n")

# ── Calibration on logit scale ──────────────────────────────────────────────
test_1pp_ca[, pred_0_clamped := clamp_prob(pred_0)]
test_1pp_ca[, pred_a_clamped := clamp_prob(pred_a)]
test_1pp_b[,  pred_b_clamped := clamp_prob(pred_b)]

cal_fit_0 <- glm(outcome ~ I(qlogis(pred_0_clamped)), data = test_1pp_ca, family = binomial)
cal_fit_a <- glm(outcome ~ I(qlogis(pred_a_clamped)), data = test_1pp_ca, family = binomial)
cal_fit_b <- glm(outcome ~ I(qlogis(pred_b_clamped)), data = test_1pp_b, family = binomial)

cat("\n── Calibration (logit scale) ──\n")
cat("Model 0: intercept =", round(coef(cal_fit_0)[1], 3),
    ", slope =", round(coef(cal_fit_0)[2], 3), "\n")
cat("Model A: intercept =", round(coef(cal_fit_a)[1], 3),
    ", slope =", round(coef(cal_fit_a)[2], 3), "\n")
cat("Model B: intercept =", round(coef(cal_fit_b)[1], 3),
    ", slope =", round(coef(cal_fit_b)[2], 3), "\n")

# ── Figure 4: Calibration plot (Model A) ────────────────────────────────────
test_1pp_ca[, pred_decile := cut(pred_a, breaks = quantile(pred_a, probs = seq(0, 1, 0.1)),
                                 include.lowest = TRUE, labels = 1:10)]
cal_data <- test_1pp_ca[, .(observed = mean(outcome),
                            predicted = mean(pred_a),
                            n = .N), by = pred_decile]

p4 <- ggplot(cal_data, aes(x = predicted, y = observed)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = FALSE, color = "#1F77B4") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability", y = "Observed proportion",
       title = "Figure 4: Calibration Plot (Model A)") +
  theme_minimal(base_size = 12)
ggsave("fig4_calibration.pdf", p4, width = 6, height = 5)

# ── Figure 3: Primary ROC comparison ────────────────────────────────────────
pdf("fig3_roc_curves.pdf", width = 6, height = 5)
plot(roc_0, col = "#7F7F7F", lwd = 2,
     main = "Figure 3: ROC Curves for Primary Scientific Question")
plot(roc_a, col = "#1F77B4", lwd = 2, add = TRUE)
legend("bottomright",
       legend = c(sprintf("Model 0 (clinical-only): AUC = %.3f", auc(roc_0)),
                  sprintf("Model A (clinical + CEA): AUC = %.3f", auc(roc_a))),
       col = c("#7F7F7F", "#1F77B4"), lwd = 2, bty = "n")
dev.off()

# ── Summary Table 4 ─────────────────────────────────────────────────────────
cat("\n========== TABLE 4: Model Comparison ==========\n")
cat(sprintf("%-22s %-16s %-20s %-10s %-10s\n",
            "Model", "Eval subset", "AUC (95% CI)", "Brier", "Cal.Slope"))
cat(sprintf("%-22s %-16s %-20s %-10s %-10s\n",
            "0: GLMM clin-only",
            "clin+CEA CC",
            sprintf("%.3f (%.3f-%.3f)", auc(roc_0), ci.auc(roc_0)[1], ci.auc(roc_0)[3]),
            sprintf("%.4f", brier_0),
            sprintf("%.3f", coef(cal_fit_0)[2])))
cat(sprintf("%-22s %-16s %-20s %-10s %-10s\n",
            "A: GLMM clin+CEA",
            "clin+CEA CC",
            sprintf("%.3f (%.3f-%.3f)", auc(roc_a), ci.auc(roc_a)[1], ci.auc(roc_a)[3]),
            sprintf("%.4f", brier_a),
            sprintf("%.3f", coef(cal_fit_a)[2])))
cat(sprintf("%-22s %-16s %-20s %-10s %-10s\n",
            "B: LASSO + genomic",
            "genomic CC",
            sprintf("%.3f (%.3f-%.3f)", auc(roc_b), ci.auc(roc_b)[1], ci.auc(roc_b)[3]),
            sprintf("%.4f", brier_b),
            sprintf("%.3f", coef(cal_fit_b)[2])))

cat("\nPrimary scientific question: compare Model 0 vs Model A on the same subset.\n")
cat("Genomic extension: compare Model A vs Model B only on the shared genomic subset.\n")

# ── High-precision audit: are Model A and Model B AUCs literally identical? ─
cat("\n========== HIGH-PRECISION AUC AUDIT (Model A vs B) ==========\n")
cat(sprintf("Model A AUC on test_1pp_ca           = %.6f (n = %d)\n",
            as.numeric(auc(roc_a)), nrow(test_1pp_ca)))
cat(sprintf("Model B AUC on test_1pp_b            = %.6f (n = %d)\n",
            as.numeric(auc(roc_b)), nrow(test_1pp_b)))
cat(sprintf("Model A AUC on shared A/B subset     = %.6f (n = %d)\n",
            as.numeric(auc(roc_a_s)), nrow(test_1pp_shared_ab)))
cat(sprintf("Model B AUC on shared A/B subset     = %.6f (n = %d)\n",
            as.numeric(auc(roc_b_s)), nrow(test_1pp_shared_ab)))

# Correlation of predictions on the shared subset
b_preds_1pp <- test_b[, .SD[which.min(landmark_time)], by = PATIENT_ID][,
                       .(landmark_id, pred_b)]
audit_m <- merge(test_1pp_shared_ab[, .(landmark_id, pred_a)], b_preds_1pp,
                 by = "landmark_id")
cat(sprintf("Pearson r  (pred_a, pred_b) = %.4f\n",
            cor(audit_m$pred_a, audit_m$pred_b)))
cat(sprintf("Spearman r (pred_a, pred_b) = %.4f  [AUC is rank-based]\n",
            cor(audit_m$pred_a, audit_m$pred_b, method = "spearman")))
# Audit interpretation: if Spearman r ~ 1, the two models induce nearly the
# same patient ordering and identical AUC is mathematically expected; the
# genomic features have not moved any progressor / non-progressor pair across
# the decision threshold of any other pair. The DeLong p-value (0.977) then
# reflects the size of the within-pair perturbation, not the absence of any
# effect.

# Preserve names used later by the conformal section.
train_a <- train_ca
test_a <- test_ca

###############################################################################
## STEP 8: CONFORMAL PREDICTION SETS
##         FIX: Use glm instead of glmer for speed
###############################################################################

cat("\n========== CONFORMAL PREDICTION SETS ==========\n")

# Split training set into proper-train (50% of full) and calibration (20% of full)
set.seed(629)
train_pids_shuffle <- sample(train_pids)
n_proper <- round(0.714 * length(train_pids))  # 50/70 of training
proper_pids <- train_pids_shuffle[1:n_proper]
cal_pids    <- train_pids_shuffle[(n_proper + 1):length(train_pids_shuffle)]

conf_vars <- vars_a  # clinical variables only

proper_train <- train[PATIENT_ID %in% proper_pids]
proper_train <- proper_train[complete.cases(proper_train[, ..conf_vars])]

# For calibration: use 1-per-patient (first landmark) to satisfy exchangeability
cal_set_all  <- train[PATIENT_ID %in% cal_pids]
cal_set_all  <- cal_set_all[complete.cases(cal_set_all[, ..conf_vars])]
cal_set      <- cal_set_all[, .SD[which.min(landmark_time)], by = PATIENT_ID]

cat("Proper train:", nrow(proper_train), "landmarks from",
    uniqueN(proper_train$PATIENT_ID), "patients\n")
cat("Calibration (1-per-patient):", nrow(cal_set), "patients\n")

# Use glm (not glmer) for conformal — faster, and exchangeability
# assumption is cleaner without random effects
fml_conf <- outcome ~ log_cea_current + cea_slope_90d + cea_cv +
  time_since_dx_yr + tx_type + age + sex + stage + msi_status + prior_med

model_conf <- glm(fml_conf, data = proper_train, family = binomial)

cat("\n── Conformal base model (GLM) ──\n")
cat("AIC:", AIC(model_conf), "\n")

# Nonconformity scores on calibration set (1-per-patient)
cal_set[, p_hat := predict(model_conf, newdata = cal_set, type = "response")]
# Score: 1 - p_hat(true class)
cal_set[, score := ifelse(outcome == 1, 1 - p_hat, p_hat)]

# Conformal quantile at alpha = 0.10 (90% coverage)
alpha <- 0.10
n_cal <- nrow(cal_set)
q_level <- ceiling((1 - alpha) * (n_cal + 1)) / n_cal
q_hat <- quantile(cal_set$score, probs = min(q_level, 1))
cat("\nConformal quantile (alpha=0.10):", round(q_hat, 4), "\n")

# Build prediction sets on test data (also 1-per-patient for fair evaluation)
test_conf <- copy(test_a)
test_conf[, p_hat_conf := predict(model_conf, newdata = test_conf, type = "response")]

# For each test point, include label y if 1 - p_hat(y|x) <= q_hat
# score for y=0 is p_hat; score for y=1 is 1-p_hat
test_conf[, include_0 := p_hat_conf <= q_hat]
test_conf[, include_1 := (1 - p_hat_conf) <= q_hat]

test_conf[, pred_set_size := as.integer(include_0) + as.integer(include_1)]

# Coverage: true label is in the prediction set
test_conf[, covered := ifelse(outcome == 1, include_1, include_0)]

# Use 1-per-patient for coverage reporting
test_conf_1pp <- test_conf[, .SD[which.min(landmark_time)], by = PATIENT_ID]

cat("\n── Conformal Results (1 per patient) ──\n")
cat("Empirical coverage:", round(mean(test_conf_1pp$covered), 3),
    "(target >= ", 1 - alpha, ")\n")
cat("Singleton sets (size=1):", round(mean(test_conf_1pp$pred_set_size == 1), 3), "\n")
cat("Empty sets (size=0):", round(mean(test_conf_1pp$pred_set_size == 0), 3), "\n")
cat("Ambiguous sets (size=2):", round(mean(test_conf_1pp$pred_set_size == 2), 3), "\n")

# Singleton fraction by outcome
cat("\nSingleton fraction among progressors:",
    round(mean(test_conf_1pp[outcome == 1]$pred_set_size == 1), 3), "\n")
cat("Singleton fraction among non-progressors:",
    round(mean(test_conf_1pp[outcome == 0]$pred_set_size == 1), 3), "\n")

# ── Figure 5: Singleton fraction by predicted risk decile ────────────────────
test_conf_1pp[, risk_decile := cut(p_hat_conf,
                                    breaks = quantile(p_hat_conf, probs = seq(0, 1, 0.1)),
                                    include.lowest = TRUE, labels = 1:10)]
conf_plot_data <- test_conf_1pp[, .(singleton_frac = mean(pred_set_size == 1),
                                     mean_risk = mean(p_hat_conf)),
                                 by = risk_decile]

p5 <- ggplot(conf_plot_data, aes(x = as.numeric(risk_decile), y = singleton_frac)) +
  geom_bar(stat = "identity", fill = "#1F77B4", alpha = 0.7) +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  labs(x = "Predicted risk decile", y = "Fraction of singleton prediction sets",
       title = "Figure 5: Conformal Informativeness by Risk Decile") +
  theme_minimal(base_size = 12)
ggsave("fig5_conformal_singletons.pdf", p5, width = 7, height = 5)

###############################################################################
## DONE
###############################################################################

cat("\n\n========== ANALYSIS COMPLETE ==========\n")
cat("Output files:\n")
cat("  fig1_cea_trajectories.pdf\n")
cat("  fig2_cea_slope_boxplot.pdf\n")
cat("  fig3_roc_curves.pdf\n")
cat("  fig4_calibration.pdf\n")
cat("  fig5_conformal_singletons.pdf\n")

# ── Reproducibility footprint: capture R version + package versions ──────────
# Writes a sessionInfo.txt that records exactly which package versions were
# loaded for this run. Anyone re-running the script can compare against this
# file to detect environment drift.
capture.output(sessionInfo(), file = "sessionInfo.txt")
cat("  sessionInfo.txt\n")
