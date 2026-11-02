################################################################################
#  This script extracts the actual parameters for
#  calibrating the MC simulation from real data
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(tidyr)
  library(here)
})
source("R/00_visit_time.R")

cat("===============================================\n")
cat("  PARAMETER ESTIMATION FOR MC SIMULATION\n")
cat("===============================================\n\n")


CLIN <- paste0(here(), "/raw/als-metadata/answerals/clinical/")

df <- read.csv(paste0(here(), "/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"),
               check.names = FALSE)
mapper <- read.csv(paste0(here(), "/raw/als-metadata/answerals/aals_dataportal_datatable.csv"),
                   check.names = FALSE)
alsfrs <- read.csv(paste0(CLIN, "ALSFRS_R.csv"))
annotation <- read.csv(paste0(here(), "/raw/somascan/documents/annotation.csv"))
dem <- read.csv(paste0(CLIN, "Demographics.csv"), check.names = FALSE)


misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
seq_cols <- annotation %>%
  filter(!is.na(UniProt.ID), UniProt.ID != "") %>%
  pull(SeqId) %>% as.character() %>%
  intersect(setdiff(colnames(df), misc_columns))
df <- df[, c(misc_columns, seq_cols)]
cat("SomaScan aptamers with a UniProt accession:", length(seq_cols), "\n\n")


meta_cols <- c("Participant_ID", "GUID", "Number_of_Visits", "SEX",
               "SUBJECT_GROUP", "Site_of_Onset", "AGE_AT_SYMPTOM_ONSET")

# ===============================================
# 1. SAMPLE SIZE & VISIT STRUCTURE
# ===============================================

cat("   1. SAMPLE SIZE & VISIT STRUCTURE\n\n")

# Missing totals and dates are removed before counting, so a patient's visit count
# is the number of usable assessments.
alsfrs_proc <- alsfrs %>%
  select(SubjectUID, Visit_Date, alsfrst) %>%
  rename(GUID = SubjectUID) %>%
  filter(!is.na(alsfrst), !is.na(Visit_Date)) %>%
  group_by(GUID) %>%
  filter(n() >= 2) %>%
  arrange(Visit_Date, .by_group = TRUE) %>%
  mutate(time_years = Visit_Date / 365.25) %>%
  ungroup()


als_guids <- mapper %>%
  filter(SUBJECT_GROUP == "ALS") %>%
  pull(GUID)

alsfrs_als <- alsfrs_proc %>%
  filter(GUID %in% als_guids)

visit_counts <- alsfrs_als %>%
  group_by(GUID) %>%
  summarise(n_visits = n(), .groups = "drop")

N_patients <- nrow(visit_counts)

cat("N patients (ALS, >=2 ALSFRS visits):", N_patients, "\n")
cat("Visit count distribution:\n")
print(table(visit_counts$n_visits))
cat("Mean visits per patient:", round(mean(visit_counts$n_visits), 2), "\n")
cat("Median visits per patient:", median(visit_counts$n_visits), "\n\n")

# Inter-visit gaps
visit_gaps <- alsfrs_als %>%
  group_by(GUID) %>%
  arrange(time_years, .by_group = TRUE) %>%
  mutate(gap = time_years - lag(time_years)) %>%
  filter(!is.na(gap))

cat("Inter-visit gap (years):\n")
cat("  Mean:", round(mean(visit_gaps$gap), 3), "\n")
cat("  SD:  ", round(sd(visit_gaps$gap), 3), "\n")
cat("  Min: ", round(min(visit_gaps$gap), 3), "\n")
cat("  Max: ", round(max(visit_gaps$gap), 3), "\n")
cat("  Median:", round(median(visit_gaps$gap), 3), "\n\n")

# Follow-up duration
followup <- alsfrs_als %>%
  group_by(GUID) %>%
  summarise(
    fu_years = max(time_years) - min(time_years),
    .groups = "drop"
  )

cat("Total follow-up per patient (years):\n")
cat("  Mean:", round(mean(followup$fu_years, na.rm = TRUE), 3), "\n")
cat("  SD:  ", round(sd(followup$fu_years, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(followup$fu_years, na.rm = TRUE), 3), "\n\n")


# ===============================================
# 2. ALSFRS-R MIXED MODEL VARIANCE COMPONENTS
# ===============================================

cat("--- 2. ALSFRS-R MIXED MODEL---\n\n")

# Baseline ALSFRS
baseline_als <- alsfrs_als %>%
  group_by(GUID) %>%
  filter(time_years == min(time_years)) %>%
  slice(1) %>%
  ungroup()

cat("Baseline ALSFRS-R:\n")
cat("  Mean:", round(mean(baseline_als$alsfrst, na.rm = TRUE), 2), "\n")
cat("  SD:  ", round(sd(baseline_als$alsfrst, na.rm = TRUE), 2), "\n")
cat("  Median:", median(baseline_als$alsfrst, na.rm = TRUE), "\n")
cat("  Range:", range(baseline_als$alsfrst, na.rm = TRUE), "\n\n")


m_alsfrs <- lmer(
  alsfrst ~ time_years + (time_years | GUID),
  data = alsfrs_als,
  REML = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

cat("Model summary:\n")
print(summary(m_alsfrs))

vc <- VarCorr(m_alsfrs)
vc_df <- as.data.frame(vc)

cat("\n\nVariance components (for simulation):\n")
cat("─────────────────────────────────────\n")

intercept_var <- vc_df$vcov[which(vc_df$var1 == "(Intercept)" & is.na(vc_df$var2))]
slope_var <- vc_df$vcov[which(vc_df$var1 == "time_years"  & is.na(vc_df$var2))]
cov_int_slope <- vc_df$vcov[which(vc_df$var1 == "(Intercept)" & vc_df$var2 == "time_years")]
residual_var <- vc_df$vcov[which(vc_df$grp == "Residual")]

cat("  Random intercept SD:", round(sqrt(intercept_var), 4), "\n")
cat("  Random slope SD:    ", round(sqrt(slope_var), 4), "\n")
cat("  Intercept-slope cor:", round(cov_int_slope / (sqrt(intercept_var) * sqrt(slope_var)), 4), "\n")
cat("  Residual SD:        ", round(sqrt(residual_var), 4), "\n")

# Fixed effects
fe <- fixef(m_alsfrs)
fe_alsfrs <- fe

cat("\nFixed effects:\n")
cat("  Intercept (mean baseline):", round(fe["(Intercept)"], 4), "\n")
cat("  Slope (mean decline/year):", round(fe["time_years"], 4), "\n\n")


# ===============================================
# 3. ALSFRS-R OLS SLOPE DISTRIBUTION
# ===============================================

cat("--- 3. EMPIRICAL OLS SLOPES ---\n\n")

ols_slopes <- alsfrs_als %>%
  filter(!is.na(alsfrst), !is.na(time_years)) %>%
  group_by(GUID) %>%
  filter(n() >= 2, var(time_years) > 0) %>%
  summarise(
    ols_slope = coef(lm(alsfrst ~ time_years))[2],
    ols_intercept = coef(lm(alsfrst ~ time_years))[1],
    n_visits = n(),
    .groups = "drop"
  ) %>%
  filter(is.finite(ols_slope))

n_dropped <- nrow(visit_counts) - nrow(ols_slopes)
if (n_dropped > 0)
  cat("  NOTE:", n_dropped, "of", nrow(visit_counts),
      "patients dropped (no time spread)\n")

cat("OLS slopes (pts/year):\n")
cat("  Mean:  ", round(mean(ols_slopes$ols_slope), 4), "\n")
cat("  SD:    ", round(sd(ols_slopes$ols_slope), 4), "\n")
cat("  Median:", round(median(ols_slopes$ols_slope), 4), "\n")
cat("  Q25:   ", round(quantile(ols_slopes$ols_slope, 0.25), 4), "\n")
cat("  Q75:   ", round(quantile(ols_slopes$ols_slope, 0.75), 4), "\n")
cat("  Min:   ", round(min(ols_slopes$ols_slope), 4), "\n")
cat("  Max:   ", round(max(ols_slopes$ols_slope), 4), "\n\n")

# BLUP slopes for comparison
blup_slopes <- ranef(m_alsfrs)$GUID %>%
  tibble::rownames_to_column("GUID") %>%
  mutate(blup_slope = time_years + fe_alsfrs["time_years"]) %>%
  select(GUID, blup_slope)

# Compare OLS vs BLUP
merged_slopes <- inner_join(ols_slopes, blup_slopes, by = "GUID")

cat("OLS vs BLUP slope comparison (sanity check):\n")
cat("  Correlation:", round(cor(merged_slopes$ols_slope, merged_slopes$blup_slope,
                                use = "complete.obs"), 4), "\n")
cat("  OLS SD:    ", round(sd(merged_slopes$ols_slope, na.rm = TRUE), 4), "\n")
cat("  BLUP SD:   ", round(sd(merged_slopes$blup_slope, na.rm = TRUE), 4), "\n")
cat("  Shrinkage ratio (BLUP_SD / OLS_SD):",
    round(sd(merged_slopes$blup_slope, na.rm = TRUE) /
          sd(merged_slopes$ols_slope, na.rm = TRUE), 4), "\n")
cat("  Cohort-level OLS reliability, slope variance / var(OLS slope):",
    round(slope_var / var(ols_slopes$ols_slope), 3), "\n\n")


# ===============================================
# 4. SOMASCAN PROTEIN TRAJECTORY PARAMETERS
# ===============================================

cat("--- 4. SOMASCAN PROTEIN TRAJECTORIES ---\n\n")

all_prot_cols <- seq_cols

vt <- soma_visit_time(paste0(CLIN, "Plasma_Sample.csv"), paste0(CLIN, "ALSFRS_R.csv"))

prot_with_time <- df %>%
  filter(SubjectUID %in% als_guids) %>%
  attach_soma_time(vt) %>%
  filter(GUID %in% visit_counts$GUID) %>%
  select(GUID, TimePoint, time_years, all_of(all_prot_cols))

cat("  Time source for these draws:\n")
print(table((df %>% filter(SubjectUID %in% als_guids) %>%
               attach_soma_time(vt, verbose = FALSE) %>%
               filter(GUID %in% visit_counts$GUID))$time_source))

# Count how many patients have ≥2 protein timepoints
prot_visit_counts <- prot_with_time %>%
  group_by(GUID) %>%
  summarise(n_prot_visits = n(), .groups = "drop")

cat("Patients with >=2 SomaScan timepoints:", sum(prot_visit_counts$n_prot_visits >= 2), "\n")
cat("Patients with >=3 SomaScan timepoints:", sum(prot_visit_counts$n_prot_visits >= 3), "\n\n")

# ---- SomaScan timepoint distribution ---------
# this drives protein-slope precision.
cat("SomaScan timepoint distribution (drives protein slopes):\n")
print(table(prot_visit_counts$n_prot_visits))
cat("Mean SomaScan timepoints per patient:", round(mean(prot_visit_counts$n_prot_visits), 2), "\n")
cat("Median:", median(prot_visit_counts$n_prot_visits), "\n\n")

# Restricted to the >=2-draw patients the two-stage analysis actually uses
pvc2 <- prot_visit_counts %>% filter(n_prot_visits >= 2)
cat("SomaScan timepoint distribution, >=2-draw patients only (the sim's N):\n")
print(table(pvc2$n_prot_visits))
cat("  N =", nrow(pvc2), "  Mean:", round(mean(pvc2$n_prot_visits), 2),
    "  Median:", median(pvc2$n_prot_visits), "\n\n")

# Inter-draw gaps and follow-up span for SomaScan (protein slope precision ~ visit spread)
prot_gaps <- prot_with_time %>%
  group_by(GUID) %>% filter(n() >= 2) %>%
  arrange(time_years, .by_group = TRUE) %>%
  mutate(gap = time_years - lag(time_years)) %>%
  filter(!is.na(gap)) %>% ungroup()
cat("SomaScan inter-draw gap (years):\n")
cat("  Mean:", round(mean(prot_gaps$gap), 3),
    "  SD:", round(sd(prot_gaps$gap), 3),
    "  Median:", round(median(prot_gaps$gap), 3), "\n")

prot_span <- prot_with_time %>%
  group_by(GUID) %>% filter(n() >= 2) %>%
  summarise(span = max(time_years) - min(time_years), .groups = "drop")
cat("SomaScan follow-up span (years): Mean:", round(mean(prot_span$span), 3),
    " Median:", round(median(prot_span$span), 3), "\n\n")

# How many ALSFRS visits do the SAME patients have?
visits_both <- visit_counts %>%
  inner_join(prot_visit_counts, by = "GUID") %>%
  filter(n_prot_visits >= 2)
cat("For the >=2-draw patients: ALSFRS visits vs SomaScan draws\n")
cat("  Mean ALSFRS visits:", round(mean(visits_both$n_visits), 2),
    "  Mean SomaScan draws:", round(mean(visits_both$n_prot_visits), 2), "\n")
cat("  Correlation:", round(cor(visits_both$n_visits, visits_both$n_prot_visits), 3), "\n")

cat("\nALSFRS visit distribution for the >=2-draw patients (the sim's N):\n")
print(table(visits_both$n_visits))
cat("  <=3 ALSFRS visits:", sum(visits_both$n_visits <= 3),
    "of", nrow(visits_both),
    sprintf("(%.1f%%)", 100 * mean(visits_both$n_visits <= 3)), "\n")
cat("  vs full ALSFRS cohort <=3 visits:", sum(visit_counts$n_visits <= 3),
    "of", nrow(visit_counts),
    sprintf("(%.1f%%)", 100 * mean(visit_counts$n_visits <= 3)), "\n")

# ALSFRS follow-up span for the same patients (spread, not just count, drives it)
fu_both <- followup %>% filter(GUID %in% visits_both$GUID)
cat("  ALSFRS follow-up span, >=2-draw patients: mean",
    round(mean(fu_both$fu_years, na.rm = TRUE), 3),
    " median", round(median(fu_both$fu_years, na.rm = TRUE), 3), "\n\n")

# Restrict to patients with ≥2 protein visits
guids_longitudinal <- prot_visit_counts %>%
  filter(n_prot_visits >= 2) %>%
  pull(GUID)

prot_long_filt <- prot_with_time %>%
  filter(GUID %in% guids_longitudinal)

# Sample a set of proteins for variance component estimation
# (fitting lmer to all 7000 would take ages to converge).
set.seed(42)
n_sample_prots <- 100
sampled_prots <- sample(all_prot_cols, min(n_sample_prots, length(all_prot_cols)))

cat("Estimating variance components for", length(sampled_prots),
    "randomly sampled proteins...\n\n")

prot_vc_list <- list()
prot_ols_sd_list <- numeric()

for (i in seq_along(sampled_prots)) {
  p <- sampled_prots[i]

  d <- prot_long_filt %>%
    select(GUID, time_years, all_of(p)) %>%
    rename(y = all_of(p)) %>%
    mutate(y = log2(y)) %>%
    filter(is.finite(y))

  # Mixed model approach
  m <- tryCatch(
    lmer(y ~ time_years + (time_years | GUID), data = d, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL,
    warning = function(w) {
      suppressWarnings(
        lmer(y ~ time_years + (time_years | GUID), data = d, REML = TRUE,
             control = lmerControl(optimizer = "bobyqa"))
      )
    }
  )

  if (!is.null(m)) {
    vc <- as.data.frame(VarCorr(m))
    fe <- fixef(m)

    int_var <- vc$vcov[which(vc$var1 == "(Intercept)" & is.na(vc$var2))]
    slp_var <- vc$vcov[which(vc$var1 == "time_years"  & is.na(vc$var2))]
    res_var <- vc$vcov[which(vc$grp == "Residual")]

    if (length(int_var) == 1L && length(slp_var) == 1L && length(res_var) == 1L) {
      prot_vc_list[[p]] <- data.frame(
        protein = p,
        fixed_int = unname(fe["(Intercept)"]),
        fixed_slope = unname(fe["time_years"]),
        intercept_sd = sqrt(int_var),
        slope_sd = sqrt(slp_var),
        residual_sd = sqrt(res_var),
        stringsAsFactors = FALSE
      )
    }
  }

  # OLS slopes for this protein (for empirical effect size estimation)
  ols_p <- d %>%
    group_by(GUID) %>%
    filter(n() >= 2, var(time_years) > 0) %>%
    summarise(prot_slope = coef(lm(y ~ time_years))[2], .groups = "drop")

  prot_ols_sd_list[p] <- sd(ols_p$prot_slope, na.rm = TRUE)
}

prot_vc_df <- do.call(rbind, prot_vc_list)
rownames(prot_vc_df) <- NULL

# Drop any protein whose fit produced a non-finite component.
n_prot_fit <- nrow(prot_vc_df)
prot_vc_df <- prot_vc_df[is.finite(prot_vc_df$intercept_sd) &
                         is.finite(prot_vc_df$slope_sd) &
                         is.finite(prot_vc_df$residual_sd) &
                         is.finite(prot_vc_df$fixed_int) &
                         is.finite(prot_vc_df$fixed_slope), ]
if (nrow(prot_vc_df) < n_prot_fit)
  cat("  NOTE:", n_prot_fit - nrow(prot_vc_df), "of", n_prot_fit,
      "protein fits dropped (non-finite variance components)\n")
cat("  Proteins attempted:", length(sampled_prots),
    "  usable fits:", nrow(prot_vc_df), "\n")

cat("Protein variance components (across", nrow(prot_vc_df), "proteins):\n")
cat("=====================================\n")

cat("\nFixed intercept (log2 RFU):\n")
cat("  Mean:", round(mean(prot_vc_df$fixed_int, na.rm = TRUE), 4), "\n")
cat("  SD:  ", round(sd(prot_vc_df$fixed_int, na.rm = TRUE), 4), "\n")
cat("  Range:", round(range(prot_vc_df$fixed_int, na.rm = TRUE), 4), "\n")

cat("\nFixed slope (log2/year):\n")
cat("  Mean:", round(mean(prot_vc_df$fixed_slope, na.rm = TRUE), 4), "\n")
cat("  SD:  ", round(sd(prot_vc_df$fixed_slope, na.rm = TRUE), 4), "\n")

cat("\nRandom intercept SD (between-patient):\n")
cat("  Median:", round(median(prot_vc_df$intercept_sd, na.rm = TRUE), 4), "\n")
cat("  IQR:   ", round(quantile(prot_vc_df$intercept_sd, c(0.25, 0.75), na.rm = TRUE), 4), "\n")

cat("\nRandom slope SD (between-patient):\n")
cat("  Median:", round(median(prot_vc_df$slope_sd, na.rm = TRUE), 4), "\n")
cat("  IQR:   ", round(quantile(prot_vc_df$slope_sd, c(0.25, 0.75), na.rm = TRUE), 4), "\n")

cat("\nResidual SD (within-patient):\n")
cat("  Median:", round(median(prot_vc_df$residual_sd, na.rm = TRUE), 4), "\n")
cat("  IQR:   ", round(quantile(prot_vc_df$residual_sd, c(0.25, 0.75), na.rm = TRUE), 4), "\n")

cat("\nOLS slope SD across patients (per protein):\n")
cat("  Median:", round(median(prot_ols_sd_list, na.rm = TRUE), 4), "\n")
cat("  IQR:   ", round(quantile(prot_ols_sd_list, c(0.25, 0.75), na.rm = TRUE), 4), "\n\n")


# =======================================================
# 5. EMPIRICAL EFFECT SIZES: PROTEIN SLOPE ~ ALSFRS SLOPE
# =======================================================

cat("-------- 5. EMPIRICAL EFFECT SIZES ----------------------\n\n")

# Per-patient OLS slopes for the sampled proteins. OLS on purpose: beta here is
# compared with the OLS `est` that the simulation's effect sizes are set against.
prot_slopes_all <- prot_long_filt %>%
  filter(GUID %in% ols_slopes$GUID) %>%
  group_by(GUID) %>%
  filter(n() >= 2) %>%
  group_modify(~ {
    dat <- .x
    slopes <- sapply(sampled_prots, function(p) {
      y  <- log2(dat[[p]])
      ok <- is.finite(y) & is.finite(dat$time_years)
      if (sum(ok) >= 2 && var(dat$time_years[ok]) > 0) {
        unname(coef(lm(y[ok] ~ dat$time_years[ok]))[2])
      } else {
        NA_real_
      }
    })
    as.data.frame(t(slopes), check.names = FALSE)
  }) %>%
  ungroup()

missing_cols <- setdiff(sampled_prots, colnames(prot_slopes_all))
if (length(missing_cols) > 0)
  stop("prot_slopes_all is missing ", length(missing_cols),
       " protein columns (name mangling?). First: ", missing_cols[1])

# Merge with ALSFRS slopes
effect_df <- prot_slopes_all %>%
  inner_join(ols_slopes %>% select(GUID, ols_slope), by = "GUID")

# Compute beta (protein_slope ~ alsfrs_slope) for each protein
betas <- sapply(sampled_prots, function(p) {
  y <- effect_df[[p]]
  x <- effect_df$ols_slope
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) >= 10) {
    unname(coef(lm(y[ok] ~ x[ok]))[2])
  } else {
    NA_real_
  }
})

cat("Proteins with an estimable beta:", sum(!is.na(betas)), "of",
    length(sampled_prots), "  (patients per fit ~", nrow(effect_df), ")\n")

correlations <- sapply(sampled_prots, function(p) {
  y <- effect_df[[p]]
  x <- effect_df$ols_slope
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) >= 10) {
    unname(cor(y[ok], x[ok]))
  } else {
    NA_real_
  }
})

# ---- BASELINE (gamma) --------------------------------
# beta above is the slope and answers if protein rate tracks decline rate.
# gamma here is the baseline: does baseline protein level predict decline rate.
# It is what Dichotomise, CrudeRate_DFS and Baseline_Pred read.
prot_baseline <- prot_long_filt %>%
  filter(GUID %in% ols_slopes$GUID) %>%
  group_by(GUID) %>%
  arrange(time_years, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

gamma_df <- prot_baseline %>%
  inner_join(ols_slopes %>% select(GUID, ols_slope), by = "GUID")

gammas <- sapply(sampled_prots, function(p) {
  y <- log2(gamma_df[[p]])
  x <- as.numeric(scale(gamma_df$ols_slope))   # SD units, to match the generator
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) >= 10) unname(coef(lm(y[ok] ~ x[ok]))[2]) else NA_real_
})

gamma_cors <- sapply(sampled_prots, function(p) {
  y  <- log2(gamma_df[[p]])
  x  <- gamma_df$ols_slope
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) >= 10) unname(cor(y[ok], x[ok])) else NA_real_
})

# ols_slope is a noisy proxy for the true decline rate.
atten <- sqrt(slope_var) / sd(ols_slopes$ols_slope, na.rm = TRUE)

cat("\nBASELINE channel (correlation: baseline log2 protein vs ALSFRS slope):\n")
cat("  This is for gamma_true_values in functions.R.\n")
cat("  Estimable:", sum(!is.na(gamma_cors)), "of", length(sampled_prots),
    "  (patients:", nrow(gamma_df), ")\n")
cat("  Observed |r|:  median", round(median(abs(gamma_cors), na.rm = TRUE), 4),
    "  75th", round(quantile(abs(gamma_cors), 0.75, na.rm = TRUE), 4),
    "  max", round(max(abs(gamma_cors), na.rm = TRUE), 4), "\n")
cat("  attenuation factor:", round(atten, 3), "\n")
cat("  TRUE (deattenuated) |r|:  median",
    round(median(abs(gamma_cors), na.rm = TRUE) / atten, 4),
    "  75th", round(quantile(abs(gamma_cors), 0.75, na.rm = TRUE) / atten, 4),
    "  max", round(max(abs(gamma_cors), na.rm = TRUE) / atten, 4), "\n")
cat("  Proteins with observed |r| > 0.2:", sum(abs(gamma_cors) > 0.2, na.rm = TRUE),
    "of", sum(!is.na(gamma_cors)), "\n")

cat("\n  (absolute gamma, log2 per SD of slope\n")
cat("   median", round(median(gammas, na.rm = TRUE), 5),
    "  75th |gamma|", round(quantile(abs(gammas), 0.75, na.rm = TRUE), 5),
    "  max |gamma|", round(max(abs(gammas), na.rm = TRUE), 5), "\n")

cat("Effect sizes (beta: protein_slope ~ ALSFRS_slope):\n")
cat("  Mean:  ", round(mean(betas, na.rm = TRUE), 6), "\n")
cat("  SD:    ", round(sd(betas, na.rm = TRUE), 6), "\n")
cat("  Median:", round(median(betas, na.rm = TRUE), 6), "\n")
cat("  Q25:   ", round(quantile(betas, 0.25, na.rm = TRUE), 6), "\n")
cat("  Q75:   ", round(quantile(betas, 0.75, na.rm = TRUE), 6), "\n")
cat("  Range: ", round(range(betas, na.rm = TRUE), 6), "\n")

cat("\nCorrelations (protein_slope vs ALSFRS_slope):\n")
cat("  Mean:  ", round(mean(correlations, na.rm = TRUE), 4), "\n")
cat("  SD:    ", round(sd(correlations, na.rm = TRUE), 4), "\n")
cat("  Median:", round(median(correlations, na.rm = TRUE), 4), "\n")

# Top/bottom effect sizes (potential true hits)
beta_sorted <- sort(betas[!is.na(betas)])
cat("\nStrongest positive betas (top 5):\n")
print(round(tail(beta_sorted, 5), 6))
cat("Strongest negative betas (top 5):\n")
print(round(head(beta_sorted, 5), 6))

# What fraction might be "real" signals?
n_potential_hits <- sum(abs(correlations) > 0.2, na.rm = TRUE)
cat("\nProteins with |r| > 0.2:", n_potential_hits, "out of",
    sum(!is.na(correlations)), "\n")
cat("Proportion:", round(n_potential_hits / sum(!is.na(correlations)), 4), "\n\n")


# ========================================
# 6. INFORMATIVE DROPOUT CHECK
# ========================================

cat("---- 6. INFORMATIVE DROPOUT ----------------\n\n")

# ols_slopes carries n_visits (the visits that fed the slope fit). visit_counts carries
# a column of the same name, so the raw count is renamed to keep both.
dropout_check <- ols_slopes %>%
  inner_join(visit_counts %>% rename(n_visits_raw = n_visits), by = "GUID")

cat("Correlation between OLS slope and number of visits:\n")
ct <- cor.test(dropout_check$ols_slope, dropout_check$n_visits)
cat("  r =", round(ct$estimate, 4), " p =", signif(ct$p.value, 4), "\n")

ct_raw <- cor.test(dropout_check$ols_slope, dropout_check$n_visits_raw)
cat("  vs raw visit count: r =", round(ct_raw$estimate, 4),
    " p =", signif(ct_raw$p.value, 4), "\n")

dropout_check <- dropout_check %>%
  mutate(fast = ols_slope < median(ols_slope))

cat("\nMean visits — fast progressors:",
    round(mean(dropout_check$n_visits[dropout_check$fast]), 2), "\n")
cat("Mean visits — slow progressors:",
    round(mean(dropout_check$n_visits[!dropout_check$fast]), 2), "\n")

wt <- wilcox.test(n_visits ~ fast, data = dropout_check)
cat("Wilcoxon p-value:", signif(wt$p.value, 4), "\n")
sub_check <- dropout_check %>% filter(GUID %in% visits_both$GUID)
ct_sub <- cor.test(sub_check$ols_slope, sub_check$n_visits)
cat("  >=2-draw patients (n =", nrow(sub_check), "): r =", round(ct_sub$estimate, 4),
    " p =", signif(ct_sub$p.value, 4), "| mean ALSFRS-R visits =",
    round(mean(visits_both$n_visits), 2), "\n\n")


# =====================================
# 6b. SYMPTOM ONSET AND DELTA FS
# =====================================

cat("----- 6b. SYMPTOM ONSET AND DELTA FS -------------\n\n")

onset_df <- mapper %>%
  filter(SUBJECT_GROUP == "ALS") %>%
  transmute(GUID, Participant_ID,
            age_onset = suppressWarnings(as.numeric(AGE_AT_SYMPTOM_ONSET))) %>%
  left_join(dem %>%
              transmute(Participant_ID, age_screen = suppressWarnings(as.numeric(age))) %>%
              distinct(Participant_ID, .keep_all = TRUE),
            by = "Participant_ID") %>%
  inner_join(alsfrs_als %>% group_by(GUID) %>%
               slice_min(Visit_Date, n = 1, with_ties = FALSE) %>% ungroup() %>%
               select(GUID, alsfrs_bl = alsfrst, bl_day = Visit_Date),
             by = "GUID") %>%
  mutate(duration_years = pmax(0.25, (age_screen - age_onset) + bl_day / 365.25),
         delta_fs = (48 - alsfrs_bl) / (duration_years * 12)) %>%
  filter(is.finite(duration_years), age_screen >= age_onset)

dur_sub <- onset_df %>% filter(GUID %in% visits_both$GUID)
dfs_sub <- dur_sub %>%
  filter(is.finite(delta_fs), delta_fs > 0) %>%
  inner_join(ols_slopes %>% select(GUID, ols_slope), by = "GUID")
ct_dfs <- suppressWarnings(cor.test(dfs_sub$delta_fs, dfs_sub$ols_slope, method = "spearman"))

cat("Time from onset to first ALSFRS-R, >=2-draw patients (n =", nrow(dur_sub), "):\n")
cat("  Median", round(median(dur_sub$duration_years), 2), "years | IQR",
    round(quantile(dur_sub$duration_years, 0.25), 2), "-",
    round(quantile(dur_sub$duration_years, 0.75), 2), "\n")
cat("Spearman r(Delta FS, OLS ALSFRS-R slope), >=2-draw patients (n =", nrow(dfs_sub), "): r =",
    round(ct_dfs$estimate, 3), " p =", signif(ct_dfs$p.value, 3), "\n\n")


# =======================================
# 7. SUMMARY: SIMULATION PARAMETERS
# =======================================

cat("==========================================\n")
cat("  RECOMMENDED SIMULATION PARAMETERS\n")
cat("==========================================\n\n")

params <- data.frame(
  Parameter = c(
    "N (patients with >=2 draws)",
    "Visits per patient (median)",
    "Visits per patient (range)",
    "Inter-visit gap (mean years)",
    "Inter-visit gap (SD years)",
    "Follow-up duration (mean years)",
    "",
    "ALSFRS intercept (mean)",
    "ALSFRS slope (mean pts/year)",
    "ALSFRS intercept SD",
    "ALSFRS slope SD",
    "ALSFRS intercept-slope cor",
    "ALSFRS residual SD",
    "",
    "Protein intercept (mean log2)",
    "Protein intercept SD (median)",
    "Protein slope SD (median)",
    "Protein residual SD (median)",
    "",
    "Effect size beta (median)",
    "Effect size beta (75th pctl)",
    "Effect size beta (strongest)",
    "Proportion |r| > 0.2"
  ),
  Value = c(
    nrow(visits_both),
    median(visits_both$n_visits),
    paste(range(visits_both$n_visits), collapse = "-"),
    round(mean(visit_gaps$gap, na.rm = TRUE), 3),
    round(sd(visit_gaps$gap, na.rm = TRUE), 3),
    round(mean(fu_both$fu_years, na.rm = TRUE), 3),
    "",
    round(fe_alsfrs["(Intercept)"], 2),
    round(fe_alsfrs["time_years"], 2),
    round(sqrt(intercept_var), 2),
    round(sqrt(slope_var), 2),
    round(cov_int_slope / (sqrt(intercept_var) * sqrt(slope_var)), 3),
    round(sqrt(residual_var), 2),
    "",
    round(mean(prot_vc_df$fixed_int, na.rm = TRUE), 2),
    round(median(prot_vc_df$intercept_sd, na.rm = TRUE), 3),
    round(median(prot_vc_df$slope_sd, na.rm = TRUE), 3),
    round(median(prot_vc_df$residual_sd, na.rm = TRUE), 3),
    "",
    round(median(betas, na.rm = TRUE), 5),
    round(quantile(abs(betas), 0.75, na.rm = TRUE), 5),
    round(max(abs(betas), na.rm = TRUE), 5),
    round(n_potential_hits / sum(!is.na(correlations)), 3)
  ),
  stringsAsFactors = FALSE
)

print(params, right = FALSE, row.names = FALSE)


# =========================================
# 8. SIM_PARAMS
#
# only realistically measuered quantities go in here. Design choices (n_proteins, n_true_hits,
# the true effect sizes) and assumptions stay in functions.R, where they are visible choices.
# =========================================

.lnfit <- function(x) {
  x <- x[is.finite(x) & x > 0]
  q <- unname(quantile(x, c(0.25, 0.5, 0.75)))
  list(meanlog = log(q[2]), sdlog = (log(q[3]) - log(q[1])) / (2 * qnorm(0.75)))
}
.as_probs <- function(tab) { v <- as.numeric(tab); names(v) <- names(tab); v / sum(v) }

.slope_fit <- .lnfit(prot_vc_df$slope_sd)
.int_fit <- .lnfit(prot_vc_df$intercept_sd)

sim_params <- list(
  N = nrow(visits_both),

  #  ALSFRS mixed model
  alsfrs_intercept_mean = unname(fe_alsfrs["(Intercept)"]),
  alsfrs_slope_mean = unname(fe_alsfrs["time_years"]),
  alsfrs_intercept_sd = sqrt(intercept_var),
  alsfrs_slope_sd = sqrt(slope_var),
  alsfrs_rho = cov_int_slope / (sqrt(intercept_var) * sqrt(slope_var)),
  alsfrs_residual_sd = sqrt(residual_var),

  # ALSFRS visit structure, >=2-DRAW patients only
  visit_count_probs = .as_probs(table(visits_both$n_visits)),
  visit_gap_mean = mean(visit_gaps$gap, na.rm = TRUE),
  visit_gap_sd = sd(visit_gaps$gap, na.rm = TRUE),

  # ALSFRS visit structure, FULL ALS cohort
  visit_count_probs_cohort = .as_probs(table(visit_counts$n_visits)),
  N_cohort = nrow(visit_counts),

  # SomaScan draw counts, same patients
  draw_count_probs = .as_probs(table(pvc2$n_prot_visits)),

  # protein variance components
  prot_intercept_mean = mean(prot_vc_df$fixed_int, na.rm = TRUE),
  prot_residual_sd = median(prot_vc_df$residual_sd, na.rm = TRUE),
  prot_slope_sd_meanlog = .slope_fit$meanlog,
  prot_slope_sd_sdlog = .slope_fit$sdlog,
  prot_intercept_sd_meanlog = .int_fit$meanlog,
  prot_intercept_sd_sdlog = .int_fit$sdlog,

  # effect-size distributions
  beta_75  = unname(quantile(abs(betas), 0.75, na.rm = TRUE)),
  beta_max = max(abs(betas), na.rm = TRUE),
  gamma_deattenuated_75  = unname(quantile(abs(gamma_cors), 0.75, na.rm = TRUE)) / atten,
  gamma_deattenuated_max = max(abs(gamma_cors), na.rm = TRUE) / atten,

  # informative dropout targets
  dropout_r_slope_visits = unname(ct$estimate),
  dropout_r_slope_visits_subset = unname(ct_sub$estimate),
  mean_visits_subset = mean(visits_both$n_visits),

  # onset and Delta FS calibration inputs
  onset_duration_meanlog = mean(log(dur_sub$duration_years)),
  onset_duration_sdlog = sd(log(dur_sub$duration_years)),
  dfs_r_slope_subset = unname(ct_dfs$estimate),

  provenance = list(generated = Sys.time(),
                    script = "01_estimate_parameters.R",
                    draw_timing = "plasma collection date, R/00_visit_time.R",
                    protein_key = "SomaScan SeqId",
                    n_proteins_sampled = nrow(prot_vc_df),
                    R_version = R.version.string)
)

cat("\n========================================================\n")
cat("  SIM_PARAMS -> functions.R  (measured quantities only)\n")
cat("==========================================================\n\n")
for (nm in setdiff(names(sim_params),
                   c("provenance", "visit_count_probs", "draw_count_probs",
                     "visit_count_probs_cohort")))
  cat(sprintf("  %-26s %s\n", nm, format(sim_params[[nm]], digits = 4)))
cat("  visit_count_probs         ",
    paste0(names(sim_params$visit_count_probs), ":",
           round(sim_params$visit_count_probs, 3), collapse = " "), "\n")
cat("  visit_count_probs_cohort  ",
    paste0(names(sim_params$visit_count_probs_cohort), ":",
           round(sim_params$visit_count_probs_cohort, 3), collapse = " "), "\n")
cat("  draw_count_probs          ",
    paste0(names(sim_params$draw_count_probs), ":",
           round(sim_params$draw_count_probs, 3), collapse = " "), "\n")

save(
  params, prot_vc_df, betas, correlations,
  ols_slopes, merged_slopes, visit_counts, visit_gaps,
  prot_visit_counts, prot_gaps, prot_span, visits_both,
  sim_params,
  file = paste0(here(), "/output/simulation_parameters.RData")
)
cat("\nSaved to output/simulation_parameters.RData\n")



cat("\n===============================================\n")
cat("  9. CALIBRATING DROPOUT AND ONSET RECALL\n")
cat("=================================================\n\n")

source("R/functions.R")   # loads the parameters saved above
set.seed(2025)

dropout_cal <- calibrate_dropout(target_r = sim_params$dropout_r_slope_visits_subset)
sim_params$dropout_par <- dropout_cal

recall_cal <- calibrate_recall_cv(target_r = sim_params$dfs_r_slope_subset,
                                  dropout_par = dropout_cal)
sim_params$onset_recall_cv <- recall_cal$cv
sim_params$onset_recall_calibration <- recall_cal

save(
  params, prot_vc_df, betas, correlations,
  ols_slopes, merged_slopes, visit_counts, visit_gaps,
  prot_visit_counts, prot_gaps, prot_span, visits_both,
  sim_params,
  file = paste0(here(), "/output/simulation_parameters.RData")
)
cat("\nCalibrations added to output/simulation_parameters.RData\n")
