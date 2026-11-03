################################################################################
# Validation in real data. 02 script tells how the methods behaved in simulations. 
# this one puts them in answer als
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4)
  library(lmerTest)
  library(ggplot2); library(parallel); library(here)
})
source(paste0(here(), "/R/00_visit_time.R"))

RNGkind("L'Ecuyer-CMRG")
set.seed(2025)
n_cores <- max(1, detectCores() - 2)
top_n   <- 20
min_pat <- 10   # patients per protein in the second stage

# Parallel apply that works on both platforms. mclapply forks, which Windows cannot
# do, so Windows gets a socket cluster instead.
run_par <- function(X, FUN, export = character(0)) {
  if (.Platform$OS.type == "windows") {
    caller <- parent.frame()
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, suppressPackageStartupMessages({
      library(dplyr); library(tidyr); library(lme4); library(lmerTest)
    }))
    if (length(export)) clusterExport(cl, export, envir = caller)
    clusterSetRNGStream(cl, 2025)
    parLapply(cl, X, FUN)
  } else {
    mclapply(X, FUN, mc.cores = n_cores)
  }
}

# load&prep
cat("Loading data...\n")
CLIN <- paste0(here(), "/raw/als-metadata/answerals/clinical/")

df <- read.csv(paste0(here(), "/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"),
               check.names = FALSE)
mapper <- read.csv(paste0(here(), "/raw/als-metadata/answerals/aals_dataportal_datatable.csv"),
                   check.names = FALSE)
alsfrs  <- read.csv(paste0(CLIN, "ALSFRS_R.csv"))
cbs_raw <- read.csv(paste0(CLIN, "ALS_CBS.csv"))
lns_raw <- read.csv(paste0(CLIN, "CNS_Lability_Scale.csv"))
dem <- read.csv(paste0(CLIN, "Demographics.csv"), check.names = FALSE)
annotation <- read.csv(paste0(here(), "/raw/somascan/documents/annotation.csv"))

# SeqId is the key throughout. 1,125 aptamers share a UniProt accession with another
# aptamer, so UniProt is an annotation column and never a name.
misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
seq_map <- annotation %>%
  transmute(SeqId = as.character(SeqId), uniprot = as.character(UniProt.ID),
            Target.Name, Target.Full.Name) %>%
  filter(SeqId %in% colnames(df), !is.na(uniprot), uniprot != "") %>%
  distinct(SeqId, .keep_all = TRUE)
all_seq <- seq_map$SeqId
df <- df[, c(misc_columns, all_seq)]
cat("Aptamers with a UniProt accession:", length(all_seq),
    "| unique accessions:", n_distinct(seq_map$uniprot), "\n")

als_guids <- mapper %>% filter(SUBJECT_GROUP == "ALS") %>% pull(GUID)

prep_outcome <- function(d, val_col) {
  d$visit_no <- if ("Visit_Name" %in% names(d)) visit_number(d$Visit_Name) else NA_integer_
  d %>%
    transmute(GUID = SubjectUID, Visit_Date = suppressWarnings(as.numeric(Visit_Date)),
              visit_no, value = suppressWarnings(as.numeric(.data[[val_col]]))) %>%
    filter(GUID %in% als_guids, !is.na(value), !is.na(Visit_Date)) %>%
    group_by(GUID) %>% filter(n() >= 2) %>%
    arrange(Visit_Date, .by_group = TRUE) %>%
    mutate(time_years = Visit_Date / 365.25) %>%
    ungroup()
}

alsfrs_proc <- prep_outcome(alsfrs, "alsfrst")
cbs_proc <- prep_outcome(cbs_raw %>% filter(Child_Name == "ALS CBS", cbstot != 0), "cbstot")
lns_proc <- prep_outcome(lns_raw, "cnslstot")
cat("ALS with >= 2 ALSFRS-R:", n_distinct(alsfrs_proc$GUID), "(paper: 663)\n")

# Proteomics with time from the plasma collection date.
vt <- soma_visit_time(paste0(CLIN, "Plasma_Sample.csv"), paste0(CLIN, "ALSFRS_R.csv"))

prot_with_time <- df %>%
  filter(SubjectUID %in% als_guids) %>%
  attach_soma_time(vt) %>%
  filter(GUID %in% alsfrs_proc$GUID) %>%       # needs an ALSFRS-R slope
  group_by(GUID) %>% filter(n() >= 2) %>% ungroup()

longitudinal_guids <- unique(prot_with_time$GUID)
cat("Longitudinal proteomics patients:", length(longitudinal_guids), "(paper: 229)\n")
cat("Draws:", nrow(prot_with_time), "| time source:\n")
print(table(prot_with_time$time_source))


# ======================================
# 1. PRECOMPUTE ALL SLOPES (ONE TIME)
# ======================================

cat("\nPrecomputing slopes...\n")
t0 <- Sys.time()

# Long format, one row per patient x draw x aptamer, log2 RFU.
prot_long <- prot_with_time %>%
  select(GUID, TimePoint, time_years, all_of(all_seq)) %>%
  pivot_longer(all_of(all_seq), names_to = "SeqId", values_to = "y") %>%
  mutate(y = log2(as.numeric(y))) %>%
  filter(is.finite(y))

prot_slopes_ols <- prot_long %>%
  group_by(SeqId, GUID) %>% filter(n() >= 2) %>%
  summarise(prot_slope = cov(y, time_years) / var(time_years), .groups = "drop") %>%
  filter(is.finite(prot_slope))
cat("  OLS protein slopes:", nrow(prot_slopes_ols), "patient x aptamer\n")


blup_one <- function(d) {
  if (sum(table(d$GUID) >= 2) < 20) return(NULL)
  m <- tryCatch(suppressMessages(suppressWarnings(
    lmer(y ~ time_years + (time_years | GUID), data = d,
         control = lmerControl(calc.derivs = FALSE)))),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  re <- ranef(m)$GUID
  data.frame(GUID = rownames(re),
             prot_slope = unname(fixef(m)["time_years"]) + re[["time_years"]],
             stringsAsFactors = FALSE)
}

long_split <- split(prot_long[, c("GUID", "time_years", "y")], prot_long$SeqId)
blup_list <- run_par(names(long_split), function(s) {
  out <- blup_one(long_split[[s]])
  if (!is.null(out)) out$SeqId <- s
  out
}, export = c("blup_one", "long_split"))
prot_slopes_blup <- bind_rows(blup_list)
rm(long_split, blup_list); gc()
cat("  BLUP protein slopes:", n_distinct(prot_slopes_blup$SeqId), "of", length(all_seq),
    "aptamers fitted\n")
saveRDS(prot_slopes_blup, paste0(here(), "/output/validation/protein_blup_slopes.rds"))


first_draw <- prot_with_time %>% group_by(GUID) %>%
  summarise(t0 = min(time_years), .groups = "drop")
prot_baseline <- prot_long %>%
  inner_join(first_draw, by = "GUID") %>%
  filter(time_years == t0) %>%
  distinct(GUID, SeqId, .keep_all = TRUE) %>%
  select(GUID, SeqId, prot_bl = y)


prot_first_rows <- prot_long %>%
  inner_join(first_draw, by = "GUID") %>%
  filter(time_years == t0) %>%
  distinct(GUID, SeqId, .keep_all = TRUE) %>%
  select(GUID, TimePoint, SeqId, y)

# Outcome slopes
compute_outcome_slopes <- function(outcome_df) {
  outcome_df %>%
    group_by(GUID) %>% filter(n() >= 2, var(time_years) > 0) %>%
    summarise(outcome_slope = cov(value, time_years) / var(time_years), .groups = "drop") %>%
    filter(is.finite(outcome_slope))
}

compute_blup_slopes <- function(outcome_df) {
  m <- tryCatch(suppressWarnings(
    lmer(value ~ time_years + (time_years | GUID), data = outcome_df, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa"))),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  fe <- fixef(m)
  ranef(m)$GUID %>% tibble::rownames_to_column("GUID") %>%
    mutate(outcome_slope_blup = time_years + fe["time_years"]) %>%
    select(GUID, outcome_slope_blup)
}

alsfrs_slopes <- compute_outcome_slopes(alsfrs_proc)
alsfrs_blups <- compute_blup_slopes(alsfrs_proc)


onset <- mapper %>%
  filter(SUBJECT_GROUP == "ALS") %>%
  transmute(GUID, Participant_ID,
            age_onset = suppressWarnings(as.numeric(AGE_AT_SYMPTOM_ONSET))) %>%
  left_join(dem %>%
              transmute(Participant_ID, age_screen = suppressWarnings(as.numeric(age))) %>%
              distinct(Participant_ID, .keep_all = TRUE),
            by = "Participant_ID")

alsfrs_bl <- alsfrs_proc %>% group_by(GUID) %>%
  slice_min(Visit_Date, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(GUID, alsfrs_bl = value, bl_day = Visit_Date)

delta_fs <- alsfrs_bl %>%
  inner_join(onset, by = "GUID") %>%
  mutate(months_from_onset = pmax(3, (age_screen - age_onset) * 12 + bl_day / 30.4375),
         delta_fs = (48 - alsfrs_bl) / months_from_onset) %>%
  filter(is.finite(delta_fs), age_screen >= age_onset) %>%
  select(GUID, delta_fs)
cat("  Delta FS from recorded onset:", nrow(delta_fs), "ALS patients |",
    sum(delta_fs$GUID %in% longitudinal_guids), "in the longitudinal set\n")

cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n\n")


# =====================================================
# 2. METHODS ON REAL DATA (vectorized second stage)
# =====================================================

finish <- function(x) {
  x %>%
    mutate(t = r * sqrt((n - 2) / (1 - r^2)), pval = 2 * pt(-abs(t), n - 2),
           p_fdr = p.adjust(pval, "fdr")) %>%
    arrange(pval)
}

# Protein slope against an outcome slope, per aptamer
stage2 <- function(ps, os, xcol) {
  ps %>% inner_join(os, by = "GUID") %>%
    group_by(SeqId) %>% filter(n() >= min_pat) %>%
    summarise(r = cor(prot_slope, .data[[xcol]]), n = n(),
              est = cov(prot_slope, .data[[xcol]]) / var(.data[[xcol]]),
              .groups = "drop") %>%
    finish()
}

# Baseline protein against a per-patient outcome, per aptamer
baseline_vs <- function(bl, ydf, ycol) {
  bl %>% inner_join(ydf, by = "GUID") %>%
    group_by(SeqId) %>% filter(n() >= min_pat) %>%
    summarise(r = cor(prot_bl, .data[[ycol]]), n = n(), .groups = "drop") %>%
    mutate(est = r) %>%
    finish()
}

run_all_real <- function(guid_subset, ps_ols, ps_blup, prot_bl, prot_rows,
                         outcome_sl, outcome_blup, dfs_df, outcome_long) {
  po <- ps_ols  %>% filter(GUID %in% guid_subset)
  pb <- ps_blup %>% filter(GUID %in% guid_subset)
  os <- outcome_sl %>% filter(GUID %in% guid_subset)
  ob <- if (!is.null(outcome_blup)) outcome_blup %>% filter(GUID %in% guid_subset) else NULL
  bl <- prot_bl %>% filter(GUID %in% guid_subset)

  results <- list()

  # M1: Two-stage OLS (OLS slopes on both sides)
  results$OLS <- stage2(po, os, "outcome_slope")

  if (!is.null(ob)) {
    # M2: Two-stage BLUP (BLUP slopes on both sides, as simulated and as in the paper)
    results$BLUP <- stage2(pb, ob, "outcome_slope_blup")
    # M2b: BLUP outcome slope against OLS protein slope (the July "BLUP")
    results$BLUP_ALSFRSonly <- stage2(po, ob, "outcome_slope_blup")
    # M8: Hybrid. Select on the BLUP p-value, report the OLS coefficient.
    results$Hybrid <- results$OLS %>% select(SeqId, r, n, est) %>%
      inner_join(results$BLUP %>% select(SeqId, pval, p_fdr), by = "SeqId") %>%
      arrange(pval)
  }

  ol <- outcome_long %>% filter(GUID %in% guid_subset, !is.na(visit_no)) %>%
    group_by(GUID, visit_no) %>%
    summarise(outcome_val = mean(value), .groups = "drop")
  if (nrow(ol) > 0) {
    pcs <- prot_rows %>% filter(GUID %in% guid_subset) %>%
      inner_join(ol, by = c("GUID", "TimePoint" = "visit_no"))
    if (nrow(pcs) > 0)
      results$CrossSectional <- pcs %>% group_by(SeqId) %>% filter(n() >= 20) %>%
        summarise(r = cor(y, outcome_val), n = n(), .groups = "drop") %>%
        mutate(est = r) %>% finish()
  }

  # M4: Crude rate Delta FS (recorded onset)
  dfs <- dfs_df %>% filter(GUID %in% guid_subset)
  if (nrow(dfs) > 0) results$CrudeRate <- baseline_vs(bl, dfs, "delta_fs")

  # M5: Dichotomise at the median OLS outcome slope
  med_slope <- median(os$outcome_slope)
  groups <- os %>% mutate(group = ifelse(outcome_slope < med_slope, "fast", "slow"))
  results$Dichotomise <- bl %>% inner_join(groups %>% select(GUID, group), by = "GUID") %>%
    group_by(SeqId) %>%
    filter(sum(group == "fast") >= 5, sum(group == "slow") >= 5) %>%
    summarise({
      tt <- tryCatch(t.test(prot_bl[group == "fast"], prot_bl[group == "slow"]),
                     error = function(e) NULL)
      if (!is.null(tt)) data.frame(est = unname(tt$statistic), pval = tt$p.value)
      else data.frame(est = NA_real_, pval = NA_real_)
    }, .groups = "drop") %>%
    filter(!is.na(pval)) %>%
    mutate(p_fdr = p.adjust(pval, "fdr")) %>% arrange(pval)

  if (!is.null(ob)) results$Baseline_BLUP <- baseline_vs(bl, ob, "outcome_slope_blup")

  # M7: Regression calibration. Pooled residual variance from patients with residual
  rc <- outcome_long %>% filter(GUID %in% guid_subset) %>%
    group_by(GUID) %>% filter(n() >= 2, var(time_years) > 0) %>%
    summarise(slope = cov(value, time_years) / var(time_years),
              ssx   = sum((time_years - mean(time_years))^2),
              rss   = sum((value - mean(value) - slope * (time_years - mean(time_years)))^2),
              df_resid = n() - 2, .groups = "drop")
  has_df <- rc$df_resid > 0
  sigma2_pooled <- sum(rc$rss[has_df]) / sum(rc$df_resid[has_df])
  lambda <- max(0.01, 1 - mean(sigma2_pooled / rc$ssx) / var(rc$slope))
  results$RegCalibration <- results$OLS %>% mutate(est = est / lambda)

  top_hits <- lapply(results, function(r) if (!is.null(r)) head(r$SeqId, top_n) else character(0))
  list(results = results, top_hits = top_hits, lambda = lambda)
}


# ==============================================
# 3. FULL-DATA RESULTS
# ==============================================

cat("===================================\n")
cat("  FULL DATA RESULTS\n")
cat("===================================\n\n")

res_alsfrs <- run_all_real(longitudinal_guids, prot_slopes_ols, prot_slopes_blup,
                           prot_baseline, prot_first_rows, alsfrs_slopes, alsfrs_blups,
                           delta_fs, alsfrs_proc)

cat("Regression-calibration lambda (ALSFRS-R):", round(res_alsfrs$lambda, 3), "\n\n")
cat("FDR-significant per method (aptamers | unique proteins):\n")
for (m in names(res_alsfrs$results)) {
  sig <- res_alsfrs$results[[m]] %>% filter(p_fdr < 0.05)
  cat(sprintf("  %-16s %4d | %4d\n", m, nrow(sig),
              n_distinct(seq_map$uniprot[match(sig$SeqId, seq_map$SeqId)])))
}

cat("\n")
for (m in names(res_alsfrs$results)) {
  cat(" Method:", m, "\n")
  top <- res_alsfrs$results[[m]] %>% head(top_n) %>%
    left_join(seq_map, by = "SeqId")
  print(top %>% select(SeqId, uniprot, Target.Name, est, pval, p_fdr) %>%
          as.data.frame(), digits = 4, row.names = FALSE)
  cat("\n")
  genes <- top %>% pull(Target.Name) %>% na.omit() %>% unique()
  writeLines(genes, paste0(here(), "/output/validation/top_genes_", m, ".txt"))
}

# --- CHECK AGAINST THE PAPER --------------
# The BLUP row here and the paper's script 02 use the same data, the same draw
# timing and the same two models, so their hit lists should agree. hopefully
paper_candidates <- c(paste0(here(), "/output/paper/progression_with_grounding.csv"),
                      "D:/AZ/als_long/output/grounding/progression_with_grounding.csv")
paper_file <- paper_candidates[file.exists(paper_candidates)][1]
paper_check <- NULL
if (!is.na(paper_file)) {
  pp <- read.csv(paper_file, stringsAsFactors = FALSE) %>% mutate(SeqId = as.character(SeqId))
  ours <- res_alsfrs$results$BLUP
  paper_hits <- pp$SeqId[pp$p_fdr < 0.05]
  our_hits   <- ours$SeqId[ours$p_fdr < 0.05]
  cmp <- inner_join(ours %>% select(SeqId, r_here = r),
                    pp %>% select(SeqId, r_paper = cor_r), by = "SeqId")
  paper_check <- list(n_paper = length(paper_hits), n_here = length(our_hits),
                      n_shared = length(intersect(paper_hits, our_hits)),
                      r_agreement = cor(cmp$r_here, cmp$r_paper),
                      max_abs_diff = max(abs(cmp$r_here - cmp$r_paper)))
  cat("======================================\n")
  cat("  BLUP HITS AGAINST THE PAPER (", basename(dirname(paper_file)), ")\n")
  cat("======================================\n")
  cat("  Paper FDR hits:", paper_check$n_paper, "| here:", paper_check$n_here,
      "| shared:", paper_check$n_shared, "\n")
  cat("  Correlation of per-aptamer r across", nrow(cmp), "aptamers:",
      round(paper_check$r_agreement, 4), "| largest |difference|:",
      signif(paper_check$max_abs_diff, 3), "\n\n")
} else {
  cat("Paper progression file not found. Copy progression_with_grounding.csv from\n")
}


# ===================================================================
# 3b. BASELINE PREDICTOR (full data only) + SLOPE-vs-BASELINE VERDICT
# ===================================================================

cat("\n==================================================\n")
cat("  BASELINE PREDICTOR (mixed model, full data)\n")
cat("====================================================\n\n")

run_baseline_pred <- function(guid_subset, prot_bl, outcome_long) {
  bl  <- prot_bl %>% filter(GUID %in% guid_subset)
  als <- outcome_long %>% filter(GUID %in% guid_subset) %>% select(GUID, time_years, value)
  bl_split <- split(bl[, c("GUID", "prot_bl")], bl$SeqId)

  out <- run_par(names(bl_split), function(s) {
    suppressPackageStartupMessages({library(lme4); library(lmerTest); library(dplyr)})
    d <- als %>% inner_join(bl_split[[s]], by = "GUID")
    if (length(unique(d$GUID)) < 10) return(NULL)
    d$prot_bl <- as.numeric(scale(d$prot_bl))   # z-score: est comparable across aptamers
    m <- tryCatch(suppressWarnings(
      lmer(value ~ time_years * prot_bl + (time_years | GUID), data = d,
           REML = TRUE, control = lmerControl(optimizer = "bobyqa"))),
      error = function(e) NULL)
    if (is.null(m)) return(NULL)
    sm <- summary(m)$coefficients
    rn <- "time_years:prot_bl"
    if (!rn %in% rownames(sm)) return(NULL)
    if (!"Pr(>|t|)" %in% colnames(sm))
      stop("lmerTest is not attached -- summary() has no Pr(>|t|) column.")
    data.frame(SeqId = s, n = length(unique(d$GUID)),
               est = sm[rn, "Estimate"], se = sm[rn, "Std. Error"],
               pval = sm[rn, "Pr(>|t|)"], stringsAsFactors = FALSE)
  }, export = c("bl_split", "als"))

  bind_rows(out) %>% mutate(p_fdr = p.adjust(pval, "fdr")) %>% arrange(pval)
}

t0 <- Sys.time()
res_blpred <- run_baseline_pred(longitudinal_guids, prot_baseline, alsfrs_proc)
cat("Baseline_Pred done in", round(difftime(Sys.time(), t0, units = "mins"), 1),
    "min |", nrow(res_blpred), "aptamers fitted\n")
cat("FDR-significant:", sum(res_blpred$p_fdr < 0.05, na.rm = TRUE), "\n\n")

if (nrow(res_blpred) > 0) {
  print(res_blpred %>% head(top_n) %>% left_join(seq_map, by = "SeqId") %>%
          select(SeqId, uniprot, Target.Name, est, pval, p_fdr) %>%
          as.data.frame(), digits = 4, row.names = FALSE)
}

# ---- REPEATED-MEASURES MIXED MODEL (field-standard comparator, full data only) ----
run_repeated_measures <- function(guid_subset, prot_rows, outcome_long) {
  ol <- outcome_long %>% filter(GUID %in% guid_subset, !is.na(visit_no)) %>%
    group_by(GUID, visit_no) %>% summarise(alsfrs = mean(value), .groups = "drop")
  d_all <- prot_rows %>% filter(GUID %in% guid_subset) %>%
    inner_join(ol, by = c("GUID", "TimePoint" = "visit_no"))
  rm_split <- split(d_all[, c("GUID", "y", "alsfrs")], d_all$SeqId)
  out <- run_par(names(rm_split), function(s) {
    d <- rm_split[[s]]
    if (length(unique(d$GUID)) < 20) return(NULL)
    m <- tryCatch(suppressWarnings(suppressMessages(
      lmer(alsfrs ~ y + (1 | GUID), data = d, REML = TRUE,
           control = lmerControl(calc.derivs = FALSE)))), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    sm <- summary(m)$coefficients
    if (!"y" %in% rownames(sm) || !"Pr(>|t|)" %in% colnames(sm)) return(NULL)
    data.frame(SeqId = s, n = length(unique(d$GUID)), n_obs = nrow(d),
               est = sm["y", "Estimate"], se = sm["y", "Std. Error"],
               pval = sm["y", "Pr(>|t|)"], stringsAsFactors = FALSE)
  }, export = "rm_split")
  bind_rows(out) %>% mutate(p_fdr = p.adjust(pval, "fdr")) %>% arrange(pval)
}

cat("\n=================================================\n")
cat("  REPEATED-MEASURES MIXED MODEL (full data)\n")
cat("===================================================\n\n")
t0 <- Sys.time()
res_rm <- run_repeated_measures(longitudinal_guids, prot_long, alsfrs_proc)
cat("RepeatedMeasures_LMM done in", round(difftime(Sys.time(), t0, units = "mins"), 1),
    "min |", nrow(res_rm), "aptamers fitted | median draws per aptamer:",
    median(res_rm$n_obs), "\n")
cat("FDR-significant:", sum(res_rm$p_fdr < 0.05, na.rm = TRUE), "\n")
rm_hits <- res_rm$SeqId[res_rm$p_fdr < 0.05]
blup_hits <- res_alsfrs$results$BLUP$SeqId[res_alsfrs$results$BLUP$p_fdr < 0.05]
cat("Overlap with the BLUP hits:", length(intersect(rm_hits, blup_hits)), "of",
    length(blup_hits), "BLUP hits are also repeated-measures hits\n\n")
print(res_rm %>% head(top_n) %>% left_join(seq_map, by = "SeqId") %>%
        select(SeqId, uniprot, Target.Name, est, pval, p_fdr) %>%
        as.data.frame(), digits = 4, row.names = FALSE)

# THE VERDICT
n_hits <- function(x) if (is.null(x)) NA_integer_ else sum(x$p_fdr < 0.05, na.rm = TRUE)

verdict <- data.frame(
  uses_protein = c("SLOPE (needs >=2 draws)", "SLOPE (needs >=2 draws)",
                   "SLOPE (needs >=2 draws)",
                   "BASELINE (1 draw)", "BASELINE (1 draw)",
                   "BASELINE (1 draw)", "BASELINE (1 draw)"),
  method = c("TwoStage_BLUP", "TwoStage_OLS", "BLUP_ALSFRSonly",
             "Baseline_BLUP", "Baseline_Pred", "Dichotomise", "CrudeRate"),
  fair_comparator = c("-", "-", "-", "YES: identical to BLUP but baseline",
                      "YES: principled mixed model",
                      "no: median split loses power", "no: crude rate from recorded onset"),
  hits = c(n_hits(res_alsfrs$results$BLUP), n_hits(res_alsfrs$results$OLS),
           n_hits(res_alsfrs$results$BLUP_ALSFRSonly),
           n_hits(res_alsfrs$results$Baseline_BLUP), n_hits(res_blpred),
           n_hits(res_alsfrs$results$Dichotomise), n_hits(res_alsfrs$results$CrudeRate)),
  stringsAsFactors = FALSE
)
verdict <- rbind(verdict, data.frame(
  uses_protein = "LEVELS (all draws, same-visit ALSFRS-R)",
  method = "RepeatedMeasures_LMM",
  fair_comparator = "YES: the field-standard longitudinal model",
  hits = n_hits(res_rm), stringsAsFactors = FALSE))

cat("\n===================================================\n")
cat("  DOES THE LONGITUDINAL PROTEIN ARM EARN ITS COST?\n")
cat("  (same", length(longitudinal_guids), "patients, same outcome, same FDR)\n")
cat("=====================================================\n\n")
print(verdict, row.names = FALSE)

if (!is.null(res_alsfrs$results$Baseline_BLUP)) {
  sl_hits <- res_alsfrs$results$BLUP %>% filter(p_fdr < 0.05) %>% pull(SeqId)
  bl_hits <- res_alsfrs$results$Baseline_BLUP %>% filter(p_fdr < 0.05) %>% pull(SeqId)
  bp_hits <- res_blpred %>% filter(p_fdr < 0.05) %>% pull(SeqId)
  cat("\nHit overlap:\n")
  cat("  BLUP only     :", length(setdiff(sl_hits, union(bl_hits, bp_hits))), "\n")
  cat("  baseline only :", length(setdiff(union(bl_hits, bp_hits), sl_hits)), "\n")
  cat("  both          :", length(intersect(sl_hits, union(bl_hits, bp_hits))), "\n")
}

save(res_blpred, res_rm, verdict,
     file = paste0(here(), "/output/validation/baseline_pred.RData"))


# =================================
# 4. SPLIT-HALF STABILITY
# =================================

cat("==================================\n")
cat("  SPLIT-HALF STABILITY\n")
cat("==================================\n\n")

n_splits <- 100
t0 <- Sys.time()

split_results <- run_par(seq_len(n_splits), function(i) {
  guids <- longitudinal_guids
  n_half <- floor(length(guids) / 2)
  idx <- sample(length(guids))
  disc <- guids[idx[1:n_half]]
  val <- guids[idx[(n_half + 1):length(guids)]]

  d <- run_all_real(disc, prot_slopes_ols, prot_slopes_blup, prot_baseline, prot_first_rows,
                    alsfrs_slopes, alsfrs_blups, delta_fs, alsfrs_proc)
  v <- run_all_real(val,  prot_slopes_ols, prot_slopes_blup, prot_baseline, prot_first_rows,
                    alsfrs_slopes, alsfrs_blups, delta_fs, alsfrs_proc)

  sapply(names(d$top_hits), function(m)
    length(intersect(d$top_hits[[m]], v$top_hits[[m]])) / top_n)
}, export = c("run_all_real", "stage2", "finish", "baseline_vs", "min_pat", "top_n",
              "longitudinal_guids", "prot_slopes_ols", "prot_slopes_blup",
              "prot_baseline", "prot_first_rows", "alsfrs_slopes", "alsfrs_blups",
              "delta_fs", "alsfrs_proc"))

cat("Split-half done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n\n")

split_df <- do.call(rbind, split_results) %>% as.data.frame() %>%
  pivot_longer(everything(), names_to = "method", values_to = "replication_rate")

split_summary <- split_df %>% group_by(method) %>%
  summarise(mean_repl = mean(replication_rate),
            sd_repl = sd(replication_rate), .groups = "drop")

cat("Replication rate (top", top_n, "):\n")
print(as.data.frame(split_summary), digits = 3, row.names = FALSE)

save(res_alsfrs, res_rm, split_df, split_summary, seq_map, longitudinal_guids, paper_check,
     file = paste0(here(), "/output/validation/validation_results.RData"))
cat("Saved interim output/validation/validation_results.RData\n")

make_figures <- TRUE

if (make_figures) {
fig_split <- split_df %>%
  ggplot(aes(x = method, y = replication_rate, fill = method)) +
  geom_boxplot(alpha = 0.6) +
  labs(title = "Split-half replication rate", x = NULL, y = "Fraction replicated") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  coord_cartesian(ylim = c(0, 1))
split_png <- paste0(here(), "/output/figures/fig_val_splithalf.png")
tryCatch({ png(split_png, width = 8, height = 5, units = "in", res = 200)
           print(fig_split); dev.off()
           cat("  Saved:", split_png, "\n") },
         error = function(e) { try(dev.off(), silent = TRUE)
                               cat("  figure not saved:", conditionMessage(e), "\n") })
} else cat("Split-half figure skipped (make_figures = FALSE)\n")


# =====================================
# 5. CROSS-OUTCOME CONSISTENCY
# =====================================

cat("\n==================================\n")
cat("  CROSS-OUTCOME CONSISTENCY\n")
cat("====================================\n\n")


run_cross_outcome <- FALSE
res_cbs <- res_lns <- cross_outcome <- NULL
if (!run_cross_outcome) cat("Skipped (run_cross_outcome = FALSE)\n")
if (run_cross_outcome) {
cbs_slopes <- compute_outcome_slopes(cbs_proc)
lns_slopes <- compute_outcome_slopes(lns_proc)
cbs_guids <- intersect(unique(cbs_slopes$GUID), longitudinal_guids)
lns_guids <- intersect(unique(lns_slopes$GUID), longitudinal_guids)
cat("CBS n=", length(cbs_guids), " LNS n=", length(lns_guids), "\n")

cbs_blups <- compute_blup_slopes(cbs_proc)
lns_blups <- compute_blup_slopes(lns_proc)

empty_dfs <- data.frame(GUID = character(0), delta_fs = numeric(0))

res_cbs <- run_all_real(cbs_guids, prot_slopes_ols, prot_slopes_blup, prot_baseline, prot_first_rows,
                        cbs_slopes, cbs_blups, empty_dfs, cbs_proc)
res_lns <- run_all_real(lns_guids, prot_slopes_ols, prot_slopes_blup, prot_baseline, prot_first_rows,
                        lns_slopes, lns_blups, empty_dfs, lns_proc)

jaccard <- function(a, b) {
  u <- length(union(a, b)); if (u == 0) NA else length(intersect(a, b)) / u
}

cross_outcome <- data.frame()
for (m in intersect(names(res_alsfrs$top_hits),
                    intersect(names(res_cbs$top_hits), names(res_lns$top_hits)))) {
  cross_outcome <- rbind(cross_outcome, data.frame(
    method = m,
    ALS_CBS = jaccard(res_alsfrs$top_hits[[m]], res_cbs$top_hits[[m]]),
    ALS_LNS = jaccard(res_alsfrs$top_hits[[m]], res_lns$top_hits[[m]]),
    CBS_LNS = jaccard(res_cbs$top_hits[[m]], res_lns$top_hits[[m]])))
}
cross_outcome$mean_j <- rowMeans(cross_outcome[, 2:4], na.rm = TRUE)
cat("Cross-outcome Jaccard:\n")
print(as.data.frame(cross_outcome), digits = 3, row.names = FALSE)
}  # end run_cross_outcome


# ============================
# 6. SAVE
# ===========================

save(res_alsfrs, res_cbs, res_lns, split_df, split_summary, cross_outcome,
     seq_map, longitudinal_guids, paper_check, res_rm,
     file = paste0(here(), "/output/validation/validation_results.RData"))

cat("\n=====================================\n")
cat("  VALIDATION COMPLETE\n")
cat("=======================================\n")
