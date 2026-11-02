################################################################################
#  Simulation comparison of the screening methods
#  Part A realistic design, Part B visit-count sweep, Part D mechanisms, Part C figures.
# Note: I ran on both Windows and Linux. there's a setup to accommodate both. Linux is preferable.
# This is one of the time-consuming scripts
################################################################################

source("R/functions.R")
source("R/design_check.R"); design_check()

RNGkind("L'Ecuyer-CMRG")
set.seed(2025)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

if (file.exists("output/simulation_parameters.RData"))
  load("output/simulation_parameters.RData")

n_cores <- get_cores()
cat("Cores:", n_cores, "\n\n")


.cl <- commandArgs(trailingOnly = TRUE)
run_parts <- if (length(.cl)) toupper(trimws(strsplit(.cl[1], ",")[[1]])) else c("A", "B", "D", "C")
cat("Running parts:", paste(run_parts, collapse = ", "), "\n\n")

# ========================================
# PART A. MAIN SIMULATION
# ========================================


n_sim_main <- 100

if ("A" %in% run_parts) {
cat("==========================================================\n")
cat("  PART A: MAIN SIMULATION (", n_sim_main, "sims, 9 methods)\n")
cat("==========================================================\n\n")

run_main_sim <- function(sim_id) {
  suppressPackageStartupMessages({library(lme4); library(lmerTest); library(dplyr)})

  dat <- suppressWarnings(
    generate_data(N, n_proteins, n_true_hits, beta_true_values))

  args <- list(als_df = dat$als_df, prot_df = dat$prot_df,
               n_proteins = n_proteins, onset_offset = dat$onset_offset)

  r1 <- tryCatch(suppressWarnings(do.call(method_twostage_ols, args)), error=function(e) NULL)
  r2 <- tryCatch(suppressWarnings(do.call(method_twostage_blup, args)), error=function(e) NULL)
  r3 <- tryCatch(suppressWarnings(do.call(method_cross_sectional, args)), error=function(e) NULL)
  r4 <- tryCatch(suppressWarnings(do.call(method_crude_rate, args)), error=function(e) NULL)
  r5 <- tryCatch(suppressWarnings(do.call(method_dichotomise, args)), error=function(e) NULL)
  r6 <- tryCatch(suppressWarnings(do.call(method_baseline_predictor, args)), error=function(e) NULL)
  r7 <- tryCatch(suppressWarnings(do.call(method_regression_calibration, args)), error=function(e) NULL)
  r8 <- method_hybrid(r1, r2)
  r9 <- tryCatch(suppressWarnings(do.call(method_repeated_measures, args)), error=function(e) NULL)

  bind_rows(r1, r2, r3, r4, r5, r6, r7, r8, r9) %>%
    mutate(sim = sim_id,
           is_true_hit = protein_id <= n_true_hits,
           protein_type = dat$protein_type[protein_id],
           true_beta = dat$true_betas[protein_id])
}

t0 <- Sys.time()
main_results <- bind_rows(par_apply(seq_len(n_sim_main), run_main_sim,
                                   export = "run_main_sim"))
cat("Part A done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Methods:\n"); print(table(main_results$method)); cat("\n")

# Metrics
main_pf <- compute_power_fdr(main_results)
main_bias <- compute_bias(main_results)
main_coverage <- compute_coverage(main_results)
main_bias_sel <- compute_bias_selected(main_results)

cat("---POWER & FDR ------\n")
print(as.data.frame(main_pf$summary), digits = 3, row.names = FALSE)
cat("\n------ BIAS (slope-estimating methods)----\n")
print(as.data.frame(main_bias), digits = 4, row.names = FALSE)
cat("\n------ COVERAGE --------\n")
print(as.data.frame(main_coverage), digits = 4, row.names = FALSE)
cat("\n------ BIAS AMONG SELECTED TRUE HITS ----\n")
print(as.data.frame(main_bias_sel), digits = 4, row.names = FALSE)

save(main_results, main_pf, main_bias, main_coverage, main_bias_sel,
     file = "output/simulation_results.RData")
gc()
}  # end Part A


# =============================================
# PART B. SENSITIVITY: VISIT COUNT
# =============================================

visit_scenarios <- c(2, 3, 4, 6, 8, 12)
n_sim_sens <- 50


if ("B" %in% run_parts) {
cohort_counts <- if (!is.null(sim_params$visit_count_probs_cohort)) {
  sim_params$visit_count_probs_cohort
} else stop("sim_params has no visit_count_probs_cohort. Rerun 01.")
cat("Cohort scenario visit distribution:",
    paste0(names(cohort_counts), ":", round(cohort_counts / sum(cohort_counts), 3), collapse = " "), "\n")
cat("\n=============================================\n")
cat("  PART B: SENSITIVITY (", n_sim_sens, " sims x ",
    length(visit_scenarios) + 1, " scenarios)\n")
cat("==============================================\n\n")

# Faster subset: skip baseline_predictor (lmer per protein) in sensitivity
run_sens_sim <- function(sim_id, nv) {
  suppressPackageStartupMessages({library(lme4); library(lmerTest); library(dplyr)})

  if (identical(nv, "cohort")) {
    visit_count_probs <<- cohort_counts / sum(cohort_counts)
    visit_count_vals <<- as.integer(names(cohort_counts))
    dat <- suppressWarnings(
      generate_data(N, n_proteins, n_true_hits, beta_true_values,
                    dropout = FALSE, visit_override = NULL, decouple_draws = FALSE))
  } else {
    dat <- suppressWarnings(
      generate_data(N, n_proteins, n_true_hits, beta_true_values,
                    dropout = FALSE, visit_override = nv))
  }

  args <- list(als_df = dat$als_df, prot_df = dat$prot_df,
               n_proteins = n_proteins, onset_offset = dat$onset_offset)

  r1 <- tryCatch(suppressWarnings(do.call(method_twostage_ols, args)), error=function(e) NULL)
  r2 <- tryCatch(suppressWarnings(do.call(method_twostage_blup, args)), error=function(e) NULL)
  r3 <- tryCatch(suppressWarnings(do.call(method_cross_sectional, args)), error=function(e) NULL)
  r4 <- tryCatch(suppressWarnings(do.call(method_crude_rate, args)), error=function(e) NULL)
  r5 <- tryCatch(suppressWarnings(do.call(method_dichotomise, args)), error=function(e) NULL)
  r7 <- tryCatch(suppressWarnings(do.call(method_regression_calibration, args)), error=function(e) NULL)
  r8 <- method_hybrid(r1, r2)
  r9 <- tryCatch(suppressWarnings(do.call(method_repeated_measures, args)), error=function(e) NULL)

  rel <- tryCatch({
    s <- dat$als_df %>% group_by(id) %>%
      summarise(sxx = sum((time - mean(time))^2), .groups = "drop")
    median(alsfrs_slope_sd^2 / (alsfrs_slope_sd^2 + alsfrs_residual_sd^2 / s$sxx))
  }, error = function(e) NA_real_)

  bind_rows(r1, r2, r3, r4, r5, r7, r8, r9) %>%
    mutate(sim = sim_id, n_visits = as.character(nv),
           median_visits = median(dat$visit_counts), ols_rel_median = rel,
           is_true_hit  = protein_id <= n_true_hits,
           protein_type = dat$protein_type[protein_id],
           true_beta    = dat$true_betas[protein_id])
}

sens_results <- list()
for (nv in c(as.list(visit_scenarios), list("cohort"))) {
  t0 <- Sys.time()
  cat("  ", nv, "visits...")
  res <- bind_rows(par_apply(seq_len(n_sim_sens), function(s) run_sens_sim(s, nv),
                             export = c("run_sens_sim", "nv", "cohort_counts")))
  sens_results[[as.character(nv)]] <- res
  cat(" done (", round(difftime(Sys.time(), t0, units = "mins"), 1), "min)\n")
  gc()
}

sens_df <- bind_rows(sens_results)
sens_df$n_visits <- factor(sens_df$n_visits, levels = c(as.character(visit_scenarios), "cohort"))
cat("\nScenario design (median visits, median per-patient OLS slope reliability):\n")
print(as.data.frame(sens_df %>% distinct(n_visits, sim, median_visits, ols_rel_median) %>%
  group_by(n_visits) %>% summarise(median_visits = median(median_visits),
                                   ols_rel_median = median(ols_rel_median), .groups = "drop")),
  digits = 3, row.names = FALSE)
cat("\nMethods by scenario:\n")
print(table(sens_df$n_visits, sens_df$method))

# Sensitivity metrics
sens_power <- sens_df %>%
  filter(!is.na(pval)) %>%
  group_by(method, n_visits, sim) %>%
  mutate(p_fdr = p.adjust(pval, "fdr")) %>% ungroup() %>%
  mutate(rejected = p_fdr < 0.05) %>%
  group_by(method, n_visits, sim) %>%
  summarise(TP = sum(rejected & is_true_hit), FP = sum(rejected & !is_true_hit),
            FN = sum(!rejected & is_true_hit),
            power = TP / pmax(TP + FN, 1),
            fdr = ifelse(TP+FP>0, FP/(TP+FP), 0), .groups = "drop") %>%
  group_by(method, n_visits) %>%
  summarise(mean_power = mean(power), mean_fdr = mean(fdr),
            mean_TP = mean(TP), mean_FP = mean(FP), .groups = "drop")

sens_bias <- sens_df %>%
  filter(is_true_hit, method %in% c("TwoStage_OLS","TwoStage_BLUP","RegCalibration","Hybrid_BLUPsel_OLSest")) %>%
  group_by(method, n_visits, true_beta) %>%
  summarise(bias = mean(est - true_beta, na.rm = TRUE),
            rel_bias = mean((est - true_beta)/true_beta, na.rm = TRUE),
            .groups = "drop")

sens_coverage <- sens_df %>%
  filter(is_true_hit, method %in% c("TwoStage_OLS","TwoStage_BLUP","RegCalibration","Hybrid_BLUPsel_OLSest"),
         !is.na(se), se > 0) %>%
  mutate(covers = (est-1.96*se <= true_beta) & (true_beta <= est+1.96*se)) %>%
  group_by(method, n_visits) %>%
  summarise(coverage = mean(covers, na.rm = TRUE), .groups = "drop")

cat("\n--- POWER BY VISITS ----\n")
print(as.data.frame(sens_power %>% arrange(n_visits, method)), digits = 3, row.names = FALSE)
cat("\n--- BIAS BY VISITS ----\n")
print(as.data.frame(sens_bias), digits = 4, row.names = FALSE)
cat("\n---- COVERAGE BY VISITS ----\n")
print(as.data.frame(sens_coverage), digits = 4, row.names = FALSE)

sens_bias_sel <- sens_df %>%
  filter(!is.na(pval), method %in% c("TwoStage_OLS","TwoStage_BLUP","RegCalibration","Hybrid_BLUPsel_OLSest")) %>%
  group_by(method, n_visits, sim) %>%
  mutate(p_fdr = p.adjust(pval, "fdr")) %>% ungroup() %>%
  filter(is_true_hit, p_fdr < 0.05) %>%
  group_by(method, n_visits) %>%
  summarise(n_selected = n(), rel_bias = mean((est - true_beta)/true_beta, na.rm = TRUE), .groups = "drop")
cat("\n---- RELATIVE BIAS AMONG SELECTED TRUE HITS, BY VISITS -----\n")
print(as.data.frame(sens_bias_sel %>% arrange(n_visits, method)), digits = 4, row.names = FALSE)

save(sens_df, sens_power, sens_bias, sens_coverage, sens_bias_sel,
     file = "output/sensitivity_results.RData")
gc()
}  # end Part B


# ==============================================
# PART D. MECHANISM ISOLATION
# Which feature breaks two-stage OLS: visit heterogeneity,
# informative dropout, or ALSFRS floor censoring? Vary one at a time.
# ==============================================

# 50, not 200: scales with n_proteins (now 1000).
n_sim_iso <- 50

if ("D" %in% run_parts) {
.chk <- generate_data(N, n_proteins, n_true_hits, beta_true_values)
cat("\nDesign check | ALSFRS visits mean", round(mean(.chk$visit_counts), 2),
    "(real", round(sim_params$mean_visits_subset %||% NA, 2), ")  | SomaScan draws mean",
    round(mean(.chk$draw_counts), 2), "(real",
    round(sum(draw_count_vals * draw_count_probs), 2), ")\n")
rm(.chk)

cat("\n===========================================================\n")
cat("  PART D: MECHANISM ISOLATION (", n_sim_iso, "sims/condition)\n")
cat("=============================================================\n\n")

run_iso_sim <- function(sim_id, informative, nv, censor, dec) {
  suppressPackageStartupMessages({library(lme4); library(lmerTest); library(dplyr)})
  dpar <- if (informative) sim_dropout_par else list(rho = 0)
  dat <- suppressWarnings(
    generate_data(N, n_proteins, n_true_hits, beta_true_values,
                  dropout = TRUE, visit_override = nv, censor = censor,
                  decouple_draws = dec, dropout_par = dpar))
  args <- list(als_df = dat$als_df, prot_df = dat$prot_df, n_proteins = n_proteins)
  r1 <- tryCatch(suppressWarnings(do.call(method_twostage_ols, args)), error=function(e) NULL)
  r2 <- tryCatch(suppressWarnings(do.call(method_twostage_blup, args)), error=function(e) NULL)
  r7 <- tryCatch(suppressWarnings(do.call(method_regression_calibration, args)), error=function(e) NULL)
  r8 <- method_hybrid(r1, r2)
  r9 <- tryCatch(suppressWarnings(do.call(method_repeated_measures, args)), error=function(e) NULL)
  bind_rows(r1, r2, r7, r8, r9) %>%
    mutate(sim = sim_id,
           is_true_hit  = protein_id <= n_true_hits,
           true_beta    = dat$true_betas[protein_id])
}


iso_conditions <- list(
  list(tag = "1 HET + informative dropout, shared",   informative = TRUE,  nv = NULL, censor = TRUE,  dec = FALSE),
  list(tag = "2 HET + uninformative dropout, shared", informative = FALSE, nv = NULL, censor = TRUE,  dec = FALSE),
  list(tag = "3 HOMO-4",                              informative = FALSE, nv = 4L,   censor = TRUE,  dec = FALSE),
  list(tag = "4 HOMO-12 + censor",                    informative = FALSE, nv = 12L,  censor = TRUE,  dec = FALSE),
  list(tag = "5 HOMO-12 + NO censor",                 informative = FALSE, nv = 12L,  censor = FALSE, dec = FALSE),
  list(tag = "6 HET + informative dropout, DECOUPLED (=Part A)",
                                                      informative = TRUE,  nv = NULL, censor = TRUE,  dec = TRUE)
)

iso_all <- list()
for (cond in iso_conditions) {
  t0 <- Sys.time(); cat("  ", cond$tag, "...")
  res <- bind_rows(par_apply(seq_len(n_sim_iso),
    function(s) run_iso_sim(s, cond$informative, cond$nv, cond$censor, cond$dec),
    export = c("run_iso_sim", "cond")))
  res$condition <- cond$tag
  iso_all[[cond$tag]] <- res
  cat(" done (", round(difftime(Sys.time(), t0, units = "mins"), 1), "min)\n"); gc()
}
iso_df <- bind_rows(iso_all)

iso_fdr <- iso_df %>% filter(!is.na(pval)) %>%
  group_by(condition, method, sim) %>%
  mutate(p_fdr = p.adjust(pval, "fdr")) %>% ungroup() %>%
  mutate(rejected = p_fdr < 0.05) %>%
  group_by(condition, method, sim) %>%
  summarise(TP = sum(rejected & is_true_hit), FP = sum(rejected & !is_true_hit),
            FN = sum(!rejected & is_true_hit),
            power = TP/pmax(TP+FN,1), fdr = ifelse(TP+FP>0, FP/(TP+FP), 0), .groups="drop") %>%
  group_by(condition, method) %>%
  summarise(mean_power = mean(power), mean_fdr = mean(fdr), mean_FP = mean(FP), .groups="drop")

iso_bias <- iso_df %>%
  filter(is_true_hit, method %in% c("TwoStage_OLS","TwoStage_BLUP","RegCalibration","Hybrid_BLUPsel_OLSest")) %>%
  group_by(condition, method, true_beta) %>%
  summarise(rel_bias = mean((est - true_beta)/true_beta, na.rm = TRUE), .groups="drop")

cat("\n--- FDR / POWER BY CONDITION (dropout & heterogeneity) ---\n")
print(as.data.frame(iso_fdr %>% arrange(condition, method)), digits = 3, row.names = FALSE)
cat("\n--- RELATIVE BIAS BY CONDITION (censoring: compare 4 vs 5) ---\n")
print(as.data.frame(iso_bias %>% arrange(condition, method)), digits = 4, row.names = FALSE)

dropout_par_used <- sim_dropout_par
save(iso_df, iso_fdr, iso_bias, dropout_par_used, file = "output/isolation_results.RData")
gc()
}  # end Part D


# ==================================================
# PART C FIGURES
# ==================================================

if ("C" %in% run_parts) {
if (!exists("main_pf") && file.exists("output/simulation_results.RData"))
  load("output/simulation_results.RData")
if (!exists("sens_power") && file.exists("output/sensitivity_results.RData"))
  load("output/sensitivity_results.RData")
as_visits <- function(d) d %>% filter(as.character(n_visits) != "cohort") %>%
  mutate(n_visits = as.numeric(as.character(n_visits)))
sens_power <- as_visits(sens_power)
sens_bias <- as_visits(sens_bias)
sens_coverage <- as_visits(sens_coverage)
cat("\n  Generating figures...\n")
theme_pub <- theme_bw(base_size = 12) + theme(legend.position = "bottom")

# Main summary
fig1 <- main_pf$summary %>%
  select(method, mean_power, mean_fdr) %>%
  pivot_longer(-method, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("mean_power","mean_fdr"),
                         labels = c("Power","FDR"))) %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_col(alpha = 0.7) +
  geom_hline(data = data.frame(metric = factor(c("Power","FDR"),
    levels = c("Power","FDR")), target = c(0.8, 0.05)),
    aes(yintercept = target), linetype = "dashed", colour = "red") +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "Methods under Answer ALS conditions",
       x = NULL, y = NULL) +
  theme_pub + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
safe_ggsave("output/figures/fig1_main_summary.png", fig1, width = 10, height = 5, dpi = 200)

# Power x visits
fig2 <- sens_power %>%
  ggplot(aes(x = n_visits, y = mean_power, colour = method, group = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.4) +
  scale_x_continuous(breaks = visit_scenarios) +
  labs(title = "Power by visit count", x = "Visits per patient",
       y = "Mean power (FDR < 0.05)") +
  theme_pub + coord_cartesian(ylim = c(0,1))
safe_ggsave("output/figures/fig2_power_visits.png", fig2, width = 9, height = 6, dpi = 200)

# FDR x visits
fig3 <- sens_power %>%
  ggplot(aes(x = n_visits, y = mean_fdr, colour = method, group = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
  scale_x_continuous(breaks = visit_scenarios) +
  labs(title = "FDR by visit count", x = "Visits per patient",
       y = "Mean FDR") +
  theme_pub + coord_cartesian(ylim = c(0, 1))
safe_ggsave("output/figures/fig3_fdr_visits.png", fig3, width = 9, height = 6, dpi = 200)

# Bias x visits
fig4 <- sens_bias %>%
  group_by(method, n_visits) %>%
  summarise(rel_bias = mean(rel_bias), .groups = "drop") %>%
  ggplot(aes(x = n_visits, y = rel_bias, colour = method, group = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = visit_scenarios) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Relative bias (slope-estimating methods)",
       x = "Visits per patient", y = "Relative bias") +
  theme_pub
safe_ggsave("output/figures/fig4_bias_visits.png", fig4, width = 9, height = 6, dpi = 200)

# Coverage x visits
fig5 <- sens_coverage %>%
  ggplot(aes(x = n_visits, y = coverage, colour = method, group = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "red") +
  scale_x_continuous(breaks = visit_scenarios) +
  labs(title = "95% CI coverage", x = "Visits per patient", y = "Coverage") +
  theme_pub + coord_cartesian(ylim = c(0, 1))
safe_ggsave("output/figures/fig5_coverage_visits.png", fig5, width = 9, height = 6, dpi = 200)
}  # end Part C

cat("\n════════════════════════════════════════════════════════\n")
cat("  DONE. Results in output/\n")
cat("════════════════════════════════════════════════════════\n")
