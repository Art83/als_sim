################################################################################
#  BLUP power precision + variance-source diagnostic (essentially doing part A of script 02 better)
################################################################################

source("R/functions.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(2025)
source("R/design_check.R"); design_check()
n_cores <- get_cores()

n_sim_blup <- 300

cat("Cores:", n_cores, "\n\n")
cat("====================================================================\n")
cat("  BLUP PRECISION (", n_sim_blup, "sims, 7 methods, no Baseline_Pred)\n")
cat("====================================================================\n\n")

run_blup_sim <- function(sim_id) {
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
  r7 <- tryCatch(suppressWarnings(do.call(method_regression_calibration, args)), error=function(e) NULL)
  r9 <- tryCatch(suppressWarnings(do.call(method_repeated_measures, args)), error=function(e) NULL)

  hs <- dat$slope_sds[seq_len(n_true_hits)]

  bind_rows(r1, r2, r3, r4, r5, r7, r9) %>%
    mutate(sim = sim_id,
           is_true_hit  = protein_id <= n_true_hits,
           protein_type = dat$protein_type[protein_id],
           true_beta = dat$true_betas[protein_id],
           hit_sd_mean = mean(hs),
           hit_sd_min = min(hs))
}

t0 <- Sys.time()
blup_results <- bind_rows(par_apply(seq_len(n_sim_blup), run_blup_sim,
                                   export = "run_blup_sim"))
cat("Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Methods:\n"); print(table(blup_results$method)); cat("\n")


# Headline metrics as in 02 so the numbers are comparable
pf <- compute_power_fdr(blup_results)
bias <- compute_bias(blup_results)
coverage <- compute_coverage(blup_results)

cat("POWER & FDR (", n_sim_blup, "sims)\n")
print(as.data.frame(pf$summary), digits = 3, row.names = FALSE)
cat("\n BIAS\n")
print(as.data.frame(bias), digits = 4, row.names = FALSE)
cat("\nCOVERAGE\n")
print(as.data.frame(coverage), digits = 4, row.names = FALSE)


per_sim <- blup_results %>%
  filter(!is.na(pval)) %>%
  group_by(method, sim) %>%
  mutate(p_fdr = p.adjust(pval, "fdr")) %>%
  summarise(power = sum(p_fdr < 0.05 & is_true_hit) / n_true_hits,
            hit_sd_mean = first(hit_sd_mean),
            hit_sd_min = first(hit_sd_min),
            .groups = "drop")

bl <- per_sim %>% filter(method == "TwoStage_BLUP")

cat("\n====================================================\n")
cat("  BLUP power and FDR with Monte Carlo intervals\n")
cat("======================================================\n\n")

bl_fdr <- pf$per_sim %>% filter(method == "TwoStage_BLUP")
m  <- mean(bl$power); s <- sd(bl$power); n <- nrow(bl)
mf <- mean(bl_fdr$fdr); sf <- sd(bl_fdr$fdr)
cat("  Power over", n, "sims: mean", round(m, 4), "| 95% CI",
    round(m - 1.96 * s / sqrt(n), 4), "-", round(m + 1.96 * s / sqrt(n), 4), "\n")
cat("  FDR over", nrow(bl_fdr), "sims: mean", round(mf, 4), "| 95% CI",
    round(mf - 1.96 * sf / sqrt(nrow(bl_fdr)), 4), "-",
    round(mf + 1.96 * sf / sqrt(nrow(bl_fdr)), 4), "\n")


cat("\n========================================================\n")
cat("  Does power track the true markers' drawn slope SD?\n")
cat("==========================================================\n\n")

ct <- suppressWarnings(cor.test(bl$power, bl$hit_sd_mean, method = "spearman"))
cat("  Spearman rho(power, mean hit slope_sd) =", round(ct$estimate, 3),
    "  p =", signif(ct$p.value, 3), "\n\n")

tert <- bl %>%
  mutate(sd_tertile = ntile(hit_sd_mean, 3)) %>%
  group_by(sd_tertile) %>%
  summarise(n = n(),
            mean_hit_sd = mean(hit_sd_mean),
            mean_power = mean(power),
            sd_power = sd(power), .groups = "drop")
print(as.data.frame(tert), digits = 4, row.names = FALSE)

cat("\n  Signal for a beta=0.025 hit is beta * als_slope_sd =",
    round(0.025 * alsfrs_slope_sd, 4), "(log2/yr)\n")
cat("  Hit slope_sd across sims: range",
    round(min(bl$hit_sd_mean), 4), "-", round(max(bl$hit_sd_mean), 4), "\n\n")

save(blup_results, per_sim, pf, bias, coverage,
     file = "output/blup_precision.RData")
cat("\nSaved to output/blup_precision.RData\n")
