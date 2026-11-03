################################################################################
#  Two-stage estimators against a one-stage bivariate model or Why not a joint model?
################################################################################

source("R/functions.R")
if (!requireNamespace("glmmTMB", quietly = TRUE)) install.packages("glmmTMB")

RNGkind("L'Ecuyer-CMRG")
set.seed(2025)
if (file.exists("output/simulation_parameters.RData"))
  load("output/simulation_parameters.RData")
source("R/design_check.R"); design_check()

n_sim_bv  <- 100   # simulated datasets
n_null_bv <- 20    # null proteins fitted per dataset

cat("=================================================\n")
cat("  ONE-STAGE vs TWO-STAGE (", n_sim_bv, "sims,", n_true_hits, "hits +",
    n_null_bv, "nulls each)\n")
cat("=================================================\n\n")

run_bv_sim <- function(sim_id) {
  suppressPackageStartupMessages({library(lme4); library(lmerTest); library(dplyr)})

  dat <- suppressWarnings(generate_data(N, n_proteins, n_true_hits, beta_true_values))

  ids <- c(seq_len(n_true_hits), sample((n_true_hits + 1):n_proteins, n_null_bv))
  pd <- dat$prot_df[dat$prot_df$protein_id %in% ids, ]
  pd$protein_id <- match(pd$protein_id, ids)
  k <- length(ids)
  args <- list(als_df = dat$als_df, prot_df = pd, n_proteins = k)

  r1 <- tryCatch(suppressWarnings(do.call(method_twostage_ols, args)), error = function(e) NULL)
  r2 <- tryCatch(suppressWarnings(do.call(method_twostage_blup, args)), error = function(e) NULL)
  r7 <- tryCatch(suppressWarnings(do.call(method_regression_calibration, args)), error = function(e) NULL)
  r8 <- method_hybrid(r1, r2)

  r9 <- tryCatch(method_bivariate(dat$als_df, pd, seq_len(k)),
                 error = function(e) data.frame(
                   protein_id = seq_len(k), est = NA_real_, r_latent = NA_real_,
                   converged = FALSE, error = paste("method_bivariate:", conditionMessage(e)),
                   se = NA_real_, pval = NA_real_, method = "OneStage_Bivariate",
                   stringsAsFactors = FALSE))

  # True latent slope correlation: b1_prot = N(0, slope_sd) + beta * b1_als.
  true_beta <- dat$true_betas[ids]
  true_r <- true_beta * alsfrs_slope_sd /
            sqrt(dat$slope_sds[ids]^2 + true_beta^2 * alsfrs_slope_sd^2)

  bind_rows(r1, r2, r7, r8, r9) %>%
    mutate(sim = sim_id,
           orig_id = ids[protein_id],
           is_true_hit = orig_id <= n_true_hits,
           true_beta = true_beta[protein_id],
           true_r = true_r[protein_id])
}

t0 <- Sys.time()
bv <- bind_rows(par_apply(seq_len(n_sim_bv), run_bv_sim,
                          export = c("run_bv_sim", "n_null_bv")))
cat("Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Rows by method:\n"); print(table(bv$method)); cat("\n")

# ------ Convergence of the bivariate fits ------------
bv_one <- bv %>% filter(method == "OneStage_Bivariate")
expected <- n_sim_bv * (n_true_hits + n_null_bv)
if (nrow(bv_one) == 0 || !"error" %in% names(bv_one)) {
  stop("No bivariate rows came back at all.")
}
bv_conv <- bv_one %>%
  summarise(n_rows = n(), n_estimated = sum(is.finite(est)),
            pct_pdHess = round(100 * mean(converged), 1))
cat("Bivariate fits:", bv_conv$n_estimated, "estimated of", expected,
    "| positive-definite Hessian:", bv_conv$pct_pdHess, "%\n")
if (any(!is.na(bv_one$error))) {
  cat("Most common failure messages:\n")
  print(head(sort(table(bv_one$error), decreasing = TRUE), 5))
}
if (bv_conv$n_estimated == 0) stop("Every bivariate fit failed.")
cat("\n")

# --- Bias and RMSE on the true hits, every method against the same estimand ----
bv_bias <- bv %>% filter(is_true_hit, is.finite(est)) %>%
  group_by(method, true_beta) %>%
  summarise(n_est = n(), mean_est = mean(est),
            bias = mean(est - true_beta),
            rel_bias = mean((est - true_beta) / true_beta),
            rmse = sqrt(mean((est - true_beta)^2)), .groups = "drop") %>%
  arrange(true_beta, abs(rel_bias))
cat("BIAS ON TRUE HITS (est targets true_beta for every method)\n")
print(as.data.frame(bv_bias), digits = 4, row.names = FALSE)

# --- Coverage for the methods with a standard error ----------
bv_cov <- compute_coverage(bv)
cat("\n95% CI COVERAGE (two-stage methods; the bivariate fit has no SE here)\n")
print(as.data.frame(bv_cov), digits = 4, row.names = FALSE)

# --- Behaviour on nulls (true_beta = 0) ------------
bv_null <- bv %>% filter(!is_true_hit, is.finite(est)) %>%
  group_by(method) %>%
  summarise(n = n(), mean_est = mean(est), sd_est = sd(est),
            mean_abs_est = mean(abs(est)), .groups = "drop")
cat("\nNULL PROTEINS (true_beta = 0)\n")
print(as.data.frame(bv_null), digits = 4, row.names = FALSE)

# --- Latent correlation recovered by the bivariate model ------------
bv_r <- bv %>% filter(method == "OneStage_Bivariate", is.finite(r_latent)) %>%
  group_by(is_true_hit) %>%
  summarise(n = n(), mean_true_r = mean(true_r), mean_r_latent = mean(r_latent),
            median_abs_err = median(abs(r_latent - true_r)),
            cor_with_truth = if (n() > 2 && sd(true_r) > 0) cor(r_latent, true_r) else NA_real_,
            .groups = "drop")
cat("\n LATENT SLOPE CORRELATION, bivariate model against truth ──\n")
print(as.data.frame(bv_r), digits = 4, row.names = FALSE)

save(bv, bv_conv, bv_bias, bv_cov, bv_null, bv_r,
     file = "output/onestage_results.RData")
cat("\nSaved to output/onestage_results.RData\n")
