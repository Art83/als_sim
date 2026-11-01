################################################################################
#  functions.R Data generation, screening methods
################################################################################

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. Measured params. What we can realistically measure in data ------------

.pfile <- if (requireNamespace("here", quietly = TRUE))
  file.path(here::here(), "output", "simulation_parameters.RData") else
  "output/simulation_parameters.RData"

sim_params <- NULL
if (file.exists(.pfile)) {
  .e <- new.env(); load(.pfile, envir = .e); sim_params <- .e$sim_params; rm(.e)
}
if (is.null(sim_params)) {
  stop("functions.R: no sim_params in ", .pfile, ". Run 01_estimate_parameters.R first.",
       call. = FALSE)
} else {
  message("functions.R: parameters loaded from ", basename(.pfile),
          " (01 run ", format(sim_params$provenance$generated), ", ",
          sim_params$provenance$n_proteins_sampled, " proteins sampled)")
}

N <- sim_params$N
alsfrs_intercept_mean <- sim_params$alsfrs_intercept_mean
alsfrs_slope_mean <- sim_params$alsfrs_slope_mean
alsfrs_intercept_sd <- sim_params$alsfrs_intercept_sd
alsfrs_slope_sd <- sim_params$alsfrs_slope_sd
alsfrs_rho <- sim_params$alsfrs_rho
alsfrs_residual_sd <- sim_params$alsfrs_residual_sd

visit_count_probs <- sim_params$visit_count_probs / sum(sim_params$visit_count_probs)
visit_count_vals <- as.integer(names(visit_count_probs))
visit_gap_mean <- sim_params$visit_gap_mean
visit_gap_sd <- sim_params$visit_gap_sd

draw_count_probs <- sim_params$draw_count_probs / sum(sim_params$draw_count_probs)
draw_count_vals <- as.integer(names(draw_count_probs))

prot_intercept_mean <- sim_params$prot_intercept_mean
prot_residual_sd <- sim_params$prot_residual_sd
prot_slope_sd_meanlog <- sim_params$prot_slope_sd_meanlog
prot_slope_sd_sdlog <- sim_params$prot_slope_sd_sdlog
prot_intercept_sd_meanlog <- sim_params$prot_intercept_sd_meanlog
prot_intercept_sd_sdlog <- sim_params$prot_intercept_sd_sdlog

# --- 2. CHOSEN -  decisions informed by Measured but not directly estimated -------
n_proteins  <- 1000
n_true_hits <- 5

beta_true_values  <- c(0.015, 0.025)   

gamma_true_values <- c(0.15, 0.25)     

# --- 3. CALIBRATED in 01 against summary of the cohort-----
onset_recall_cv <- sim_params$onset_recall_cv
sim_dropout_par <- sim_params$dropout_par
onset_dur_meanlog <- sim_params$onset_duration_meanlog
onset_dur_sdlog <- sim_params$onset_duration_sdlog
if (is.null(onset_recall_cv) || is.null(sim_dropout_par))
  message("functions.R: dropout and onset calibrations are not in sim_params yet. ")


# DATA GENERATION

generate_visits <- function(n) {
  nv <- sample(visit_count_vals, n, replace = TRUE, prob = visit_count_probs)
  shape <- (visit_gap_mean / visit_gap_sd)^2
  rate  <- visit_gap_mean / visit_gap_sd^2
  lapply(nv, function(k) {
    if (k == 1) return(0)
    gaps <- rgamma(k - 1, shape = shape, rate = rate)
    gaps <- pmax(gaps, 0.08)
    sort(c(0, cumsum(gaps)))
  })
}

generate_data <- function(N, n_proteins, n_true_hits, beta_values,
                          gamma_values = gamma_true_values,
                          dropout = TRUE, visit_override = NULL,
                          censor = TRUE, recall_cv = onset_recall_cv,
                          decouple_draws = TRUE, dropout_par = sim_dropout_par) {
  Sigma_als <- matrix(c(
    alsfrs_intercept_sd^2,
    alsfrs_rho * alsfrs_intercept_sd * alsfrs_slope_sd,
    alsfrs_rho * alsfrs_intercept_sd * alsfrs_slope_sd,
    alsfrs_slope_sd^2), 2, 2)

  re_als <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_als)
  b0_als <- alsfrs_intercept_mean + re_als[, 1]
  b1_als <- alsfrs_slope_mean     + re_als[, 2]

  onset <- simulate_onset(b0_als, b1_als, recall_cv)
  onset_true <- onset$true; onset_offset <- onset$recalled

  if (!is.null(visit_override)) {
    shape <- (visit_gap_mean / visit_gap_sd)^2
    rate  <- visit_gap_mean / visit_gap_sd^2
    visit_times <- lapply(seq_len(N), function(i) {
      if (visit_override == 1) return(0)
      gaps <- pmax(rgamma(visit_override - 1, shape = shape, rate = rate), 0.08)
      sort(c(0, cumsum(gaps)))
    })
  } else if (dropout) {
    if (is.null(dropout_par))
      stop("No calibrated dropout in sim_params.")
    visit_times <- monotone_schedule(b1_als, dropout_par)
  } else {
    visit_times <- generate_visits(N)
  }

  # Somascan blood sampling 

  if (decouple_draws && is.null(visit_override)) {
    k_draw <- sample(draw_count_vals, N, replace = TRUE, prob = draw_count_probs)
    draw_times <- lapply(seq_len(N), function(i)
      visit_times[[i]][seq_len(min(k_draw[i], length(visit_times[[i]])))])
  } else {
    draw_times <- visit_times
  }

  als_list <- lapply(seq_len(N), function(i) {
    t <- visit_times[[i]]
    y <- b0_als[i] + b1_als[i] * t + rnorm(length(t), 0, alsfrs_residual_sd)
    if (censor) y <- pmin(pmax(y, 0), 48) 
    data.frame(id = i, time = t, alsfrs = y)
  })
  als_df <- do.call(rbind, als_list)

  slope_sd_vec <- rlnorm(n_proteins, prot_slope_sd_meanlog, prot_slope_sd_sdlog)
  int_sd_vec <- rlnorm(n_proteins, prot_intercept_sd_meanlog, prot_intercept_sd_sdlog)
  protein_type <- ifelse(seq_len(n_proteins) <= n_true_hits, "hit",
                  ifelse(slope_sd_vec > exp(prot_slope_sd_meanlog), "variable", "null"))

  beta_vec <- rep(0, n_proteins)
  for (h in seq_len(n_true_hits))
    beta_vec[h] <- beta_values[((h - 1) %% length(beta_values)) + 1]

  gamma_vec <- rep(0, n_proteins)
  for (h in seq_len(n_true_hits))
    gamma_vec[h] <- gamma_values[((h - 1) %% length(gamma_values)) + 1]

  z_als <- as.numeric(scale(b1_als))

  prot_data_list <- list()
  for (p in seq_len(n_proteins)) {
    b0_prot <- prot_intercept_mean +
               int_sd_vec[p] * (gamma_vec[p] * z_als +
                                sqrt(1 - gamma_vec[p]^2) * rnorm(N))
    b1_prot <- rnorm(N, 0, slope_sd_vec[p]) + beta_vec[p] * b1_als
    prot_list <- lapply(seq_len(N), function(i) {
      t <- draw_times[[i]]
      y <- b0_prot[i] + b1_prot[i] * t + rnorm(length(t), 0, prot_residual_sd)
      data.frame(id = i, time = t, prot = y, protein_id = p)
    })
    prot_data_list[[p]] <- do.call(rbind, prot_list)
  }
  prot_df <- do.call(rbind, prot_data_list)

  list(als_df = als_df, prot_df = prot_df,
       true_slopes = b1_als, true_betas = beta_vec, true_gammas = gamma_vec,
       protein_type = protein_type,
       slope_sds = slope_sd_vec,
       int_sds = int_sd_vec,
       onset_offset = onset_offset, onset_true = onset_true,
       visit_counts = sapply(visit_times, length),
       draw_counts = sapply(draw_times, length))
}


# Methods for comparisons
# 1: tow stage OLS

method_twostage_ols <- function(als_df, prot_df, n_proteins, ...) {
  als_slopes <- als_df %>% group_by(id) %>% filter(n() >= 2) %>%
    summarise(als_slope = cov(alsfrs, time) / var(time), .groups = "drop")

  prot_slopes <- prot_df %>% group_by(protein_id, id) %>% filter(n() >= 2) %>%
    summarise(prot_slope = cov(prot, time) / var(time), .groups = "drop")

  merged <- inner_join(prot_slopes, als_slopes, by = "id")

  merged %>% group_by(protein_id) %>% filter(n() >= 10) %>%
    summarise(
      est = cov(prot_slope, als_slope) / var(als_slope),
      r = cor(prot_slope, als_slope),
      pval = {n <- n(); t <- r * sqrt((n-2)/(1-r^2)); 2*pt(-abs(t), n-2)},
      se = {
        fit <- lm(prot_slope ~ als_slope)
        summary(fit)$coefficients[2, 2]
      },
      .groups = "drop"
    ) %>%
    mutate(method = "TwoStage_OLS") %>%
    select(protein_id, est, se, pval, method)
}


# 2: two stage BLUP

method_twostage_blup <- function(als_df, prot_df, n_proteins, ...) {
  m_als <- tryCatch(suppressWarnings(
    lmer(alsfrs ~ time + (time | id), data = als_df, REML = TRUE,
         # lme4 refuses to fit when n_obs <= n_random_effects
         control = lmerControl(optimizer = "bobyqa", check.nobs.vs.nRE = "ignore"))),
    error = function(e) NULL)
  if (is.null(m_als)) return(NULL)

  fe_als <- fixef(m_als)
  als_slopes <- ranef(m_als)$id %>%
    tibble::rownames_to_column("id") %>%
    mutate(id = as.integer(id), als_slope = time + fe_als["time"]) %>%
    select(id, als_slope)

  prot_slopes_list <- list()
  for (p in seq_len(n_proteins)) {
    d <- prot_df[prot_df$protein_id == p, ]
    m_p <- tryCatch(suppressWarnings(
      lmer(prot ~ time + (time | id), data = d, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa", check.nobs.vs.nRE = "ignore"))),
      error = function(e) NULL)
    if (is.null(m_p)) next
    fe_p <- fixef(m_p)
    sl <- ranef(m_p)$id %>%
      tibble::rownames_to_column("id") %>%
      mutate(id = as.integer(id), prot_slope = time + fe_p["time"]) %>%
      select(id, prot_slope)
    sl$protein_id <- p
    prot_slopes_list[[p]] <- sl
  }
  prot_slopes <- do.call(rbind, prot_slopes_list)
  merged <- inner_join(prot_slopes, als_slopes, by = "id")

  merged %>% group_by(protein_id) %>% filter(n() >= 10) %>%
    summarise(
      est = cov(prot_slope, als_slope) / var(als_slope),
      r = cor(prot_slope, als_slope),
      pval = {n <- n(); t <- r*sqrt((n-2)/(1-r^2)); 2*pt(-abs(t), n-2)},
      se = tryCatch(summary(lm(prot_slope ~ als_slope))$coefficients[2,2],
                       error = function(e) NA),
      .groups = "drop"
    ) %>%
    mutate(method = "TwoStage_BLUP") %>%
    select(protein_id, est, se, pval, method)
}


# 3: Cross-sectional correlation, first draw per patient

method_cross_sectional <- function(als_df, prot_df, n_proteins, ...) {
  first <- prot_df %>% group_by(protein_id, id) %>%
    slice_min(time, n = 1, with_ties = FALSE) %>% ungroup()
  merged <- inner_join(first, als_df, by = c("id", "time"))

  merged %>% group_by(protein_id) %>% filter(n() >= 20) %>%
    summarise(
      r = cor(prot, alsfrs),
      n = n(),
      t = r * sqrt((n - 2) / (1 - r^2)),
      pval = 2 * pt(-abs(t), df = n - 2),
      est = r,
      se = sqrt((1 - r^2) / (n - 2)),
      .groups = "drop"
    ) %>%
    mutate(method = "CrossSectional") %>%
    select(protein_id, est, se, pval, method)
}



# 4: Cruse rate deltaFS


method_crude_rate <- function(als_df, prot_df, n_proteins, onset_offset, ...) {
  als_bl <- als_df %>% group_by(id) %>%
    filter(time == min(time)) %>% slice(1) %>% ungroup() %>%
    rename(alsfrs_bl = alsfrs)

  # delta FS per patient
  dfs <- data.frame(id = seq_along(onset_offset),
                    delta_fs = (48 - als_bl$alsfrs_bl[match(seq_along(onset_offset), als_bl$id)]) /
                               (onset_offset * 12)) 
  dfs <- dfs %>% filter(is.finite(delta_fs), delta_fs > 0)

  prot_bl <- prot_df %>% group_by(protein_id, id) %>%
    filter(time == min(time)) %>% slice(1) %>% ungroup() %>%
    rename(prot_bl = prot)

  merged <- inner_join(prot_bl, dfs, by = "id")

  merged %>% group_by(protein_id) %>% filter(n() >= 10) %>%
    summarise(
      r = cor(prot_bl, delta_fs),
      n = n(),
      t = r * sqrt((n - 2) / (1 - r^2)),
      pval = 2 * pt(-abs(t), df = n - 2),
      est = r,
      se = sqrt((1 - r^2) / (n - 2)),
      .groups = "drop"
    ) %>%
    mutate(method = "CrudeRate_DFS") %>%
    select(protein_id, est, se, pval, method)
}


# 5: dichotomise FAST vs SLOW

method_dichotomise <- function(als_df, prot_df, n_proteins, ...) {
  als_slopes <- als_df %>% group_by(id) %>% filter(n() >= 2) %>%
    summarise(als_slope = cov(alsfrs, time) / var(time), .groups = "drop")

  med <- median(als_slopes$als_slope)
  als_slopes$group <- ifelse(als_slopes$als_slope < med, "fast", "slow")

  prot_bl <- prot_df %>% group_by(protein_id, id) %>%
    filter(time == min(time)) %>% slice(1) %>% ungroup() %>%
    rename(prot_bl = prot)

  merged <- inner_join(prot_bl, als_slopes %>% select(id, group), by = "id")

  merged %>% group_by(protein_id) %>%
    filter(sum(group == "fast") >= 5, sum(group == "slow") >= 5) %>%
    summarise({
      fast_vals <- prot_bl[group == "fast"]
      slow_vals <- prot_bl[group == "slow"]
      tt <- tryCatch(t.test(fast_vals, slow_vals), error = function(e) NULL)
      if (!is.null(tt)) {
        data.frame(est = tt$statistic, se = NA_real_, pval = tt$p.value)
      } else {
        data.frame(est = NA_real_, se = NA_real_, pval = NA_real_)
      }
    }, .groups = "drop") %>%
    filter(!is.na(pval)) %>%
    mutate(method = "Dichotomise") %>%
    select(protein_id, est, se, pval, method)
}


# 6: baseline ~ rate
# lmer(alsfrs ~ time * baseline_protein + (time | id))
# Target: time:baseline_protein interaction

method_baseline_predictor <- function(als_df, prot_df, n_proteins, ...) {
  # Baseline protein per patient per protein
  prot_bl <- prot_df %>% group_by(protein_id, id) %>%
    filter(time == min(time)) %>% slice(1) %>% ungroup() %>%
    rename(prot_bl = prot)

  results_list <- list()
  for (p in seq_len(n_proteins)) {
    bl <- prot_bl %>% filter(protein_id == p) %>% select(id, prot_bl)
    if (nrow(bl) < 30) next

    bl$prot_bl <- as.numeric(scale(bl$prot_bl))
    d <- inner_join(als_df, bl, by = "id")
    if (n_distinct(d$id) < 30) next

    m <- tryCatch(suppressWarnings(
      lmer(alsfrs ~ time * prot_bl + (time | id), data = d, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa"))),
      error = function(e) NULL)

    if (!is.null(m)) {
      s <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
      if (!is.null(s) && "time:prot_bl" %in% rownames(s)) {
        results_list[[p]] <- data.frame(
          protein_id = p,
          est  = s["time:prot_bl", "Estimate"],
          se   = s["time:prot_bl", "Std. Error"],
          pval = s["time:prot_bl", "Pr(>|t|)"])
      }
    }
  }
  res <- do.call(rbind, results_list)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$method <- "Baseline_Pred"
  res %>% select(protein_id, est, se, pval, method)
}


# 7: regression calibration
# Two-stage OLS with attenuation correction.
# beta_corrected = beta_naive / lambda, where lambda = 1 - mean(SE^2) / var(slopes)


method_regression_calibration <- function(als_df, prot_df, n_proteins, ...) {
  als_sl <- als_df %>% group_by(id) %>% filter(n() >= 2) %>%
    summarise(
      als_slope = cov(alsfrs, time) / var(time),
      ssx = sum((time - mean(time))^2),
      rss = {
        slope <- cov(alsfrs, time) / var(time)
        int <- mean(alsfrs) - slope * mean(time)
        sum((alsfrs - (int + slope * time))^2)
      },
      df_resid  = n() - 2,
      .groups = "drop"
    ) %>% filter(is.finite(als_slope), ssx > 0)
  has_df <- als_sl$df_resid > 0
  sigma2_pooled <- if (any(has_df))
    sum(als_sl$rss[has_df]) / sum(als_sl$df_resid[has_df]) else NA_real_

  S2_als <- var(als_sl$als_slope, na.rm = TRUE)
  lambda_als <- if (is.finite(sigma2_pooled) && is.finite(S2_als) && S2_als > 0) {
    sigma2_als <- mean(sigma2_pooled / als_sl$ssx)   # marginal slope error variance, ALL patients
    max(0.01, 1 - sigma2_als / S2_als)
  } else 1   # no reliability info -> no correction

  prot_sl <- prot_df %>% group_by(protein_id, id) %>% filter(n() >= 2) %>%
    summarise(prot_slope = cov(prot, time) / var(time), .groups = "drop")

  merged <- inner_join(prot_sl, als_sl %>% select(id, als_slope), by = "id")

  merged %>% group_by(protein_id) %>% filter(n() >= 10) %>%
    summarise({
      n <- n()
      r <- cor(prot_slope, als_slope)
      s <- summary(lm(prot_slope ~ als_slope))$coefficients
      # Divide beta and se by lambda: corrects attenuation, leaves t (hence p) = naive OLS.
      p_ols <- 2 * pt(-abs(r * sqrt((n - 2)/(1 - r^2))), n - 2)
      data.frame(est = s[2, 1] / lambda_als, se = s[2, 2] / lambda_als, pval = p_ols)
    }, .groups = "drop") %>%
    mutate(method = "RegCalibration") %>%
    select(protein_id, est, se, pval, method)
}


# 8: hybrid. BLUP for selection, OLS for the reported coefficient.
# Built from the OLS and BLUP results already computed

method_hybrid <- function(r_ols, r_blup, ...) {
  if (is.null(r_ols) || is.null(r_blup)) return(NULL)
  inner_join(r_ols  %>% select(protein_id, est, se),
             r_blup %>% select(protein_id, pval), by = "protein_id") %>%
    mutate(method = "Hybrid_BLUPsel_OLSest") %>%
    select(protein_id, est, se, pval, method)
}


# 9: repeated measure mixed model

method_repeated_measures <- function(als_df, prot_df, n_proteins, ...) {
  merged <- inner_join(prot_df, als_df, by = c("id", "time"))
  res <- vector("list", n_proteins)
  for (p in seq_len(n_proteins)) {
    d <- merged[merged$protein_id == p, ]
    if (length(unique(d$id)) < 20) next
    m <- tryCatch(suppressWarnings(suppressMessages(
      lmer(alsfrs ~ prot + (1 | id), data = d, REML = TRUE,
           control = lmerControl(calc.derivs = FALSE)))),
      error = function(e) NULL)
    if (is.null(m)) next
    s <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
    if (is.null(s) || !"prot" %in% rownames(s) || !"Pr(>|t|)" %in% colnames(s)) next
    res[[p]] <- data.frame(protein_id = p,
                           est = s["prot", "Estimate"],
                           se = s["prot", "Std. Error"],
                           pval = s["prot", "Pr(>|t|)"])
  }
  out <- bind_rows(res)
  if (nrow(out) == 0) return(NULL)
  out %>% mutate(method = "RepeatedMeasures_LMM") %>%
    select(protein_id, est, se, pval, method)
}



# metrics

compute_power_fdr <- function(results_df) {
  per_sim <- results_df %>%
    filter(!is.na(pval)) %>%
    group_by(method, sim) %>%
    mutate(p_fdr = p.adjust(pval, method = "fdr")) %>%
    ungroup() %>%
    mutate(rejected = p_fdr < 0.05) %>%
    group_by(method, sim) %>%
    summarise(
      TP = sum(rejected & is_true_hit, na.rm = TRUE),
      FP = sum(rejected & !is_true_hit, na.rm = TRUE),
      FN = sum(!rejected & is_true_hit, na.rm = TRUE),
      power = TP / pmax(TP + FN, 1),
      fdr = ifelse(TP + FP > 0, FP / (TP + FP), 0),
      .groups = "drop")

  summary <- per_sim %>% group_by(method) %>%
    summarise(mean_power = mean(power), sd_power = sd(power),
              mean_fdr = mean(fdr), sd_fdr = sd(fdr),
              mean_TP = mean(TP), mean_FP = mean(FP), .groups = "drop")
  list(per_sim = per_sim, summary = summary)
}

compute_bias <- function(results_df, beta_methods = c("TwoStage_OLS","TwoStage_BLUP","RegCalibration","Hybrid_BLUPsel_OLSest")) {
  results_df %>%
    filter(is_true_hit, method %in% beta_methods) %>%
    group_by(method, true_beta) %>%
    summarise(
      n_est = sum(!is.na(est)), mean_est = mean(est, na.rm = TRUE),
      bias = mean(est - true_beta, na.rm = TRUE),
      rel_bias = mean((est - true_beta)/true_beta, na.rm = TRUE),
      rmse = sqrt(mean((est - true_beta)^2, na.rm = TRUE)),
      .groups = "drop")
}

compute_bias_selected <- function(results_df, beta_methods = c("TwoStage_OLS","TwoStage_BLUP","RegCalibration","Hybrid_BLUPsel_OLSest")) {
  results_df %>%
    filter(!is.na(pval), method %in% beta_methods) %>%
    group_by(method, sim) %>%
    mutate(p_fdr = p.adjust(pval, method = "fdr")) %>% ungroup() %>%
    filter(is_true_hit, p_fdr < 0.05) %>%
    group_by(method, true_beta) %>%
    summarise(
      n_selected = n(), mean_est = mean(est, na.rm = TRUE),
      bias = mean(est - true_beta, na.rm = TRUE),
      rel_bias = mean((est - true_beta)/true_beta, na.rm = TRUE),
      .groups = "drop")
}

compute_coverage <- function(results_df, beta_methods = c("TwoStage_OLS","TwoStage_BLUP","RegCalibration","Hybrid_BLUPsel_OLSest")) {
  results_df %>%
    filter(is_true_hit, method %in% beta_methods, !is.na(se), se > 0) %>%
    mutate(covers = (est - 1.96*se <= true_beta) & (true_beta <= est + 1.96*se)) %>%
    group_by(method) %>%
    summarise(n = n(), coverage = mean(covers, na.rm = TRUE),
              mean_ci_width = mean(2*1.96*se, na.rm = TRUE), .groups = "drop")
}


get_cores <- function() max(1, detectCores() - 2)

# ggsave crashes my windows, inclear why. this is the safe version
safe_ggsave <- function(filename, plot, width = 8, height = 6, dpi = 200, ...) {
  tryCatch({
    png(filename, width = width, height = height, units = "in", res = dpi)
    print(plot); dev.off()
    cat("  Saved:", filename, "\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("  FAILED:", filename, "-", conditionMessage(e), "\n")
  })
}


par_apply <- function(X, FUN, export = character(0), envir = parent.frame(),
                      n_cores = get_cores()) {
  if (.Platform$OS.type == "windows") {
    seed <- sample.int(.Machine$integer.max, 1)
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(wd) {
      setwd(wd)
      suppressMessages(suppressWarnings(source("R/functions.R")))
      NULL
    }, getwd())
    if (length(export)) parallel::clusterExport(cl, export, envir = envir)
    parallel::clusterSetRNGStream(cl, seed)
    parallel::parLapply(cl, X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = n_cores)
  }
}


monotone_schedule <- function(b1, par, max_visits = max(visit_count_vals)) {
  shape <- (visit_gap_mean / visit_gap_sd)^2
  rate  <- visit_gap_mean / visit_gap_sd^2
  z <- -as.numeric(scale(b1))
  h <- plogis(par$base + par$k * z)
  lapply(seq_along(b1), function(i) {
    gaps <- pmax(rgamma(max_visits - 1, shape = shape, rate = rate), 0.08)
    t <- c(0, cumsum(gaps))
    n <- 2
    while (n < max_visits && runif(1) > h[i]) n <- n + 1
    t[seq_len(n)]
  })
}

calibrate_monotone_dropout <- function(target_mean, target_r, n = N, reps = 5,
                                       base_grid = seq(-4, 0, by = 0.2),
                                       k_grid = seq(0, 2, by = 0.1)) {
  Sigma_als <- matrix(c(alsfrs_intercept_sd^2,
                        alsfrs_rho * alsfrs_intercept_sd * alsfrs_slope_sd,
                        alsfrs_rho * alsfrs_intercept_sd * alsfrs_slope_sd,
                        alsfrs_slope_sd^2), 2, 2)
  one <- function(base, k) {
    re <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_als)
    b0 <- alsfrs_intercept_mean + re[, 1]
    b1 <- alsfrs_slope_mean     + re[, 2]
    vt <- monotone_schedule(b1, list(base = base, k = k))
    nv <- lengths(vt)
    sl <- vapply(seq_len(n), function(i) {
      t <- vt[[i]]
      y <- pmin(pmax(b0[i] + b1[i] * t + rnorm(length(t), 0, alsfrs_residual_sd), 0), 48)
      cov(y, t) / var(t)
    }, numeric(1))
    c(mean_visits = mean(nv), r = suppressWarnings(cor(sl, nv)))
  }
  grid <- expand.grid(base = base_grid, k = k_grid)
  res <- t(vapply(seq_len(nrow(grid)), function(g)
    rowMeans(replicate(reps, one(grid$base[g], grid$k[g]))), numeric(2)))
  grid$mean_visits <- res[, 1]; grid$r <- res[, 2]
  grid$dist <- ((grid$mean_visits - target_mean) / target_mean)^2 +
               ((grid$r - target_r) / max(abs(target_r), 0.05))^2
  best <- grid[which.min(grid$dist), ]
  cat(sprintf("  Monotone dropout calibration: base = %.2f, k = %.2f | mean visits %.2f (target %.2f) | r(OLS slope, visits) %.3f (target %.3f)\n",
              best$base, best$k, best$mean_visits, target_mean, best$r, target_r))
  list(base = best$base, k = best$k, achieved_mean = best$mean_visits,
       achieved_r = best$r, target_mean = target_mean, target_r = target_r)
}


simulate_onset <- function(b0, b1, recall_cv) {
  if (is.null(recall_cv)) stop("No calibrated onset recall CV in sim_params. Run 01 to completion.")
  if (is.null(onset_dur_meanlog) || is.null(onset_dur_sdlog))
    stop("No onset duration distribution in sim_params. Rerun 01.")
  n <- length(b0)
  dur_obs <- rlnorm(n, onset_dur_meanlog, onset_dur_sdlog)
  true <- ifelse(b1 < -0.5, pmax(0.25, (48 - b0) / (-b1)), dur_obs)
  list(true = true, recalled = true * exp(rnorm(n, 0, recall_cv)))
}


calibrate_recall_cv <- function(target_r, n = N, reps = 20,
                                cv_grid = seq(0, 2, by = 0.05),
                                dropout_par = sim_dropout_par) {
  Sigma_als <- matrix(c(alsfrs_intercept_sd^2,
                        alsfrs_rho * alsfrs_intercept_sd * alsfrs_slope_sd,
                        alsfrs_rho * alsfrs_intercept_sd * alsfrs_slope_sd,
                        alsfrs_slope_sd^2), 2, 2)
  one <- function(cv) {
    re <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_als)
    b0 <- alsfrs_intercept_mean + re[, 1]
    b1 <- alsfrs_slope_mean     + re[, 2]
    vt <- if (!is.null(dropout_par)) monotone_schedule(b1, dropout_par) else generate_visits(n)
    y_first <- numeric(n); sl <- numeric(n)
    for (i in seq_len(n)) {
      t <- vt[[i]]
      y <- pmin(pmax(b0[i] + b1[i] * t + rnorm(length(t), 0, alsfrs_residual_sd), 0), 48)
      y_first[i] <- y[1]
      sl[i] <- cov(y, t) / var(t)
    }
    on <- simulate_onset(b0, b1, cv)
    dfs <- (48 - y_first) / (on$recalled * 12)
    ok <- is.finite(dfs) & dfs > 0 & is.finite(sl)
    cor(dfs[ok], sl[ok], method = "spearman")
  }
  r <- vapply(cv_grid, function(cv) mean(replicate(reps, one(cv))), numeric(1))
  best <- which.min(abs(r - target_r))
  cat(sprintf("  Delta FS calibration: recall CV = %.2f | Spearman r(Delta FS, OLS slope) %.3f (target %.3f)\n",
              cv_grid[best], r[best], target_r))
  if (abs(r[best] - target_r) > 0.05)
    warning("Delta FS calibration misses its target by more than 0.05; the grid may not reach it.")
  list(cv = cv_grid[best], achieved_r = r[best], target_r = target_r,
       grid = data.frame(cv = cv_grid, r = r))
}



# one stage bivariate mm

fit_bivariate <- function(als_d, prot_d) {
  if (!requireNamespace("glmmTMB", quietly = TRUE))
    stop("glmmTMB is not installed")
  sa <- sd(als_d$alsfrs); sp <- sd(prot_d$prot)
  d <- rbind(
    data.frame(id = as.character(als_d$id), time = als_d$time,
               y = (als_d$alsfrs - mean(als_d$alsfrs)) / sa, outcome = "als"),
    data.frame(id = as.character(prot_d$id), time = prot_d$time,
               y = (prot_d$prot - mean(prot_d$prot)) / sp, outcome = "prot"))
  d$outcome <- factor(d$outcome, levels = c("als", "prot"))
  d$id <- factor(d$id)
  fail <- function(msg) data.frame(est = NA_real_, r_latent = NA_real_,
                                   converged = FALSE, error = msg,
                                   stringsAsFactors = FALSE)
  m <- tryCatch(suppressWarnings(suppressMessages(
    glmmTMB::glmmTMB(y ~ 0 + outcome + outcome:time + (0 + outcome + outcome:time | id),
                     dispformula = ~ 0 + outcome, data = d, REML = TRUE))),
    error = function(e) e)
  if (inherits(m, "error")) return(fail(paste("fit:", conditionMessage(m))))
  S <- tryCatch(as.matrix(getS3method("VarCorr", "glmmTMB")(m)$cond$id),
                error = function(e) e)
  if (inherits(S, "error")) return(fail(paste("VarCorr:", conditionMessage(S))))
  nm <- rownames(S)
  ia <- which(grepl("als", nm) & grepl("time", nm))
  ip <- which(grepl("prot", nm) & grepl("time", nm))
  if (length(ia) != 1 || length(ip) != 1)
    return(fail(paste("term names:", paste(nm, collapse = ", "))))
  data.frame(est       = S[ia, ip] * sp / (S[ia, ia] * sa),
             r_latent  = S[ia, ip] / sqrt(S[ia, ia] * S[ip, ip]),
             converged = isTRUE(m$sdr$pdHess),
             error     = NA_character_, stringsAsFactors = FALSE)
}


method_bivariate <- function(als_df, prot_df, protein_ids, ...) {
  rows <- lapply(protein_ids, function(p) {
    f <- fit_bivariate(als_df, prot_df[prot_df$protein_id == p, ])
    if (is.null(f)) NULL else cbind(protein_id = p, f)
  })
  res <- bind_rows(rows)
  if (nrow(res) == 0) return(NULL)
  res %>% mutate(se = NA_real_, pval = NA_real_, method = "OneStage_Bivariate")
}
