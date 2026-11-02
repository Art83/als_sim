# Compares the simulated design with the cohort before any long run, and stops if it
# does not match. 

design_check <- function(n_rep = 20, tol = 0.10) {
  real <- c(
    visits = sim_params$mean_visits_subset,
    short = unname(sum(sim_params$visit_count_probs[c("2", "3")]) / sum(sim_params$visit_count_probs)),
    draws = sum(draw_count_vals * draw_count_probs),
    r = sim_params$dropout_r_slope_visits_subset)

  s <- replicate(n_rep, {
    d <- suppressWarnings(generate_data(N, 1, 1, beta_true_values))
    sl <- d$als_df %>% group_by(id) %>%
      summarise(s = cov(alsfrs, time) / var(time), n = n(), .groups = "drop")
    c(visits = mean(d$visit_counts), short = mean(d$visit_counts <= 3),
      draws = mean(d$draw_counts), r = cor(sl$s, sl$n))
  })
  sim <- rowMeans(s)

  tab <- data.frame(quantity = c("mean ALSFRS-R visits", "share with <=3 visits",
                                 "mean plasma draws", "r(OLS slope, visits)"),
                    real = real, simulated = sim)
  cat("\nDesign check (", n_rep, "simulated cohorts):\n")
  print(tab, digits = 3, row.names = FALSE)
  off <- abs(tab$simulated - tab$real) > pmax(tol * abs(tab$real), 0.03)
  if (any(off))
    stop("Simulated design does not match the cohort (", paste(tab$quantity[off], collapse = ", "),
         "). Stopping before the long run.", call. = FALSE)
  cat("Design check passed.\n\n")
}
