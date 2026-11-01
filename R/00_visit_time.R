################################################################################
#  00. VISIT TIME HELPER
################################################################################

visit_number <- function(visit_name) {
  v <- as.character(visit_name)
  out <- suppressWarnings(as.integer(sub(".*Visit (\\d+)$", "\\1", v)))
  out[grepl("Screening", v)] <- 1L
  out[!grepl("Screening|Visit \\d+$", v)] <- NA_integer_
  out
}

# Returns one row per (GUID, TimePoint) with time_years and time_source.
soma_visit_time <- function(plasma_path, alsfrs_path) {
  pl <- read.csv(plasma_path, stringsAsFactors = FALSE) %>%
    dplyr::transmute(GUID = SubjectUID, TimePoint = visit_number(Visit_Name),
                     t_plasma = as.numeric(Visit_Date)) %>%
    dplyr::filter(!is.na(TimePoint), !is.na(t_plasma)) %>%
    dplyr::distinct(GUID, TimePoint, .keep_all = TRUE)

  fx <- read.csv(alsfrs_path, stringsAsFactors = FALSE) %>%
    dplyr::transmute(GUID = SubjectUID, TimePoint = visit_number(Visit_Name),
                     t_alsfrs = as.numeric(Visit_Date)) %>%
    dplyr::filter(!is.na(TimePoint), !is.na(t_alsfrs)) %>%
    dplyr::distinct(GUID, TimePoint, .keep_all = TRUE)

  vt <- dplyr::full_join(pl, fx, by = c("GUID", "TimePoint")) %>%
    dplyr::mutate(
      time_days = dplyr::coalesce(t_plasma, t_alsfrs),
      time_source = dplyr::case_when(!is.na(t_plasma) ~ "plasma_draw_date",
                                     !is.na(t_alsfrs) ~ "alsfrs_visit_date",
                                     TRUE ~ "none"),
      lag_days = t_plasma - t_alsfrs,
      time_years = time_days / 365.25
    ) %>%
    dplyr::filter(time_source != "none")

  cat("  Visit-time table:", nrow(vt), "GUID x visit rows |",
      sum(vt$time_source == "plasma_draw_date"), "plasma-dated |",
      sum(vt$time_source == "alsfrs_visit_date"), "ALSFRS-dated fallback\n")
  ld <- vt$lag_days[!is.na(vt$lag_days)]
  cat("  Plasma draw vs ALSFRS-R same-visit lag: n =", length(ld),
      "| identical:", sprintf("%.1f%%", 100 * mean(ld == 0)),
      "| median |lag| =", median(abs(ld)), "d | 95th pct =", quantile(abs(ld), 0.95), "d\n")
  vt %>% dplyr::select(GUID, TimePoint, time_years, time_source, lag_days)
}

# Attach time to a SomaScan participant-sample data frame (SubjectUID, TimePoint, ...).
attach_soma_time <- function(df, vt, verbose = TRUE) {
  out <- df %>%
    dplyr::rename(GUID = SubjectUID) %>%
    dplyr::inner_join(vt, by = c("GUID", "TimePoint"))
  if (verbose) cat("  SomaScan draws with a usable time:", nrow(out), "of", nrow(df),
                   "(dropped", nrow(df) - nrow(out), "with no draw date)\n")
  out
}
