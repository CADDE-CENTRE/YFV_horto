# scripts/genus_drivers_simple.R
# Genus-specific NB models of monthly abundance vs mean temperature & rainfall
# Inputs (searched):
#   data/climate/SP_met_and_haemagogus_data_MONTHLY.csv (or .CSV) or SP_met.csv
#   data/supplementary/SupplementaryData1.csv (or "Supplementary Data 1.csv")
# Outputs:
#   outputs/supp_fig/SI_Fig_GenusDrivers.(pdf|png)
#   outputs/stats/SI_Table_GenusDrivers.(html|csv)
#   outputs/stats/SI_GenusDrivers_Haemagogus_standardised.csv

# ---------- packages ----------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(MASS)        # glm.nb
  library(broom)       # tidy()
  library(gt)          # SI table
  library(patchwork)   # panel figure
  library(corrplot)
  library(janitor)
  library(scales)
})

# ---------- helpers ----------
pseudo_r2 <- function(m) 1 - m$deviance / m$null.deviance

ok_for_glm <- function(d, y = "n", x = c("tmean_c","rain_cm")) {
  d <- d %>% dplyr::select(dplyr::all_of(c(y, x))) %>% tidyr::drop_na()
  nrow(d) >= 5 && var(d[[y]]) > 0 &&
    sum(purrr::map_dbl(d[x], ~ var(.x, na.rm = TRUE))) > 0
}

fit_nb_safe <- function(d) {
  d <- d %>% tidyr::drop_na(n, tmean_c, rain_cm)
  if (!ok_for_glm(d)) return(NULL)
  suppressWarnings( tryCatch(glm.nb(n ~ tmean_c + rain_cm, data = d),
                             error = function(e) NULL) )
}

tidy_nb <- function(m, genus) {
  if (is.null(m)) return(NULL)
  broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term %in% c("tmean_c","rain_cm")) %>%
    dplyr::mutate(
      genus = genus,
      term  = dplyr::recode(term,
                            tmean_c = "Mean temperature (°C)",
                            rain_cm = "Rainfall (cm)"),
      pseudo_r2 = pseudo_r2(m),
      aic       = AIC(m),
      n_obs     = nobs(m)
    ) %>%
    dplyr::select(genus, term, estimate, conf.low, conf.high, p.value, pseudo_r2, aic, n_obs)
}

# ---------- resolve inputs ----------
clim_cand <- c(
  file.path(DATA_CLIM, "SP_met_and_haemagogus_data_MONTHLY.csv"),
  file.path(DATA_CLIM, "SP_met_and_haemagogus_data_MONTHLY.CSV"),
  file.path(DATA_CLIM, "SP_met.csv")
)
clim_file <- clim_cand[file.exists(clim_cand)][1]
if (is.na(clim_file)) stop("Could not find climate file. Tried:\n", paste(clim_cand, collapse = "\n"))

mosq_cand <- c(
  file.path(DATA_SUPP, "SupplementaryData1.csv"),
  file.path(DATA_SUPP, "Supplementary Data 1.csv")
)
mosq_file <- mosq_cand[file.exists(mosq_cand)][1]
if (is.na(mosq_file)) stop("Could not find mosquito file. Tried:\n", paste(mosq_cand, collapse = "\n"))

dir.create(file.path(OUT, "supp_fig"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "stats"),    showWarnings = FALSE, recursive = TRUE)

# ===================== Data ======================
# 1) Meteorology (first of month dates; cum_rain to cm)
met <- readr::read_csv(clim_file, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  mutate(
    date     = make_date(year, month_original, 1),
    cum_rain = cum_rain * 100
  ) %>%
  transmute(
    date,
    rain_cm   = cum_rain,
    tmin_c    = min_temp,
    tmean_c   = mean_temp,
    tmax_c    = max_temp,
    hmin_pct  = min_hum,
    hmean_pct = mean_hum,
    hmax_pct  = max_hum
  )

# Weather correlation (numeric cols only)
# pdf(file.path(OUT, "supp_fig", "SI_Fig_weather_correlations.pdf"), width = 5, height = 5)
# num_mat <- met %>% dplyr::select(-date) %>% as.matrix()
# corrplot::corrplot(stats::cor(num_mat, use = "pairwise.complete.obs"), method = "circle")
# dev.off()

# 2) Mosquito pools -> monthly counts by genus
mosq_raw <- readr::read_csv(mosq_file, show_col_types = FALSE)
names(mosq_raw) <- tolower(names(mosq_raw))

# robust column picks
date_col <- intersect(c("date", "collection_date"), names(mosq_raw))[1]
genus_col <- intersect(c("genus"), names(mosq_raw))[1]
qty_col <- intersect(c("quantity", "count", "n"), names(mosq_raw))[1]
stopifnot(!is.na(date_col), !is.na(genus_col), !is.na(qty_col))

# parse date robustly
parse_any_date <- function(x){
  d <- suppressWarnings(lubridate::ymd(x))
  if (all(is.na(d))) d <- suppressWarnings(lubridate::dmy(x))
  if (all(is.na(d))) d <- suppressWarnings(lubridate::mdy(x))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x))
  d
}

mosq_mo <- mosq_raw %>%
  transmute(
    date_raw = .data[[date_col]],
    date     = floor_date(parse_any_date(date_raw), "month"),
    genus    = .data[[genus_col]],
    quantity = suppressWarnings(as.numeric(.data[[qty_col]]))
  ) %>%
  mutate(
    genus = case_when(
      str_detect(genus,  regex("^haemagogus", ignore_case = TRUE)) ~ "Haemagogus",
      str_detect(genus,  regex("^aedes",      ignore_case = TRUE)) ~ "Aedes",
      str_detect(genus,  regex("^culex",      ignore_case = TRUE)) ~ "Culex",
      str_detect(genus,  regex("^limatus",    ignore_case = TRUE)) ~ "Limatus",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(date, genus) %>%
  summarise(n = sum(quantity, na.rm = TRUE), .groups = "drop")

# 3) Merge + keep four genera of interest
df <- mosq_mo %>%
  filter(genus %in% c("Haemagogus","Aedes","Limatus","Culex")) %>%
  left_join(met, by = "date") %>%
  arrange(genus, date)

# ------------------- sanity checks -------------------
message("Sanity check: covariates non-NA by genus")
df %>%
  group_by(genus) %>%
  summarise(
    months      = dplyr::n(),
    nonNA_tmean = sum(!is.na(tmean_c)),
    nonNA_rain  = sum(!is.na(rain_cm)),
    var_n       = var(n, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

# ===================== Models ====================
genus_list <- c("Haemagogus","Aedes","Limatus","Culex")

models <- purrr::map(genus_list, ~{
  d <- dplyr::filter(df, genus == .x)
  m <- fit_nb_safe(d)
  list(genus = .x, data = d, model = m)
})

# SI table rows
tidy_rows <- purrr::map_dfr(models, ~ tidy_nb(.x$model, .x$genus))

# Keep failed genera as NA rows so table stays complete
failed <- setdiff(genus_list, unique(tidy_rows$genus))
if (length(failed)) {
  tidy_rows <- bind_rows(
    tidy_rows,
    tibble(
      genus = rep(failed, each = 2),
      term  = rep(c("Mean temperature (°C)", "Rainfall (cm)"), times = length(failed)),
      estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_,
      pseudo_r2 = NA_real_, aic = NA_real_, n_obs = NA_integer_
    )
  )
}

# ===================== Supplementary Table 3  ==================
si_table <- tidy_rows %>%
  mutate(
    `Rate ratio` = estimate,
    `95% CI`     = dplyr::if_else(is.na(conf.low), NA_character_,
                                  sprintf("%.2f–%.2f", conf.low, conf.high)),
    p            = dplyr::if_else(is.na(p.value), NA_character_,
                                  scales::pvalue(p.value, accuracy = 1e-4)),
    `Pseudo-R²`  = dplyr::if_else(is.na(pseudo_r2), NA_character_,
                                  sprintf("%.2f", pseudo_r2)),
    `N (months)` = n_obs
  ) %>%
  dplyr::select(Genus = genus, Term = term, `Rate ratio`, `95% CI`, p, `Pseudo-R²`, AIC = aic, `N (months)`)

gt_tbl <- gt(si_table) %>%
  tab_header(
    title = md("**Genus-specific negative-binomial models of abundance**"),
    subtitle = "Dependent variable: monthly mosquito counts"
  ) %>%
  fmt_number(columns = c(`Rate ratio`, AIC, `N (months)`), decimals = 2)

#gtsave(gt_tbl, file.path(OUT, "stats", "SM_Table3.html"))
readr::write_csv(si_table, file.path(OUT, "stats", "SM_Table3.csv"))

# ===================== 4-panel figure (sanity check, not included in the manuscript) =============
panel_plot <- function(g, m, d) {
  if (is.null(m)) {
    ggplot(d, aes(x = date)) +
      geom_col(aes(y = n), alpha = 0.6) +
      labs(title = paste0(g, " (model not converged)"),
           y = "Mosquitoes captured", x = NULL) +
      theme_bw()
  } else {
    d2 <- d %>% tidyr::drop_na(n, tmean_c, rain_cm) %>% mutate(fit = predict(m, type = "response"))
    ggplot(d2, aes(x = date)) +
      geom_col(aes(y = n), alpha = 0.6) +
      geom_line(aes(y = fit), linewidth = 1) +
      labs(title = g, y = "Mosquitoes captured", x = NULL) +
      theme_bw()
  }
}

plots <- purrr::map(models, ~ panel_plot(.x$genus, .x$model, .x$data))
fig4 <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])

ggsave(file.path(OUT, "supp_fig", "Sanity_check_genus_drivers.pdf"), fig4, width = 10, height = 6)

message("✓ Wrote figures to outputs/supp_fig/ and tables to outputs/stats/")
