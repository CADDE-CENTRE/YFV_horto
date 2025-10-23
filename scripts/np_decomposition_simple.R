# scripts/np_decomposition_simple.R
# NP carcass decomposition (Early vs Later) — medians BEFORE any 70% filter
# Inputs: data/supplementary/SupplementaryData4.csv
# Outputs:
#   - outputs/stats/np_decomposition_medians_all.csv  (medians over ALL NP per stage)
#   - outputs/stats/np_decomposition_summary.csv      (medians + %>70 + Wilcoxon p for Ct/RPM/N50)

source("scripts/00_setup_simple.R")  # sets DATA_SUPP, OUT, fmt_p(), fmt_num(), rng()

# ---- load ----
infile <- file.path(DATA_SUPP, "SupplementaryData4.csv")
stopifnot(file.exists(infile))
dat <- readr::read_csv(infile, show_col_types = FALSE)
names(dat) <- tolower(names(dat))

# ---- required columns ----
need <- c("host", "decomposition_simple", "ct", "smart9n_cov10x_pct")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "))

# Flexible names
rpm_nm <- intersect(c("rpm","rpm_log10","viral_load_rpmillion"), names(dat))
n50_nm <- intersect(c("n50","smart9n_n50","smart_9n_n50","smart-9n_n50"), names(dat))
if (!length(rpm_nm)) stop("Missing RPM (log10) column.")
if (!length(n50_nm)) stop("Missing N50 column.")

# ---- NP only + stage mapping ----
np <- dat %>%
  mutate(
    host_lc = tolower(host),
    dec_raw = tolower(decomposition_simple)
  ) %>%
  filter(host_lc == "neotropicalprimate") %>%
  mutate(
    stage = dplyr::case_when(
      dec_raw == "intact_live_2d"     ~ "Early (intact/alive)",
      dec_raw == "medium_or_advanced" ~ "Later (medium/advanced)",
      TRUE                            ~ NA_character_
    ),
    ct     = suppressWarnings(as.numeric(ct)),
    rpm    = suppressWarnings(as.numeric(.data[[rpm_nm[1]]])),
    n50    = suppressWarnings(as.numeric(.data[[n50_nm[1]]])),
    cov10x = suppressWarnings(as.numeric(smart9n_cov10x_pct))
  )

# ---- medians BEFORE any 70% filter (NP only, by stage) ----
meds_all <- np %>%
  filter(stage %in% c("Early (intact/alive)", "Later (medium/advanced)")) %>%
  group_by(stage) %>%
  summarise(
    n_samples         = dplyr::n(),
    median_ct         = median(ct,    na.rm = TRUE),
    median_rpm_log10  = median(rpm,   na.rm = TRUE),
    median_n50        = median(n50,   na.rm = TRUE),
    median_cov10x_pct = median(cov10x,na.rm = TRUE),
    # range (optional, also pre-threshold)
    range_ct          = rng(ct,    d = 2),
    range_rpm_log10   = rng(rpm,   d = 2),
    range_n50         = rng(n50,   d = 0),
    range_cov10x_pct  = rng(cov10x,d = 1),
    .groups = "drop"
  )

# ---- % > 70% coverage (does NOT affect medians) ----
prop70 <- np %>%
  filter(stage %in% c("Early (intact/alive)", "Later (medium/advanced)"), is.finite(cov10x)) %>%
  group_by(stage) %>%
  summarise(pct_over_70 = mean(cov10x > 70, na.rm = TRUE) * 100, .groups = "drop")

# ---- Wilcoxon tests (Early vs Later) for Ct, RPM, N50 (no 70% filter) ----
wilc <- function(df, var){
  tmp <- df %>%
    filter(stage %in% c("Early (intact/alive)", "Later (medium/advanced)")) %>%
    select(stage, value = dplyr::all_of(var)) %>%
    filter(is.finite(value))
  if (nrow(tmp) < 3 || dplyr::n_distinct(tmp$stage) < 2) return(NA_real_)
  stats::wilcox.test(value ~ stage, data = tmp, exact = FALSE)$p.value
}
p_ct  <- wilc(np, "ct")
p_rpm <- wilc(np, "rpm")
p_n50 <- wilc(np, "n50")

# ---- write outputs ----
out_dir <- file.path(OUT, "stats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Medians only (pre-threshold)
readr::write_csv(meds_all, file.path(out_dir, "np_decomposition_medians_all.csv"))

# Combined summary (medians + %>70 + p-values)
summary_tbl <- meds_all %>%
  left_join(prop70, by = "stage") %>%
  mutate(
    wilcoxon_p_ct        = fmt_p(p_ct),
    wilcoxon_p_rpm_log10 = fmt_p(p_rpm),
    wilcoxon_p_n50       = fmt_p(p_n50)
  )

readr::write_csv(summary_tbl, file.path(out_dir, "np_decomposition_summary.csv"))

message("✓ wrote: ",
        file.path(out_dir, "np_decomposition_medians_all.csv"), " and ",
        file.path(out_dir, "np_decomposition_summary.csv"))
