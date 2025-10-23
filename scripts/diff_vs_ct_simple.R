# scripts/diff_vs_ct_simple.R
# Coverage difference (metagenomics – tiled-amplicon) vs Ct
# Uses the simple project setup from scripts/00_setup_simple.R

source("scripts/00_setup_simple.R")  # loads packages, sets DATA_SUPP and OUT, defines fmt_p, fmt_num, rng

# ---- load & standardise ----
infile <- file.path(DATA_SUPP, "SupplementaryData4.csv")
stopifnot(file.exists(infile))

dat <- readr::read_csv(infile, show_col_types = FALSE)
names(dat) <- tolower(names(dat))

# ---- required columns ----
req <- c("amplicon_cov_pct", "smart9n_cov10x_pct", "ct")
missing <- setdiff(req, names(dat))
if (length(missing)) stop("Missing required column(s): ", paste(missing, collapse = ", "))

# ---- paired subset + difference ----
paired <- dat %>%
  mutate(
    amplicon_cov_pct   = as.numeric(amplicon_cov_pct),
    smart9n_cov10x_pct = as.numeric(smart9n_cov10x_pct),
    ct                 = as.numeric(ct),
    diff_pp            = smart9n_cov10x_pct - amplicon_cov_pct
  ) %>%
  filter(is.finite(amplicon_cov_pct), is.finite(smart9n_cov10x_pct), is.finite(ct))

n_pairs <- nrow(paired)
if (n_pairs < 3) stop("Fewer than 3 paired samples with Ct available.")

# ---- Spearman correlation ----
ct_diff_cor <- suppressWarnings(cor.test(paired$ct, paired$diff_pp, method = "spearman", exact = FALSE))
rho  <- unname(ct_diff_cor$estimate)
pval <- ct_diff_cor$p.value

fmt_r <- function(x) formatC(x, format = "f", digits = 2)

# ---- plot: single global LM line; points coloured by host if present ----
has_host <- "host" %in% names(paired)

p <- ggplot(paired, aes(x = ct, y = diff_pp)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.7)

if (has_host) {
  p <- p + geom_point(aes(color = host), size = 2, alpha = 0.85, na.rm = TRUE)
} else {
  p <- p + geom_point(size = 2, alpha = 0.85, na.rm = TRUE)
}

p <- p +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "black") +
  labs(
    x = "Ct",
    y = "Coverage difference (metagenomics – tiled-amplicon, percentage points)",
    color = if (has_host) "Host" else NULL
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.4,
    label = glue::glue("Spearman r = {fmt_r(rho)}\np = {fmt_p(pval)}\nn = {n_pairs}")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = if (has_host) "right" else "none",
    panel.grid.minor = element_blank()
  )

# ---- save outputs ----
out_fig_dir  <- file.path(OUT, "supp_fig")
out_stat_dir <- file.path(OUT, "stats")
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_stat_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_fig_dir,  "diff_vs_ct_global.pdf"), p, width = 6.5, height = 4.5, device = cairo_pdf)

readr::write_csv(
  tibble::tibble(n_pairs = n_pairs, spearman_r = rho, p_value = pval),
  file.path(out_stat_dir, "diff_vs_ct_stats.csv")
)

message("✓ wrote: ",
        file.path(out_fig_dir,  "diff_vs_ct_global.(pdf|png)"),
        " and ", file.path(out_stat_dir, "diff_vs_ct_stats.csv"))
