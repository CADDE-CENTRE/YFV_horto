# scripts/panel_d_np_decomposition_simple.R
# Panel 
# Fig 2d: NP decomposition stage vs genome coverage (SMART-9N ≥10×)
# - Groups on x by stage, ordering samples within each stage from highest to lowest coverage.

source("scripts/00_setup_simple.R")  # sets DATA_SUPP, OUT and loads pkgs

# ---- I/O ----
infile_std <- file.path(DATA_SUPP, "SupplementaryData4.csv")
infile_alt <- file.path(DATA_SUPP, "SupplementaryData4_finalissimo_98samples.csv")
infile <- if (file.exists(infile_std)) infile_std else infile_alt
stopifnot(file.exists(infile))

# ---- Load & standardise ----
dat <- readr::read_csv(infile, show_col_types = FALSE)
names(dat) <- tolower(names(dat))

need <- c("host", "decomposition_simple", "smart9n_cov10x_pct")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "))

# ---- NP only + stage labels ----
np <- dat %>%
  dplyr::filter(tolower(host) == "neotropicalprimate") %>%
  dplyr::mutate(
    stage = dplyr::case_when(
      decomposition_simple == "intact_live_2d"     ~ "Intact/Alive",
      decomposition_simple == "medium_or_advanced" ~ "Medium/Advanced",
      TRUE                                         ~ "Not available"
    ),
    cov10x = suppressWarnings(as.numeric(smart9n_cov10x_pct))
  ) %>%
  dplyr::filter(is.finite(cov10x))

stage_levels <- c("Intact/Alive", "Medium/Advanced", "Not available")
np$stage <- factor(np$stage, levels = stage_levels)

# ---- Order samples within each stage by coverage (desc) and place stages sequentially ----
# counts per stage (including zeros for missing stages)
stage_counts <- np %>%
  dplyr::count(stage, name = "n") %>%
  tidyr::complete(stage = factor(stage_levels, levels = stage_levels), fill = list(n = 0)) %>%
  dplyr::arrange(stage)

# block starts (no gap between blocks; set 'gap' > 0 if desired)
gap <- 0L
block_start <- c(0L,
                 stage_counts$n[1] + gap,
                 stage_counts$n[1] + stage_counts$n[2] + 2L * gap)
block_map <- tibble::tibble(stage = factor(stage_levels, levels = stage_levels),
                            block_start = block_start)

np_ord <- np %>%
  dplyr::arrange(stage, dplyr::desc(cov10x)) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(idx_in_group = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(block_map, by = "stage") %>%
  dplyr::mutate(sample_idx = idx_in_group + block_start) %>%
  dplyr::arrange(sample_idx)

# ---- Plot ----
pal <- c("Intact/Alive"    = "#40627A",
         "Medium/Advanced" = "#A7BBD1",
         "Not available"   = "#D9D9D9")

p <- ggplot(np_ord, aes(x = sample_idx, y = cov10x, fill = stage)) +
  geom_col(width = 0.9, color = "#333333", linewidth = 0.2) +
  geom_hline(yintercept = 70, linetype = "dashed", linewidth = 0.6, color = "#D55E00") +
  scale_fill_manual(values = pal, name = "NP decomposition stage:") +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Sample NPs (ordered within stage)", y = "YFV genome coverage ≥10× (%)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# ---- Save figure ----
out_fig_dir  <- file.path(OUT, "fig2")
out_stat_dir <- file.path(OUT, "stats")
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_stat_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_fig_dir, "panel_d_np_decomposition_coverage.pdf"),
       p, width = 7.2, height = 3.8, device = cairo_pdf)

# ---- Save small stats table (counts + medians/range + %>70) ----
stats_tbl <- np %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(
    n_samples   = dplyr::n(),
    median_cov  = stats::median(cov10x, na.rm = TRUE),
    min_cov     = min(cov10x, na.rm = TRUE),
    max_cov     = max(cov10x, na.rm = TRUE),
    pct_over70  = mean(cov10x > 70, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(stage, levels = stage_levels))

readr::write_csv(stats_tbl, file.path(out_stat_dir, "panel_d_np_decomposition_coverage_stats.csv"))

message("✓ wrote: ",
        file.path(out_fig_dir, "panel_d_np_decomposition_coverage.(png|pdf)"),
        " and ", file.path(out_stat_dir, "panel_d_np_decomposition_coverage_stats.csv"))
