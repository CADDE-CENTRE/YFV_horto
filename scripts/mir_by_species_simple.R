# scripts/mir_by_species_from_suppl_noagg.R
# MIR by species from SupplementaryData3.csv, tested_pools >= 10, NO aggregation.

source("scripts/00_setup_simple.R")  # provides DATA_SUPP, OUT and loads tidyverse/ggplot2

suppressPackageStartupMessages({
  library(ggpubr)
})

# --- I/O ---
infile <- file.path(DATA_SUPP, "SupplementaryData3.csv")
stopifnot(file.exists(infile))
dir.create(file.path(OUT, "supp_fig"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "stats"),    showWarnings = FALSE, recursive = TRUE)

# --- Load & standardise ---
dat <- readr::read_csv(infile, show_col_types = FALSE)
names(dat) <- tolower(names(dat))

need <- c("species","municipality","mir","tested_pools")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "))

df <- dat %>%
  mutate(
    include_flag = if ("include" %in% names(.)) tolower(include) %in% c("yes","y","1") else TRUE,
    species_raw  = tolower(species),
    MIR          = suppressWarnings(as.numeric(mir)),
    tested_pools = suppressWarnings(as.numeric(tested_pools))
  ) %>%
  filter(include_flag, is.finite(MIR), is.finite(tested_pools), tested_pools >= 10) %>%
  filter(species_raw %in% c("leucocelaenus","janthinomys","other")) %>%
  mutate(
    species_grp = dplyr::recode(species_raw, "other" = "Other"),
    species_grp = factor(species_grp, levels = c("leucocelaenus","janthinomys","Other"))
  )

stopifnot(nrow(df) > 0, dplyr::n_distinct(df$species_grp) >= 2)

# --- Stats (no aggregation; one row per record) ---
kruskal_test <- ggpubr::compare_means(MIR ~ species_grp, data = df, method = "kruskal.test")
pairwise     <- ggpubr::compare_means(MIR ~ species_grp, data = df,
                                      method = "wilcox.test", p.adjust.method = "fdr",
                                      exact = FALSE)

readr::write_csv(
  dplyr::bind_rows(kruskal_test %>% dplyr::mutate(test = "kruskal"),
                   pairwise %>% dplyr::mutate(test = "wilcoxon_pairwise")),
  file.path(OUT, "stats", "MIR_by_species_tests_noagg.csv")
)
readr::write_csv(df, file.path(OUT, "stats", "MIR_by_species_dataset_used_noagg.csv"))

# --- Write stats CSVs ---
dir.create(file.path(OUT, "stats"), showWarnings = FALSE, recursive = TRUE)

# inferential tests
tests_tbl <- dplyr::bind_rows(
  kruskal_test %>% dplyr::mutate(test = "kruskal"),
  pairwise %>% dplyr::mutate(test = "wilcoxon_pairwise")
)
readr::write_csv(tests_tbl, file.path(OUT, "stats", "MIR_by_species_tests.csv"))

# descriptive summary by species group
summary_tbl <- df %>%
  dplyr::group_by(species_grp) %>%
  dplyr::summarise(
    n           = dplyr::n(),
    median_MIR  = median(MIR, na.rm = TRUE),
    min_MIR     = min(MIR, na.rm = TRUE),
    max_MIR     = max(MIR, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(summary_tbl, file.path(OUT, "stats", "MIR_by_species_group_summary.csv"))

# pairwise p-values as a wide matrix
pw_mat <- pairwise.wilcox.test(df$MIR, df$species_grp, p.adjust.method = "fdr", exact = FALSE)$p.value
pw_long <- as.data.frame(as.table(pw_mat))
names(pw_long) <- c("group1","group2","p_adj")
readr::write_csv(pw_long, file.path(OUT, "stats", "MIR_by_species_pairwise_matrix.csv"))


# --- Plot ---
my_comparisons <- list(
  c("leucocelaenus","janthinomys"),
  c("leucocelaenus","Other"),
  c("janthinomys","Other")
)

p <- ggpubr::ggboxplot(
  df, x = "species_grp", y = "MIR",
  color = "species_grp", palette = "jco",
  add = "jitter", add.params = list(size = 2, alpha = 0.7)
) +
  ggpubr::stat_compare_means(
    method = "kruskal.test",
    label.y = max(df$MIR, na.rm = TRUE) + 0.05 * max(df$MIR, na.rm = TRUE),
    label.x.npc = "centre"
  ) +
  ggpubr::stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    label = "p.signif",
    step.increase = 0.07
  ) +
  labs(x = NULL, y = "MIR (per 1,000 mosquitoes)",
       title = "Minimum Infection Rate by Species (tested_pools ≥ 10)") +
  ggpubr::theme_pubr()

ggplot2::ggsave(file.path(OUT, "supp_fig", "MIR_by_species_boxplot.pdf"),
                p, width = 6.5, height = 4.2, device = cairo_pdf)

message("✓ wrote stats tables to outputs/stats/ and figures to outputs/supp_fig/")

