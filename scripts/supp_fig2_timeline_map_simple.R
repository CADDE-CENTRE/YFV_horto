# scripts/supp_fig2_timeline_map_simple.R
# Supplementary Figure 2: Timeline (stacked bars) + Map (pie charts)
# Inputs: data/supplementary/SupplementaryData1.csv
# Outputs:
#   - outputs/supp_fig/supp_fig2_timeline_map.(pdf)
#   - outputs/stats/supp_fig2_timeline_by_genus.csv
#   - outputs/stats/supp_fig2_site_pies.csv

source("scripts/00_setup_simple.R")  # provides DATA_SUPP, OUT; loads tidyverse/ggplot2 if you set it there

# scripts/supp_fig2_timeline_map_simple.R
# Supplementary Figure 2: Timeline (stacked bars) + Map (pie charts)

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(scatterpie)
  library(patchwork)
})

# -------------------------------
# 1) Load data (from simple-repo paths)
# -------------------------------
cand <- c(
  file.path(DATA_SUPP, "SupplementaryData1.csv"),
  file.path(DATA_SUPP, "Supplementary Data 1.csv")
)
infile <- cand[file.exists(cand)][1]
stopifnot(!is.na(infile))

dir.create(file.path(OUT, "supp_fig"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "stats"),    showWarnings = FALSE, recursive = TRUE)

df <- read_csv(infile, show_col_types = FALSE)

# Coerce types (keep original column names/casing)
df <- df %>%
  mutate(
    Date     = as.Date(Date),
    Quantity = suppressWarnings(as.numeric(Quantity)),
    Genus    = as.factor(Genus),
    Latitude = suppressWarnings(as.numeric(Latitude)),
    Longitude= suppressWarnings(as.numeric(Longitude))
  ) %>%
  # Safety: drop rows without coords/date/quantity
  filter(!is.na(Date), is.finite(Latitude), is.finite(Longitude), is.finite(Quantity))

# -------------------------------
# 2) Shared colours (consistent across plots)
# -------------------------------
genus_levels <- sort(unique(df$Genus))
base_pal <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
colours  <- setNames(colorRampPalette(base_pal)(length(genus_levels)), genus_levels)

# -------------------------------
# 3) Plot 1: Stacked bar timeline (by Genus)
# -------------------------------
timeline_genus <- df %>%
  group_by(Date, Genus) %>%
  summarise(Total = sum(Quantity, na.rm = TRUE), .groups = "drop") %>%
  arrange(Date, Genus)

p1 <- ggplot(timeline_genus, aes(x = Date, y = Total, fill = Genus)) +
  geom_col() +
  labs(x = "Date of collection", y = "Number of mosquitoes collected", fill = "Genus") +
  scale_fill_manual(values = colours, drop = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8),
    plot.title      = element_blank()
  )

# -------------------------------
# 4) Plot 2: Map with pie charts at sampling sites
# -------------------------------
site_genus <- df %>%
  group_by(Latitude, Longitude, Genus) %>%
  summarise(Total = sum(Quantity, na.rm = TRUE), .groups = "drop")

site_wide <- site_genus %>%
  tidyr::pivot_wider(
    names_from  = Genus,
    values_from = Total,
    values_fill = 0
  ) %>%
  # total per site for radius scaling
  mutate(Total_sum = rowSums(across(all_of(as.character(genus_levels))), na.rm = TRUE))

# --- Dynamic radius: make the 95th percentile site ≈ 0.15° ---
q95 <- stats::quantile(site_wide$Total_sum[site_wide$Total_sum > 0], probs = 0.95, na.rm = TRUE)
target_deg <- 0.0005
scale_r <- if (is.finite(q95) && q95 > 0) target_deg / sqrt(q95) else 0.01
# (Optional) hard-cap largest pies at target radius:
# site_wide <- site_wide %>% mutate(radius_deg = pmin(sqrt(Total_sum) * scale_r, target_deg))
site_wide <- site_wide %>% mutate(radius_deg = sqrt(Total_sum) * scale_r)

p2 <- ggplot() +
  scatterpie::geom_scatterpie(
    data = site_wide,
    aes(
      x = Longitude,
      y = Latitude,
      r = radius_deg
    ),
    cols = as.character(genus_levels)
  ) +
  coord_equal() +
  scale_fill_manual(values = colours, name = "Genus", drop = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_blank())

# -------------------------------
# 5) Combine & save
# -------------------------------
combined <- p1 / p2 + patchwork::plot_layout(heights = c(2, 1))

ggsave(file.path(OUT, "supp_fig", "supp_fig2_timeline_map.pdf"),
       combined, width = 10, height = 10, units = "in")


# -------------------------------
# 6) Save data used
# -------------------------------
#readr::write_csv(timeline_genus, file.path(OUT, "stats", "supp_fig2_timeline_by_genus.csv"))

message("✓ wrote: ",
        file.path(OUT, "supp_fig", "supp_fig2_timeline_map.(pdf)"),
        " and stats CSVs in outputs/stats/")
