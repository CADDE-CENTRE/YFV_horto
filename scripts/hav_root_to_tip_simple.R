# scripts/hav_root_to_tip_simple.R
# Root-to-tip regression plot for HAV (TempEst output table)

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

# --- I/O (simple-repo defaults) ---
if (!exists("DATA_PHY")) DATA_PHY <- file.path("data", "phylo", "trees")
if (!exists("OUT"))      OUT      <- "outputs"
infile <- file.path(DATA_PHY, "HAV_tempest.csv")

stopifnot(file.exists(infile))
dir.create(file.path(OUT, "supp_fig"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "stats"),    showWarnings = FALSE, recursive = TRUE)

# --- Flexible reader (handles varied column names) ---
pick_first <- function(nm, choices) {
  hit <- intersect(tolower(choices), tolower(nm))
  if (length(hit)) nm[match(hit[1], tolower(nm))] else NA_character_
}

raw <- readr::read_csv(infile, show_col_types = FALSE)
nm  <- names(raw)

date_col <- pick_first(nm, c("date","decimal_date","decimalyear","decimal_year","year"))
dist_col <- pick_first(nm, c("distance","root_to_tip","root2tip","root_to_tip_distance","divergence"))
res_col  <- pick_first(nm, c("residual","resid"))
loc_col  <- pick_first(nm, c("location","region","country","state","area"))

if (any(is.na(c(date_col, dist_col)))) {
  stop("Couldn't find required columns. Need decimal date and root-to-tip distance.\n",
       "Found names: ", paste(nm, collapse = ", "))
}

hepA <- raw %>%
  transmute(
    tip       = if ("tip" %in% nm) .data[["tip"]] else NA_character_,
    date      = suppressWarnings(as.numeric(.data[[date_col]])),
    distance  = suppressWarnings(as.numeric(.data[[dist_col]])),
    residual  = if (!is.na(res_col)) suppressWarnings(as.numeric(.data[[res_col]])) else NA_real_,
    location  = if (!is.na(loc_col)) as.character(.data[[loc_col]]) else "unknown"
  )

# --- Optional: inspect rows that would be excluded from plotting limits ---
bad <- subset(hepA, is.na(date) | is.na(distance) | distance < 0 | date < 1960 | date > 2025)
if (nrow(bad) > 0) print(bad[, c("tip","date","distance","location")])

# --- Clean for plotting (keep full range; we zoom with coord_cartesian) ---
hepA_clean <- subset(hepA, !is.na(date) & !is.na(distance) & distance >= 0)

# --- Fit simple LM for annotation / export ---
lm_fit <- lm(distance ~ date, data = hepA_clean)
lm_sum <- summary(lm_fit)
stats_out <- tibble::tibble(
  n            = nrow(hepA_clean),
  intercept    = unname(coef(lm_fit)[1]),
  slope        = unname(coef(lm_fit)[2]),
  r_squared    = unname(lm_sum$r.squared),
  adj_r_squared= unname(lm_sum$adj.r.squared),
  p_value      = unname(coef(lm_sum)[2, "Pr(>|t|)"])
)
#readr::write_csv(stats_out, file.path(OUT, "stats", "hav_tempest_root_to_tip_lm.csv"))

# --- Plot ---
ymax <- max(hepA_clean$distance, na.rm = TRUE) * 1.05

p <- ggplot(hepA_clean, aes(x = date, y = distance, colour = location)) +
  geom_point(alpha = 0.25, size = 3, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.6, na.rm = TRUE) +
  scale_x_continuous(breaks = seq(1960, 2025, by = 10)) +
  coord_cartesian(xlim = c(1960, 2025), ylim = c(0, ymax)) +
  labs(x = "Time in years", y = "Root-to-tip divergence", colour = "Region") +
  theme_minimal(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

ggsave(file.path(OUT, "supp_fig", "hav_tempest_root_to_tip.pdf"),
       p, width = 6.8, height = 4.6)

message("✓ wrote: ", file.path(OUT, "supp_fig", "hav_tempest_root_to_tip.(pdf)"))
