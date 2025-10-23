# scripts/fig_ct_rpm_ridges_simple.R
# Panels:
# Fig2 (a) Ct vs RPM(log10) scatter + single global LM
# Fig2 (b) Ridgelines for Ct, RPM(log10), and N50 by Group ("I (MO)", "II (NP)", "III (NP)")

source("scripts/00_setup_simple.R")  # sets DATA_SUPP, OUT and loads pkgs

# ---- I/O ----
infile_std <- file.path(DATA_SUPP, "SupplementaryData4.csv")
#infile_alt <- file.path(DATA_SUPP, "SupplementaryData4_finalissimo_98samples.csv")
infile <- if (file.exists(infile_std)) infile_std else infile_alt
stopifnot(file.exists(infile))

# ---- Load & standardise ----
dat <- readr::read_csv(infile, show_col_types = FALSE)
names(dat) <- tolower(names(dat))

pick_first <- function(nm, choices) {
  ch <- intersect(choices, nm)
  if (length(ch)) ch[1] else NA_character_
}

ct_col   <- pick_first(names(dat), c("ct"))
rpm_col  <- pick_first(names(dat), c("rpm","rpm_log10","viral_load_rpmillion","log10_rpm","yfv_rpm_log10"))
n50_col  <- pick_first(names(dat), c("n50","smart9n_n50","smart_9n_n50","smart-9n_n50"))
host_col <- pick_first(names(dat), c("host","name_host"))
grp_col  <- pick_first(names(dat), c("group","summary_metadata"))

if (is.na(ct_col))  stop("Couldn't find a 'ct' column.")
if (is.na(rpm_col)) stop("Couldn't find a log10 RPM column (e.g., 'rpm').")

dat <- dat %>%
  mutate(
    ct   = suppressWarnings(as.numeric(.data[[ct_col]])),
    rpm  = suppressWarnings(as.numeric(.data[[rpm_col]])),  # already log10
    n50  = if (!is.na(n50_col)) suppressWarnings(as.numeric(.data[[n50_col]])) else NA_real_,
    host_raw = tolower(as.character(if (!is.na(host_col)) .data[[host_col]] else NA_character_)),
    grp_raw  = tolower(as.character(if (!is.na(grp_col))  .data[[grp_col]]  else NA_character_))
  )

# ---- Group labels: "I (MO)", "II (NP)", "III (NP)" ----
grp_num <- dplyr::case_when(
  stringr::str_detect(dat$grp_raw, "group[_ -]*i\\b")   ~ "I",
  stringr::str_detect(dat$grp_raw, "group[_ -]*ii\\b")  ~ "II",
  stringr::str_detect(dat$grp_raw, "group[_ -]*iii\\b") ~ "III",
  TRUE ~ NA_character_
)
host_tag <- dplyr::case_when(
  stringr::str_detect(dat$host_raw, "mosq") ~ "MO",
  stringr::str_detect(dat$host_raw, "neotropical") | stringr::str_detect(dat$host_raw, "primate") ~ "NP",
  TRUE ~ NA_character_
)
# Fallback if host missing
host_tag <- ifelse(is.na(host_tag) & grp_num == "I", "MO", host_tag)
host_tag <- ifelse(is.na(host_tag) & grp_num %in% c("II","III"), "NP", host_tag)

dat$group_lab <- ifelse(!is.na(grp_num) & !is.na(host_tag),
                        paste0(grp_num, " (", host_tag, ")"),
                        NA_character_)
plot_df <- dat %>% filter(is.finite(ct), is.finite(rpm), !is.na(group_lab))
plot_df$group_lab <- factor(plot_df$group_lab, levels = c("I (MO)","II (NP)","III (NP)"))

# ---- Stats (global Spearman) ----
ct_rpm_spear <- suppressWarnings(cor.test(plot_df$ct, plot_df$rpm, method = "spearman", exact = FALSE))
rho  <- unname(ct_rpm_spear$estimate)
pval <- ct_rpm_spear$p.value
fmt_r <- function(x) formatC(x, format = "f", digits = 2)
fmt_p <- function(p) if (p < 1e-3) sub("e(-?)(\\d+)", " x 10^\\1\\2", format(p, digits = 2, scientific = TRUE)) else formatC(p, format = "f", digits = 3)

# ---- Colours ----
pal <- c("I (MO)" = "#e76f51", "II (NP)" = "#8ecae6", "III (NP)" = "#2a6f9e")

# ---- Panel (a): scatter + global LM ----
pa <- ggplot(plot_df, aes(x = ct, y = rpm)) +
  geom_point(aes(color = group_lab), size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "black") +
  scale_color_manual(values = pal, name = "Group:") +
  labs(x = "Ct", y = expression(paste("Viral load (log"[10], " RPM)"))) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.6,
           label = glue::glue("Spearman's rho = {fmt_r(rho)}\n",
                              "p-value = {fmt_p(pval)}"),
           size = 3.5) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

# ---- Panel (b): ridgelines ----
ridge_theme <- ggridges::theme_ridges(font_size = 12, grid = TRUE, center_axis_labels = FALSE) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

pb_ct <- ggplot(plot_df, aes(x = ct, y = group_lab, fill = group_lab)) +
  ggridges::stat_density_ridges(geom = "density_ridges", scale = 1.1, rel_min_height = 0.01,
                                jittered_points = TRUE, point_shape = "|", point_size = 1.2, alpha = 0.9) +
  scale_fill_manual(values = pal) + labs(x = "Ct", y = NULL) + ridge_theme

pb_rpm <- ggplot(plot_df, aes(x = rpm, y = group_lab, fill = group_lab)) +
  ggridges::stat_density_ridges(geom = "density_ridges", scale = 1.1, rel_min_height = 0.01,
                                jittered_points = TRUE, point_shape = "|", point_size = 1.2, alpha = 0.9) +
  scale_fill_manual(values = pal) + labs(x = "RPM (log10)", y = NULL) + ridge_theme

if (any(is.finite(plot_df$n50))) {
  plot_df_n50 <- plot_df %>% filter(is.finite(n50))
  pb_n50 <- ggplot(plot_df_n50, aes(x = n50, y = group_lab, fill = group_lab)) +
    ggridges::stat_density_ridges(geom = "density_ridges", scale = 1.1, rel_min_height = 0.01,
                                  jittered_points = TRUE, point_shape = "|", point_size = 1.2, alpha = 0.9) +
    scale_fill_manual(values = pal) + labs(x = "N50", y = NULL) + ridge_theme
  pb <- pb_ct / pb_rpm / pb_n50
} else {
  pb <- pb_ct / pb_rpm
}

# ---- Assemble & save ----
fig <- pa + pb + patchwork::plot_layout(widths = c(2, 1)) + patchwork::plot_annotation(tag_levels = "a")

out_fig_dir  <- file.path(OUT, "fig2")
out_stat_dir <- file.path(OUT, "stats")
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_stat_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_fig_dir, "fig_ct_rpm_ridges.pdf"),  fig, width = 9, height = 4.7, device = cairo_pdf)

readr::write_csv(
  tibble::tibble(n = nrow(plot_df), spearman_r = as.numeric(rho), p_value = as.numeric(pval)),
  file.path(out_stat_dir, "fig_ct_rpm_ridges_stats.csv")
)

message("✓ wrote: ",
        file.path(out_fig_dir, "fig_ct_rpm_ridges.(png|pdf)"),
        " and ", file.path(out_stat_dir, "fig_ct_rpm_ridges_stats.csv"))
