# ==============================================================================
# Figure 1: Epidemiological Context and Surveillance (Integrated Version)
# ==============================================================================

# 1. Load Libraries
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)
library(readr)
library(tidyr)

# 2. Load Data
# Note: Ensure these files are in your working directory
data_np_b <- read.csv("TableS3_PEAL_NP_short.csv") 
df_met    <- read_csv("SP_met_and_haemagogus_data_MONTHLY.CSV", show_col_types = FALSE)
df_mosq_raw <- read_csv("SupplementaryData1.csv", show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 3. Data Preparation for Panel B (NHP Timeline)
# ------------------------------------------------------------------------------

# Clean dates
data_np_b$Date <- as.Date(data_np_b$Date_death, format = "%d/%m/%Y")

# Color mapping for bars
status_map <- c(
  "negative"                 = "#89b6d8", 
  "not_collected_decomp"     = "#f8e69a", 
  "positive_not_seq_low_cov" = "#e65539", 
  "Sequenced_high_coverage"  = "#b31b1b"  
)

# Process cumulative data for the line and dots (All but negative)
data_np_b_line <- data_np_b %>%
  filter(!is.na(Date), Zone_PEAL_Status_simple != "negative") %>%
  arrange(Date) %>%
  mutate(
    cum_n = row_number(),
    cum_pct = (cum_n / n()) * 100
  )

# Data for bars (All cases)
data_np_b_bars <- data_np_b %>%
  filter(!is.na(Date)) %>%
  mutate(
    Status = factor(Zone_PEAL_Status_simple, levels = names(status_map)),
    epiweek = floor_date(Date, unit = "week", week_start = 7)
  )

# Calculate Weekly Midpoint Breaks for x-axis ticks
breaks_b <- seq(as.Date("2017-10-01"), as.Date("2018-01-14"), by = "1 week") + 3.5

# Function to show labels only for the first week of each month
get_monthly_labels <- function(breaks) {
  labels <- format(breaks, "%b\n%Y")
  labels[duplicated(format(breaks, "%m-%Y"))] <- ""
  return(labels)
}

# ------------------------------------------------------------------------------
# 4. PLOT PANEL B
# ------------------------------------------------------------------------------

scale_factor_b <- 2 # Max 50 deaths aligns with 100%

p_panel_b <- ggplot() +
  geom_bar(data = data_np_b_bars, aes(x = epiweek + 3.5, fill = Status), 
           color = "black", linewidth = 0.2, width = 6) +
  geom_line(data = data_np_b_line, aes(x = Date, y = cum_pct / scale_factor_b), 
            color = "black", linewidth = 0.8) +
  # Circles plotted in grey for all non-negative cases
  geom_point(data = data_np_b_line, aes(x = Date, y = cum_pct / scale_factor_b), 
             shape = 21, size = 1.5, fill = "grey80", color = "black") +
  scale_fill_manual(values = status_map, name = "NP YFV result:") +
  scale_y_continuous(
    name = "No. NP deaths per epiweek",
    limits = c(0, 50),      
    breaks = seq(0, 50, 10),
    sec.axis = sec_axis(~ . * scale_factor_b, name = "Cumulative % confirmed NP deaths", 
                        breaks = seq(0, 100, 25))
  ) +
  scale_x_date(breaks = breaks_b, labels = get_monthly_labels(breaks_b)) +
  coord_cartesian(xlim = as.Date(c("2017-10-01", "2018-01-07"))) +
  theme_classic(base_size = 12) +
  labs(tag = "b", x = "") +
  theme(legend.position = "right", plot.tag = element_text(face = "bold"))

# ------------------------------------------------------------------------------
# 5. PLOT PANEL C (User Provided Script)
# ------------------------------------------------------------------------------

# Prepare climate and mosquito data
df_met <- df_met %>% mutate(month_date = as.Date(paste0(Date, "-01")))
df_mosq_raw <- df_mosq_raw %>% mutate(Date = as.Date(Date))

hg_combined <- df_mosq_raw %>%
  filter(Species == "leucocelaenus") %>%
  filter(PCR_YFV %in% c("Positive", "Negative")) %>%
  mutate(
    epiweek = floor_date(Date, unit = "week", week_start = 7),
    Status = ifelse(PCR_YFV == "Positive", "Positive and sequenced", "Negative")
  ) %>%
  mutate(Status = factor(Status, levels = c("Negative", "Positive and sequenced")))

# Alignment Logic: 0 pools = 10°C (Offset = 10)
offset <- 10
xlims_c <- as.Date(c("2017-11-01", "2018-10-15"))

p_panel_c <- ggplot() +
  # 1. Temperature Background
  geom_ribbon(data = df_met, 
              aes(x = month_date, ymin = min_temp - offset, ymax = max_temp - offset), 
              fill = "grey90", alpha = 0.8) +
  geom_line(data = df_met, 
            aes(x = month_date, y = mean_temp - offset), 
            color = "grey30", linewidth = 0.8) +
  
  # 2. Mosquito Bars (Epi Weeks Midpoint Ticks)
  geom_bar(data = hg_combined, aes(x = epiweek + 3.5, fill = Status), 
           color = "black", linewidth = 0.2, position = "stack", width = 5.5) +
  
  # 3. Scales and Labels
  scale_fill_manual(values = c("Negative" = "#89b6d8", 
                               "Positive and sequenced" = "#EB0000"),
                    name = "MO YFV result:") +
  
  scale_x_date(date_labels = "%b\n%Y", date_breaks = "2 months") +
  
  scale_y_continuous(
    name = "No. Haemagogus\nleucocelaenus pools",
    limits = c(0, 20), 
    breaks = seq(0, 20, 5),
    sec.axis = sec_axis(~ . + offset, name = "Temperature (C)", breaks = seq(10, 35, 5)) 
  ) +
  
  coord_cartesian(xlim = xlims_c) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.title.y.right = element_text(color = "grey30"),
    plot.tag = element_text(face = "bold")
  ) +
  labs(tag = "c", x = "")

# ------------------------------------------------------------------------------
# 6. Combine and Save
# ------------------------------------------------------------------------------

final_plot <- p_panel_b / p_panel_c + plot_layout(heights = c(1, 1))
ggsave("Figure1_BC_Final.pdf", final_plot, width = 10, height = 10)

# Statistics
plotted_circle_count <- nrow(data_np_b_line)
cat("--- Final Statistics ---\n")
cat("Total number of grey circles plotted in Panel B:", plotted_circle_count, "\n")