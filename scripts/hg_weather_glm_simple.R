# scripts/hg_weather_glm_simple.R
# Haemagogus abundance vs meteorological variables (monthly); negative binomial GLMs
# Inputs (searched in order): 
#   data/climate/SP_met_and_haemagogus_data_MONTHLY.csv
#   data/climate/SP_met_and_haemagogus_data_MONTHLY.CSV
#   data/climate/SP_met.csv
# Outputs:
#   - outputs/stats/lagged.html
#   - outputs/stats/contemporaneous.html
#   - outputs/stats/contemporaneous_funcs.html
#   - outputs/fig1/correlation_weather.pdf
#   - outputs/fig1/variables_time.pdf

source("scripts/00_setup_simple.R")  # sets DATA_CLIM, OUT, loads tidyverse/ggplot2/etc.

suppressPackageStartupMessages({
  library(stargazer)
  library(MASS)        # glm.nb
  library(corrplot)
  library(lubridate)
  library(forcats)
})

# ---------- resolve input ----------
cand <- c(
  file.path(DATA_CLIM, "SP_met_and_haemagogus_data_MONTHLY.csv"),
  file.path(DATA_CLIM, "SP_met.csv")
)
infile <- cand[file.exists(cand)][1]
if (is.na(infile)) stop("Could not find input in data/climate/. Tried:\n  - ", paste(cand, collapse = "\n  - "))

# ---------- read & harmonise ----------
df <- readr::read_csv(infile, show_col_types = FALSE) |>
  rename_with(tolower) |>
  # expected original columns: haemagogus_mosquitoes, cum_rain, min_temp, max_temp, mean_temp, min_hum, max_hum, mean_hum, year, month_original
  rename(n = haemagogus_mosquitoes) |>
  mutate(
    cum_rain = cum_rain * 100,  # to cm
    `rain (cm)`        = cum_rain,
    `min temp. (C)`    = min_temp,
    `max temp. (C)`    = max_temp,
    `mean temp. (C)`   = mean_temp,
    `min humidity (%)` = min_hum,
    `max humidity (%)` = max_hum,
    `mean humidity (%)`= mean_hum
  ) |>
  mutate(
    # build a mid-month date (assumes month_original = 1..12)
    date = make_date(year, month_original, 15)
  ) |>
  arrange(year, month_original)

# quick sanity on required columns
req <- c("n","rain (cm)","min temp. (C)","max temp. (C)","mean temp. (C)","min humidity (%)","max humidity (%)","mean humidity (%)")
miss <- setdiff(req, names(df))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

# ---------- helpers ----------
pseudo_r2 <- function(model) {
  r2 <- 1 - model$deviance / model$null.deviance
  format(round(r2, 3), nsmall = 3)
}
pseudo_r2_list_models <- function(models) vapply(models, pseudo_r2, character(1))
output_stargazer <- function(models, output_name, ...) {
  r2 <- pseudo_r2_list_models(models)
  stargazer(models,
            out = file.path(OUT, "stats", output_name),
            add.lines = list(c("Pseudo-R²", r2)),
            dep.var.caption = "Dependent variable: Mosquitoes captured",
            dep.var.labels = "",
            apply.coef = exp,         # multiplicative interpretation
            omit = "Constant",
            ci = TRUE,
            t.auto = FALSE, p.auto = FALSE,
            type = "html",
            ...
  )
}

# ensure output dirs
dir.create(file.path(OUT, "stats"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "fig1c"),  showWarnings = FALSE, recursive = TRUE)

# ---------- LAGGED MODELS (lag = 1 month) ----------
model_r       <- glm.nb(n ~ dplyr::lag(`rain (cm)`, 1),        data = df)
model_temp_1  <- glm.nb(n ~ dplyr::lag(`min temp. (C)`, 1),    data = df)
model_temp_2  <- glm.nb(n ~ dplyr::lag(`mean temp. (C)`, 1),   data = df)
model_temp_3  <- glm.nb(n ~ dplyr::lag(`max temp. (C)`, 1),    data = df)
model_hum_1   <- glm.nb(n ~ dplyr::lag(`min humidity (%)`, 1), data = df)
model_hum_2   <- glm.nb(n ~ dplyr::lag(`mean humidity (%)`, 1),data = df)
model_hum_3   <- glm.nb(n ~ dplyr::lag(`max humidity (%)`, 1), data = df)

model_all     <- glm.nb(n ~ dplyr::lag(`min humidity (%)`, 1) + dplyr::lag(`mean temp. (C)`, 1) + dplyr::lag(`rain (cm)`, 1), data = df)
model_all_1   <- glm.nb(n ~ dplyr::lag(`mean temp. (C)`,1) + dplyr::lag(`rain (cm)`,1), data = df)
model_all_2   <- glm.nb(n ~ dplyr::lag(`min temp. (C)`,1)  + dplyr::lag(`rain (cm)`,1), data = df)
model_all_3   <- glm.nb(n ~ dplyr::lag(`max temp. (C)`, 1) + dplyr::lag(`rain (cm)`, 1), data = df)
model_all_1a  <- glm.nb(n ~ dplyr::lag(`mean temp. (C)`, 1) + dplyr::lag(`min humidity (%)`, 1), data = df)
model_all_2a  <- glm.nb(n ~ dplyr::lag(`mean temp. (C)`, 1) + dplyr::lag(`mean humidity (%)`, 1), data = df)
model_all_3a  <- glm.nb(n ~ dplyr::lag(`mean temp. (C)`, 1) + dplyr::lag(`max humidity (%)`, 1), data = df)

models_lag <- list(
  model_r, model_temp_1, model_temp_2, model_temp_3,
  model_hum_1, model_hum_2, model_hum_3,
  model_all, model_all_1, model_all_2, model_all_3,
  model_all_1a, model_all_2a, model_all_3a
)
output_stargazer(models_lag, "lagged.html")

# ---------- CONTEMPORANEOUS MODELS ----------
model_r       <- glm.nb(n ~ `rain (cm)`,        data = df)
model_temp_1  <- glm.nb(n ~ `min temp. (C)`,    data = df)
model_temp_2  <- glm.nb(n ~ `mean temp. (C)`,   data = df)
model_temp_3  <- glm.nb(n ~ `max temp. (C)`,    data = df)
model_hum_1   <- glm.nb(n ~ `min humidity (%)`, data = df)
model_hum_2   <- glm.nb(n ~ `mean humidity (%)`,data = df)
model_hum_3   <- glm.nb(n ~ `max humidity (%)`, data = df)

model_all     <- glm.nb(n ~ `mean temp. (C)` + `min humidity (%)` + `rain (cm)`, data = df)
model_all_1   <- glm.nb(n ~ `mean temp. (C)` + `rain (cm)`,       data = df)
model_all_2   <- glm.nb(n ~ `min temp. (C)`  + `rain (cm)`,       data = df)
model_all_3   <- glm.nb(n ~ `max temp. (C)`  + `rain (cm)`,       data = df)
model_all_1a  <- glm.nb(n ~ `mean temp. (C)` + `mean humidity (%)`,data = df)
model_all_2a  <- glm.nb(n ~ `mean temp. (C)` + `min humidity (%)`, data = df)
model_all_3a  <- glm.nb(n ~ `mean temp. (C)` + `max humidity (%)`, data = df)

models_contemp <- list(
  model_r, model_temp_1, model_temp_2, model_temp_3,
  model_hum_1, model_hum_2, model_hum_3,
  model_all, model_all_1, model_all_2, model_all_3,
  model_all_1a, model_all_2a, model_all_3a
)
output_stargazer(models_contemp, "contemporaneous.html")

# ANOVA on the 3-var contemporaneous model
anova(model_all)

# ---------- OTHER FUNCTIONAL FORMS ----------
# guard logs against non-positive values
rain_pos <- pmax(df$`rain (cm)`, .Machine$double.eps)
mint_pos <- pmax(df$`min temp. (C)`, .Machine$double.eps)
meant_pos<- pmax(df$`mean temp. (C)`, .Machine$double.eps)
maxt_pos <- pmax(df$`max temp. (C)`, .Machine$double.eps)

model_all_1  <- glm.nb(n ~ `mean temp. (C)` + `rain (cm)`, data = df)
model_all_2  <- glm.nb(n ~ `min temp. (C)`  + `rain (cm)`, data = df)
model_all_3  <- glm.nb(n ~ `max temp. (C)`  + `rain (cm)`, data = df)

model_all_1a <- glm.nb(n ~ log(meant_pos) + log(rain_pos), data = df)
model_all_2a <- glm.nb(n ~ log(mint_pos)  + log(rain_pos), data = df)
model_all_3a <- glm.nb(n ~ log(maxt_pos)  + log(rain_pos), data = df)

model_all_1b <- glm.nb(n ~ `mean temp. (C)` + I(`mean temp. (C)`^2) + `rain (cm)` + I(`rain (cm)`^2), data = df)
model_all_2b <- glm.nb(n ~ `min temp. (C)`  + I(`min temp. (C)`^2)  + `rain (cm)` + I(`rain (cm)`^2), data = df)
model_all_3b <- glm.nb(n ~ `max temp. (C)`  + I(`max temp. (C)`^2)  + `rain (cm)` + I(`rain (cm)`^2), data = df)

model_all_1c <- glm.nb(n ~ `mean temp. (C)` + `rain (cm)` + `mean temp. (C)`:`rain (cm)`, data = df)
model_all_2c <- glm.nb(n ~ `min temp. (C)`  + `rain (cm)` + `min temp. (C)`:`rain (cm)`,  data = df)
model_all_3c <- glm.nb(n ~ `max temp. (C)`  + `rain (cm)` + `max temp. (C)`:`rain (cm)`,  data = df)

models_funcs <- list(
  model_all_1, model_all_2, model_all_3,
  model_all_1a, model_all_2a, model_all_3a,
  model_all_1b, model_all_2b, model_all_3b,
  model_all_1c, model_all_2c, model_all_3c
)
output_stargazer(models_funcs, "contemporaneous_funcs.html")

# ---------- Correlation plot of weather vars (robust; includes rain) ----------
# Works whether your columns are pretty-labeled (e.g., "rain (cm)") or raw (e.g., cum_rain)

pick <- function(df, choices) {
  hit <- intersect(tolower(choices), tolower(names(df)))
  if (length(hit)) names(df)[match(hit[1], tolower(names(df)))] else NA_character_
}

vars <- c(
  min_temp = pick(df, c("min temp. (c)", "min_temp")),
  mean_temp= pick(df, c("mean temp. (c)", "mean_temp")),
  max_temp = pick(df, c("max temp. (c)", "max_temp")),
  min_hum  = pick(df, c("min humidity (%)", "min_hum")),
  mean_hum = pick(df, c("mean humidity (%)", "mean_hum")),
  max_hum  = pick(df, c("max humidity (%)", "max_hum")),
  rain_cm  = pick(df, c("rain (cm)", "cum_rain", "rain_cm", "rain"))
)

# keep only those found
vars <- vars[!is.na(vars)]

x <- df[, unname(vars), drop = FALSE]
colnames(x) <- names(vars)

m_cor <- stats::cor(x, use = "pairwise.complete.obs", method = "pearson")

pdf(file.path(OUT, "supp_fig", "correlation_weather.pdf"))
corrplot::corrplot(
  m_cor,
  method   = "circle",
  type     = "lower",
  addCoef.col = "black",
  tl.col   = "black",
  tl.srt   = 45,
  number.cex = 0.7
)
dev.off()


# ---------- Time-variation plot (robust) ----------
lookup <- tibble::tribble(
  ~name,              ~category,
  "rain (cm)",        "Rain",
  "min temp. (C)",    "Temperature",
  "max temp. (C)",    "Temperature",
  "mean temp. (C)",   "Temperature",
  "min humidity (%)", "Humidity",
  "max humidity (%)", "Humidity",
  "mean humidity (%)","Humidity"
)

# pick variables by name (avoid range selection)
var_names <- c(
  "min temp. (C)", "mean temp. (C)", "max temp. (C)",
  "min humidity (%)", "mean humidity (%)", "max humidity (%)",
  "rain (cm)"
)

df_short <- df |>
  dplyr::select(date, n, dplyr::any_of(var_names)) |>
  tidyr::pivot_longer(cols = dplyr::any_of(var_names),
                      names_to = "name", values_to = "value") |>
  dplyr::left_join(lookup, by = "name")

# fill NA counts to avoid dropped bars
df_short <- df_short |>
  dplyr::mutate(n_fill = tidyr::replace_na(n, 0))

# compute present categories & names for clean releveling
cats_present  <- unique(df_short$category)
names_present <- unique(df_short$name)

g <- df_short |>
  dplyr::group_by(name) |>
  dplyr::mutate(value = value / max(value, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    value    = value * max(n_fill, na.rm = TRUE),
    category = forcats::fct_relevel(
      factor(category),
      intersect(c("Rain", "Temperature", "Humidity"), cats_present)
    ),
    name     = forcats::fct_relevel(
      factor(name),
      intersect(c(
        "rain (cm)",
        "min temp. (C)", "mean temp. (C)", "max temp. (C)",
        "min humidity (%)", "mean humidity (%)", "max humidity (%)"
      ), names_present)
    )
  ) |>
  ggplot(aes(x = date, y = value, group = name)) +
  geom_col(aes(y = n_fill), fill = "grey85") +
  geom_line(aes(colour = name)) +
  scale_color_brewer("Variable", palette = "Dark2") +
  xlab("Date") + ylab("Mosquitoes captured") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text()) +
  facet_wrap(~ category, nrow = 1, scales = "free_y")

ggsave(file.path(OUT, "supp_fig", "variables_time.pdf"), g, width = 10, height = 4)

message("✓ wrote stargazer tables to outputs/stats/ and figures to outputs/fig1/ and outputs/supp_fig/")
