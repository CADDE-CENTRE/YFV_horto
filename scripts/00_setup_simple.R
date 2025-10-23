# 1) scripts/00_setup_simple.R (one-time per session)
# scripts/00_setup_simple.R
# Loads packages, sets paths, creates output dirs, and defines a tiny helper.

# install once if missing
pkgs <- c("here","readr","dplyr","stringr","ggplot2","ggridges","patchwork")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, quiet = TRUE)

suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(stringr)
  library(ggplot2); library(ggridges); library(patchwork)
})

# paths (all relative to project root)
DATA_SUPP <- here::here("data","supplementary")
DATA_CLIM <- here::here("data","climate")
OUT       <- here::here("outputs")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT,"fig2"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT,"stats"),  showWarnings = FALSE, recursive = TRUE)

# quick check
stopifnot(
  file.exists(file.path(DATA_SUPP,"SupplementaryData4.csv")),
  file.exists(file.path(DATA_SUPP,"SupplementaryData2.csv"))
)

# tiny util
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-3) sub("e(-?)(\\d+)", " x 10^\\1\\2", format(p, digits=2, scientific=TRUE))
  else formatC(p, format="f", digits=3)
}
fmt_num <- function(x, d=2) formatC(x, format="f", digits=d)
rng <- function(x, d=2) { x <- x[is.finite(x)]; if(!length(x)) "NA–NA" else paste0(fmt_num(min(x),d),"–",fmt_num(max(x),d)) }

message("setup: OK ✓")