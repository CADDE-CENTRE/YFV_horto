getwd()                    # should end with .../_ZENODO
file.exists("data/supplementary/SupplementaryData4.csv")  # TRUE

# 1) scripts/ct_rpm_simple.R (stats only; reads SupplementaryData4.csv) 

# scripts/ct_rpm_simple.R

dat <- readr::read_csv(file.path(DATA_SUPP,"SupplementaryData4.csv"), show_col_types = FALSE)
names(dat) <- tolower(names(dat))

# columns
stopifnot("ct" %in% names(dat))
rpm_col <- intersect(c("rpm","rpm_log10","viral_load_rpmillion"), names(dat))
stopifnot(length(rpm_col) >= 1)
dat <- dat %>% mutate(
  ct  = suppressWarnings(as.numeric(ct)),
  rpm = suppressWarnings(as.numeric(.data[[rpm_col[1]]])),
  n50 = suppressWarnings(as.numeric(dplyr::coalesce(!!!rlang::syms(intersect(
    c("n50","smart9n_n50","smart_9n_n50","smart-9n_n50"), names(dat)
  )))))
)

# groups
host_str    <- tolower(as.character(dplyr::coalesce(dat$host, dat$host)))
species_str <- tolower(as.character(dplyr::coalesce(dat$species)))
is_np <- (!is.na(host_str) & str_detect(host_str, "neotropical\\s*primate")) |
  (!is.na(species_str) & str_detect(species_str, "alouatta|bugio"))
is_hg <- !is.na(species_str) & str_detect(species_str, "leucocel|leucocelaenus|haemagogus\\s*leuc")
dat$group2 <- dplyr::case_when(is_np ~ "NP", is_hg ~ "Hg.leucocelaenus", TRUE ~ NA_character_)

# (i) Spearman Ct vs RPM (log10)
corr_df <- dat %>% filter(is.finite(ct), is.finite(rpm))
ct_rpm_spear <- suppressWarnings(cor.test(corr_df$ct, corr_df$rpm, method="spearman", exact=FALSE))
rho <- unname(ct_rpm_spear$estimate); p_rho <- ct_rpm_spear$p.value; n_corr <- nrow(corr_df)
cat(sprintf("\nCt vs RPM (log10) Spearman (n=%d): rho=%.2f, p=%s\n", n_corr, rho, fmt_p(p_rho)))

# helper
wtest <- function(df, var, digits = 2){
  sub <- df %>%
    filter(group2 %in% c("NP","Hg.leucocelaenus"), is.finite(.data[[var]])) %>%
    select(group2, value = all_of(var))   # <- no .data in tidyselect
  
  v1 <- sub$value[sub$group2 == "NP"]; v2 <- sub$value[sub$group2 == "Hg.leucocelaenus"]
  if (!length(v1) || !length(v2))
    return(list(n1=length(v1), n2=length(v2), med1=NA, med2=NA, r1="NA–NA", r2="NA–NA", p=NA))
  p <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE)$p.value)
  list(n1=length(v1), n2=length(v2),
       med1=median(v1, na.rm=TRUE), med2=median(v2, na.rm=TRUE),
       r1=rng(v1, d=digits), r2=rng(v2, d=digits), p=p)
}

res_ct  <- wtest(dat,"ct", 2)
res_rpm <- wtest(dat,"rpm",2)
res_n50 <- wtest(dat,"n50",0)

cat(sprintf("\nCt NP vs Hg.leucocelaenus: NP median=%.2f (%s, n=%d); HG median=%.2f (%s, n=%d); p=%s\n",
            res_ct$med1, res_ct$r1, res_ct$n1, res_ct$med2, res_ct$r2, res_ct$n2, fmt_p(res_ct$p)))
cat(sprintf("\nRPM(log10) NP vs Hg.leucocelaenus: NP median=%.2f (%s, n=%d); HG median=%.2f (%s, n=%d); p=%s\n",
            res_rpm$med1, res_rpm$r1, res_rpm$n1, res_rpm$med2, res_rpm$r2, res_rpm$n2, fmt_p(res_rpm$p)))
cat(sprintf("\nN50 NP vs Hg.leucocelaenus: NP median=%.0f (%s, n=%d); HG median=%.0f (%s, n=%d); p=%s\n",
            res_n50$med1, res_n50$r1, res_n50$n1, res_n50$med2, res_n50$r2, res_n50$n2, fmt_p(res_n50$p)))

# save a small summary table
summ <- tibble::tibble(
  metric    = c("spearman_ct_vs_rpm","ct_np_vs_hg","rpmlog10_np_vs_hg","n50_np_vs_hg"),
  n_total   = c(n_corr, NA, NA, NA),
  rho       = c(rho, NA, NA, NA),
  p_value   = c(p_rho, res_ct$p, res_rpm$p, res_n50$p),
  n_np      = c(NA, res_ct$n1, res_rpm$n1, res_n50$n1),
  n_hg      = c(NA, res_ct$n2, res_rpm$n2, res_n50$n2),
  median_np = c(NA, res_ct$med1, res_rpm$med1, res_n50$med1),
  median_hg = c(NA, res_ct$med2, res_rpm$med2, res_n50$med2),
  range_np  = c(NA, res_ct$r1,   res_rpm$r1,   res_n50$r1),
  range_hg  = c(NA, res_ct$r2,   res_rpm$r2,   res_n50$r2)
)
readr::write_csv(summ, file.path(OUT,"stats","ct_rpm_stats_summary.csv"))
message("wrote: ", file.path(OUT,"stats","ct_rpm_stats_summary.csv"))
