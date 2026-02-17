# Supplementary Data Bundle

Contains the key input tables used to reproduce the main analyses and figures.

**License:** CC BY 4.0 (data tables). Cite the manuscript and this bundle's DOI.

## Contents

- `SupplementaryData1.csv` — rows: 818, columns: 13. Intended use: Figure 1 maps / site locations (meteorological & geographical context)
- `SupplementaryData2.csv` — rows: 118, columns: 13. Intended use: General analysis inputs (see Data Dictionary)
- `SupplementaryData3.csv` — rows: 53, columns: 13. Intended use: General analysis inputs (see Data Dictionary)
- `SupplementaryData4.csv` — rows: 98, columns: 24. Intended use: Figure 2 & sequencing performance (Ct, RPM[log10], N50, coverage); Figure 1 maps / site locations (meteorological & geographical context)

## Mapping to scripts

- `scripts/fig1_meteo_geo.R`: expects site coordinates columns (`lat`/`latitude`, `lon`/`longitude`) and optional site names.
- `scripts/glm_hg_abundance.R`: expects trap counts joined to climate (`date`, `site`, `species`, `count`, and lagged climate).
- `scripts/mir_estimation.R`: expects pooled PCR results (`pool_id`, `pool_size`, `yfv_pcr_positive`).
- `scripts/fig2.R` & `scripts/ct_rpm_simple.R`: expect sequencing stats (`ct`, `rpm` [log10], `n50`, `smart9n_cov10x_pct`, `amplicon_cov_pct`, `host`, `group`, `decomposition_simple`).
- `scripts/phylo/align_iqtree.sh`: uses FASTA alignments; if a metadata CSV exists, include `accession`, `country`, `year`.