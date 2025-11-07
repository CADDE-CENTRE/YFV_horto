**Code** [![DOI](https://zenodo.org/badge/1002171612.svg)](https://zenodo.org/badge/latestdoi/1002171612) **Extended & Source Data** [![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17427038.svg)](https://doi.org/10.5281/zenodo.17427038)

# Evolution and spillover dynamics of yellow fever at the forest–urban interface

**Code and Data to Reproduce Analyses**

This repository contains the scripts and input tables used to reproduce the main and supplementary analyses for the study entitled *“Evolution and spillover dynamics of yellow fever at the forest–urban interface”* (preprint: https://doi.org/10.21203/rs.3.rs-7048803/v1). The deposit includes curated input data, parameter files, and one-click R scripts that regenerate all figures and statistics referenced in the manuscript and response letter.

---

## Quick start

1. Open the R project: `YFV_horto.Rproj`.
2. In R, run:
   ```r
   source("scripts/00_setup_simple.R")
   ```
3. Then run the scripts listed under **How to reproduce** below.

---

## Repository layout

```
data/
  climate/
    SP_met_and_haemagogus_data_MONTHLY.CSV
    SP_met.RDS
  modelling/
    Febre Amarela Horto.csv
    processed_HortoData.rds
  phylo/
    alignments/           # FASTA alignments used to infer trees
    trees/                # IQ-TREE and BEAST outputs; HAV_tempest.csv
  supplementary/
    SupplementaryData1.csv  # mosquito collections
    SupplementaryData2.csv
    SupplementaryData3.csv  # MIR/MLE meta-analysis
    SupplementaryData4.csv  # sequencing stats; RPM already in log10 units
    DATA_DICTIONARY.md
    SUPPLEMENTARY_DATA_README.md

scripts/
  00_setup_simple.R                # sets paths DATA_* and OUT, loads packages
  ct_rpm_simple.R
  fig_ct_rpm_ridges_simple.R
  diff_vs_ct_simple.R
  np_decomposition_simple.R
  panel_d_np_decomposition_simple.R
  hg_weather_glm_simple.R
  genus_drivers_simple.R
  mir_by_species_simple.R
  supp_fig2_timeline_map_simple.R
  hav_root_to_tip_simple.R
  bioinformatics/                 # one-line command files for consensus, mafft, iqtree
  beast_xmls/                     # BEAST XMLs used in time-scaled analyses
  modelling/                      # contains code required to implement IBM-based modelling of transmission dynamics

outputs/
  fig2/       # Figure 2 panels
  supp_fig/   # Supplementary figures
  stats/      # CSV and HTML stats exports
  modelling/  # IBM modelling & R0 estimation outputs
```

---

## Software

- **R** ≥ 4.2 (tested on macOS 13–14).

**Required packages:** `tidyverse`, `ggplot2`, `ggridges`, `patchwork`, `MASS`, `ggpubr`, `corrplot`, `gt`, `janitor`, `lubridate`, `scatterpie`, `glue`, `scales`.

Install once in R (example):
```r
pkgs <- c(
  "tidyverse","ggplot2","ggridges","patchwork","MASS","ggpubr",
  "corrplot","gt","janitor","lubridate","scatterpie","glue","scales"
)
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

---

## Data notes

- `SupplementaryData4.csv`: column **`rpm`** is already **log10**; do **not** transform again.
- Carcass decomposition: column **`decomposition_simple`** has values `intact_live_2d`, `medium_or_advanced`, or `NA`.
- Phylogenetics: alignments under `data/phylo/alignments`; trees under `data/phylo/trees`; root-to-tip table at `data/phylo/trees/HAV_tempest.csv`.

---

## How to reproduce

> Always start your R session with:
> ```r
> source("scripts/00_setup_simple.R")
> ```

### Figure 2a–b: Ct vs RPM (log10) scatter + ridgelines
- **Script:** `scripts/fig_ct_rpm_ridges_simple.R`  
- **Input:** `supplementary/SupplementaryData4.csv`  
- **Output:** `outputs/fig2/fig_ct_rpm_ridges.pdf`  
- **Console:** prints Spearman ρ and *n*.

### Figure 2d panel: NP decomposition vs genome coverage (ordered bars)
- **Script:** `scripts/panel_d_np_decomposition_simple.R`  
- **Input:** `supplementary/SupplementaryData4.csv`  
- **Output:** `outputs/fig2/panel_d_np_decomposition_coverage.pdf`

### NP decomposition statistics (Ct, RPM log10, N50, coverage)
- **Script:** `scripts/np_decomposition_simple.R`  
- **Input:** `supplementary/SupplementaryData4.csv`  
- **Output:** `outputs/stats/np_decomposition_summary.csv` (medians, ranges, % >70%)

### Between-method coverage difference vs Ct (metagenomics – amplicon)
- **Script:** `scripts/diff_vs_ct_simple.R`  
- **Input:** `supplementary/SupplementaryData4.csv`  
- **Outputs:**  
  - `outputs/supp_fig/diff_vs_ct_global.pdf`  
  - `outputs/stats/diff_vs_ct_stats.csv`

### Ct vs RPM (log10) + group contrasts (NP vs *Hg. leucocelaenus*; N50)
- **Script:** `scripts/ct_rpm_simple.R`  
- **Input:** `supplementary/SupplementaryData4.csv`  
- **Output:** console summaries and CSVs in `outputs/stats/`

### Supplementary Figure 2: timeline of collections by genus + map (genus pies)
- **Script:** `scripts/supp_fig2_timeline_map_simple.R`  
- **Input:** `supplementary/SupplementaryData1.csv`  
- **Output:** `outputs/supp_fig/supp_fig2_timeline_map.pdf`

### Weather GLMs for *Hg. leucocelaenus* (Figure 1c + SI)
- **Script:** `scripts/hg_weather_glm_simple.R`  
- **Input:** `data/climate/SP_met_and_haemagogus_data_MONTHLY.CSV`  
- **Outputs:**  
  - `outputs/fig1c/`  
  - `outputs/supp_fig/variables_time.pdf`  
  - `outputs/supp_fig/correlation_weather.pdf`  
  - Model tables in `outputs/stats/`

### Genus-specific negative-binomial models & weather correlations (SI)
- **Script:** `scripts/genus_drivers_simple.R`  
- **Inputs:** `supplementary/SupplementaryData1.csv`, climate CSV above  
- **Outputs:**  
  - `outputs/supp_fig/SI_Fig_GenusDrivers.pdf`  
  - `outputs/supp_fig/SI_Fig_weather_correlations.pdf`  
  - `outputs/stats/SI_Table_GenusDrivers.csv` and `.html`

### HAV root-to-tip regression (TempEst export)
- **Script:** `scripts/hav_root_to_tip_simple.R`  
- **Input:** `data/phylo/trees/HAV_tempest.csv`  
- **Outputs:**  
  - `outputs/supp_fig/hav_tempest_root_to_tip.pdf`  
  - `outputs/stats/hav_tempest_root_to_tip_lm.csv`

---

## Phylogenetics

- **Alignments (FASTA):** `data/phylo/alignments/`  
- **Trees (IQ-TREE `.tre`, BEAST `.MCC.tre`):** `data/phylo/trees/`  
- **BEAST XMLs:** `scripts/beast_xmls/`  
- **Command files (MAFFT, IQ-TREE, consensus):** `scripts/bioinformatics/`  
- Phylogenetic figures are regenerated from the included trees.

---

## Modelling & R0 estimation

- **Data:** `data/modelling/processed_HortoData.rds`
- **Model:** `functions/IBM_model.R`
- **Script:** `scripts/modelling/3_horto_YFV_R0_estimation/2_outbreakInference_parameterScan_yesImportations_yesStartDate.R`

## Licence

- **Code:** MIT License (© Authors)  
- **Data (tables, trees, alignments):** CC BY 4.0

---

## How to cite

Please cite the preprint (https://doi.org/10.21203/rs.3.rs-7048803/v1) or the manuscript (when available), and this dataset/software deposit (Zenodo DOI when minted).

**Example:** Author list. Title. *Journal* (Year). Code and data: Zenodo DOI: **[TBD]**.

---

## Contact

**Corresponding author:** Nuno R. Faria — <nfaria@ic.ac.uk>

For issues or requests, please use the contact above.
