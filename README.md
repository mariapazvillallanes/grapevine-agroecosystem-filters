# Grapevine agroecosystem filters

This repository contains the R scripts and datasets used to analyse culturable
bacterial and fungal communities associated with grapevine agroecosystems under
different management systems, culture media and plant compartments.

## Repository structure

- `data/raw/`  
  Original isolate-level datasets.

- `data/processed/`  
  Processed datasets with unique SampleID used for all analyses.

- `scripts/`  
  R scripts for data processing, statistical analyses and figure generation:
  - `01_add_sampleid.R`
  - `02_bacteria_analysis.R`
  - `03_fungi_analysis.R`
  - `04_figures_main_bacteria.R`
  - `04_figures_main_fungi.R`
  - `05_figures_supplementary_bacteria.R`
  - `06_figures_supplementary_fungi.R`

- `outputs/`  
  Automatically generated figures, tables and R objects (ignored by git).

## Reproducibility

Analyses should be run in the following order:

1. `01_add_sampleid.R`
2. `02_bacteria_analysis.R`
3. `03_fungi_analysis.R`
4. Figure scripts (`04–06`)

All analyses were performed in R using packages from CRAN.

## License

This project is licensed under the MIT License.
