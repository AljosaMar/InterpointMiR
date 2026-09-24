# InterpointMiR

R code accompanying the article **“Distance based feature selection for microRNA-cancer classification.”**

The repository contains the real-data analysis and simulation code for the distance-distribution feature-selection procedures studied in the article.

## Repository contents

- `GEO_download.R` — downloads and prepares the 11 GEO datasets used in the article.
- `Utilities.R` — core feature-selection, evaluation, plotting, and grouped train/test split utilities.
- `data_analysis.R` — real-data comparison of DE, ReliefF, mRMR, LASSO, Elastic Net, Bottom-up, and Top-down procedures.
- `Simulations.R` — simulation experiments reported in the article and Supplement.

## GEO datasets

The real-data analysis uses the following GEO studies:

- GSE10694
- GSE25508
- GSE34535
- GSE34536
- GSE41655
- GSE45666
- GSE53870
- GSE54751
- GSE60978
- GSE76260
- GSE102286

`GEO_download.R` downloads the GEO Series Matrix data through `GEOquery`, maps platform probes to miRNA identifiers, collapses multiple probes for the same miRNA by their median expression, assigns the binary labels `normal` and `tumor`, and writes the processed matrices expected by `data_analysis.R`.

The processed files are stored as

```text
GEO_Download/<GSE_ID>_miRNA_by_samples.csv
```

with samples in rows, miRNAs in columns, and a final `label` column.

The script also creates

```text
GEO_Download/successful_gse_list.txt
GEO_Download/failed_gse_list.txt
GEO_Download/download_report.csv
GEO_Download/GEO_download.log
```

`successful_gse_list.txt` contains one successfully prepared GEO accession per line and is read directly by `data_analysis.R`.

## Running the analysis

Run the scripts from the repository root.

First prepare the GEO data:

```r
source("GEO_download.R")
```

Then run the real-data analysis:

```r
source("data_analysis.R")
```

The real-data results are written to the output directory specified in `data_analysis.R`.

Run the simulations separately:

```r
source("Simulations.R")
```

The simulation script writes its tables, raw simulation results, and figures to the output directories defined inside `Simulations.R`.

## Real-data analysis settings

For reproduction of the article, the real-data analysis uses:

- 250 repeated stratified 70/30 train/test splits;
- patient-level grouping for paired or repeated specimens;
- feature selection and prefiltering performed on the training data only;
- a top-250 limma differential-expression candidate prefilter for the distance-based procedures;
- 200 stability-selection repetitions;
- subsampling fraction 0.5;
- Bottom-up stability threshold 0.5;
- Top-down stability threshold 0.8;
- maximum Bottom-up panel size 20;
- minimum Top-down panel size 20;
- \(l_2\)-distance, implemented through `p = 2`;
- unpenalized logistic regression for the final AUC comparison, with standardization fitted on the training set only.


## Parallel computation

The real-data analysis uses parallel workers through the `future`/`furrr` framework. The number of workers can be controlled with the `N_WORKERS` environment variable. The simulations also use parallel computation and can be computationally intensive.

## R packages

The scripts use packages from CRAN and Bioconductor, including:

- GEOquery
- Biobase
- limma
- dplyr
- tibble
- ggplot2
- ggrepel
- glmnet
- pROC
- furrr
- future
- progressr
- R.utils
- infotheo
- filelock
- foreach
- doParallel
- doRNG
- patchwork

The scripts install missing packages where appropriate.
