

################################################################################
#  DE vs ReliefF vs Bottom-up vs Top-down vs ElasticNet vs LASSO
# - Parallel + overall progress bar
################################################################################

library(GEOquery)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(glmnet)
library(pROC)
library(limma)
library(furrr)
library(progressr)
library(R.utils)
library(infotheo)
library(filelock)



options(stringsAsFactors = FALSE)

source("Utilities.R")

# -------------------- Package setup --------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("dplyr", "tibble", "ggplot2", "ggrepel", "glmnet", "pROC",
               "furrr", "progressr", "R.utils", "infotheo", "filelock")
bioc_pkgs <- c("GEOquery", "limma")

for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)




# Tempdir guard (Windows + long runs)
# Some Windows setups (cleanup tools / antivirus / parallel workers) can delete
# the active R session temp directory mid-run. Many base and package utilities
# use tempfile() internally


.ensure_tempdir()

# ------ Parallel + progress ----------------
n_workers <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", unset = "4")))
if (is.na(n_workers) || n_workers < 1) n_workers <- 6
plan(multisession, workers = n_workers)    # set N_WORKERS env var (or edit default "4")
handlers("txtprogressbar")

# Paths and study list

download_dir <- file.path(getwd(), "GEO_Download")
success_file <- file.path(download_dir, "successful_gse_list.txt")

if (!file.exists(success_file)) {
  stop(
    paste0(
      "Cannot find successful_gse_list.txt in: ",
      download_dir,
      "\nRun the data-download/preparation script first."
    )
  )
}

available_gse <- trimws(readLines(success_file, warn = FALSE))
available_gse <- available_gse[nzchar(available_gse)]

# The 11 studies

gse_list <- c(
  "GSE10694",
  "GSE25508",
  "GSE34535",
  "GSE34536",
  "GSE41655",
  "GSE45666",
  "GSE53870",
  "GSE54751",
  "GSE60978",
  "GSE76260",
  "GSE102286"
)

missing_from_cache_index <- setdiff(gse_list, available_gse)
if (length(missing_from_cache_index) > 0L) {
  warning(
    paste0(
      "The following article datasets are not listed in successful_gse_list.txt: ",
      paste(missing_from_cache_index, collapse = ", "),
      ". The script will still try to find their cached CSV files."
    )
  )
}

# Output directory 

results_dir <- file.path(getwd(), "full_safe")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Rolling global progress log 

global_log_file     <- file.path(results_dir, "GLOBAL_PROGRESS.log")
global_log_lockfile <- paste0(global_log_file, ".lock")
completed_file      <- file.path(results_dir, "_completed_gse_ids.txt")
started_file        <- file.path(results_dir, "_started_gse_ids.txt")
GLOBAL_LOG_TAIL_N   <- 60L



# Parameters 
timeout_sec_per_gse <- 800 * 60

k_panel      <- 20
n_splits_auc <- 250L

# Distance-based selection parameters

dist_prefilter_topN <- 250
dist_B      <- 200
dist_sub    <- 0.5
dist_k      <- 20
dist_p      <- 2
dist_pi_thr <- 0.6
dist_eps    <- 0
dist_min_n  <- 3L


# Top-down (overlap) selection parameters 

td_B      <- dist_B
td_sub    <- dist_sub
td_q      <- dist_k      # keep q features (analogous to k for bottom-up)
td_p <- 2
td_pi_thr <- 0.8
td_eps    <- 0
td_min_n  <- dist_min_n


# Initialize global log

tryCatch({ global_log_event("PIPELINE", "START run") }, error = function(e) {})
# -------------------- Run all with progress bar --------------------
with_progress({
  p <- progressor(steps = length(gse_list))
  furrr::future_walk(gse_list, function(gse_id) {
    safe_process_one(gse_id)
    p()
  })
})


.safe_set_single_thread()

