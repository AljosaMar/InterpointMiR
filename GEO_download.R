################################################################################
# GEO_download.R
#
# Download and preprocessing pipeline for the 11 GEO datasets used in:
# "Distance based feature selection for microRNA-cancer classification"
#
# Output expected by data_analysis.R:
#   GEO_Download/<GSE_ID>_miRNA_by_samples.csv
#   GEO_Download/successful_gse_list.txt
#
# Additional bookkeeping:
#   GEO_Download/failed_gse_list.txt
#   GEO_Download/download_report.csv
#   GEO_Download/GEO_download.log
#
# Run this script from the repository root.
################################################################################

options(stringsAsFactors = FALSE)
options(timeout = max(600, getOption("timeout")))

# ==============================================================================
# Package setup
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c(
  "dplyr",
  "tibble"
)

bioc_pkgs <- c(
  "GEOquery",
  "Biobase"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(GEOquery)
library(Biobase)
library(dplyr)
library(tibble)

# ==============================================================================
# Paths and article datasets
# ==============================================================================

download_dir <- file.path(getwd(), "GEO_Download")
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

success_file <- file.path(download_dir, "successful_gse_list.txt")
failed_file <- file.path(download_dir, "failed_gse_list.txt")
report_file <- file.path(download_dir, "download_report.csv")
log_file <- file.path(download_dir, "GEO_download.log")

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

# Sample counts reported in the article.  The two class counts are compared
# without assuming an order here; the phenotype rules below determine which
# class is "normal" and which is "tumor".
expected_class_sizes <- list(
  GSE10694  = c(78L, 78L),
  GSE25508  = c(26L, 34L),
  GSE34535  = c(7L, 7L),
  GSE34536  = c(7L, 7L),
  GSE41655  = c(33L, 15L),
  GSE45666  = c(101L, 15L),
  GSE53870  = c(63L, 9L),
  GSE54751  = c(10L, 10L),
  GSE60978  = c(83L, 6L),
  GSE76260  = c(32L, 32L),
  GSE102286 = c(91L, 88L)
)

# ==============================================================================
# Logging
# ==============================================================================

log_message <- function(gse_id, message) {
  line <- sprintf(
    "%s | %s | %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    gse_id,
    message
  )

  cat(line, "\n")
  write(line, file = log_file, append = TRUE)
  invisible(line)
}

# ==============================================================================
# Sample labels
# ==============================================================================

get_sample_labels <- function(pheno_df, gse_id) {
  if (identical(gse_id, "GSE53870")) {
    if (!("title" %in% colnames(pheno_df))) {
      stop("GSE53870 phenotype data has no 'title' column.")
    }

    titles <- tolower(as.character(pheno_df$title))

    labels <- rep(
      NA_character_,
      length(titles)
    )

    # Tumor samples: e.g. "microRNA profile of case ICC1"
    labels[
      grepl("\\bicc[0-9]+\\b", titles)
    ] <- "tumor"

    # Normal samples: e.g. "microRNA profile of case N1"
    labels[
      grepl("\\bcase n[0-9]+\\b", titles)
    ] <- "normal"

    return(labels)
  }

  tumor_keywords <- paste(
    "\\bcancer\\b",
    "\\btumou?r\\b",
    "\\bcarcinoma\\b",
    "\\bmalignant\\b",
    "\\bneoplasm\\b",
    "\\blesion\\b",
    "\\badenocarcinoma\\b",
    sep = "|"
  )

  normal_keywords <- paste(
    "\\bnormal\\b",
    "\\bcontrol\\b",
    "\\bhealthy\\b",
    "\\bnon-?tumou?r\\b",
    "\\bnon-?neoplastic\\b",
    "\\bbenign\\b",
    "\\badjacent\\b",
    sep = "|"
  )

  cols <- grep(
    "characteristics|source_name|title",
    colnames(pheno_df),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(cols) == 0L) {
    warning("No phenotype text columns matched characteristics/source_name/title.")
    return(rep(NA_character_, nrow(pheno_df)))
  }

  text_blob <- apply(
    pheno_df[, cols, drop = FALSE],
    1,
    function(row) tolower(paste(na.omit(row), collapse = " "))
  )

  labels <- rep(NA_character_, nrow(pheno_df))

  is_tumor <- grepl(
    tumor_keywords,
    text_blob,
    ignore.case = TRUE
  )

  is_normal <- grepl(
    normal_keywords,
    text_blob,
    ignore.case = TRUE
  )

  labels[is_tumor] <- "tumor"
  labels[is_normal] <- "normal"

  if (any(is.na(labels))) {
    warning("Could not assign labels for all samples. Unlabelled samples will be retained with NA labels.")
  }

  labels
}

# ==============================================================================
# Expression preprocessing and probe-to-miRNA mapping
# ==============================================================================

preprocess_expression_data <- function(expr_matrix) {
  expr_vals <- as.vector(expr_matrix)
  expr_vals <- expr_vals[!is.na(expr_vals) & expr_vals > 0]

  if (length(expr_vals) == 0L) {
    return(expr_matrix)
  }

  is_whole_numbers <- all(
    abs(expr_vals - round(expr_vals)) < 1e-10,
    na.rm = TRUE
  )

  has_large_values <- any(
    expr_vals > 50,
    na.rm = TRUE
  )

  if (is_whole_numbers && has_large_values) {
    cat("Detected raw count-like data, applying log2(count + 1)\n")
    expr_matrix <- log2(expr_matrix + 1)
  } else {
    has_negative <- any(
      expr_matrix < 0,
      na.rm = TRUE
    )

    if (
      !has_negative &&
      max(expr_vals, na.rm = TRUE) > 20
    ) {
      cat("Applying log2(x + 1) to expression data\n")
      expr_matrix <- log2(expr_matrix + 1)
    }
  }

  expr_matrix
}

infer_mirna_col <- function(anno_df) {
  prefs <- c(
    "miRNA_ID",
    "miRNA",
    "MIRNA",
    "mirna",
    "miRNA.id",
    "miRNA_ID_REF"
  )

  hit <- prefs[prefs %in% colnames(anno_df)]

  if (length(hit) > 0L) {
    return(hit[1])
  }

  cand <- grep(
    "mirna|microRNA|miR|MIR",
    colnames(anno_df),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(cand) > 0L) {
    return(cand[1])
  }

  NA_character_
}

build_mirna_matrix <- function(eset, gse_id) {
  expr <- Biobase::exprs(eset)
  pheno <- Biobase::pData(eset)

  expr <- preprocess_expression_data(expr)
  pheno$label <- get_sample_labels(pheno, gse_id)

  label_counts <- table(pheno$label)

  if (
    length(label_counts) < 2L ||
    any(label_counts < 3L)
  ) {
    stop(
      "Insufficient labelled samples after phenotype classification; ",
      paste(
        names(label_counts),
        as.integer(label_counts),
        sep = "=",
        collapse = ", "
      )
    )
  }

  gpl_id <- Biobase::annotation(eset)

  if (
    length(gpl_id) != 1L ||
    is.na(gpl_id) ||
    !nzchar(gpl_id)
  ) {
    stop("ExpressionSet has no usable GPL annotation identifier.")
  }

  gpl <- GEOquery::getGEO(gpl_id)
  anno <- GEOquery::Table(gpl)

  if (!("ID" %in% colnames(anno))) {
    stop("GPL annotation has no 'ID' column.")
  }

  mir_col <- infer_mirna_col(anno)

  if (is.na(mir_col)) {
    stop("Could not infer a miRNA column from the GPL annotation table.")
  }

  acc_hits <- grep(
    "MIMAT|Accession",
    colnames(anno),
    value = TRUE
  )

  seq_hits <- grep(
    "Sequence|SEQUENCE",
    colnames(anno),
    value = TRUE
  )

  acc_col <- if (length(acc_hits) > 0L) acc_hits[1] else NA_character_
  seq_col <- if (length(seq_hits) > 0L) seq_hits[1] else NA_character_

  keep_cols <- unique(
    na.omit(
      c(
        "ID",
        mir_col,
        acc_col,
        seq_col
      )
    )
  )

  map <- anno[
    ,
    keep_cols,
    drop = FALSE
  ]

  colnames(map)[colnames(map) == "ID"] <- "probe_id"
  colnames(map)[colnames(map) == mir_col] <- "miRNA"

  map$probe_id <- as.character(map$probe_id)

  expr_df <- tibble::rownames_to_column(
    as.data.frame(expr),
    var = "probe_id"
  )

  expr_df$probe_id <- as.character(expr_df$probe_id)

  expr_annot <- suppressMessages(
    dplyr::left_join(
      map,
      expr_df,
      by = "probe_id"
    )
  )

  expr_by_miRNA <- expr_annot |>
    dplyr::filter(!is.na(miRNA) & miRNA != "") |>
    dplyr::group_by(miRNA) |>
    dplyr::summarise(
      dplyr::across(
        where(is.numeric),
        median,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  if (nrow(expr_by_miRNA) < 1L) {
    stop("No miRNA features remained after GPL annotation and aggregation.")
  }

  mat <- as.data.frame(
    t(
      as.matrix(
        expr_by_miRNA[
          ,
          -1,
          drop = FALSE
        ]
      )
    )
  )

  colnames(mat) <- expr_by_miRNA$miRNA

  mat$label <- pheno$label[
    match(
      rownames(mat),
      rownames(pheno)
    )
  ]

  list(
    mat = mat,
    pheno = pheno,
    gpl_id = gpl_id
  )
}

# ==============================================================================
# Validation
# ==============================================================================

validate_matrix <- function(mat, gse_id) {
  if (!is.data.frame(mat)) {
    stop("Processed object is not a data.frame.")
  }

  if (!("label" %in% colnames(mat))) {
    stop("Processed matrix has no 'label' column.")
  }

  if (ncol(mat) < 2L) {
    stop("Processed matrix contains no miRNA features.")
  }

  labels <- factor(
    mat$label,
    levels = c("normal", "tumor")
  )

  n_normal <- sum(
    labels == "normal",
    na.rm = TRUE
  )

  n_tumor <- sum(
    labels == "tumor",
    na.rm = TRUE
  )

  if (
    n_normal < 3L ||
    n_tumor < 3L
  ) {
    stop(
      sprintf(
        "Insufficient labelled samples: normal=%d, tumor=%d.",
        n_normal,
        n_tumor
      )
    )
  }

  expected <- expected_class_sizes[[gse_id]]

  if (!is.null(expected)) {
    observed <- sort(
      c(
        as.integer(n_normal),
        as.integer(n_tumor)
      )
    )

    if (!identical(observed, sort(as.integer(expected)))) {
      stop(
        sprintf(
          paste0(
            "Labelled class sizes do not match the article dataset: ",
            "observed normal=%d, tumor=%d; expected class sizes %d/%d."
          ),
          n_normal,
          n_tumor,
          expected[1],
          expected[2]
        )
      )
    }
  }

  feature_cols <- setdiff(
    colnames(mat),
    "label"
  )

  numeric_ok <- vapply(
    mat[
      ,
      feature_cols,
      drop = FALSE
    ],
    is.numeric,
    logical(1)
  )

  if (!all(numeric_ok)) {
    stop("At least one miRNA expression column is not numeric.")
  }

  list(
    n_normal = n_normal,
    n_tumor = n_tumor,
    n_features = length(feature_cols),
    n_unlabelled = sum(is.na(labels))
  )
}

read_and_validate_cache <- function(cached_file, gse_id) {
  mat <- read.csv(
    cached_file,
    row.names = 1,
    check.names = FALSE
  )

  validation <- validate_matrix(
    mat,
    gse_id
  )

  list(
    mat = mat,
    validation = validation
  )
}

# ==============================================================================
# Download one GEO study
# ==============================================================================

download_one_gse <- function(gse_id) {
  cached_file <- file.path(
    download_dir,
    paste0(
      gse_id,
      "_miRNA_by_samples.csv"
    )
  )

  if (file.exists(cached_file)) {
    cached <- tryCatch(
      read_and_validate_cache(
        cached_file,
        gse_id
      ),
      error = function(e) {
        log_message(
          gse_id,
          paste0(
            "Existing cache failed validation; rebuilding: ",
            conditionMessage(e)
          )
        )

        NULL
      }
    )

    if (!is.null(cached)) {
      log_message(
        gse_id,
        sprintf(
          paste0(
            "CACHE OK: normal=%d; tumor=%d; ",
            "unlabelled=%d; features=%d"
          ),
          cached$validation$n_normal,
          cached$validation$n_tumor,
          cached$validation$n_unlabelled,
          cached$validation$n_features
        )
      )

      return(
        data.frame(
          GSE = gse_id,
          status = "success_cached",
          ExpressionSet = NA_integer_,
          GPL = NA_character_,
          normal = cached$validation$n_normal,
          tumor = cached$validation$n_tumor,
          unlabelled = cached$validation$n_unlabelled,
          features = cached$validation$n_features,
          file = cached_file,
          error = "",
          stringsAsFactors = FALSE
        )
      )
    }
  }

  log_message(
    gse_id,
    "DOWNLOAD: start"
  )

  gse_obj <- GEOquery::getGEO(
    gse_id,
    GSEMatrix = TRUE
  )

  if (!is.list(gse_obj)) {
    gse_obj <- list(gse_obj)
  }

  if (length(gse_obj) < 1L) {
    stop("GEOquery returned no ExpressionSet objects.")
  }

  candidate_errors <- character(0)

  for (i in seq_along(gse_obj)) {
    log_message(
      gse_id,
      sprintf(
        "Trying ExpressionSet %d/%d",
        i,
        length(gse_obj)
      )
    )

    built <- tryCatch(
      build_mirna_matrix(
        gse_obj[[i]],
        gse_id
      ),
      error = function(e) {
        candidate_errors <<- c(
          candidate_errors,
          paste0(
            "ExpressionSet ",
            i,
            ": ",
            conditionMessage(e)
          )
        )

        NULL
      }
    )

    if (is.null(built)) {
      next
    }

    validation <- tryCatch(
      validate_matrix(
        built$mat,
        gse_id
      ),
      error = function(e) {
        candidate_errors <<- c(
          candidate_errors,
          paste0(
            "ExpressionSet ",
            i,
            ": ",
            conditionMessage(e)
          )
        )

        NULL
      }
    )

    if (is.null(validation)) {
      next
    }

    write.csv(
      built$mat,
      cached_file,
      row.names = TRUE
    )

    log_message(
      gse_id,
      sprintf(
        paste0(
          "DOWNLOAD OK: ExpressionSet=%d; GPL=%s; ",
          "normal=%d; tumor=%d; unlabelled=%d; features=%d"
        ),
        i,
        built$gpl_id,
        validation$n_normal,
        validation$n_tumor,
        validation$n_unlabelled,
        validation$n_features
      )
    )

    return(
      data.frame(
        GSE = gse_id,
        status = "success_downloaded",
        ExpressionSet = i,
        GPL = built$gpl_id,
        normal = validation$n_normal,
        tumor = validation$n_tumor,
        unlabelled = validation$n_unlabelled,
        features = validation$n_features,
        file = cached_file,
        error = "",
        stringsAsFactors = FALSE
      )
    )
  }

  stop(
    paste0(
      "No usable GEO ExpressionSet was found.\n",
      paste(
        candidate_errors,
        collapse = "\n"
      )
    )
  )
}

# ==============================================================================
# Run all 11 article datasets
# ==============================================================================

if (file.exists(log_file)) {
  file.remove(log_file)
}

writeLines(character(0), success_file)
writeLines(character(0), failed_file)

successful_gse <- character(0)
failed_gse <- character(0)
report_rows <- vector(
  "list",
  length(gse_list)
)

for (i in seq_along(gse_list)) {
  gse_id <- gse_list[i]

  result <- tryCatch(
    download_one_gse(gse_id),
    error = function(e) {
      msg <- conditionMessage(e)

      log_message(
        gse_id,
        paste0(
          "FAILED: ",
          msg
        )
      )

      data.frame(
        GSE = gse_id,
        status = "failed",
        ExpressionSet = NA_integer_,
        GPL = NA_character_,
        normal = NA_integer_,
        tumor = NA_integer_,
        unlabelled = NA_integer_,
        features = NA_integer_,
        file = file.path(
          download_dir,
          paste0(
            gse_id,
            "_miRNA_by_samples.csv"
          )
        ),
        error = msg,
        stringsAsFactors = FALSE
      )
    }
  )

  report_rows[[i]] <- result

  if (identical(result$status[1], "failed")) {
    failed_gse <- c(
      failed_gse,
      gse_id
    )
  } else {
    successful_gse <- c(
      successful_gse,
      gse_id
    )
  }

  writeLines(
    successful_gse,
    success_file
  )

  writeLines(
    failed_gse,
    failed_file
  )

  current_report <- do.call(
    rbind,
    report_rows[
      seq_len(i)
    ]
  )

  write.csv(
    current_report,
    report_file,
    row.names = FALSE
  )
}

report <- do.call(
  rbind,
  report_rows
)

rownames(report) <- NULL

write.csv(
  report,
  report_file,
  row.names = FALSE
)

writeLines(
  successful_gse,
  success_file
)

writeLines(
  failed_gse,
  failed_file
)

cat(
  "\nGEO preparation complete.\n",
  "Successful: ",
  length(successful_gse),
  "/",
  length(gse_list),
  "\n",
  "Output directory: ",
  download_dir,
  "\n",
  sep = ""
)

if (length(failed_gse) > 0L) {
  stop(
    paste0(
      "The following GEO datasets failed: ",
      paste(
        failed_gse,
        collapse = ", "
      ),
      ". See ",
      report_file,
      " and ",
      log_file,
      "."
    )
  )
}
