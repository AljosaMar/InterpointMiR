

#
# Simulations.R
# 
# Same data generation for the first two supplement tables:
# n = 50 per group, d = 100, relevant features = 86:100.


source("Utilities.R")
library(ggplot2)
library(patchwork)

set.seed(123)

n_mc <- 100L
n_per_group <- 50L
d <- 100L
delta_grid <- c(0, 0.5, 1, 3, 5)
relevant_features <- 86:100

methods <- c(
  "Uncorrected",
  "B-H",
  "B-Y",
  "Bonferroni",
  "Bottom-up (P_A)",
  "Top-down (O_A)",
  "Bottom-up stability",
  "Top-down stability"
)

output_dir <- file.path(getwd(), "simulation_outputs_first_example")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# Location experiment

n_cores <- 4L
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(20260326)

clusterEvalQ(cl, {
  source("Utilities.R")
  NULL
})

clusterExport(cl, varlist = c("methods", "relevant_features", "n_per_group", "d"),
              envir = environment())

acc_mean <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                   dimnames = list(methods, as.character(delta_grid)))
acc_sd   <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                   dimnames = list(methods, as.character(delta_grid)))
acc_se   <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                   dimnames = list(methods, as.character(delta_grid)))

cmp_mean <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                   dimnames = list(methods, as.character(delta_grid)))
cmp_sd   <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                   dimnames = list(methods, as.character(delta_grid)))
cmp_se   <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                   dimnames = list(methods, as.character(delta_grid)))

exact_prob <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                     dimnames = list(methods, as.character(delta_grid)))
size_mean  <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                     dimnames = list(methods, as.character(delta_grid)))
size_se    <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                     dimnames = list(methods, as.character(delta_grid)))
jacc_mean  <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                     dimnames = list(methods, as.character(delta_grid)))
jacc_se    <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                     dimnames = list(methods, as.character(delta_grid)))
fp_prob    <- matrix(NA_real_, nrow = length(methods), ncol = length(delta_grid),
                     dimnames = list(methods, as.character(delta_grid)))

for (g in seq_along(delta_grid)) {
  delta_mu <- delta_grid[g]
  
  cat(sprintf("[%s] Starting Delta = %s with %d parallel workers\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              delta_mu, n_cores))
  flush.console()
  
  mc_out <- foreach(
    b = seq_len(n_mc),
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = "stats"
  ) %dorng% {
    
    X <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1), nrow = n_per_group, ncol = d)
    Y <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1), nrow = n_per_group, ncol = d)
    Y[, relevant_features] <- Y[, relevant_features, drop = FALSE] + delta_mu
    
    colnames(X) <- paste0("V", seq_len(d))
    colnames(Y) <- paste0("V", seq_len(d))
    
    pvals <- rep(NA_real_, d)
    for (j in seq_len(d)) {
      pvals[j] <- t.test(X[, j], Y[, j], var.equal = TRUE, alternative = "two.sided")$p.value
    }
    
    bu_idx <- bottom_up_pa(X, Y, k = d, p = 1, eps = 0)$idx
    td_idx <- top_down_oa(X, Y, q = 1, p = 1, eps = 0)$idx
    
    bu_stab_idx <- stability_bottom_up(
      X, Y,
      B = 200,
      sub = 0.7,
      pi_thr = 0.5,
      seed = 1000L + 100L * g + b,
      progress_every = 0L
    )$selected_idx
    
    td_stab_idx <- stability_top_down(
      X, Y,
      B = 200,
      sub = 0.7,
      pi_thr = 0.5,
      seed = 2000L + 100L * g + b,
      progress_every = 0L
    )$selected_idx
    
    selected_list <- list(
      "Uncorrected"         = which(pvals < 0.05),
      "B-H"                 = which(p.adjust(pvals, method = "BH") < 0.05),
      "B-Y"                 = which(p.adjust(pvals, method = "BY") < 0.05),
      "Bonferroni"          = which(p.adjust(pvals, method = "bonferroni") < 0.05),
      "Bottom-up (P_A)"     = bu_idx,
      "Top-down (O_A)"      = td_idx,
      "Bottom-up stability" = bu_stab_idx,
      "Top-down stability"  = td_stab_idx
    )
    
    if (delta_mu == 0) {
      true_set <- integer(0)
    } else {
      true_set <- relevant_features
    }
    
    out_list <- vector("list", length(methods))
    
    for (i in seq_along(methods)) {
      m <- methods[i]
      A <- sort(unique(as.integer(selected_list[[m]])))
      tp <- length(intersect(A, true_set))
      fp <- length(setdiff(A, true_set))
      union_size <- length(union(A, true_set))
      
      acc_val <- if (length(A) == 0L) 1 else tp / length(A)
      cmp_val <- if (delta_mu == 0) fp else tp / length(true_set)
      exact_val <- as.integer(setequal(A, true_set))
      size_val <- length(A)
      jacc_val <- if (union_size == 0L) 1 else tp / union_size
      fp_val <- as.integer(fp > 0L)
      
      out_list[[i]] <- data.frame(
        Replicate = b,
        Delta = delta_mu,
        Method = m,
        Accuracy = acc_val,
        Completeness = cmp_val,
        ExactSupportRecovery = exact_val,
        SelectedPanelSize = size_val,
        JaccardSimilarity = jacc_val,
        PFalsePositive = fp_val,
        stringsAsFactors = FALSE
      )
    }
    
    do.call(rbind, out_list)
  }
  
  for (m in methods) {
    idx <- mc_out$Method == m
    
    acc_mean[m, g] <- mean(mc_out$Accuracy[idx])
    acc_sd[m, g]   <- sd(mc_out$Accuracy[idx])
    acc_se[m, g]   <- acc_sd[m, g] / sqrt(n_mc)
    
    cmp_mean[m, g] <- mean(mc_out$Completeness[idx])
    cmp_sd[m, g]   <- sd(mc_out$Completeness[idx])
    cmp_se[m, g]   <- cmp_sd[m, g] / sqrt(n_mc)
    
    exact_prob[m, g] <- mean(mc_out$ExactSupportRecovery[idx])
    size_mean[m, g]  <- mean(mc_out$SelectedPanelSize[idx])
    size_se[m, g]    <- sd(mc_out$SelectedPanelSize[idx]) / sqrt(n_mc)
    jacc_mean[m, g]  <- mean(mc_out$JaccardSimilarity[idx])
    jacc_se[m, g]    <- sd(mc_out$JaccardSimilarity[idx]) / sqrt(n_mc)
    fp_prob[m, g]    <- mean(mc_out$PFalsePositive[idx])
  }
  
  cat(sprintf("[%s] Finished Delta = %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              delta_mu))
  flush.console()
}

stopCluster(cl)
registerDoSEQ()

result_rows <- list()
row_counter <- 1L

for (m in methods) {
  row_accuracy <- c(Method = m, Metric = "Accuracy")
  row_completeness <- c(Method = m, Metric = "Completeness / #FP at null")
  row_exact <- c(Method = m, Metric = "Exact support recovery")
  row_size <- c(Method = m, Metric = "Selected panel size")
  row_jacc <- c(Method = m, Metric = "Jaccard similarity")
  row_fp <- c(Method = m, Metric = "P(false positive)")
  
  for (g in seq_along(delta_grid)) {
    delta_name <- paste0("Delta_", delta_grid[g])
    
    row_accuracy[delta_name] <- sprintf("%.3f ± %.3f", acc_mean[m, g], acc_se[m, g])
    row_completeness[delta_name] <- sprintf("%.3f ± %.3f", cmp_mean[m, g], cmp_se[m, g])
    row_exact[delta_name] <- sprintf("%.3f ± %.3f", exact_prob[m, g], sqrt(exact_prob[m, g] * (1 - exact_prob[m, g]) / (n_mc - 1)))
    row_size[delta_name]  <- sprintf("%.3f ± %.3f", size_mean[m, g], size_se[m, g])
    row_jacc[delta_name]  <- sprintf("%.3f ± %.3f", jacc_mean[m, g], jacc_se[m, g])
    row_fp[delta_name]    <- sprintf("%.3f ± %.3f", fp_prob[m, g], sqrt(fp_prob[m, g] * (1 - fp_prob[m, g]) / (n_mc - 1)))
  }
  
  result_rows[[row_counter]] <- as.data.frame(as.list(row_accuracy), stringsAsFactors = FALSE)
  row_counter <- row_counter + 1L
  result_rows[[row_counter]] <- as.data.frame(as.list(row_completeness), stringsAsFactors = FALSE)
  row_counter <- row_counter + 1L
  result_rows[[row_counter]] <- as.data.frame(as.list(row_exact), stringsAsFactors = FALSE)
  row_counter <- row_counter + 1L
  result_rows[[row_counter]] <- as.data.frame(as.list(row_size), stringsAsFactors = FALSE)
  row_counter <- row_counter + 1L
  result_rows[[row_counter]] <- as.data.frame(as.list(row_jacc), stringsAsFactors = FALSE)
  row_counter <- row_counter + 1L
  result_rows[[row_counter]] <- as.data.frame(as.list(row_fp), stringsAsFactors = FALSE)
  row_counter <- row_counter + 1L
}

results_table <- do.call(rbind, result_rows)
row.names(results_table) <- NULL

write.csv(results_table,
          file = file.path(output_dir, "simulation_table_first_example.csv"),
          row.names = FALSE,
          quote = TRUE)

write.table(results_table,
            file = file.path(output_dir, "simulation_table_first_example.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

capture.output(print(results_table, row.names = FALSE, right = FALSE),
               file = file.path(output_dir, "simulation_table_first_example.txt"))

plot_rows <- list()
plot_counter <- 1L

for (m in methods) {
  for (g in seq_along(delta_grid)) {
    if (delta_grid[g] > 0) {
      plot_rows[[plot_counter]] <- data.frame(
        Method = m,
        Delta = delta_grid[g],
        Metric = "Accuracy",
        Mean = acc_mean[m, g],
        SE = acc_se[m, g],
        stringsAsFactors = FALSE
      )
      plot_counter <- plot_counter + 1L
      
      plot_rows[[plot_counter]] <- data.frame(
        Method = m,
        Delta = delta_grid[g],
        Metric = "Completeness",
        Mean = cmp_mean[m, g],
        SE = cmp_se[m, g],
        stringsAsFactors = FALSE
      )
      plot_counter <- plot_counter + 1L
      
      plot_rows[[plot_counter]] <- data.frame(
        Method = m,
        Delta = delta_grid[g],
        Metric = "Selected panel size",
        Mean = size_mean[m, g],
        SE = size_se[m, g],
        stringsAsFactors = FALSE
      )
      plot_counter <- plot_counter + 1L
      
      plot_rows[[plot_counter]] <- data.frame(
        Method = m,
        Delta = delta_grid[g],
        Metric = "Exact support recovery",
        Mean = exact_prob[m, g],
        SE = sqrt(pmax(exact_prob[m, g] * (1 - exact_prob[m, g]) / n_mc, 0)),
        stringsAsFactors = FALSE
      )
      plot_counter <- plot_counter + 1L
      
      plot_rows[[plot_counter]] <- data.frame(
        Method = m,
        Delta = delta_grid[g],
        Metric = "Jaccard similarity",
        Mean = jacc_mean[m, g],
        SE = jacc_se[m, g],
        stringsAsFactors = FALSE
      )
      plot_counter <- plot_counter + 1L
      
      plot_rows[[plot_counter]] <- data.frame(
        Method = m,
        Delta = delta_grid[g],
        Metric = "P(false positive)",
        Mean = fp_prob[m, g],
        SE = sqrt(pmax(fp_prob[m, g] * (1 - fp_prob[m, g]) / n_mc, 0)),
        stringsAsFactors = FALSE
      )
      plot_counter <- plot_counter + 1L
    }
  }
}

plot_df <- do.call(rbind, plot_rows)

write.csv(plot_df,
          file = file.path(output_dir, "simulation_plot_data_first_example.csv"),
          row.names = FALSE)

plot_df$Method <- factor(plot_df$Method, levels = methods)

method_colors <- c(
  "Uncorrected"         = "black",
  "B-H"                 = "blue",
  "B-Y"                 = "purple",
  "Bonferroni"          = "orange",
  "Bottom-up (P_A)"     = "red",
  "Top-down (O_A)"      = "green3",
  "Bottom-up stability" = "red",
  "Top-down stability"  = "green3"
)

method_linetypes <- c(
  "Uncorrected"         = "dashed",
  "B-H"                 = "dotted",
  "B-Y"                 = "dotdash",
  "Bonferroni"          = "longdash",
  "Bottom-up (P_A)"     = "dashed",
  "Top-down (O_A)"      = "dashed",
  "Bottom-up stability" = "solid",
  "Top-down stability"  = "solid"
)

method_shapes <- c(
  "Uncorrected"         = 16,
  "B-H"                 = 17,
  "B-Y"                 = 15,
  "Bonferroni"          = 18,
  "Bottom-up (P_A)"     = 19,
  "Top-down (O_A)"      = 8,
  "Bottom-up stability" = 1,
  "Top-down stability"  = 2
)

base_theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

p1 <- ggplot(
  plot_df[plot_df$Metric == "Accuracy", ],
  aes(x = Delta, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.08, linewidth = 0.5) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(Delta * mu), y = "Accuracy", title = "1. Accuracy") +
  base_theme

p2 <- ggplot(
  plot_df[plot_df$Metric == "Completeness", ],
  aes(x = Delta, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.08, linewidth = 0.5) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(Delta * mu), y = "Completeness", title = "2. Completeness") +
  base_theme

p3 <- ggplot(
  plot_df[plot_df$Metric == "Selected panel size", ],
  aes(x = Delta, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.08, linewidth = 0.5) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(Delta * mu), y = "Selected panel size", title = "3. Selected panel size") +
  base_theme

p4 <- ggplot(
  plot_df[plot_df$Metric == "Exact support recovery", ],
  aes(x = Delta, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.08, linewidth = 0.5) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(Delta * mu), y = "Exact support recovery", title = "4. Exact support recovery") +
  base_theme

p5 <- ggplot(
  plot_df[plot_df$Metric == "Jaccard similarity", ],
  aes(x = Delta, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.08, linewidth = 0.5) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(Delta * mu), y = "Jaccard similarity", title = "5. Jaccard similarity") +
  base_theme

p6 <- ggplot(
  plot_df[plot_df$Metric == "P(false positive)", ],
  aes(x = Delta, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.08, linewidth = 0.5) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(Delta * mu), y = "P(false positive)", title = "6. P(false positive)") +
  base_theme

combined_plot <- ((p1 | p2) / (p3 | p4) / (p5 | p6)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, "simulation_metrics_first_example_6panel.png"),
       plot = combined_plot, width = 12, height = 13.5, dpi = 300)

ggsave(filename = file.path(output_dir, "simulation_metrics_first_example_6panel.pdf"),
       plot = combined_plot, width = 12, height = 13.5)

print(results_table, row.names = FALSE, right = FALSE)
print(combined_plot)
cat("\nSaved outputs in:\n", output_dir, "\n", sep = "")



# Scale experiment

n_cores_scale <- 8L
cl <- makeCluster(n_cores_scale)
registerDoParallel(cl)
registerDoRNG(20260326)

scale_var_grid <- c(0.9, 0.8, 0.7, 0.6)

scale_acc_mean <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                         dimnames = list(methods, as.character(scale_var_grid)))
scale_acc_sd   <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                         dimnames = list(methods, as.character(scale_var_grid)))
scale_acc_se   <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                         dimnames = list(methods, as.character(scale_var_grid)))

scale_cmp_mean <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                         dimnames = list(methods, as.character(scale_var_grid)))
scale_cmp_sd   <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                         dimnames = list(methods, as.character(scale_var_grid)))
scale_cmp_se   <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                         dimnames = list(methods, as.character(scale_var_grid)))

scale_exact_prob <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                           dimnames = list(methods, as.character(scale_var_grid)))
scale_size_mean  <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                           dimnames = list(methods, as.character(scale_var_grid)))
scale_size_se    <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                           dimnames = list(methods, as.character(scale_var_grid)))
scale_jacc_mean  <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                           dimnames = list(methods, as.character(scale_var_grid)))
scale_jacc_se    <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                           dimnames = list(methods, as.character(scale_var_grid)))
scale_fp_prob    <- matrix(NA_real_, nrow = length(methods), ncol = length(scale_var_grid),
                           dimnames = list(methods, as.character(scale_var_grid)))

scale_mc_results_list <- vector("list", length(scale_var_grid))

for (g in seq_along(scale_var_grid)) {
  scale_var <- scale_var_grid[g]
  
  cat(sprintf("[%s] Starting scale variance = %.1f with %d parallel workers\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              scale_var, n_cores_scale))
  flush.console()
  
  mc_out <- foreach(
    b = seq_len(n_mc),
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = "stats",
    .export = c("methods", "relevant_features", "n_per_group", "d",
                "bottom_up_pa", "top_down_oa",
                "stability_bottom_up", "stability_top_down")
  ) %dorng% {
    
    X <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1), nrow = n_per_group, ncol = d)
    Y <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1), nrow = n_per_group, ncol = d)
    
    # Relevant features differ only in variance
    Y[, relevant_features] <- sqrt(scale_var) * Y[, relevant_features, drop = FALSE]
    
    colnames(X) <- paste0("V", seq_len(d))
    colnames(Y) <- paste0("V", seq_len(d))
    
    pvals <- rep(NA_real_, d)
    for (j in seq_len(d)) {
      pvals[j] <- var.test(X[, j], Y[, j], alternative = "two.sided")$p.value
    }
    
    bu_idx <- bottom_up_pa(X, Y, k = d, p = 1, eps = 0)$idx
    td_idx <- top_down_oa(X, Y, q = 1, p = 1, eps = 0)$idx
    
    bu_stab_idx <- stability_bottom_up(
      X, Y,
      B = 200,
      sub = 0.7,
      pi_thr = 0.5,
      seed = 1000L + 100L * g + b,
      progress_every = 0L
    )$selected_idx
    
    td_stab_idx <- stability_top_down(
      X, Y,
      B = 200,
      sub = 0.7,
      pi_thr = 0.5,
      seed = 2000L + 100L * g + b,
      progress_every = 0L
    )$selected_idx
    
    selected_list <- list(
      "Uncorrected"         = which(pvals < 0.05),
      "B-H"                 = which(p.adjust(pvals, method = "BH") < 0.05),
      "B-Y"                 = which(p.adjust(pvals, method = "BY") < 0.05),
      "Bonferroni"          = which(p.adjust(pvals, method = "bonferroni") < 0.05),
      "Bottom-up (P_A)"     = bu_idx,
      "Top-down (O_A)"      = td_idx,
      "Bottom-up stability" = bu_stab_idx,
      "Top-down stability"  = td_stab_idx
    )
    
    true_set <- relevant_features
    
    out_list <- vector("list", length(methods))
    
    for (i in seq_along(methods)) {
      m <- methods[i]
      A <- sort(unique(as.integer(selected_list[[m]])))
      tp <- length(intersect(A, true_set))
      fp <- length(setdiff(A, true_set))
      union_size <- length(union(A, true_set))
      
      acc_val <- if (length(A) == 0L) 1 else tp / length(A)
      cmp_val <- tp / length(true_set)
      exact_val <- as.integer(setequal(A, true_set))
      size_val <- length(A)
      jacc_val <- if (union_size == 0L) 1 else tp / union_size
      fp_val <- as.integer(fp > 0L)
      
      out_list[[i]] <- data.frame(
        Replicate = b,
        ScaleVar = scale_var,
        Method = m,
        Accuracy = acc_val,
        Completeness = cmp_val,
        ExactSupportRecovery = exact_val,
        SelectedPanelSize = size_val,
        JaccardSimilarity = jacc_val,
        PFalsePositive = fp_val,
        stringsAsFactors = FALSE
      )
    }
    
    do.call(rbind, out_list)
  }
  
  scale_mc_results_list[[g]] <- mc_out
  
  for (m in methods) {
    idx <- mc_out$Method == m
    
    scale_acc_mean[m, g] <- mean(mc_out$Accuracy[idx])
    scale_acc_sd[m, g]   <- sd(mc_out$Accuracy[idx])
    scale_acc_se[m, g]   <- scale_acc_sd[m, g] / sqrt(n_mc)
    
    scale_cmp_mean[m, g] <- mean(mc_out$Completeness[idx])
    scale_cmp_sd[m, g]   <- sd(mc_out$Completeness[idx])
    scale_cmp_se[m, g]   <- scale_cmp_sd[m, g] / sqrt(n_mc)
    
    scale_exact_prob[m, g] <- mean(mc_out$ExactSupportRecovery[idx])
    scale_size_mean[m, g]  <- mean(mc_out$SelectedPanelSize[idx])
    scale_size_se[m, g]    <- sd(mc_out$SelectedPanelSize[idx]) / sqrt(n_mc)
    scale_jacc_mean[m, g]  <- mean(mc_out$JaccardSimilarity[idx])
    scale_jacc_se[m, g]    <- sd(mc_out$JaccardSimilarity[idx]) / sqrt(n_mc)
    scale_fp_prob[m, g]    <- mean(mc_out$PFalsePositive[idx])
  }
  
  cat(sprintf("[%s] Finished scale variance = %.1f\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              scale_var))
  flush.console()
}

scale_mc_results <- do.call(rbind, scale_mc_results_list)
row.names(scale_mc_results) <- NULL

scale_result_rows <- list()
scale_row_counter <- 1L

for (m in methods) {
  row_accuracy <- c(Method = m, Metric = "Accuracy")
  row_completeness <- c(Method = m, Metric = "Completeness")
  row_exact <- c(Method = m, Metric = "Exact support recovery")
  row_size <- c(Method = m, Metric = "Selected panel size")
  row_jacc <- c(Method = m, Metric = "Jaccard similarity")
  row_fp <- c(Method = m, Metric = "P(false positive)")
  
  for (g in seq_along(scale_var_grid)) {
    scale_name <- paste0("Var_", format(scale_var_grid[g], nsmall = 1, trim = TRUE))
    
    row_accuracy[scale_name]     <- sprintf("%.3f ± %.3f", scale_acc_mean[m, g], scale_acc_se[m, g])
    row_completeness[scale_name] <- sprintf("%.3f ± %.3f", scale_cmp_mean[m, g], scale_cmp_se[m, g])
    row_exact[scale_name]        <- sprintf("%.3f ± %.3f", scale_exact_prob[m, g], sqrt(scale_exact_prob[m, g] * (1 - scale_exact_prob[m, g]) / (n_mc - 1)))
    row_size[scale_name]         <- sprintf("%.3f ± %.3f", scale_size_mean[m, g], scale_size_se[m, g])
    row_jacc[scale_name]         <- sprintf("%.3f ± %.3f", scale_jacc_mean[m, g], scale_jacc_se[m, g])
    row_fp[scale_name]           <- sprintf("%.3f ± %.3f", scale_fp_prob[m, g], sqrt(scale_fp_prob[m, g] * (1 - scale_fp_prob[m, g]) / (n_mc - 1)))
  }
  
  scale_result_rows[[scale_row_counter]] <- as.data.frame(as.list(row_accuracy), stringsAsFactors = FALSE)
  scale_row_counter <- scale_row_counter + 1L
  scale_result_rows[[scale_row_counter]] <- as.data.frame(as.list(row_completeness), stringsAsFactors = FALSE)
  scale_row_counter <- scale_row_counter + 1L
  scale_result_rows[[scale_row_counter]] <- as.data.frame(as.list(row_exact), stringsAsFactors = FALSE)
  scale_row_counter <- scale_row_counter + 1L
  scale_result_rows[[scale_row_counter]] <- as.data.frame(as.list(row_size), stringsAsFactors = FALSE)
  scale_row_counter <- scale_row_counter + 1L
  scale_result_rows[[scale_row_counter]] <- as.data.frame(as.list(row_jacc), stringsAsFactors = FALSE)
  scale_row_counter <- scale_row_counter + 1L
  scale_result_rows[[scale_row_counter]] <- as.data.frame(as.list(row_fp), stringsAsFactors = FALSE)
  scale_row_counter <- scale_row_counter + 1L
}

scale_results_table <- do.call(rbind, scale_result_rows)
row.names(scale_results_table) <- NULL

write.csv(scale_results_table,
          file = file.path(output_dir, "simulation_table_scale_example_parallel.csv"),
          row.names = FALSE,
          quote = TRUE)

write.table(scale_results_table,
            file = file.path(output_dir, "simulation_table_scale_example_parallel.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

capture.output(print(scale_results_table, row.names = FALSE, right = FALSE),
               file = file.path(output_dir, "simulation_table_scale_example_parallel.txt"))

write.csv(scale_mc_results,
          file = file.path(output_dir, "simulation_raw_scale_example_parallel.csv"),
          row.names = FALSE)

saveRDS(
  list(
    scale_var_grid = scale_var_grid,
    methods = methods,
    mc_results = scale_mc_results,
    results_table = scale_results_table,
    scale_acc_mean = scale_acc_mean,
    scale_acc_sd = scale_acc_sd,
    scale_acc_se = scale_acc_se,
    scale_cmp_mean = scale_cmp_mean,
    scale_cmp_sd = scale_cmp_sd,
    scale_cmp_se = scale_cmp_se,
    scale_exact_prob = scale_exact_prob,
    scale_size_mean = scale_size_mean,
    scale_size_se = scale_size_se,
    scale_jacc_mean = scale_jacc_mean,
    scale_jacc_se = scale_jacc_se,
    scale_fp_prob = scale_fp_prob
  ),
  file = file.path(output_dir, "simulation_scale_example_parallel.rds")
)

stopCluster(cl)
registerDoSEQ()

print(scale_results_table, row.names = FALSE, right = FALSE)
cat("\nSaved parallel scale-experiment outputs in:\n", output_dir, "\n", sep = "")



# Plots for the scale experiment


scale_plot_rows <- list()
scale_plot_counter <- 1L

for (m in methods) {
  for (g in seq_along(scale_var_grid)) {
    scale_plot_rows[[scale_plot_counter]] <- data.frame(
      Method = m,
      ScaleVar = scale_var_grid[g],
      Metric = "Accuracy",
      Mean = scale_acc_mean[m, g],
      SE = scale_acc_se[m, g],
      stringsAsFactors = FALSE
    )
    scale_plot_counter <- scale_plot_counter + 1L
    
    scale_plot_rows[[scale_plot_counter]] <- data.frame(
      Method = m,
      ScaleVar = scale_var_grid[g],
      Metric = "Completeness",
      Mean = scale_cmp_mean[m, g],
      SE = scale_cmp_se[m, g],
      stringsAsFactors = FALSE
    )
    scale_plot_counter <- scale_plot_counter + 1L
    
    scale_plot_rows[[scale_plot_counter]] <- data.frame(
      Method = m,
      ScaleVar = scale_var_grid[g],
      Metric = "Selected panel size",
      Mean = scale_size_mean[m, g],
      SE = scale_size_se[m, g],
      stringsAsFactors = FALSE
    )
    scale_plot_counter <- scale_plot_counter + 1L
    
    scale_plot_rows[[scale_plot_counter]] <- data.frame(
      Method = m,
      ScaleVar = scale_var_grid[g],
      Metric = "Exact support recovery",
      Mean = scale_exact_prob[m, g],
      SE = sqrt(pmax(scale_exact_prob[m, g] * (1 - scale_exact_prob[m, g]) / n_mc, 0)),
      stringsAsFactors = FALSE
    )
    scale_plot_counter <- scale_plot_counter + 1L
    
    scale_plot_rows[[scale_plot_counter]] <- data.frame(
      Method = m,
      ScaleVar = scale_var_grid[g],
      Metric = "Jaccard similarity",
      Mean = scale_jacc_mean[m, g],
      SE = scale_jacc_se[m, g],
      stringsAsFactors = FALSE
    )
    scale_plot_counter <- scale_plot_counter + 1L
    
    scale_plot_rows[[scale_plot_counter]] <- data.frame(
      Method = m,
      ScaleVar = scale_var_grid[g],
      Metric = "P(false positive)",
      Mean = scale_fp_prob[m, g],
      SE = sqrt(pmax(scale_fp_prob[m, g] * (1 - scale_fp_prob[m, g]) / n_mc, 0)),
      stringsAsFactors = FALSE
    )
    scale_plot_counter <- scale_plot_counter + 1L
  }
}

scale_plot_df <- do.call(rbind, scale_plot_rows)

write.csv(scale_plot_df,
          file = file.path(output_dir, "simulation_plot_data_scale_example.csv"),
          row.names = FALSE)

scale_plot_df$Method <- factor(scale_plot_df$Method, levels = methods)

method_colors <- c(
  "Uncorrected"         = "black",
  "B-H"                 = "blue",
  "B-Y"                 = "purple",
  "Bonferroni"          = "orange",
  "Bottom-up (P_A)"     = "red",
  "Top-down (O_A)"      = "green3",
  "Bottom-up stability" = "red",
  "Top-down stability"  = "green3"
)

method_linetypes <- c(
  "Uncorrected"         = "dashed",
  "B-H"                 = "dotted",
  "B-Y"                 = "dotdash",
  "Bonferroni"          = "longdash",
  "Bottom-up (P_A)"     = "dashed",
  "Top-down (O_A)"      = "dashed",
  "Bottom-up stability" = "solid",
  "Top-down stability"  = "solid"
)

method_shapes <- c(
  "Uncorrected"         = 16,
  "B-H"                 = 17,
  "B-Y"                 = 15,
  "Bonferroni"          = 18,
  "Bottom-up (P_A)"     = 19,
  "Top-down (O_A)"      = 8,
  "Bottom-up stability" = 1,
  "Top-down stability"  = 2
)

base_theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

scale_p1 <- ggplot(
  scale_plot_df[scale_plot_df$Metric == "Accuracy", ],
  aes(x = ScaleVar, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.015, linewidth = 0.5) +
  scale_x_reverse(breaks = scale_var_grid) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(sigma[Y]^2), y = "Accuracy", title = "1. Accuracy") +
  base_theme

scale_p2 <- ggplot(
  scale_plot_df[scale_plot_df$Metric == "Completeness", ],
  aes(x = ScaleVar, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.015, linewidth = 0.5) +
  scale_x_reverse(breaks = scale_var_grid) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(sigma[Y]^2), y = "Completeness", title = "2. Completeness") +
  base_theme

scale_p3 <- ggplot(
  scale_plot_df[scale_plot_df$Metric == "Selected panel size", ],
  aes(x = ScaleVar, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.015, linewidth = 0.5) +
  scale_x_reverse(breaks = scale_var_grid) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(sigma[Y]^2), y = "Selected panel size", title = "3. Selected panel size") +
  base_theme

scale_p4 <- ggplot(
  scale_plot_df[scale_plot_df$Metric == "Exact support recovery", ],
  aes(x = ScaleVar, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.015, linewidth = 0.5) +
  scale_x_reverse(breaks = scale_var_grid) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(sigma[Y]^2), y = "Exact support recovery", title = "4. Exact support recovery") +
  base_theme

scale_p5 <- ggplot(
  scale_plot_df[scale_plot_df$Metric == "Jaccard similarity", ],
  aes(x = ScaleVar, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.015, linewidth = 0.5) +
  scale_x_reverse(breaks = scale_var_grid) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(sigma[Y]^2), y = "Jaccard similarity", title = "5. Jaccard similarity") +
  base_theme

scale_p6 <- ggplot(
  scale_plot_df[scale_plot_df$Metric == "P(false positive)", ],
  aes(x = ScaleVar, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.015, linewidth = 0.5) +
  scale_x_reverse(breaks = scale_var_grid) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  labs(x = expression(sigma[Y]^2), y = "P(false positive)", title = "6. P(false positive)") +
  base_theme

scale_combined_plot <- ((scale_p1 | scale_p2) / (scale_p3 | scale_p4) / (scale_p5 | scale_p6)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, "simulation_metrics_scale_example_6panel.png"),
       plot = scale_combined_plot, width = 12, height = 13.5, dpi = 300)

ggsave(filename = file.path(output_dir, "simulation_metrics_scale_example_6panel.pdf"),
       plot = scale_combined_plot, width = 12, height = 13.5)

print(scale_results_table, row.names = FALSE, right = FALSE)
print(scale_combined_plot)
cat("\nSaved scale-experiment outputs in:\n", output_dir, "\n", sep = "")



# eps-sensitivity experiment for the location problem
# eps here is the already implemented stopping tolerance in Utilities.R


n_mc <- 100L
n_grid <- c(50L, 100L)
d <- 100L
delta_grid <- c(0, 0.5, 1, 3, 5)
eps_grid <- c(0, 1e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2)
eps_labels <- c("0", "1e-5", "1e-4", "2e-4", "5e-4", "1e-3", "2e-3", "5e-3", "1e-2")
relevant_features <- 86:100

methods_eps <- c(
  "Bottom-up (P_A)",
  "Top-down (O_A)"
)

output_dir_eps <- file.path(getwd(), "simulation_outputs_eps_location")
dir.create(output_dir_eps, showWarnings = FALSE, recursive = TRUE)



n_cores <- 8L
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(20260331)

clusterEvalQ(cl, {
  source("Utilities.R")
  NULL
})

clusterExport(
  cl,
  varlist = c("methods_eps", "relevant_features", "d"),
  envir = environment()
)

eps_mc_results_list <- list()
eps_result_counter <- 1L

for (n_idx in seq_along(n_grid)) {
  n_per_group <- n_grid[n_idx]
  
  for (g in seq_along(delta_grid)) {
    delta_mu <- delta_grid[g]
    
    for (e in seq_along(eps_grid)) {
      eps_val <- eps_grid[e]
      
      cat(sprintf("[%s] Starting n = %d, Delta = %s, eps = %s with %d parallel workers\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  n_per_group, delta_mu, eps_labels[e], n_cores))
      flush.console()
      
      mc_out <- foreach(
        b = seq_len(n_mc),
        .combine = rbind,
        .multicombine = TRUE,
        .inorder = FALSE,
        .packages = "stats"
      ) %dorng% {
        
        X <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1),
                    nrow = n_per_group, ncol = d)
        Y <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1),
                    nrow = n_per_group, ncol = d)
        
        Y[, relevant_features] <- Y[, relevant_features, drop = FALSE] + delta_mu
        
        colnames(X) <- paste0("V", seq_len(d))
        colnames(Y) <- paste0("V", seq_len(d))
        
        bu_idx <- bottom_up_pa(X, Y, k = d, p = 1, eps = eps_val)$idx
        td_idx <- top_down_oa(X, Y, q = 1, p = 1, eps = eps_val)$idx
        
        selected_list <- list(
          "Bottom-up (P_A)" = bu_idx,
          "Top-down (O_A)"  = td_idx
        )
        
        if (delta_mu == 0) {
          true_set <- integer(0)
        } else {
          true_set <- relevant_features
        }
        
        out_list <- vector("list", length(methods_eps))
        
        for (i in seq_along(methods_eps)) {
          m <- methods_eps[i]
          A <- sort(unique(as.integer(selected_list[[m]])))
          tp <- length(intersect(A, true_set))
          fp <- length(setdiff(A, true_set))
          union_size <- length(union(A, true_set))
          
          acc_val <- if (length(A) == 0L) 1 else tp / length(A)
          cmp_val <- if (delta_mu == 0) fp else tp / length(true_set)
          exact_val <- as.integer(setequal(A, true_set))
          size_val <- length(A)
          jacc_val <- if (union_size == 0L) 1 else tp / union_size
          fp_val <- as.integer(fp > 0L)
          
          out_list[[i]] <- data.frame(
            Replicate = b,
            NPerGroup = n_per_group,
            Delta = delta_mu,
            Eps = eps_val,
            Method = m,
            Accuracy = acc_val,
            Completeness = cmp_val,
            ExactSupportRecovery = exact_val,
            SelectedPanelSize = size_val,
            JaccardSimilarity = jacc_val,
            PFalsePositive = fp_val,
            stringsAsFactors = FALSE
          )
        }
        
        do.call(rbind, out_list)
      }
      
      eps_mc_results_list[[eps_result_counter]] <- mc_out
      eps_result_counter <- eps_result_counter + 1L
      
      cat(sprintf("[%s] Finished n = %d, Delta = %s, eps = %s\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  n_per_group, delta_mu, eps_labels[e]))
      flush.console()
    }
  }
}

stopCluster(cl)
registerDoSEQ()

eps_mc_results <- do.call(rbind, eps_mc_results_list)
row.names(eps_mc_results) <- NULL

write.csv(
  eps_mc_results,
  file = file.path(output_dir_eps, "eps_location_mc_results.csv"),
  row.names = FALSE
)

# Summary table

summary_rows <- list()
summary_counter <- 1L

for (n_idx in seq_along(n_grid)) {
  n_per_group <- n_grid[n_idx]
  
  for (g in seq_along(delta_grid)) {
    delta_mu <- delta_grid[g]
    
    for (e in seq_along(eps_grid)) {
      eps_val <- eps_grid[e]
      
      for (m in methods_eps) {
        idx <- eps_mc_results$NPerGroup == n_per_group &
          eps_mc_results$Delta == delta_mu &
          eps_mc_results$Eps == eps_val &
          eps_mc_results$Method == m
        
        n_used <- sum(idx)
        
        acc_mean <- mean(eps_mc_results$Accuracy[idx])
        acc_se <- sd(eps_mc_results$Accuracy[idx]) / sqrt(n_used)
        
        cmp_mean <- mean(eps_mc_results$Completeness[idx])
        cmp_se <- sd(eps_mc_results$Completeness[idx]) / sqrt(n_used)
        
        exact_mean <- mean(eps_mc_results$ExactSupportRecovery[idx])
        exact_se <- sd(eps_mc_results$ExactSupportRecovery[idx]) / sqrt(n_used)
        
        size_mean <- mean(eps_mc_results$SelectedPanelSize[idx])
        size_se <- sd(eps_mc_results$SelectedPanelSize[idx]) / sqrt(n_used)
        
        jacc_mean <- mean(eps_mc_results$JaccardSimilarity[idx])
        jacc_se <- sd(eps_mc_results$JaccardSimilarity[idx]) / sqrt(n_used)
        
        fp_mean <- mean(eps_mc_results$PFalsePositive[idx])
        fp_se <- sd(eps_mc_results$PFalsePositive[idx]) / sqrt(n_used)
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          Delta = delta_mu,
          Eps = eps_val,
          Method = m,
          Metric = "Accuracy",
          Mean = acc_mean,
          SE = acc_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          Delta = delta_mu,
          Eps = eps_val,
          Method = m,
          Metric = "Completeness / #FP at null",
          Mean = cmp_mean,
          SE = cmp_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          Delta = delta_mu,
          Eps = eps_val,
          Method = m,
          Metric = "Exact support recovery",
          Mean = exact_mean,
          SE = exact_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          Delta = delta_mu,
          Eps = eps_val,
          Method = m,
          Metric = "Selected panel size",
          Mean = size_mean,
          SE = size_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          Delta = delta_mu,
          Eps = eps_val,
          Method = m,
          Metric = "Jaccard similarity",
          Mean = jacc_mean,
          SE = jacc_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          Delta = delta_mu,
          Eps = eps_val,
          Method = m,
          Metric = "P(false positive)",
          Mean = fp_mean,
          SE = fp_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
      }
    }
  }
}

eps_summary_df <- do.call(rbind, summary_rows)
row.names(eps_summary_df) <- NULL

eps_summary_df$EpsLabel <- factor(
  eps_labels[match(eps_summary_df$Eps, eps_grid)],
  levels = eps_labels
)

eps_plot_grid <- c(min(eps_grid[eps_grid > 0]) / 3, eps_grid[eps_grid > 0])

eps_summary_df$EpsPlot <- eps_plot_grid[match(eps_summary_df$Eps, eps_grid)]
eps_summary_df$DeltaLabel <- factor(
  paste0("Delta = ", eps_summary_df$Delta),
  levels = paste0("Delta = ", delta_grid)
)
eps_summary_df$NLabel <- factor(
  paste0("n = ", eps_summary_df$NPerGroup),
  levels = paste0("n = ", n_grid)
)
eps_summary_df$Method <- factor(
  eps_summary_df$Method,
  levels = methods_eps
)

write.csv(
  eps_summary_df,
  file = file.path(output_dir_eps, "eps_location_summary.csv"),
  row.names = FALSE
)

# Heatmaps


base_theme_hm <- theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

hm_jacc <- ggplot(
  eps_summary_df[eps_summary_df$Metric == "Jaccard similarity", ],
  aes(x = EpsLabel, y = DeltaLabel, fill = Mean)
) +
  geom_tile(color = "grey85") +
  facet_grid(NLabel ~ Method) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  labs(
    x = "eps",
    y = expression(Delta * mu),
    title = "Heatmap: Jaccard similarity"
  ) +
  base_theme_hm

hm_size <- ggplot(
  eps_summary_df[eps_summary_df$Metric == "Selected panel size", ],
  aes(x = EpsLabel, y = DeltaLabel, fill = Mean)
) +
  geom_tile(color = "grey85") +
  facet_grid(NLabel ~ Method) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(
    x = "eps",
    y = expression(Delta * mu),
    title = "Heatmap: Selected panel size"
  ) +
  base_theme_hm

hm_fp <- ggplot(
  eps_summary_df[eps_summary_df$Metric == "P(false positive)", ],
  aes(x = EpsLabel, y = DeltaLabel, fill = Mean)
) +
  geom_tile(color = "grey85") +
  facet_grid(NLabel ~ Method) +
  scale_fill_gradient(low = "white", high = "darkgreen", limits = c(0, 1)) +
  labs(
    x = "eps",
    y = expression(Delta * mu),
    title = "Heatmap: P(false positive)"
  ) +
  base_theme_hm

heatmap_plot <- (hm_jacc / hm_size / hm_fp)

ggsave(
  filename = file.path(output_dir_eps, "eps_location_heatmaps.png"),
  plot = heatmap_plot,
  width = 12,
  height = 14,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir_eps, "eps_location_heatmaps.pdf"),
  plot = heatmap_plot,
  width = 12,
  height = 14
)


# Paired line plots


method_colors_eps <- c(
  "Bottom-up (P_A)" = "red",
  "Top-down (O_A)"  = "green3"
)

method_linetypes_eps <- c(
  "Bottom-up (P_A)" = "solid",
  "Top-down (O_A)"  = "solid"
)

method_shapes_eps <- c(
  "Bottom-up (P_A)" = 16,
  "Top-down (O_A)"  = 17
)

base_theme_lines <- theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

p_acc <- ggplot(
  eps_summary_df[eps_summary_df$Metric == "Accuracy", ],
  aes(x = EpsPlot, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0, linewidth = 0.4) +
  scale_x_log10(breaks = eps_plot_grid, labels = eps_labels) +
  scale_color_manual(values = method_colors_eps) +
  scale_linetype_manual(values = method_linetypes_eps) +
  scale_shape_manual(values = method_shapes_eps) +
  facet_grid(NLabel ~ DeltaLabel) +
  labs(x = "eps", y = "Accuracy", title = "1. Accuracy") +
  base_theme_lines

p_cmp <- ggplot(
  eps_summary_df[eps_summary_df$Metric == "Completeness / #FP at null", ],
  aes(x = EpsPlot, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0, linewidth = 0.4) +
  scale_x_log10(breaks = eps_plot_grid, labels = eps_labels) +
  scale_color_manual(values = method_colors_eps) +
  scale_linetype_manual(values = method_linetypes_eps) +
  scale_shape_manual(values = method_shapes_eps) +
  facet_grid(NLabel ~ DeltaLabel, scales = "free_y") +
  labs(x = "eps", y = "Completeness / #FP at null", title = "2. Completeness / #FP at null") +
  base_theme_lines

p_size <- ggplot(
  eps_summary_df[eps_summary_df$Metric == "Selected panel size", ],
  aes(x = EpsPlot, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0, linewidth = 0.4) +
  scale_x_log10(breaks = eps_plot_grid, labels = eps_labels) +
  scale_color_manual(values = method_colors_eps) +
  scale_linetype_manual(values = method_linetypes_eps) +
  scale_shape_manual(values = method_shapes_eps) +
  facet_grid(NLabel ~ DeltaLabel) +
  labs(x = "eps", y = "Selected panel size", title = "3. Selected panel size") +
  base_theme_lines

p_jacc <- ggplot(
  eps_summary_df[eps_summary_df$Metric == "Jaccard similarity", ],
  aes(x = EpsPlot, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0, linewidth = 0.4) +
  scale_x_log10(breaks = eps_plot_grid, labels = eps_labels) +
  scale_color_manual(values = method_colors_eps) +
  scale_linetype_manual(values = method_linetypes_eps) +
  scale_shape_manual(values = method_shapes_eps) +
  facet_grid(NLabel ~ DeltaLabel) +
  labs(x = "eps", y = "Jaccard similarity", title = "4. Jaccard similarity") +
  base_theme_lines

line_plot <- ((p_acc / p_cmp) / (p_size / p_jacc)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(output_dir_eps, "eps_location_lineplots.png"),
  plot = line_plot,
  width = 18,
  height = 14,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir_eps, "eps_location_lineplots.pdf"),
  plot = line_plot,
  width = 18,
  height = 14
)

print(heatmap_plot)
print(line_plot)

cat("\nSaved eps-sensitivity outputs in:\n", output_dir_eps, "\n", sep = "")


# Scale experiment, epsilon examination

n_mc <- 100L
n_grid <- c(50L, 100L)
d <- 100L
scale_var_grid <- c(1.0, 0.9, 0.8, 0.7, 0.6)
scale_var_labels <- c("1.0", "0.9", "0.8", "0.7", "0.6")
eps_grid <- c(0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
eps_labels <- c("0", "1e-12", "1e-10", "1e-8", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2")
relevant_features <- 86:100

methods_eps <- c(
  "Bottom-up (P_A)",
  "Top-down (O_A)"
)

output_dir_eps_scale <- file.path(getwd(), "simulation_outputs_eps_scale")
dir.create(output_dir_eps_scale, showWarnings = FALSE, recursive = TRUE)

n_cores <- 8L
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(20260401)

clusterEvalQ(cl, {
  source("Utilities.R")
  NULL
})

clusterExport(
  cl,
  varlist = c("methods_eps", "relevant_features", "d"),
  envir = environment()
)

eps_scale_mc_results_list <- list()
eps_scale_result_counter <- 1L

for (n_idx in seq_along(n_grid)) {
  n_per_group <- n_grid[n_idx]
  
  for (s in seq_along(scale_var_grid)) {
    scale_var <- scale_var_grid[s]
    
    for (e in seq_along(eps_grid)) {
      eps_val <- eps_grid[e]
      
      cat(sprintf("[%s] Starting n = %d, sigmaY^2 = %s, eps = %s with %d parallel workers\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  n_per_group, scale_var_labels[s], eps_labels[e], n_cores))
      flush.console()
      
      mc_out <- foreach(
        b = seq_len(n_mc),
        .combine = rbind,
        .multicombine = TRUE,
        .inorder = FALSE,
        .packages = "stats"
      ) %dorng% {
        
        X <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1),
                    nrow = n_per_group, ncol = d)
        Y <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1),
                    nrow = n_per_group, ncol = d)
        
        Y[, relevant_features] <- sqrt(scale_var) * Y[, relevant_features, drop = FALSE]
        
        colnames(X) <- paste0("V", seq_len(d))
        colnames(Y) <- paste0("V", seq_len(d))
        
        bu_idx <- bottom_up_pa(X, Y, k = d, p = 1, eps = eps_val)$idx
        td_idx <- top_down_oa(X, Y, q = 1, p = 1, eps = eps_val)$idx
        
        selected_list <- list(
          "Bottom-up (P_A)" = bu_idx,
          "Top-down (O_A)"  = td_idx
        )
        
        if (scale_var == 1) {
          true_set <- integer(0)
        } else {
          true_set <- relevant_features
        }
        
        out_list <- vector("list", length(methods_eps))
        
        for (i in seq_along(methods_eps)) {
          m <- methods_eps[i]
          A <- sort(unique(as.integer(selected_list[[m]])))
          tp <- length(intersect(A, true_set))
          fp <- length(setdiff(A, true_set))
          union_size <- length(union(A, true_set))
          
          acc_val <- if (length(A) == 0L) 1 else tp / length(A)
          cmp_val <- if (scale_var == 1) fp else tp / length(true_set)
          exact_val <- as.integer(setequal(A, true_set))
          size_val <- length(A)
          jacc_val <- if (union_size == 0L) 1 else tp / union_size
          fp_val <- as.integer(fp > 0L)
          
          out_list[[i]] <- data.frame(
            Replicate = b,
            NPerGroup = n_per_group,
            ScaleVar = scale_var,
            Eps = eps_val,
            Method = m,
            Accuracy = acc_val,
            Completeness = cmp_val,
            ExactSupportRecovery = exact_val,
            SelectedPanelSize = size_val,
            JaccardSimilarity = jacc_val,
            PFalsePositive = fp_val,
            stringsAsFactors = FALSE
          )
        }
        
        do.call(rbind, out_list)
      }
      
      eps_scale_mc_results_list[[eps_scale_result_counter]] <- mc_out
      eps_scale_result_counter <- eps_scale_result_counter + 1L
      
      cat(sprintf("[%s] Finished n = %d, sigmaY^2 = %s, eps = %s\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  n_per_group, scale_var_labels[s], eps_labels[e]))
      flush.console()
    }
  }
}

stopCluster(cl)
registerDoSEQ()

eps_scale_mc_results <- do.call(rbind, eps_scale_mc_results_list)
row.names(eps_scale_mc_results) <- NULL

write.csv(
  eps_scale_mc_results,
  file = file.path(output_dir_eps_scale, "eps_scale_mc_results.csv"),
  row.names = FALSE
)


# Summary table

summary_rows <- list()
summary_counter <- 1L

for (n_idx in seq_along(n_grid)) {
  n_per_group <- n_grid[n_idx]
  
  for (s in seq_along(scale_var_grid)) {
    scale_var <- scale_var_grid[s]
    
    for (e in seq_along(eps_grid)) {
      eps_val <- eps_grid[e]
      
      for (m in methods_eps) {
        idx <- eps_scale_mc_results$NPerGroup == n_per_group &
          eps_scale_mc_results$ScaleVar == scale_var &
          eps_scale_mc_results$Eps == eps_val &
          eps_scale_mc_results$Method == m
        
        n_used <- sum(idx)
        
        acc_mean <- mean(eps_scale_mc_results$Accuracy[idx])
        acc_se <- sd(eps_scale_mc_results$Accuracy[idx]) / sqrt(n_used)
        
        cmp_mean <- mean(eps_scale_mc_results$Completeness[idx])
        cmp_se <- sd(eps_scale_mc_results$Completeness[idx]) / sqrt(n_used)
        
        exact_mean <- mean(eps_scale_mc_results$ExactSupportRecovery[idx])
        exact_se <- sd(eps_scale_mc_results$ExactSupportRecovery[idx]) / sqrt(n_used)
        
        size_mean <- mean(eps_scale_mc_results$SelectedPanelSize[idx])
        size_se <- sd(eps_scale_mc_results$SelectedPanelSize[idx]) / sqrt(n_used)
        
        jacc_mean <- mean(eps_scale_mc_results$JaccardSimilarity[idx])
        jacc_se <- sd(eps_scale_mc_results$JaccardSimilarity[idx]) / sqrt(n_used)
        
        fp_mean <- mean(eps_scale_mc_results$PFalsePositive[idx])
        fp_se <- sd(eps_scale_mc_results$PFalsePositive[idx]) / sqrt(n_used)
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          ScaleVar = scale_var,
          Eps = eps_val,
          Method = m,
          Metric = "Accuracy",
          Mean = acc_mean,
          SE = acc_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          ScaleVar = scale_var,
          Eps = eps_val,
          Method = m,
          Metric = "Completeness / #FP at null",
          Mean = cmp_mean,
          SE = cmp_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          ScaleVar = scale_var,
          Eps = eps_val,
          Method = m,
          Metric = "Exact support recovery",
          Mean = exact_mean,
          SE = exact_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          ScaleVar = scale_var,
          Eps = eps_val,
          Method = m,
          Metric = "Selected panel size",
          Mean = size_mean,
          SE = size_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          ScaleVar = scale_var,
          Eps = eps_val,
          Method = m,
          Metric = "Jaccard similarity",
          Mean = jacc_mean,
          SE = jacc_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
        
        summary_rows[[summary_counter]] <- data.frame(
          NPerGroup = n_per_group,
          ScaleVar = scale_var,
          Eps = eps_val,
          Method = m,
          Metric = "P(false positive)",
          Mean = fp_mean,
          SE = fp_se,
          stringsAsFactors = FALSE
        )
        summary_counter <- summary_counter + 1L
      }
    }
  }
}

eps_scale_summary_df <- do.call(rbind, summary_rows)
row.names(eps_scale_summary_df) <- NULL

eps_scale_summary_df$EpsLabel <- factor(
  eps_labels[match(eps_scale_summary_df$Eps, eps_grid)],
  levels = eps_labels
)

eps_plot_grid <- c(min(eps_grid[eps_grid > 0]) / 3, eps_grid[eps_grid > 0])

eps_scale_summary_df$EpsPlot <- eps_plot_grid[match(eps_scale_summary_df$Eps, eps_grid)]
eps_scale_summary_df$ScaleLabel <- factor(
  paste0("Var = ", scale_var_labels[match(eps_scale_summary_df$ScaleVar, scale_var_grid)]),
  levels = paste0("Var = ", scale_var_labels)
)
eps_scale_summary_df$NLabel <- factor(
  paste0("n = ", eps_scale_summary_df$NPerGroup),
  levels = paste0("n = ", n_grid)
)
eps_scale_summary_df$Method <- factor(
  eps_scale_summary_df$Method,
  levels = methods_eps
)

write.csv(
  eps_scale_summary_df,
  file = file.path(output_dir_eps_scale, "eps_scale_summary.csv"),
  row.names = FALSE
)



# Scale/epsilon plots

base_theme_hm <- theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

hm_jacc <- ggplot(
  eps_scale_summary_df[eps_scale_summary_df$Metric == "Jaccard similarity", ],
  aes(x = EpsLabel, y = ScaleLabel, fill = Mean)
) +
  geom_tile(color = "grey85") +
  facet_grid(NLabel ~ Method) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  labs(
    x = "eps",
    y = expression(sigma[Y]^2),
    title = "Heatmap: Jaccard similarity"
  ) +
  base_theme_hm

hm_size <- ggplot(
  eps_scale_summary_df[eps_scale_summary_df$Metric == "Selected panel size", ],
  aes(x = EpsLabel, y = ScaleLabel, fill = Mean)
) +
  geom_tile(color = "grey85") +
  facet_grid(NLabel ~ Method) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(
    x = "eps",
    y = expression(sigma[Y]^2),
    title = "Heatmap: Selected panel size"
  ) +
  base_theme_hm

hm_fp <- ggplot(
  eps_scale_summary_df[eps_scale_summary_df$Metric == "P(false positive)", ],
  aes(x = EpsLabel, y = ScaleLabel, fill = Mean)
) +
  geom_tile(color = "grey85") +
  facet_grid(NLabel ~ Method) +
  scale_fill_gradient(low = "white", high = "darkgreen", limits = c(0, 1)) +
  labs(
    x = "eps",
    y = expression(sigma[Y]^2),
    title = "Heatmap: P(false positive)"
  ) +
  base_theme_hm

heatmap_plot <- (hm_jacc / hm_size / hm_fp)

ggsave(
  filename = file.path(output_dir_eps_scale, "eps_scale_heatmaps.png"),
  plot = heatmap_plot,
  width = 12,
  height = 14,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir_eps_scale, "eps_scale_heatmaps.pdf"),
  plot = heatmap_plot,
  width = 12,
  height = 14
)

method_colors_eps <- c(
  "Bottom-up (P_A)" = "red",
  "Top-down (O_A)"  = "green3"
)

method_linetypes_eps <- c(
  "Bottom-up (P_A)" = "solid",
  "Top-down (O_A)"  = "solid"
)

method_shapes_eps <- c(
  "Bottom-up (P_A)" = 16,
  "Top-down (O_A)"  = 17
)

base_theme_lines <- theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

p_acc <- ggplot(
  eps_scale_summary_df[eps_scale_summary_df$Metric == "Accuracy", ],
  aes(x = EpsPlot, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0, linewidth = 0.4) +
  scale_x_log10(breaks = eps_plot_grid, labels = eps_labels) +
  scale_color_manual(values = method_colors_eps) +
  scale_linetype_manual(values = method_linetypes_eps) +
  scale_shape_manual(values = method_shapes_eps) +
  facet_grid(NLabel ~ ScaleLabel) +
  labs(x = "eps", y = "Accuracy", title = "1. Accuracy") +
  base_theme_lines

p_cmp <- ggplot(
  eps_scale_summary_df[eps_scale_summary_df$Metric == "Completeness / #FP at null", ],
  aes(x = EpsPlot, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0, linewidth = 0.4) +
  scale_x_log10(breaks = eps_plot_grid, labels = eps_labels) +
  scale_color_manual(values = method_colors_eps) +
  scale_linetype_manual(values = method_linetypes_eps) +
  scale_shape_manual(values = method_shapes_eps) +
  facet_grid(NLabel ~ ScaleLabel, scales = "free_y") +
  labs(x = "eps", y = "Completeness / #FP at null", title = "2. Completeness / #FP at null") +
  base_theme_lines

p_size <- ggplot(
  eps_scale_summary_df[eps_scale_summary_df$Metric == "Selected panel size", ],
  aes(x = EpsPlot, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0, linewidth = 0.4) +
  scale_x_log10(breaks = eps_plot_grid, labels = eps_labels) +
  scale_color_manual(values = method_colors_eps) +
  scale_linetype_manual(values = method_linetypes_eps) +
  scale_shape_manual(values = method_shapes_eps) +
  facet_grid(NLabel ~ ScaleLabel) +
  labs(x = "eps", y = "Selected panel size", title = "3. Selected panel size") +
  base_theme_lines

p_jacc <- ggplot(
  eps_scale_summary_df[eps_scale_summary_df$Metric == "Jaccard similarity", ],
  aes(x = EpsPlot, y = Mean, group = Method, color = Method, linetype = Method, shape = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0, linewidth = 0.4) +
  scale_x_log10(breaks = eps_plot_grid, labels = eps_labels) +
  scale_color_manual(values = method_colors_eps) +
  scale_linetype_manual(values = method_linetypes_eps) +
  scale_shape_manual(values = method_shapes_eps) +
  facet_grid(NLabel ~ ScaleLabel) +
  labs(x = "eps", y = "Jaccard similarity", title = "4. Jaccard similarity") +
  base_theme_lines

line_plot <- ((p_acc / p_cmp) / (p_size / p_jacc)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(output_dir_eps_scale, "eps_scale_lineplots.png"),
  plot = line_plot,
  width = 18,
  height = 14,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir_eps_scale, "eps_scale_lineplots.pdf"),
  plot = line_plot,
  width = 18,
  height = 14
)

print(heatmap_plot)
print(line_plot)

# Examination of objectives (saturation phenomenon)

n_mc <- 200L
n_per_group <- 50L
d <- 100L
delta_mu <- 3
p_val <- 1

relevant_features <- 86:100
irrelevant_features <- setdiff(seq_len(d), relevant_features)

signal_grid <- c(0L, 1L, 3L, 5L, 10L, 15L)
noise_grid  <- c(0L, 1L, 2L, 3L, 5L, 8L, 10L, 15L, 20L, 30L, 40L, 50L)

grid_df <- expand.grid(
  Replicate = seq_len(n_mc),
  Signal = signal_grid,
  Noise = noise_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

res_list <- vector("list", nrow(grid_df))

for (b in seq_len(n_mc)) {
  X <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1), nrow = n_per_group, ncol = d)
  Y <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1), nrow = n_per_group, ncol = d)
  Y[, relevant_features] <- Y[, relevant_features, drop = FALSE] + delta_mu
  
  colnames(X) <- paste0("V", seq_len(d))
  colnames(Y) <- paste0("V", seq_len(d))
  
  idx_b <- which(grid_df$Replicate == b)
  
  for (ii in idx_b) {
    s <- grid_df$Signal[ii]
    r <- grid_df$Noise[ii]
    
    A_signal <- if (s > 0L) relevant_features[seq_len(s)] else integer(0)
    A_noise  <- if (r > 0L) irrelevant_features[seq_len(r)] else integer(0)
    A <- c(A_signal, A_noise)
    
    res_list[[ii]] <- data.frame(
      Replicate = b,
      Signal = s,
      Noise = r,
      PanelSize = length(A),
      P_A = pa_hat(X, Y, A = A, p = p_val),
      O_A = oa_hat(X, Y, A = A, p = p_val),
      stringsAsFactors = FALSE
    )
  }
  
  cat(sprintf("[%s] finished replicate %d / %d\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), b, n_mc))
  flush.console()
}

res_df <- do.call(rbind, res_list)
row.names(res_df) <- NULL

pa_mean <- aggregate(P_A ~ Signal + Noise, data = res_df, FUN = mean)
pa_sd   <- aggregate(P_A ~ Signal + Noise, data = res_df, FUN = sd)
oa_mean <- aggregate(O_A ~ Signal + Noise, data = res_df, FUN = mean)
oa_sd   <- aggregate(O_A ~ Signal + Noise, data = res_df, FUN = sd)

summary_df <- merge(pa_mean, pa_sd, by = c("Signal", "Noise"), suffixes = c("_Mean", "_SD"))
summary_df <- merge(summary_df, oa_mean, by = c("Signal", "Noise"))
summary_df <- merge(summary_df, oa_sd, by = c("Signal", "Noise"), suffixes = c("_Mean", "_SD"))

names(summary_df)[names(summary_df) == "O_A_Mean"] <- "O_A_Mean"
names(summary_df)[names(summary_df) == "O_A_SD"]   <- "O_A_SD"

summary_df$P_A_SE <- summary_df$P_A_SD / sqrt(n_mc)
summary_df$O_A_SE <- summary_df$O_A_SD / sqrt(n_mc)

summary_df$SignalF <- factor(summary_df$Signal, levels = signal_grid)
summary_df$NoiseF  <- factor(summary_df$Noise, levels = noise_grid)

long_plot_df <- rbind(
  data.frame(
    Signal = summary_df$Signal,
    Noise = summary_df$Noise,
    SignalF = summary_df$SignalF,
    NoiseF = summary_df$NoiseF,
    Objective = "P_A",
    Mean = summary_df$P_A_Mean,
    SE = summary_df$P_A_SE,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Signal = summary_df$Signal,
    Noise = summary_df$Noise,
    SignalF = summary_df$SignalF,
    NoiseF = summary_df$NoiseF,
    Objective = "O_A",
    Mean = summary_df$O_A_Mean,
    SE = summary_df$O_A_SE,
    stringsAsFactors = FALSE
  )
)

write.csv(res_df, "objective_contamination_raw.csv", row.names = FALSE)
write.csv(summary_df, "objective_contamination_summary.csv", row.names = FALSE)

line_plot_df <- subset(long_plot_df, !(Signal == 0 & Noise == 0))

pa_heat_df <- long_plot_df[long_plot_df$Objective == "P_A", ]
oa_heat_df <- long_plot_df[long_plot_df$Objective == "O_A", ]

pa_limits <- range(
  pa_heat_df$Mean[!(pa_heat_df$Signal == 0 & pa_heat_df$Noise == 0)],
  na.rm = TRUE
)

p_pa <- ggplot(
  line_plot_df[line_plot_df$Objective == "P_A", ],
  aes(x = Noise, y = Mean, group = SignalF, color = SignalF)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.8, linewidth = 0.45) +
  scale_x_continuous(breaks = noise_grid) +
  labs(
    x = "Number of noise features",
    y = expression(P[A]),
    color = "Signal features",
    title = expression("Contamination curves for " * P[A])
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_oa <- ggplot(
  line_plot_df[line_plot_df$Objective == "O_A", ],
  aes(x = Noise, y = Mean, group = SignalF, color = SignalF)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.8, linewidth = 0.45) +
  scale_x_continuous(breaks = noise_grid) +
  labs(
    x = "Number of noise features",
    y = expression(O[A]),
    color = "Signal features",
    title = expression("Contamination curves for " * O[A])
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

pa_breaks <- signif(seq(pa_limits[1], pa_limits[2], length.out = 3), 4)

heat_pa <- ggplot(
  pa_heat_df,
  aes(x = NoiseF, y = SignalF, fill = Mean)
) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = hcl.colors(100, "Inferno"),
    limits = pa_limits,
    oob = scales::squish,
    breaks = pa_breaks,
    labels = function(x) formatC(x, digits = 4, format = "f")
  ) +
  labs(
    x = "Number of noise features",
    y = "Number of signal features",
    fill = expression(P[A]),
    title = expression("Heatmap of mean " * P[A])
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


heat_oa <- ggplot(
  oa_heat_df,
  aes(x = NoiseF, y = SignalF, fill = Mean)
) +
  geom_tile(color = "white") +
  labs(
    x = "Number of noise features",
    y = "Number of signal features",
    fill = expression(O[A]),
    title = expression("Heatmap of mean " * O[A])
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

combined_plot <- ((p_pa | p_oa) / (heat_pa | heat_oa)) +
  plot_annotation(
    title = "Objective behaviour under mixtures of signal and noise",
    subtitle = paste0(
      "Location model: n = ", n_per_group,
      " per group, d = ", d,
      ", Delta_mu = ", delta_mu,
      ", Monte Carlo = ", n_mc,
      ", p = ", p_val
    )
  ) &
  theme(legend.position = "bottom")

ggsave(
  filename = "objective_contamination_plots.png",
  plot = combined_plot,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  filename = "objective_contamination_plots.pdf",
  plot = combined_plot,
  width = 14,
  height = 10
)

print(combined_plot)

n_mc <- 200L
n_per_group <- 50L
d <- 100L
delta_mu <- 1
p_val <- 1

relevant_features <- 86:100
irrelevant_features <- setdiff(seq_len(d), relevant_features)

signal_grid <- c(0L, 1L, 3L, 5L, 10L, 15L)
noise_grid  <- c(0L, 1L, 2L, 3L, 5L, 8L, 10L, 15L, 20L, 30L, 40L, 50L)

grid_df <- expand.grid(
  Replicate = seq_len(n_mc),
  Signal = signal_grid,
  Noise = noise_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

res_list <- vector("list", nrow(grid_df))

for (b in seq_len(n_mc)) {
  X <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1), nrow = n_per_group, ncol = d)
  Y <- matrix(rnorm(n_per_group * d, mean = 0, sd = 1), nrow = n_per_group, ncol = d)
  Y[, relevant_features] <- Y[, relevant_features, drop = FALSE] + delta_mu
  
  colnames(X) <- paste0("V", seq_len(d))
  colnames(Y) <- paste0("V", seq_len(d))
  
  idx_b <- which(grid_df$Replicate == b)
  
  for (ii in idx_b) {
    s <- grid_df$Signal[ii]
    r <- grid_df$Noise[ii]
    
    A_signal <- if (s > 0L) relevant_features[seq_len(s)] else integer(0)
    A_noise  <- if (r > 0L) irrelevant_features[seq_len(r)] else integer(0)
    A <- c(A_signal, A_noise)
    
    res_list[[ii]] <- data.frame(
      Replicate = b,
      Signal = s,
      Noise = r,
      PanelSize = length(A),
      P_A = pa_hat(X, Y, A = A, p = p_val),
      O_A = oa_hat(X, Y, A = A, p = p_val),
      stringsAsFactors = FALSE
    )
  }
  
  cat(sprintf("[%s] finished replicate %d / %d\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), b, n_mc))
  flush.console()
}

res_df <- do.call(rbind, res_list)
row.names(res_df) <- NULL

pa_mean <- aggregate(P_A ~ Signal + Noise, data = res_df, FUN = mean)
pa_sd   <- aggregate(P_A ~ Signal + Noise, data = res_df, FUN = sd)
oa_mean <- aggregate(O_A ~ Signal + Noise, data = res_df, FUN = mean)
oa_sd   <- aggregate(O_A ~ Signal + Noise, data = res_df, FUN = sd)

summary_df <- merge(pa_mean, pa_sd, by = c("Signal", "Noise"), suffixes = c("_Mean", "_SD"))
summary_df <- merge(summary_df, oa_mean, by = c("Signal", "Noise"))
summary_df <- merge(summary_df, oa_sd, by = c("Signal", "Noise"), suffixes = c("_Mean", "_SD"))

names(summary_df)[names(summary_df) == "O_A_Mean"] <- "O_A_Mean"
names(summary_df)[names(summary_df) == "O_A_SD"]   <- "O_A_SD"

summary_df$P_A_SE <- summary_df$P_A_SD / sqrt(n_mc)
summary_df$O_A_SE <- summary_df$O_A_SD / sqrt(n_mc)

summary_df$SignalF <- factor(summary_df$Signal, levels = signal_grid)
summary_df$NoiseF  <- factor(summary_df$Noise, levels = noise_grid)

long_plot_df <- rbind(
  data.frame(
    Signal = summary_df$Signal,
    Noise = summary_df$Noise,
    SignalF = summary_df$SignalF,
    NoiseF = summary_df$NoiseF,
    Objective = "P_A",
    Mean = summary_df$P_A_Mean,
    SE = summary_df$P_A_SE,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Signal = summary_df$Signal,
    Noise = summary_df$Noise,
    SignalF = summary_df$SignalF,
    NoiseF = summary_df$NoiseF,
    Objective = "O_A",
    Mean = summary_df$O_A_Mean,
    SE = summary_df$O_A_SE,
    stringsAsFactors = FALSE
  )
)

write.csv(res_df, "objective_contamination_delta1_raw.csv", row.names = FALSE)
write.csv(summary_df, "objective_contamination_delta1_summary.csv", row.names = FALSE)

line_plot_df <- subset(long_plot_df, !(Signal == 0 & Noise == 0))

pa_heat_df <- long_plot_df[long_plot_df$Objective == "P_A", ]
oa_heat_df <- long_plot_df[long_plot_df$Objective == "O_A", ]

pa_limits <- range(
  pa_heat_df$Mean[!(pa_heat_df$Signal == 0 & pa_heat_df$Noise == 0)],
  na.rm = TRUE
)

p_pa <- ggplot(
  line_plot_df[line_plot_df$Objective == "P_A", ],
  aes(x = Noise, y = Mean, group = SignalF, color = SignalF)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.8, linewidth = 0.45) +
  scale_x_continuous(breaks = noise_grid) +
  labs(
    x = "Number of noise features",
    y = expression(P[A]),
    color = "Signal features",
    title = expression("Contamination curves for " * P[A])
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_oa <- ggplot(
  line_plot_df[line_plot_df$Objective == "O_A", ],
  aes(x = Noise, y = Mean, group = SignalF, color = SignalF)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.8, linewidth = 0.45) +
  scale_x_continuous(breaks = noise_grid) +
  labs(
    x = "Number of noise features",
    y = expression(O[A]),
    color = "Signal features",
    title = expression("Contamination curves for " * O[A])
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

pa_breaks <- signif(seq(pa_limits[1], pa_limits[2], length.out = 3), 4)

heat_pa <- ggplot(
  pa_heat_df,
  aes(x = NoiseF, y = SignalF, fill = Mean)
) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = hcl.colors(100, "Inferno"),
    limits = pa_limits,
    oob = scales::squish,
    breaks = pa_breaks,
    labels = function(x) formatC(x, digits = 4, format = "f")
  ) +
  labs(
    x = "Number of noise features",
    y = "Number of signal features",
    fill = expression(P[A]),
    title = expression("Heatmap of mean " * P[A])
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


heat_oa <- ggplot(
  oa_heat_df,
  aes(x = NoiseF, y = SignalF, fill = Mean)
) +
  geom_tile(color = "white") +
  labs(
    x = "Number of noise features",
    y = "Number of signal features",
    fill = expression(O[A]),
    title = expression("Heatmap of mean " * O[A])
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

combined_plot <- ((p_pa | p_oa) / (heat_pa | heat_oa)) +
  plot_annotation(
    title = "Objective behaviour under mixtures of signal and noise",
    subtitle = paste0(
      "Location model: n = ", n_per_group,
      " per group, d = ", d,
      ", Delta_mu = ", delta_mu,
      ", Monte Carlo = ", n_mc,
      ", p = ", p_val
    )
  ) &
  theme(legend.position = "bottom")

ggsave(
  filename = "objective_contamination_delta1_plots.png",
  plot = combined_plot,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  filename = "objective_contamination_delta1_plots.pdf",
  plot = combined_plot,
  width = 14,
  height = 10
)

print(combined_plot)