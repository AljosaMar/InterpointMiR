
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)

# -------------------- Distance-based selection --------------------

pa_hat <- function(X, Y, A = NULL, p = 1) {
  if (is.null(A)) A <- seq_len(ncol(X))
  if (length(A) == 0L) return(3 / 4)
  
  safe_var <- function(x) {
    v <- stats::var(x)
    if (is.na(v) || !is.finite(v)) 0 else v
  }
  
  f <- function(v) {
    d <- abs(outer(v, v, "-"))^p
    u <- d[upper.tri(d)]
    c(mean = mean(u), var = safe_var(u))
  }
  
  g <- function(a, b) {
    d <- abs(outer(a, b, "-"))^p
    u <- as.vector(d)
    c(mean = mean(u), var = safe_var(u))
  }
  
  mx <- vx <- my <- vy <- mxy <- vxy <- 0
  
  for (j in A) {
    t1 <- f(X[, j]);      mx  <- mx  + t1[1]; vx  <- vx  + t1[2]
    t2 <- f(Y[, j]);      my  <- my  + t2[1]; vy  <- vy  + t2[2]
    t3 <- g(X[, j], Y[, j]); mxy <- mxy + t3[1]; vxy <- vxy + t3[2]
  }
  
  sd1 <- sqrt(pmax(vxy + vx, .Machine$double.eps))
  sd2 <- sqrt(pmax(vxy + vy, .Machine$double.eps))
  sd3 <- sqrt(pmax(vy  + vx, .Machine$double.eps))
  
  # P1 = P(dX <= dXY)
  P1 <- stats::pnorm(0, mean = mxy - mx, sd = sd1, lower.tail = FALSE)
  
  # P2 = P(dY <= dXY)
  P2 <- stats::pnorm(0, mean = mxy - my, sd = sd2, lower.tail = FALSE)
  
  # P3 = P(dX <= dY)
  P3 <- stats::pnorm(0, mean = my - mx, sd = sd3, lower.tail = FALSE)
  
  sum((c(P1, P2, P3) - 0.5)^2)
}

bottom_up_pa <- function(X, Y, k = 20, p = 1, eps = 0) {
  d <- ncol(X)
  
  if (is.null(colnames(X)) || any(colnames(X) == "")) {
    colnames(X) <- paste0("V", seq_len(d))
    colnames(Y) <- colnames(X)
  }
  
  S <- integer(0); path <- numeric(0)
  for (t in seq_len(min(k, d))) {
    cand <- setdiff(seq_len(d), S)
    sc <- sapply(cand, function(j) pa_hat(X, Y, c(S, j), p))
    j <- cand[which.max(sc)]; best <- max(sc)
    if (length(path) && best - tail(path, 1) < eps) break
    S <- c(S, j); path <- c(path, best)
  }
  list(idx = S, feat = colnames(X)[S], pa_path = path)
}

stability_bottom_up <- function(X, Y,
                                B = 100, sub = 0.5, k = 20, p = 1,
                                pi_thr = 0.6, eps = 0, min_n = 3L, seed = NULL,
                                prefix = "", progress_every = 10L) {
  if (!is.null(seed)) set.seed(seed)
  d <- ncol(X)
  
  if (is.null(colnames(X)) || any(colnames(X) == "")) {
    colnames(X) <- paste0("V", seq_len(d))
  }
  if (is.null(colnames(Y)) || any(colnames(Y) == "")) {
    colnames(Y) <- colnames(X)
  }
  
  nx <- nrow(X); ny <- nrow(Y)
  cnt <- integer(d)
  
  for (b in seq_len(B)) {
    nx_sub <- max(min_n, floor(sub * nx))
    ny_sub <- max(min_n, floor(sub * ny))
    ix <- sample.int(nx, nx_sub, replace = FALSE)
    iy <- sample.int(ny, ny_sub, replace = FALSE)
    
    S <- bottom_up_pa(X[ix, , drop = FALSE], Y[iy, , drop = FALSE], k = k, p = p, eps = eps)$idx
    if (length(S)) cnt[S] <- cnt[S] + 1L
    
    if (progress_every > 0 && (b %% progress_every == 0 || b == B)) {
      cat(sprintf("%sDistance stability: %d/%d (%.0f%%)\n", prefix, b, B, 100 * b / B))
      flush.console()
    }
  }
  
  pi <- cnt / B
  names(pi) <- colnames(X)
  
  ord <- order(pi, decreasing = TRUE)
  sel <- which(pi >= pi_thr)
  if (length(sel) == 0) sel <- head(ord, k)
  
  list(pi = pi, order = ord, selected_idx = sel, selected_feat = names(pi)[sel])
}


# -------------------- Top-down distance-based selection --------------------

.safe_var <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 2L) return(0)
  v <- stats::var(x)
  if (is.na(v) || !is.finite(v)) 0 else v
}

# Per-feature moments of |diff|^p for within/between distances.
# Returns vectors of length d with means/variances for each feature.
.feature_moments <- function(X, Y, p = 1) {
  d <- ncol(X)
  mu_dx <- mu_dy <- mu_dxy <- numeric(d)
  var_dx <- var_dy <- var_dxy <- numeric(d)
  
  f_within <- function(v) {
    D <- abs(outer(v, v, "-"))^p
    u <- D[upper.tri(D)]
    c(mean = mean(u), var = .safe_var(u))
  }
  
  f_between <- function(a, b) {
    D <- abs(outer(a, b, "-"))^p
    u <- as.vector(D)
    c(mean = mean(u), var = .safe_var(u))
  }
  
  for (j in seq_len(d)) {
    t1 <- f_within(X[, j]); mu_dx[j]  <- t1[1]; var_dx[j]  <- t1[2]
    t2 <- f_within(Y[, j]); mu_dy[j]  <- t2[1]; var_dy[j]  <- t2[2]
    t3 <- f_between(X[, j], Y[, j]); mu_dxy[j] <- t3[1]; var_dxy[j] <- t3[2]
  }
  
  list(mu_dx = mu_dx, var_dx = var_dx,
       mu_dy = mu_dy, var_dy = var_dy,
       mu_dxy = mu_dxy, var_dxy = var_dxy)
}

# Numerically evaluate the overlap between
#   0.5 N(mu_dx, var_dx) + 0.5 N(mu_dy, var_dy)
# and
#   N(mu_dxy, var_dxy).

.overlap_normal_mixture <- function(mu_dx, var_dx,
                                    mu_dy, var_dy,
                                    mu_dxy, var_dxy,
                                    rel_tol = 1e-7,
                                    abs_tol = 1e-9,
                                    subdivisions = 200L) {
  pars <- c(mu_dx, var_dx, mu_dy, var_dy, mu_dxy, var_dxy)
  if (length(pars) != 6L || any(!is.finite(pars))) {
    stop("All normal-approximation means and variances must be finite scalars.")
  }
  
  # Variances estimated from finite samples can be exactly zero. A common,
  # scale-adaptive floor gives a well-defined normal approximation while
  # preserving equality when all three approximations coincide.
  variance_scale <- max(1, abs(mu_dx)^2, abs(mu_dy)^2, abs(mu_dxy)^2,
                        var_dx, var_dy, var_dxy)
  variance_floor <- .Machine$double.eps * variance_scale
  var_dx  <- max(var_dx,  variance_floor)
  var_dy  <- max(var_dy,  variance_floor)
  var_dxy <- max(var_dxy, variance_floor)
  
  sd_dx  <- sqrt(var_dx)
  sd_dy  <- sqrt(var_dy)
  sd_dxy <- sqrt(var_dxy)
  
  # Under the null, the three normal approximations coincide and the overlap is 1.
  same_dx <- isTRUE(all.equal(mu_dx, mu_dxy, tolerance = 1e-12)) &&
    isTRUE(all.equal(var_dx, var_dxy, tolerance = 1e-12))
  same_dy <- isTRUE(all.equal(mu_dy, mu_dxy, tolerance = 1e-12)) &&
    isTRUE(all.equal(var_dy, var_dxy, tolerance = 1e-12))
  if (same_dx && same_dy) return(1)
  
  log_half <- log(0.5)
  logspace_add <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  
  integrand <- function(t) {
    z <- mu_dxy + sd_dxy * t
    log_fx <- stats::dnorm(z, mean = mu_dx,  sd = sd_dx,  log = TRUE)
    log_fy <- stats::dnorm(z, mean = mu_dy,  sd = sd_dy,  log = TRUE)
    log_fb <- stats::dnorm(z, mean = mu_dxy, sd = sd_dxy, log = TRUE)
    log_mix <- logspace_add(log_half + log_fx, log_half + log_fy)
    ratio_truncated <- exp(pmin(log_mix - log_fb, 0))
    stats::dnorm(t) * ratio_truncated
  }
  
  # The omitted standard-normal mass outside [-8, 8] is below 1.3e-15.
  ans <- suppressWarnings(tryCatch(
    stats::integrate(
      integrand,
      lower = -8,
      upper = 8,
      subdivisions = subdivisions,
      rel.tol = rel_tol,
      abs.tol = abs_tol,
      stop.on.error = FALSE
    ),
    error = function(e) NULL
  ))
  
  if (!is.null(ans) && is.finite(ans$value)) {
    return(pmin(pmax(as.numeric(ans$value), 0), 1))
  }
  
  # Deterministic fallback if adaptive integration fails.
  grid <- seq(-8, 8, length.out = 4001L)
  vals <- integrand(grid)
  h <- grid[2L] - grid[1L]
  overlap <- h * (sum(vals) - 0.5 * vals[1L] - 0.5 * vals[length(vals)])
  pmin(pmax(as.numeric(overlap), 0), 1)
}

oa_hat <- function(X, Y, A = NULL, p = 1) {
  if (is.null(A)) A <- seq_len(ncol(X))
  if (length(A) == 0L) return(1)
  
  mom <- .feature_moments(X, Y, p = p)
  
  .overlap_normal_mixture(
    mu_dx  = sum(mom$mu_dx[A]),
    var_dx = sum(mom$var_dx[A]),
    mu_dy  = sum(mom$mu_dy[A]),
    var_dy = sum(mom$var_dy[A]),
    mu_dxy  = sum(mom$mu_dxy[A]),
    var_dxy = sum(mom$var_dxy[A])
  )
}

# Top-down (backward) greedy elimination based on O_A:
# Start with all features and iteratively remove one feature if it decreases O_A.

top_down_oa <- function(X, Y, q = 20, p = 1, eps = 0) {
  d <- ncol(X)
  A <- seq_len(d)
  
  if (is.null(colnames(X)) || any(colnames(X) == "")) {
    colnames(X) <- paste0("V", seq_len(d))
    colnames(Y) <- colnames(X)
  }
  
  mom <- .feature_moments(X, Y, p = p)
  mu_dx <- mom$mu_dx; var_dx <- mom$var_dx
  mu_dy <- mom$mu_dy; var_dy <- mom$var_dy
  mu_dxy <- mom$mu_dxy; var_dxy <- mom$var_dxy
  
  sum_mu_dx  <- sum(mu_dx);  sum_var_dx  <- sum(var_dx)
  sum_mu_dy  <- sum(mu_dy);  sum_var_dy  <- sum(var_dy)
  sum_mu_dxy <- sum(mu_dxy); sum_var_dxy <- sum(var_dxy)
  
  current <- .overlap_normal_mixture(
    mu_dx = sum_mu_dx,
    var_dx = sum_var_dx,
    mu_dy = sum_mu_dy,
    var_dy = sum_var_dy,
    mu_dxy = sum_mu_dxy,
    var_dxy = sum_var_dxy
  )
  
  path <- current
  
  while (length(A) > q) {
    cand <- A
    
    # Evaluate O_{A\{k}} for all k in the current set A.
    ov <- sapply(cand, function(k) {
      .overlap_normal_mixture(
        mu_dx = sum_mu_dx - mu_dx[k],
        var_dx = sum_var_dx - var_dx[k],
        mu_dy = sum_mu_dy - mu_dy[k],
        var_dy = sum_var_dy - var_dy[k],
        mu_dxy = sum_mu_dxy - mu_dxy[k],
        var_dxy = sum_var_dxy - var_dxy[k]
      )
    })
    
    best <- min(ov)
    kbest <- cand[which.min(ov)]
    
    # Stop if no meaningful improvement is available.
    if ((current - best) < eps || best > current) break
    
    # Remove kbest and update cached sums.
    A <- setdiff(A, kbest)
    sum_mu_dx  <- sum_mu_dx  - mu_dx[kbest];  sum_var_dx  <- sum_var_dx  - var_dx[kbest]
    sum_mu_dy  <- sum_mu_dy  - mu_dy[kbest];  sum_var_dy  <- sum_var_dy  - var_dy[kbest]
    sum_mu_dxy <- sum_mu_dxy - mu_dxy[kbest]; sum_var_dxy <- sum_var_dxy - var_dxy[kbest]
    
    current <- best
    path <- c(path, current)
  }
  
  list(idx = A, feat = colnames(X)[A], oa_path = path)
}

# Stability selection wrapper for the top-down overlap procedure
# Repeat top-down on subsamples and retain features with selection frequency >= pi_thr

stability_top_down <- function(X, Y,
                               B = 100, sub = 0.5, q = 20, p = 1,
                               pi_thr = 0.6, eps = 0, min_n = 3L, seed = NULL,
                               prefix = "", progress_every = 10L) {
  if (!is.null(seed)) set.seed(seed)
  d <- ncol(X)
  
  if (is.null(colnames(X)) || any(colnames(X) == "")) {
    colnames(X) <- paste0("V", seq_len(d))
  }
  if (is.null(colnames(Y)) || any(colnames(Y) == "")) {
    colnames(Y) <- colnames(X)
  }
  
  nx <- nrow(X); ny <- nrow(Y)
  cnt <- integer(d)
  
  for (b in seq_len(B)) {
    nx_sub <- max(min_n, floor(sub * nx))
    ny_sub <- max(min_n, floor(sub * ny))
    ix <- sample.int(nx, nx_sub, replace = FALSE)
    iy <- sample.int(ny, ny_sub, replace = FALSE)
    
    S <- top_down_oa(X[ix, , drop = FALSE], Y[iy, , drop = FALSE], q = q, p = p, eps = eps)$idx
    if (length(S)) cnt[S] <- cnt[S] + 1L
    
    if (progress_every > 0 && (b %% progress_every == 0 || b == B)) {
      cat(sprintf("%sTop-down stability: %d/%d (%.0f%%)\n", prefix, b, B, 100 * b / B))
      flush.console()
    }
  }
  
  pi <- cnt / B
  names(pi) <- colnames(X)
  
  ord <- order(pi, decreasing = TRUE)
  sel <- which(pi >= pi_thr)
  if (length(sel) == 0) sel <- head(ord, q)
  
  list(pi = pi, order = ord, selected_idx = sel, selected_feat = names(pi)[sel])
}




# -------------------- Differential expression via limma --------------------

run_limma_de <- function(mat) {
  stopifnot("label" %in% colnames(mat))
  
  empty_de <- function() {
    data.frame(
      logFC = numeric(0),
      AveExpr = numeric(0),
      t = numeric(0),
      P.Value = numeric(0),
      adj.P.Val = numeric(0),
      B = numeric(0),
      feature = character(0)
    )
  }
  
  grp <- factor(mat$label, levels = c("normal", "tumor"))
  keep <- !is.na(grp)
  X <- mat[keep, , drop = FALSE]
  grp <- droplevels(grp[keep])
  
  # Need both classes
  if (length(levels(grp)) < 2) return(empty_de())
  
  E0 <- t(as.matrix(X[, setdiff(colnames(X), "label"), drop = FALSE]))
  # Force numeric; non-numeric becomes NA
  E <- suppressWarnings(matrix(as.numeric(E0), nrow = nrow(E0), dimnames = dimnames(E0)))
  E[!is.finite(E)] <- NA_real_
  
  # Drop features with <2 finite values
  keep_row <- rowSums(is.finite(E)) >= 2
  E <- E[keep_row, , drop = FALSE]
  if (nrow(E) < 2) return(empty_de())
  
  # Impute remaining NAs per-feature (row) with row median (or 0 if undefined)
  if (anyNA(E)) {
    for (r in seq_len(nrow(E))) {
      if (anyNA(E[r, ])) {
        rr <- E[r, ]
        med <- median(rr[is.finite(rr)], na.rm = TRUE)
        if (!is.finite(med)) med <- 0
        rr[!is.finite(rr)] <- med
        E[r, ] <- rr
      }
    }
  }
  
  # Remove near-zero-variance features (after imputation)
  rv <- apply(E, 1, function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2) return(0)
    var(z)
  })
  E <- E[rv > 1e-10, , drop = FALSE]
  if (nrow(E) < 2) return(empty_de())
  
  design <- model.matrix(~0 + grp)
  colnames(design) <- levels(grp)
  
  fit <- tryCatch(lmFit(E, design), error = function(e) NULL)
  if (is.null(fit)) return(empty_de())
  
  cont <- makeContrasts(tumor - normal, levels = design)
  fitc <- tryCatch(contrasts.fit(fit, cont), error = function(e) NULL)
  if (is.null(fitc)) return(empty_de())
  
  fit2 <- tryCatch(
    eBayes(fitc, trend = TRUE, robust = TRUE),
    error = function(e) {
      tryCatch(eBayes(fitc, trend = FALSE, robust = TRUE), error = function(e2) NULL)
    }
  )
  if (is.null(fit2)) return(empty_de())
  
  de <- tryCatch(topTable(fit2, coef = 1, number = Inf, sort.by = "P"), error = function(e) NULL)
  if (is.null(de)) return(empty_de())
  
  de$feature <- rownames(de)
  de
}


# -------------------- ReliefF (binary) on matrices ---------------------------

reliefF_rank_features <- function(mat, k_panel = 20, m = NULL, k_nn = 10, seed = 123) {
  set.seed(seed)
  
  feat_cols <- setdiff(colnames(mat), "label")
  y <- factor(mat$label, levels = c("normal", "tumor"))
  keep_samp <- !is.na(y)
  
  X <- data.matrix(mat[keep_samp, feat_cols, drop = FALSE])
  y <- droplevels(y[keep_samp])
  
  # Impute NAs
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[, j])) {
      X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)
    }
  }
  
  # Normalize features
  sds <- apply(X, 2, sd, na.rm = TRUE)
  sds[sds == 0 | is.na(sds)] <- 1
  Xn <- sweep(X, 2, colMeans(X, na.rm = TRUE), "-")
  Xn <- sweep(Xn, 2, sds, "/")
  Xn[is.na(Xn)] <- 0
  
  n <- nrow(Xn)
  d <- ncol(Xn)
  
  if (is.null(m)) m <- n
  m <- min(m, n)
  k_eff <- min(k_nn, n - 1)
  
  W <- rep(0, d)
  feat_kept <- colnames(Xn)
  
  for (i in seq_len(m)) {
    yi <- y[i]
    xi <- Xn[i, , drop = FALSE]
    
    diffs <- sweep(Xn, 2, as.numeric(xi), "-")
    dist2 <- rowSums(diffs * diffs)
    dist2[i] <- Inf
    
    hit_idx  <- which(y == yi)
    miss_idx <- which(y != yi)
    
    hit_ord  <- hit_idx[order(dist2[hit_idx])]
    miss_ord <- miss_idx[order(dist2[miss_idx])]
    
    k_hit  <- min(k_eff, length(hit_ord))
    k_miss <- min(k_eff, length(miss_ord))
    
    hit_nn  <- head(hit_ord,  k_hit)
    miss_nn <- head(miss_ord, k_miss)
    
    dh <- if (k_hit > 0) {
      colMeans(abs(sweep(Xn[hit_nn, , drop = FALSE], 2, as.numeric(xi), "-")))
    } else {
      rep(0, d)
    }
    dm <- if (k_miss > 0) {
      colMeans(abs(sweep(Xn[miss_nn, , drop = FALSE], 2, as.numeric(xi), "-")))
    } else {
      rep(0, d)
    }
    
    W <- W - dh / m + dm / m
  }
  
  names(W) <- feat_kept
  ord <- order(W, decreasing = TRUE, na.last = NA)
  ranked <- data.frame(feature = names(W)[ord], weight = as.numeric(W[ord]))
  panel <- head(ranked$feature, k_panel)
  
  list(ranked = ranked, panel = panel)
}

# -------------------- ElasticNet feature ranking (glmnet) ---------------

elasticnet_rank_features <- function(mat, k_panel = 20, alpha = 0.5, seed = 123) {
  set.seed(seed)
  
  feat_cols <- setdiff(colnames(mat), "label")
  yfac <- factor(mat$label, levels = c("normal", "tumor"))
  keep_samp <- !is.na(yfac)
  
  X <- data.matrix(mat[keep_samp, feat_cols, drop = FALSE])
  y <- as.numeric(yfac[keep_samp] == "tumor")
  
  if (anyNA(X)) {
    for (j in seq_len(ncol(X))) {
      if (anyNA(X[, j])) {
        med <- median(X[, j], na.rm = TRUE)
        X[is.na(X[, j]), j] <- med
      }
    }
  }
  
  if (length(unique(y)) < 2) stop("ElasticNet skipped: only one class after labeling.")
  
  nfolds <- min(5, nrow(X))
  if (nfolds < 3) nfolds <- 3
  
  fit <- suppressWarnings(tryCatch({
    cv.glmnet(X, y, family = "binomial", alpha = alpha, nfolds = nfolds)
  }, error = function(e) NULL))
  
  if (is.null(fit)) {
    v <- apply(X, 2, var, na.rm = TRUE)
    v[is.na(v)] <- -Inf
    ord <- order(v, decreasing = TRUE)
    ranked <- data.frame(
      feature = colnames(X)[ord],
      coef_abs = as.numeric(v[ord]),
      coef = NA_real_
    )
    panel <- head(ranked$feature, k_panel)
    return(list(ranked = ranked, panel = panel, fit = NULL))
  }
  
  cf <- as.matrix(coef(fit, s = "lambda.min"))
  rn <- rownames(cf)
  if (!is.null(rn) && "(Intercept)" %in% rn) {
    cf <- cf[rn != "(Intercept)", , drop = FALSE]
  }
  beta <- as.numeric(cf[, 1])
  names(beta) <- rownames(cf)
  
  ord <- order(abs(beta), decreasing = TRUE, na.last = NA)
  ranked <- data.frame(
    feature = names(beta)[ord],
    coef_abs = abs(beta[ord]),
    coef = beta[ord]
  )
  
  nz <- ranked$coef_abs > 0
  panel <- head(ranked$feature[nz], k_panel)
  if (length(panel) < k_panel) {
    panel <- unique(c(panel, head(ranked$feature, k_panel)))
    panel <- head(panel, k_panel)
  }
  
  list(ranked = ranked, panel = panel, fit = fit)
}


# -------------------- LASSO feature ranking (glmnet; classical L1)
# Implemented as ElasticNet with alpha = 1

lasso_rank_features <- function(mat, k_panel = 20, seed = 123) {
  elasticnet_rank_features(mat, k_panel = k_panel, alpha = 1, seed = seed)
}




# ------------ mRMR (mutual information filter; greedy difference criterion) -
# Discretization-based MI via infotheo; suitable as a model-agnostic filter in omics

mrmr_rank_features <- function(mat, k_panel = 20, nbins = 5, seed = 123) {
  set.seed(seed)
  
  .empty <- function() {
    list(
      ranked = data.frame(feature = character(0), step = integer(0), score = numeric(0), relevance_mi = numeric(0)),
      panel  = character(0)
    )
  }
  
  feat_cols <- setdiff(colnames(mat), "label")
  yfac <- factor(mat$label, levels = c("normal", "tumor"))
  keep <- !is.na(yfac)
  
  if (sum(keep) < 6) return(.empty())
  if (length(feat_cols) < 1) return(.empty())
  
  X <- data.matrix(mat[keep, feat_cols, drop = FALSE])
  storage.mode(X) <- "double"
  X[is.infinite(X)] <- NA_real_
  
  y <- as.integer(yfac[keep] == "tumor")
  if (sum(y == 0) < 2 || sum(y == 1) < 2) return(.empty())
  
  # Drop features with too few finite values
  
  finite_n <- colSums(is.finite(X))
  keep_feat <- finite_n >= 3
  if (!any(keep_feat)) return(.empty())
  X <- X[, keep_feat, drop = FALSE]
  feat_cols <- feat_cols[keep_feat]
  
  # Median impute per feature, drop features with undefined median
  
  med <- apply(X, 2, function(v) median(v, na.rm = TRUE))
  good_med <- is.finite(med)
  if (!any(good_med)) return(.empty())
  X <- X[, good_med, drop = FALSE]
  feat_cols <- feat_cols[good_med]
  med <- med[good_med]
  
  for (j in seq_len(ncol(X))) {
    bad <- !is.finite(X[, j])
    if (any(bad)) X[bad, j] <- med[j]
  }
  
  # Drop constant / near-constant features (mutual information becomes degenerate)
  
  uniq_n <- apply(X, 2, function(v) length(unique(v)))
  sds    <- apply(X, 2, stats::sd)
  keep2  <- uniq_n >= 2 & is.finite(sds) & sds > 0
  if (!any(keep2)) return(.empty())
  X <- X[, keep2, drop = FALSE]
  feat_cols <- feat_cols[keep2]
  
  # Safe discretization per feature (equal-frequency; falls back to pretty() breaks)
  
  safe_disc <- function(v, nb) {
    u <- length(unique(v))
    nbj <- min(nb, u)
    if (nbj < 2) return(factor(rep(1L, length(v))))
    
    qs <- unique(stats::quantile(v, probs = seq(0, 1, length.out = nbj + 1), na.rm = TRUE, type = 7))
    if (length(qs) <= 2) return(factor(rep(1L, length(v))))
    
    b <- tryCatch(cut(v, breaks = qs, include.lowest = TRUE, labels = FALSE),
                  error = function(e) NULL)
    if (is.null(b) || anyNA(b)) {
      br <- pretty(v, n = nbj)
      b <- cut(v, breaks = br, include.lowest = TRUE, labels = FALSE)
    }
    factor(b)
  }
  
  Xd <- as.data.frame(lapply(seq_len(ncol(X)), function(j) safe_disc(X[, j], nbins)))
  names(Xd) <- feat_cols
  yd <- factor(y)
  
  safe_mi <- function(a, b) {
    out <- tryCatch(infotheo::mutinformation(a, b), error = function(e) NA_real_)
    if (!is.finite(out)) 0 else out
  }
  
  # Relevance I(X_j;Y)
  
  mi_y <- sapply(seq_len(ncol(Xd)), function(j) safe_mi(Xd[[j]], yd))
  names(mi_y) <- feat_cols
  if (all(mi_y <= 0)) {
    # Still return something deterministic, top k_panel by (zero) relevance
    
    k_eff <- min(k_panel, length(feat_cols))
    sel <- head(feat_cols, k_eff)
    ranked <- data.frame(
      feature = sel,
      step = seq_along(sel),
      score = rep(0, length(sel)),
      relevance_mi = rep(0, length(sel))
    )
    return(list(ranked = ranked, panel = sel))
  }
  
  p <- length(feat_cols)
  k_eff <- min(k_panel, p)
  
  # first feature: max relevance
  
  S <- integer(0)
  score_path <- numeric(0)
  
  j1 <- which.max(mi_y)
  S <- c(S, j1)
  score_path <- c(score_path, mi_y[j1])
  
  if (k_eff >= 2) {
    for (t in 2:k_eff) {
      cand <- setdiff(seq_len(p), S)
      if (length(cand) == 0) break
      
      red <- sapply(cand, function(j) {
        mean(sapply(S, function(k) safe_mi(Xd[[j]], Xd[[k]])), na.rm = TRUE)
      })
      rel <- mi_y[cand]
      sc <- rel - red
      
      jbest <- cand[which.max(sc)]
      S <- c(S, jbest)
      score_path <- c(score_path, max(sc, na.rm = TRUE))
    }
  }
  
  selected <- feat_cols[S]
  ranked <- data.frame(
    feature = selected,
    step = seq_along(selected),
    score = score_path,
    relevance_mi = as.numeric(mi_y[S])
  )
  
  list(ranked = ranked, panel = selected)
}



# -------------------- AUC comparison helpers ---------

compute_auc_quiet <- function(y_true, y_score) {
  y_true <- as.numeric(y_true)
  y_score <- as.numeric(y_score)
  
  # Drop NA pairs 
  
  keep <- !(is.na(y_true) | is.na(y_score))
  y_true <- y_true[keep]
  y_score <- y_score[keep]
  
  # If only one class is present, AUC is not identifiable; return 0.5 (no-skill)
  
  if (length(unique(y_true)) < 2) return(0.5)
  
  # Degenerate predictor (constant score) => no-skill
  
  if (length(unique(y_score)) < 2) return(0.5)
  
  roc_obj <- suppressMessages(suppressWarnings(
    pROC::roc(response = y_true, predictor = y_score,
              levels = c(0, 1), direction = "<", quiet = TRUE)
  ))
  
  out <- suppressMessages(suppressWarnings(as.numeric(pROC::auc(roc_obj))))
  if (is.na(out) || !is.finite(out)) 0.5 else out
}


# -------------------- AUC comparison (glmnet) ----------------

compare_panels_auc_glmnet <- function(mat, panels, n_splits = 100, train_frac = 0.7, seed = 123) {
  set.seed(seed)
  
  yfac <- factor(mat$label, levels = c("normal", "tumor"))
  ybin <- as.numeric(yfac == "tumor")
  
  feat_cols <- setdiff(colnames(mat), "label")
  results <- list()
  
  for (panel_name in names(panels)) {
    feats <- intersect(panels[[panel_name]], feat_cols)
    
    # only skip if 0 features
    
    if (length(feats) < 1) {
      results[[panel_name]] <- data.frame(
        Panel = panel_name,
        AUC_mean = NA_real_,
        AUC_sd = NA_real_,
        AUC_se = NA_real_,
        N_features = length(feats),
        N_splits = 0L
      )
      next
    }
    
    X <- data.matrix(mat[, feats, drop = FALSE])  # drop=FALSE keeps 1-feature as matrix
    y <- ybin
    keep <- !is.na(y)
    X <- X[keep, , drop = FALSE]
    y <- y[keep]
    
    aucs <- numeric(0)
    
    for (i in seq_len(n_splits)) {
      n <- length(y)
      
      # Stratified split
      idx0 <- which(y == 0)
      idx1 <- which(y == 1)
      
      tr0 <- sample(idx0, size = floor(train_frac * length(idx0)), replace = FALSE)
      tr1 <- sample(idx1, size = floor(train_frac * length(idx1)), replace = FALSE)
      
      trn <- sort(c(tr0, tr1))
      tst <- setdiff(seq_len(n), trn)
      
      if (length(tst) == 0) next
      if (length(unique(y[trn])) < 2) next
      if (length(unique(y[tst])) < 2) next
      
      # Train-only scaling
      Xtr <- X[trn, , drop = FALSE]
      Xte <- X[tst, , drop = FALSE]
      
      mu <- colMeans(Xtr)
      sdv <- apply(Xtr, 2, sd)
      sdv[is.na(sdv) | sdv == 0] <- 1
      
      Xtr <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/")
      Xte <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")
      
      Xtr[is.na(Xtr)] <- 0
      Xte[is.na(Xte)] <- 0
      
      fit <- suppressWarnings(tryCatch({
        nfolds <- min(5, length(trn))
        nfolds <- max(3, nfolds)
        cv.glmnet(Xtr, y[trn],
                  family = "binomial", alpha = 0.5, nfolds = nfolds,
                  standardize = FALSE)
      }, error = function(e) NULL))
      
      if (is.null(fit)) next
      
      pred <- tryCatch({
        as.numeric(predict(fit, newx = Xte, s = "lambda.min"))
      }, error = function(e) NULL)
      
      if (is.null(pred)) next
      
      auc_val <- compute_auc_quiet(y[tst], pred)
      if (!is.na(auc_val)) aucs <- c(aucs, auc_val)
    }
    
    if (length(aucs) == 0) {
      Xall <- X
      yall <- y
      
      mu <- colMeans(Xall)
      sdv <- apply(Xall, 2, sd)
      sdv[is.na(sdv) | sdv == 0] <- 1
      
      Xall <- sweep(sweep(Xall, 2, mu, "-"), 2, sdv, "/")
      Xall[is.na(Xall)] <- 0
      
      pred_fb <- NULL
      
      fit_cv <- suppressWarnings(tryCatch({
        nfolds <- min(5L, length(yall))
        nfolds <- max(3L, nfolds)
        cv.glmnet(Xall, yall, family = "binomial", alpha = 0.5, nfolds = nfolds, standardize = FALSE)
      }, error = function(e) NULL))
      
      if (!is.null(fit_cv)) {
        pred_fb <- tryCatch(as.numeric(predict(fit_cv, newx = Xall, s = "lambda.min")), error = function(e) NULL)
      }
      
      if (is.null(pred_fb)) {
        fit_ridge <- suppressWarnings(tryCatch({
          glmnet(Xall, yall, family = "binomial", alpha = 0, standardize = FALSE)
        }, error = function(e) NULL))
        if (!is.null(fit_ridge)) {
          lambda_use <- tail(fit_ridge$lambda, 1)
          pred_fb <- tryCatch(as.numeric(predict(fit_ridge, newx = Xall, s = lambda_use)), error = function(e) NULL)
        }
      }
      
      if (is.null(pred_fb)) pred_fb <- rep(0, length(yall))
      
      aucs <- c(compute_auc_quiet(yall, pred_fb))
    }
    
    auc_mean <- mean(aucs)
    if (length(aucs) == 1) {
      auc_sd <- 0
      auc_se <- 0
    } else {
      auc_sd <- sd(aucs)
      auc_se <- auc_sd / sqrt(length(aucs))
    }
    
    results[[panel_name]] <- data.frame(
      Panel = panel_name,
      AUC_mean = auc_mean,
      AUC_sd = auc_sd,
      AUC_se = auc_se,
      N_features = length(feats),
      N_splits = length(aucs)
    )
  }
  
  do.call(rbind, results)
}

# --------------- Leakage-free AUC comparison --------------------

compare_methods_auc_glmnet_noleak <- function(
    mat,
    group_ids = NULL,
    k_panel,
    n_splits = 1,
    train_frac = 0.7,
    seed = 123,
    dist_prefilter_topN = 100,
    dist_B = 200,
    dist_sub = 0.5,
    dist_k = 20,
    dist_p = 2,
    dist_pi_thr = 0.6,
    dist_eps = 0,
    dist_min_n = 3L,
    td_B = dist_B,
    td_sub = dist_sub,
    td_q = dist_k,
    td_p = 1,
    td_pi_thr = dist_pi_thr,
    td_eps = 0,
    td_min_n = dist_min_n
) {
  # Keep only samples with non-missing labels
  yfac_all <- factor(mat$label, levels = c("normal", "tumor"))
  ybin_all <- as.numeric(yfac_all == "tumor")
  keep_all <- !is.na(ybin_all)
  
  mat2 <- mat[keep_all, , drop = FALSE]
  ybin <- ybin_all[keep_all]
  
  if (is.null(group_ids)) {
    group_ids <- rownames(mat)
  }
  if (length(group_ids) != nrow(mat)) {
    stop("group_ids must have one entry per row of mat.")
  }
  group_ids <- as.character(group_ids)[keep_all]
  
  feat_cols_all <- setdiff(colnames(mat2), "label")
  
  # Pre-generate patient/group-stratified splits once
  
  splits <- make_grouped_stratified_splits(
    ybin = ybin,
    group_ids = group_ids,
    n_splits = n_splits,
    train_frac = train_frac,
    seed = seed
  )
  
  aucs <- list(DE = numeric(0), ReliefF = numeric(0), mRMR = numeric(0), BottomUp = numeric(0), TopDown = numeric(0), ElasticNet = numeric(0), LASSO = numeric(0))
  nfeats <- list(DE = integer(0), ReliefF = integer(0), mRMR = integer(0), BottomUp = integer(0), TopDown = integer(0), ElasticNet = integer(0), LASSO = integer(0))
  
  make_stratified_foldid <- function(y, nfolds, seed_local = 1L) {
    set.seed(seed_local)
    foldid <- integer(length(y))
    for (cls in c(0, 1)) {
      idx <- which(y == cls)
      foldid[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
    }
    foldid
  }
  
  eval_panel_auc <- function(mat_tr, mat_te, panel_feats) {
    feats <- intersect(panel_feats, feat_cols_all)
    
    ytr <- as.numeric(factor(mat_tr$label, levels = c("normal", "tumor")) == "tumor")
    yte <- as.numeric(factor(mat_te$label, levels = c("normal", "tumor")) == "tumor")
    
    # Intercept-only fallback: constant score -> AUC=0.5 
    
    if (length(feats) < 1) {
      return(compute_auc_quiet(yte, rep(0, length(yte))))
    }
    
    Xtr <- data.matrix(mat_tr[, feats, drop = FALSE])
    Xte <- data.matrix(mat_te[, feats, drop = FALSE])
    
    keeptr <- !is.na(ytr)
    keepte <- !is.na(yte)
    Xtr <- Xtr[keeptr, , drop = FALSE]
    Xte <- Xte[keepte, , drop = FALSE]
    ytr <- ytr[keeptr]
    yte <- yte[keepte]
    
    if (nrow(Xtr) < 2 || nrow(Xte) < 2) return(0.5)
    if (length(unique(ytr)) < 2) return(0.5)
    if (length(unique(yte)) < 2) return(0.5)
    
    # Train-only scaling 
    
    mu <- colMeans(Xtr)
    sdv <- apply(Xtr, 2, sd)
    sdv[is.na(sdv) | sdv == 0] <- 1
    
    Xtr <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/")
    Xte <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")
    
    Xtr[is.na(Xtr)] <- 0
    Xte[is.na(Xte)] <- 0
    
    # Plain (unpenalized) logistic regression on all selected features 
    
    fit <- suppressWarnings(tryCatch({
      stats::glm.fit(x = cbind(1, Xtr), y = ytr, family = stats::binomial())
    }, error = function(e) NULL))
    
    if (is.null(fit) || is.null(fit$coefficients)) {
      return(compute_auc_quiet(yte, rep(0, length(yte))))
    }
    
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0
    
    eta <- as.numeric(cbind(1, Xte) %*% beta)
    pred <- 1 / (1 + exp(-eta))
    
    compute_auc_quiet(yte, pred)
  }
  
  
  for (i in seq_len(n_splits)) {
    trn <- splits[[i]]$trn
    tst <- splits[[i]]$tst
    
    if (length(tst) == 0) next
    if (length(unique(ybin[trn])) < 2) next
    if (length(unique(ybin[tst])) < 2) next
    
    mat_tr <- mat2[trn, , drop = FALSE]
    mat_te <- mat2[tst, , drop = FALSE]
    
    # -------------------- Feature selection on train only --------------------
    
    de_tr <- tryCatch(run_limma_de(mat_tr), error = function(e) data.frame(feature = character(0)))
    
    # ---- BottomUp panel --------
    
    panel_bu_tr <- character(0)
    tryCatch({
      top_de_features_tr <- head(de_tr$feature, dist_prefilter_topN)
      available_top_tr <- intersect(top_de_features_tr, feat_cols_all)
      
      if (length(available_top_tr) >= 2) {
        grp <- factor(mat_tr$label, levels = c("normal", "tumor"))
        keepg <- !is.na(grp)
        grp <- droplevels(grp[keepg])
        
        E_sub <- data.matrix(mat_tr[keepg, available_top_tr, drop = FALSE])
        colnames(E_sub) <- available_top_tr
        
        Xs <- E_sub[grp == "normal", , drop = FALSE]
        Ys <- E_sub[grp == "tumor",  , drop = FALSE]
        
        if (nrow(Xs) >= dist_min_n && nrow(Ys) >= dist_min_n) {
          out_dist <- stability_bottom_up(
            Xs, Ys,
            B = dist_B, sub = dist_sub, k = dist_k, p = dist_p,
            pi_thr = dist_pi_thr, eps = dist_eps, min_n = dist_min_n,
            seed = seed + i,
            prefix = paste0("[AUC:", i, "] "),
            progress_every = max(1L, floor(dist_B / 20L))
          )
          
          panel_bu_tr <- out_dist$selected_feat
          if (is.null(panel_bu_tr)) panel_bu_tr <- character(0)
        }
      }
    }, error = function(e) {
      panel_bu_tr <<- character(0)
    })
    
    # Compute BottomUp AUC 
    
    # Fair sizing: all other (non-TopDown) methods use the realized BottomUp panel size
    
    k_eff <- length(intersect(panel_bu_tr, feat_cols_all))
    
    if (k_eff < 1) {
      panel_de_tr <- character(0)
      panel_relief_tr <- character(0)
      panel_en_tr <- character(0)
      panel_lasso_tr <- character(0)
      panel_mrmr_tr <- character(0)
    } else {
      panel_de_tr <- head(de_tr$feature, k_eff)
      
      rel_tr <- tryCatch(reliefF_rank_features(mat_tr, k_panel = k_eff, m = NULL, k_nn = 10, seed = seed + i),
                         error = function(e) list(panel = character(0)))
      panel_relief_tr <- rel_tr$panel
      if (is.null(panel_relief_tr)) panel_relief_tr <- character(0)
      
      en_tr <- tryCatch(elasticnet_rank_features(mat_tr, k_panel = k_eff, alpha = 0.5, seed = seed + i),
                        error = function(e) list(panel = character(0)))
      panel_en_tr <- en_tr$panel
      if (is.null(panel_en_tr)) panel_en_tr <- character(0)
      
      la_tr <- tryCatch(lasso_rank_features(mat_tr, k_panel = k_eff, seed = seed + i),
                        error = function(e) list(panel = character(0)))
      panel_lasso_tr <- la_tr$panel
      if (is.null(panel_lasso_tr)) panel_lasso_tr <- character(0)
      
      mrmr_tr <- tryCatch(mrmr_rank_features(mat_tr, k_panel = k_eff, nbins = 5, seed = seed + i),
                          error = function(e) list(panel = character(0)))
      panel_mrmr_tr <- mrmr_tr$panel
      if (is.null(panel_mrmr_tr)) panel_mrmr_tr <- character(0)
    }
    
    
    panel_td_tr <- character(0)
    tryCatch({
      top_de_features_tr <- head(de_tr$feature, dist_prefilter_topN)
      available_top_tr <- intersect(top_de_features_tr, feat_cols_all)
      
      if (length(available_top_tr) >= 2) {
        grp <- factor(mat_tr$label, levels = c("normal", "tumor"))
        keepg <- !is.na(grp)
        grp <- droplevels(grp[keepg])
        
        E_sub <- data.matrix(mat_tr[keepg, available_top_tr, drop = FALSE])
        colnames(E_sub) <- available_top_tr
        
        Xs <- E_sub[grp == "normal", , drop = FALSE]
        Ys <- E_sub[grp == "tumor",  , drop = FALSE]
        
        if (nrow(Xs) >= td_min_n && nrow(Ys) >= td_min_n) {
          out_td <- stability_top_down(
            Xs, Ys,
            B = td_B, sub = td_sub, q = td_q, p = td_p,
            pi_thr = td_pi_thr, eps = td_eps, min_n = td_min_n,
            seed = seed + i,
            prefix = paste0("[AUC:", i, "] "),
            progress_every = max(1L, floor(td_B / 20L))
          )
          
          panel_td_tr <- out_td$selected_feat
          if (is.null(panel_td_tr)) panel_td_tr <- character(0)
        }
      }
    }, error = function(e) {
      panel_td_tr <<- character(0)
    })
    
    # -------------------- Evaluate on test data --------------------
    
    a_de <- eval_panel_auc(mat_tr, mat_te, panel_de_tr)
    a_rf <- eval_panel_auc(mat_tr, mat_te, panel_relief_tr)
    a_en <- eval_panel_auc(mat_tr, mat_te, panel_en_tr)
    a_la <- eval_panel_auc(mat_tr, mat_te, panel_lasso_tr)
    a_mr <- eval_panel_auc(mat_tr, mat_te, panel_mrmr_tr)
    a_bu <- eval_panel_auc(mat_tr, mat_te, panel_bu_tr)
    a_td <- eval_panel_auc(mat_tr, mat_te, panel_td_tr)
    
    aucs$DE        <- c(aucs$DE, a_de)
    aucs$ReliefF   <- c(aucs$ReliefF, a_rf)
    aucs$ElasticNet<- c(aucs$ElasticNet, a_en)
    aucs$LASSO     <- c(aucs$LASSO, a_la)
    aucs$mRMR      <- c(aucs$mRMR, a_mr)
    aucs$BottomUp  <- c(aucs$BottomUp, a_bu)
    aucs$TopDown   <- c(aucs$TopDown, a_td)
    
    nfeats$DE         <- c(nfeats$DE, length(intersect(panel_de_tr, feat_cols_all)))
    nfeats$ReliefF    <- c(nfeats$ReliefF, length(intersect(panel_relief_tr, feat_cols_all)))
    nfeats$ElasticNet <- c(nfeats$ElasticNet, length(intersect(panel_en_tr, feat_cols_all)))
    nfeats$LASSO      <- c(nfeats$LASSO, length(intersect(panel_lasso_tr, feat_cols_all)))
    nfeats$mRMR       <- c(nfeats$mRMR, length(intersect(panel_mrmr_tr, feat_cols_all)))
    nfeats$BottomUp   <- c(nfeats$BottomUp, length(intersect(panel_bu_tr, feat_cols_all)))
    nfeats$TopDown    <- c(nfeats$TopDown, length(intersect(panel_td_tr, feat_cols_all)))
  }
  
  summarize_one <- function(method_name) {
    v <- aucs[[method_name]]
    v <- v[!is.na(v)]
    
    auc_mean <- if (length(v)) mean(v) else NA_real_
    if (length(v) == 1) {
      auc_sd <- 0
      auc_se <- 0
    } else if (length(v) >= 2) {
      auc_sd <- sd(v)
      auc_se <- auc_sd / sqrt(length(v))
    } else {
      auc_sd <- NA_real_
      auc_se <- NA_real_
    }
    
    nf <- nfeats[[method_name]]
    nf <- nf[!is.na(nf)]
    nfeat_rep <- if (length(nf)) as.integer(round(median(nf))) else NA_integer_
    
    data.frame(
      Panel = method_name,
      AUC_mean = auc_mean,
      AUC_sd = auc_sd,
      AUC_se = auc_se,
      N_features = nfeat_rep,
      N_splits = length(v)
    )
  }
  
  do.call(rbind, lapply(c("DE", "ReliefF", "mRMR", "BottomUp", "TopDown", "ElasticNet", "LASSO"), summarize_one))
}


# -------------------- Volcano plot with panel colorings -----------------

make_volcano_plot <- function(de, panels, out_pdf) {
  df <- de
  
  # Ensure required columns exist even if DE is empty
  if (is.null(df) || nrow(df) == 0) {
    df <- data.frame(
      logFC = numeric(0),
      P.Value = numeric(0),
      feature = character(0)
    )
  }
  
  df$P.Value <- pmax(df$P.Value, .Machine$double.xmin)
  df$negLogP <- -log10(df$P.Value)
  df$logFC_plot <- pmax(pmin(df$logFC, 10), -10)
  
  # Membership across all panels
  panel_names <- names(panels)
  if (is.null(panel_names) || length(panel_names) == 0) panel_names <- character(0)
  
  if (length(panel_names) > 0) {
    for (nm in panel_names) {
      feats <- panels[[nm]]
      if (is.null(feats)) feats <- character(0)
      df[[paste0("in_", nm)]] <- df$feature %in% feats
    }
    in_cols <- paste0("in_", panel_names)
    n_in <- rowSums(as.matrix(df[, in_cols, drop = FALSE]))
  } else {
    n_in <- rep(0L, nrow(df))
  }
  
  df$group <- "None"
  if (length(panel_names) > 0) {
    for (nm in panel_names) {
      coln <- paste0("in_", nm)
      df$group[df[[coln]] & (n_in == 1)] <- nm
    }
  }
  df$group[n_in > 1] <- "Multiple"
  df$label_me <- (n_in > 0)
  
  cols <- c(
    "None"       = "grey70",
    "DE"         = "dodgerblue3",
    "ReliefF"    = "green3",
    "mRMR"       = "orange2",
    "BottomUp"   = "red3",
    "TopDown"    = "purple3",
    "ElasticNet" = "brown3",
    "LASSO"      = "cyan4",
    "Multiple"   = "black"
  )
  
  # Only keep colors for groups that exist in this dataset 
  
  cols_use <- cols[names(cols) %in% unique(df$group)]
  
  p <- ggplot(df, aes(x = logFC_plot, y = negLogP)) +
    geom_point(aes(color = group), alpha = 0.85) +
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_text_repel(
      data = df[df$label_me, , drop = FALSE],
      aes(label = feature),
      size = 3, max.overlaps = Inf
    ) +
    scale_color_manual(values = cols_use, breaks = names(cols_use)) +
    coord_cartesian(xlim = c(-10, 10)) +
    theme_minimal() +
    labs(x = "logFC (tumor - normal)", y = "-log10(P)",
         title = "Volcano plot with all panels (single vs multiple overlap)")
  
  ggsave(out_pdf, p, width = 8, height = 5.5)
  invisible(p)
}

# -------------------- One-dataset processing wrapper -----------

safe_process_one <- function(gse_id) {
  on.exit({ gc() }, add = TRUE)
  log_file <- file.path(results_dir, paste0(gse_id, "_log.txt"))
  msg_prefix <- paste0("[", gse_id, "] ")
  
  # Global progress bookkeeping (rolling tail log)
  tryCatch(.update_id_list(started_file, gse_id), error = function(e) {})
  tryCatch(global_log_event(gse_id, "START dataset"), error = function(e) {})
  
  tryCatch({
    withTimeout({
      
      # Recreate tempdir() if it was deleted mid-run (Windows cleanup / AV)
      
      .ensure_tempdir()
      
      cat(Sys.time(), msg_prefix, "Processing dataset...\n")      
      cached_file <- file.path(download_dir, paste0(gse_id, "_miRNA_by_samples.csv"))
      if (!file.exists(cached_file)) stop(paste0("Cached dataset not found: ", cached_file))
      
      mat <- read.csv(cached_file, row.names = 1, check.names = FALSE)
      if (!("label" %in% colnames(mat))) stop("Cached dataset has no 'label' column.")
      group_ids <- get_subject_ids_for_gse(gse_id, rownames(mat))
      
      write.csv(mat, file.path(results_dir, paste0(gse_id, "_miRNA_by_samples.csv")),
                row.names = TRUE)
      
      global_log_event(gse_id, "DE: start")
      # ---- DE ----
      de <- run_limma_de(mat)
      global_log_event(gse_id, "DE: done")
      write.csv(de, file.path(results_dir, paste0(gse_id, "_DE.csv")), row.names = FALSE)
      panel_de <- head(de$feature, k_panel)
      
      global_log_event(gse_id, "ReliefF: start")
      # ---- ReliefF ----
      rel <- reliefF_rank_features(mat, k_panel = k_panel, m = NULL, k_nn = 10, seed = 123)
      panel_relief <- rel$panel
      write.csv(rel$ranked, file.path(results_dir, paste0(gse_id, "_ReliefF_weights.csv")), row.names = FALSE)
      write.csv(data.frame(feature = panel_relief),
                file.path(results_dir, paste0(gse_id, "_ReliefF_panel.csv")),
                row.names = FALSE)
      global_log_event(gse_id, "ReliefF: done")
      
      # ---- ElasticNet ----
      
      en <- elasticnet_rank_features(mat, k_panel = k_panel, alpha = 0.5, seed = 123)
      panel_en <- en$panel
      write.csv(en$ranked, file.path(results_dir, paste0(gse_id, "_ElasticNet_coefficients.csv")), row.names = FALSE)
      write.csv(data.frame(feature = panel_en),
                file.path(results_dir, paste0(gse_id, "_ElasticNet_panel.csv")),
                row.names = FALSE)
      global_log_event(gse_id, "ElasticNet: done")
      
      
      global_log_event(gse_id, "LASSO: start")
      
      # ---- LASSO (classical) ----
      
      la <- lasso_rank_features(mat, k_panel = k_panel, seed = 123)
      panel_lasso <- la$panel
      write.csv(la$ranked, file.path(results_dir, paste0(gse_id, "_LASSO_coefficients.csv")), row.names = FALSE)
      write.csv(data.frame(feature = panel_lasso),
                file.path(results_dir, paste0(gse_id, "_LASSO_panel.csv")),
                row.names = FALSE)
      global_log_event(gse_id, "LASSO: done")
      
      # ---- mRMR (mutual information filter) ----
      
      panel_mrmr <- character(0)
      mr_ranked <- data.frame(feature = character(0), step = integer(0), score = numeric(0), relevance_mi = numeric(0))
      
      global_log_event(gse_id, "mRMR: start")
      tryCatch({
        mr <- mrmr_rank_features(mat, k_panel = k_panel, nbins = 5, seed = 123)
        panel_mrmr <- mr$panel
        mr_ranked <- mr$ranked
        
        write.csv(mr_ranked, file.path(results_dir, paste0(gse_id, "_mRMR_ranked.csv")), row.names = FALSE)
        write.csv(data.frame(feature = panel_mrmr),
                  file.path(results_dir, paste0(gse_id, "_mRMR_panel.csv")),
                  row.names = FALSE)
        
        global_log_event(gse_id, "mRMR: done")
      }, error = function(e) {
        global_log_event(gse_id, paste0("mRMR: ERROR - ", conditionMessage(e)))
        
        # Write empty outputs so downstream plotting/evaluation doesn't fail
        write.csv(mr_ranked, file.path(results_dir, paste0(gse_id, "_mRMR_ranked.csv")), row.names = FALSE)
        write.csv(data.frame(feature = panel_mrmr),
                  file.path(results_dir, paste0(gse_id, "_mRMR_panel.csv")),
                  row.names = FALSE)
      })
      
      
      
      global_log_event(gse_id, "Distance selection: start")
      
      # ---- Distance-based ----
      # Bottom-up (P_A) and Top-down (overlap) use the SAME DE prefiltering (top dist_prefilter_topN DE features)
      
      panel_bu <- character(0)
      panel_td <- character(0)
      
      tryCatch({
        feat_cols <- setdiff(colnames(mat), "label")
        top_de_features <- head(de$feature, dist_prefilter_topN)
        available_top <- intersect(top_de_features, feat_cols)
        
        if (length(available_top) < 2) {
          cat(msg_prefix, "Distance selection skipped: <2 available prefiltered features.\n")
          panel_bu <- character(0)
          panel_td <- character(0)
        } else {
          grp <- factor(mat$label, levels = c("normal", "tumor"))
          keep <- !is.na(grp)
          grp <- droplevels(grp[keep])
          
          E_sub <- data.matrix(mat[keep, available_top, drop = FALSE])
          colnames(E_sub) <- available_top
          
          Xs <- E_sub[grp == "normal", , drop = FALSE]
          Ys <- E_sub[grp == "tumor",  , drop = FALSE]
          
          if (nrow(Xs) < dist_min_n || nrow(Ys) < dist_min_n) {
            cat(msg_prefix, "Distance selection skipped: too few samples after labeling.\n")
            panel_bu <- character(0)
            panel_td <- character(0)
          } else {
            cat(msg_prefix, sprintf("Distance selection on %d prefiltered features...\n", ncol(E_sub)))
            
            # Bottom-up 
            
            tryCatch({
              out_bu <- stability_bottom_up(
                Xs, Ys,
                B = dist_B, sub = dist_sub, k = dist_k, p = dist_p,
                pi_thr = dist_pi_thr, eps = dist_eps, min_n = dist_min_n,
                seed = 123,
                prefix = msg_prefix,
                progress_every = max(1L, floor(dist_B / 20L))
              )
              
              panel_bu <- out_bu$selected_feat
              if (is.null(panel_bu)) panel_bu <- character(0)
              
              write.csv(
                data.frame(feature = panel_bu),
                file.path(results_dir, paste0(gse_id, "_bottom_up_panel.csv")),
                row.names = FALSE
              )
              
              pi_vec <- out_bu$pi
              if (is.null(pi_vec)) pi_vec <- numeric(0)
              
              feat_names <- names(pi_vec)
              if (is.null(feat_names) || length(feat_names) == 0) feat_names <- colnames(Xs)
              if (is.null(feat_names) || length(feat_names) == 0) feat_names <- paste0("feat_", seq_along(pi_vec))
              if (length(feat_names) != length(pi_vec)) {
                m <- min(length(feat_names), length(pi_vec))
                feat_names <- feat_names[seq_len(m)]
                pi_vec <- pi_vec[seq_len(m)]
              }
              
              pi_df <- data.frame(feature = feat_names, pi = as.numeric(pi_vec))
              write.csv(
                pi_df,
                file.path(results_dir, paste0(gse_id, "_bottom_up_pi.csv")),
                row.names = FALSE
              )
            }, error = function(e) {
              cat(Sys.time(), msg_prefix, "Bottom-up selection ERROR (continuing):", e$message, "\n")
              write(paste0(Sys.time(), " ", msg_prefix, "Bottom-up selection ERROR: ", e$message, "\n"),
                    file = log_file, append = TRUE)
              panel_bu <<- character(0)
            })
            
            # Top-down
            
            tryCatch({
              out_td <- stability_top_down(
                Xs, Ys,
                B = td_B, sub = td_sub, q = td_q, p = td_p,
                pi_thr = td_pi_thr, eps = td_eps, min_n = td_min_n,
                seed = 123,
                prefix = msg_prefix,
                progress_every = max(1L, floor(td_B / 20L))
              )
              
              panel_td <- out_td$selected_feat
              if (is.null(panel_td)) panel_td <- character(0)
              
              write.csv(
                data.frame(feature = panel_td),
                file.path(results_dir, paste0(gse_id, "_top_down_panel.csv")),
                row.names = FALSE
              )
              
              pi_vec_td <- out_td$pi
              if (is.null(pi_vec_td)) pi_vec_td <- numeric(0)
              
              feat_names_td <- names(pi_vec_td)
              if (is.null(feat_names_td) || length(feat_names_td) == 0) feat_names_td <- colnames(Xs)
              if (is.null(feat_names_td) || length(feat_names_td) == 0) feat_names_td <- paste0("feat_", seq_along(pi_vec_td))
              if (length(feat_names_td) != length(pi_vec_td)) {
                m <- min(length(feat_names_td), length(pi_vec_td))
                feat_names_td <- feat_names_td[seq_len(m)]
                pi_vec_td <- pi_vec_td[seq_len(m)]
              }
              
              pi_df_td <- data.frame(feature = feat_names_td, pi = as.numeric(pi_vec_td))
              write.csv(
                pi_df_td,
                file.path(results_dir, paste0(gse_id, "_top_down_pi.csv")),
                row.names = FALSE
              )
            }, error = function(e) {
              cat(Sys.time(), msg_prefix, "Top-down selection ERROR (continuing):", e$message, "\n")
              write(paste0(Sys.time(), " ", msg_prefix, "Top-down selection ERROR: ", e$message, "\n"),
                    file = log_file, append = TRUE)
              panel_td <<- character(0)
            })
          }
        }
      }, error = function(e) {
        cat(Sys.time(), msg_prefix, "Distance selection ERROR (continuing):", e$message, "\n")
        write(paste0(Sys.time(), " ", msg_prefix, "Distance selection ERROR: ", e$message, "\n"),
              file = log_file, append = TRUE)
        panel_bu <<- character(0)
        panel_td <<- character(0)
      })
      
      global_log_event(gse_id, "Distance selection: done")
      
      # AUC comparison 
      
      global_log_event(gse_id, sprintf("AUC comparison: start (n_splits=%d)", n_splits_auc))
      res_auc <- compare_methods_auc_glmnet_noleak(
        mat,
        group_ids = group_ids,
        k_panel = k_panel,
        n_splits = n_splits_auc,
        train_frac = 0.7,
        seed = 123,
        dist_prefilter_topN = dist_prefilter_topN,
        dist_B = dist_B,
        dist_sub = dist_sub,
        dist_k = dist_k,
        dist_p = dist_p,
        dist_pi_thr = dist_pi_thr,
        dist_eps = dist_eps,
        dist_min_n = dist_min_n,
        td_B = td_B,
        td_sub = td_sub,
        td_q = td_q,
        td_p = td_p,
        td_pi_thr = td_pi_thr,
        td_eps = td_eps,
        td_min_n = td_min_n
      )
      write.csv(res_auc,
                file.path(results_dir, paste0(gse_id, "_AUC_comparison_DE_vs_ReliefF_vs_mRMR_vs_BottomUp_vs_TopDown_vs_ElasticNet_vs_LASSO.csv")),
                row.names = FALSE)
      global_log_event(gse_id, "AUC comparison: done")
      
      # ---- Volcano plot ---
      
      out_vol <- file.path(results_dir, paste0(gse_id, "_volcano_panels_ALL_methods.pdf"))
      tryCatch({
        make_volcano_plot(de, list(DE = panel_de, ReliefF = panel_relief, mRMR = panel_mrmr, BottomUp = panel_bu, TopDown = panel_td, ElasticNet = panel_en, LASSO = panel_lasso), out_vol)
        global_log_event(gse_id, "Volcano plot: done")
      }, error = function(e) {
        cat(Sys.time(), msg_prefix, "Volcano plot ERROR:", e$message, "\n")
        write(paste0(Sys.time(), " ", msg_prefix, "Volcano plot ERROR: ", e$message, "\n"),
              file = log_file, append = TRUE)
      })
      
      cat(Sys.time(), msg_prefix, "Done successfully.\n")
      tryCatch(.update_id_list(completed_file, gse_id), error = function(e) {})
      tryCatch(global_log_event(gse_id, "DONE dataset"), error = function(e) {})
      
    }, timeout = timeout_sec_per_gse)
    
  }, error = function(e) {
    errmsg <- paste0(Sys.time(), " ", msg_prefix, "ERROR: ", e$message, "\n")
    cat(errmsg)
    write(errmsg, file = log_file, append = TRUE)
    tryCatch(global_log_event(gse_id, paste0("ERROR dataset: ", e$message)), error = function(e2) {})
    tryCatch(.update_id_list(completed_file, gse_id), error = function(e2) {})
  })
}

#  Safety: limit implicit threading 

.safe_set_single_thread <- function() {
  
  # Limit common BLAS/OpenMP thread envs (works even if RhpcBLASctl is not installed)
  
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    suppressWarnings(try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE))
    suppressWarnings(try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE))
  }
}


.safe_read_lines <- function(fp) {
  if (!file.exists(fp)) return(character(0))
  x <- tryCatch(readLines(fp, warn = FALSE), error = function(e) character(0))
  x <- trimws(x)
  x[nzchar(x)]
}


.safe_write_lines <- function(lines, fp) {
  .ensure_tempdir()
  tmp <- paste0(fp, ".tmp")
  writeLines(lines, tmp, useBytes = TRUE)
  if (file.exists(fp)) file.remove(fp)
  file.rename(tmp, fp)
}

.with_global_lock <- function(expr) {
  lk <- NULL
  lk <- filelock::lock(global_log_lockfile, timeout = 60000)
  on.exit({ if (!is.null(lk)) filelock::unlock(lk) }, add = TRUE)
  force(expr)
}

.update_id_list <- function(fp, id) {
  .with_global_lock({
    ids <- .safe_read_lines(fp)
    if (!(id %in% ids)) {
      ids <- c(ids, id)
      .safe_write_lines(ids, fp)
    }
  })
}

.get_progress <- function() {
  done <- .safe_read_lines(completed_file)
  ndone <- length(unique(done))
  ntot  <- length(gse_list)
  pct   <- if (ntot > 0) round(100 * ndone / ntot, 1) else NA_real_
  list(done = ndone, total = ntot, pct = pct)
}

global_log_event <- function(gse_id, event) {
  .with_global_lock({
    pr <- .get_progress()
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- sprintf("%s | %d/%d (%.1f%%) | %s | %s", ts, pr$done, pr$total, pr$pct, gse_id, event)
    
    prev <- .safe_read_lines(global_log_file)
    out <- c(prev, line)
    if (length(out) > GLOBAL_LOG_TAIL_N) out <- tail(out, GLOBAL_LOG_TAIL_N)
    
    .safe_write_lines(out, global_log_file)
  })
}


.ensure_tempdir <- function() {
  td <- tempdir()
  if (!dir.exists(td)) {
    dir.create(td, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(td)
}

# -------------------- Robust sample labeling --------------------

get_sample_labels <- function(pheno_df) {
  tumor_keywords  <- "\\bcancer\\b|\\btumou?r\\b|\\bcarcinoma\\b|\\bmalignant\\b|\\bneoplasm\\b|\\blesion\\b|\\badenocarcinoma\\b"
  normal_keywords <- "\\bnormal\\b|\\bcontrol\\b|\\bhealthy\\b|\\bnon-?tumou?r\\b|\\bnon-?neoplastic\\b|\\bbenign\\b|\\badjacent\\b"
  
  cols <- grep("characteristics|source_name|title", colnames(pheno_df),
               value = TRUE, ignore.case = TRUE)
  
  if (length(cols) == 0) {
    warning("No phenotype text columns matched (characteristics/source_name/title).")
    return(rep(NA_character_, nrow(pheno_df)))
  }
  
  text_blob <- apply(
    pheno_df[, cols, drop = FALSE],
    1,
    function(row) tolower(paste(na.omit(row), collapse = " "))
  )
  
  labels <- rep(NA_character_, nrow(pheno_df))
  is_tumor  <- grepl(tumor_keywords,  text_blob, ignore.case = TRUE)
  is_normal <- grepl(normal_keywords, text_blob, ignore.case = TRUE)
  
  labels[is_tumor]  <- "tumor"
  labels[is_normal] <- "normal"  # normal overwrites tumor if both match
  
  if (any(is.na(labels))) warning("Could not assign labels for all samples. Check phenotype keywords.")
  labels
}

# -------------------- Expression preprocessing --------------------

preprocess_expression_data <- function(expr_matrix) {
  expr_vals <- as.vector(expr_matrix)
  expr_vals <- expr_vals[!is.na(expr_vals) & expr_vals > 0]
  if (length(expr_vals) == 0) return(expr_matrix)
  
  is_whole_numbers <- all(abs(expr_vals - round(expr_vals)) < 1e-10, na.rm = TRUE)
  has_large_values <- any(expr_vals > 50, na.rm = TRUE)
  
  if (is_whole_numbers && has_large_values) {
    cat("Detected raw count-like data, applying log2(count + 1)\n")
    expr_matrix <- log2(expr_matrix + 1)
  } else {
    has_negative <- any(expr_matrix < 0, na.rm = TRUE)
    if (!has_negative && max(expr_vals, na.rm = TRUE) > 20) {
      cat("Applying log2(x + 1) to expression data\n")
      expr_matrix <- log2(expr_matrix + 1)
    }
  }
  expr_matrix
}

# -------------------- Utility: infer miRNA column from GPL table -------

infer_mirna_col <- function(anno_df) {
  prefs <- c("miRNA_ID", "miRNA", "MIRNA", "mirna", "miRNA.id", "miRNA_ID_REF")
  hit <- prefs[prefs %in% colnames(anno_df)]
  if (length(hit) > 0) return(hit[1])
  
  cand <- grep("mirna|microRNA|miR|MIR", colnames(anno_df), value = TRUE, ignore.case = TRUE)
  if (length(cand) > 0) return(cand[1])
  
  NA_character_
}

# -------------------- Build miRNA-by-sample matrix --------------------

build_mirna_matrix <- function(eset) {
  expr <- exprs(eset)
  pheno <- pData(eset)
  
  expr <- preprocess_expression_data(expr)
  pheno$label <- get_sample_labels(pheno)
  
  label_counts <- table(pheno$label)
  if (length(label_counts) < 2 || any(label_counts < 3)) {
    stop("Insufficient labeled samples (need >= 3 per class).")
  }
  
  gpl_id <- annotation(eset)
  gpl <- getGEO(gpl_id)
  anno <- Table(gpl)
  
  if (!("ID" %in% colnames(anno))) stop("GPL annotation has no 'ID' column.")
  
  mir_col <- infer_mirna_col(anno)
  if (is.na(mir_col)) stop("Could not infer miRNA column from GPL annotation table.")
  
  acc_col <- grep("MIMAT|Accession", colnames(anno), value = TRUE)[1]
  seq_col <- grep("Sequence|SEQUENCE", colnames(anno), value = TRUE)[1]
  
  keep_cols <- unique(na.omit(c("ID", mir_col, acc_col, seq_col)))
  map <- anno[, keep_cols, drop = FALSE]
  colnames(map)[colnames(map) == "ID"] <- "probe_id"
  colnames(map)[colnames(map) == mir_col] <- "miRNA"
  map$probe_id <- as.character(map$probe_id)
  
  expr_df <- rownames_to_column(as.data.frame(expr), var = "probe_id")
  expr_df$probe_id <- as.character(expr_df$probe_id)
  
  expr_annot <- suppressMessages(left_join(map, expr_df, by = "probe_id"))
  
  expr_by_miRNA <- expr_annot |>
    filter(!is.na(miRNA) & miRNA != "") |>
    group_by(miRNA) |>
    summarise(across(where(is.numeric), median, na.rm = TRUE), .groups = "drop")
  
  mat <- as.data.frame(t(as.matrix(expr_by_miRNA[, -1, drop = FALSE])))
  colnames(mat) <- expr_by_miRNA$miRNA
  
  mat$label <- pheno$label[match(rownames(mat), rownames(pheno))]
  
  list(mat = mat, pheno = pheno, gpl_id = gpl_id)
}


get_subject_ids_for_gse <- function(gse_id, sample_ids) {
  sample_ids <- as.character(sample_ids)
  
  paired_gse <- c(
    "GSE10694", "GSE25508", "GSE34535", "GSE34536",
    "GSE45666", "GSE54751", "GSE76260", "GSE102286"
  )
  
  # Unpaired datasets: each sample is its own independent group.
  
  if (!(gse_id %in% paired_gse)) {
    return(sample_ids)
  }
  
  gse_obj <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  if (!is.list(gse_obj)) gse_obj <- list(gse_obj)
  
  overlaps <- vapply(
    gse_obj,
    function(eset) {
      ph <- Biobase::pData(eset)
      sum(sample_ids %in% rownames(ph))
    },
    numeric(1)
  )
  
  eset <- gse_obj[[which.max(overlaps)]]
  pheno <- Biobase::pData(eset)
  
  if (!("title" %in% colnames(pheno))) {
    stop("Cannot construct patient groups for ", gse_id, ": GEO phenotype data has no title column.")
  }
  
  titles <- as.character(pheno$title)
  names(titles) <- rownames(pheno)
  titles <- titles[match(sample_ids, names(titles))]
  
  if (anyNA(titles)) {
    stop("Cannot construct patient groups for ", gse_id, ": some cached samples were not found in GEO phenotype metadata.")
  }
  
  # Default: unique sample-level groups. Dataset-specific rules below only join
  # samples when GEO metadata identifies the same patient/specimen.
  
  ids <- sample_ids
  
  if (gse_id == "GSE10694") {
    m <- regexec("patient\\s*([0-9]+)", titles, ignore.case = TRUE)
    z <- regmatches(titles, m)
    ok <- lengths(z) >= 2L
    ids[ok] <- paste0("patient_", vapply(z[ok], `[`, character(1), 2L))
  }
  
  if (gse_id == "GSE25508") {
    
    # Paired tumor/normal samples share the same title apart from the tissue word.
    
    key <- tolower(titles)
    key <- gsub("_(tumou?r|normal|norma)", "", key, ignore.case = TRUE)
    key <- gsub("\\s+", " ", trimws(key))
    ids <- paste0("subject_", key)
  }
  
  if (gse_id %in% c("GSE34535", "GSE34536")) {
    # Titles are 1_BCC ... 7_BCC / 1_control ... 7_control, or SCC analogues.
    m <- regexec("^([0-9]+)_", titles)
    z <- regmatches(titles, m)
    ok <- lengths(z) >= 2L
    ids[ok] <- paste0("patient_", vapply(z[ok], `[`, character(1), 2L))
  }
  
  if (gse_id == "GSE45666") {
    # Join adjacent-normal/tumor samples and technical replicates by sample code.
    # Examples include BreastTumor-S941-Rep1 / BreastTumor-S941-Rep2.
    key <- tolower(titles)
    key <- sub("^breast(tumou?r|adjacentnormal|normal)-", "", key, ignore.case = TRUE)
    key <- sub("-rep[0-9]+$", "", key, ignore.case = TRUE)
    ids <- paste0("subject_", key)
  }
  
  if (gse_id == "GSE54751") {
    # GEO lists the 20 specimens as 10 consecutive tumor/non-tumor pairs.
    ord <- match(sample_ids, rownames(pheno))
    ids <- paste0("pair_", ceiling(ord / 2))
  }
  
  if (gse_id == "GSE76260") {
    # Titles contain e.g. "tumor tissue patient p2" / "normal tissue patient p2".
    m <- regexec("patient\\s+(p[0-9]+)", titles, ignore.case = TRUE)
    z <- regmatches(titles, m)
    ok <- lengths(z) >= 2L
    ids[ok] <- paste0("patient_", tolower(vapply(z[ok], `[`, character(1), 2L)))
  }
  
  if (gse_id == "GSE102286") {
    # Titles contain patient codes such as [10097N] and [10097T].
    m <- regexec("\\[([0-9]+)[NT]\\]", titles, ignore.case = TRUE)
    z <- regmatches(titles, m)
    ok <- lengths(z) >= 2L
    ids[ok] <- paste0("patient_", vapply(z[ok], `[`, character(1), 2L))
  }
  
  if (anyNA(ids) || any(!nzchar(ids))) {
    stop("Failed to construct complete patient grouping for ", gse_id, ".")
  }
  
  ids
}


make_grouped_stratified_splits <- function(
    ybin,
    group_ids,
    n_splits,
    train_frac,
    seed
) {
  if (length(group_ids) != length(ybin)) {
    stop("group_ids must have the same length as ybin.")
  }
  
  set.seed(seed)
  group_ids <- as.character(group_ids)
  groups <- split(seq_along(ybin), group_ids)
  
  # Stratify groups by the class composition they contain. Thus ordinary
  # unpaired samples are stratified by class, while matched pairs form their
  # own stratum and are always allocated together.
  group_stratum <- vapply(
    groups,
    function(ii) paste(sort(unique(ybin[ii])), collapse = ""),
    character(1)
  )
  
  strata <- split(names(groups), group_stratum)
  splits <- vector("list", n_splits)
  
  for (i in seq_len(n_splits)) {
    train_groups <- character(0)
    
    for (s in names(strata)) {
      gs <- strata[[s]]
      ng <- length(gs)
      
      if (ng == 1L) {
        n_train <- if (train_frac >= 0.5) 1L else 0L
      } else {
        n_train <- floor(train_frac * ng)
        n_train <- max(1L, min(ng - 1L, n_train))
      }
      
      if (n_train > 0L) {
        train_groups <- c(
          train_groups,
          sample(gs, size = n_train, replace = FALSE)
        )
      }
    }
    
    trn <- sort(unlist(groups[train_groups], use.names = FALSE))
    tst <- setdiff(seq_along(ybin), trn)
    
    if (
      length(unique(ybin[trn])) < 2L ||
      length(unique(ybin[tst])) < 2L
    ) {
      stop("Grouped split produced a train or test set lacking one class.")
    }
    
    splits[[i]] <- list(trn = trn, tst = tst, train = trn, test = tst)
  }
  
  splits
}

