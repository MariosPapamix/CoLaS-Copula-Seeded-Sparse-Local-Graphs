## ============================================================
## E5: Robustness & model-selection sensitivity (SELF-CONTAINED)
## Includes FIX for calibrate_lambda signature mismatch.
## Also avoids "double LCC conditioning":
##   - extract observed LCC once
##   - set use_lcc=FALSE throughout E5 evaluation
##
## Outputs:
##   - observed_edgelist.txt
##   - E5_results.csv
## ============================================================

## --------------------------
## 0) Packages
## --------------------------
required_pkgs <- c("igraph", "copula", "Matrix")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(igraph)
library(copula)
library(Matrix)

set.seed(123)

## --------------------------
## 1) Helpers: LCC + torus geometry
## --------------------------
largest_cc <- function(g) {
  g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
  if (vcount(g) == 0) return(g)
  comp <- components(g)
  keep <- which(comp$membership == which.max(comp$csize))
  induced_subgraph(g, keep)
}

torus_diff <- function(x, y) {
  d <- x - y
  d - round(d)  # wrap to [-1/2,1/2)
}

## --------------------------
## 2) Kernels (support radius 1)
## --------------------------
k_ball <- function(u) as.numeric(rowSums(u^2) <= 1)
k_epanechnikov <- function(u) pmax(1 - rowSums(u^2), 0)
k_triangular <- function(u) pmax(1 - sqrt(rowSums(u^2)), 0)

get_kernel <- function(name) {
  switch(name,
         ball = k_ball,
         epanechnikov = k_epanechnikov,
         triangular = k_triangular,
         stop("Unknown kernel: ", name))
}

## --------------------------
## 3) Copula (bivariate) + empirical W quantile map
## --------------------------
make_bicop <- function(family, tau) {
  family <- match.arg(family, c("gaussian", "frank", "clayton", "gumbel"))
  
  # clip tau slightly away from +/-1 for stability
  if (family %in% c("gaussian", "frank")) {
    tau <- max(min(tau, 0.95), -0.95)
  } else {
    tau <- max(min(tau, 0.95), 0.0)
  }
  
  cop <- switch(
    family,
    gaussian = normalCopula(param = 0.0, dim = 2, dispstr = "un"),
    frank    = frankCopula(param = 1.0, dim = 2),
    clayton  = claytonCopula(param = 1.0, dim = 2),
    gumbel   = gumbelCopula(param = 2.0, dim = 2)
  )
  
  cop@parameters <- iTau(cop, tau)
  cop
}

emp_q_vec <- function(u, x_sorted) {
  n <- length(x_sorted)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  idx <- ceiling(u * n)
  idx[idx < 1] <- 1
  idx[idx > n] <- n
  x_sorted[idx]
}

sample_latents <- function(n, d, cop_family, tau_dep, qW, tau_trunc = Inf) {
  cop <- make_bicop(cop_family, tau_dep)
  UV <- rCopula(n, cop)
  U  <- UV[, 1]
  V1 <- UV[, 2]
  
  W <- qW(U)
  if (is.finite(tau_trunc)) W <- pmin(W, tau_trunc)
  W <- pmax(W, 1e-12)
  
  if (d == 1L) {
    X <- matrix(V1, ncol = 1)
  } else {
    Xrest <- matrix(runif(n * (d - 1L)), ncol = d - 1L)
    X <- cbind(V1, Xrest)
  }
  
  list(W = W, X = X)
}

## --------------------------
## 4) Simulator: CoLaS (fixed) and CoLaS-HT (HT)
## Uses cell hashing on the torus for speed.
## rho_n = n * eps_n^d -> rho; eps_n = (rho / n)^(1/d)
## --------------------------
simulate_colaS <- function(W, X, d,
                           regime = c("fixed", "HT"),
                           kernel_name = "ball",
                           rho = 1,
                           lambda = 1,
                           seed = NULL) {
  regime <- match.arg(regime)
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(W)
  k_fn <- get_kernel(kernel_name)
  
  eps_n <- (rho / n)^(1 / d)
  eps_n <- min(eps_n, 0.5)   # safety cap on torus scale
  
  # Max interaction radius for hashing
  if (regime == "fixed") {
    r_max <- eps_n
  } else {
    maxW <- max(W)
    r_max <- eps_n * (maxW^2)^(1 / d)
    r_max <- min(r_max, 1)
  }
  
  cell_size <- max(r_max, 1e-8)
  m <- max(1L, ceiling(1 / cell_size))     # cells per dimension
  base <- m^(0:(d - 1L))
  
  idx <- floor(X * m)
  idx[idx == m] <- m - 1L
  cell_id <- as.integer(idx %*% base) + 1L
  
  cells <- split(seq_len(n), cell_id)
  cell_map <- vector("list", m^d)
  cell_map[as.integer(names(cells))] <- cells
  nonempty <- sort(as.integer(names(cells)))
  
  decode_cell <- function(cid) {
    cid0 <- cid - 1L
    out <- integer(d)
    for (k in seq_len(d)) {
      out[k] <- cid0 %% m
      cid0 <- cid0 %/% m
    }
    out
  }
  encode_cell <- function(idx_vec) as.integer(sum((idx_vec %% m) * base)) + 1L
  
  offsets <- as.matrix(expand.grid(rep(list(-1:1), d)))
  
  from <- integer(0)
  to   <- integer(0)
  
  for (c in nonempty) {
    nodes_i <- cell_map[[c]]
    if (length(nodes_i) == 0) next
    idx_c <- decode_cell(c)
    
    for (r in seq_len(nrow(offsets))) {
      idx2 <- (idx_c + offsets[r, ]) %% m
      c2 <- encode_cell(idx2)
      if (c2 < c) next
      
      nodes_j <- cell_map[[c2]]
      if (length(nodes_j) == 0) next
      
      if (c2 == c) {
        if (length(nodes_i) < 2) next
        comb <- combn(nodes_i, 2)
        ii <- comb[1, ]
        jj <- comb[2, ]
      } else {
        ii <- rep(nodes_i, times = length(nodes_j))
        jj <- rep(nodes_j, each  = length(nodes_i))
      }
      
      dx <- torus_diff(X[ii, , drop = FALSE], X[jj, , drop = FALSE])
      dist <- sqrt(rowSums(dx^2))
      
      if (regime == "fixed") {
        keep <- dist <= eps_n
        if (!any(keep)) next
        ii <- ii[keep]; jj <- jj[keep]
        dx <- dx[keep, , drop = FALSE]
        
        u <- dx / eps_n
        k_val <- k_fn(u)
        p <- 1 - exp(-(lambda / rho) * W[ii] * W[jj] * k_val)
        
      } else {
        rad <- eps_n * (W[ii] * W[jj])^(1 / d)
        keep <- dist <= rad
        if (!any(keep)) next
        ii <- ii[keep]; jj <- jj[keep]
        dx <- dx[keep, , drop = FALSE]
        rad <- rad[keep]
        
        u <- dx / rad
        k_val <- k_fn(u)
        p <- 1 - exp(-(lambda / rho) * k_val)
      }
      
      keep_edge <- runif(length(p)) < p
      if (any(keep_edge)) {
        from <- c(from, ii[keep_edge])
        to   <- c(to,   jj[keep_edge])
      }
    }
  }
  
  g <- make_empty_graph(n, directed = FALSE)
  if (length(from) > 0) {
    g <- add_edges(g, as.vector(rbind(from, to)))
  }
  simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
}

## --------------------------
## 5) Stats (matched + held-out)
## --------------------------
safe_transitivity <- function(g) {
  x <- suppressWarnings(transitivity(g, type = "global"))
  if (is.na(x) || is.nan(x)) 0 else x
}
safe_assort <- function(g) {
  x <- suppressWarnings(assortativity_degree(g, directed = FALSE))
  if (is.na(x) || is.nan(x)) 0 else x
}

degree_clust_curve <- function(g) {
  deg <- degree(g)
  cl  <- transitivity(g, type = "local", isolates = "zero")
  ck  <- tapply(cl, deg, mean, na.rm = TRUE)
  ck[is.na(ck)] <- 0
  ck
}

binned_jdd <- function(g, n_bins = 20) {
  deg <- degree(g)
  n <- length(deg)
  if (ecount(g) == 0) return(matrix(0, n_bins, n_bins))
  
  r <- rank(deg, ties.method = "random")
  bins <- ceiling(r / n * n_bins)
  bins[bins < 1] <- 1
  bins[bins > n_bins] <- n_bins
  
  E <- ends(g, E(g), names = FALSE)
  b1 <- bins[E[, 1]]
  b2 <- bins[E[, 2]]
  
  M <- matrix(0, n_bins, n_bins)
  for (k in seq_len(nrow(E))) {
    i <- b1[k]; j <- b2[k]
    M[i, j] <- M[i, j] + 1
    M[j, i] <- M[j, i] + 1
  }
  s <- sum(M)
  if (s > 0) M / s else M
}

core_hist <- function(g) {
  c <- coreness(g)
  tab <- table(c)
  as.numeric(tab) / sum(tab)
}

spectral_radius_power <- function(g, iters = 50) {
  if (ecount(g) == 0) return(0)
  A <- as_adjacency_matrix(g, sparse = TRUE)
  n <- vcount(g)
  x <- runif(n)
  x <- x / sqrt(sum(x^2))
  for (t in seq_len(iters)) {
    x_new <- as.vector(A %*% x)
    nrm <- sqrt(sum(x_new^2))
    if (nrm == 0) return(0)
    x <- x_new / nrm
  }
  as.numeric(crossprod(x, A %*% x))
}

prep_for_stats <- function(g, use_lcc = TRUE) {
  g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
  if (use_lcc) g <- largest_cc(g)
  g
}

extract_stats <- function(g, n_bins = 20, use_lcc = TRUE) {
  g2 <- prep_for_stats(g, use_lcc = use_lcc)
  deg <- degree(g2)
  list(
    mean_deg = if (vcount(g2) > 0) mean(deg) else 0,
    trans    = safe_transitivity(g2),
    assort   = safe_assort(g2),
    
    # held-out:
    deg      = deg,
    ck       = degree_clust_curve(g2),
    jdd      = binned_jdd(g2, n_bins = n_bins),
    core     = core_hist(g2),
    spec_rad = spectral_radius_power(g2)
  )
}

## --------------------------
## 6) GOF distances
## --------------------------
ks_dist_deg <- function(d1, d2) {
  x <- sort(d1); y <- sort(d2)
  allv <- sort(unique(c(x, y)))
  Fx <- ecdf(x); Fy <- ecdf(y)
  max(abs(Fx(allv) - Fy(allv)))
}
ck_dist <- function(ck1, ck2) {
  k <- sort(unique(c(as.integer(names(ck1)), as.integer(names(ck2)))))
  v1 <- ck1[as.character(k)]
  v2 <- ck2[as.character(k)]
  v1[is.na(v1)] <- 0
  v2[is.na(v2)] <- 0
  sqrt(mean((v1 - v2)^2))
}
l1_dist_mat <- function(A, B) sum(abs(A - B))
l1_dist_hist <- function(h1, h2) {
  m <- max(length(h1), length(h2))
  a <- rep(0, m); b <- rep(0, m)
  a[seq_along(h1)] <- h1
  b[seq_along(h2)] <- h2
  sum(abs(a - b))
}

gof_from_sims <- function(obs, sims,
                          w = list(degKS = 1, ck = 1, jdd = 1, core = 1, spec = 1)) {
  per_sim <- lapply(sims, function(st) {
    degKS <- ks_dist_deg(obs$deg, st$deg)
    ckE   <- ck_dist(obs$ck, st$ck)
    jddE  <- l1_dist_mat(obs$jdd, st$jdd)
    coreE <- l1_dist_hist(obs$core, st$core)
    specE <- if (is.na(obs$spec_rad) || obs$spec_rad == 0) {
      abs(st$spec_rad - obs$spec_rad)
    } else {
      abs(st$spec_rad - obs$spec_rad) / obs$spec_rad
    }
    score <- w$degKS*degKS + w$ck*ckE + w$jdd*jddE + w$core*coreE + w$spec*specE
    c(degKS = degKS, ck = ckE, jdd = jddE, core = coreE, specE = specE, score = score)
  })
  M <- do.call(rbind, per_sim)
  colMeans(M)
}

## --------------------------
## 7) Calibration: lambda by bisection, and rho candidates
## (FIX: signature includes lambda0/lambda_max, plus ... for robustness)
## --------------------------
calibrate_lambda <- function(target_kbar, sim_mean_deg_fn,
                             lambda0 = 1,
                             lambda_max = 1e6,
                             B = 1,
                             iters = 8,
                             ...) {
  mean_deg_at <- function(lam) {
    ks <- replicate(B, sim_mean_deg_fn(lam))
    mean(ks)
  }
  
  lo <- 1e-10
  hi <- lambda0
  k_hi <- mean_deg_at(hi)
  
  while (k_hi < target_kbar && hi < lambda_max) {
    hi <- hi * 2
    k_hi <- mean_deg_at(hi)
  }
  
  if (k_hi < target_kbar) {
    return(list(feasible = FALSE, lambda = NA_real_, k_hi = k_hi))
  }
  
  for (t in seq_len(iters)) {
    mid <- sqrt(lo * hi)
    k_mid <- mean_deg_at(mid)
    if (k_mid < target_kbar) lo <- mid else hi <- mid
  }
  
  list(feasible = TRUE, lambda = hi, k_hi = k_hi)
}

vol_unit_ball <- function(d) {
  pi^(d/2) / gamma(d/2 + 1)
}

rho_candidates_auto <- function(d, regime, kbar_target,
                                safety = 1.25,
                                rho_floor = 0.25,
                                rho_cap = 128) {
  v <- vol_unit_ball(d)
  
  if (regime == "fixed") {
    rho0 <- max(rho_floor, safety * kbar_target / max(v, 1e-12))
  } else {
    rho0 <- max(rho_floor, 0.75 * kbar_target / max(v, 1e-12))
  }
  
  rho0 <- min(rho0, rho_cap)
  unique(pmin(rho_cap, pmax(rho_floor, rho0 * c(0.5, 1, 2, 4))))
}

tau_grid_default <- function(family) {
  if (family %in% c("clayton", "gumbel")) {
    unique(c(0, 0.2, 0.4, 0.6, 0.8))
  } else {
    unique(c(-0.6, -0.3, 0, 0.3, 0.6))
  }
}

## --------------------------
## 8) Composite-likelihood IC (optional)
## FIX: d=1 uses 2D FR then angle -> 1D coordinate
## --------------------------
rank_to_unit <- function(v) {
  n <- length(v)
  (rank(v, ties.method = "random") - 0.5) / n
}

estimate_positions_fr_unit <- function(g, d) {
  if (vcount(g) == 0) return(matrix(numeric(0), ncol = d))
  
  if (d == 1L) {
    coords2 <- layout_with_fr(g, dim = 2)
    ang <- atan2(coords2[, 2], coords2[, 1])     # (-pi, pi]
    x1  <- (ang + pi) / (2 * pi)                 # [0,1)
    return(matrix(rank_to_unit(x1), ncol = 1))
  }
  
  coords <- layout_with_fr(g, dim = d)           # d is 2 or 3
  X <- apply(coords, 2, rank_to_unit)
  X
}

pseudo_cl_ic <- function(g_obs, d, regime, kernel,
                         rho_hat, lambda_hat,
                         cop_family, tau_hat,
                         tau_trunc,
                         nonedge_multiplier = 5L,
                         seed = 1,
                         use_lcc = TRUE) {
  set.seed(seed)
  
  g0 <- prep_for_stats(g_obs, use_lcc = use_lcc)
  n <- vcount(g0)
  m <- ecount(g0)
  if (n < 3 || m == 0) {
    return(list(cl_loglik = NA_real_, cl_aic = NA_real_, cl_bic = NA_real_))
  }
  
  k_fn <- get_kernel(kernel)
  eps_n <- (rho_hat / n)^(1 / d)
  eps_n <- min(eps_n, 0.5)
  
  deg <- degree(g0)
  w_hat <- (deg + 1) / mean(deg + 1)
  if (is.finite(tau_trunc)) w_hat <- pmin(w_hat, tau_trunc)
  w_hat <- pmax(w_hat, 1e-12)
  
  X_hat <- estimate_positions_fr_unit(g0, d)
  
  # copula density term for (U_hat, V1_hat)
  U_hat  <- rank_to_unit(w_hat)
  V1_hat <- rank_to_unit(X_hat[, 1])
  cop <- make_bicop(cop_family, tau_hat)
  dens <- dCopula(cbind(U_hat, V1_hat), cop)
  ll_cop <- sum(log(pmax(dens, 1e-12)))
  
  A <- as_adjacency_matrix(g0, sparse = TRUE)
  E <- ends(g0, E(g0), names = FALSE)
  ii_e <- E[, 1]; jj_e <- E[, 2]
  
  Npairs <- n * (n - 1) / 2
  Nne <- Npairs - m
  target_non <- min(as.integer(nonedge_multiplier * m), as.integer(Nne))
  
  ii_ne <- integer(0); jj_ne <- integer(0)
  tries <- 0L
  while (length(ii_ne) < target_non && tries < target_non * 100L) {
    tries <- tries + 1L
    i <- sample.int(n, 1L)
    j <- sample.int(n, 1L)
    if (i == j) next
    if (i > j) { tmp <- i; i <- j; j <- tmp }
    if (A[i, j] != 0) next
    ii_ne <- c(ii_ne, i); jj_ne <- c(jj_ne, j)
  }
  
  pair_prob <- function(ii, jj) {
    dx <- torus_diff(X_hat[ii, , drop = FALSE], X_hat[jj, , drop = FALSE])
    dist <- sqrt(rowSums(dx^2))
    
    if (regime == "fixed") {
      keep <- dist <= eps_n
      p <- rep(0, length(ii))
      if (any(keep)) {
        u <- dx[keep, , drop = FALSE] / eps_n
        k_val <- k_fn(u)
        p[keep] <- 1 - exp(-(lambda_hat / rho_hat) * w_hat[ii[keep]] * w_hat[jj[keep]] * k_val)
      }
      p
    } else {
      rad <- eps_n * (w_hat[ii] * w_hat[jj])^(1 / d)
      keep <- dist <= rad
      p <- rep(0, length(ii))
      if (any(keep)) {
        u <- dx[keep, , drop = FALSE] / rad[keep]
        k_val <- k_fn(u)
        p[keep] <- 1 - exp(-(lambda_hat / rho_hat) * k_val)
      }
      p
    }
  }
  
  p_e <- pmax(pmin(pair_prob(ii_e, jj_e), 1 - 1e-12), 1e-12)
  ll_e <- sum(log(p_e))
  
  ll_ne <- 0
  if (length(ii_ne) > 0) {
    p_ne <- pmax(pmin(pair_prob(ii_ne, jj_ne), 1 - 1e-12), 1e-12)
    w_ne <- Nne / length(p_ne)
    ll_ne <- w_ne * sum(log(1 - p_ne))
  }
  
  cl <- ll_cop + ll_e + ll_ne
  k_par <- 3  # (tau_hat, rho_hat, lambda_hat)
  aic <- -2 * cl + 2 * k_par
  bic <- -2 * cl + log(Npairs) * k_par
  
  list(cl_loglik = cl, cl_aic = aic, cl_bic = bic)
}

## --------------------------
## 9) Fit (tau, rho, lambda) for one spec
## Matches: mean degree (via rho+lambda), transitivity + assortativity (via tau)
## --------------------------
fit_params_for_spec <- function(g_obs,
                                d, regime, kernel, cop_family, tau_trunc,
                                tau_grid = NULL,
                                B_fit = 3,
                                B_lambda = 2,
                                use_lcc = TRUE,
                                seed = 1) {
  set.seed(seed)
  obs <- extract_stats(g_obs, use_lcc = use_lcc)
  kbar_target <- obs$mean_deg
  
  # empirical proxy for W (degree-based), apply truncation
  deg_obs <- degree(prep_for_stats(g_obs, use_lcc = use_lcc))
  w_emp <- (deg_obs + 1) / mean(deg_obs + 1)
  if (is.finite(tau_trunc)) w_emp <- pmin(w_emp, tau_trunc)
  w_sorted <- sort(w_emp)
  qW_emp <- function(u) emp_q_vec(u, w_sorted)
  
  if (is.null(tau_grid)) tau_grid <- tau_grid_default(cop_family)
  rho_grid <- rho_candidates_auto(d, regime, kbar_target)
  
  best <- list(loss = Inf, tau = NA_real_, rho = NA_real_, lambda = NA_real_)
  
  for (rho_try in rho_grid) {
    for (tau_try in tau_grid) {
      
      sim_mean_deg <- function(lambda) {
        n_sim <- vcount(prep_for_stats(g_obs, use_lcc = use_lcc))
        lat <- sample_latents(n_sim, d, cop_family, tau_try,
                              qW = qW_emp, tau_trunc = tau_trunc)
        g_sim <- simulate_colaS(lat$W, lat$X, d,
                                regime = regime,
                                kernel_name = kernel,
                                rho = rho_try,
                                lambda = lambda,
                                seed = seed + sample.int(1e9, 1))
        st <- extract_stats(g_sim, use_lcc = use_lcc)
        st$mean_deg
      }
      
      lam_fit <- calibrate_lambda(kbar_target, sim_mean_deg,
                                  lambda0 = 1,
                                  lambda_max = 1e6,
                                  B = B_lambda,
                                  iters = 8)
      if (!lam_fit$feasible) next
      
      # simulate B_fit graphs to estimate trans & assort under fitted (rho, lambda, tau)
      trans_vals <- numeric(B_fit)
      assort_vals <- numeric(B_fit)
      
      for (b in seq_len(B_fit)) {
        n_sim <- vcount(prep_for_stats(g_obs, use_lcc = use_lcc))
        lat <- sample_latents(n_sim, d, cop_family, tau_try,
                              qW = qW_emp, tau_trunc = tau_trunc)
        g_sim <- simulate_colaS(lat$W, lat$X, d,
                                regime = regime,
                                kernel_name = kernel,
                                rho = rho_try,
                                lambda = lam_fit$lambda,
                                seed = seed + b + sample.int(1e9, 1))
        st <- extract_stats(g_sim, use_lcc = use_lcc)
        trans_vals[b] <- st$trans
        assort_vals[b] <- st$assort
      }
      
      trans_bar <- mean(trans_vals)
      assort_bar <- mean(assort_vals)
      
      loss <- (trans_bar - obs$trans)^2 + (assort_bar - obs$assort)^2
      
      if (is.finite(loss) && loss < best$loss) {
        best <- list(loss = loss,
                     tau = tau_try,
                     rho = rho_try,
                     lambda = lam_fit$lambda)
      }
    }
  }
  
  best
}

## --------------------------
## 10) Fit + GOF + optional IC for one E5 row
## --------------------------
fit_and_eval_spec <- function(g_obs, spec,
                              B_fit = 3,
                              B_gof = 30,
                              B_lambda = 2,
                              n_bins = 20,
                              use_lcc = TRUE,
                              use_clic = TRUE,
                              seed = 1) {
  
  d <- spec$d
  regime <- spec$regime
  kernel <- spec$kernel
  cop_family <- spec$copula
  tau_trunc <- spec$tau_trunc
  
  fit <- fit_params_for_spec(g_obs,
                             d = d, regime = regime, kernel = kernel,
                             cop_family = cop_family, tau_trunc = tau_trunc,
                             tau_grid = NULL,
                             B_fit = B_fit,
                             B_lambda = B_lambda,
                             use_lcc = use_lcc,
                             seed = seed)
  
  if (!is.finite(fit$loss) || is.na(fit$tau) || is.na(fit$rho) || is.na(fit$lambda)) {
    return(data.frame(spec,
                      tau_hat = NA_real_, rho_hat = NA_real_, lambda_hat = NA_real_,
                      fit_loss = NA_real_,
                      gof_score = NA_real_,
                      degKS = NA_real_, ck = NA_real_, jdd = NA_real_, core = NA_real_, specE = NA_real_,
                      cl_loglik = NA_real_, cl_aic = NA_real_, cl_bic = NA_real_,
                      stringsAsFactors = FALSE))
  }
  
  obs_stats <- extract_stats(g_obs, n_bins = n_bins, use_lcc = use_lcc)
  
  # rebuild empirical W quantile
  deg_obs <- degree(prep_for_stats(g_obs, use_lcc = use_lcc))
  w_emp <- (deg_obs + 1) / mean(deg_obs + 1)
  if (is.finite(tau_trunc)) w_emp <- pmin(w_emp, tau_trunc)
  w_sorted <- sort(w_emp)
  qW_emp <- function(u) emp_q_vec(u, w_sorted)
  
  sims_stats <- vector("list", B_gof)
  for (b in seq_len(B_gof)) {
    n_sim <- vcount(prep_for_stats(g_obs, use_lcc = use_lcc))
    lat <- sample_latents(n_sim, d, cop_family, fit$tau,
                          qW = qW_emp, tau_trunc = tau_trunc)
    g_sim <- simulate_colaS(lat$W, lat$X, d,
                            regime = regime,
                            kernel_name = kernel,
                            rho = fit$rho,
                            lambda = fit$lambda,
                            seed = seed + 10000 + b + sample.int(1e9, 1))
    sims_stats[[b]] <- extract_stats(g_sim, n_bins = n_bins, use_lcc = use_lcc)
  }
  
  gof <- gof_from_sims(obs_stats, sims_stats)
  
  ic <- list(cl_loglik = NA_real_, cl_aic = NA_real_, cl_bic = NA_real_)
  if (use_clic) {
    ic <- pseudo_cl_ic(g_obs, d, regime, kernel,
                       rho_hat = fit$rho,
                       lambda_hat = fit$lambda,
                       cop_family = cop_family,
                       tau_hat = fit$tau,
                       tau_trunc = tau_trunc,
                       seed = seed,
                       use_lcc = use_lcc)
  }
  
  data.frame(spec,
             tau_hat = fit$tau,
             rho_hat = fit$rho,
             lambda_hat = fit$lambda,
             fit_loss = fit$loss,
             gof_score = unname(gof["score"]),
             degKS = unname(gof["degKS"]),
             ck    = unname(gof["ck"]),
             jdd   = unname(gof["jdd"]),
             core  = unname(gof["core"]),
             specE = unname(gof["specE"]),
             cl_loglik = ic$cl_loglik,
             cl_aic    = ic$cl_aic,
             cl_bic    = ic$cl_bic,
             stringsAsFactors = FALSE)
}

## --------------------------
## 11) E5 runner (full grid)
## --------------------------
run_E5 <- function(g_obs,
                   copula_families = c("gaussian", "frank", "clayton", "gumbel"),
                   ds = c(1L, 2L, 3L),
                   kernels = c("ball", "epanechnikov", "triangular"),
                   regimes = c("fixed", "HT"),
                   tau_trunc = NULL,
                   tau_probs = c(0.99, 0.995),
                   B_fit = 3,
                   B_gof = 30,
                   B_lambda = 2,
                   n_bins = 20,
                   use_lcc = TRUE,
                   use_clic = TRUE,
                   seed = 1) {
  
  g0 <- prep_for_stats(g_obs, use_lcc = use_lcc)
  
  if (is.null(tau_trunc)) {
    deg <- degree(g0)
    w_emp <- (deg + 1) / mean(deg + 1)
    qs <- as.numeric(quantile(w_emp, probs = tau_probs, na.rm = TRUE))
    tau_trunc <- unique(c(Inf, qs))
  }
  
  grid <- expand.grid(
    copula = copula_families,
    d = as.integer(ds),
    kernel = kernels,
    regime = regimes,
    tau_trunc = tau_trunc,
    stringsAsFactors = FALSE
  )
  
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    spec <- grid[i, , drop = FALSE]
    message(sprintf("E5 [%d/%d] copula=%s d=%d kernel=%s regime=%s tau_trunc=%s",
                    i, nrow(grid),
                    spec$copula, spec$d, spec$kernel, spec$regime,
                    ifelse(is.infinite(spec$tau_trunc), "Inf", format(spec$tau_trunc, digits = 6))))
    
    out[[i]] <- fit_and_eval_spec(g0, spec,
                                  B_fit = B_fit,
                                  B_gof = B_gof,
                                  B_lambda = B_lambda,
                                  n_bins = n_bins,
                                  use_lcc = use_lcc,
                                  use_clic = use_clic,
                                  seed = seed + i * 100L)
  }
  
  res <- do.call(rbind, out)
  res[order(res$gof_score), ]
}

## ============================================================
## 12) Load real dataset OR generate synthetic observed network
## ============================================================
edgelist_path <- "observed_edgelist.txt"

if (file.exists(edgelist_path)) {
  message("Loading observed graph from: ", edgelist_path)
  g_obs <- read_graph(edgelist_path, format = "ncol", directed = FALSE)
} else {
  message("No observed_edgelist.txt found; generating a synthetic observed graph...")
  
  # Pareto quantile for W (heavy-ish) and HT regime generator
  qpareto <- function(u, alpha = 3, xmin = 1) {
    u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
    xmin * (1 - u)^(-1 / alpha)
  }
  
  n_obs <- 220
  alpha_true <- 3
  xmin_true <- 1
  mean_W <- alpha_true * xmin_true / (alpha_true - 1)
  qW_true <- function(u) qpareto(u, alpha_true, xmin_true) / mean_W
  
  true_spec <- list(d = 2L, regime = "HT", kernel = "epanechnikov",
                    copula = "frank", tau = 0.4,
                    rho = 1.5)
  
  target_kbar_obs <- 10
  
  sim_true_graph <- function(lambda) {
    lat <- sample_latents(n_obs, true_spec$d, true_spec$copula, true_spec$tau,
                          qW = qW_true, tau_trunc = Inf)
    simulate_colaS(lat$W, lat$X, true_spec$d,
                   regime = true_spec$regime,
                   kernel_name = true_spec$kernel,
                   rho = true_spec$rho,
                   lambda = lambda,
                   seed = sample.int(1e9, 1))
  }
  
  # calibrate lambda for synthetic observed graph (rho fixed at true rho)
  get_kbar <- function(lam) {
    g <- sim_true_graph(lam)
    mean(degree(largest_cc(g)))
  }
  lam_fit <- calibrate_lambda(target_kbar_obs, get_kbar, lambda0 = 1, B = 1, iters = 10)
  lambda_true <- ifelse(lam_fit$feasible, lam_fit$lambda, 10)
  
  g_obs_full <- sim_true_graph(lambda_true)
  g_obs <- g_obs_full
}

## --------------------------
## KEY: Extract observed LCC ONCE, then do NOT LCC-filter sims
## --------------------------
g_obs <- largest_cc(g_obs)
USE_LCC <- FALSE

message("Observed graph summary (LCC already extracted):")
st0 <- extract_stats(g_obs, use_lcc = USE_LCC)
message(sprintf("  n=%d, m=%d, mean degree=%.3f, transitivity=%.3f, assortativity=%.3f",
                vcount(g_obs), ecount(g_obs),
                st0$mean_deg, st0$trans, st0$assort))

# Write out edge list (so your run is reproducible)
write.table(as_edgelist(g_obs, names = FALSE),
            file = edgelist_path, row.names = FALSE, col.names = FALSE)

## ============================================================
## 13) Run E5
## ============================================================
E5_results <- run_E5(
  g_obs,
  copula_families = c("gaussian", "frank", "clayton", "gumbel"),
  ds = c(1L, 2L, 3L),
  kernels = c("ball", "epanechnikov", "triangular"),
  regimes = c("fixed", "HT"),
  tau_trunc = NULL,
  B_fit = 3,
  B_gof = 30,
  B_lambda = 2,
  n_bins = 20,
  use_lcc = USE_LCC,   # IMPORTANT: FALSE (no simulated-LCC conditioning)
  use_clic = TRUE,
  seed = 999
)

write.csv(E5_results, "E5_results.csv", row.names = FALSE)

message("Top 10 (lowest GOF score):")
print(head(E5_results, 10))

message("Saved: observed_edgelist.txt and E5_results.csv")
