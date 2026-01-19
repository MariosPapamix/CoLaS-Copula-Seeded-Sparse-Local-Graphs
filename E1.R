# ============================================================
# Paper E2: Separate knobs / partial-effect matrix  (your E1.R)
# REVISED: adds CoLaS-HT (tail-inheriting) and runs the "degree tail"
#          knob (vary F_W via Pareto alpha) in CoLaS-HT as the paper intends.
#
# Default here: run ALL three knobs in CoLaS-HT so tail-index inference is meaningful.
# (You can switch geo/theta back to "fixed" by editing cfg$regime_geo / cfg$regime_theta.)
#
# Outputs (in cfg$out_dir):
#   figures/
#     E2_knob_sanity.png
#     E2_mean_degree_flatness.png
#     E2_lambda_calibration.png
#     E2_alpha_degree_CCDF.png
#     E2_partial_effects_heatmap.png
#   tables/
#     E2_replicates.csv
#     E2_summaries_mean_sd.csv
#     E2_partial_effects.csv
# ============================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Config (FAST defaults)
# -----------------------------
cfg <- list(
  out_dir = "colas_PaperE2_knobs_outputs",
  seed = 1,
  
  # graph size + replicates
  n = 2500,
  B = 4,                   # raise to 8–10 for paper-quality
  target_mean_deg = 12,
  
  # calibration budget
  cal_iters = 7,
  cal_B = 2,
  
  # base model
  d = 2,
  rho0 = 20,               # in d=2: candidates per node ~ pi*rho0 ≈ 63
  kernel = "triangular",
  kernel_beta = 1.0,
  
  # popularity marginal (Pareto)
  pareto_wmin = 1.0,
  pareto_alpha0 = 2.8,
  
  # copula / dependence knob (Gaussian copula by default)
  copula_family = "gaussian",
  copula_theta0 = 0.0,
  
  # HT neighbor-search “heavy” threshold (speed hack; exact once combined with brute-force heavy list)
  heavy_q = 0.99,
  
  # knob grids
  alpha_grid = c(2.0, 2.3, 2.8, 3.4, 4.2),            # degree-tail knob (run in HT)
  rho_grid   = c(0.5, 0.75, 1.0, 1.5, 2.0) * 20,      # clustering knob (locality)
  theta_grid = c(0.0, 0.25, 0.5, 0.75, 0.95),         # assortativity knob (one-sided)
  
  # regime per knob (paper E2 wants alpha in HT)
  regime_alpha = "ht",
  regime_geo   = "ht",    # set to "fixed" if you want the old E1.R behavior for the geo knob
  regime_theta = "ht"     # set to "fixed" if you want the old E1.R behavior for the theta knob
)

set.seed(cfg$seed)

# -----------------------------
# 1) Packages
# -----------------------------
need_pkgs <- c("igraph", "copula", "data.table", "ggplot2")
install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    utils::install.packages(missing, repos = "https://cloud.r-project.org")
  }
}
install_if_missing(need_pkgs)

suppressWarnings(suppressMessages({
  library(igraph)
  library(copula)
  library(data.table)
  library(ggplot2)
}))

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2) Helpers
# -----------------------------
kernel_val <- function(dist_scaled, family = c("ball","triangular","epanechnikov"), beta = 1.0) {
  family <- match.arg(family)
  ds <- dist_scaled
  if (family == "ball") return(as.numeric(ds <= 1))
  if (family == "triangular") return(pmax(0, 1 - ds)^beta)
  if (family == "epanechnikov") return(pmax(0, 1 - ds^2)^beta)
  stop("Unknown kernel family")
}

qpareto1 <- function(p, alpha, wmin = 1.0) wmin * ((1 - p)^(-1/alpha))

hill_alpha <- function(x, tail_frac = 0.20, k_min = 15) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n < 2 * k_min) return(NA_real_)
  x <- sort(x, decreasing = TRUE)
  k <- max(k_min, floor(tail_frac * n))
  k <- min(k, n - 1)
  xk1 <- x[k + 1]
  if (!is.finite(xk1) || xk1 <= 0) return(NA_real_)
  gamma_hat <- mean(log(x[1:k] / xk1))
  if (!is.finite(gamma_hat) || gamma_hat <= 0) return(NA_real_)
  1 / gamma_hat
}

copula_bounds <- function(family) {
  family <- tolower(family)
  if (family == "gaussian") return(c(-0.99, 0.99))
  if (family == "frank")    return(c(-20, 20))
  if (family == "clayton")  return(c(0, 10))
  if (family == "gumbel")   return(c(1, 10))
  stop("Unknown copula family")
}

make_copula <- function(family, theta) {
  family <- tolower(family)
  if (family == "gaussian") return(normalCopula(param = theta, dim = 2, dispstr = "un"))
  if (family == "frank")    return(frankCopula(param = theta, dim = 2))
  if (family == "clayton")  return(claytonCopula(param = theta, dim = 2))
  if (family == "gumbel")   return(gumbelCopula(param = theta, dim = 2))
  stop("Unknown copula family")
}

torus_delta <- function(xi, Xj) {
  delta <- abs(sweep(Xj, 2, xi, "-"))
  pmin(delta, 1 - delta)
}

# -----------------------------
# 3) Latents (W, X) with copula seeding
# -----------------------------
sample_latents <- function(n, d,
                           pareto_alpha, pareto_wmin,
                           copula_family, copula_theta) {
  cop <- make_copula(copula_family, copula_theta)
  UV <- rCopula(n, cop)
  U  <- UV[, 1]
  V1 <- UV[, 2]
  
  W <- qpareto1(U, alpha = pareto_alpha, wmin = pareto_wmin)
  
  # positions in [0,1]^d; couple W to first coordinate via V1
  if (d == 1) {
    X <- matrix(V1, ncol = 1)
  } else {
    X <- cbind(V1, matrix(runif(n * (d - 1)), nrow = n, ncol = d - 1))
  }
  list(W = W, X = X)
}

# -----------------------------
# 4) Cell grid (torus) helper
#    Use cell_size = r_n (as in your E0.R spirit)
# -----------------------------
build_cell_grid <- function(X, cell_size) {
  n <- nrow(X); d <- ncol(X)
  ncell <- max(1L, as.integer(ceiling(1 / cell_size)))
  cell_idx <- floor(X / cell_size)
  cell_idx <- cell_idx %% ncell
  mult <- ncell^(0:(d - 1))
  cell_id <- as.integer(cell_idx %*% mult)
  
  cells <- split(seq_len(n), as.character(cell_id))
  cell_env <- list2env(cells, hash = TRUE, parent = emptyenv())
  
  list(ncell = ncell, cell_size = cell_size, cell_idx = cell_idx, mult = mult, cell_env = cell_env)
}

# -----------------------------
# 5) Edge builders: fixed-range vs HT
# -----------------------------
# Fixed-range CoLaS (baseline):
#   p_ij = 1 - exp(-(lambda/rho) * W_i W_j * k(dist/r_n)), with dist <= r_n support
colas_build_fixed <- function(W, X, rho, lambda, kernel, kernel_beta) {
  n <- length(W); d <- ncol(X)
  r_n <- (rho / n)^(1 / d)
  grid <- build_cell_grid(X, cell_size = r_n)
  
  ncell <- grid$ncell
  cell_idx <- grid$cell_idx
  mult <- grid$mult
  cell_env <- grid$cell_env
  rad <- r_n
  
  off <- as.matrix(expand.grid(rep(list(-1:1), d)))
  n_off <- nrow(off)
  
  scale_fac <- lambda / rho
  
  from_list <- vector("list", n)
  to_list   <- vector("list", n)
  
  for (i in seq_len(n - 1)) {
    xi <- X[i, ]
    wi <- W[i]
    ci <- cell_idx[i, ]
    
    nbr_ids <- integer(0)
    for (t in seq_len(n_off)) {
      nbr <- (ci + off[t, ]) %% ncell
      nbr_ids <- c(nbr_ids, as.integer(sum(nbr * mult)))
    }
    nbr_ids <- unique(nbr_ids)
    
    cand_all <- integer(0)
    for (nbr_id in nbr_ids) {
      key <- as.character(nbr_id)
      if (exists(key, envir = cell_env, inherits = FALSE)) {
        cand_all <- c(cand_all, get(key, envir = cell_env, inherits = FALSE))
      }
    }
    
    cand <- cand_all[cand_all > i]
    if (!length(cand)) next
    
    Xj <- X[cand, , drop = FALSE]
    dist <- sqrt(rowSums(torus_delta(xi, Xj)^2))
    
    ok <- dist <= rad
    if (!any(ok)) next
    
    cand2 <- cand[ok]
    dist2 <- dist[ok]
    wj <- W[cand2]
    
    ds <- dist2 / rad
    kval <- kernel_val(ds, family = kernel, beta = kernel_beta)
    
    p <- 1 - exp(-scale_fac * wi * wj * kval)
    sel <- runif(length(p)) < p
    if (!any(sel)) next
    
    to_list[[i]]   <- cand2[sel]
    from_list[[i]] <- rep.int(i, sum(sel))
  }
  
  from <- unlist(from_list, use.names = FALSE)
  to   <- unlist(to_list,   use.names = FALSE)
  
  g <- make_empty_graph(n = n, directed = FALSE)
  if (length(from) > 0) g <- add_edges(g, as.vector(rbind(from, to)))
  simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

# CoLaS-HT (tail-inheriting):
#   p_ij = 1 - exp(-(lambda/rho) * k(dist / (r_n*(W_i W_j)^{1/d})))
# (weights move the range; no W_i W_j multiplicative intensity factor)
colas_build_ht <- function(W, X, rho, lambda, kernel, kernel_beta, heavy_q = 0.99) {
  n <- length(W); d <- ncol(X)
  r_n <- (rho / n)^(1 / d)
  grid <- build_cell_grid(X, cell_size = r_n)
  
  ncell <- grid$ncell
  cell_size <- grid$cell_size
  cell_idx <- grid$cell_idx
  mult <- grid$mult
  cell_env <- grid$cell_env
  
  w_cap <- as.numeric(stats::quantile(W, probs = heavy_q, names = FALSE, type = 1))
  if (!is.finite(w_cap) || w_cap <= 0) w_cap <- max(W[is.finite(W)], na.rm = TRUE)
  
  is_heavy <- (W > w_cap)
  heavy_idx <- which(is_heavy)
  
  # Precompute offsets for h=0..max_h (for non-heavy neighbor search)
  max_h <- as.integer(ceiling((w_cap * w_cap)^(1 / d)))  # since cell_size=r_n cancels
  if (!is.finite(max_h) || max_h < 0) max_h <- 0L
  
  off_list <- vector("list", max_h + 1L)
  for (h in 0:max_h) off_list[[h + 1L]] <- as.matrix(expand.grid(rep(list(-h:h), d)))
  
  scale_fac <- lambda / rho
  
  from_list <- vector("list", n)
  to_list   <- vector("list", n)
  
  for (i in seq_len(n - 1)) {
    xi <- X[i, ]
    wi <- W[i]
    
    # helper to add sampled edges i -> cand (j>i)
    add_edges_from_candidates <- function(cand) {
      if (!length(cand)) return(NULL)
      
      Xj <- X[cand, , drop = FALSE]
      dist <- sqrt(rowSums(torus_delta(xi, Xj)^2))
      
      eff_rad <- r_n * (wi * W[cand])^(1 / d)
      ok <- is.finite(eff_rad) & (eff_rad > 0) & (dist <= eff_rad)
      if (!any(ok)) return(NULL)
      
      cand2 <- cand[ok]
      dist2 <- dist[ok]
      eff2  <- eff_rad[ok]
      
      ds <- dist2 / eff2
      kval <- kernel_val(ds, family = kernel, beta = kernel_beta)
      
      p <- 1 - exp(-scale_fac * kval)
      sel <- runif(length(p)) < p
      if (!any(sel)) return(NULL)
      
      list(to = cand2[sel], from = rep.int(i, sum(sel)))
    }
    
    if (is_heavy[i]) {
      # heavy node: brute force over all j>i
      cand <- (i + 1):n
      out <- add_edges_from_candidates(cand)
      if (!is.null(out)) {
        to_list[[i]]   <- out$to
        from_list[[i]] <- out$from
      }
      next
    }
    
    # (a) non-heavy candidates via cell grid (bounded by w_cap)
    h_i <- as.integer(ceiling((wi * w_cap)^(1 / d) * (r_n / cell_size)))  # r_n/cell_size = 1 here
    if (!is.finite(h_i) || h_i < 0) h_i <- 0L
    if (h_i > max_h) h_i <- max_h
    
    # If h_i is so large it essentially covers the whole torus grid, just brute force.
    if (h_i >= (ncell - 1L)) {
      cand_all <- (i + 1):n
      cand_nh <- cand_all[!is_heavy[cand_all]]
      outA <- add_edges_from_candidates(cand_nh)
      if (!is.null(outA)) {
        to_list[[i]]   <- c(to_list[[i]], outA$to)
        from_list[[i]] <- c(from_list[[i]], outA$from)
      }
    } else {
      off <- off_list[[h_i + 1L]]
      nbr_ids <- integer(0)
      for (t in seq_len(nrow(off))) {
        nbr <- (cell_idx[i, ] + off[t, ]) %% ncell
        nbr_ids <- c(nbr_ids, as.integer(sum(nbr * mult)))
      }
      nbr_ids <- unique(nbr_ids)
      
      cand_all <- integer(0)
      for (nbr_id in nbr_ids) {
        key <- as.character(nbr_id)
        if (exists(key, envir = cell_env, inherits = FALSE)) {
          cand_all <- c(cand_all, get(key, envir = cell_env, inherits = FALSE))
        }
      }
      
      cand_nh <- cand_all[cand_all > i & !is_heavy[cand_all]]
      outA <- add_edges_from_candidates(cand_nh)
      if (!is.null(outA)) {
        to_list[[i]]   <- c(to_list[[i]], outA$to)
        from_list[[i]] <- c(from_list[[i]], outA$from)
      }
    }
    
    # (b) edges to heavy nodes (few): brute force those j>i
    candH <- heavy_idx[heavy_idx > i]
    outB <- add_edges_from_candidates(candH)
    if (!is.null(outB)) {
      to_list[[i]]   <- c(to_list[[i]], outB$to)
      from_list[[i]] <- c(from_list[[i]], outB$from)
    }
  }
  
  from <- unlist(from_list, use.names = FALSE)
  to   <- unlist(to_list,   use.names = FALSE)
  
  g <- make_empty_graph(n = n, directed = FALSE)
  if (length(from) > 0) g <- add_edges(g, as.vector(rbind(from, to)))
  simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

# Unified generator
colas_generate <- function(n, d, rho, lambda,
                           kernel, kernel_beta,
                           pareto_alpha, pareto_wmin,
                           copula_family, copula_theta,
                           regime = c("fixed", "ht"),
                           heavy_q = 0.99) {
  regime <- match.arg(regime)
  lat <- sample_latents(n, d, pareto_alpha, pareto_wmin, copula_family, copula_theta)
  W <- lat$W; X <- lat$X
  
  if (regime == "fixed") {
    colas_build_fixed(W, X, rho, lambda, kernel, kernel_beta)
  } else {
    colas_build_ht(W, X, rho, lambda, kernel, kernel_beta, heavy_q = heavy_q)
  }
}

# -----------------------------
# 6) Stats + CCDF
# -----------------------------
compute_stats_E2 <- function(g) {
  deg <- degree(g)
  mean_deg <- mean(deg)
  
  Cg <- suppressWarnings(transitivity(g, type = "global"))
  if (!is.finite(Cg)) Cg <- 0
  
  r <- suppressWarnings(assortativity_degree(g, directed = FALSE))
  if (!is.finite(r)) r <- 0
  
  alpha_hat <- hill_alpha(deg, tail_frac = 0.20, k_min = 15)
  q99 <- as.numeric(stats::quantile(deg, 0.99, names = FALSE, type = 1))
  
  list(
    mean_deg = mean_deg,
    clustering = Cg,
    assort = r,
    alpha_hat = alpha_hat,
    deg_q99 = q99,
    deg = deg
  )
}

degree_ccdf_dt <- function(deg_vec, label) {
  deg_vec <- deg_vec[is.finite(deg_vec) & deg_vec >= 0]
  if (length(deg_vec) == 0) return(data.table(k = integer(0), ccdf = numeric(0), label = character(0)))
  
  tab <- table(deg_vec)
  k <- as.integer(names(tab))
  ord <- order(k)
  k <- k[ord]
  cnt <- as.numeric(tab)[ord]
  p <- cnt / sum(cnt)
  ccdf <- rev(cumsum(rev(p)))
  
  data.table(k = k, ccdf = ccdf, label = label)[k >= 1 & is.finite(ccdf) & ccdf > 0]
}

# -----------------------------
# 7) Lambda calibration (multiplicative)
# -----------------------------
calibrate_lambda <- function(pars, target_k, lambda_init = 1.0, iters = 7, B = 2,
                             clamp = c(1e-6, 1e6)) {
  lam <- as.numeric(lambda_init)
  
  for (it in seq_len(iters)) {
    k_hat <- numeric(B)
    for (b in seq_len(B)) {
      g <- colas_generate(
        n = pars$n, d = pars$d, rho = pars$rho, lambda = lam,
        kernel = pars$kernel, kernel_beta = pars$kernel_beta,
        pareto_alpha = pars$pareto_alpha, pareto_wmin = pars$pareto_wmin,
        copula_family = pars$copula_family, copula_theta = pars$copula_theta,
        regime = pars$regime, heavy_q = pars$heavy_q
      )
      k_hat[b] <- mean(degree(g))
    }
    
    k_bar <- mean(k_hat, na.rm = TRUE)
    if (!is.finite(k_bar) || k_bar <= 0) break
    
    lam <- lam * (target_k / (k_bar + 1e-12))
    lam <- min(max(lam, clamp[1]), clamp[2])
  }
  
  lam
}

# -----------------------------
# 8) Runner: Paper E2 (effect matrix)
# -----------------------------
run_Paper_E2 <- function(cfg) {
  message("\n=== Paper E2: Separate knobs / partial-effect matrix (synthetic) ===")
  
  # base parameters (knob-specific regimes will override pars$regime)
  base <- list(
    n = cfg$n,
    d = cfg$d,
    rho = cfg$rho0,
    lambda = 1.0,
    kernel = cfg$kernel,
    kernel_beta = cfg$kernel_beta,
    pareto_alpha = cfg$pareto_alpha0,
    pareto_wmin = cfg$pareto_wmin,
    copula_family = cfg$copula_family,
    copula_theta = cfg$copula_theta0,
    heavy_q = cfg$heavy_q,
    regime = "ht"
  )
  
  # Calibrate a baseline lambda per regime used (good initial values)
  used_regimes <- unique(c(cfg$regime_alpha, cfg$regime_geo, cfg$regime_theta))
  base_lambda <- list()
  
  for (rg in used_regimes) {
    tmp <- base
    tmp$regime <- rg
    lam0 <- calibrate_lambda(
      tmp, target_k = cfg$target_mean_deg,
      lambda_init = 1.0, iters = cfg$cal_iters, B = cfg$cal_B
    )
    base_lambda[[rg]] <- lam0
    message("Baseline lambda (regime=", rg, ") to hit mean degree ~", cfg$target_mean_deg, ": ", signif(lam0, 4))
  }
  
  reps_all <- list()
  alpha_deg_pool <- list()
  
  simulate_level <- function(knob_name, knob_value, pars_mod, regime) {
    pars <- modifyList(base, pars_mod)
    pars$regime <- regime
    
    lam_init <- base_lambda[[regime]]
    if (is.null(lam_init) || !is.finite(lam_init) || lam_init <= 0) lam_init <- 1.0
    
    # calibrate lambda at THIS knob level (holding mean degree ~constant)
    lam_hat <- calibrate_lambda(
      pars, target_k = cfg$target_mean_deg,
      lambda_init = lam_init, iters = cfg$cal_iters, B = cfg$cal_B
    )
    pars$lambda <- lam_hat
    
    dt_rep <- vector("list", cfg$B)
    deg_pool <- integer(0)
    
    for (b in seq_len(cfg$B)) {
      g <- colas_generate(
        n = pars$n, d = pars$d, rho = pars$rho, lambda = pars$lambda,
        kernel = pars$kernel, kernel_beta = pars$kernel_beta,
        pareto_alpha = pars$pareto_alpha, pareto_wmin = pars$pareto_wmin,
        copula_family = pars$copula_family, copula_theta = pars$copula_theta,
        regime = pars$regime, heavy_q = pars$heavy_q
      )
      st <- compute_stats_E2(g)
      deg_pool <- c(deg_pool, st$deg)
      
      dt_rep[[b]] <- data.table(
        knob = knob_name,
        level = knob_value,
        rep = b,
        regime = pars$regime,
        lambda_cal = pars$lambda,
        mean_deg = st$mean_deg,
        clustering = st$clustering,
        assort = st$assort,
        alpha_hat = st$alpha_hat,
        deg_q99 = st$deg_q99
      )
    }
    
    list(dt = rbindlist(dt_rep), deg_pool = deg_pool)
  }
  
  # --- degree-tail knob: vary Pareto alpha (paper wants this in CoLaS-HT) ---
  for (a in cfg$alpha_grid) {
    out <- simulate_level("alpha", a, list(pareto_alpha = a), regime = cfg$regime_alpha)
    reps_all[[length(reps_all) + 1]] <- out$dt
    alpha_deg_pool[[as.character(a)]] <- out$deg_pool
  }
  
  # --- clustering knob: vary rho (locality), hold (F_W, theta) fixed ---
  for (rr in cfg$rho_grid) {
    out <- simulate_level("geo", rr, list(rho = rr), regime = cfg$regime_geo)
    reps_all[[length(reps_all) + 1]] <- out$dt
  }
  
  # --- assortativity knob: vary theta, hold (F_W, rho) fixed ---
  bnd <- copula_bounds(cfg$copula_family)
  theta_grid <- pmin(pmax(cfg$theta_grid, bnd[1]), bnd[2])
  for (th in theta_grid) {
    out <- simulate_level("theta", th, list(copula_theta = th), regime = cfg$regime_theta)
    reps_all[[length(reps_all) + 1]] <- out$dt
  }
  
  dt_reps <- rbindlist(reps_all, fill = TRUE)
  fwrite(dt_reps, file.path(cfg$out_dir, "tables", "E2_replicates.csv"))
  
  # --- mean degree flatness plot ---
  mean_deg_summ <- dt_reps[, .(
    mean = mean(mean_deg, na.rm = TRUE),
    sd   = sd(mean_deg,   na.rm = TRUE)
  ), by = .(knob, level)][, level_num := as.numeric(level)]
  
  p_md <- ggplot(mean_deg_summ, aes(x = level_num, y = mean)) +
    geom_hline(yintercept = cfg$target_mean_deg, linetype = "dashed", linewidth = 0.5) +
    geom_line(aes(group = 1), linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, na.rm = TRUE) +
    facet_wrap(~knob, scales = "free_x") +
    theme_minimal() +
    labs(
      title = "Paper E2: Mean degree held ~constant via per-level lambda calibration",
      subtitle = paste0("Dashed line = target mean degree = ", cfg$target_mean_deg),
      x = "knob level",
      y = "mean degree"
    )
  
  ggsave(file.path(cfg$out_dir, "figures", "E2_mean_degree_flatness.png"),
         p_md, width = 9, height = 3.6, dpi = 220)
  
  # --- lambda calibration plot ---
  lambda_summ <- dt_reps[, .(
    mean = mean(lambda_cal, na.rm = TRUE),
    sd   = sd(lambda_cal,   na.rm = TRUE)
  ), by = .(knob, level)][, level_num := as.numeric(level)]
  
  p_lam <- ggplot(lambda_summ, aes(x = level_num, y = mean)) +
    geom_line(aes(group = 1), linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, na.rm = TRUE) +
    facet_wrap(~knob, scales = "free_x") +
    theme_minimal() +
    labs(
      title = "Paper E2: Calibrated lambda per knob level (mean ± sd across replicates)",
      x = "knob level",
      y = "lambda_cal"
    )
  
  ggsave(file.path(cfg$out_dir, "figures", "E2_lambda_calibration.png"),
         p_lam, width = 9, height = 3.6, dpi = 220)
  
  # --- summary mean±sd by knob level ---
  stats <- c("mean_deg", "clustering", "assort", "alpha_hat", "deg_q99")
  summ_long <- rbindlist(lapply(stats, function(s) {
    dt_reps[, .(
      mean = mean(get(s), na.rm = TRUE),
      sd   = sd(get(s),   na.rm = TRUE)
    ), by = .(knob, level)][, stat := s][]
  }))
  summ_long[, level_num := as.numeric(level)]
  fwrite(summ_long, file.path(cfg$out_dir, "tables", "E2_summaries_mean_sd.csv"))
  
  # --- knob sanity plot ---
  p_sanity <- ggplot(summ_long, aes(x = level_num, y = mean)) +
    geom_line(aes(group = 1), linewidth = 0.6, na.rm = TRUE) +
    geom_point(size = 1.8, na.rm = TRUE) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, na.rm = TRUE) +
    facet_grid(knob ~ stat, scales = "free_y") +
    theme_minimal() +
    labs(
      title = "Paper E2: knob sanity check (mean ± sd across replicates)\n(lambda calibrated per level to hold mean degree ~constant)",
      x = "knob level (numeric)",
      y = "statistic"
    )
  
  ggsave(file.path(cfg$out_dir, "figures", "E2_knob_sanity.png"),
         p_sanity, width = 12, height = 6, dpi = 220)
  
  # --- alpha knob: degree CCDF overlay ---
  ccdf_dt <- rbindlist(lapply(names(alpha_deg_pool), function(a) {
    degree_ccdf_dt(alpha_deg_pool[[a]], label = a)
  }), fill = TRUE)
  ccdf_dt[, alpha := label]
  
  p_ccdf <- ggplot(ccdf_dt, aes(x = k, y = ccdf, color = alpha)) +
    geom_line(aes(group = alpha), linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 1.2, na.rm = TRUE) +
    scale_x_log10() + scale_y_log10() +
    theme_minimal() +
    labs(
      title = "Paper E2: Degree CCDF overlay for alpha knob\n(mean degree held ~constant via per-level lambda calibration)",
      x = "degree k (log)",
      y = "CCDF P(D ≥ k) (log)",
      color = "alpha"
    )
  
  ggsave(file.path(cfg$out_dir, "figures", "E2_alpha_degree_CCDF.png"),
         p_ccdf, width = 9, height = 5.2, dpi = 220)
  
  # --- partial-effect matrix: standardized slopes (replicate-level regression) ---
  eff <- rbindlist(lapply(unique(dt_reps$knob), function(kn) {
    dtk <- dt_reps[knob == kn]
    x <- as.numeric(dtk$level)
    sx <- sd(x, na.rm = TRUE)
    
    rbindlist(lapply(stats, function(s) {
      y <- dtk[[s]]
      sy <- sd(y, na.rm = TRUE)
      
      if (!is.finite(sx) || sx == 0 || !is.finite(sy) || sy == 0) {
        return(data.table(knob = kn, stat = s, slope = 0, std_effect = 0))
      }
      
      slope <- tryCatch(coef(lm(y ~ x))[2], error = function(e) NA_real_)
      std_eff <- as.numeric(slope) * sx / sy  # SD(stat) per SD(knob)
      data.table(knob = kn, stat = s, slope = slope, std_effect = std_eff)
    }))
  }))
  
  fwrite(eff, file.path(cfg$out_dir, "tables", "E2_partial_effects.csv"))
  
  p_eff <- ggplot(eff, aes(x = knob, y = stat, fill = std_effect)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", std_effect)), size = 3) +
    scale_fill_gradient2(midpoint = 0) +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = "Paper E2: Partial-effect matrix (standardized slopes; SD(stat) per SD(knob))",
      x = "knob varied",
      y = "statistic",
      fill = "std_effect"
    )
  
  ggsave(file.path(cfg$out_dir, "figures", "E2_partial_effects_heatmap.png"),
         p_eff, width = 10, height = 6, dpi = 240)
  
  message("\nDONE. Paper E2 outputs in: ", normalizePath(cfg$out_dir))
  invisible(list(replicates = dt_reps, summaries = summ_long, effects = eff))
}

# Execute when pasted / sourced:
if (sys.nframe() == 0) {
  run_Paper_E2(cfg)
}
