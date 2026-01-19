# ============================================================
# Paper E3 / Script E2.R: One-graph recovery of theta
# CoLaS + finite-moment regime E[W^6] < Inf
#
# Key paper-alignment changes vs your working E2.R:
#   (1) Enforce finite moments for W:
#       - Default: Pareto alpha_true > 6 (so E[W^6] < Inf)
#       - Optional: truncation W <- pmin(W, tau_n)
#   (2) Remove Hill tail-index estimation from degrees (not meaningful in fixed-range)
#   (3) Focus evaluation on theta_hat (and lambda_hat)
# ============================================================

options(stringsAsFactors = FALSE)

cfg <- list(
  out_dir = "colas_PaperE3_onegraph_outputs",
  seed = 1,
  
  # experiment grid
  n_grid = c(600, 1200, 2500, 5000),
  R = 10,                 # replicates per n (increase to 20–30 for smoother MSE)
  
  # fitting budget
  B_fit = 2,              # sims per theta candidate (5–10 for paper-quality)
  cal_iters = 7,
  cal_B = 2,
  
  # true params
  d = 2,
  rho0 = 20,              # ensure target mean degree is feasible (d=2 => candidates ~ pi*rho0)
  
  kernel = "triangular",
  kernel_beta = 1.0,
  
  # ---- Weight marginal F_W: finite-moment regime ----
  # Pareto alpha>6 ensures E[W^6] < Inf
  pareto_wmin = 1.0,
  alpha_true = 8.0,        # <-- PAPER E3 finite-moment default (CHANGE THIS if you want)
  
  # Optional truncation W -> W ∧ tau_n (paper allows this too)
  truncate_W = FALSE,      # set TRUE if you insist on alpha_true <= 6
  tau0 = 50,               # tau_n = tau0 * n^tau_exp
  tau_exp = 0.0,           # 0.0 => constant truncation (strongest “finite-moment” enforcement)
  
  copula_family = "gaussian",
  theta_true = 0.5,
  
  target_mean_deg = 12,
  
  # theta grid searched by estimator
  theta_grid = seq(-0.95, 0.95, length.out = 9)
)

set.seed(cfg$seed)

# --------- sanity check for paper regime ---------
if (!cfg$truncate_W && cfg$alpha_true <= 6) {
  stop("Paper E3 finite-moment regime needs E[W^6]<Inf. For Pareto this requires alpha_true>6, ",
       "OR set truncate_W=TRUE (W <- W ∧ tau_n).")
}

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

# ---------- helpers ----------
kernel_val <- function(dist_scaled, family = c("ball","triangular","epanechnikov"), beta = 1.0) {
  family <- match.arg(family)
  ds <- dist_scaled
  if (family == "ball") return(as.numeric(ds <= 1))
  if (family == "triangular") return(pmax(0, 1 - ds)^beta)
  if (family == "epanechnikov") return(pmax(0, 1 - ds^2)^beta)
  stop("Unknown kernel family")
}

# Pareto quantile (with clamp for numerical safety)
qpareto1 <- function(p, alpha, wmin = 1.0) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  wmin * ((1 - p)^(-1/alpha))
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

# sample W from copula-uniform U, with optional truncation
draw_W <- function(U, n, pareto_alpha, pareto_wmin, truncate_W, tau0, tau_exp) {
  W <- qpareto1(U, alpha = pareto_alpha, wmin = pareto_wmin)
  if (isTRUE(truncate_W)) {
    tau_n <- tau0 * (n^tau_exp)
    W <- pmin(W, tau_n)
  }
  W
}

# ---------- generator: fixed-range CoLaS ----------
colas_generate <- function(n, d, rho, lambda,
                           kernel, kernel_beta,
                           pareto_alpha, pareto_wmin,
                           copula_family, copula_theta,
                           truncate_W = FALSE, tau0 = 50, tau_exp = 0) {
  
  r_n <- (rho / n)^(1/d)
  if (!is.finite(r_n) || r_n <= 0) stop("Bad r_n; check rho,n,d")
  ncell <- max(1L, as.integer(floor(1 / r_n)))
  cell_size <- 1 / ncell
  rad <- r_n
  
  cop <- make_copula(copula_family, copula_theta)
  UV <- rCopula(n, cop)
  U <- UV[, 1]; V <- UV[, 2]
  
  W <- draw_W(U, n,
              pareto_alpha = pareto_alpha,
              pareto_wmin  = pareto_wmin,
              truncate_W   = truncate_W,
              tau0         = tau0,
              tau_exp      = tau_exp)
  
  # positions: use V as 1st coordinate to couple W-X via copula, rest uniform
  if (d == 1) {
    X <- matrix(V, ncol = 1)
  } else {
    X <- cbind(V, matrix(runif(n * (d - 1)), nrow = n, ncol = d - 1))
  }
  
  cell_idx <- floor(X / cell_size)
  cell_idx[cell_idx >= ncell] <- ncell - 1
  cell_idx[cell_idx < 0] <- 0
  mult <- ncell^(0:(d-1))
  cell_id <- as.integer(cell_idx %*% mult)
  
  cells <- split(seq_len(n), as.character(cell_id))
  cell_env <- list2env(cells, hash = TRUE, parent = emptyenv())
  
  off <- as.matrix(expand.grid(rep(list(-1:1), d)))
  n_off <- nrow(off)
  
  from_list <- vector("list", n)
  to_list   <- vector("list", n)
  
  scale_fac <- lambda / rho
  
  for (i in seq_len(n)) {
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
    if (length(cand_all) == 0) next
    cand <- cand_all[cand_all > i]
    if (length(cand) == 0) next
    
    Xj <- X[cand, , drop = FALSE]
    dj <- torus_delta(xi, Xj)
    dist <- sqrt(rowSums(dj^2))
    
    in_ball <- dist <= rad
    if (!any(in_ball)) next
    
    cand2 <- cand[in_ball]
    dist2 <- dist[in_ball]
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

compute_stats <- function(g) {
  deg <- degree(g)
  mean_deg <- mean(deg)
  Cg <- suppressWarnings(transitivity(g, type = "global"))
  if (!is.finite(Cg)) Cg <- 0
  r <- suppressWarnings(assortativity_degree(g, directed = FALSE))
  if (!is.finite(r)) r <- 0
  list(mean_deg = mean_deg, clustering = Cg, assort = r)
}

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
        truncate_W = pars$truncate_W, tau0 = pars$tau0, tau_exp = pars$tau_exp
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

# ---------- one-graph fit: estimate theta by matching (C_hat, r_hat) ----------
fit_colas_onegraph <- function(g_obs, pars_fixed, theta_grid, B_fit = 2, cal_iters = 7, cal_B = 2) {
  st_obs <- compute_stats(g_obs)
  
  # clamp theta to family bounds
  bnd <- copula_bounds(pars_fixed$copula_family)
  theta_grid <- pmin(pmax(theta_grid, bnd[1]), bnd[2])
  
  dist <- rep(NA_real_, length(theta_grid))
  lam_vec <- rep(NA_real_, length(theta_grid))
  
  for (i in seq_along(theta_grid)) {
    th <- theta_grid[i]
    pars <- pars_fixed
    pars$copula_theta <- th
    
    # match mean degree by calibrating lambda (as in your E2.R)
    lam <- calibrate_lambda(
      pars, target_k = st_obs$mean_deg,
      lambda_init = pars_fixed$lambda,
      iters = cal_iters, B = cal_B
    )
    pars$lambda <- lam
    lam_vec[i] <- lam
    
    # simulate and compare clustering + assort
    Csim <- numeric(B_fit)
    rsim <- numeric(B_fit)
    for (b in seq_len(B_fit)) {
      g <- colas_generate(
        n = pars$n, d = pars$d, rho = pars$rho, lambda = pars$lambda,
        kernel = pars$kernel, kernel_beta = pars$kernel_beta,
        pareto_alpha = pars$pareto_alpha, pareto_wmin = pars$pareto_wmin,
        copula_family = pars$copula_family, copula_theta = pars$copula_theta,
        truncate_W = pars$truncate_W, tau0 = pars$tau0, tau_exp = pars$tau_exp
      )
      st <- compute_stats(g)
      Csim[b] <- st$clustering
      rsim[b] <- st$assort
    }
    
    Cbar <- mean(Csim, na.rm = TRUE)
    rbar <- mean(rsim, na.rm = TRUE)
    
    # scale by observed magnitude (same idea as your script)
    scC <- max(0.05, abs(st_obs$clustering))
    scr <- max(0.05, abs(st_obs$assort))
    
    dist[i] <- ((Cbar - st_obs$clustering) / scC)^2 + ((rbar - st_obs$assort) / scr)^2
  }
  
  ok <- which(is.finite(dist))
  if (length(ok) == 0) {
    return(list(theta_hat = 0, lambda_hat = pars_fixed$lambda))
  }
  best <- ok[which.min(dist[ok])]
  list(theta_hat = theta_grid[best], lambda_hat = lam_vec[best])
}

# ---------- calibrate lambda_true(n) per n ----------
lambda_true_by_n <- list()
for (n in cfg$n_grid) {
  pars_true_n <- list(
    n = n, d = cfg$d, rho = cfg$rho0, lambda = 1.0,
    kernel = cfg$kernel, kernel_beta = cfg$kernel_beta,
    pareto_alpha = cfg$alpha_true, pareto_wmin = cfg$pareto_wmin,
    copula_family = cfg$copula_family, copula_theta = cfg$theta_true,
    truncate_W = cfg$truncate_W, tau0 = cfg$tau0, tau_exp = cfg$tau_exp
  )
  lam_n <- calibrate_lambda(
    pars_true_n, target_k = cfg$target_mean_deg,
    lambda_init = 1.0, iters = cfg$cal_iters, B = cfg$cal_B
  )
  lambda_true_by_n[[as.character(n)]] <- lam_n
  message("lambda_true calibrated at n=", n, " -> ", signif(lam_n, 4))
}

# ---------- run Paper E3 experiment ----------
rows <- list()
iter <- 0L
total <- length(cfg$n_grid) * cfg$R
pb <- utils::txtProgressBar(min = 0, max = total, style = 3)

for (n in cfg$n_grid) {
  lam_true_n <- lambda_true_by_n[[as.character(n)]]
  
  for (r in seq_len(cfg$R)) {
    iter <- iter + 1L
    utils::setTxtProgressBar(pb, iter)
    
    set.seed(cfg$seed + 10000L * n + r)
    
    pars_true <- list(
      n = n, d = cfg$d, rho = cfg$rho0, lambda = lam_true_n,
      kernel = cfg$kernel, kernel_beta = cfg$kernel_beta,
      pareto_alpha = cfg$alpha_true, pareto_wmin = cfg$pareto_wmin,
      copula_family = cfg$copula_family, copula_theta = cfg$theta_true,
      truncate_W = cfg$truncate_W, tau0 = cfg$tau0, tau_exp = cfg$tau_exp
    )
    
    g <- colas_generate(
      n = pars_true$n, d = pars_true$d, rho = pars_true$rho, lambda = pars_true$lambda,
      kernel = pars_true$kernel, kernel_beta = pars_true$kernel_beta,
      pareto_alpha = pars_true$pareto_alpha, pareto_wmin = pars_true$pareto_wmin,
      copula_family = pars_true$copula_family, copula_theta = pars_true$copula_theta,
      truncate_W = pars_true$truncate_W, tau0 = pars_true$tau0, tau_exp = pars_true$tau_exp
    )
    
    pars_fixed <- pars_true
    
    fit <- fit_colas_onegraph(
      g_obs = g,
      pars_fixed = pars_fixed,
      theta_grid = cfg$theta_grid,
      B_fit = cfg$B_fit,
      cal_iters = cfg$cal_iters,
      cal_B = cfg$cal_B
    )
    
    rows[[length(rows) + 1]] <- data.table(
      n = n, rep = r,
      theta_true = cfg$theta_true, theta_hat = fit$theta_hat,
      lambda_true = lam_true_n, lambda_hat = fit$lambda_hat
    )
  }
}
close(pb)

dt <- rbindlist(rows)
fwrite(dt, file.path(cfg$out_dir, "tables", "PaperE3_recovery_raw.csv"))

# squared errors + MSE
dt[, se_theta  := (theta_hat  - theta_true)^2]
dt[, se_lambda := (lambda_hat - lambda_true)^2]
dt[, se_loglambda := (log(lambda_hat + 1e-12) - log(lambda_true + 1e-12))^2]

mse <- dt[, .(
  MSE_theta  = mean(se_theta,  na.rm = TRUE),
  MSE_lambda = mean(se_lambda, na.rm = TRUE),
  MSE_loglambda = mean(se_loglambda, na.rm = TRUE)
), by = n][order(n)]

fwrite(mse, file.path(cfg$out_dir, "tables", "PaperE3_recovery_mse.csv"))

# ---- CI + slope + plots (kept from your E2.R structure) ----
mse_with_ci <- function(dt, se_col) {
  out <- dt[, .(
    MSE = mean(get(se_col), na.rm = TRUE),
    SD  = sd(get(se_col),   na.rm = TRUE),
    N   = sum(is.finite(get(se_col)))
  ), by = n][order(n)]
  
  out[, SE := SD / sqrt(pmax(N, 1L))]
  out[, `:=`(
    CI_low  = pmax(MSE - 1.96 * SE, 1e-12),
    CI_high = MSE + 1.96 * SE
  )]
  out
}

fit_loglog_rate <- function(mse_ci_dt) {
  d <- mse_ci_dt[is.finite(MSE) & MSE > 0]
  fit <- lm(log10(MSE) ~ log10(n), data = d)
  ci  <- suppressMessages(confint(fit))[2, ]
  list(
    slope   = unname(coef(fit)[2]),
    ci_low  = unname(ci[1]),
    ci_high = unname(ci[2]),
    r2      = summary(fit)$r.squared
  )
}

plot_mse_loglog_ci <- function(dt_raw, mse_ci_dt, se_col, title, out_path) {
  raw <- data.table::copy(dt_raw)
  raw[, se_plot := pmax(get(se_col), 1e-12)]
  
  summ <- data.table::copy(mse_ci_dt)
  summ[, `:=`(
    mse_plot = pmax(MSE, 1e-12),
    lo_plot  = pmax(CI_low, 1e-12),
    hi_plot  = pmax(CI_high, 1e-12)
  )]
  
  rate <- fit_loglog_rate(summ)
  subtitle <- sprintf(
    "Fit: log10(MSE)=a+b log10(n);  b=%.2f [%.2f, %.2f],  R^2=%.2f",
    rate$slope, rate$ci_low, rate$ci_high, rate$r2
  )
  
  p <- ggplot() +
    geom_point(data = raw, aes(x = n, y = se_plot),
               alpha = 0.25, size = 1) +
    geom_line(data = summ, aes(x = n, y = mse_plot, group = 1),
              linewidth = 0.8) +
    geom_point(data = summ, aes(x = n, y = mse_plot),
               size = 2.2) +
    geom_errorbar(data = summ, aes(x = n, ymin = lo_plot, ymax = hi_plot),
                  width = 0) +
    scale_x_log10() + scale_y_log10() +
    theme_minimal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "n (log10 scale)",
      y = "MSE (log10 scale)"
    )
  
  ggsave(out_path, p, width = 8, height = 4.8, dpi = 220)
  invisible(rate)
}

mse_theta_ci     <- mse_with_ci(dt, "se_theta")
mse_lambda_ci    <- mse_with_ci(dt, "se_lambda")
mse_loglambda_ci <- mse_with_ci(dt, "se_loglambda")

fwrite(mse_theta_ci,     file.path(cfg$out_dir, "tables", "PaperE3_mse_ci_theta.csv"))
fwrite(mse_lambda_ci,    file.path(cfg$out_dir, "tables", "PaperE3_mse_ci_lambda.csv"))
fwrite(mse_loglambda_ci, file.path(cfg$out_dir, "tables", "PaperE3_mse_ci_loglambda.csv"))

rate_theta     <- fit_loglog_rate(mse_theta_ci)
rate_lambda    <- fit_loglog_rate(mse_lambda_ci)
rate_loglambda <- fit_loglog_rate(mse_loglambda_ci)

rates <- rbindlist(list(
  data.table(param = "theta_hat",     rate_theta),
  data.table(param = "lambda_hat",    rate_lambda),
  data.table(param = "loglambda_hat", rate_loglambda)
))
fwrite(rates, file.path(cfg$out_dir, "tables", "PaperE3_loglog_slopes.csv"))

plot_mse_loglog_ci(
  dt_raw = dt, mse_ci_dt = mse_theta_ci, se_col = "se_theta",
  title = "Paper E3: MSE(theta_hat) vs n (log-log) with uncertainty",
  out_path = file.path(cfg$out_dir, "figures", "PaperE3_MSE_theta_withCI.png")
)

plot_mse_loglog_ci(
  dt_raw = dt, mse_ci_dt = mse_lambda_ci, se_col = "se_lambda",
  title = "Paper E3: MSE(lambda_hat) vs n (log-log) with uncertainty",
  out_path = file.path(cfg$out_dir, "figures", "PaperE3_MSE_lambda_withCI.png")
)

plot_mse_loglog_ci(
  dt_raw = dt, mse_ci_dt = mse_loglambda_ci, se_col = "se_loglambda",
  title = "Paper E3: MSE(log lambda_hat) vs n (log-log) with uncertainty",
  out_path = file.path(cfg$out_dir, "figures", "PaperE3_MSE_loglambda_withCI.png")
)

message("\nDONE. Paper E3 outputs in: ", normalizePath(cfg$out_dir))
