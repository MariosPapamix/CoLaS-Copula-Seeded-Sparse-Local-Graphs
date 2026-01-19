# ============================================================
# E4: Generalization beyond matched moments
# ============================================================

# ----------------------------
# 0) Packages (auto-install)
# ----------------------------
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
ensure_pkg("igraph")
ensure_pkg("Matrix")

# ----------------------------
# 1) Graph hygiene helpers
# ----------------------------
as_simple_undirected <- function(g) {
  g <- igraph::as.undirected(g, mode = "collapse")
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

giant_component <- function(g) {
  g <- as_simple_undirected(g)
  comps <- igraph::components(g)
  keep <- which(comps$membership == which.max(comps$csize))
  igraph::induced_subgraph(g, vids = keep)
}

safe_transitivity <- function(g) {
  x <- tryCatch(igraph::transitivity(g, type = "globalundirected"),
                error = function(e) NA_real_)
  if (!is.finite(x)) 0 else x
}

safe_assort <- function(g) {
  x <- tryCatch(igraph::assortativity_degree(g, directed = FALSE),
                error = function(e) NA_real_)
  if (!is.finite(x)) 0 else x
}

mean_degree <- function(g) mean(igraph::degree(g))
global_trans <- safe_transitivity
edge_assort  <- safe_assort

# ----------------------------
# 2) Tail-index estimator (Hill) for CoLaS-HT tail fit
# ----------------------------
hill_alpha <- function(x, top_frac = 0.05, xmin = NULL) {
  x <- x[x > 0]
  if (length(x) < 50) return(NA_real_)
  x <- sort(x, decreasing = TRUE)
  m <- max(10L, floor(length(x) * top_frac))
  x_top <- x[seq_len(m)]
  if (is.null(xmin)) xmin <- min(x_top)
  if (!is.finite(xmin) || xmin <= 0) return(NA_real_)
  a <- 1 / mean(log(x_top / xmin))
  if (!is.finite(a)) NA_real_ else a
}

calibration_targets <- function(g) {
  c(
    mean_deg = mean_degree(g),
    trans    = global_trans(g),
    assort   = edge_assort(g)
  )
}

# ----------------------------
# 3) Held-out summaries for E4
# ----------------------------
jdd_prob <- function(g) {
  g <- as_simple_undirected(g)
  if (igraph::ecount(g) == 0) return(setNames(numeric(0), character(0)))
  d <- igraph::degree(g)
  el <- igraph::ends(g, igraph::E(g))
  a <- pmin(d[el[, 1]], d[el[, 2]])
  b <- pmax(d[el[, 1]], d[el[, 2]])
  key <- paste(a, b, sep = "|")
  tab <- table(key)
  p <- as.numeric(tab) / sum(tab)
  names(p) <- names(tab)
  p
}

tv_dist_named <- function(p, q) {
  keys <- base::union(names(p), names(q))
  if (length(keys) == 0) return(0)
  p2 <- setNames(numeric(length(keys)), keys)
  q2 <- p2
  p2[names(p)] <- p
  q2[names(q)] <- q
  0.5 * sum(abs(p2 - q2))
}

ck_curve <- function(g, min_count = 10L) {
  g <- as_simple_undirected(g)
  k <- igraph::degree(g)
  c_loc <- tryCatch(
    igraph::transitivity(g, type = "localundirected", isolates = "zero"),
    error = function(e) rep(0, igraph::vcount(g))
  )
  df <- data.frame(k = k, c = c_loc)
  nk <- as.integer(table(df$k))
  kvals <- as.integer(names(table(df$k)))
  keep <- kvals[nk >= min_count]
  df <- df[df$k %in% keep, , drop = FALSE]
  if (nrow(df) == 0) return(data.frame(k = integer(), Ck = numeric(), n = integer()))
  aggC <- tapply(df$c, df$k, function(x) mean(x, na.rm = TRUE))
  aggn <- tapply(df$c, df$k, length)
  data.frame(k = as.integer(names(aggC)), Ck = as.numeric(aggC), n = as.integer(aggn))
}

ck_rmse <- function(ck_obs, ck_sim) {
  if (nrow(ck_obs) == 0 || nrow(ck_sim) == 0) return(NA_real_)
  ks <- ck_obs$k
  sim_map <- setNames(ck_sim$Ck, ck_sim$k)
  sim_vals <- sim_map[as.character(ks)]
  sim_vals[is.na(sim_vals)] <- 0
  sqrt(mean((ck_obs$Ck - sim_vals)^2))
}

triangles_count <- function(g) {
  g <- as_simple_undirected(g)
  if (igraph::ecount(g) == 0) return(0)
  sum(igraph::count_triangles(g)) / 3
}

four_cycles_count <- function(g) {
  g <- as_simple_undirected(g)
  if (igraph::ecount(g) == 0) return(0)
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  C <- A %*% A
  s <- summary(C)
  keep <- s$i < s$j
  x <- s$x[keep]
  sum(x * (x - 1) / 2) / 2
}

core_dist <- function(g) igraph::coreness(as_simple_undirected(g))

ks_dist <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(as.numeric(stats::ks.test(x, y)$statistic))
}

distance_samples <- function(g, n_sources = 150L) {
  gg <- giant_component(g)
  Vn <- igraph::vcount(gg)
  if (Vn < 2) return(numeric())
  src <- sample(igraph::V(gg), size = min(n_sources, Vn))
  D <- igraph::distances(gg, v = src, to = igraph::V(gg))
  d <- as.vector(D)
  d <- d[is.finite(d) & d > 0]
  d
}

spectral_radius <- function(g, iters = 150L, tol = 1e-7) {
  gg <- giant_component(g)
  n <- igraph::vcount(gg)
  if (n < 2) return(NA_real_)
  A <- igraph::as_adjacency_matrix(gg, sparse = TRUE)
  x <- rep(1, n); x <- x / sqrt(sum(x^2))
  lam_old <- 0
  for (t in seq_len(iters)) {
    y <- as.numeric(A %*% x)
    ny <- sqrt(sum(y^2))
    if (!is.finite(ny) || ny == 0) return(NA_real_)
    x <- y / ny
    lam <- sum(x * as.numeric(A %*% x))
    if (is.finite(lam_old) && abs(lam - lam_old) < tol) break
    lam_old <- lam
  }
  if (!is.finite(lam)) NA_real_ else lam
}

rel_err <- function(obs, sim, eps = 1e-12) abs(sim - obs) / (abs(obs) + eps)

heldout_stats <- function(g, ck_min_count = 10L, dist_sources = 150L) {
  g <- as_simple_undirected(g)
  list(
    jdd      = jdd_prob(g),
    ck       = ck_curve(g, min_count = ck_min_count),
    motifs   = c(triangles = triangles_count(g), four_cycles = four_cycles_count(g)),
    cores    = core_dist(g),
    dists    = distance_samples(g, n_sources = dist_sources),
    spec_rad = spectral_radius(g)
  )
}

score_heldout <- function(obs, sim) {
  c(
    JDD_TV       = tv_dist_named(obs$jdd, sim$jdd),
    Ck_RMSE      = ck_rmse(obs$ck, sim$ck),
    Tri_RE       = unname(rel_err(obs$motifs["triangles"],   sim$motifs["triangles"])),
    FourCycle_RE = unname(rel_err(obs$motifs["four_cycles"], sim$motifs["four_cycles"])),
    Core_KS      = ks_dist(obs$cores, sim$cores),
    Dist_KS      = ks_dist(obs$dists, sim$dists),
    SpecRad_RE   = unname(rel_err(obs$spec_rad, sim$spec_rad))
  )
}

# ----------------------------
# 4) Robust copula + toy CoLaS simulators
# ----------------------------
gaussian_copula_uv <- function(n, d = 2L, r = 0.2) {
  r_max <- (1 / sqrt(d)) - 1e-3
  r <- max(min(r, r_max), -r_max)
  
  p <- d + 1L
  Sigma <- diag(p)
  Sigma[1, 2:p] <- r
  Sigma[2:p, 1] <- r
  
  L <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(L)) return(matrix(stats::runif(n * p), n, p))
  
  Z <- matrix(stats::rnorm(n * p), n, p) %*% L
  stats::pnorm(Z)
}

r_pareto <- function(u, alpha = 8, xmin = 1, tau = Inf) {
  w <- xmin * (pmax(1e-12, 1 - u))^(-1 / alpha)
  pmin(w, tau)
}

torus_dist <- function(xi, xj) {
  diff <- abs(xi - xj)
  diff <- pmin(diff, 1 - diff)
  sqrt(sum(diff^2))
}

make_graph_from_edges <- function(n, edges) {
  if (is.null(edges) || nrow(edges) == 0) return(igraph::make_empty_graph(n = n, directed = FALSE))
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  if (igraph::vcount(g) < n) g <- igraph::add_vertices(g, n - igraph::vcount(g))
  as_simple_undirected(g)
}

simulate_colas_fixed <- function(n, d = 2L, eps = 0.07, lambda = 3,
                                 cop_r = 0.2, alpha = 8, tau = 10,
                                 seed = 1L) {
  set.seed(seed)
  UV <- gaussian_copula_uv(n, d = d, r = cop_r)
  U0 <- UV[, 1]
  X  <- UV[, 2:(d + 1), drop = FALSE]
  W  <- r_pareto(U0, alpha = alpha, xmin = 1, tau = tau)
  
  rho_n <- n * (eps^d)
  edges <- matrix(0L, 0, 2)
  
  for (i in 1:(n - 1)) {
    xi <- X[i, ]
    wi <- W[i]
    for (j in (i + 1):n) {
      if (torus_dist(xi, X[j, ]) <= eps) {
        pij <- 1 - exp(-(lambda / rho_n) * wi * W[j])
        if (stats::runif(1) < pij) edges <- rbind(edges, c(i, j))
      }
    }
  }
  make_graph_from_edges(n, edges)
}

simulate_colas_ht <- function(n, d = 2L, eps = 0.07, lambda = 6,
                              cop_r = 0.2, alpha = 8, tau = 10,
                              seed = 1L) {
  set.seed(seed)
  UV <- gaussian_copula_uv(n, d = d, r = cop_r)
  U0 <- UV[, 1]
  X  <- UV[, 2:(d + 1), drop = FALSE]
  W  <- r_pareto(U0, alpha = alpha, xmin = 1, tau = tau)
  
  rho_n <- n * (eps^d)
  p0 <- 1 - exp(-(lambda / rho_n))
  
  edges <- matrix(0L, 0, 2)
  for (i in 1:(n - 1)) {
    xi <- X[i, ]
    wi <- W[i]
    for (j in (i + 1):n) {
      rad <- eps * (wi * W[j])^(1 / d)
      rad <- min(rad, 0.75)
      if (torus_dist(xi, X[j, ]) <= rad) {
        if (stats::runif(1) < p0) edges <- rbind(edges, c(i, j))
      }
    }
  }
  make_graph_from_edges(n, edges)
}

# ----------------------------
# 5) Robust simulation-based moment matching (minimal calibration)
# ----------------------------
fit_by_sim_moments <- function(sim_fn, n, targets,
                               par0, lower, upper,
                               B_fit = 2L, seed = 1L) {
  stopifnot(length(par0) == length(lower), length(lower) == length(upper))
  
  penalty <- function(par, base = 1e12) base + 1e6 * sum(abs(par))
  
  obj <- function(par) {
    if (any(!is.finite(par))) return(penalty(par))
    sims <- vector("list", B_fit)
    for (b in seq_len(B_fit)) {
      g <- tryCatch(sim_fn(par, n, seed + 1000L * b), error = function(e) NULL)
      if (is.null(g)) return(penalty(par))
      sims[[b]] <- as_simple_undirected(g)
    }
    
    sim_stats <- vapply(sims, function(g) {
      c(mean_deg = mean_degree(g), trans = global_trans(g), assort = edge_assort(g))
    }, FUN.VALUE = targets)
    
    avg <- rowMeans(sim_stats)
    if (any(!is.finite(avg))) return(penalty(par))
    
    z <- (avg - targets) / (abs(targets) + 1e-9)
    val <- sum(z^2)
    if (!is.finite(val)) penalty(par) else val
  }
  
  optim(par0, obj, method = "L-BFGS-B", lower = lower, upper = upper)
}

fit_ER_minimal <- function(g_obs) {
  n <- igraph::vcount(g_obs)
  md <- mean_degree(g_obs)
  p <- min(1, max(0, md / (n - 1)))
  list(
    name = "ER",
    simulate_one = function(seed) {
      set.seed(seed)
      igraph::sample_gnp(n, p, directed = FALSE, loops = FALSE)
    }
  )
}

fit_CoLaS_fixed_minimal <- function(g_obs) {
  n <- igraph::vcount(g_obs)
  targets <- calibration_targets(g_obs)
  
  par0  <- c(log(3), log(0.07), 0.1)
  lower <- c(log(0.3), log(0.02), -0.60)
  upper <- c(log(20),  log(0.20),  0.60)
  
  sim_fn <- function(par, n, seed) {
    simulate_colas_fixed(
      n, d = 2L,
      eps = exp(par[2]),
      lambda = exp(par[1]),
      cop_r = par[3],
      alpha = 8, tau = 10,
      seed = seed
    )
  }
  
  opt <- fit_by_sim_moments(sim_fn, n, targets, par0, lower, upper, B_fit = 2L, seed = 1L)
  
  ph <- opt$par
  list(
    name = "CoLaS_fixed",
    par_hat = c(lambda = exp(ph[1]), eps = exp(ph[2]), cop_r = ph[3]),
    simulate_one = function(seed) {
      simulate_colas_fixed(n, d = 2L, eps = exp(ph[2]), lambda = exp(ph[1]),
                           cop_r = ph[3], alpha = 8, tau = 10, seed = seed)
    }
  )
}

fit_CoLaSHT_minimal <- function(g_obs) {
  n <- igraph::vcount(g_obs)
  
  alpha_hat <- hill_alpha(igraph::degree(g_obs), top_frac = 0.05)
  if (!is.finite(alpha_hat)) alpha_hat <- 8
  
  targets <- calibration_targets(g_obs)
  
  par0  <- c(log(6), log(0.07), 0.1)
  lower <- c(log(0.3), log(0.02), -0.60)
  upper <- c(log(30),  log(0.20),  0.60)
  
  sim_fn <- function(par, n, seed) {
    simulate_colas_ht(
      n, d = 2L,
      eps = exp(par[2]),
      lambda = exp(par[1]),
      cop_r = par[3],
      alpha = alpha_hat, tau = 10,
      seed = seed
    )
  }
  
  opt <- fit_by_sim_moments(sim_fn, n, targets, par0, lower, upper, B_fit = 2L, seed = 2L)
  
  ph <- opt$par
  list(
    name = "CoLaS_HT",
    par_hat = c(lambda = exp(ph[1]), eps = exp(ph[2]), cop_r = ph[3], alpha = alpha_hat),
    simulate_one = function(seed) {
      simulate_colas_ht(n, d = 2L, eps = exp(ph[2]), lambda = exp(ph[1]),
                        cop_r = ph[3], alpha = alpha_hat, tau = 10, seed = seed)
    }
  )
}

# ----------------------------
# 6) E4 collection + tables
# ----------------------------
summarize_metrics <- function(metrics_mat, model_name) {
  data.frame(
    model  = model_name,
    metric = colnames(metrics_mat),
    mean   = colMeans(metrics_mat, na.rm = TRUE),
    median = apply(metrics_mat, 2, median, na.rm = TRUE),
    sd     = apply(metrics_mat, 2, sd, na.rm = TRUE),
    row.names = NULL
  )
}

wide_table <- function(e4_table) {
  u_models <- unique(e4_table$model)
  u_metrics <- unique(e4_table$metric)
  out <- data.frame(model = u_models)
  for (m in u_metrics) {
    out[[m]] <- NA_real_
    for (i in seq_along(u_models)) {
      val <- e4_table$mean[e4_table$model == u_models[i] & e4_table$metric == m]
      if (length(val) == 1) out[[m]][i] <- val
    }
  }
  out
}

collect_E4 <- function(g_obs, fitted_models, B = 20L,
                       ck_min_count = 10L, dist_sources = 150L,
                       seed = 999) {
  set.seed(seed)
  g_obs <- as_simple_undirected(g_obs)
  obs <- heldout_stats(g_obs, ck_min_count = ck_min_count, dist_sources = dist_sources)
  
  per_model <- list()
  table_rows <- list()
  
  for (m in fitted_models) {
    cat("\nSimulating + scoring model:", m$name, "\n")
    sims <- vector("list", B)
    metrics_list <- vector("list", B)
    
    for (b in seq_len(B)) {
      g_sim <- as_simple_undirected(m$simulate_one(seed + 1000L * b))
      sim_stats <- heldout_stats(g_sim, ck_min_count = ck_min_count, dist_sources = dist_sources)
      sims[[b]] <- sim_stats
      metrics_list[[b]] <- score_heldout(obs, sim_stats)
    }
    
    metrics_mat <- do.call(rbind, metrics_list)
    per_model[[m$name]] <- list(sims = sims, metrics = metrics_mat)
    table_rows[[m$name]] <- summarize_metrics(metrics_mat, m$name)
  }
  
  e4_table <- do.call(rbind, table_rows)
  e4_wide  <- wide_table(e4_table)
  
  list(obs = obs, per_model = per_model, e4_table = e4_table, e4_wide = e4_wide)
}

# ----------------------------
# 7) Plot helpers -> multi-page PDF
# (ASCII-only labels, no Unicode warnings)
# ----------------------------
mean_named_prob <- function(vec_list) {
  keys <- unique(unlist(lapply(vec_list, names)))
  if (length(keys) == 0) return(setNames(numeric(0), character(0)))
  mat <- vapply(vec_list, function(v) {
    out <- setNames(numeric(length(keys)), keys)
    if (length(v) > 0) out[names(v)] <- v
    out
  }, FUN.VALUE = numeric(length(keys)))
  setNames(rowMeans(mat), keys)
}

ck_summary <- function(ck_list) {
  allk <- sort(unique(unlist(lapply(ck_list, function(df) df$k))))
  if (length(allk) == 0) return(data.frame(k=integer(), mean=numeric(), lo=numeric(), hi=numeric()))
  M <- matrix(NA_real_, nrow = length(allk), ncol = length(ck_list),
              dimnames = list(as.character(allk), NULL))
  for (b in seq_along(ck_list)) {
    df <- ck_list[[b]]
    if (nrow(df) == 0) next
    M[as.character(df$k), b] <- df$Ck
  }
  data.frame(
    k    = allk,
    mean = apply(M, 1, function(x) mean(x, na.rm = TRUE)),
    lo   = apply(M, 1, function(x) stats::quantile(x, 0.10, na.rm = TRUE, names = FALSE)),
    hi   = apply(M, 1, function(x) stats::quantile(x, 0.90, na.rm = TRUE, names = FALSE))
  )
}

plot_ck_panel <- function(obs_ck, sim_ck_list, main = "C(k): observed vs simulated") {
  sim_sum <- ck_summary(sim_ck_list)
  if (nrow(obs_ck) == 0 || nrow(sim_sum) == 0) {
    plot.new(); title(main); text(0.5, 0.5, "C(k) not available"); return(invisible())
  }
  
  ks <- sort(unique(obs_ck$k))
  sim_map_mean <- setNames(sim_sum$mean, sim_sum$k)
  sim_map_lo   <- setNames(sim_sum$lo,   sim_sum$k)
  sim_map_hi   <- setNames(sim_sum$hi,   sim_sum$k)
  
  y_obs <- obs_ck$Ck[match(ks, obs_ck$k)]
  y_mu  <- sim_map_mean[as.character(ks)]
  y_lo  <- sim_map_lo[as.character(ks)]
  y_hi  <- sim_map_hi[as.character(ks)]
  
  ylim <- range(c(y_obs, y_lo, y_hi), finite = TRUE)
  plot(ks, y_obs, type = "p", xlab = "k", ylab = "C(k)", main = main, ylim = ylim)
  lines(ks, y_mu, type = "l")
  lines(ks, y_lo, lty = 2)
  lines(ks, y_hi, lty = 2)
  legend("topright", legend = c("Observed", "Sim mean", "Sim 10-90%"),
         lty = c(NA, 1, 2), pch = c(1, NA, NA), bty = "n")
}

plot_jdd_scatter <- function(obs_jdd, sim_jdd_list, top_m = 250, main = "JDD: obs vs sim mean") {
  sim_mean <- mean_named_prob(sim_jdd_list)
  if (length(obs_jdd) == 0 || length(sim_mean) == 0) {
    plot.new(); title(main); text(0.5, 0.5, "JDD not available"); return(invisible())
  }
  
  ord <- order(obs_jdd, decreasing = TRUE)
  keep <- names(obs_jdd)[ord][seq_len(min(top_m, length(obs_jdd)))]
  x <- obs_jdd[keep]
  y <- sim_mean[keep]
  y[is.na(y)] <- 0
  
  eps <- 1e-10
  plot(log10(x + eps), log10(y + eps),
       xlab = "log10 Obs JDD prob", ylab = "log10 Sim mean JDD prob",
       main = main)
  abline(0, 1, lty = 2)
}

plot_motifs_bar <- function(obs_motifs, sim_motifs_mat, main = "Motifs: observed vs sim mean") {
  mu <- colMeans(sim_motifs_mat, na.rm = TRUE)
  sd <- apply(sim_motifs_mat, 2, sd, na.rm = TRUE)
  
  vals <- rbind(Observed = obs_motifs, SimMean = mu)
  bp <- barplot(vals, beside = TRUE, main = main, ylab = "count")
  
  x1 <- bp[2, ]
  y1 <- mu - sd
  y2 <- mu + sd
  segments(x1, y1, x1, y2)
  segments(x1 - 0.05, y1, x1 + 0.05, y1)
  segments(x1 - 0.05, y2, x1 + 0.05, y2)
  legend("topright", legend = c("Observed", "Sim mean +/- sd"), bty = "n")
}

plot_ecdf_overlay <- function(obs_vec, sim_vec, main = "ECDF overlay", xlab = "value") {
  if (length(obs_vec) < 2 || length(sim_vec) < 2) {
    plot.new(); title(main); text(0.5, 0.5, "Not enough data"); return(invisible())
  }
  e1 <- stats::ecdf(obs_vec)
  e2 <- stats::ecdf(sim_vec)
  xs <- sort(unique(c(obs_vec, sim_vec)))
  plot(xs, e1(xs), type = "l", main = main, xlab = xlab, ylab = "ECDF")
  lines(xs, e2(xs), lty = 2)
  legend("bottomright", legend = c("Observed", "Sim pooled"), lty = c(1,2), bty = "n")
}

plot_specrad_box <- function(obs_sr, sim_sr, main = "Spectral radius across sims") {
  boxplot(sim_sr, main = main, ylab = "spectral radius")
  abline(h = obs_sr, lty = 2)
  legend("topright", legend = c("Observed"), lty = 2, bty = "n")
}

make_figures_pdf <- function(res, pdf_file = "E4_results.pdf") {
  obs <- res$obs
  
  # Use cairo_pdf if available (better font/encoding); fallback to pdf()
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(pdf_file, width = 8.5, height = 6.5)
  } else {
    pdf(pdf_file, width = 8.5, height = 6.5)
  }
  on.exit(dev.off(), add = TRUE)
  
  for (model_name in names(res$per_model)) {
    pm <- res$per_model[[model_name]]
    sims <- pm$sims
    
    sim_ck_list        <- lapply(sims, `[[`, "ck")
    sim_jdd_list       <- lapply(sims, `[[`, "jdd")
    sim_motifs_mat     <- do.call(rbind, lapply(sims, `[[`, "motifs"))
    sim_core_pooled    <- unlist(lapply(sims, `[[`, "cores"))
    sim_dist_pooled    <- unlist(lapply(sims, `[[`, "dists"))
    sim_sr             <- unlist(lapply(sims, `[[`, "spec_rad"))
    
    plot_ck_panel(obs$ck, sim_ck_list, main = paste0(model_name, " - C(k) (held-out)"))
    plot_jdd_scatter(obs$jdd, sim_jdd_list, main = paste0(model_name, " - JDD (held-out)"))
    plot_motifs_bar(obs$motifs, sim_motifs_mat, main = paste0(model_name, " - Motifs (held-out)"))
    plot_ecdf_overlay(obs$cores, sim_core_pooled,
                      main = paste0(model_name, " - Core distribution (held-out)"),
                      xlab = "coreness")
    plot_ecdf_overlay(obs$dists, sim_dist_pooled,
                      main = paste0(model_name, " - Distances on LCC (held-out)"),
                      xlab = "shortest-path length")
    plot_specrad_box(obs$spec_rad, sim_sr, main = paste0(model_name, " - Spectral radius (held-out)"))
  }
  
  plot.new()
  title("E4 mean-error table (wide; lower is better)")
  txt <- capture.output(print(res$e4_wide, row.names = FALSE))
  text(0, 1, paste(txt, collapse = "\n"), adj = c(0,1), family = "mono", cex = 0.8)
  
  invisible(pdf_file)
}

# ============================================================
# 8) Self-contained data: generate an observed graph
# ============================================================
set.seed(123)
n_obs <- 350

g_obs_true <- simulate_colas_ht(
  n = n_obs, d = 2L,
  eps = 0.07, lambda = 6,
  cop_r = 0.45, alpha = 8, tau = 10,
  seed = 123
)
g_obs_true <- as_simple_undirected(g_obs_true)

# Optional IO demo (igraph edgelist)
tmp_file <- tempfile(fileext = ".igraph_edgelist.txt")
igraph::write_graph(g_obs_true, tmp_file, format = "edgelist")
g_obs <- igraph::read_graph(tmp_file, format = "edgelist", directed = FALSE)
g_obs <- as_simple_undirected(g_obs)

cat("\nObserved graph: v=", igraph::vcount(g_obs), " e=", igraph::ecount(g_obs), "\n")
cat("Calibration targets (minimal set):\n")
print(calibration_targets(g_obs))
cat("Tail index (Hill) used for HT tail-fit:\n")
print(hill_alpha(igraph::degree(g_obs), top_frac = 0.05))

# ============================================================
# 9) Fit minimal-calibration models (E4 requirement)
# ============================================================
m_er  <- fit_ER_minimal(g_obs)
m_fix <- fit_CoLaS_fixed_minimal(g_obs)
m_ht  <- fit_CoLaSHT_minimal(g_obs)

cat("\nFitted parameters:\n")
cat("CoLaS_fixed par_hat:\n"); print(m_fix$par_hat)
cat("CoLaS_HT par_hat:\n");    print(m_ht$par_hat)

# ============================================================
# 10) Run E4 + save tables
# ============================================================
res <- collect_E4(
  g_obs,
  fitted_models = list(m_er, m_fix, m_ht),
  B = 20,
  ck_min_count = 10,
  dist_sources = 150,
  seed = 999
)

cat("\nE4 held-out error table (lower is better):\n")
print(res$e4_table)

cat("\nE4 wide table (mean errors):\n")
print(res$e4_wide)

utils::write.csv(res$e4_table, file = "E4_summary_long.csv", row.names = FALSE)
utils::write.csv(res$e4_wide,  file = "E4_summary_wide.csv", row.names = FALSE)
cat("\nWrote: E4_summary_long.csv and E4_summary_wide.csv\n")

ranked <- aggregate(mean ~ model, data = res$e4_table, FUN = function(x) mean(x, na.rm = TRUE))
ranked <- ranked[order(ranked$mean), ]
cat("\nNaive overall ranking (avg of metric means):\n")
print(ranked)

# ============================================================
# 11) Produce figures PDF
# ============================================================
pdf_file <- make_figures_pdf(res, pdf_file = "E4_results.pdf")
cat("\nWrote figures to:", pdf_file, "\n")

# ============================================================
# Notes for YOUR OWN FILE:
# - 2 columns only: read_graph(..., format="ncol")
# - igraph edgelist w/ header "n m": format="edgelist"
# ============================================================
