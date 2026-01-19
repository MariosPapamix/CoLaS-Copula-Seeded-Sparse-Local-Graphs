# ============================================================
# CoLaS E3: Real-data goodness-of-fit across domains
# ------------------------------------------------------------
# Paper-readiness additions:
#   (A) OBS point shown on (assort, clustering) scatter + legend
#   (B) Compact GOF tables (mean ± sd) per dataset + combined:
#       deg_KS, ck_rmse, jdd_L1, core_L1, spec_rel, mean_dist_rel
#   (C) Optional tail-inheriting regime ("tail_inheriting") so
#       heavy degree tails are actually achievable when claimed
# ============================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Config
# -----------------------------
cfg <- list(
  out_dir = "colas_E3_outputs",
  seed = 1,
  
  # evaluation budget (paper runs: increase further if runtime allows)
  B_eval = 80,   # replicate graphs per fitted model (suggest 50–200)
  B_fit  = 8,    # replicate graphs per candidate (rho,theta) (suggest 5–10)
  
  # safety for huge real datasets
  n_max = 5000,  # downsample to at most this many nodes (set Inf to disable)
  
  # ---- CoLaS regime ----
  # "fixed"          : fixed-range kernel (light-tail at extremes is expected)
  # "tail_inheriting": weight-dependent radii (needed if you claim heavy-tail match)
  colas_regime = "fixed",
  
  # CoLaS generator settings
  d = 2,
  kernel = "triangular",
  kernel_beta = 1.0,
  copula_family = "gaussian",   # gaussian / frank / clayton / gumbel
  pareto_wmin = 1.0,
  
  # Hill tail estimate on degrees (proxy for weight tail)
  hill_tail_frac = 0.10,
  hill_k_min = 50,
  
  # lambda calibration budget
  cal_iters = 10,
  cal_B = 2,
  lambda_clamp = c(1e-6, 1e6),
  
  # fitting grids
  rho_grid = c(1, 2, 4, 8, 16, 32),
  theta_grid = seq(-0.9, 0.9, length.out = 9),
  
  # held-out stats settings
  ck_bins = 12,      # bins for C(k)
  jdd_bins = 12,     # bins for joint-degree distribution
  dist_sources = 80, # number of BFS sources for distance sampling
  
  # connectivity enforcement for simulated graphs
  ensure_connected = TRUE,
  connect_tries = 30,
  
  # GOF compact metrics for paper tables
  gof_compact_metrics = c("deg_KS","ck_rmse","jdd_L1","core_L1","spec_rel","mean_dist_rel")
)

set.seed(cfg$seed)

# -----------------------------
# 1) Packages
# -----------------------------
need_pkgs <- c("igraph", "copula", "data.table", "ggplot2", "R.utils")
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
  library(R.utils)
}))

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# 2) igraph compatibility wrappers
# -----------------------------
ig_is_connected <- function(g) {
  if (exists("is_connected", where = asNamespace("igraph"), inherits = FALSE)) {
    return(igraph::is_connected(g))
  }
  if (exists("is.connected", where = asNamespace("igraph"), inherits = FALSE)) {
    return(igraph::is.connected(g))
  }
  comps <- igraph::components(g)
  comps$no == 1
}

ig_distances <- function(g, v, to) {
  if (exists("distances", where = asNamespace("igraph"), inherits = FALSE)) {
    return(igraph::distances(g, v = v, to = to, mode = "all"))
  }
  return(igraph::shortest.paths(g, v = v, to = to, mode = "all"))
}

ig_random_walk <- function(g, start, steps) {
  if (exists("random_walk", where = asNamespace("igraph"), inherits = FALSE)) {
    return(igraph::random_walk(g, start = start, steps = steps, mode = "all"))
  }
  return(igraph::random.walk(g, start = start, steps = steps, mode = "all"))
}

ig_assortativity_degree <- function(g, directed = FALSE) {
  if (exists("assortativity_degree", where = asNamespace("igraph"), inherits = FALSE)) {
    return(igraph::assortativity_degree(g, directed = directed))
  }
  if (exists("assortativity.degree", where = asNamespace("igraph"), inherits = FALSE)) {
    return(igraph::assortativity.degree(g, directed = directed))
  }
  # fallback via edge-endpoint degree correlation
  if (ecount(g) == 0) return(NA_real_)
  deg <- degree(g)
  e <- ends(g, E(g), names = FALSE)
  stats::cor(deg[e[,1]], deg[e[,2]], use = "complete.obs")
}

ig_eigen_centrality_value <- function(g) {
  if (ecount(g) == 0) return(0)
  if (exists("eigen_centrality", where = asNamespace("igraph"), inherits = FALSE)) {
    out <- tryCatch(igraph::eigen_centrality(g, directed = FALSE, scale = FALSE),
                    error = function(e) NULL)
    if (!is.null(out) && is.finite(out$value)) return(as.numeric(out$value))
  }
  if (exists("evcent", where = asNamespace("igraph"), inherits = FALSE)) {
    out <- tryCatch(igraph::evcent(g, directed = FALSE, scale = FALSE),
                    error = function(e) NULL)
    if (!is.null(out) && is.finite(out$value)) return(as.numeric(out$value))
  }
  NA_real_
}

# -----------------------------
# 3) Utility helpers
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

hill_alpha <- function(x, tail_frac = 0.10, k_min = 50) {
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
  if (family == "gaussian") return(c(-0.95, 0.95))
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

# Morton/Z-order mapping: scalar -> ~uniform points on [0,1]^d
scalar_to_unitcube_morton <- function(v, d, bits = 30L) {
  v <- v %% 1
  n <- length(v)
  if (d == 1) return(matrix(v, ncol = 1))
  bits <- as.integer(bits)
  bits_adj <- bits - (bits %% d)
  bpc <- as.integer(bits_adj / d)
  scale <- 2^bits_adj
  I <- as.integer(floor(v * scale))
  X <- matrix(0, nrow = n, ncol = d)
  for (j in 0:(d-1)) {
    coord <- integer(n)
    for (t in 0:(bpc-1)) {
      bitpos <- as.integer(j + d*t)
      bit <- bitwAnd(bitwShiftR(I, bitpos), 1L)
      coord <- coord + bitwShiftL(bit, t)
    }
    X[, j+1] <- coord / 2^bpc
  }
  X
}

torus_delta <- function(xi, Xj) {
  delta <- abs(sweep(Xj, 2, xi, "-"))
  pmin(delta, 1 - delta)
}

# -----------------------------
# 4) CoLaS generator (fixed-range AND tail-inheriting option)
# -----------------------------
OFF_CACHE <- new.env(parent = emptyenv())

get_offsets <- function(reach, d) {
  key <- paste0("d", d, "_r", reach)
  if (exists(key, envir = OFF_CACHE, inherits = FALSE)) {
    return(get(key, envir = OFF_CACHE, inherits = FALSE))
  }
  off <- as.matrix(expand.grid(rep(list(-reach:reach), d)))
  assign(key, off, envir = OFF_CACHE)
  off
}

colas_generate <- function(n,
                           d = 2,
                           rho = 4,
                           lambda = 1.0,
                           kernel = "triangular",
                           kernel_beta = 1.0,
                           pareto_alpha = 2.5,
                           pareto_wmin = 1.0,
                           copula_family = "gaussian",
                           copula_theta = 0.0,
                           regime = c("fixed","tail_inheriting")) {
  regime <- match.arg(regime)
  stopifnot(n >= 2, d >= 1)
  
  rho <- as.numeric(rho)
  lambda <- as.numeric(lambda)
  
  # baseline scale r_n so that n * r_n^d = rho
  base_rad <- (rho / n)^(1/d)
  if (!is.finite(base_rad) || base_rad <= 0) stop("Bad base_rad; check rho,n,d")
  
  # grid hashing on cell_size = base_rad (works well for fixed; acceptable for HT with downsample)
  cell_size <- base_rad
  ncell <- max(1L, as.integer(ceiling(1 / cell_size)))
  cell_size <- 1 / ncell  # snap to exact grid
  
  # copula-seeded weights and positions
  cop <- make_copula(copula_family, copula_theta)
  UV <- rCopula(n, cop)
  U <- UV[, 1]; V <- UV[, 2]
  W <- qpareto1(U, alpha = pareto_alpha, wmin = pareto_wmin)
  X <- scalar_to_unitcube_morton(V, d = d, bits = 30L)
  
  # cell indexing
  cell_idx <- floor(X / cell_size)
  cell_idx[cell_idx >= ncell] <- ncell - 1
  cell_idx[cell_idx < 0] <- 0
  mult <- ncell^(0:(d-1))
  cell_id <- as.integer(cell_idx %*% mult)
  
  cells <- split(seq_len(n), as.character(cell_id))
  cell_env <- list2env(cells, hash = TRUE, parent = emptyenv())
  
  from_list <- vector("list", n)
  to_list   <- vector("list", n)
  
  scale_fac <- (lambda / rho)
  
  # for tail-inheriting: gamma = 1/d (as in the writeup)
  gamma <- 1 / d
  w_max <- max(W)
  
  for (i in seq_len(n)) {
    xi <- X[i, ]
    wi <- W[i]
    ci <- cell_idx[i, ]
    
    # candidate cell reach
    if (regime == "fixed") {
      reach <- 1L
    } else {
      # ensure candidate set covers radius up to base_rad*(wi*w_max)^(1/d)
      rad_i_search <- base_rad * (wi * w_max)^gamma
      reach <- as.integer(ceiling(rad_i_search / cell_size))
      reach <- max(1L, min(reach, ncell))  # clamp
    }
    
    # if reach covers all cells, just consider all nodes > i
    if (reach >= ncell) {
      cand <- (i+1):n
      if (length(cand) == 0) next
      Xj <- X[cand, , drop = FALSE]
      dj <- torus_delta(xi, Xj)
      dist <- sqrt(rowSums(dj^2))
      
      if (regime == "fixed") {
        in_ball <- dist <= base_rad
        if (!any(in_ball)) next
        cand2 <- cand[in_ball]
        dist2 <- dist[in_ball]
        wj <- W[cand2]
        ds <- dist2 / base_rad
        kval <- kernel_val(ds, family = kernel, beta = kernel_beta)
        p <- 1 - exp(-scale_fac * wi * wj * kval)
      } else {
        wj <- W[cand]
        rad_ij <- base_rad * (wi * wj)^gamma
        ds <- dist / rad_ij
        kval <- kernel_val(ds, family = kernel, beta = kernel_beta)
        p <- 1 - exp(-scale_fac * kval)
        cand2 <- cand
      }
      
      sel <- runif(length(p)) < p
      if (!any(sel)) next
      to_list[[i]]   <- cand2[sel]
      from_list[[i]] <- rep.int(i, sum(sel))
      next
    }
    
    # otherwise, gather candidates from nearby cells
    off <- get_offsets(reach, d)
    nbr_ids <- integer(nrow(off))
    for (t in seq_len(nrow(off))) {
      nbr <- (ci + off[t, ]) %% ncell
      nbr_ids[t] <- as.integer(sum(nbr * mult))
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
    cand <- unique(cand_all[cand_all > i])
    if (length(cand) == 0) next
    
    Xj <- X[cand, , drop = FALSE]
    dj <- torus_delta(xi, Xj)
    dist <- sqrt(rowSums(dj^2))
    
    if (regime == "fixed") {
      in_ball <- dist <= base_rad
      if (!any(in_ball)) next
      cand2 <- cand[in_ball]
      dist2 <- dist[in_ball]
      wj <- W[cand2]
      
      ds <- dist2 / base_rad
      kval <- kernel_val(ds, family = kernel, beta = kernel_beta)
      p <- 1 - exp(-scale_fac * wi * wj * kval)
    } else {
      # tail-inheriting: radius depends on (wi*wj)^(1/d); no multiplicative wi*wj in exponent
      wj <- W[cand]
      rad_ij <- base_rad * (wi * wj)^gamma
      ds <- dist / rad_ij
      kval <- kernel_val(ds, family = kernel, beta = kernel_beta)
      p <- 1 - exp(-scale_fac * kval)
      cand2 <- cand
    }
    
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

# -----------------------------
# 5) Lambda calibration to match mean degree
# -----------------------------
calibrate_lambda <- function(pars, target_k,
                             lambda_init = 1.0,
                             iters = 10, B = 2,
                             clamp = c(1e-6, 1e6),
                             regime = c("fixed","tail_inheriting")) {
  regime <- match.arg(regime)
  lam <- as.numeric(lambda_init)
  for (it in seq_len(iters)) {
    k_hat <- numeric(B)
    for (b in seq_len(B)) {
      g <- colas_generate(
        n = pars$n, d = pars$d, rho = pars$rho, lambda = lam,
        kernel = pars$kernel, kernel_beta = pars$kernel_beta,
        pareto_alpha = pars$pareto_alpha, pareto_wmin = pars$pareto_wmin,
        copula_family = pars$copula_family, copula_theta = pars$copula_theta,
        regime = regime
      )
      k_hat[b] <- mean(degree(g))
    }
    k_bar <- mean(k_hat, na.rm = TRUE)
    if (!is.finite(k_bar) || k_bar <= 0) break
    
    lam <- lam * (target_k / (k_bar + 1e-12))
    lam <- min(max(lam, clamp[1]), clamp[2])
    
    if (abs(k_bar - target_k) / max(1e-12, target_k) < 0.02) break
  }
  lam
}

# -----------------------------
# 6) Robust SNAP I/O (data.table comment compat)
# -----------------------------
download_to <- function(url, dest_dir) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  fname <- basename(url)
  dest <- file.path(dest_dir, fname)
  if (!file.exists(dest)) {
    message("Downloading: ", url)
    utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
  }
  dest
}

maybe_gunzip <- function(path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(path)
  out <- sub("\\.gz$", "", path, ignore.case = TRUE)
  if (!file.exists(out)) {
    message("Decompressing: ", basename(path))
    R.utils::gunzip(path, destname = out, remove = FALSE, overwrite = TRUE)
  }
  out
}

fread_two_col_skip_comments <- function(path) {
  f <- data.table::fread
  formals_names <- names(formals(f))
  
  if ("comment.char" %in% formals_names) {
    return(f(path,
             header = FALSE,
             select = c(1, 2),
             col.names = c("from", "to"),
             comment.char = "#",
             showProgress = FALSE))
  }
  
  if (nzchar(Sys.which("grep"))) {
    cmd <- paste("grep -v '^#' ", shQuote(path))
    return(f(cmd = cmd,
             header = FALSE,
             select = c(1, 2),
             col.names = c("from", "to"),
             showProgress = FALSE))
  }
  
  tmp <- tempfile(fileext = ".txt")
  con <- file(path, open = "r")
  out <- file(tmp, open = "w")
  on.exit({
    try(close(con), silent = TRUE)
    try(close(out), silent = TRUE)
  }, add = TRUE)
  
  repeat {
    lines <- readLines(con, n = 500000)
    if (length(lines) == 0) break
    lines <- lines[!grepl("^\\s*#", lines)]
    if (length(lines) > 0) writeLines(lines, out, sep = "\n")
  }
  close(out); close(con)
  
  f(tmp,
    header = FALSE,
    select = c(1, 2),
    col.names = c("from", "to"),
    showProgress = FALSE)
}

read_snap_edgelist <- function(path, directed = FALSE) {
  dt <- fread_two_col_skip_comments(path)
  g <- graph_from_data_frame(dt, directed = directed)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  if (directed) g <- as.undirected(g, mode = "collapse")
  g
}

largest_cc <- function(g) {
  if (vcount(g) == 0) return(g)
  comps <- components(g)
  giant <- which.max(comps$csize)
  vids <- which(comps$membership == giant)
  induced_subgraph(g, vids)
}

subsample_graph <- function(g, n_max, seed = 1) {
  if (!is.finite(n_max)) return(g)
  n <- vcount(g)
  if (n <= n_max) return(g)
  
  set.seed(seed)
  start <- sample(V(g), 1)
  walk <- ig_random_walk(g, start = start, steps = max(10L, 10L * n_max))
  nodes <- unique(as.integer(walk))
  if (length(nodes) < n_max) {
    nodes <- unique(c(nodes, sample(seq_len(n), n_max - length(nodes))))
  }
  nodes <- nodes[seq_len(min(n_max, length(nodes)))]
  g2 <- induced_subgraph(g, nodes)
  largest_cc(simplify(g2, remove.multiple = TRUE, remove.loops = TRUE))
}

# -----------------------------
# 7) Feature extraction (targets + held-out)
# -----------------------------
make_degree_breaks <- function(deg, nbins = 12) {
  k <- deg[is.finite(deg) & deg > 0]
  if (length(k) == 0) return(c(-0.5, 0.5, 1.5))
  kmin <- max(1, min(k))
  kmax <- max(k)
  raw <- exp(seq(log(kmin), log(kmax + 1), length.out = nbins + 1)) - 1
  br <- sort(unique(floor(raw)))
  br <- unique(c(min(k), br, max(k) + 1))
  br <- sort(unique(br))
  br <- unique(pmax(br, 0))
  if (length(br) < 3) br <- unique(c(min(k), floor((min(k)+max(k))/2), max(k)+1))
  br
}

ks_dist <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  vals <- sort(unique(c(x, y)))
  Fx <- ecdf(x)(vals)
  Fy <- ecdf(y)(vals)
  max(abs(Fx - Fy))
}

ck_curve <- function(g, breaks) {
  deg <- degree(g)
  cl  <- transitivity(g, type = "local", isolates = "zero")
  bin <- cut(deg, breaks = breaks, include.lowest = TRUE, right = FALSE)
  dt <- data.table(bin = bin, deg = deg, cl = cl)
  dt <- dt[!is.na(bin)]
  dt[, .(k_mid = mean(deg), Ck = mean(cl, na.rm = TRUE)), by = bin][order(k_mid)]
}

ck_dist <- function(ck_obs, ck_sim) {
  if (nrow(ck_obs) < 2 || nrow(ck_sim) < 2) return(NA_real_)
  x <- ck_obs$k_mid
  y_obs <- ck_obs$Ck
  y_sim <- approx(ck_sim$k_mid, ck_sim$Ck, xout = x, rule = 2)$y
  sqrt(mean((y_sim - y_obs)^2, na.rm = TRUE))
}

jdd_matrix <- function(g, breaks) {
  deg <- degree(g)
  e <- ends(g, E(g), names = FALSE)
  nb <- length(breaks) - 1
  if (nrow(e) == 0) return(matrix(0, nb, nb))
  
  k1 <- deg[e[,1]]
  k2 <- deg[e[,2]]
  b1 <- cut(k1, breaks = breaks, include.lowest = TRUE, right = FALSE)
  b2 <- cut(k2, breaks = breaks, include.lowest = TRUE, right = FALSE)
  
  lv <- levels(b1)
  b1 <- factor(b1, levels = lv)
  b2 <- factor(b2, levels = lv)
  
  mat <- matrix(0, nb, nb, dimnames = list(lv, lv))
  for (i in seq_along(b1)) {
    ii <- as.integer(b1[i])
    jj <- as.integer(b2[i])
    if (is.na(ii) || is.na(jj)) next
    mat[ii, jj] <- mat[ii, jj] + 1
    mat[jj, ii] <- mat[jj, ii] + 1
  }
  s <- sum(mat)
  if (s > 0) mat / s else mat
}

jdd_dist <- function(j_obs, j_sim) {
  if (any(dim(j_obs) != dim(j_sim))) return(NA_real_)
  sum(abs(j_obs - j_sim))
}

core_hist <- function(g) {
  c <- coreness(g)
  tab <- table(c)
  data.table(core = as.integer(names(tab)), p = as.numeric(tab) / sum(tab))
}

core_dist <- function(h_obs, h_sim) {
  dt <- merge(h_obs, h_sim, by = "core", all = TRUE, suffixes = c("_obs","_sim"))
  dt[is.na(p_obs), p_obs := 0]
  dt[is.na(p_sim), p_sim := 0]
  sum(abs(dt$p_sim - dt$p_obs))
}

mean_dist_sample <- function(g, n_src = 80, seed = 1) {
  n <- vcount(g)
  if (n <= 2) return(NA_real_)
  set.seed(seed)
  src <- sample(seq_len(n), min(n_src, n))
  dmat <- ig_distances(g, v = src, to = V(g))
  vals <- as.numeric(dmat)
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

safe_scalar <- function(x, fallback = 0) {
  if (!is.finite(x)) fallback else x
}

compute_features <- function(g, breaks_deg = NULL, seed = 1) {
  deg <- degree(g)
  if (is.null(breaks_deg)) breaks_deg <- make_degree_breaks(deg, nbins = cfg$ck_bins)
  
  cl_global <- safe_scalar(suppressWarnings(transitivity(g, type = "global")), 0)
  assort    <- safe_scalar(suppressWarnings(ig_assortativity_degree(g, directed = FALSE)), 0)
  
  # Triangle count consistent with "global transitivity" definition: C = 3T/W
  Wedges <- sum(choose(deg, 2))
  tri <- if (is.finite(cl_global) && is.finite(Wedges) && Wedges > 0) (cl_global * Wedges / 3) else NA_real_
  
  list(
    scalars = list(
      n = vcount(g),
      m = ecount(g),
      mean_deg = mean(deg),
      clustering = cl_global,
      assort = assort,
      alpha_hat = hill_alpha(deg, tail_frac = cfg$hill_tail_frac, k_min = cfg$hill_k_min),
      tri = safe_scalar(tri, 0),
      spec = safe_scalar(ig_eigen_centrality_value(g), NA_real_),
      mean_dist = mean_dist_sample(g, n_src = cfg$dist_sources, seed = seed)
    ),
    deg = deg,
    breaks = breaks_deg,
    ck = ck_curve(g, breaks_deg),
    jdd = jdd_matrix(g, breaks_deg),
    core = core_hist(g)
  )
}

# -----------------------------
# 8) Baseline generators
# -----------------------------
gen_config_model_connected <- function(deg_seq, max_tries = 30) {
  best_g <- NULL
  best_cc <- -Inf
  
  for (t in seq_len(max_tries)) {
    g <- tryCatch(sample_degseq(deg_seq, method = "vl"),
                  error = function(e) NULL)
    if (is.null(g)) next
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    
    comps <- components(g)
    cc <- max(comps$csize)
    if (cc > best_cc) {
      best_cc <- cc
      best_g <- g
    }
    if (ig_is_connected(g)) return(g)
  }
  
  if (is.null(best_g)) {
    g <- sample_degseq(deg_seq, method = "simple")
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    return(largest_cc(g))
  }
  largest_cc(best_g)
}

gen_rgg <- function(n, target_mean_deg, torus = TRUE) {
  r <- sqrt(max(target_mean_deg, 1e-9) / (max(1, n - 1) * pi))
  if (exists("sample_grg", where = asNamespace("igraph"), inherits = FALSE)) {
    g <- igraph::sample_grg(n = n, radius = r, torus = torus)
    return(simplify(g, remove.multiple = TRUE, remove.loops = TRUE))
  }
  # brute force fallback
  X <- matrix(runif(2 * n), ncol = 2)
  A <- matrix(FALSE, n, n)
  for (i in 1:(n-1)) {
    dx <- abs(X[i,1] - X[(i+1):n, 1])
    dy <- abs(X[i,2] - X[(i+1):n, 2])
    if (torus) {
      dx <- pmin(dx, 1 - dx)
      dy <- pmin(dy, 1 - dy)
    }
    dist <- sqrt(dx^2 + dy^2)
    sel <- dist <= r
    if (any(sel)) {
      jj <- (i+1):n
      A[i, jj[sel]] <- TRUE
    }
  }
  idx <- which(A, arr.ind = TRUE)
  g <- make_empty_graph(n = n, directed = FALSE)
  if (nrow(idx) > 0) g <- add_edges(g, as.vector(t(idx)))
  simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

ensure_connected_graph <- function(gen_fun, tries = 20) {
  g_last <- NULL
  for (t in seq_len(tries)) {
    g <- gen_fun()
    g <- simplify(as.undirected(g, mode = "collapse"),
                  remove.multiple = TRUE, remove.loops = TRUE)
    g_last <- g
    if (ig_is_connected(g)) return(g)
  }
  largest_cc(g_last)
}

# -----------------------------
# 9) Fit CoLaS (and Geo+tail baseline)
# -----------------------------
fit_colas <- function(obs_feat) {
  obs <- obs_feat$scalars
  n <- obs$n
  
  alpha_fit <- obs$alpha_hat
  if (!is.finite(alpha_fit) || alpha_fit <= 1.5) alpha_fit <- 2.5
  
  bnd <- copula_bounds(cfg$copula_family)
  theta_grid <- pmin(pmax(cfg$theta_grid, bnd[1]), bnd[2])
  
  best <- list(score = Inf, pars = NULL, fit_mom = NULL)
  
  for (rho in cfg$rho_grid) {
    for (theta in theta_grid) {
      pars <- list(
        n = n, d = cfg$d, rho = rho, lambda = 1.0,
        kernel = cfg$kernel, kernel_beta = cfg$kernel_beta,
        pareto_alpha = alpha_fit, pareto_wmin = cfg$pareto_wmin,
        copula_family = cfg$copula_family, copula_theta = theta
      )
      
      pars$lambda <- calibrate_lambda(
        pars, target_k = obs$mean_deg,
        lambda_init = 1.0,
        iters = cfg$cal_iters, B = cfg$cal_B,
        clamp = cfg$lambda_clamp,
        regime = cfg$colas_regime
      )
      
      Csim <- numeric(cfg$B_fit)
      rsim <- numeric(cfg$B_fit)
      for (b in seq_len(cfg$B_fit)) {
        g <- colas_generate(
          n = pars$n, d = pars$d, rho = pars$rho, lambda = pars$lambda,
          kernel = pars$kernel, kernel_beta = pars$kernel_beta,
          pareto_alpha = pars$pareto_alpha, pareto_wmin = pars$pareto_wmin,
          copula_family = pars$copula_family, copula_theta = pars$copula_theta,
          regime = cfg$colas_regime
        )
        Csim[b] <- safe_scalar(suppressWarnings(transitivity(g, type = "global")), 0)
        rsim[b] <- safe_scalar(suppressWarnings(ig_assortativity_degree(g, directed = FALSE)), 0)
      }
      
      Cbar <- mean(Csim, na.rm = TRUE)
      rbar <- mean(rsim, na.rm = TRUE)
      
      scC <- max(0.02, abs(obs$clustering))
      scr <- max(0.02, abs(obs$assort))
      score <- ((Cbar - obs$clustering) / scC)^2 + ((rbar - obs$assort) / scr)^2
      
      if (is.finite(score) && score < best$score) {
        best <- list(
          score = score,
          pars = pars,
          fit_mom = list(Cbar = Cbar, rbar = rbar, alpha_fit = alpha_fit)
        )
      }
    }
  }
  
  if (is.null(best$pars)) {
    rho <- max(cfg$rho_grid)
    theta <- 0
    pars <- list(
      n = n, d = cfg$d, rho = rho, lambda = 1.0,
      kernel = cfg$kernel, kernel_beta = cfg$kernel_beta,
      pareto_alpha = alpha_fit, pareto_wmin = cfg$pareto_wmin,
      copula_family = cfg$copula_family, copula_theta = theta
    )
    pars$lambda <- calibrate_lambda(pars, target_k = obs$mean_deg,
                                    lambda_init = 1.0, iters = cfg$cal_iters, B = cfg$cal_B,
                                    clamp = cfg$lambda_clamp,
                                    regime = cfg$colas_regime)
    best <- list(score = Inf, pars = pars,
                 fit_mom = list(Cbar = NA, rbar = NA, alpha_fit = alpha_fit))
  }
  
  best
}

fit_geotail_indep <- function(obs_feat, alpha_fit, rho_grid = NULL) {
  obs <- obs_feat$scalars
  n <- obs$n
  rho_grid <- rho_grid %||% cfg$rho_grid
  
  best <- list(score = Inf, pars = NULL)
  
  for (rho in rho_grid) {
    pars <- list(
      n = n, d = cfg$d, rho = rho, lambda = 1.0,
      kernel = cfg$kernel, kernel_beta = cfg$kernel_beta,
      pareto_alpha = alpha_fit, pareto_wmin = cfg$pareto_wmin,
      copula_family = cfg$copula_family, copula_theta = 0.0
    )
    pars$lambda <- calibrate_lambda(pars, target_k = obs$mean_deg,
                                    lambda_init = 1.0, iters = cfg$cal_iters, B = cfg$cal_B,
                                    clamp = cfg$lambda_clamp,
                                    regime = cfg$colas_regime)
    
    Csim <- numeric(cfg$B_fit)
    for (b in seq_len(cfg$B_fit)) {
      g <- colas_generate(
        n = pars$n, d = pars$d, rho = pars$rho, lambda = pars$lambda,
        kernel = pars$kernel, kernel_beta = pars$kernel_beta,
        pareto_alpha = pars$pareto_alpha, pareto_wmin = pars$pareto_wmin,
        copula_family = pars$copula_family, copula_theta = pars$copula_theta,
        regime = cfg$colas_regime
      )
      Csim[b] <- safe_scalar(suppressWarnings(transitivity(g, type = "global")), 0)
    }
    Cbar <- mean(Csim, na.rm = TRUE)
    scC <- max(0.02, abs(obs$clustering))
    score <- ((Cbar - obs$clustering) / scC)^2
    if (is.finite(score) && score < best$score) best <- list(score = score, pars = pars)
  }
  
  if (is.null(best$pars)) {
    rho <- max(rho_grid)
    pars <- list(
      n = n, d = cfg$d, rho = rho, lambda = 1.0,
      kernel = cfg$kernel, kernel_beta = cfg$kernel_beta,
      pareto_alpha = alpha_fit, pareto_wmin = cfg$pareto_wmin,
      copula_family = cfg$copula_family, copula_theta = 0.0
    )
    pars$lambda <- calibrate_lambda(pars, target_k = obs$mean_deg,
                                    lambda_init = 1.0, iters = cfg$cal_iters, B = cfg$cal_B,
                                    clamp = cfg$lambda_clamp,
                                    regime = cfg$colas_regime)
    best <- list(score = Inf, pars = pars)
  }
  
  best
}

# -----------------------------
# 10) GOF distances
# -----------------------------
feature_to_row <- function(feat, model, rep) {
  s <- feat$scalars
  data.table(
    model = model, rep = rep,
    n = s$n, m = s$m,
    mean_deg = s$mean_deg,
    clustering = s$clustering,
    assort = s$assort,
    alpha_hat = s$alpha_hat,
    tri = s$tri,
    spec = s$spec,
    mean_dist = s$mean_dist
  )
}

gof_distances <- function(obs_feat, sim_feat) {
  obs <- obs_feat$scalars
  sim <- sim_feat$scalars
  
  data.table(
    err_mean_deg = abs(sim$mean_deg - obs$mean_deg),
    err_clustering = abs(sim$clustering - obs$clustering),
    err_assort = abs(sim$assort - obs$assort),
    err_alpha = abs(sim$alpha_hat - obs$alpha_hat),
    deg_KS = ks_dist(sim_feat$deg, obs_feat$deg),
    ck_rmse = ck_dist(obs_feat$ck, sim_feat$ck),
    jdd_L1 = jdd_dist(obs_feat$jdd, sim_feat$jdd),
    core_L1 = core_dist(obs_feat$core, sim_feat$core),
    spec_rel = ifelse(is.finite(obs$spec) && obs$spec > 0,
                      abs(sim$spec - obs$spec) / obs$spec, NA_real_),
    mean_dist_rel = ifelse(is.finite(obs$mean_dist) && obs$mean_dist > 0,
                           abs(sim$mean_dist - obs$mean_dist) / obs$mean_dist, NA_real_)
  )
}

# compact GOF summary (paper table)
summarize_compact_gof <- function(dist_dt, metrics = cfg$gof_compact_metrics) {
  # numeric means/sds
  num <- dist_dt[, {
    m <- sapply(.SD, mean, na.rm = TRUE)
    s <- sapply(.SD, sd, na.rm = TRUE)
    out <- c(as.list(setNames(as.numeric(m), paste0(names(m), "_mean"))),
             as.list(setNames(as.numeric(s), paste0(names(s), "_sd"))))
    out
  }, by = model, .SDcols = metrics]
  
  # formatted mean ± sd
  fmt <- num[, {
    out <- list(model = model)
    for (mm in metrics) {
      mu <- get(paste0(mm, "_mean"))
      sd <- get(paste0(mm, "_sd"))
      out[[mm]] <- sprintf("%.4g \u00B1 %.4g", mu, sd)
    }
    out
  }, by = model]
  
  list(numeric = num, formatted = fmt)
}

# -----------------------------
# 11) Degree CCDF overlay
# -----------------------------
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

plot_degree_ccdf_overlay <- function(obs_deg, sim_deg_list, out_png, title) {
  dt <- degree_ccdf_dt(obs_deg, "OBS")
  for (nm in names(sim_deg_list)) {
    dt <- rbind(dt, degree_ccdf_dt(sim_deg_list[[nm]], nm), fill = TRUE)
  }
  p <- ggplot(dt, aes(x = k, y = ccdf, color = label)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    scale_x_log10() + scale_y_log10() +
    theme_minimal() +
    labs(title = title, x = "degree k (log)", y = "CCDF P(D >= k) (log)", color = "")
  ggsave(out_png, p, width = 9, height = 5.2, dpi = 200)
}

# -----------------------------
# 12) Dataset list (Representative datasets)
# -----------------------------
datasets <- list(
  list(name = "facebook_ego", source = "snap",
       url = "https://snap.stanford.edu/data/facebook_combined.txt.gz", directed = FALSE, enabled = TRUE),
  list(name = "livejournal", source = "snap",
       url = "https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz", directed = FALSE, enabled = TRUE),
  
  list(name = "arxiv_AstroPh", source = "snap",
       url = "https://snap.stanford.edu/data/ca-AstroPh.txt.gz", directed = FALSE, enabled = TRUE),
  list(name = "arxiv_CondMat", source = "snap",
       url = "https://snap.stanford.edu/data/ca-CondMat.txt.gz", directed = FALSE, enabled = TRUE),
  
  list(name = "as_caida20071105", source = "snap",
       url = "https://snap.stanford.edu/data/as-caida20071105.txt.gz", directed = FALSE, enabled = TRUE),
  
  list(name = "amazon_copurchase", source = "snap",
       url = "https://snap.stanford.edu/data/com-Amazon.ungraph.txt.gz", directed = FALSE, enabled = TRUE),
  
  list(name = "biogrid_ppi", source = "local",
       file = "BIOGRID-ALL.tab2.txt", directed = FALSE, enabled = file.exists("BIOGRID-ALL.tab2.txt"))
)

datasets <- Filter(function(x) isTRUE(x$enabled), datasets)

# -----------------------------
# 13) Main E3 runner
# -----------------------------
run_E3_one_dataset <- function(ds) {
  message("\n=== E3: ", ds$name, " (regime=", cfg$colas_regime, ") ===")
  
  ds_dir <- file.path(cfg$out_dir, ds$name)
  dir.create(file.path(ds_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(ds_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
  
  # ---- load graph ----
  if (ds$source == "snap") {
    raw <- download_to(ds$url, dest_dir = file.path(cfg$out_dir, "_downloads"))
    path <- maybe_gunzip(raw)
    g0 <- read_snap_edgelist(path, directed = ds$directed %||% FALSE)
  } else {
    g0 <- read_snap_edgelist(ds$file, directed = ds$directed %||% FALSE)
  }
  
  # preprocess: undirected + simplify + LCC + optional downsample
  g0 <- simplify(as.undirected(g0, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
  g0 <- largest_cc(g0)
  g0 <- subsample_graph(g0, n_max = cfg$n_max, seed = cfg$seed)
  g0 <- simplify(g0, remove.multiple = TRUE, remove.loops = TRUE)
  
  # observed features
  obs_feat <- compute_features(g0, seed = cfg$seed)
  data.table::fwrite(as.data.table(obs_feat$scalars), file.path(ds_dir, "tables", "obs_scalars.csv"))
  
  # ---- fit models ----
  fitC <- fit_colas(obs_feat)
  pars_colas <- fitC$pars
  alpha_fit <- fitC$fit_mom$alpha_fit %||% pars_colas$pareto_alpha
  
  fitGT <- fit_geotail_indep(obs_feat, alpha_fit = alpha_fit)
  deg_obs <- obs_feat$deg
  n <- obs_feat$scalars$n
  
  fits_dt <- rbindlist(list(
    data.table(model = "CoLaS",
               regime = cfg$colas_regime,
               rho = pars_colas$rho, theta = pars_colas$copula_theta,
               lambda = pars_colas$lambda, alpha = pars_colas$pareto_alpha),
    data.table(model = "GeoTail_indep",
               regime = cfg$colas_regime,
               rho = fitGT$pars$rho, theta = 0.0,
               lambda = fitGT$pars$lambda, alpha = fitGT$pars$pareto_alpha),
    data.table(model = "RGG",
               regime = NA_character_, rho = NA_real_, theta = NA_real_,
               lambda = NA_real_, alpha = NA_real_),
    data.table(model = "Config",
               regime = NA_character_, rho = NA_real_, theta = NA_real_,
               lambda = NA_real_, alpha = NA_real_)
  ), fill = TRUE)
  data.table::fwrite(fits_dt, file.path(ds_dir, "tables", "fits.csv"))
  
  # ---- generators ----
  gen_colas <- function() colas_generate(
    n = n, d = pars_colas$d, rho = pars_colas$rho, lambda = pars_colas$lambda,
    kernel = pars_colas$kernel, kernel_beta = pars_colas$kernel_beta,
    pareto_alpha = pars_colas$pareto_alpha, pareto_wmin = pars_colas$pareto_wmin,
    copula_family = pars_colas$copula_family, copula_theta = pars_colas$copula_theta,
    regime = cfg$colas_regime
  )
  
  gen_geotail <- function() colas_generate(
    n = n, d = fitGT$pars$d, rho = fitGT$pars$rho, lambda = fitGT$pars$lambda,
    kernel = fitGT$pars$kernel, kernel_beta = fitGT$pars$kernel_beta,
    pareto_alpha = fitGT$pars$pareto_alpha, pareto_wmin = fitGT$pars$pareto_wmin,
    copula_family = fitGT$pars$copula_family, copula_theta = 0.0,
    regime = cfg$colas_regime
  )
  
  gen_rgg_fun <- function() gen_rgg(n = n, target_mean_deg = obs_feat$scalars$mean_deg, torus = TRUE)
  gen_cfg_fun <- function() gen_config_model_connected(deg_seq = deg_obs, max_tries = cfg$connect_tries)
  
  gens <- list(
    CoLaS = gen_colas,
    GeoTail_indep = gen_geotail,
    RGG = gen_rgg_fun,
    Config = gen_cfg_fun
  )
  
  # ---- simulate + score ----
  sim_rows <- list()
  dist_rows <- list()
  sim_deg_pool <- list()
  
  for (model in names(gens)) {
    message("Simulating model: ", model)
    deg_pool <- integer(0)
    
    for (b in seq_len(cfg$B_eval)) {
      set.seed(cfg$seed + 100000L + b + 1000L * match(model, names(gens)))
      
      gen0 <- function() gens[[model]]()
      
      gsim <- if (isTRUE(cfg$ensure_connected)) {
        ensure_connected_graph(gen0, tries = cfg$connect_tries)
      } else {
        g <- gen0()
        simplify(as.undirected(g, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
      }
      
      # Evaluate on LCC like OBS
      gsim <- largest_cc(gsim)
      
      feat_sim <- compute_features(gsim, breaks_deg = obs_feat$breaks, seed = cfg$seed + b)
      
      sim_rows[[length(sim_rows) + 1]] <- feature_to_row(feat_sim, model = model, rep = b)
      
      dist_dt <- gof_distances(obs_feat, feat_sim)
      dist_dt[, `:=`(model = model, rep = b)]
      dist_rows[[length(dist_rows) + 1]] <- dist_dt
      
      deg_pool <- c(deg_pool, feat_sim$deg)
    }
    sim_deg_pool[[model]] <- deg_pool
  }
  
  sim_dt <- rbindlist(sim_rows, fill = TRUE)
  dist_dt <- rbindlist(dist_rows, fill = TRUE)
  
  data.table::fwrite(sim_dt,  file.path(ds_dir, "tables", "sim_scalars_by_rep.csv"))
  data.table::fwrite(dist_dt, file.path(ds_dir, "tables", "gof_distances_by_rep.csv"))
  
  # aggregate full summaries (kept)
  sim_sum <- sim_dt[, lapply(.SD, function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))),
                    by = model,
                    .SDcols = c("mean_deg","clustering","assort","alpha_hat","tri","spec","mean_dist")]
  sim_sum_flat <- sim_sum[, {
    out <- list(model = model)
    for (nm in setdiff(names(sim_sum), "model")) {
      out[[paste0(nm, "_mean")]] <- get(nm)[1]
      out[[paste0(nm, "_sd")]]   <- get(nm)[2]
    }
    out
  }, by = model]
  
  dist_sum <- dist_dt[, lapply(.SD, function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))),
                      by = model,
                      .SDcols = setdiff(names(dist_dt), c("model","rep"))]
  dist_sum_flat <- dist_sum[, {
    out <- list(model = model)
    for (nm in setdiff(names(dist_sum), "model")) {
      out[[paste0(nm, "_mean")]] <- get(nm)[1]
      out[[paste0(nm, "_sd")]]   <- get(nm)[2]
    }
    out
  }, by = model]
  
  data.table::fwrite(sim_sum_flat,  file.path(ds_dir, "tables", "sim_scalars_summary.csv"))
  data.table::fwrite(dist_sum_flat, file.path(ds_dir, "tables", "gof_distances_summary.csv"))
  
  # ---- NEW: compact GOF tables (paper) ----
  comp <- summarize_compact_gof(dist_dt, metrics = cfg$gof_compact_metrics)
  data.table::fwrite(comp$numeric,   file.path(ds_dir, "tables", "gof_compact_summary_numeric.csv"))
  data.table::fwrite(comp$formatted, file.path(ds_dir, "tables", "gof_compact_summary_fmt.csv"))
  
  # ---- plots ----
  plot_degree_ccdf_overlay(
    obs_deg = obs_feat$deg,
    sim_deg_list = sim_deg_pool,
    out_png = file.path(ds_dir, "figures", "degree_ccdf_overlay.png"),
    title = paste0("E3: Degree CCDF (", ds$name, ") — OBS vs fitted models (", cfg$colas_regime, ")")
  )
  
  # scatter: clustering vs assort per replicate + OBS marker + legend
  plot_dt <- sim_dt[, .(assort, clustering, model)]
  plot_dt <- rbind(plot_dt,
                   data.table(assort = obs_feat$scalars$assort,
                              clustering = obs_feat$scalars$clustering,
                              model = "OBS"),
                   fill = TRUE)
  
  plot_dt[, model := factor(model, levels = c("OBS", names(gens)))]
  
  p_sc <- ggplot() +
    geom_point(data = plot_dt[model != "OBS"],
               aes(x = assort, y = clustering, shape = model),
               alpha = 0.8) +
    geom_point(data = plot_dt[model == "OBS"],
               aes(x = assort, y = clustering, shape = model),
               size = 3.3, stroke = 1.2) +
    theme_minimal() +
    labs(title = paste0("E3: (assort, clustering) across replicates — ", ds$name,
                        " (", cfg$colas_regime, ")"),
         x = "degree assortativity", y = "global clustering", shape = "model")
  
  ggsave(file.path(ds_dir, "figures", "assort_vs_clustering.png"),
         p_sc, width = 7.5, height = 5, dpi = 200)
  
  invisible(list(obs = obs_feat, fits = fits_dt, sim = sim_dt, gof = dist_dt))
}

run_E3 <- function() {
  all_full <- list()
  all_comp_num <- list()
  all_comp_fmt <- list()
  
  for (ds in datasets) {
    run_E3_one_dataset(ds)
    ds_dir <- file.path(cfg$out_dir, ds$name)
    
    # full summary (existing)
    dist_sum <- data.table::fread(file.path(ds_dir, "tables", "gof_distances_summary.csv"))
    dist_sum[, dataset := ds$name]
    all_full[[length(all_full) + 1]] <- dist_sum
    
    # compact summaries (new)
    cnum <- data.table::fread(file.path(ds_dir, "tables", "gof_compact_summary_numeric.csv"))
    cnum[, dataset := ds$name]
    all_comp_num[[length(all_comp_num) + 1]] <- cnum
    
    cfmt <- data.table::fread(file.path(ds_dir, "tables", "gof_compact_summary_fmt.csv"))
    cfmt[, dataset := ds$name]
    all_comp_fmt[[length(all_comp_fmt) + 1]] <- cfmt
  }
  
  all_dt <- rbindlist(all_full, fill = TRUE)
  data.table::fwrite(all_dt, file.path(cfg$out_dir, "E3_all_datasets_gof_summary.csv"))
  
  comp_num <- rbindlist(all_comp_num, fill = TRUE)
  data.table::fwrite(comp_num, file.path(cfg$out_dir, "E3_all_datasets_gof_compact_summary_numeric.csv"))
  
  comp_fmt <- rbindlist(all_comp_fmt, fill = TRUE)
  data.table::fwrite(comp_fmt, file.path(cfg$out_dir, "E3_all_datasets_gof_compact_summary_fmt.csv"))
  
  message("\nDONE. E3 outputs in: ", normalizePath(cfg$out_dir))
  message("Wrote combined compact GOF tables:")
  message("  - E3_all_datasets_gof_compact_summary_numeric.csv")
  message("  - E3_all_datasets_gof_compact_summary_fmt.csv")
  
  invisible(list(full = all_dt, compact_numeric = comp_num, compact_fmt = comp_fmt))
}

# Execute when pasted / sourced:
if (sys.nframe() == 0) {
  run_E3()
}
