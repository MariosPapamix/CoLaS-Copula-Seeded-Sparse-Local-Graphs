###############################################
## E1: Light-tail impossibility vs tail inheritance
## - Fixed-range CoLaS: edge rule (11)
## - CoLaS-HT:          edge rule (12)
## Using the same heavy-tailed F_W and same geometric kernel
###############################################

set.seed(1)

## --------------------------
## 1) Parameters
## --------------------------
n   <- 5000
d   <- 2

alpha <- 2.5   # Pareto tail index for W (heavy-tailed)
xm    <- 1

rho   <- 10                     # sparsity scale
eps   <- (rho / n)^(1/d)         # makes rho_n = n*eps^d = rho (exact)

target_mean_deg <- 20            # optional: roughly match mean degree across regimes

## --------------------------
## 2) Heavy-tailed weights + latent positions on T^2
## (Here: independence copula for simplicity; E0 is about tail regime, not assortativity.)
## --------------------------
rpareto <- function(n, alpha, xm=1) xm * (runif(n)^(-1/alpha))
W <- rpareto(n, alpha, xm)

X <- cbind(runif(n), runif(n))   # positions on [0,1)^2 (torus via wrap-around distance)

## --------------------------
## 3) Geometric kernel k(u)=1{||u||<=1}
## (bounded + compact support, which is the regime where fixed-range degrees are provably light-tailed)
## --------------------------
k_indicator <- function(r) as.numeric(r <= 1)  # r = ||u||

torus_dx <- function(a, b){
  dx <- abs(a - b)
  pmin(dx, 1 - dx)
}

## --------------------------
## 4) Rough lambda calibration (optional)
## - HT: uses small-eps approximation
## - Fixed-range: Monte Carlo solve for mean degree target
## --------------------------
calibrate_lambdas <- function(W, rho, target_mean_deg, m_mc=200000){
  if(target_mean_deg >= pi*rho){
    stop("For indicator kernel in fixed-range, mean degree is bounded by ~pi*rho; choose target_mean_deg < pi*rho.")
  }
  
  Ew <- mean(W)
  
  ## HT: mean degree ≈ (1-exp(-λ/ρ)) * pi*rho*E[W]^2
  p0 <- target_mean_deg / (pi * rho * Ew^2)
  p0 <- max(min(p0, 0.95), 1e-6)
  lambda_HT <- -rho * log1p(-p0)
  
  ## Fixed-range: mean degree ≈ pi*rho * E[1 - exp(-(λ/ρ) W1 W2)]
  a <- sample(W, m_mc, replace=TRUE)
  b <- sample(W, m_mc, replace=TRUE)
  rhs <- target_mean_deg / (pi * rho)
  f <- function(lam) mean(1 - exp(-(lam/rho)*a*b)) - rhs
  
  lo <- 1e-6; hi <- 100
  while(f(hi) < 0) hi <- hi * 2
  lambda_FR <- uniroot(f, c(lo, hi))$root
  
  list(lambda_FR=lambda_FR, lambda_HT=lambda_HT)
}

lams <- calibrate_lambdas(W, rho, target_mean_deg)
lambda_FR <- lams$lambda_FR
lambda_HT <- lams$lambda_HT

cat(sprintf("Using eps=%.4f (rho_n=rho=%.1f)\n", eps, rho))
cat(sprintf("lambda_FR=%.3f, lambda_HT=%.3f (target mean degree ≈ %g)\n\n",
            lambda_FR, lambda_HT, target_mean_deg))

## --------------------------
## 5) Fast degree-only simulators via cell grid
## --------------------------
make_grid <- function(X, cell_size){
  m  <- ceiling(1/cell_size)
  cx <- floor(X[,1]/cell_size) %% m
  cy <- floor(X[,2]/cell_size) %% m
  id <- cx + m*cy + 1L
  list(m=m, cx=cx, cy=cy, grid=split(seq_len(nrow(X)), id))
}

cell_id <- function(cx, cy, m) ((cx %% m) + m*(cy %% m) + 1L)

## Fixed-range CoLaS (11) with indicator kernel:
## p_ij = 1 - exp(-(λ/ρ) W_i W_j) if ||X_i-X_j|| <= eps; else 0.
simulate_fixed_range_deg <- function(W, X, eps, rho, lambda){
  n <- length(W)
  g <- make_grid(X, cell_size=eps)
  m <- g$m; cx <- g$cx; cy <- g$cy; grid <- g$grid
  
  deg <- integer(n)
  
  for(i in 1:(n-1)){
    for(dx in -1:1){
      for(dy in -1:1){
        cid  <- cell_id(cx[i]+dx, cy[i]+dy, m)
        cand <- grid[[as.character(cid)]]
        if(is.null(cand)) next
        cand <- cand[cand > i]
        if(!length(cand)) next
        
        dxv  <- torus_dx(X[cand,1], X[i,1])
        dyv  <- torus_dx(X[cand,2], X[i,2])
        dist <- sqrt(dxv*dxv + dyv*dyv)
        
        ok <- (dist <= eps)
        if(!any(ok)) next
        j <- cand[ok]
        
        p <- 1 - exp(-(lambda/rho) * W[i] * W[j])
        keep <- runif(length(j)) < p
        if(any(keep)){
          jkeep <- j[keep]
          deg[i]     <- deg[i] + length(jkeep)
          deg[jkeep] <- deg[jkeep] + 1L
        }
      }
    }
  }
  deg
}

## CoLaS-HT (12) with indicator kernel:
## p_ij = 1 - exp(-(λ/ρ)) if ||X_i-X_j|| <= eps*(W_i W_j)^(1/d); else 0.
simulate_ht_deg <- function(W, X, eps, rho, lambda, heavy_q=0.99){
  n <- length(W)
  g <- make_grid(X, cell_size=eps)
  m <- g$m; cx <- g$cx; cy <- g$cy; grid <- g$grid
  
  w_cap    <- as.numeric(quantile(W, probs=heavy_q))
  heavy    <- which(W > w_cap)
  is_heavy <- logical(n); is_heavy[heavy] <- TRUE
  
  p0 <- 1 - exp(-(lambda/rho))   # constant inside the weight-scaled ball
  
  deg <- integer(n)
  
  for(i in 1:(n-1)){
    if(is_heavy[i]){
      ## brute force over all j>i (heavy nodes can have large radii)
      cand <- (i+1):n
      dxv  <- torus_dx(X[cand,1], X[i,1])
      dyv  <- torus_dx(X[cand,2], X[i,2])
      dist <- sqrt(dxv*dxv + dyv*dyv)
      
      thresh <- eps * sqrt(W[i] * W[cand])  # since d=2
      ok <- (dist <= thresh)
      if(!any(ok)) next
      j <- cand[ok]
      
      keep <- runif(length(j)) < p0
      if(any(keep)){
        jkeep <- j[keep]
        deg[i]     <- deg[i] + length(jkeep)
        deg[jkeep] <- deg[jkeep] + 1L
      }
      
    } else {
      ## (a) exact local search among non-heavy nodes (since W_j <= w_cap bounds the radius)
      Ri <- eps * sqrt(W[i] * w_cap)
      h  <- ceiling(Ri/eps)
      
      for(dx in -h:h){
        for(dy in -h:h){
          cid  <- cell_id(cx[i]+dx, cy[i]+dy, m)
          cand <- grid[[as.character(cid)]]
          if(is.null(cand)) next
          cand <- cand[cand > i & !is_heavy[cand]]
          if(!length(cand)) next
          
          dxv  <- torus_dx(X[cand,1], X[i,1])
          dyv  <- torus_dx(X[cand,2], X[i,2])
          dist <- sqrt(dxv*dxv + dyv*dyv)
          
          thresh <- eps * sqrt(W[i] * W[cand])
          ok <- (dist <= thresh)
          if(!any(ok)) next
          j <- cand[ok]
          
          keep <- runif(length(j)) < p0
          if(any(keep)){
            jkeep <- j[keep]
            deg[i]     <- deg[i] + length(jkeep)
            deg[jkeep] <- deg[jkeep] + 1L
          }
        }
      }
      
      ## (b) edges from i to heavy nodes with j>i (few, so brute force)
      cand <- heavy[heavy > i]
      if(length(cand)){
        dxv  <- torus_dx(X[cand,1], X[i,1])
        dyv  <- torus_dx(X[cand,2], X[i,2])
        dist <- sqrt(dxv*dxv + dyv*dyv)
        
        thresh <- eps * sqrt(W[i] * W[cand])
        ok <- (dist <= thresh)
        if(any(ok)){
          j <- cand[ok]
          keep <- runif(length(j)) < p0
          if(any(keep)){
            jkeep <- j[keep]
            deg[i]     <- deg[i] + length(jkeep)
            deg[jkeep] <- deg[jkeep] + 1L
          }
        }
      }
    }
  }
  deg
}

## Simulate the two graphs (same W, same X)
deg_FR <- simulate_fixed_range_deg(W, X, eps, rho, lambda_FR)
deg_HT <- simulate_ht_deg(W, X, eps, rho, lambda_HT, heavy_q=0.99)

cat(sprintf("Realized mean degree: fixed-range=%.2f, HT=%.2f\n\n",
            mean(deg_FR), mean(deg_HT)))

## --------------------------
## 6) E0 diagnostics: log–log CCDF and Hill tail-index plots
## --------------------------
ccdf <- function(x){
  x <- x[is.finite(x) & x > 0]
  x <- sort(x)
  n <- length(x)
  ux <- unique(x)
  idx <- match(ux, x)
  surv <- (n - idx + 1) / n
  data.frame(x=ux, surv=surv)
}

hill_alpha <- function(x, k_min=20, k_max=400){
  x <- x[is.finite(x) & x > 0]
  x <- sort(x, decreasing=TRUE)
  n <- length(x)
  k_max <- min(k_max, n-1)
  ks <- k_min:k_max
  xi_hat <- sapply(ks, function(k){
    mean(log(x[1:k])) - log(x[k+1])
  })
  data.frame(k=ks, alpha=1/xi_hat)
}

S_W  <- ccdf(W)
S_FR <- ccdf(deg_FR)
S_HT <- ccdf(deg_HT)

par(mfrow=c(1,3))
plot(log10(S_W$x),  log10(S_W$surv),  type="l",
     xlab="log10(x)", ylab="log10(P(X>=x))",
     main="Weights W (heavy-tailed)")
plot(log10(S_FR$x), log10(S_FR$surv), type="l",
     xlab="log10(k)", ylab="log10(P(D>=k))",
     main="Degrees: fixed-range CoLaS (11)")
plot(log10(S_HT$x), log10(S_HT$surv), type="l",
     xlab="log10(k)", ylab="log10(P(D>=k))",
     main="Degrees: CoLaS-HT (12)")

H_W  <- hill_alpha(W,      k_min=50, k_max=800)
H_FR <- hill_alpha(deg_FR, k_min=20, k_max=400)
H_HT <- hill_alpha(deg_HT, k_min=20, k_max=400)

par(mfrow=c(1,3))
plot(H_W$k,  H_W$alpha,  type="l", xlab="k (top order stats)", ylab="Hill α-hat",
     main="Hill: W")
abline(h=alpha, lty=2)
plot(H_FR$k, H_FR$alpha, type="l", xlab="k", ylab="Hill α-hat",
     main="Hill: fixed-range degrees")
abline(h=alpha, lty=2)
plot(H_HT$k, H_HT$alpha, type="l", xlab="k", ylab="Hill α-hat",
     main="Hill: HT degrees")
abline(h=alpha, lty=2)

cat(sprintf("True alpha = %.2f\n", alpha))
cat(sprintf("Median Hill(alpha-hat): W=%.2f, fixed-range D=%.2f, HT D=%.2f\n",
            median(H_W$alpha), median(H_FR$alpha), median(H_HT$alpha)))
