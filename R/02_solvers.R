#' @title Compute Multi-Scale Metagenomic Differences
#' @description Calculates OT distances at specified levels.
#' @export
asmot_calculate <- function(obj, level = "joint", datatype = "ra", mode = "UOT", rho = NULL) {

  if (datatype == "ra") {
    R <- obj@real_ra; S <- obj@synth_ra
  } else {
    R <- obj@real_counts; S <- obj@synth_counts
  }

  M <- obj@cost_matrix
  
  # Heuristic for rho and lambda
  if (is.null(rho)) rho <- quantile(M, 0.5)
  lambda <- 0.05 * quantile(M, 0.5)

  # --- Level 1: Univariate (1D Wasserstein) ---
  if (level == "univariate") {
    dists <- sapply(seq_along(obj@taxa_names), function(i) {
      transport::wasserstein1d(R[,i], S[,i])
    })
    return(data.frame(Taxon = obj@taxa_names, Distance = dists))
  }

  # --- Level 2: Joint (Sample Geometry) ---
  else if (level == "joint") {
    # 1. Compute Cost Matrix between ALL samples (Real vs Synthetic)
    # We essentially append S to R to compute the cross-distance
    D_samples <- as.matrix(dist(rbind(R, S)))
    n_r <- nrow(R); n_s <- nrow(S)
    
    # Extract the Real-vs-Synthetic block
    Cost_RS <- D_samples[1:n_r, (n_r + 1):(n_r + n_s)]

    # Weights (Uniform)
    w_r <- rep(1/n_r, n_r)
    w_s <- rep(1/n_s, n_s)

    if (mode == "OT") {
      # Exact OT (Balanced) using 'transport' package
      res <- transport::transport(w_r, w_s, Cost_RS)
      return(sum(res$mass * res$cost))
    } else {
      # Unbalanced Sinkhorn (Robust Internal Solver)
      res <- internal_uot_sinkhorn(w_r, w_s, Cost_RS, lambda = lambda, rho = rho)
      return(res)
    }
  }

  # --- Level 3: Structural (Network Topology) ---
  else if (level == "structural") {
    # Since taxa are ALIGNED (same columns), we can compare correlation structures directly.
    # This is equivalent to Gromov-Wasserstein assuming the Identity mapping is optimal.
    D_R <- cor_to_dist(R)
    D_S <- cor_to_dist(S)
    
    # Frobenius Norm of the difference between distance matrices
    # Measures how much the "web of connections" has distorted
    diff <- norm(D_R - D_S, type = "F")
    return(diff)
  }
}

#' @title Audit High-Order Interactions
#' @export
asmot_k_scan <- function(obj, k = 3, n_subsamples = 50, mode = "UOT", rho = 1.0, seed = 42) {
  set.seed(seed)
  R <- obj@real_ra; S <- obj@synth_ra
  n_taxa <- ncol(R)
  distances <- numeric(n_subsamples)
  lambda <- 0.05 * median(obj@cost_matrix)

  for (i in 1:n_subsamples) {
    idx <- sample(1:n_taxa, k)
    R_sub <- R[, idx, drop=FALSE]; S_sub <- S[, idx, drop=FALSE]
    M_sub <- obj@cost_matrix[idx, idx]
    
    a <- normalize_vec(colMeans(R_sub))
    b <- normalize_vec(colMeans(S_sub))

    # Use internal solver for consistency
    distances[i] <- internal_uot_sinkhorn(a, b, M_sub, lambda = lambda, rho = rho)
  }
  
  res <- list(k = k, mean_distance = mean(distances), all_distances = distances)
  class(res) <- "asmot_k_result"
  return(res)
}

# --- ROBUST INTERNAL SOLVERS (No external dependencies) ---

#' Internal Unbalanced Sinkhorn Solver
#' @keywords internal
internal_uot_sinkhorn <- function(a, b, cost, lambda, rho, n_iter = 100) {
  # Add epsilon to avoid division by zero
  a <- pmax(a, 1e-10)
  b <- pmax(b, 1e-10)
  
  # Gibbs Kernel
  K <- exp(-cost / lambda)
  
  # Scaling power (Soft marginal constraint)
  fi <- rho / (rho + lambda)
  
  # Init potentials
  u <- rep(1, length(a))
  v <- rep(1, length(b))
  
  for(i in 1:n_iter) {
    # Sinkhorn iterations
    Kv <- as.vector(K %*% v)
    u <- (a / Kv)^fi
    
    Ktu <- as.vector(t(K) %*% u)
    v <- (b / Ktu)^fi
  }
  
  # Compute Transport Plan
  gamma <- diag(u) %*% K %*% diag(v)
  
  # Return Cost
  return(sum(gamma * cost))
}
