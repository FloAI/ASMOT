#' @title Compute Multi-Scale Metagenomic Differences
#' @description Calculates OT distances at specified levels.
#' @export
asmot_calculate <- function(obj, level = "joint", datatype = "ra", mode = "UOT", rho = NULL) {

  if (datatype == "ra") {
    R <- obj@real_ra; S <- obj@synth_ra
  } else {
    R <- obj@real_counts; S <- obj@synth_counts
    if (level %in% c("joint", "structural")) {
      warning("Using raw counts for Joint/Structural analysis is not recommended. Consider 'ra' or preprocess first.")
    }
  }

  M <- obj@cost_matrix
  # Default rho estimation if NULL
  if (is.null(rho) && mode == "UOT") rho <- estimate_rho(M)
  # Default lambda (entropy) for Sinkhorn
  lambda <- 0.1 * median(M)

  # --- Level 1: Univariate ---
  if (level == "univariate") {
    # Use transport::wasserstein1d for reliable 1D computation
    dists <- sapply(seq_along(obj@taxa_names), function(i) {
      transport::wasserstein1d(R[,i], S[,i])
    })
    return(data.frame(Taxon = obj@taxa_names, Distance = dists))
  }

  # --- Level 2: Joint ---
  else if (level == "joint") {
    D_samples <- as.matrix(dist(rbind(R, S)))
    n_r <- nrow(R); n_s <- nrow(S)
    Cost_SS <- D_samples[1:n_r, (n_r + 1):(n_r + n_s)]

    w_r <- rep(1/n_r, n_r); w_s <- rep(1/n_s, n_s)

    if (mode == "OT") {
      # Exact OT using transport package
      res <- transport::transport(w_r, w_s, Cost_SS)
      return(sum(res$mass * res$cost))
    } else {
      # Unbalanced OT using our internal solver
      res <- internal_uot_sinkhorn(w_r, w_s, Cost_SS, lambda = lambda, rho = rho)
      return(res)
    }
  }

  # --- Level 3: Structural ---
  else if (level == "structural") {
    D_R <- cor_to_dist(R)
    D_S <- cor_to_dist(S)

    if (mode == "OT") {
      # Standard GW (Balanced)
      res <- T4transport::gw(D_R, D_S)
      return(res$distance)
    } else {
      # Unbalanced GW
      res <- T4transport::ugw(D_R, D_S, rho = rho)
      return(res$distance)
    }
  }
}

#' @title Audit High-Order (k-Dimensional) Interactions
#' @export
asmot_k_scan <- function(obj, k = 3, n_subsamples = 50, mode = "UOT", rho = 1.0, seed = 42) {
  
  set.seed(seed)
  R <- obj@real_ra; S <- obj@synth_ra
  n_taxa <- ncol(R)
  if (k > n_taxa) stop("Dimension k cannot exceed total number of taxa.")

  distances <- numeric(n_subsamples)
  lambda <- 0.1 * median(obj@cost_matrix)

  for (i in 1:n_subsamples) {
    idx <- sample(1:n_taxa, k)
    R_sub <- R[, idx, drop=FALSE]; S_sub <- S[, idx, drop=FALSE]
    M_sub <- obj@cost_matrix[idx, idx]
    
    a <- normalize_vec(colMeans(R_sub))
    b <- normalize_vec(colMeans(S_sub))

    if (mode == "OT") {
      # T4transport::sinkhorn for balanced case
      res <- T4transport::sinkhorn(a, b, M_sub, lambda = lambda)
      distances[i] <- res$distance
    } else {
      # Use internal UOT solver
      distances[i] <- internal_uot_sinkhorn(a, b, M_sub, lambda = lambda, rho = rho)
    }
  }
  
  res <- list(k = k, mean_distance = mean(distances), all_distances = distances)
  class(res) <- "asmot_k_result"
  return(res)
}

# --- Internal Solver (The Fix) ---

#' Internal Unbalanced Sinkhorn Solver
#' Implements Sinkhorn-Knopp algorithm with marginal relaxation (KL penalty).
#' @keywords internal
internal_uot_sinkhorn <- function(a, b, cost, lambda, rho, n_iter = 100) {
  # Avoid log(0)
  a <- pmax(a, 1e-10)
  b <- pmax(b, 1e-10)
  
  # Kernel
  K <- exp(-cost / lambda)
  
  # Scaling factor for unbalanced update (rho / (rho + epsilon))
  fi <- rho / (rho + lambda)
  
  # Initialization
  u <- rep(1, length(a))
  v <- rep(1, length(b))
  
  for(i in 1:n_iter) {
    # Update u
    Kv <- as.vector(K %*% v)
    u <- (a / Kv)^fi
    
    # Update v
    Ktu <- as.vector(t(K) %*% u)
    v <- (b / Ktu)^fi
  }
  
  # Calculate Transport Plan Gamma
  gamma <- diag(u) %*% K %*% diag(v)
  
  # Return Total Transport Cost (Primal Objective approximation)
  return(sum(gamma * cost))
}
