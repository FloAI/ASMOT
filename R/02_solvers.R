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
  if (is.null(rho) && mode == "UOT") rho <- estimate_rho(M)

  # --- Level 1: Univariate ---
  if (level == "univariate") {
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
      res <- transport::transport(w_r, w_s, Cost_SS)
      return(sum(res$mass * res$cost))
    } else {
      res <- T4transport::uot(w_r, w_s, Cost_SS, rho = rho)
      return(res$distance)
    }
  }

  # --- Level 3: Structural ---
  else if (level == "structural") {
    D_R <- cor_to_dist(R)
    D_S <- cor_to_dist(S)

    if (mode == "OT") {
      res <- T4transport::gw(D_R, D_S)
    } else {
      res <- T4transport::ugw(D_R, D_S, rho = rho)
    }
    return(res$distance)
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

  for (i in 1:n_subsamples) {
    idx <- sample(1:n_taxa, k)
    R_sub <- R[, idx, drop=FALSE]; S_sub <- S[, idx, drop=FALSE]
    M_sub <- obj@cost_matrix[idx, idx]
    
    a <- normalize_vec(colMeans(R_sub))
    b <- normalize_vec(colMeans(S_sub))

    if (mode == "OT") {
      res <- T4transport::sinkhorn(a, b, M_sub, lambda=0.05)
      distances[i] <- res$distance
    } else {
      res <- T4transport::uot(a, b, M_sub, rho = rho)
      distances[i] <- res$distance
    }
  }
  
  res <- list(k = k, mean_distance = mean(distances), all_distances = distances)
  class(res) <- "asmot_k_result"
  return(res)
}
