#' @title Compute Realism Score with Adaptive Weights
#' @param weights Either "auto" (Inverse Null Variance) or a numeric vector of length 3.
#' @export
asmot_audit <- function(obj, B = 100, weights = "auto", rho = NULL, fast_structural = TRUE) {

  # 1. Observed Distances
  u_res <- asmot_calculate(obj, "univariate", mode = "UOT", rho = rho)
  j_dist <- asmot_calculate(obj, "joint", mode = "UOT", rho = rho)
  s_dist <- asmot_calculate(obj, "structural", mode = "UOT", rho = rho)
  
  u_dist_mean <- mean(u_res$Distance)

  # 2. Bootstrap Null Distributions
  null_univ <- numeric(B)
  null_joint <- numeric(B)
  null_struct <- numeric(B)
  
  M <- obj@cost_matrix
  if(is.null(rho)) rho <- estimate_rho(M)

  for(i in 1:B) {
    idx <- sample(1:nrow(obj@real_ra), floor(nrow(obj@real_ra)/2))
    R1 <- obj@real_ra[idx, ]; R2 <- obj@real_ra[-idx, ]

    # Univariate Null
    d_u_null <- sapply(1:ncol(R1), function(j) transport::wasserstein1d(R1[,j], R2[,j]))
    null_univ[i] <- mean(d_u_null)

    # Joint Null
    null_joint[i] <- T4transport::uot(normalize_vec(colMeans(R1)), normalize_vec(colMeans(R2)), M, rho = rho)$distance

    # Structural Null
    C1 <- cor_to_dist(R1); C2 <- cor_to_dist(R2)
    if(fast_structural) {
      null_struct[i] <- norm(C1 - C2, type = "F") 
    } else {
      null_struct[i] <- T4transport::ugw(C1, C2, rho = rho)$distance
    }
  }

  # 3. Compute Adaptive Weights
  if (identical(weights, "auto")) {
    vars <- c(
      Univariate = var(null_univ),
      Joint      = var(null_joint),
      Structural = var(null_struct)
    )
    vars[vars < 1e-9] <- 1e-9
    inv_vars <- 1 / vars
    final_weights <- inv_vars / sum(inv_vars)
  } else {
    if (length(weights) != 3) stop("Weights must be 'auto' or a numeric vector of length 3.")
    final_weights <- weights / sum(weights)
    names(final_weights) <- c("Univariate", "Joint", "Structural")
  }

  # 4. P-Value Calculation
  p_val_j <- mean(null_joint >= j_dist)
  p_val_s <- if(fast_structural) {
    obs_proxy <- norm(cor_to_dist(obj@real_ra) - cor_to_dist(obj@synth_ra), type = "F")
    mean(null_struct >= obs_proxy)
  } else {
    mean(null_struct >= s_dist)
  }
  final_p <- min(p_val_j, p_val_s)

  # 5. Realism Score Calculation
  s_u <- exp(-u_dist_mean)
  s_j <- exp(-j_dist)
  s_s <- exp(-s_dist)
  
  score <- (final_weights["Univariate"] * s_u + 
            final_weights["Joint"] * s_j + 
            final_weights["Structural"] * s_s) * 100

  res <- list(
    realism_score = round(score, 2),
    p_value = round(final_p, 4),
    weights_used = final_weights,
    null_dist = null_joint, 
    obs_dist = j_dist,
    components = c(Univariate=u_dist_mean, Joint=j_dist, Structural=s_dist)
  )
  class(res) <- "asmot_audit_res"
  return(res)
}

#' @title Classifier-Based Audit
#' @export
asmot_classifier <- function(obj) {
  real_center <- normalize_vec(colMeans(obj@real_ra))
  
  calc_dist <- function(mat) {
    sapply(1:nrow(mat), function(i) {
      T4transport::uot(normalize_vec(mat[i,]), real_center, obj@cost_matrix)$distance
    })
  }

  feat_real <- calc_dist(obj@real_ra)
  feat_synth <- calc_dist(obj@synth_ra)
  
  df <- data.frame(
    OT_Distance = c(feat_real, feat_synth),
    Label = factor(c(rep("Real", length(feat_real)), rep("Synthetic", length(feat_synth))))
  )

  model <- randomForest::randomForest(Label ~ OT_Distance, data = df)
  accuracy <- mean(predict(model) == df$Label)

  return(list(Accuracy = accuracy, Data = df))
}

#' @title Compare Counts vs Relative Abundance
#' @export
asmot_bias_check <- function(obj) {
  dist_ra <- asmot_calculate(obj, "joint", datatype = "ra", mode = "UOT")
  dist_count <- asmot_calculate(obj, "joint", datatype = "counts", mode = "UOT")
  ratio <- dist_count / dist_ra
  return(list(Dist_RA = dist_ra, Dist_Counts = dist_count, Ratio = ratio))
}
