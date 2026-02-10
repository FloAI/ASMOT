#' @title Perform Full ASMOT Audit
#' @export
asmot_audit <- function(obj, B = 100, weights = "auto", rho = NULL, fast_structural = TRUE) {
  
  # 1. Observed Distances
  u_res <- asmot_calculate(obj, "univariate")
  j_dist <- asmot_calculate(obj, "joint", rho = rho)
  s_dist <- asmot_calculate(obj, "structural")
  
  u_dist_mean <- mean(u_res$Distance)

  # 2. Bootstrap Null Distributions
  null_univ <- numeric(B)
  null_joint <- numeric(B)
  null_struct <- numeric(B)
  
  M <- obj@cost_matrix
  if(is.null(rho)) rho <- quantile(M, 0.5)
  lambda <- 0.05 * quantile(M, 0.5)

  for(i in 1:B) {
    # Split real data
    idx <- sample(1:nrow(obj@real_ra), floor(nrow(obj@real_ra)/2))
    R1 <- obj@real_ra[idx, ]; R2 <- obj@real_ra[-idx, ]

    # Univariate Null
    d_u_null <- sapply(1:ncol(R1), function(j) transport::wasserstein1d(R1[,j], R2[,j]))
    null_univ[i] <- mean(d_u_null)

    # Joint Null (Using our internal solver directly for speed)
    a <- normalize_vec(colMeans(R1))
    b <- normalize_vec(colMeans(R2))
    null_joint[i] <- internal_uot_sinkhorn(a, b, M, lambda = lambda, rho = rho)

    # Structural Null
    C1 <- cor_to_dist(R1); C2 <- cor_to_dist(R2)
    null_struct[i] <- norm(C1 - C2, type = "F") 
  }

  # 3. Adaptive Weights
  if (identical(weights, "auto")) {
    vars <- c(var(null_univ), var(null_joint), var(null_struct))
    vars[vars < 1e-9] <- 1e-9
    inv_vars <- 1 / vars
    final_weights <- inv_vars / sum(inv_vars)
  } else {
    final_weights <- weights / sum(weights)
  }
  names(final_weights) <- c("Univariate", "Joint", "Structural")

  # 4. Scores
  p_val_j <- mean(null_joint >= j_dist)
  p_val_s <- mean(null_struct >= s_dist) # Direct comparison since we use Frobenius for both
  final_p <- min(p_val_j, p_val_s)

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

#' @title Run Adversarial Classifier Audit
#' @export
asmot_classifier <- function(obj) {
  real_center <- normalize_vec(colMeans(obj@real_ra))
  lambda <- 0.05 * quantile(obj@cost_matrix, 0.5)
  rho <- quantile(obj@cost_matrix, 0.5)

  calc_dist <- function(mat) {
    sapply(1:nrow(mat), function(i) {
      internal_uot_sinkhorn(normalize_vec(mat[i,]), real_center, obj@cost_matrix, lambda, rho)
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
