#' @title Perform Full ASMOT Audit
#' @description 
#' Runs the complete Multi-Scale Audit (Univariate, Joint, Structural). Computes a global "Realism Score"
#' using Adaptive Weighting based on the variance of the null distribution.
#'
#' @param obj An \code{ASMOT} object.
#' @param B Integer. Number of bootstrap permutations for the null distribution (default = 100).
#' @param weights Character or Numeric. If \code{"auto"} (default), weights are inversely proportional
#'        to the null variance. If numeric vector of length 3, weights are fixed.
#' @param rho Numeric. Mass penalty parameter for UOT.
#' @param fast_structural Logical. If TRUE (default), uses Frobenius Norm as a fast proxy for 
#'        Gromov-Wasserstein distance during the bootstrap loop to speed up computation.
#'
#' @return An object of class \code{asmot_audit_res} containing scores, p-values, and weights.
#' @export
asmot_audit <- function(obj, B = 100, weights = "auto", rho = NULL, fast_structural = TRUE) {
  # (Implementation as provided in previous steps)
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

#' @title Run Adversarial Classifier Audit
#' @description 
#' Trains a Random Forest classifier to distinguish Real from Synthetic samples based on their
#' Optimal Transport distance to the Real centroid. Serves as a "Turing Test" for the generator.
#'
#' @param obj An \code{ASMOT} object.
#' @return A list containing the classification Accuracy (0.5 = Indistinguishable) and the data used for plotting.
#' @export
asmot_classifier <- function(obj) {
  # (Implementation as provided in previous steps)
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
