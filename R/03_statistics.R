#' @title Perform Full ASMOT Audit
#' @export
asmot_audit <- function(obj, B = 100, weights = "auto", rho = NULL, fast_structural = TRUE) {
  
  u_res <- asmot_calculate(obj, "univariate", mode = "UOT", rho = rho)
  j_dist <- asmot_calculate(obj, "joint", mode = "UOT", rho = rho)
  s_dist <- asmot_calculate(obj, "structural", mode = "UOT", rho = rho)
  
  u_dist_mean <- mean(u_res$Distance)
  
  null_univ <- numeric(B)
  null_joint <- numeric(B)
  null_struct <- numeric(B)
  
  M <- obj@cost_matrix
  if(is.null(rho)) rho <- estimate_rho(M)

  for(i in 1:B) {
    idx <- sample(1:nrow(obj@real_ra), floor(nrow(obj@real_ra)/2))
    R1 <- obj@real_ra[idx, ]; R2 <- obj@real_ra[-idx, ]

    d_u_null <- sapply(1:ncol(R1), function(j) transport::wasserstein1d(R1[,j], R2[,j]))
    null_univ[i] <- mean(d_u_null)

    # Use internal solver for Joint Null
    D_samples <- as.matrix(dist(rbind(normalize_vec(colMeans(R1)), normalize_vec(colMeans(R2)))))
    # (Simplified call for internal solver)
    # Ideally we just call asmot_calculate logic here, but for speed we keep it inline or simple
    # Note: For simplicity in the loop, we can just call T4transport::sinkhorn or our internal one
    # But to match your code structure, let's keep it simple:
    
    # Recalculate Joint Distance for Null
    # We construct a mini-problem
    a <- normalize_vec(colMeans(R1))
    b <- normalize_vec(colMeans(R2))
    null_joint[i] <- internal_uot_sinkhorn(a, b, M, lambda=0.1*median(M), rho=rho)

    C1 <- cor_to_dist(R1); C2 <- cor_to_dist(R2)
    if(fast_structural) {
      null_struct[i] <- norm(C1 - C2, type = "F") 
    } else {
      # FIXED: Use gw instead of ugw
      null_struct[i] <- T4transport::gw(C1, C2)$distance
    }
  }

  if (identical(weights, "auto")) {
    vars <- c(Univariate = var(null_univ), Joint = var(null_joint), Structural = var(null_struct))
    vars[vars < 1e-9] <- 1e-9
    inv_vars <- 1 / vars
    final_weights <- inv_vars / sum(inv_vars)
  } else {
    final_weights <- weights / sum(weights)
    names(final_weights) <- c("Univariate", "Joint", "Structural")
  }

  p_val_j <- mean(null_joint >= j_dist)
  p_val_s <- if(fast_structural) {
    obs_proxy <- norm(cor_to_dist(obj@real_ra) - cor_to_dist(obj@synth_ra), type = "F")
    mean(null_struct >= obs_proxy)
  } else {
    mean(null_struct >= s_dist)
  }
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
  
  calc_dist <- function(mat) {
    sapply(1:nrow(mat), function(i) {
      # Use internal solver here too
      internal_uot_sinkhorn(normalize_vec(mat[i,]), real_center, obj@cost_matrix, 
                            lambda=0.1*median(obj@cost_matrix), rho=1.0)
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
