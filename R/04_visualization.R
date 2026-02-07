#' @title Plot Audit Dashboard
#' @description 
#' Visualizes the "Realism Gap" by plotting the null distribution of real-vs-real distances
#' against the observed real-vs-synthetic distance.
#'
#' @param audit_res Result object from \code{asmot_audit}.
#' @param bias_res (Optional) Result object from \code{asmot_bias_check}.
#' @return A ggplot object.
#' @export
asmot_plot_dashboard <- function(audit_res, bias_res = NULL) {
  
  # (Implementation as provided in previous steps)
  p1 <- ggplot(data.frame(d = audit_res$null_dist), aes(d)) +
    geom_density(fill = "#69b3a2", alpha = 0.6) +
    geom_vline(xintercept = audit_res$obs_dist, color = "red", linetype = "dashed", size=1) +
    labs(title = paste0("Realism Score: ", audit_res$realism_score, "%"),
         subtitle = paste0("P-value: ", audit_res$p_value),
         x = "Joint OT Distance", y = "Density") +
    theme_minimal()

  if(!is.null(bias_res)) {
    df_bias <- data.frame(Type = c("Relative Abundance", "Raw Counts"),
                          Distance = c(bias_res$Dist_RA, bias_res$Dist_Counts))
    p2 <- ggplot(df_bias, aes(x=Type, y=Distance, fill=Type)) +
      geom_bar(stat="identity", width=0.5) +
      theme_minimal() + theme(legend.position="none")
    return(list(Realism_Plot = p1, Bias_Plot = p2))
  }
  return(p1)
}

#' @title Plot High-Order Interaction Scan
#' @description 
#' Wrapper that runs \code{asmot_k_scan} for dimensions 1 to \code{max_k} and plots the error profile.
#' Rising error indicates failure to capture complex ecological motifs.
#'
#' @param obj An \code{ASMOT} object.
#' @param max_k Integer. Maximum dimension to scan (default = 8).
#' @param seed Integer. Random seed.
#' @return An object of class \code{asmot_k_scan} containing the plot and raw data.
#' @export
asmot_plot_scan <- function(obj, max_k = 8, seed = 42) {
  
  # (Implementation as provided in previous steps)
  results <- data.frame()
  for (k in 1:max_k) {
    res <- asmot_k_scan(obj, k = k, n_subsamples = 30, seed = seed)
    results <- rbind(results, data.frame(k = k, Error = res$all_distances))
  }
  
  p <- ggplot(results, aes(x = factor(k), y = Error, fill = factor(k))) +
    geom_boxplot() +
    labs(title = "High-Order Interaction Audit", x = "Dimension (k)", y = "OT Distance") +
    theme_minimal() + theme(legend.position = "none")
  
  out <- list(plot = p, data = results)
  class(out) <- "asmot_k_scan"
  return(out)
}
