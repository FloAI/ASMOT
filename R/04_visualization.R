#' @import ggplot2
#' @import reshape2

#' @title Plot Audit Dashboard
#' @export
asmot_plot_dashboard <- function(audit_res, bias_res = NULL) {
  
  
  
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

#' @title Plot Dimension Scan
#' @export
asmot_plot_scan <- function(obj, max_k = 8, seed = 42) {
  
  
  
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

#' @export
asmot_plot_classifier <- function(class_res) {
  
  
  
  ggplot(class_res$Data, aes(x = OT_Distance, fill = Label)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Classifier Accuracy:", round(class_res$Accuracy, 2)),
         subtitle = "Overlap = Realism. Separation = Artifacts.",
         x = "Distance to Real Center", y = "Density") +
    theme_minimal()
}

#' @export
asmot_plot_sparsity <- function(obj) {
  calc_stats <- function(mat, label) {
    data.frame(Mean = colMeans(mat), Sparsity = colMeans(mat == 0), Dataset = label)
  }
  df <- rbind(calc_stats(obj@real_ra, "Real"), calc_stats(obj@synth_ra, "Synthetic"))
  
  
  
  ggplot(df, aes(x = Mean, y = Sparsity, color = Dataset)) +
    geom_point(alpha = 0.5) + scale_x_log10() +
    labs(title = "Sparsity vs Abundance", x = "Log10 Mean Abundance", y = "Frequency of Zeros") +
    theme_minimal()
}
