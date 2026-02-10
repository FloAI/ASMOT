library(ASMOT)
library(ggplot2)
library(dplyr)
library(gridExtra)

set.seed(123) # Ensure reproducibility

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

interpret_asmot <- function(res, label) {
  cat("\n========================================\n")
  cat(paste("REPORT CARD:", label, "\n"))
  cat("========================================\n")
  cat(sprintf("Realism Score:   %.2f%%\n", res$realism_score))
  cat(sprintf("P-value:         %.4f\n", res$p_value))

  comps <- res$components
  # Normalize structural error for display if needed
  if (comps["Structural"] > 1.0) comps["Structural"] <- comps["Structural"] / 100

  print(round(comps, 4))
  cat("----------------------------------------\n")
}

# --- Robust Scan Function (Returns Data for Unified Plotting) ---
compute_scan_data <- function(obj, label, max_k = 10) {
  n_taxa <- ncol(obj@real_ra)
  run_k <- min(max_k, n_taxa)

  if (run_k < 2) return(NULL)

  all_results <- data.frame()

  for (k in 1:run_k) {
    errors <- numeric(20) # 20 reps per k
    for (i in 1:20) {
      idx <- sample(seq_len(n_taxa), k)
      R_sub <- obj@real_ra[, idx, drop=FALSE]
      S_sub <- obj@synth_ra[, idx, drop=FALSE]

      # Error Metric: Mean Diff + Scaled Covariance Diff
      diff_mean <- sum(abs(colMeans(R_sub) - colMeans(S_sub)))
      diff_cov <- if (k > 1) sum(abs(cov(R_sub) - cov(S_sub))) / (k*k) else 0

      errors[i] <- diff_mean + diff_cov
    }
    # Store data for Boxplot
    batch <- data.frame(k = as.factor(k), Error = errors, Scenario = label)
    all_results <- rbind(all_results, batch)
  }
  return(all_results)
}

# ------------------------------------------------------------------------------
# 2. LOAD & AUDIT
# ------------------------------------------------------------------------------
message("Loading Datasets...")
load("data/amgut1.filt.rda")
load("data/amgut1.filt.cs.rda")
load("data/synth_zinb.rda")
load("data/synth_htln.rda")

# Create Objects
obj_filt_zinb <- create_asmot(amgut1.filt, synth_zinb)
obj_filt_htln <- create_asmot(amgut1.filt, synth_htln)
obj_cs_zinb   <- create_asmot(amgut1.filt.cs, synth_zinb)
obj_cs_htln   <- create_asmot(amgut1.filt.cs, synth_htln)

# Run Audits (B=2000, Equal Weights)
message("\n--- Running Audits ---")
B_reps <- 2000
weights <- c(1, 1, 1)

res_filt_zinb <- asmot_audit(obj_filt_zinb, B = B_reps, weights = weights)
res_filt_htln <- asmot_audit(obj_filt_htln, B = B_reps, weights = weights)
res_cs_zinb   <- asmot_audit(obj_cs_zinb, B = B_reps, weights = weights)
res_cs_htln   <- asmot_audit(obj_cs_htln, B = B_reps, weights = weights)

# Print Reports
interpret_asmot(res_filt_zinb, "Standard Data vs ZINB")
interpret_asmot(res_filt_htln, "Standard Data vs HTLN")
interpret_asmot(res_cs_zinb,   "CS Data vs ZINB")
interpret_asmot(res_cs_htln,   "CS Data vs HTLN")

# ------------------------------------------------------------------------------
# 3. PLOT 1: SUMMARY BAR CHART
# ------------------------------------------------------------------------------
results_df <- data.frame(
  Ground_Truth = c(rep("Standard", 2), rep("CS (Complex)", 2)),
  Model        = c("ZINB", "HTLN", "ZINB", "HTLN"),
  Realism      = c(res_filt_zinb$realism_score, res_filt_htln$realism_score,
                   res_cs_zinb$realism_score,   res_cs_htln$realism_score)
)

p1 <- ggplot(results_df, aes(x = Ground_Truth, y = Realism, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Realism, 1)), position = position_dodge(0.8), vjust = -0.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 105) +
  labs(title = "Benchmark Results: Realism Score", y = "ASMOT Score (%)") +
  theme_minimal()
print(p1)

# ------------------------------------------------------------------------------
# 4. PLOT 2: HIGH-ORDER SCAN (Regular Scales)
# ------------------------------------------------------------------------------
message("\n--- Running k-Scans (Box Plot / Unified Scale) ---")

# Compute data individually, then combine
all_scan_data <- rbind(
  compute_scan_data(obj_filt_zinb, "Standard vs ZINB"),
  compute_scan_data(obj_filt_htln, "Standard vs HTLN"),
  compute_scan_data(obj_cs_zinb,   "CS vs ZINB"),
  compute_scan_data(obj_cs_htln,   "CS vs HTLN")
)

# Plot with Facets -> Automatically enforces fixed scales
p2 <- ggplot(all_scan_data, aes(x = k, y = Error, fill = Scenario)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~Scenario) +  # Fixed scales by default
  scale_fill_brewer(palette = "Set2") +
  labs(title = "High-Order Consistency Scan",
       subtitle = "Boxplots show error stability (Unified Scale)",
       x = "Complexity (k Taxa)", y = "Geometric Error") +
  theme_minimal() +
  theme(legend.position = "none")

print(p2)

# ------------------------------------------------------------------------------
# 5. PLOT 3: TURING TEST (Regular Scales)
# ------------------------------------------------------------------------------
message("\n--- Running Classifiers (Unified Scale) ---")

rf_filt_zinb <- asmot_classifier(obj_filt_zinb)
rf_filt_htln <- asmot_classifier(obj_filt_htln)
rf_cs_zinb   <- asmot_classifier(obj_cs_zinb)
rf_cs_htln   <- asmot_classifier(obj_cs_htln)

combine_rf <- function(res, name) { df <- res$Data; df$Scenario <- name; df }

all_rf_data <- rbind(
  combine_rf(rf_filt_zinb, "Standard vs ZINB"),
  combine_rf(rf_filt_htln, "Standard vs HTLN"),
  combine_rf(rf_cs_zinb,   "CS vs ZINB"),
  combine_rf(rf_cs_htln,   "CS vs HTLN")
)

# Print Accuracy Table
acc_df <- data.frame(
  Scenario = c("Standard vs ZINB", "Standard vs HTLN", "CS vs ZINB", "CS vs HTLN"),
  Accuracy = c(rf_filt_zinb$Accuracy, rf_filt_htln$Accuracy, rf_cs_zinb$Accuracy, rf_cs_htln$Accuracy)
)
print(acc_df)

# Plot with Facets -> Enforces fixed scales (Removed 'scales="free"')
p3 <- ggplot(all_rf_data, aes(x = OT_Distance, fill = Label)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Scenario) + # Fixed Scales (Regular)
  scale_fill_manual(values = c("Real" = "black", "Synthetic" = "#E69F00")) +
  labs(title = "The Turing Test: Classifier Separation",
       subtitle = "X-axis fixed to show absolute magnitude of deviation",
       x = "Distance to Centroid",
       y = "Density") +
  theme_minimal()

print(p3)
