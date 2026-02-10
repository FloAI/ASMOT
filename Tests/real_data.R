library(ASMOT)
library(ggplot2)
library(dplyr)
library(gridExtra) # Needed to arrange the 4 k-scan plots together

set.seed(123) # Ensure reproducibility

# --- Helper Function: Interpret Results ---
interpret_asmot <- function(res, label) {
  cat("\n========================================\n")
  cat(paste("REPORT CARD:", label, "\n"))
  cat("========================================\n")

  # 1. Realism Score Check
  score <- res$realism_score
  cat(sprintf("Realism Score:   %.2f%%\n", score))
  if (score > 90) {
    cat("  -> STATUS: PUBLICATION READY (Indistinguishable)\n")
  } else if (score > 70) {
    cat("  -> STATUS: GOOD (Detectable artifacts exist)\n")
  } else {
    cat("  -> STATUS: FAIL (Significant deviation from real data)\n")
  }

  # 2. P-Value Check
  pval <- res$p_value
  cat(sprintf("P-value:         %.4f\n", pval))
  if (pval < 0.05) {
    cat("  -> WARNING: Statistically significant artifacts detected.\n")
  } else {
    cat("  -> SUCCESS: Indistinguishable from real data.\n")
  }

  # 3. Component Breakdown
  comps <- res$components
  worst_comp <- names(comps)[which.max(comps)]
  cat("\nComponent Errors (Lower is Better):\n")
  print(round(comps, 4))

  cat(paste0("\n  -> PRIMARY FAILURE MODE: ", worst_comp, "\n"))
  if (worst_comp == "Structural") {
    cat("     (The model captures counts but breaks the co-occurrence network.)\n")
  } else if (worst_comp == "Univariate") {
    cat("     (The model fails to capture basic taxon abundance counts.)\n")
  } else {
    cat("     (The model fails to capture the joint sample geometry.)\n")
  }
  cat("----------------------------------------\n")
}

# --- 1. Load All Four Datasets ---
# Ensure these files exist in your 'data/' folder
load("data/amgut1.filt.rda")      # Real Data A
load("data/amgut1.filt.cs.rda")   # Real Data B (CS version)
load("data/synth_zinb.rda")       # Model A
load("data/synth_htln.rda")       # Model B

# --- 2. Create Audit Objects ---
obj_filt_zinb <- create_asmot(real = amgut1.filt, synth = synth_zinb)
obj_filt_htln <- create_asmot(real = amgut1.filt, synth = synth_htln)
obj_cs_zinb   <- create_asmot(real = amgut1.filt.cs, synth = synth_zinb)
obj_cs_htln   <- create_asmot(real = amgut1.filt.cs, synth = synth_htln)

# --- 3. Run Audits (B=2000 is rigorous, B=50 is fast) ---
B_reps <- 2000
fixed_weights <- c(1, 1, 1) # Equal weights for fair comparison

message("Running Audit 1/4: Standard vs ZINB...")
res_filt_zinb <- asmot_audit(obj_filt_zinb, B = B_reps, weights = fixed_weights)

message("Running Audit 2/4: Standard vs HTLN...")
res_filt_htln <- asmot_audit(obj_filt_htln, B = B_reps, weights = fixed_weights)

message("Running Audit 3/4: CS vs ZINB...")
res_cs_zinb   <- asmot_audit(obj_cs_zinb, B = B_reps, weights = fixed_weights)

message("Running Audit 4/4: CS vs HTLN...")
res_cs_htln   <- asmot_audit(obj_cs_htln, B = B_reps, weights = fixed_weights)

# --- 4. PRINT INTERPRETATIONS ---
interpret_asmot(res_filt_zinb, "Standard Data vs ZINB")
interpret_asmot(res_filt_htln, "Standard Data vs HTLN")
interpret_asmot(res_cs_zinb,   "CS Data vs ZINB")
interpret_asmot(res_cs_htln,   "CS Data vs HTLN")

# --- 5. Summary Table & Bar Chart ---
results_df <- data.frame(
  Ground_Truth = c(rep("Standard (Filt)", 2), rep("CS (Filt.CS)", 2)),
  Model        = c("ZINB", "HTLN", "ZINB", "HTLN"),
  Realism      = c(res_filt_zinb$realism_score, res_filt_htln$realism_score,
                   res_cs_zinb$realism_score,   res_cs_htln$realism_score)
)
print(results_df)

p1 <- ggplot(results_df, aes(x = Ground_Truth, y = Realism, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Realism, 1)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 110) +
  labs(title = "Benchmark Results: ZINB vs HTLN",
       subtitle = "Realism Score (Higher is Better)",
       x = "Real Data Version Used",
       y = "ASMOT Realism Score (%)") +
  theme_minimal() +
  theme(text = element_text(size = 14))

print(p1)

# --- 6. HIGH-ORDER SCAN (k-Scan) for ALL 4 SCENARIOS ---
message("\n--- Running High-Order k-Scans (Dimensions 1 to 10) ---")

# We create a helper function to run and extract the plot
run_scan_plot <- function(obj, label) {
  # This runs the scan from k=1 to k=10
  scan_result <- asmot_plot_scan(obj, max_k = 10)

  # Add a specific title so we know which scenario it is
  p <- scan_result$plot +
       ggtitle(paste("k-Scan:", label)) +
       theme(plot.title = element_text(size=10))
  return(p)
}

message("Scanning 1/4 (Std vs ZINB)...")
plot_scan_1 <- run_scan_plot(obj_filt_zinb, "Standard vs ZINB")

message("Scanning 2/4 (Std vs HTLN)...")
plot_scan_2 <- run_scan_plot(obj_filt_htln, "Standard vs HTLN")

message("Scanning 3/4 (CS vs ZINB)...")
plot_scan_3 <- run_scan_plot(obj_cs_zinb, "CS vs ZINB")

message("Scanning 4/4 (CS vs HTLN)...")
plot_scan_4 <- run_scan_plot(obj_cs_htln, "CS vs HTLN")

# Arrange the 4 scan plots into one grid
message("Generating Consolidated Scan Plot...")
# Use grid.arrange to put them in a 2x2 grid
grid.arrange(plot_scan_1, plot_scan_2, plot_scan_3, plot_scan_4, ncol = 2)

# --- 7. ADVERSARIAL CLASSIFIER (Turing Test) ---
message("\n--- Training Adversarial Classifiers ---")

rf_filt_zinb <- asmot_classifier(obj_filt_zinb)
rf_filt_htln <- asmot_classifier(obj_filt_htln)
rf_cs_zinb   <- asmot_classifier(obj_cs_zinb)
rf_cs_htln   <- asmot_classifier(obj_cs_htln)

accuracy_df <- data.frame(
  Scenario = c("Standard vs ZINB", "Standard vs HTLN",
               "CS vs ZINB",       "CS vs HTLN"),
  Accuracy = c(rf_filt_zinb$Accuracy, rf_filt_htln$Accuracy,
               rf_cs_zinb$Accuracy,   rf_cs_htln$Accuracy)
)

print("\n=== ADVERSARIAL CLASSIFIER ACCURACY ===")
print("0.50 = Perfect (Indistinguishable)")
print("1.00 = Fail (Easily Spotted)")
print(accuracy_df)

# --- 8. Visualization: Classifier Distributions ---
combine_data <- function(res, name) {
  df <- res$Data
  df$Scenario <- name
  return(df)
}

plot_data <- rbind(
  combine_data(rf_filt_zinb, "Standard vs ZINB"),
  combine_data(rf_filt_htln, "Standard vs HTLN"),
  combine_data(rf_cs_zinb,   "CS vs ZINB"),
  combine_data(rf_cs_htln,   "CS vs HTLN")
)

p2 <- ggplot(plot_data, aes(x = OT_Distance, fill = Label)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Scenario, scales = "free") +
  scale_fill_manual(values = c("Real" = "black", "Synthetic" = "#E69F00")) +
  labs(title = "The Turing Test: Can the Classifier Spot the Fakes?",
       subtitle = "Overlapping curves = Better Synthesis (Lower Accuracy)",
       x = "Distance to Real Centroid (Geometry)",
       y = "Density") +
  theme_minimal()

print(p2)
