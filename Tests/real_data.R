library(ASMOT)
library(ggplot2)

set.seed(123) # Ensure reproducibility

# --- Helper Function: Interpret Results ---
# This prints the specific "Metric Result Interpretation" you asked for
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

  # 3. Component Breakdown (Why did it fail?)
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
load("data/amgut1.filt.rda")      # Real Data A
load("data/amgut1.filt.cs.rda")   # Real Data B (CS version)
load("data/synth_zinb.rda")       # Model A
load("data/synth_htln.rda")       # Model B

# --- 2. Create Audit Objects ---
obj_filt_zinb <- create_asmot(real = amgut1.filt, synth = synth_zinb)
obj_filt_htln <- create_asmot(real = amgut1.filt, synth = synth_htln)
obj_cs_zinb   <- create_asmot(real = amgut1.filt.cs, synth = synth_zinb)
obj_cs_htln   <- create_asmot(real = amgut1.filt.cs, synth = synth_htln)

# --- 3. Run Audits (B=50 for speed, use 100+ for paper) ---
B_reps <- 50
fixed_weights <- c(0.2, 0.4, 0.4)
message("Running Audit 1/4: Standard vs ZINB...")
res_filt_zinb <- asmot_audit(obj_filt_zinb, B = B_reps, weights = fixed_weights)

message("Running Audit 2/4: Standard vs HTLN...")
res_filt_htln <- asmot_audit(obj_filt_htln, B = B_reps, weights = fixed_weights)

message("Running Audit 3/4: CS vs ZINB...")
res_cs_zinb   <- asmot_audit(obj_cs_zinb, B = B_reps, weights = fixed_weights)

message("Running Audit 4/4: CS vs HTLN...")
res_cs_htln   <- asmot_audit(obj_cs_htln, B = B_reps, weights = fixed_weights)

# --- 4. PRINT INTERPRETATIONS ---
# This triggers the detailed text output for each
interpret_asmot(res_filt_zinb, "Standard Data vs ZINB")
interpret_asmot(res_filt_htln, "Standard Data vs HTLN")
interpret_asmot(res_cs_zinb,   "CS Data vs ZINB")
interpret_asmot(res_cs_htln,   "CS Data vs HTLN")

# --- 5. Summary Table ---
results_df <- data.frame(
  Ground_Truth = c(rep("Standard (Filt)", 2), rep("CS (Filt.CS)", 2)),
  Model        = c("ZINB", "HTLN", "ZINB", "HTLN"),
  Realism      = c(res_filt_zinb$realism_score, res_filt_htln$realism_score,
                   res_cs_zinb$realism_score,   res_cs_htln$realism_score)
)
print(results_df)

# --- 6. Visualization: Grouped Bar Chart ---
# Example Scan (High-order interactions)
message("Running k-Scan for CS vs HTLN...")
scan_res <- asmot_plot_scan(obj_cs_htln, max_k = 5)
print(scan_res$plot)

p1 <- ggplot(results_df, aes(x = Ground_Truth, y = Realism, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Realism, 1)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 110) +
  labs(title = "Comprehensive Benchmark: ZINB vs HTLN",
       subtitle = "Comparing Geometric Fidelity against two Ground Truths",
       x = "Real Data Version Used",
       y = "ASMOT Realism Score (%)") +
  theme_minimal() +
  theme(text = element_text(size = 14))

print(p1)
