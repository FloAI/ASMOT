# --- 1. Simulation Setup ---
set.seed(42)
# library(ASMOT) 

# Generate Dummy Data
n_taxa <- 20; n_samples <- 50
real_mat <- abs(matrix(rnorm(n_samples*n_taxa, mean=10, sd=2), nrow=n_samples))
synth_mat <- abs(matrix(rnorm(n_samples*n_taxa, mean=10, sd=2), nrow=n_samples))
colnames(real_mat) <- colnames(synth_mat) <- paste0("Taxon_", 1:n_taxa)

# --- 2. Initialize ---
obj <- create_asmot(real = real_mat, synth = synth_mat)

# --- 3. Run Audit (Adaptive Weights) ---
res <- asmot_audit(obj, B = 50, weights = "auto")
print(res)

# --- 4. High-Order Scan ---
scan_res <- asmot_plot_scan(obj, max_k = 5)
print(scan_res)

# --- 5. Dashboard ---
dash <- asmot_plot_dashboard(res)
print(dash)
