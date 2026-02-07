#' @import methods
NULL

# --- Class Definition ---
#' @title ASMOT Object
#' @description Stores real and synthetic metagenomic data for auditing.
#' @slot real_counts Matrix of raw counts (Samples x Taxa)
#' @slot synth_counts Matrix of raw counts (Samples x Taxa)
#' @slot real_ra Matrix of relative abundances (Samples x Taxa)
#' @slot synth_ra Matrix of relative abundances (Samples x Taxa)
#' @slot cost_matrix Matrix defining taxon-taxon distances
#' @slot taxa_names Character vector of taxon names
#' @export
setClass("ASMOT",
         slots = list(
           real_counts = "matrix",
           synth_counts = "matrix",
           real_ra = "matrix",
           synth_ra = "matrix",
           cost_matrix = "matrix",
           taxa_names = "character"
         ))

#' @title Initialize ASMOT Audit
#' @export
create_asmot <- function(real, synth, tree = NULL) {
  common <- intersect(colnames(real), colnames(synth))
  if(length(common) == 0) stop("No common taxa found between datasets.")

  R_counts <- as.matrix(real[, common])
  S_counts <- as.matrix(synth[, common])

  # Compute Relative Abundances (TSS)
  R_ra <- t(apply(R_counts, 1, function(x) if(sum(x)>0) x/sum(x) else x))
  S_ra <- t(apply(S_counts, 1, function(x) if(sum(x)>0) x/sum(x) else x))

  # Build Cost Matrix
  if (!is.null(tree)) {
    M <- phangorn::cophenetic.phylo(tree)
    M <- M[common, common]
  } else {
    M <- as.matrix(dist(diag(length(common))))
  }

  new("ASMOT",
      real_counts = R_counts, synth_counts = S_counts,
      real_ra = R_ra, synth_ra = S_ra,
      cost_matrix = M, taxa_names = common)
}

# --- S3 Print Methods ---

#' @export
print.asmot_audit_res <- function(x, ...) {
  cat("\n=== ASMOT Audit Results ===\n")
  cat(sprintf("Realism Score:   %.2f%%\n", x$realism_score))
  cat(sprintf("P-value:         %.4f\n", x$p_value))
  
  cat("\n--- Adaptive Weights Used ---\n")
  print(round(x$weights_used, 3))
  
  cat("\n--- Component Scores (Distance) ---\n")
  print(x$components)
  cat("\n")
}

#' @export
print.asmot_k_scan <- function(x, ...) {
  cat("\n=== ASMOT Interaction Scan (Plot Object) ===\n")
  cat(sprintf("Dimensions Scanned: %d to %d\n", min(x$data$k), max(x$data$k)))
  cat(sprintf("Failure Point: k = %d (Error: %.4f)\n", 
              x$data$k[which.max(x$data$Error)], max(x$data$Error)))
}

#' @export
print.asmot_k_result <- function(x, ...) {
  cat("\n=== Single K-Dimension Result ===\n")
  cat(sprintf("Dimension (k): %d\n", x$k))
  cat(sprintf("Mean Distance: %.4f\n", x$mean_distance))
}
