#' @title Explicit Metagenomic Preprocessing
#' @description 
#' Helper function to load and transform count data. ASMOT requires explicit preprocessing
#' to ensure transparency (no hidden transforms).
#'
#' @param input Either a file path (CSV) or a numeric matrix/data.frame.
#' @param method Character. One of:
#' \itemize{
#'   \item \code{"tss"}: Total Sum Scaling (Relative Abundance).
#'   \item \code{"clr"}: Centered Log-Ratio (for correlations).
#'   \item \code{"hellinger"}: Hellinger transformation.
#'   \item \code{"log10"}: Log-10 transformation.
#' }
#' @param pseudocount Numeric. Added to zeros for log-based methods (default 1e-6).
#' @param reference_col Integer. Column index for ALR transformation (if implemented).
#'
#' @return A numeric matrix of transformed data.
#' @export
asmot_preprocess <- function(input, method = "tss", pseudocount = 1e-6, reference_col = 1) {
  # (Implementation as provided in previous steps)
  if (is.character(input) && file.exists(input)) {
    data <- read.csv(input, row.names = 1, check.names = FALSE)
  } else {
    data <- input
  }
  mat <- as.matrix(data)
  
  if (any(mat == 0) && method %in% c("clr", "alr", "log10")) {
    mat <- mat + pseudocount
  }
  
  result <- switch(method,
    "tss" = t(apply(mat, 1, function(x) x / sum(x))),
    "clr" = t(apply(mat, 1, function(x) {
      g_mean <- exp(mean(log(x)))
      log(x / g_mean)
    })),
    "hellinger" = t(apply(mat, 1, function(x) sqrt(x / sum(x)))),
    "log10" = log10(mat),
    "alr" = {
      ref <- mat[, reference_col]
      t(apply(mat, 1, function(x) log(x / ref)))
    },
    stop("Unknown method.")
  )
  return(result)
}
