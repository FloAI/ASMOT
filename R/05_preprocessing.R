#' @title Preprocess Metagenomic Data
#' @description Reads a CSV or Matrix and applies standard metagenomic transformations.
#' @export
asmot_preprocess <- function(input, method = "tss", pseudocount = 1e-6, reference_col = 1) {
  
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
