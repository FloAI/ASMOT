#' Real Metagenomic Data (Standard Filtered)
#'
#' A matrix of real microbiome counts from the American Gut Project (filtered).
#' Rows are samples, columns are taxa.
#'
#' @format A numeric matrix with 50 rows and 20 columns.
#' @source Simulated/Derived for testing.
#' @usage data(amgut1.filt)
"amgut1.filt"

#' Real Metagenomic Data (CS Filtered)
#'
#' A variant of the real microbiome counts processed with Compressive Sensing (CS).
#' Rows are samples, columns are taxa.
#'
#' @format A numeric matrix with 50 rows and 20 columns.
#' @usage data(amgut1.filt.cs)
"amgut1.filt.cs"

#' Synthetic Data (ZINB Model)
#'
#' Synthetic data generated using a Zero-Inflated Negative Binomial (ZINB) model.
#' Intended as a baseline for benchmarking.
#'
#' @format A numeric matrix with 50 rows and 20 columns.
#' @usage data(synth_zinb)
"synth_zinb"

#' Synthetic Data (HTLN Model)
#'
#' Synthetic data generated using a Hurdle-Truncated Log-Normal(HTLN) model.
#' Intended for benchmarking against ZINB.
#'
#' @format A numeric matrix with 50 rows and 20 columns.
#' @usage data(synth_htln)
"synth_htln"
