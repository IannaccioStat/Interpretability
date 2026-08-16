# =============================================================================
# calc_polychoric.R
# =============================================================================
# DESCRIPTION:
#   MATLAB-to-R Bridge Script for Mixed / Polychoric Correlation Estimation.
#   This script is invoked automatically by the MATLAB routine `processData.m`
#   to compute polychoric, tetrachoric, or mixed (polyserial/pearson) 
#   correlation matrices from raw data vectors.
#
# INTERACTION PROTOCOL:
#   Input File  : 'temp_X_data.mat'  (Contains matrix X from MATLAB)
#   Output File : 'temp_R_poly.mat'  (Exports R_poly, num_ordinal, num_continuous)
#   Engine      : R psych package (psych::polychoric, psych::mixedCor)
#
# REFERENCES:
#
#   1. R psych Package:
#      Revelle, W. (2024). psych: Procedures for Psychological, Psychometric, 
#      and Personality Research. Northwestern University, Evanston, Illinois. 
#      R package version 2.4.3. https://CRAN.R-project.org/package=psych
#
#   2. Polychoric & Tetrachoric Estimation Theory:
#      Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation 
#      coefficient. Psychometrika, 44(4), 443–460. 
#      https://doi.org/10.1007/BF02296207
#
#   3. Polyserial Correlation (Mixed Continuous & Ordinal Data):
#      Olsson, U., Drasgow, F., & Dorans, N. J. (1982). The polyserial correlation 
#      coefficient. Psychometrika, 47(3), 337–347. 
#      https://doi.org/10.1007/BF02294164
#
#   4. Correlation Matrix Smoothing (Positive Semi-Definiteness):
#      Knol, D. L., & ten Berge, J. M. (1989). Least-squares approximation of an 
#      inproper correlation matrix by a proper one. Psychometrika, 54(1), 53–61. 
#      https://doi.org/10.1007/BF02294451
# =============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(psych)
    library(R.matlab)
  })
})

input_file  <- "temp_X_data.mat"
output_file <- "temp_R_poly.mat"

if (file.exists(input_file)) {
  
  # ---------------------------------------------------------------------------
  # 1. Load Data Exchange File from MATLAB
  # ---------------------------------------------------------------------------
  data_struct <- readMat(input_file)
  X <- as.matrix(data_struct$X)
  J <- ncol(X)
  X_df <- as.data.frame(X)
  
  ordinal_cols          <- c()
  continuous_cols       <- c()
  imputed_rounded_count <- 0
  
  # ---------------------------------------------------------------------------
  # 2. Variable Classification (Max Value Rule: <= 7 is Ordinal)
  # ---------------------------------------------------------------------------
  for (col_idx in seq_len(J)) {
    vals <- X_df[[col_idx]]
    vals <- vals[is.finite(vals)] # Filter out non-finite/NA for evaluation
    
    if (length(vals) == 0) {
      continuous_cols <- c(continuous_cols, col_idx)
      next
    }
    
    val_max <- max(vals)
    
    # RULE: If max <= 7, round floating/imputed decimals to nearest integer & treat as Ordinal.
    if (val_max <= 7) {
      if (any(vals %% 1 != 0)) {
        X_df[[col_idx]] <- round(X_df[[col_idx]])
        imputed_rounded_count <- imputed_rounded_count + 1
      }
      ordinal_cols <- c(ordinal_cols, col_idx)
    } else {
      continuous_cols <- c(continuous_cols, col_idx)
    }
  }
  
  # ---------------------------------------------------------------------------
  # 3. Console Reporting to MATLAB Output Stream
  # ---------------------------------------------------------------------------
  if (imputed_rounded_count > 0) {
    cat(sprintf("[R-Bridge] %d column(s) with continuous decimals detected on Likert items. Values rounded to nearest integer.\n", 
                imputed_rounded_count))
  }
  
  cat(sprintf("[R-Bridge] Final Classification: %d Ordinal/Likert item(s), %d Continuous item(s).\n", 
              length(ordinal_cols), length(continuous_cols)))
  
  # ---------------------------------------------------------------------------
  # 4. Correlation Matrix Estimation via psych::mixedCor / psych::polychoric
  # ---------------------------------------------------------------------------
  R_poly <- NULL
  
  tryCatch({
    if (length(continuous_cols) == 0) {
      # --- PATH A: PURE ORDINAL MATRIX ---
      # Explicitly convert columns to ordered factors for psych::polychoric
      X_factors <- as.data.frame(lapply(X_df, function(col) factor(col, ordered = TRUE)))
      poly_res  <- suppressWarnings(psych::polychoric(X_factors, smooth = TRUE, global = FALSE, na.rm = TRUE))
      R_poly    <- poly_res$rho
      
    } else if (length(ordinal_cols) == 0) {
      # --- PATH B: PURE CONTINUOUS MATRIX ---
      R_poly <- cor(X_df, use = "pairwise.complete.obs")
      
    } else {
      # --- PATH C: MIXED (ORDINAL + CONTINUOUS) MATRIX ---
      # Ensure ordinal columns are explicit integer vectors for mixedCor
      for (ord_idx in ordinal_cols) {
        X_df[[ord_idx]] <- as.integer(X_df[[ord_idx]])
      }
      mixed_res <- suppressWarnings(psych::mixedCor(data = X_df, c = continuous_cols, p = ordinal_cols, use = "pairwise"))
      R_poly    <- mixed_res$rho
    }
  }, error = function(e) {
    cat("[R-Bridge] ERROR during mixed/polychoric correlation computation:\n")
    cat(sprintf("          %s\n", e$message))
  })
  
  # ---------------------------------------------------------------------------
  # 5. Fallback Guardrail & MATLAB Export
  # ---------------------------------------------------------------------------
  if (is.null(R_poly) || !is.matrix(R_poly)) {
    cat("[R-Bridge] WARNING: Polychoric optimization failed to converge. Falling back to pairwise Pearson correlation.\n")
    R_poly  <- cor(X, use = "pairwise.complete.obs")
    num_ord <- 0
    num_con <- J
  } else {
    num_ord <- length(ordinal_cols)
    num_con <- length(continuous_cols)
  }
  
  # Write exchange MAT-file back to MATLAB directory
  writeMat(output_file, 
           R_poly         = as.matrix(R_poly), 
           num_ordinal    = as.numeric(num_ord), 
           num_continuous = as.numeric(num_con))
}