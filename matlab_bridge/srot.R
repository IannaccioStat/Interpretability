# =============================================================================
# calc_rotation.R
# =============================================================================
# DESCRIPTION:
#   MATLAB-to-R Bridge Script for Oblique Factor Rotations using GPArotation.
#   Executes specialized rotation criteria via Gradient Projection Algorithms 
#   with multi-start random initializations to find global minima.
#
# INTERACTION PROTOCOL:
#   Input File  : 'temp_A_norm.mat' (Contains A_norm, method_str, num_starts)
#   Output File : 'temp_rot_out.mat' (Exports Lambda_norm, T, f_opt)
#
# REFERENCES:
#
#   1. GPArotation Package & Framework:
#      Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms 
#      and software for arbitrary rotation criteria in factor analysis. 
#      Educational and Psychological Measurement, 65(5), 761–796. 
#      https://doi.org/10.1177/0013164404272507
#
#   2. Bentler Invariant Criterion (bentler):
#      Bentler, P. M. (1977). Factor simplicity index and transformations. 
#      Psychometrika, 42(2), 277–295. 
#      https://doi.org/10.1007/BF02294054
#
#   3. Geomin Criterion (geomin):
#      Yates, A. (1987). Multivariate Exploratory Data Analysis: A Perspective 
#      on Exploratory Factor Analysis. State University of New York Press.
#
#   4. Infomax Criterion (infomax):
#      Asuncion, A., & Jennrich, R. I. (2002). Infomax rotation in factor analysis. 
#      Psychometrika, 67(1), 111–120. 
#      https://doi.org/10.1007/BF02294712
#
#   5. McCammon Minimum Entropy Criterion (mccammon):
#      McCammon, R. B. (1966). Principal component analysis and its application 
#      in limnology. Journal of Geology, 74(5), 721–726. 
#      https://doi.org/10.1086/627202
# =============================================================================
set.seed(123)
suppressPackageStartupMessages({
  suppressWarnings({
    library(GPArotation)
    library(R.matlab)
  })
})

input_file  <- "temp_A_norm.mat"
output_file <- "temp_rot_out.mat"

if (file.exists(input_file)) {
  data_struct <- readMat(input_file)
  
  A_norm     <- as.matrix(data_struct$A.norm)
  num_starts <- as.numeric(data_struct$num.starts)
  
  # Flatten and clean string input from MATLAB
  raw_method <- data_struct$method.str
  method_str <- toupper(trimws(paste(as.character(unlist(raw_method)), collapse = "")))
  
  res <- NULL
  
  tryCatch({
    if (method_str == "BENTLER") {
      res <- bentlerQ(A_norm, randomStarts = num_starts, normalize = FALSE)
      
    } else if (method_str == "GEOMIN") {
      res <- geominQ(A_norm, randomStarts = num_starts, normalize = FALSE)
      
    } else if (method_str == "INFOMAX") {
      res <- infomaxQ(A_norm, randomStarts = num_starts, normalize = FALSE)
      
    } else if (method_str == "MCCAMMON") {
      res <- mccammon(A_norm, randomStarts = num_starts, normalize = FALSE)
      
    } else {
      # Fallback for standard built-in functions in GPArotation
      res <- GPFoblq(A_norm, method = tolower(method_str), normalize = FALSE)
    }
  }, error = function(e) {
    cat("[R-Bridge Rotation] ERROR during GPA rotation execution:\n")
    cat(sprintf("          %s\n", e$message))
  })
  
  # Export results back to MATLAB
  if (!is.null(res)) {
    writeMat(output_file,
             Lambda_norm = as.matrix(res$loadings),
             T           = as.matrix(res$Th),
             f_opt       = as.numeric(res$f))
  }
}