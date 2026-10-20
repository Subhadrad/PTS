source("function_generate_data_version_1.R")
source("function_variable_reduction_methods_version2.R")
source("Function_Subsample_Tagging_Algorithms.R")


# PTS Wrapper function ####
function_select_PTS_subdata <- function(X, y,
                                   screening_methods = "TH", initial_prescreeing_random_subsample_percent,  
                                   threshold = 0.05, top_q = NULL, clumping_threshold = 0.1, 
                                   L = 10, ld_regions = NULL, SuSie_pval_cutoff = 0.05, 
                                   pip_cutoff = 0.5, include_CS = FALSE,
                                   subsampling_methods = "IBOSS",
                                   subsample_size, seed = 1) {
  
                            X <- as.matrix(X)
                            stopifnot(is.numeric(X), is.numeric(y), is.null(dim(y)), nrow(X) == length(y),
                                      nrow(X) >= 3, ncol(X) >= 1, all(is.finite(X)), all(is.finite(y)), var(y) > 0,
                                      length(subsample_size) == 1, subsample_size >= 3, subsample_size <= nrow(X),
                                      subsample_size == floor(subsample_size), length(threshold) == 1, threshold > 0, threshold < 1,
                                      length(clumping_threshold) == 1, clumping_threshold >= 0, clumping_threshold <= 1)
                            
                            
                            if (is.null(ld_regions)) ld_regions <- list(seq_len(ncol(X)))
                            stopifnot(is.list(ld_regions))
                            
                            # method =  screening_methods[3]
                            run_PTS_instance <- function(method) {
                              
                                                  set.seed(seed )
                                                  initial_prescreeing_random_subsample_indices <- sample.int( n = nrow(X),  size = floor(initial_prescreeing_random_subsample_percent / 100 * nrow(X)), replace = FALSE )

                                                  result <- switch(method,
                                                                   TH = function_variable_reduction_Thresholding(X, y, initial_prescreeing_random_subsample_indices,
                                                                                                                 pval_cutoff = threshold, top_q),
                                                                   CT = function_variable_reduction_C_T(         X, y, initial_prescreeing_random_subsample_indices,
                                                                                                                 pval_cutoff = threshold, corr_cutoff = clumping_threshold),
                                                                   SuSiE = function_variable_reduction_SuSi_v2(  X, y, initial_prescreeing_random_subsample_indices,
                                                                                                                 ld_regions,L, pval_cutoff = SuSie_pval_cutoff, pip_cutoff,include_CS ) )
                                                  
                                                                
                                                  
                                                  # Stage 1: prescreen based on initail random subsample
                                                  stage_I_prescreened_var <- result$selected_cols
                                                  if (!length(stage_I_prescreened_var)) stop("Stage 1 retained no covariates. Try less strict screening.")
                                                  X_stage_1 <- X[, stage_I_prescreened_var, drop = FALSE]
                              
                                                  # Select rows using the requested subsampling algorithm.
                                                  stage_II_subsample_indices <- function_iboss_subsample_index(seedVal = seed, X_input = X_stage_1, subsample_size = subsample_size)$subsample_index_iboss
                              
                          
                                                  X_subdata <- X[stage_II_subsample_indices, stage_I_prescreened_var, drop = FALSE]
                                                  y_subdata <- y[stage_II_subsample_indices]
                                                  
                                                  list(screening = method, subsampling = subsampling_methods,
                                                       stage_I_prescreened_var = stage_I_prescreened_var,
                                                       stage_II_subsample_indices = stage_II_subsample_indices,
                                                       X_subdata = X_subdata, y_subdata = y_subdata
                                                       )
                                                                  }
                            
                            # if (length(screening_methods) == 1) return(run_PTS_instance(screening_methods))
                            results <- lapply(screening_methods, function(method) tryCatch(run_PTS_instance(method),
                                                                                 error = function(e) 
                                                                                   list(screening = method, subsampling = subsampling_methods, error = conditionMessage(e))))
                            
                            results <- setNames(results, paste0(screening_methods, "_", subsampling_methods))
                                                                  }

#Generate example data ####
function_generate_example_data <- function(n = 1000, LD_rho = 0.7, seed = 1) {
  set.seed(seed)
  ld_regions <- list(1:10, 11:20, 21:40)
  X_train <- matrix(NA_real_, n, 40)
  
  for (region in ld_regions) {
    Sigma <- function_generate_covariance_matrix(length(region), "AR1", LD_rho)
    X_train[, region] <- MASS::mvrnorm(n, mu = rep(0, length(region)), Sigma = Sigma)
  }
  
  beta_true <- rep(0, 40)
  beta_true[c(2, 5, 8)] <- c(3, -2, 2)
  beta_true[c(12, 15, 18)] <- c(2, -3, 1.5)
  beta_true[c(23, 30, 37)] <- c(-2, 3, 2)
  beta0 <- 1
  y_train <- drop(beta0 + X_train %*% beta_true + rnorm(n))
  
  list(X_train = X_train, y_train = y_train, ld_regions = ld_regions,
       beta_true = beta_true, beta0 = beta0)
}

