library(pracma)
library(MASS)
library(BEDMatrix)
library(data.table)




function_generate_covariance_matrix <- function(p, covariate_covariance_structure, LD_rho){
  
  noiseMatrix <- outer(1:p, 1:p, function(i,j) 0.5^abs(i-j) )

  
  R <- randortho(p, type="orthonormal")
  
  if(covariate_covariance_structure == "Low-dimensional-on 2 PC" ){
    EigenValues <- rep(0, p)
    EigenValues[1] <- 100; EigenValues[2] <- 2; 
    temp = matrix(0, ncol = p, nrow = p)
    for (i in 1:2)
    {temp= temp + EigenValues[i]  * as.matrix(R[,i]) %*% t( as.matrix(R[,i]))}
    Sigma  <-       temp + .01 * noiseMatrix
    # eigen(temp)$values
    # eigen(design_matrix_covariance)$values  
  }else if(covariate_covariance_structure == "Equal_Corelation_0.5"){
    Sigma <- matrix(0.5, nrow = p, ncol = p)
    diag(Sigma) <- 1 
  }else if(covariate_covariance_structure == "AR1"){
    Sigma <- outer(1:p, 1:p, function(i,j) LD_rho^abs(i-j) )
  }else if(covariate_covariance_structure == "iid"){
    Sigma <- diag(p)
  }
  Sigma
}
function_generate_covariates_block_mvn <- function( n, LD_block_size, n_LD_blocks, covariate_covariance_structure, LD_rho) {
  
  # Checks
  stopifnot( n >= 3, LD_block_size >= 1, n_LD_blocks >= 1)
  stopifnot(LD_rho > -1 && LD_rho < 1)
  
  
  X_blocks <- vector("list", n_LD_blocks)
  Sigma <- function_generate_covariance_matrix(LD_block_size, covariate_covariance_structure, LD_rho)
  
  for (b in seq_len(n_LD_blocks)) {
    X_blocks[[b]] <- MASS::mvrnorm(n = n, mu = rep(0, LD_block_size), Sigma = Sigma)
    #  <- scale(Xb)  # standardize columns
  }
  
  
  X <- do.call(cbind, X_blocks)
  X
}
# p = 10
function_generate_beta <- function(beta_type, LD_block_size,  n_LD_blocks, sparsity, intercept) {
  p = LD_block_size* n_LD_blocks
  if(intercept == "Yes"){
    if(beta_type == "Sparse_Unif_minus1_1"){
      beta_0 <- runif( 1 , -1, 1)  
      beta_1 <- c( runif(  floor(sparsity/100* p) , -1, 1), rep(0, p-  floor(sparsity* p/100))  )    
    }
    if(beta_type == "fixed_1_notSparse"){
      beta_0 <- 1
      beta_1 <- rep(1, p)
    }
    if(beta_type == "fixed_1_Sparse"){
      beta_0 <- 1
      beta_1 <- c(rep(1, floor(sparsity/100* p)), rep(0, p-floor(sparsity/100* p)))
    }
    if(beta_type == "LD-1"){
      beta_0 <- 1
      beta_1 <- replace(numeric(p), c(1, 2, floor(p/2) , floor(p/2)+1, p ), c(1.1, 1.2, 0.8, 0.5, 2))
    
    }
    if(beta_type == "LD-2"){
      p = LD_block_size  
      beta_0 <- 1
      beta_1_c1 <- replace(numeric(p), c(1, 2, floor(p/2) , floor(p/2)+1, p ), c(1.1, 1.2, 0.8, 0.5, 2))
      beta_1_c2 <- replace(numeric(p), c(4, 5, floor(p/2) , floor(p/2)+1 ), c(-1.1, 1.2, 2, 0.2))
      beta_1 <- c(beta_1_c1, beta_1_c2)
    }
    if(beta_type == "LD-3"){
      p = LD_block_size  
      beta_0 <- 1
      beta_1_c1 <- replace(numeric(p), c(1, 2, floor(p/2) , floor(p/2)+1, p ), c(1.1, 1.2, 0.8, 0.5, 2))
      beta_1_c2 <- replace(numeric(p), c(4, 5, floor(p/2) , floor(p/2)+1 ), c(-1.1, 1.2, 2, 1))
      beta_1_c3 <- replace(numeric(p), c(1, floor(p/2)+1 ), c(1,  0.2))
      beta_1 <- c(beta_1_c1, beta_1_c2, beta_1_c3)
    }
    if(beta_type == "LD-4"){
      p = LD_block_size  
      beta_0 <- 1
      beta_1_c1 <- replace(numeric(p), c(1, 2, floor(p/2) , floor(p/2)+1, p ), c(1.1, 1.2, 0.8, 0.5, 2))
      beta_1_c2 <- replace(numeric(p), c(4, 5, floor(p/2) , floor(p/2)+1 ), c(-1.1, 1.2, 2, 1))
      beta_1_c3 <- replace(numeric(p), c(1, floor(p/2)+1 ), c(1,  0.2))
      beta_1_c4 <- replace(numeric(p), c( floor(p/2)   ), c( 2))
      beta_1 <- c(beta_1_c1, beta_1_c2, beta_1_c3, beta_1_c4)
    }
    # "LD-5" ####
    if(beta_type == "LD-5"){
      p = LD_block_size  
      beta_0 <- 1
      beta_1_c1 <- replace(numeric(p), c(1, 2, floor(p/2) , floor(p/2)+1, p ), c(1.1, 1.2, 0.8, 0.5, 2))
      beta_1_c2 <- replace(numeric(p), c(4, 5, floor(p/2) , floor(p/2)+1 ), c(-1.1, 1.2, 2, 1))
      beta_1_c3 <- replace(numeric(p), c(1, floor(p/2)+1 ), c(1,  0.2))
      beta_1_c4 <- replace(numeric(p), c( floor(p/2)   ), c( 2))
      beta_1_c5 <- numeric(p)
      beta_1 <- c(beta_1_c1, beta_1_c2, beta_1_c3, beta_1_c4, beta_1_c5)
      
      beta_1[beta_1!=0]
    }
    
    
    # "LD-5 sparsity Type-I" ####
    if(beta_type == "LD-5 sparsity Type-I"){
      p = LD_block_size* n_LD_blocks
      
      beta_0 <- 1
      beta_1 <- numeric(p)
      
      # total nonzeros
      k_total <- floor(sparsity / 100 * p)
      signal_sites <- sample(1:p, k_total, replace = FALSE)
      beta_1[signal_sites] <- sample(c(-1, 1), length(signal_sites), replace = TRUE)
    }
    
    # "LD-5 sparsity Type-II"
    if(beta_type == "LD-5 sparsity"){
      beta_0 <- 1
      weights = c(40, 30, 20, 10, 0)
      adjacent_frac = 0.5
      stopifnot(p %% 5 == 0)
      m <- p / 5
      beta_1 <- numeric(p)
      
      # total nonzeros
      k_total <- round(sparsity / 100 * p)
      
      # allocate to blocks by weights, then fix rounding so sums match
      w <- weights / sum(weights)
      k_block <- floor(k_total * w)
      
      # # give leftover 1s to blocks with largest fractional parts
      # remainder <- k_total - sum(k_block)
      # if (remainder > 0) {
      #   frac <- k_total * w - floor(k_total * w)
      #   add_to <- order(frac, decreasing = TRUE)[1:remainder]
      #   k_block[add_to] <- k_block[add_to] + 1
      # }
      
      for (b in 1:5) {
        # b =1 
        k <- k_block[b]
        if (k <= 0) next
        
        k_adj <- round(adjacent_frac * k)
        k_far <- k - k_adj
        
        idx <- ((b - 1) * m + 1):(b * m)   # global indices for this block
        chosen <- integer(0)
        
        # adjacent run (one chunk)
        if (k_adj > 0) {
          
          # chunk sizes: 10,10,..., remainder
          n_full <- k_adj %/% 10
          rem <- k_adj %% 10
          sizes <- c(rep(10, n_full), if (rem > 0) rem)
          
          # deterministic starts with gap = k_adj: 1, 1+k_adj, 1+2*k_adj, ...
          starts <- 1 + 2*max(sizes)  * (0:(length(sizes) - 1))
          
          # safety check: must fit in the block of length m
          ends <- starts + sizes - 1
          if (max(ends) > m) stop("Adjacent chunks do not fit in block; need larger block size m or smaller k_adj.")
          
          # collect local positions for all chunks
          adj_local <- integer(0)
          for (i in seq_along(sizes)) {
            adj_local <- c(adj_local, starts[i]:(starts[i] + sizes[i] - 1))
          }
          
          # map to global indices in beta_1
          chosen <- idx[adj_local]
        }
        
        # far picks: avoid neighbors of chosen if possible
        if (k_far > 0) {
          chosen_local <- if (length(chosen) > 0) chosen - min(idx) + 1 else integer(0)
          forbidden <- unique(c(chosen_local - 1, chosen_local, chosen_local + 1))
          forbidden <- forbidden[forbidden >= 1 & forbidden <= m]
          candidates <- setdiff(1:m, forbidden)
          
          # if not enough candidates, relax and just avoid already-chosen
          if (length(candidates) < k_far) {
            candidates <- setdiff(1:m, chosen_local)
          }
          
          far_local <- sample(candidates, k_far, replace = FALSE)
          chosen <- c(chosen, idx[far_local])
        }
        
        beta_1[chosen] <- sample(c(-1, 1), length(chosen), replace = TRUE)
      }
      
      # sanity check: exact sparsity enforced
      # stopifnot(sum(beta_1 != 0) == k_total)
      
    }
    
    # "HAPNEST LD-5 sparsity Type-I" ####
    
    if(beta_type == "HAPNEST LD-5 sparsity Type-I"){
      p = LD_block_size* n_LD_blocks
      
      beta_0 <- 1
      beta_1 <- numeric(p)
      
      # total nonzeros
      k_total <- floor(sparsity / 100 * p)
      signal_sites <- sample(1:p, k_total, replace = FALSE)
      beta_1[signal_sites] <- sample(c(-1, 1), length(signal_sites), replace = TRUE)
    }
    
    
    
    }
    
    
    
    
    # length(which(beta_1 != 0 ))    
  
  list(beta_1 = beta_1, beta_0 = beta_0, beta = c(beta_0, beta_1))
  
}
function_compute_sigma_sq <- function(beta, sigma_blocks, n_LD_blocks, SNR, intercept_present_in_beta   ) {
if(intercept_present_in_beta == "Yes" )
  {beta = beta[-1]}

  # beta = beta_1 
  LD_block_size <- nrow(sigma_blocks)
  p <- LD_block_size * n_LD_blocks
  stopifnot(length(beta) == p)
  
  # length(Data$beta_1)
  # beta = Data$beta_1
  
  # reshape beta into block_size x n_LD_blocks
  beta_mat <- matrix(beta, nrow = LD_block_size, ncol = n_LD_blocks)
  
  signal_var <- sum(vapply(seq_len(n_LD_blocks), function(j) {
    bj <- beta_mat[, j, drop = FALSE]
    as.numeric(t(bj) %*% sigma_blocks %*% bj)
  }, numeric(1)))
  
  sigma_sq <- signal_var / SNR
  # list(signal_var = signal_var, sigma_sq = sigma_sq, sqrt_sigma_sq = sqrt(sigma_sq))
  sqrt(sigma_sq)
}


function_generate_data_version_1 <- function(seedvalue, DataSet_Name,  LD_block_size, n_LD_blocks, covariate_covariance_structure, LD_rho, n_train, n_test,   beta_type,  sparsity, SNR, noise_variance, intercept ){
  
  
  set.seed(seedvalue)
  # Simulated LD Structure ####
  if(DataSet_Name == "Simulated LD Structure"){
    
  p <- LD_block_size * n_LD_blocks  
  sigma_blocks <-  function_generate_covariance_matrix(LD_block_size, covariate_covariance_structure, LD_rho) 
  beta = function_generate_beta(beta_type, LD_block_size,  n_LD_blocks, sparsity, intercept)
  beta_0 = beta$beta_0
  beta_1 = beta$beta_1

  n = n_train + n_test

  
  if(noise_variance != "TBD")
  {
    sqrt_sigma_sq = sqrt(noise_variance)
  }else{
    sqrt_sigma_sq = function_compute_sigma_sq(beta_1, sigma_blocks, n_LD_blocks, SNR, "No")
  }
    
  
  

   
     
     
    
    X =  function_generate_covariates_block_mvn( n_train + n_test, LD_block_size, n_LD_blocks, covariate_covariance_structure, LD_rho)
    epsilon <- rnorm(n, mean = 0, sd = sqrt_sigma_sq)
    y <- beta_0 + X %*% beta_1 + epsilon
    
    
    

    X_train = X[1:n_train,]  
    X_test =  X[ (n_train+1): (n_train+n_test),]
    
    y_train = y[1:n_train]
    y_test = y[ (n_train+1): (n_train+n_test)]
      
    return(list(X_train = X_train, y_train = y_train,
                X_test = X_test, y_test = y_test,
                beta_0 = beta_0, beta_1 = beta_1, beta =  c(beta_0, beta_1)))  

  } 


  
}

