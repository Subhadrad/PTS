function_iboss_subsample_index <- function(seedVal, X_input, subsample_size) {
  set.seed(seedVal)
  start_time <- Sys.time()
  # X_input = Method_TH_X_train_reduced 
  # seed_parallel_runs = 1
  
  # dim(Method_TH_X_train_reduced)

  n_input = nrow(X_input)
  p_input = ncol(X_input)
  
  allocations_per_coordinate <- ceiling(subsample_size/(2*p_input)) 
  subsample_index_iboss <- NULL
  
  # Now pick the points with extreme coordinates
  for(i in 1:p_input) 
  { # i =1 i=2
    subsample_index <- NULL
    
    subsample_temp_lo <-   setdiff(order(X_input[, i]), subsample_index_iboss)[1:allocations_per_coordinate]     
    subsample_temp_hi <-   setdiff( order(X_input[, i], decreasing = TRUE), subsample_index_iboss)[1:allocations_per_coordinate]     
    subsample_index  <- unique(c(subsample_temp_lo, subsample_temp_hi)) 

    subsample_index_iboss <- c(subsample_index_iboss, subsample_index)
  }
  
  # subsample_size - length(subsample_index_iboss) 
  # NOTE: Due to the ceiling rounding there could be slightly more subsamples than requested.
  # In order to keep the comparison fair with other subsamples, we choose indices from subsample_iboss_index randomly   
  
  
  if (length(subsample_index_iboss) > subsample_size) {
    subsample_index_iboss <- sample(subsample_index_iboss, subsample_size, replace = FALSE)
    }
  
  elapsed_time <- as.numeric( Sys.time() - start_time, units="secs" )
  
  return( list(subsample_index_iboss = subsample_index_iboss, elapsed_time= elapsed_time) )

}


function_uniform_subsampling <- function(seedVal, X_train, y_train,  subsample_size){
 set.seed(seedVal)
 start_time <- Sys.time()
 subsample_index =  sample(1:nrow(X_train), subsample_size)
 elapsed_time <- as.numeric( Sys.time() - start_time, units="secs" )
      
return(list(subsample_index = subsample_index  , elapsed_time = elapsed_time) )

  
}


function_osmac_subsample_index_weights <- function(seedvalue, DataSet, X_train, y_train, subsample_size){ 

  k0_fold <- 1/5
  set.seed(seedvalue)
  
  subsample_size_I <- floor(subsample_size * k0_fold)
  subsample_size_II <- subsample_size - subsample_size_I

  if(DataSet %in% c("Synthetic and Sparse" , "Normal", "IBOSS_t2", "IBOSS_normal_vary_SNR"  )) 
     {  
       subsample_index_OSMAC_I <-   sample(1:nrow(X_train),  subsample_size_I, replace = TRUE)
       y_train_OSMAC_initial <-  y_train[subsample_index_OSMAC_I]
       X_train_OSMAC_initial <-  X_train[subsample_index_OSMAC_I, ]
       model_OSMAC_initial <- lm(y_train_OSMAC_initial ~ ., data = as.data.frame(X_train_OSMAC_initial))  
       numerator <- abs(y_train - predict(model_OSMAC_initial, newdata = as.data.frame(X_train))) * sqrt(apply(X_train^2,1,sum))
     }
  # if( DataSet ==  c("UK Biobank-Semisynthetic With Screening")){
  #   
  #   X_train_Scaled <- scale(X_train) 
  #   
  #   subsample_index_OSMAC_I <-   sample(1:nrow(X_train_Scaled),  subsample_size_I, replace = TRUE)
  #   y_train_OSMAC_initial <-  y_train[subsample_index_OSMAC_I]
  #   X_train_OSMAC_initial <-  X_train_Scaled[subsample_index_OSMAC_I, ]
  #   model_OSMAC_initial <- lm(y_train_OSMAC_initial ~ ., data = as.data.frame(X_train_OSMAC_initial))  
  #   numerator <- abs(y_train - predict(model_OSMAC_initial, newdata = as.data.frame(X_train_Scaled))) * sqrt(apply(X_train_Scaled^2,1,sum))
  #   
  #     
  # }
  
  
  # boxplot(numerator)
  pi_i <- numerator/sum(numerator)
  weights <- ((1/k0_fold-1)*k0_fold)*pi_i+(k0_fold)*(1/nrow(X_train))
  
# return(list(pi_i= pi_i , weights = weights, subsample_index_OSMAC_I = subsample_index_OSMAC_I, subsample_size_II = subsample_size_II ))


subsample_index_OSMAC_II <- sample(1:nrow(X_train),subsample_size_II, prob = pi_i,replace = TRUE)



subsample_index  <-         c(subsample_index_OSMAC_I, subsample_index_OSMAC_II)
subsample_weight <- weights[c(subsample_index_OSMAC_I, subsample_index_OSMAC_II) ]
return(list(subsample_index = subsample_index, subsample_weight = subsample_weight  ))



# X_train_OSMAC <- X_train[c(subsample_index_OSMAC_I, subsample_index_OSMAC_II), ]
# y_train_OSMAC <- y_train[c(subsample_index_OSMAC_I, subsample_index_OSMAC_II) ]
# weight_OSMAC <- weights[c(subsample_index_OSMAC_I, subsample_index_OSMAC_II) ]
# return(list(X_train_OSMAC= X_train_OSMAC , y_train_OSMAC= y_train_OSMAC,
#             weight_OSMAC = weight_OSMAC  ))
  

}   



