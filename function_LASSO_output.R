function_LASSO_output =  function(seedVal, X_train, y_train, X_test, y_test){
  set.seed(seedVal)
  start_time <- Sys.time()
  lasso_cv <- cv.glmnet(as.matrix(X_train), y_train, alpha = 1)
  elapsed_time <- as.numeric( Sys.time() - start_time, units="secs" )
  
  optimal_lambda <- lasso_cv$lambda.min # plot(cv_lasso)
  
  beta_estimate <- coef(lasso_cv, s = optimal_lambda)
  lasso_selected_variables = which(beta_estimate[-1] != 0)

  
  
  y_pred_test <- as.numeric( predict(lasso_cv, newx = as.matrix(X_test), s = optimal_lambda))
  


  LASSO_mse <- mean(( y_pred_test  - y_test)^2)
  LASSO_r_sq_test = cor(y_pred_test, y_test)^2
  
  
  return(list(mse = LASSO_mse, r_sq_test = LASSO_r_sq_test, elapsed_time = elapsed_time ,
              optimal_lambda=  optimal_lambda,
              beta_with_intercept = beta_estimate, selected_variables = lasso_selected_variables, 
              method_used = "LASSO"))
  
}


function_OLS_output <- function(X_train, y_train, X_test, y_test) {
  
  
  start_time <- Sys.time()
  ols_fit <- lm( y_train ~ ., data = as.data.frame(X_train) )
  elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
  
  
  ols_beta_with_intercept =  ols_fit$coefficients
  ols_selected_variables = which(ols_beta_with_intercept[-1] != 0)
  names(ols_selected_variables) = NULL
  
  
  
  y_pred_test <- as.numeric(predict(ols_fit, newdata = as.data.frame(X_test)))
  
  
  
  ols_mse <- mean((y_pred_test - y_test)^2)
  r_sq_test <- cor(y_pred_test, y_test)^2
  
  return(list(
    mse = ols_mse,    r_sq_test = r_sq_test,    elapsed_time = elapsed_time,
    beta_with_intercept = ols_fit$coefficients, selected_variables = ols_selected_variables,
    method_used = "OLS"
  ))
}



function_OLS_LASSO_output <- function(seedVal,X_train, y_train,
                                      X_test, y_test) {
  
  
  
  if ( ncol(X_train) < (nrow(X_train) / 2)) {
    out <- function_OLS_output(  X_train, y_train, X_test, y_test)
  } else {
    out <- function_LASSO_output( seedVal , X_train, y_train , X_test,  y_test)
  }
  
  
  return(out)
}

