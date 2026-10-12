library(glmnet)
library(MASS)
library(doParallel)
library(doRNG)
library(foreach)
library(bigsnpr)
library(dplyr)
library(pracma)


num_cores <- 5
num_simulations <- 100

cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1)
set.seed(1)

source("function_generate_data_version_1.R")
source("function_variable_reduction_methods_version2.R")
source("Function_Subsample_Tagging_Algorithms.R")
source("function_LASSO_output.R")
all_results <- NULL
all_summary <- NULL

# Simulation Parameters  ####
{
  series =  "Series 50c 1000by5 rhopt95 cases L10"
  LD_block_size = 1000
  n_LD_blocks = 5
  
  sparsity_1 = c(1) # Dummy value
  sparsity_3 = c(1)

  DataSet_Name = "Simulated LD Structure"
  covariate_covariance_structure = "AR1"
  noise_variance = "TBD"
  intercept = "Yes"
  
  n_train_list = c( 100*1000 )
  n_train = 100*1000
  LD_rho_list = c(  0.95 )
  
  beta_type_list_1 = c("LD-5")
  beta_type_list_3 = c( "LD-5 sparsity Type-I") 
  subsample_percent_variable_selection_list = c(50 )
  subsample_percent_list =  c(1.25/2, 1.25,   2.5, 5)
  SNR_list_1 = c(0.5)

  
}



params_1 = expand.grid(n_train_list, LD_rho_list,  SNR_list_1, beta_type_list_1, sparsity_1 )
params_3 = expand.grid(n_train_list, LD_rho_list,  SNR_list_1, beta_type_list_3, sparsity_3 )


params = rbind(params_1,  params_3)
dim(params)
nrow(params) * length(subsample_percent_variable_selection_list) * length(subsample_percent_list) * num_simulations





# The Simulation function ####
function_run_ensembles_folded <- function(sim_id,   n_train, LD_rho, SNR, beta_type, 
                                   subsample_percent_variable_selection_list, 
                                   subsample_percent_list ) {
  # sim_id = 1
  # sim_id = 74
  
  # subsample_percent_variable_selection = 50
  # subsample_percent_variable_selection = 100
  
  result = NULL
  ## STEP 1 ####
  ## a. Get or simulate the full data (both the test and train sets) ####
  ## b. Train Lasso full model ####
  ## c. Prescreen using the perfect variable selection ####
  
  if(  DataSet_Name %in% c("Simulated LD Structure")  ){
    n_test =  0.30 * n_train
    
    Data <- function_generate_data_version_1(sim_id, DataSet_Name,  LD_block_size, n_LD_blocks, covariate_covariance_structure, LD_rho, n_train, n_test,   beta_type,  sparsity, SNR, noise_variance, intercept )
    
    # Data <- function_generate_data_version_1(seedvalue = sim_id, DataSet_Name = DataSet_Name, LD_block_size = LD_block_size, n_LD_blocks = n_LD_blocks, covariate_covariance_structure = covariate_covariance_structure, LD_rho = LD_rho, n_train = n_train, n_test = n_test, beta_type = beta_type, sparsity = sparsity, SNR = SNR, noise_variance = noise_variance, intercept = intercept)
      
    X_train <- Data$X_train
    y_train <- Data$y_train
    
    X_test <- Data$X_test 
    y_test <- Data$y_test
    
    beta_true <- Data$beta
    beta_one_true <- Data$beta_1
    

    sparsity_corrected = ifelse(beta_type %in% c("LD-1","LD-2","LD-3","LD-4","LD-5"),
                                sum( beta_one_true!= 0)/(n_LD_blocks * LD_block_size) *100,
                                      sparsity  )
      
    
  }
  LASSO_Fullsample  =   function_LASSO_output(sim_id, X_train, y_train, X_test, y_test)
  
  
  
  Final_variable_selection_fullsample_LASSO =  function_variable_reduction_evaluate_selection( LASSO_Fullsample$selected_variables, Data$beta_1, eps = 1e-12) 
  
  
  Method_Perfect_variable_selection =  function_variable_reduction_perfect(X_train, Data$beta_1 )
  Method_Perfect_evaluate_variable_selection <- function_variable_reduction_evaluate_selection( Method_Perfect_variable_selection$selected_cols, 
                                                                                                 Data$beta_1, eps = 1e-12)
  
  
  
  ## STEP 2 Loop over list of c settings ####
  ## a. Pre-screened using TH, CT, and SuSiE  ####
      for (c_setting in subsample_percent_variable_selection_list) {
        
        # c_setting = subsample_percent_variable_selection_list[1]
        subsample_percent_variable_selection = c_setting 
        subsample_percent_variable_selection_indices =  1: floor(subsample_percent_variable_selection *  n_train/100)   
        
        Method_TH_variable_selection =  function_variable_reduction_Thresholding(X_train,  y_train, subsample_percent_variable_selection_indices,  pval_cutoff = 0.05)
        Method_CT_variable_selection =  function_variable_reduction_C_T(X_train, y_train, subsample_percent_variable_selection_indices, pval_cutoff = 0.05, corr_cutoff = 0.1)
        Method_SuSi_variable_selection =  function_variable_reduction_SuSi(X_train, y_train, subsample_percent_variable_selection_indices, LD_block_size, n_LD_blocks)
        
        # Method_Perfect_variable_selection_X_reduced_100 =  Method_Perfect_variable_selection$X_reduced
        # Method_Perfect_variable_selection_X_reduced_50 =  Method_Perfect_variable_selection$X_reduced
        # sum(Method_Perfect_variable_selection_X_reduced_100 != Method_Perfect_variable_selection_X_reduced_50)
        
        Method_TH_evaluate_variable_selection <- function_variable_reduction_evaluate_selection(Method_TH_variable_selection$selected_cols, Data$beta_1, eps = 1e-12)
        Method_CT_evaluate_variable_selection <- function_variable_reduction_evaluate_selection(Method_CT_variable_selection$selected_cols, Data$beta_1, eps = 1e-12)
        Method_SuSi_evaluate_variable_selection <- function_variable_reduction_evaluate_selection( Method_SuSi_variable_selection$selected_cols, Data$beta_1, eps = 1e-12)
        
        
        for (k_settings in subsample_percent_list) {
          # STEP 3 Loop over list of k ####
          # a. Subsample tagging based on prescreened matrix (TH, C+T, SUSIE, PVS, Uniform) ####
          # b. Train model based on the subdata ####
          # c. Record the result vector ####
          

          
          # Step 3.a Subsample tagging ####
          # k_settings = subsample_percent_list[1]
          subsample_percent = k_settings
          subsample_size <- round(subsample_percent/100 * nrow(X_train))
          
          Method_Uniform_subsampling = function_uniform_subsampling(sim_id, X_train, y_train,  subsample_size)
          Method_TH_IBOSS_subsampling = function_iboss_subsample_index(sim_id,Method_TH_variable_selection$X_reduced, subsample_size)
          Method_CT_IBOSS_subsampling = function_iboss_subsample_index(sim_id,Method_CT_variable_selection$X_reduced, subsample_size)
          Method_SuSi_IBOSS_subsampling = function_iboss_subsample_index(sim_id,Method_SuSi_variable_selection$X_reduced, subsample_size)
          Method_Perfect_IBOSS_subsampling = function_iboss_subsample_index(sim_id,Method_Perfect_variable_selection$X_reduced, subsample_size)
          
          # Step 3.b Subdata Model Training ####
          
          # i) p variate LASSO   ####
          LASSO_Subsample_Uniform  =   function_LASSO_output(sim_id, X_train[Method_Uniform_subsampling$subsample_index , ],
                                                             y_train[Method_Uniform_subsampling$subsample_index ], X_test, y_test)
          LASSO_Subsample_TH    =   function_LASSO_output(sim_id, X_train[ Method_TH_IBOSS_subsampling$subsample_index, ],
                                                          y_train[ Method_TH_IBOSS_subsampling$subsample_index ], X_test, y_test)
          LASSO_Subsample_CT    =   function_LASSO_output(sim_id,  X_train[ Method_CT_IBOSS_subsampling$subsample_index, ], 
                                                          y_train[ Method_CT_IBOSS_subsampling$subsample_index ], X_test, y_test)
          LASSO_Subsample_SuSi  =   function_LASSO_output(sim_id,X_train[ Method_SuSi_IBOSS_subsampling$subsample_index, ], 
                                                          y_train[ Method_SuSi_IBOSS_subsampling$subsample_index ] , X_test, y_test)
          LASSO_Subsample_Perfect  =   function_LASSO_output(sim_id,X_train[ Method_Perfect_IBOSS_subsampling$subsample_index, ], 
                                                             y_train[ Method_Perfect_IBOSS_subsampling$subsample_index ] , X_test, y_test)
          
          
          Final_variable_selection_performance_subsample_p_LASSO_uniform =  function_variable_reduction_evaluate_selection( LASSO_Subsample_Uniform$selected_variables, Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_p_LASSO_TH      =  function_variable_reduction_evaluate_selection( LASSO_Subsample_TH$selected_variables,      Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_p_LASSO_CT      =  function_variable_reduction_evaluate_selection( LASSO_Subsample_CT$selected_variables,      Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_p_LASSO_SuSi    =  function_variable_reduction_evaluate_selection( LASSO_Subsample_SuSi$selected_variables,    Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_p_LASSO_Perfect =  function_variable_reduction_evaluate_selection( LASSO_Subsample_Perfect$selected_variables, Data$beta_1, eps = 1e-12) 
          
          
          #   
          
          # ii) q  variate LASSO   ####
          LASSO_ScreenedVariables_Subsample_TH    =   function_LASSO_output(sim_id, X_train[ Method_TH_IBOSS_subsampling$subsample_index,Method_TH_variable_selection$selected_cols ],
                                                                            y_train[ Method_TH_IBOSS_subsampling$subsample_index ],
                                                                            X_test[,Method_TH_variable_selection$selected_cols], y_test)
          LASSO_ScreenedVariables_Subsample_CT    =   function_LASSO_output(sim_id,  X_train[ Method_CT_IBOSS_subsampling$subsample_index ,Method_CT_variable_selection$selected_cols ], 
                                                                            y_train[ Method_CT_IBOSS_subsampling$subsample_index ], 
                                                                            X_test[,Method_CT_variable_selection$selected_cols], y_test)
          LASSO_ScreenedVariables_Subsample_SuSi  =   function_LASSO_output(sim_id,X_train[ Method_SuSi_IBOSS_subsampling$subsample_index,Method_SuSi_variable_selection$selected_cols ],
                                                                            y_train[ Method_SuSi_IBOSS_subsampling$subsample_index ],
                                                                            X_test[,Method_SuSi_variable_selection$selected_cols], y_test)
          LASSO_ScreenedVariables_Subsample_Perfect  =   function_LASSO_output(sim_id,X_train[ Method_Perfect_IBOSS_subsampling$subsample_index,Method_Perfect_variable_selection$selected_cols ], 
                                                                               y_train[ Method_Perfect_IBOSS_subsampling$subsample_index ],
                                                                               X_test[,Method_Perfect_variable_selection$selected_cols], y_test)
          
          
          Final_variable_selection_performance_subsample_q_LASSO_TH      =  function_variable_reduction_evaluate_selection( Method_TH_variable_selection$selected_cols[LASSO_ScreenedVariables_Subsample_TH$selected_variables],      Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_q_LASSO_CT      =  function_variable_reduction_evaluate_selection( Method_CT_variable_selection$selected_cols[LASSO_ScreenedVariables_Subsample_CT$selected_variables],      Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_q_LASSO_SuSi    =  function_variable_reduction_evaluate_selection( Method_SuSi_variable_selection$selected_cols[LASSO_ScreenedVariables_Subsample_SuSi$selected_variables],    Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_q_LASSO_Perfect =  function_variable_reduction_evaluate_selection( Method_Perfect_variable_selection$selected_cols[LASSO_ScreenedVariables_Subsample_Perfect$selected_variables] , Data$beta_1, eps = 1e-12) 
          
          

         

          
          
          # iii) q  variate LASSO/OLS   ####
          LASSO_OLS_ScreenedVariables_Subsample_TH    =      function_OLS_LASSO_output(sim_id, X_train[ Method_TH_IBOSS_subsampling$subsample_index, Method_TH_variable_selection$selected_cols],
                                                                                       y_train[ Method_TH_IBOSS_subsampling$subsample_index ], X_test[, Method_TH_variable_selection$selected_cols], y_test)
          LASSO_OLS_ScreenedVariables_Subsample_CT    =      function_OLS_LASSO_output(sim_id, X_train[ Method_CT_IBOSS_subsampling$subsample_index, Method_CT_variable_selection$selected_cols], 
                                                                                       y_train[ Method_CT_IBOSS_subsampling$subsample_index ], X_test[,Method_CT_variable_selection$selected_cols], y_test)
          LASSO_OLS_ScreenedVariables_Subsample_SuSi  =      function_OLS_LASSO_output(sim_id, X_train[ Method_SuSi_IBOSS_subsampling$subsample_index, Method_SuSi_variable_selection$selected_cols], 
                                                                                       y_train[ Method_SuSi_IBOSS_subsampling$subsample_index ] , X_test[, Method_SuSi_variable_selection$selected_cols], y_test)
          LASSO_OLS_ScreenedVariables_Subsample_Perfect  =   function_OLS_LASSO_output(sim_id, X_train[ Method_Perfect_IBOSS_subsampling$subsample_index, Method_Perfect_variable_selection$selected_cols ], 
                                                                                       y_train[ Method_Perfect_IBOSS_subsampling$subsample_index ] , X_test[,Method_Perfect_variable_selection$selected_cols], y_test)
          
          
          Final_variable_selection_performance_subsample_q_OLSLASSO_TH      =  function_variable_reduction_evaluate_selection( Method_TH_variable_selection$selected_cols[LASSO_OLS_ScreenedVariables_Subsample_TH$selected_variables],           Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_q_OLSLASSO_CT      =  function_variable_reduction_evaluate_selection( Method_CT_variable_selection$selected_cols[LASSO_OLS_ScreenedVariables_Subsample_CT$selected_variables],           Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_q_OLSLASSO_SuSi    =  function_variable_reduction_evaluate_selection( Method_SuSi_variable_selection$selected_cols[LASSO_OLS_ScreenedVariables_Subsample_SuSi$selected_variables],       Data$beta_1, eps = 1e-12) 
          Final_variable_selection_performance_subsample_q_OLSLASSO_Perfect =  function_variable_reduction_evaluate_selection( Method_Perfect_variable_selection$selected_cols[LASSO_OLS_ScreenedVariables_Subsample_Perfect$selected_variables], Data$beta_1, eps = 1e-12) 
          
          
          # Step 3.c Record the result vector ####
          {
            result_case <- data.frame(
              n_train = n_train, 
              n_test = n_test,
              LD_block_size = LD_block_size ,
              n_LD_blocks = n_LD_blocks , 
              
              DataSet_Name = DataSet_Name,
              covariate_covariance_structure = covariate_covariance_structure,
              noise_variance = "TBD",
              intercept = "Yes",
              
              
              LD_rho = LD_rho, 
              k = subsample_size,
              subsample_percent = subsample_percent,
              # num_models = num_models,
              sparsity = sparsity,
              sparsity_corrected = sparsity_corrected, 
              SNR = SNR,
              beta_type = beta_type, 
              num_simulations = num_simulations,
              
              subsample_percent_variable_selection = subsample_percent_variable_selection,
              
              
              
              
              
              
              Eval_TH_SignalsSelected = Method_TH_evaluate_variable_selection$num_selected,
              Eval_TH_Power = Method_TH_evaluate_variable_selection$TPR,
              Eval_TH_Precision = Method_TH_evaluate_variable_selection$PPV,
              Eval_CT_SignalsSelected = Method_CT_evaluate_variable_selection$num_selected,
              Eval_CT_Power = Method_CT_evaluate_variable_selection$TPR,
              Eval_CT_Precision = Method_CT_evaluate_variable_selection$PPV,
              Eval_SuSi_SignalsSelected = Method_SuSi_evaluate_variable_selection$num_selected,
              Eval_SuSi_Power = Method_SuSi_evaluate_variable_selection$TPR,
              Eval_SuSi_Precision = Method_SuSi_evaluate_variable_selection$PPV,
              Eval_Perfect_SignalsSelected = Method_Perfect_evaluate_variable_selection$num_selected,
              Eval_Perfect_Power = Method_Perfect_evaluate_variable_selection$TPR,
              Eval_Perfect_Precision = Method_Perfect_evaluate_variable_selection$PPV  ,
              
              
              # subdata_LASSO- Elapsed time, mse,r_sq
              Elapsed_time_total_LASSO_Fullsample = LASSO_Fullsample$elapsed_time,
              Elapsed_time_total_LASSO_Subsample_Uniform =  Method_Uniform_subsampling$elapsed_time  + LASSO_Subsample_Uniform$elapsed_time,
              Elapsed_time_total_LASSO_Subsample_TH =   Method_TH_variable_selection$elapsed_time +    Method_TH_IBOSS_subsampling$elapsed_time +   LASSO_Subsample_TH$elapsed_time,
              Elapsed_time_total_LASSO_Subsample_CT =   Method_CT_variable_selection$elapsed_time +    Method_CT_IBOSS_subsampling$elapsed_time +   LASSO_Subsample_CT$elapsed_time,
              Elapsed_time_total_LASSO_Subsample_SuSi = Method_SuSi_variable_selection$elapsed_time +  Method_SuSi_IBOSS_subsampling$elapsed_time + LASSO_Subsample_SuSi$elapsed_time,
              Elapsed_time_total_LASSO_Subsample_Perfect = Method_Perfect_variable_selection$elapsed_time +  Method_Perfect_IBOSS_subsampling$elapsed_time + LASSO_Subsample_Perfect$elapsed_time,
              
              mse_LASSO_Fullsample = LASSO_Fullsample$mse,
              mse_LASSO_Subsample_Uniform = LASSO_Subsample_Uniform$mse,
              mse_LASSO_Subsample_TH = LASSO_Subsample_TH$mse,
              mse_LASSO_Subsample_CT = LASSO_Subsample_CT$mse,
              mse_LASSO_Subsample_SuSi =  LASSO_Subsample_SuSi$mse,
              mse_LASSO_Subsample_Perfect =  LASSO_Subsample_Perfect$mse,
              
              r_sq_test_LASSO_Fullsample = LASSO_Fullsample$r_sq_test,
              r_sq_test_LASSO_Subsample_Uniform = LASSO_Subsample_Uniform$r_sq_test,
              r_sq_test_LASSO_Subsample_TH = LASSO_Subsample_TH$r_sq_test,
              r_sq_test_LASSO_Subsample_CT = LASSO_Subsample_CT$r_sq_test,
              r_sq_test_LASSO_Subsample_SuSi = LASSO_Subsample_SuSi$r_sq_test,
              r_sq_test_LASSO_Subsample_Perfect = LASSO_Subsample_Perfect$r_sq_test,
              
              # subdata_LASSOScreenedVar - Elapsed time, mse,r_sq
              Elapsed_time_total_LASSO_ScreenedVariables_Subsample_TH =   Method_TH_variable_selection$elapsed_time +    Method_TH_IBOSS_subsampling$elapsed_time +   LASSO_ScreenedVariables_Subsample_TH$elapsed_time,
              Elapsed_time_total_LASSO_ScreenedVariables_Subsample_CT =   Method_CT_variable_selection$elapsed_time +    Method_CT_IBOSS_subsampling$elapsed_time +   LASSO_ScreenedVariables_Subsample_CT$elapsed_time,
              Elapsed_time_total_LASSO_ScreenedVariables_Subsample_SuSi = Method_SuSi_variable_selection$elapsed_time +  Method_SuSi_IBOSS_subsampling$elapsed_time + LASSO_ScreenedVariables_Subsample_SuSi$elapsed_time,
              Elapsed_time_total_LASSO_ScreenedVariables_Subsample_Perfect = Method_Perfect_variable_selection$elapsed_time +  Method_Perfect_IBOSS_subsampling$elapsed_time + LASSO_ScreenedVariables_Subsample_Perfect$elapsed_time,
              
              mse_LASSO_ScreenedVariables_Subsample_TH = LASSO_ScreenedVariables_Subsample_TH$mse,
              mse_LASSO_ScreenedVariables_Subsample_CT = LASSO_ScreenedVariables_Subsample_CT$mse,
              mse_LASSO_ScreenedVariables_Subsample_SuSi =  LASSO_ScreenedVariables_Subsample_SuSi$mse,
              mse_LASSO_ScreenedVariables_Subsample_Perfect =  LASSO_ScreenedVariables_Subsample_Perfect$mse,
              
              r_sq_test_LASSO_ScreenedVariables_Subsample_TH = LASSO_ScreenedVariables_Subsample_TH$r_sq_test,
              r_sq_test_LASSO_ScreenedVariables_Subsample_CT = LASSO_ScreenedVariables_Subsample_CT$r_sq_test,
              r_sq_test_LASSO_ScreenedVariables_Subsample_SuSi = LASSO_ScreenedVariables_Subsample_SuSi$r_sq_test,
              r_sq_test_LASSO_ScreenedVariables_Subsample_Perfect = LASSO_ScreenedVariables_Subsample_Perfect$r_sq_test,
              
              # subdata_LASSO/OLS_ScreenedVar - Elapsed time, mse, r_sq 
              
              
              
              Elapsed_time_total_LASSO_OLS_ScreenedVariables_Subsample_TH =      Method_TH_variable_selection$elapsed_time +       Method_TH_IBOSS_subsampling$elapsed_time +      LASSO_OLS_ScreenedVariables_Subsample_TH$elapsed_time,
              Elapsed_time_total_LASSO_OLS_ScreenedVariables_Subsample_CT =      Method_CT_variable_selection$elapsed_time +       Method_CT_IBOSS_subsampling$elapsed_time +      LASSO_OLS_ScreenedVariables_Subsample_CT$elapsed_time,
              Elapsed_time_total_LASSO_OLS_ScreenedVariables_Subsample_SuSi =    Method_SuSi_variable_selection$elapsed_time +     Method_SuSi_IBOSS_subsampling$elapsed_time +    LASSO_OLS_ScreenedVariables_Subsample_SuSi$elapsed_time,
              Elapsed_time_total_LASSO_OLS_ScreenedVariables_Subsample_Perfect = Method_Perfect_variable_selection$elapsed_time +  Method_Perfect_IBOSS_subsampling$elapsed_time + LASSO_OLS_ScreenedVariables_Subsample_Perfect$elapsed_time,
              
              
              mse_LASSO_OLS_ScreenedVariables_Subsample_TH = LASSO_OLS_ScreenedVariables_Subsample_TH$mse,
              mse_LASSO_OLS_ScreenedVariables_Subsample_CT = LASSO_OLS_ScreenedVariables_Subsample_CT$mse,
              mse_LASSO_OLS_ScreenedVariables_Subsample_SuSi =  LASSO_OLS_ScreenedVariables_Subsample_SuSi$mse,
              mse_LASSO_OLS_ScreenedVariables_Subsample_Perfect =  LASSO_OLS_ScreenedVariables_Subsample_Perfect$mse,
              
              
              
              r_sq_test_LASSO_OLS_ScreenedVariables_Subsample_TH = LASSO_OLS_ScreenedVariables_Subsample_TH$r_sq_test,
              r_sq_test_LASSO_OLS_ScreenedVariables_Subsample_CT = LASSO_OLS_ScreenedVariables_Subsample_CT$r_sq_test,
              r_sq_test_LASSO_OLS_ScreenedVariables_Subsample_SuSi = LASSO_OLS_ScreenedVariables_Subsample_SuSi$r_sq_test,
              r_sq_test_LASSO_OLS_ScreenedVariables_Subsample_Perfect = LASSO_OLS_ScreenedVariables_Subsample_Perfect$r_sq_test,
              
              
              
              # Variable selection - Elapsed time
              Elapsed_time_Method_TH_variable_selection = Method_TH_variable_selection$elapsed_time,
              Elapsed_time_Method_CT_variable_selection = Method_CT_variable_selection$elapsed_time ,
              Elapsed_time_Method_SuSi_variable_selection = Method_SuSi_variable_selection$elapsed_time,
              # SubsampleTagging selection - Elapsed time
              Elapsed_time_SS_tag_Method_Uniform_IBOSS = Method_Uniform_subsampling$elapsed_time,
              Elapsed_time_SS_IBOSS_tag_Method_TH_IBOSS = Method_TH_IBOSS_subsampling$elapsed_time,
              Elapsed_time_SS_IBOSS_tag_Method_CT_IBOSS = Method_CT_IBOSS_subsampling$elapsed_time,
              Elapsed_time_SS_IBOSS_tag_Method_SuSi_IBOSS = Method_SuSi_IBOSS_subsampling$elapsed_time,
              Elapsed_time_SS_IBOSS_tag_Method_Perfect_IBOSS = Method_Perfect_IBOSS_subsampling$elapsed_time,
              # LASSO on suddataAll variables- Elapsed time
              Elapsed_time_LASSO_Fullsample = LASSO_Fullsample$elapsed_time,
              Elapsed_time_LASSO_Subsample_Uniform = LASSO_Subsample_Uniform$elapsed_time,
              Elapsed_time_LASSO_Subsample_TH = LASSO_Subsample_TH$elapsed_time,
              Elapsed_time_LASSO_Subsample_CT = LASSO_Subsample_CT$elapsed_time,
              Elapsed_time_LASSO_Subsample_SuSi = LASSO_Subsample_SuSi$elapsed_time,
              Elapsed_time_LASSO_Subsample_Perfect = LASSO_Subsample_Perfect$elapsed_time,
              # LASSO on prescreened variables- Elapsed time?
              # OLS/LASSO on prescreened variables- Elapsed time?
              
              
              Elapsed_time_LASSO_ScreenedVariables_Subsample_TH = LASSO_ScreenedVariables_Subsample_TH$elapsed_time,
              Elapsed_time_LASSO_ScreenedVariables_Subsample_CT = LASSO_ScreenedVariables_Subsample_CT$elapsed_time,
              Elapsed_time_LASSO_ScreenedVariables_Subsample_SuSi = LASSO_ScreenedVariables_Subsample_SuSi$elapsed_time,
              Elapsed_time_LASSO_ScreenedVariables_Subsample_Perfect = LASSO_ScreenedVariables_Subsample_Perfect$elapsed_time,
            
            # Final Variable selection performance
                ## Full sample lasso
            NumSelected_Final_variable_selection_fullsample_LASSO = Final_variable_selection_fullsample_LASSO$num_selected,
            Power_Final_variable_selection_fullsample_LASSO =   Final_variable_selection_fullsample_LASSO$TPR,
            Precision_Final_variable_selection_fullsample_LASSO_Precision =   Final_variable_selection_fullsample_LASSO$PPV,
            
               ## subsample p variate lasso  
            NumSelected_Final_variable_selection_performance_subsample_p_LASSO_uniform  =  Final_variable_selection_performance_subsample_p_LASSO_uniform$num_selected ,
            Power_Final_variable_selection_performance_subsample_p_LASSO_uniform  =  Final_variable_selection_performance_subsample_p_LASSO_uniform$TPR ,
            Precision_Final_variable_selection_performance_subsample_p_LASSO_uniform = Final_variable_selection_performance_subsample_p_LASSO_uniform$PPV , 
            
            NumSelected_Final_variable_selection_performance_subsample_p_LASSO_TH = Final_variable_selection_performance_subsample_p_LASSO_TH$num_selected,
            Power_Final_variable_selection_performance_subsample_p_LASSO_TH  = Final_variable_selection_performance_subsample_p_LASSO_TH$TPR ,
            Precision_Final_variable_selection_performance_subsample_p_LASSO_TH = Final_variable_selection_performance_subsample_p_LASSO_TH$PPV,             
            
            # length(LASSO_Subsample_CT$selected_variables)
            NumSelected_Final_variable_selection_performance_subsample_p_LASSO_CT = Final_variable_selection_performance_subsample_p_LASSO_CT$num_selected,
            Power_Final_variable_selection_performance_subsample_p_LASSO_CT  = Final_variable_selection_performance_subsample_p_LASSO_CT$TPR ,
            Precision_Final_variable_selection_performance_subsample_p_LASSO_CT = Final_variable_selection_performance_subsample_p_LASSO_CT$PPV,             
            
            NumSelected_Final_variable_selection_performance_subsample_p_LASSO_SuSi = Final_variable_selection_performance_subsample_p_LASSO_SuSi$num_selected,
            Power_Final_variable_selection_performance_subsample_p_LASSO_SuSi  = Final_variable_selection_performance_subsample_p_LASSO_SuSi$TPR ,
            Precision_Final_variable_selection_performance_subsample_p_LASSO_SuSi = Final_variable_selection_performance_subsample_p_LASSO_SuSi$PPV,             
            
            NumSelected_Final_variable_selection_performance_subsample_p_LASSO_Perfect= Final_variable_selection_performance_subsample_p_LASSO_Perfect$num_selected,
            Power_Final_variable_selection_performance_subsample_p_LASSO_Perfect  = Final_variable_selection_performance_subsample_p_LASSO_Perfect$TPR ,
            Precision_Final_variable_selection_performance_subsample_p_LASSO_Perfect = Final_variable_selection_performance_subsample_p_LASSO_Perfect$PPV,             
            
               ## subsample q variate lasso  
            NumSelected_Final_variable_selection_performance_subsample_q_LASSO_TH = Final_variable_selection_performance_subsample_q_LASSO_TH$num_selected,
            Power_Final_variable_selection_performance_subsample_q_LASSO_TH  = Final_variable_selection_performance_subsample_q_LASSO_TH$TPR ,
            Precision_Final_variable_selection_performance_subsample_q_LASSO_TH = Final_variable_selection_performance_subsample_q_LASSO_TH$PPV,             
            
            NumSelected_Final_variable_selection_performance_subsample_q_LASSO_CT = Final_variable_selection_performance_subsample_q_LASSO_CT$num_selected,
            Power_Final_variable_selection_performance_subsample_q_LASSO_CT   = Final_variable_selection_performance_subsample_q_LASSO_CT$TPR ,
            Precision_Final_variable_selection_performance_subsample_q_LASSO_CT  = Final_variable_selection_performance_subsample_q_LASSO_CT$PPV,             
            
            NumSelected_Final_variable_selection_performance_subsample_q_LASSO_SuSi   = Final_variable_selection_performance_subsample_q_LASSO_SuSi$num_selected,
            Power_Final_variable_selection_performance_subsample_q_LASSO_SuSi  = Final_variable_selection_performance_subsample_q_LASSO_SuSi$TPR ,
            Precision_Final_variable_selection_performance_subsample_q_LASSO_SuSi = Final_variable_selection_performance_subsample_q_LASSO_SuSi$PPV,             
            
            NumSelected_Final_variable_selection_performance_subsample_q_LASSO_Perfect = Final_variable_selection_performance_subsample_q_LASSO_Perfect$num_selected,
            Power_Final_variable_selection_performance_subsample_q_LASSO_Perfect  = Final_variable_selection_performance_subsample_q_LASSO_Perfect$TPR ,
            Precision_Final_variable_selection_performance_subsample_q_LASSO_Perfect = Final_variable_selection_performance_subsample_q_LASSO_Perfect$PPV,             
            
              ## subsample q variate ols/lasso  
            NumSelected_Final_variable_selection_performance_subsample_q_OLSLASSO_TH = Final_variable_selection_performance_subsample_q_OLSLASSO_TH$num_selected,
            Power_Final_variable_selection_performance_subsample_q_OLSLASSO_TH   = Final_variable_selection_performance_subsample_q_OLSLASSO_TH$TPR,
            Precision_Final_variable_selection_performance_subsample_q_OLSLASSO_TH  = Final_variable_selection_performance_subsample_q_OLSLASSO_TH$PPV,             
            
            NumSelected_Final_variable_selection_performance_subsample_q_OLSLASSO_CT = Final_variable_selection_performance_subsample_q_OLSLASSO_CT$num_selected,
            Power_Final_variable_selection_performance_subsample_q_OLSLASSO_CT = Final_variable_selection_performance_subsample_q_OLSLASSO_CT$TPR ,
            Precision_Final_variable_selection_performance_subsample_q_OLSLASSO_CT =Final_variable_selection_performance_subsample_q_OLSLASSO_CT$PPV,                 
            
            NumSelected_Final_variable_selection_performance_subsample_q_OLSLASSO_SuSi = Final_variable_selection_performance_subsample_q_OLSLASSO_SuSi$num_selected,
            Power_Final_variable_selection_performance_subsample_q_OLSLASSO_SuSi  = Final_variable_selection_performance_subsample_q_OLSLASSO_SuSi$TPR ,
            Precision_Final_variable_selection_performance_subsample_q_OLSLASSO_SuSi = Final_variable_selection_performance_subsample_q_OLSLASSO_SuSi$PPV,             
            
            NumSelected_Final_variable_selection_performance_subsample_q_OLSLASSO_Perfect  = Final_variable_selection_performance_subsample_q_OLSLASSO_Perfect$num_selected,
            Power_Final_variable_selection_performance_subsample_q_OLSLASSO_Perfect  = Final_variable_selection_performance_subsample_q_OLSLASSO_Perfect$TPR ,
            Precision_Final_variable_selection_performance_subsample_q_OLSLASSO_Perfect = Final_variable_selection_performance_subsample_q_OLSLASSO_Perfect$PPV,
            
            
            
            
            
            lambda_LASSO_Fullsample = LASSO_Fullsample$optimal_lambda, 
            
            lambda_LASSO_Subsample_Uniform = LASSO_Subsample_Uniform$optimal_lambda,
            lambda_LASSO_Subsample_TH = LASSO_Subsample_TH$optimal_lambda,
            lambda_LASSO_Subsample_CT = LASSO_Subsample_CT$optimal_lambda,
            lambda_LASSO_Subsample_SuSi = LASSO_Subsample_SuSi$optimal_lambda,
            lambda_LASSO_Subsample_Perfect = LASSO_Subsample_Perfect$optimal_lambda,
            
            lambda_LASSO_ScreenedVariables_Subsample_TH = LASSO_ScreenedVariables_Subsample_TH$optimal_lambda,    
            lambda_LASSO_ScreenedVariables_Subsample_CT = LASSO_ScreenedVariables_Subsample_CT$optimal_lambda,
            lambda_LASSO_ScreenedVariables_Subsample_SuSi = LASSO_ScreenedVariables_Subsample_SuSi$optimal_lambda,
            lambda_LASSO_ScreenedVariables_Subsample_Perfect = LASSO_ScreenedVariables_Subsample_Perfect$optimal_lambda
            )
          }
          result = rbind(result, result_case)
          
        }
        
       }
  write.csv(result, paste0("Data_Temp/Temp_iter_", i,  "_sim_", sim_id, ".csv"), row.names = FALSE)
  return(result )
  
  
}





for (i in seq_len(nrow(params))){
  
  n_train = params[i,1]
  LD_rho = params[i,2]
  SNR =  params[i,3]
  beta_type = as.character(params[i,4])
  sparsity = params[i,5]

  
  results_list <- foreach(sim_id = seq_len(num_simulations), 
                          .combine = rbind,             # Combine results row-wise into a single data frame
                          .packages = c("stats","glmnet", "MASS", "pracma", "susieR", "dplyr")
  ) %dorng% {
    # To make sure worker has access to helper functions
    source("function_generate_data_version_1.R")
    source("function_variable_reduction_methods_version2.R")
    source("Function_Subsample_Tagging_Algorithms.R")
    source("function_LASSO_output.R")
    
    
    function_run_ensembles_folded(
      sim_id = sim_id,
      n_train = n_train,
      LD_rho = LD_rho,
      SNR = SNR,
      beta_type = beta_type, 
      subsample_percent_variable_selection_list, 
      subsample_percent_list    )
  }
  
  
  
  cat(paste0("The iteration for parameter setup ", i, " is complete now.\n"))
  all_results = rbind(all_results, results_list)
  write.csv(all_results, paste0( "Data/All_Results_Dataset_", DataSet_Name, "p_", LD_block_size,  "times ",  n_LD_blocks,   "_Series_", series,  ".csv"), row.names = FALSE)
  
  
  
}



stopCluster(cl)


