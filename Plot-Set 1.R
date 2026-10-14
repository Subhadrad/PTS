# Plots: elapsed time and rho_sq Vs k
source("plot_functions.R")

# Simulated data setup
{

  data <- read.csv("Data/All_Results_Dataset_Simulated LD Structurep_1000times 5_Series_Series 50 c 1000 by 5 by 5 sparsity cases L10 New.csv")

  data <- data %>%
    mutate(
      beta_type = recode(
        beta_type,
        "LD-5 sparsity Type-I" = "LD-5 sparsity I"
        #,
        #"LD-5 sparsity"        = "LD-5 sparsity II"
      )
    )
  

  set_n_train  = 1e+05 # set_n_train  =1e+05 # set_n_train = unique(data$n_train)
  # k is the variable on x axis
  set_subsample_percent_variable_selection = 50 # 100  75  50  25 # unique(data$subsample_percent_variable_selection)
  # set_p = 1000  # unique(data$LD_block_size) *   unique(data$n_LD_blocks)
  set_LD_block_size = 1000  # 200 500 1000 # unique(data$LD_block_size)
  set_n_LD_blocks = 5 # unique(data$n_LD_blocks)
  set_DataSet_Name = "Simulated LD Structure" # unique(data$DataSet_Name)
  # set_sparsity  = 0.24 # 1.2 1 # # 0.48 1 #  0.24 1.00 # # unique(data$sparsity_corrected)
  set_SNR = 0.5 #0.5  # unique(data$SNR) 
  
  
  # set_LD_rho = 0.95 # unique(data$LD_rho)
  # set_beta_type  = "LD-5"  # "LD-5" "LD-5 sparsity I" "LD-5 sparsity II"   #  unique(data$beta_type)
  # set_sparsity  = 0.24 # 1.2 1 # # 0.48 1 #  0.24 1.00 # # unique(data$sparsity_corrected)
  # plot_title = TeX("$Fixed\\; s_{a}$")
  # set_beta_type  =     "LD-5 sparsity I" # "LD-5 sparsity II"   #  unique(data$beta_type)
  # set_sparsity  =  1 # 1.2 1 # # 0.48 1 #  0.24 1.00 # # unique(data$sparsity_corrected)
  # plot_title = TeX("$Fixed\\; s$")
  
  
  # set_LD_rho = 0.95 # unique(data$LD_rho)
  # set_LD_rho = 0.1 # unique(data$LD_rho)
  set_LD_rho = 0.5 # unique(data$LD_rho)
  
  set_beta_type  = "LD-5"  # "LD-5" "LD-5 sparsity I" "LD-5 sparsity II"   #  unique(data$beta_type)
  set_sparsity  = 0.24 # 1.2 1 # # 0.48 1 #  0.24 1.00 # # unique(data$sparsity_corrected)
  plot_title = TeX("$Fixed\\; s_{a}$")
  
  set_beta_type  =     "LD-5 sparsity I" # "LD-5 sparsity II"   #  unique(data$beta_type)
  set_sparsity  =  1 # 1.2 1 # # 0.48 1 #  0.24 1.00 # # unique(data$sparsity_corrected)
  plot_title = TeX("$Fixed\\; s$")


  
}
table(data$sparsity_corrected) 



# q variables- R sqaured and time ####
# Plotting R sqaured - LASSO with q variables 
{
  variable = c("r_sq_test_LASSO_Fullsample", "r_sq_test_LASSO_Subsample_Uniform")  
  variable_like = "r_sq_test_LASSO_ScreenedVariables" 
  
  summary_like = "r_sq_test_"
  values_to_var = "R_Sq_Reduced_Var"
  
  algorithm_map = c(
    "r_sq_test_LASSO_Fullsample_mean" = "Full Sample",
    "r_sq_test_LASSO_Subsample_Uniform_mean" = "Uniform Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_TH_mean" = "TH IBOSS Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_CT_mean" = "CT IBOSS Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_SuSi_mean" = "SuSiE IBOSS Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_Perfect_mean" = "Oracle IBOSS Subsample"
  )
  # y_axis_label = expression(r^2)
  # y_axis_label =  TeX("$\\bar{r}^{ 2}_{test}$")
  y_axis_label =  TeX("${r}^{ 2}_{test}$")
  x_axis_var = "k" 
  file_name_start = "q_R_sq_RedVar"
  

  # LD-5 0.95
  function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                         set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
                         set_beta_type, set_subsample_percent_variable_selection, variable  , variable_like  ,
                         summary_like , values_to_var ,
                         algorithm_map , y_axis_label   ,
                         x_axis_var, file_name_start ,
                         legend_position = "none",
                         plot_title = "",
  y_axis_break_range = c(0.20, 0.275),
                         y_axis_break_scales = 2.5, ymin = .19)
  ## LD-5 Sparsity Type I 0.95
  # function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
  #                        set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
  #                        set_beta_type, set_subsample_percent_variable_selection, variable  , variable_like  ,
  #                        summary_like , values_to_var ,
  #                        algorithm_map ,
  #                        y_axis_label = "",
  #                        x_axis_var, file_name_start ,
  #                        legend_position = "none",
  #                        plot_title = ""
  #                        # ,
  #                        # y_axis_break_range = c(0.20, 0.29),
  #                        # y_axis_break_scales = 2.5, ymin = .19
  #                        )


  # # LD-5 0.50
  function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                         set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
                         set_beta_type, set_subsample_percent_variable_selection, variable  , variable_like  ,
                         summary_like , values_to_var ,
                         algorithm_map , y_axis_label   ,
                         x_axis_var, file_name_start ,
                         legend_position = "none",
                         plot_title = ""
                         # ,
  # y_axis_break_range = c(0.20, 0.29),
                         # y_axis_break_scales = 2.5, ymin = .19
                         )

  # # LD-5 Sparsity Type I 0.50
  # function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
  #                        set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
  #                        set_beta_type, set_subsample_percent_variable_selection, variable  , variable_like  ,
  #                        summary_like , values_to_var ,
  #                        algorithm_map ,
  #                        y_axis_label = "",
  #                        x_axis_var, file_name_start ,
  #                        legend_position = "none",
  #                        plot_title = "",
  #                        # ,
  #                        # y_axis_break_range = c(0.20, 0.29),
  #                        #  y_axis_break_scales = 2.5,
  #                        ymin = 0.0
  #                        )


  
  
  
  
}
# Plotting Time - LASSO with q variables 
{
  variable = c("Elapsed_time_total_LASSO_Fullsample",
               "Elapsed_time_total_LASSO_Subsample_Uniform")  
  variable_like = "Elapsed_time_total_LASSO_ScreenedVariables_Subsample" 

  summary_like = "Elapsed_time_total"
  values_to_var = "Elapsed_time_total" 
  

  algorithm_map = c(
    "Elapsed_time_total_LASSO_Fullsample_mean" = "Full Sample",
    "Elapsed_time_total_LASSO_Subsample_Uniform_mean" = "Uniform Subsample",
    "Elapsed_time_total_LASSO_ScreenedVariables_Subsample_TH_mean" = "TH IBOSS Subsample",
    "Elapsed_time_total_LASSO_ScreenedVariables_Subsample_CT_mean" = "CT IBOSS Subsample",
    "Elapsed_time_total_LASSO_ScreenedVariables_Subsample_SuSi_mean" = "SuSiE IBOSS Subsample",
    "Elapsed_time_total_LASSO_ScreenedVariables_Subsample_Perfect_mean" = "Oracle IBOSS Subsample"
  )
  y_axis_label =  expression(italic(T) ~ "(sec)")
  # y_axis_label = expression(bar(italic(T)) ~ "(sec)")
  
  x_axis_var = "k"
  file_name_start = "q_Time" 
  
  legend_position = "bottom"
  plot_title = "Model training time"
  
  
  
  # LD-5 rho 0.95
  function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                         set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
                         set_beta_type, set_subsample_percent_variable_selection,
                         variable  , variable_like  ,
                         summary_like , values_to_var ,
                         algorithm_map , y_axis_label = y_axis_label  ,  plot_title = "",
                         x_axis_var, file_name_start , legend_position = "none",
                         y_axis_break_range = c(55, 1200),
                         y_axis_break_scales = 0.2 )


  
  
  # LD-5 Sparsity I rho 0.95
  # function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
  #                        set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
  #                        set_beta_type, set_subsample_percent_variable_selection,
  #                        variable  , variable_like  ,
  #                        summary_like , values_to_var ,
  #                        algorithm_map , y_axis_label = " "  ,  plot_title = "",
  #                        x_axis_var, file_name_start , legend_position = "none",
  #                        y_axis_break_range = c(150, 1200),
  #                        y_axis_break_scales = 0.2 )
  # 
  
  
  
  # LD-5 rho 0.5
  # function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
  #                        set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
  #                        set_beta_type, set_subsample_percent_variable_selection,
  #                        variable  , variable_like  ,
  #                        summary_like , values_to_var ,
  #                        algorithm_map , y_axis_label = y_axis_label  ,  plot_title = "",
  #                        x_axis_var, file_name_start , legend_position = "none",
  #                        y_axis_break_range = c(60, 550),
  #                        y_axis_break_scales = 0.2 )
  # # LD-5 Sparsity I rho 0.5
  # function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
  #                        set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
  #                        set_beta_type, set_subsample_percent_variable_selection,
  #                        variable  , variable_like  ,
  #                        summary_like , values_to_var ,
  #                        algorithm_map , y_axis_label = " "  ,  plot_title = "",
  #                        x_axis_var, file_name_start , legend_position = "none",
  #                        y_axis_break_range = c(80, 500),
  #                        y_axis_break_scales = 0.2 )
  # 
  
  
}

# q variables- Power and precision ####
# Plotting power and precision - of LASSO on q variable 
{
  variable = c("Power_Final_variable_selection_fullsample_LASSO", "Power_Final_variable_selection_performance_subsample_p_LASSO_uniform")
  variable_like = "Power_Final_variable_selection_performance_subsample_q_LASSO"

  summary_like = "Power"
  values_to_var = "Power"


  algorithm_map = c(
    "Power_Final_variable_selection_fullsample_LASSO_mean" = "Full Sample",
    "Power_Final_variable_selection_performance_subsample_p_LASSO_uniform_mean" = "Uniform Subsample",
    "Power_Final_variable_selection_performance_subsample_q_LASSO_TH_mean" = "TH IBOSS Subsample",
    "Power_Final_variable_selection_performance_subsample_q_LASSO_CT_mean" = "CT IBOSS Subsample",
    "Power_Final_variable_selection_performance_subsample_q_LASSO_SuSi_mean" = "SuSiE IBOSS Subsample",
    "Power_Final_variable_selection_performance_subsample_q_LASSO_Perfect_mean" = "Oracle IBOSS Subsample"
  )
  y_axis_label = "Power"
  x_axis_var = "k"
  file_name_start = "Power_q_var"

  
  function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                         set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
                         set_beta_type, set_subsample_percent_variable_selection,
                         variable  , variable_like  ,
                         summary_like , values_to_var ,
                         algorithm_map , y_axis_label  ,
                         x_axis_var, file_name_start, legend_position = "bottom" )

  # This is for LD-5
  function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                         set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
                         set_beta_type, set_subsample_percent_variable_selection,
                         variable  , variable_like  ,
                         summary_like , values_to_var ,
                         algorithm_map , y_axis_label  ,
                         x_axis_var, file_name_start, legend_position = "none",
                         plot_title = ""   )
  # LD-5 Sparsity I
  function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                         set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
                         set_beta_type, set_subsample_percent_variable_selection,
                         variable  , variable_like  ,
                         summary_like , values_to_var ,
                         algorithm_map , y_axis_label =""  ,
                         x_axis_var, file_name_start, legend_position = "none", 
                         plot_title = "")
  
  
}
{
  variable = c("Precision_Final_variable_selection_fullsample_LASSO_Precision",
               "Precision_Final_variable_selection_performance_subsample_p_LASSO_uniform")  
  variable_like = "Precision_Final_variable_selection_performance_subsample_q_LASSO" 
  
  summary_like = "Precision"
  values_to_var = "Precision" 
  
  
  algorithm_map = c(
    "Precision_Final_variable_selection_fullsample_LASSO_Precision_mean" = "Full Sample",
    "Precision_Final_variable_selection_performance_subsample_p_LASSO_uniform_mean" = "Uniform Subsample",
    "Precision_Final_variable_selection_performance_subsample_q_LASSO_TH_mean" = "TH IBOSS Subsample",
    "Precision_Final_variable_selection_performance_subsample_q_LASSO_CT_mean" = "CT IBOSS Subsample",
    "Precision_Final_variable_selection_performance_subsample_q_LASSO_SuSi_mean" = "SuSiE IBOSS Subsample",
    "Precision_Final_variable_selection_performance_subsample_q_LASSO_Perfect_mean" = "Oracle IBOSS Subsample"
  )
  y_axis_label = "Precision"
  x_axis_var = "k"
  file_name_start = "Precision_q_var" 
  
  
  # This is for LD-5
  function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                         set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR,
                         set_beta_type, set_subsample_percent_variable_selection,
                         variable  , variable_like  ,
                         summary_like , values_to_var ,
                         algorithm_map , y_axis_label  ,
                         x_axis_var, file_name_start ,  legend_position = "none",
                         plot_title = "")

  # This is for LD-5 Sparsity I 
  function_plots_xaxis_k(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                         set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR, 
                         set_beta_type, set_subsample_percent_variable_selection,
                         variable  , variable_like  ,
                         summary_like , values_to_var , 
                         algorithm_map , y_axis_label = " "  ,
                         x_axis_var,
                         plot_title = "",
                         file_name_start ,  legend_position = "none") 
  
  
}


