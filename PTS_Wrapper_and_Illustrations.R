source("Function_PTS.R")


example_data <- function_generate_example_data()
X_train <- example_data$X_train
y_train <- example_data$y_train
ld_regions <- example_data$ld_regions
beta_true <- example_data$beta_true

initial_prescreeing_random_subsample_percent = 50
screening_methods = c("TH", "CT", "SuSiE") 
threshold = 0.05
clumping_threshold = 0.1 
L=10
subsample_size = 200
subsampling_methods = "IBOSS"
pip_cutoff = 0.5
include_CS = FALSE #TRUE
top_q = NULL # top_q = 10, corresponds SIS implementation which returns only top 10 marginally correlated covariates 
seed = 1


subdata_selection = function_select_PTS_subdata(X_train, y_train, 
                                              screening_methods, initial_prescreeing_random_subsample_percent,  
                                              threshold, top_q,
                                              clumping_threshold, 
                                              L, ld_regions, SuSie_pval_cutoff = threshold, pip_cutoff, include_CS, 
                                            subsampling_methods, subsample_size,
                                                        seed)



# Subdata Selected by PTS ####
subdata_selection$TH_IBOSS$screening
subdata_selection$TH_IBOSS$subsampling
subdata_selection$TH_IBOSS$stage_I_prescreened_var
subdata_selection$TH_IBOSS$stage_II_subsample_indices
subdata_selection$TH_IBOSS$X_subdata
subdata_selection$TH_IBOSS$y_subdata

subdata_selection$CT_IBOSS$screening
subdata_selection$CT_IBOSS$subsampling
subdata_selection$CT_IBOSS$stage_I_prescreened_var
subdata_selection$CT_IBOSS$stage_II_subsample_indices
subdata_selection$CT_IBOSS$X_subdata
subdata_selection$CT_IBOSS$y_subdata


subdata_selection$SuSiE_IBOSS$screening
subdata_selection$SuSiE_IBOSS$subsampling
subdata_selection$SuSiE_IBOSS$stage_I_prescreened_var
subdata_selection$SuSiE_IBOSS$stage_II_subsample_indices
subdata_selection$SuSiE_IBOSS$X_subdata
subdata_selection$SuSiE_IBOSS$y_subdata




which(beta_true!=0) 
sort( subdata_selection$TH_IBOSS$stage_I_prescreened_var)
sort( subdata_selection$CT_IBOSS$stage_I_prescreened_var)
sort( subdata_selection$SuSiE_IBOSS$stage_I_prescreened_var)


