library(susieR)

# Purpose:
#  Produce a reduced design matrix (for subsample tagging / IBOSS, etc.)
#  via different variable-reduction methods.


# FUNCTION 1
# TH function_variable_reduction_Thresholding ############################################################
# Marginal p-value thresholding or top-k by smallest p-value


# X = X_train 
# 
# y = y_train

pvals_univariate <- function(X, y,  intercept = "Yes") {

  if(intercept == "Yes"){
    n <- nrow(X)
    
    # center y and all columns of X
    yc <- y - mean(y)
    Xc <- sweep(X, 2, colMeans(X), "-")
    
    # plot(colMeans(Xc)) 
    # plot(colMeans(X)) 
    
    # correlation for each column: r_j = (Xc_j' yc) / sqrt((Xc_j'Xc_j)(yc'yc))
    num <- as.numeric(crossprod(yc, Xc))               # length p
    den <- sqrt(sum(yc^2) * colSums(Xc^2))             # length p
    r   <- num / den
    
    
    # r -> t -> p (df = n - 2 for univariate w/ intercept)
    t <- r * sqrt((n - 2) / (1 - r^2))
    p <- 2 * pt(-abs(t), df = n - 2)
    

    
    return(p)    
  }
    

}


function_variable_reduction_Thresholding <- function(X_train,
                                                     y_train,
                                                     subsample_percent_variable_selection_indices,
                                                     pval_cutoff = 0.05) {
  
  # dim(X_train)
  start_time <- Sys.time()
  p <- ncol(X_train)
  adj_cutoff <- pval_cutoff / p
  
  # # old 
  # pvals <- rep(NA_real_, p)
  # for (j in seq_len(p)) {
  #   fit <- lm(y_train ~ X_train[, j])
  #   pvals[j] <- summary(fit)$coefficients[2, 4]
  # }
  #   keep <- which(pvals <= adj_cutoff)

    # new
    pvals <- pvals_univariate(X_train[subsample_percent_variable_selection_indices,], y_train[subsample_percent_variable_selection_indices])
    keep <- which( pvals <= adj_cutoff)

    
    variable_selection_elapsed_time <- as.numeric( Sys.time() - start_time, units="secs" )
    
    
    
    list(
      X_reduced = X_train[, keep, drop = FALSE],
      elapsed_time = variable_selection_elapsed_time, 
      selected_cols = keep,
      method = "TH",
      pval_cutoff = pval_cutoff,
      pvals = pvals)
    
}


# FUNCTION 2
# function_variable_reduction_C_T
# Clumping + Thresholding (C+T) ############################################################

#  1) threshold by marginal p-values (or take top_k by pval)
#  2) greedy clumping by correlation among candidates

function_variable_reduction_C_T <- function(X_train,
                                            y_train,
                                            subsample_percent_variable_selection_indices,
                                            pval_cutoff = 0.05,
                                            corr_cutoff = sqrt(0.1)) {
  
  start_time <- Sys.time()
  
  p <- ncol(X_train)
  adj_cutoff <- pval_cutoff / p
  
  # # old
  # pvals <- rep(NA_real_, p)
  # for (j in seq_len(p)) {
  #   fit <- lm(y_train ~ X_train[, j])
  #   pvals[j] <- summary(fit)$coefficients[2, 4]
  # }
  # cand <- which(pvals <= adj_cutoff)


  # new
  pvals <- pvals_univariate(X_train[subsample_percent_variable_selection_indices,], y_train[subsample_percent_variable_selection_indices])
  cand <- which(pvals <= adj_cutoff)
  
  # rank candidates by strongest signal first (smallest p-value)
  cand <- cand[order(pvals[cand])]
  
  # Step 2: greedy clumping
  selected <- integer(0)
  idx = 2
  # idx = 3
  for (idx in cand) {
    if (length(selected) == 0) {
      selected <- c(selected, idx)
    } else {
      cors <- abs(cor(X_train[subsample_percent_variable_selection_indices, idx], X_train[subsample_percent_variable_selection_indices, selected, drop = FALSE]))
      if (all(cors <= corr_cutoff)) selected <- c(selected, idx)

    }
  }
  variable_selection_elapsed_time <- as.numeric( Sys.time() - start_time, units="secs" )
  
  list(
    X_reduced = X_train[, selected,  drop = FALSE],
    elapsed_time = variable_selection_elapsed_time,
    selected_cols = selected,
    pval_cutoff = pval_cutoff,
    method = "CT",
    pvals = pvals
  )
  
  
}



# FUNCTION 3
# function_variable_reduction_Susi
# Susi ############################################################

function_variable_reduction_SuSi <- function(X_train,
                                             y_train,
                                             subsample_percent_variable_selection_indices,
                                             LD_block_size, n_LD_blocks) {
  
  start_time <- Sys.time()
  pval_cutoff_block <- 0.05 / ncol(X_train)
  selected <- integer(0)
  
  for (b in 1:n_LD_blocks) {
    # b = 1
    st <- (b - 1) * LD_block_size + 1
    en <- min(b * LD_block_size, ncol(X_train))
    X_train_block <- X_train[subsample_percent_variable_selection_indices, st:en, drop = FALSE]
    
    
    pvals_block <- pvals_univariate(X_train_block, y_train[subsample_percent_variable_selection_indices])
    if (any(pvals_block <= pval_cutoff_block)) {
      keep_local <- which(pvals_block <= 0.05)
      # fitSusi <- susie( as.matrix( X_train_block[, keep_local, drop = FALSE]), y_train[subsample_percent_variable_selection_indices], 
      #                        L=  min(sum(pvals_block <= pval_cutoff_block),    5, ncol(X_train)),  max_iter = 100)
      
      
      X_susie <- as.matrix( X_train_block[, keep_local, drop = FALSE])
      storage.mode(X_susie) <- "double"
      y_susie <- as.numeric(y_train[subsample_percent_variable_selection_indices])
      
      
      fitSusi <- susie(X_susie , y_susie , 
                       L=  min(length(keep_local),  10, ncol(X_train)),  
                       max_iter = 100)
      
      selected <- c(selected, st - 1 + keep_local[which(fitSusi$pip > 0.5)] ) }
  }
  
  variable_selection_elapsed_time <- as.numeric( Sys.time() - start_time, units="secs" )
  
  list(
    X_reduced = X_train[, selected, drop = FALSE],
    elapsed_time = variable_selection_elapsed_time,
    selected_cols = selected,
    method = "SuSi")
}


# Perfect VS ############################################################


function_variable_reduction_perfect <- function(X_train, beta_one_true) {
  
  start_time <- Sys.time()

    
  selected <- which(beta_one_true != 0 )
  
  variable_selection_elapsed_time <- as.numeric( Sys.time() - start_time, units="secs" )
  
  
  
  list(
    X_reduced = X_train[, selected, drop = FALSE],
    elapsed_time = variable_selection_elapsed_time,
    selected_cols = selected,
    method = "Perfect")
}


















########################### Evaluation variable selection ########################### 
function_variable_reduction_evaluate_selection <- function(selected_cols,
                                                           beta_true,
                                                           eps = 1e-12) {
  p <- length(beta_true)
  
  truth <- abs(beta_true) > eps
  
  
  pred  <- rep(FALSE, p)
  pred[selected_cols] <- TRUE
  
  TP <- sum(pred & truth)
  FP <- sum(pred & !truth)
  TN <- sum(!pred & !truth)
  FN <- sum(!pred & truth)
  
  TPR <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)   # recall / sensitivity/ Power 
  FPR <- if ((FP + TN) == 0) NA_real_ else FP / (FP + TN)
  PPV <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)   # precision
  FDR <- if ((TP + FP) == 0) NA_real_ else FP / (TP + FP)
  F1  <- if (is.na(PPV) || is.na(TPR) || (PPV + TPR) == 0) NA_real_ else 2 * PPV * TPR / (PPV + TPR)
  
  list(
    TP = TP, FP = FP, TN = TN, FN = FN,
    TPR = TPR, FPR = FPR, PPV = PPV, FDR = FDR, F1 = F1,
    num_selected = sum(pred),
    num_true_nonzero = sum(truth)
  )
}

