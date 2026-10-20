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
                                                     pval_cutoff = 0.05, 
                                                     top_q = NULL ) {
  
  # dim(X_train)
  start_time <- Sys.time()
  p <- ncol(X_train)
  
  
  # # old 
  # pvals <- rep(NA_real_, p)
  # for (j in seq_len(p)) {
  #   fit <- lm(y_train ~ X_train[, j])
  #   pvals[j] <- summary(fit)$coefficients[2, 4]
  # }
  #   keep <- which(pvals <= adj_cutoff)
  
  # new
  pvals <- pvals_univariate(X_train[subsample_percent_variable_selection_indices,], y_train[subsample_percent_variable_selection_indices])
  
  
  if (is.null(top_q)) {
    keep <- which(pvals <= pval_cutoff / p)
  } else {
    correlations <- abs(cor(X_train[subsample_percent_variable_selection_indices,], 
                            y_train[subsample_percent_variable_selection_indices]))
    keep <- order(correlations, decreasing = TRUE)[seq_len(top_q)]
  }
  
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
  # idx = 2
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

# function_variable_reduction_SuSi <- function(X_train,
#                                             y_train,
#                                             subsample_percent_variable_selection_indices,
#                                             LD_block_size, n_LD_blocks) {
# 
# start_time <- Sys.time()
# pval_cutoff_block <- 0.05 / ncol(X_train)
# selected <- integer(0)
#     
#     for (b in 1:n_LD_blocks) {
#       # b = 1
#     st <- (b - 1) * LD_block_size + 1
#     en <- min(b * LD_block_size, ncol(X_train))
#     X_train_block <- X_train[subsample_percent_variable_selection_indices, st:en, drop = FALSE]
# 
# 
#         pvals_block <- pvals_univariate(X_train_block, y_train[subsample_percent_variable_selection_indices])
#                 if (any(pvals_block <= pval_cutoff_block)) {
#                 keep_local <- which(pvals_block <= pval_cutoff_block)
#                 # fitSusi <- susie( as.matrix( X_train_block[, keep_local, drop = FALSE]), y_train[subsample_percent_variable_selection_indices], 
#                 #                        L=  min(sum(pvals_block <= pval_cutoff_block),    5, ncol(X_train)),  max_iter = 100)
#                 
#                 
#                 X_susie <- as.matrix( X_train_block[, keep_local, drop = FALSE])
#                 storage.mode(X_susie) <- "double"
#                 y_susie <- as.numeric(y_train[subsample_percent_variable_selection_indices])
#                 
#                 
#                 fitSusi <- susie(X_susie , y_susie , 
#                                   L=  min(sum(pvals_block <= pval_cutoff_block),    10, ncol(X_train)),  
#                                     max_iter = 100)
#                 
#                 selected <- c(selected, st - 1 + keep_local[which(fitSusi$pip > 0.5)] ) }
#                              }
# 
#     variable_selection_elapsed_time <- as.numeric( Sys.time() - start_time, units="secs" )
# 
#     list(
#       X_reduced = X_train[, selected, drop = FALSE],
#       elapsed_time = variable_selection_elapsed_time,
#       selected_cols = selected,
#       method = "SuSi")
# }
# 


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



#SuSie 2 #### 
function_variable_reduction_SuSi_v2 <- function(
    X_train, y_train, subsample_percent_variable_selection_indices,
    ld_regions,
    L = 10, pval_cutoff = 0.05, pip_cutoff = 0.5, include_CS = FALSE) {
  
  start_time <- Sys.time()
  X_train <- as.matrix(X_train)
  p <- ncol(X_train)
  
  stopifnot(is.numeric(X_train), is.numeric(y_train), is.null(dim(y_train)),
            nrow(X_train) == length(y_train), p >= 1,
            length(L) == 1, is.finite(L), L >= 1, L == floor(L),
            length(pval_cutoff) == 1, pval_cutoff > 0, pval_cutoff < 1,
            length(pip_cutoff) == 1, pip_cutoff >= 0, pip_cutoff < 1, is.list(ld_regions))
  
  # Each column must belong to exactly one supplied LD region.
  region_cols <- unlist(ld_regions, use.names = FALSE)
  stopifnot(is.numeric(region_cols), all(is.finite(region_cols)),
            all(region_cols == floor(region_cols)), !anyDuplicated(region_cols),
            setequal(region_cols, seq_len(p)))
  
  # Use the supplied row indices for variable screening.
  # subsample_percent_variable_selection_indices = initial_prescreeing_random_subsample_indices
  rows <- subsample_percent_variable_selection_indices
  stopifnot(is.numeric(rows), length(rows) >= 3, all(is.finite(rows)),
            all(rows == floor(rows)), all(rows >= 1 & rows <= nrow(X_train)))
  
  a <- X_train[rows, , drop = FALSE]
  b <- y_train[rows]
  stopifnot(all(is.finite(a)), all(is.finite(b)))
  pvals <- pvals_univariate(a, b)
  selected <- integer(0)
  pip <- rep(NA_real_, p)
  
  # Explicit regions replace the calculation based on block size.
  credible_sets <- vector("list", length(ld_regions))
  
  for (region in ld_regions) {
    # region = ld_regions[[2]]
    if (!any(pvals[region] <= pval_cutoff / p)) next
    candidates <- region[pvals[region] <= pval_cutoff]
    if (!length(candidates)) next
    
    fit <- susieR::susie(a[, candidates, drop = FALSE], b,
                         L = min(L, length(candidates)), max_iter = 100)
    
    if (!isTRUE(fit$converged)) warning("A SuSiE region did not converge within 100 iterations.")
    pip[candidates] <- fit$pip
    selected <- c(selected, candidates[which(fit$pip > pip_cutoff)])
    
    if (include_CS) {
      cs <- susieR::susie_get_cs(fit, X = a[, candidates, drop = FALSE])$cs
      cs_rep <- if (length(cs)) unique(sapply(cs, function(x) x[which.max(fit$pip[x])])) else integer(0)
      selected <- union(selected, candidates[cs_rep])
    }
  }
  
  selected <- sort(selected)
  
  # Return selected columns for all rows of X_train.
  list(X_reduced = X_train[, selected, drop = FALSE], selected_cols = selected, pip = pip,
       method = "SuSiE", elapsed_time = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
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


















# FUNCTION 1b
# TH HAPNEST  ############################################################
# Marginal p-value thresholding or top-k by smallest p-value


function_variable_reduction_Thresholding_HAPNEST <- function(y,
                                                             subsample_percent_variable_selection_indices ,
                                                             pval_cutoff = 0.05) {
  
  start_time <- Sys.time()
  
  pval_cutoff_adj = pval_cutoff/p   
  
  keep_file_txt_name_location <- file.path(    geno_files_processed_location,    paste0("keep_simulation_", sim_id, ".txt")  )
  keep_df <- fam[subsample_percent_variable_selection_indices, c("FID", "IID")]
  write.table( keep_df,    file = keep_file_txt,    quote = FALSE,    row.names = FALSE,    col.names = FALSE,    sep = "\t"  )
  
  
  
  
  pheno_file_txt_name_location <- file.path(    geno_files_processed_location,    paste0("pheno_simulation_", sim_id, ".txt")  )
  assoc_out_name_location <- file.path(    geno_files_processed_location,    paste0("plink_", "simulation_", sim_id)  )
  pheno_df <- data.frame(    FID = fam$FID, IID = fam$IID,    y = y  )
  write.table(    pheno_df, file = pheno_file_txt,    quote = FALSE,    row.names = FALSE,    col.names = TRUE,    sep = "\t"  )
  # run PLINK
  # run PLINK on subset only
  system2(
    plink2_exe,
    c(
      "--bfile", HAPNEST_data_input_name_location,
      "--keep", keep_file_txt_name_location,
      "--pheno", pheno_file_txt_name_location,
      "--pheno-name", "y",
      "--glm", "allow-no-covars",
      "--out", assoc_out_name_location
    )
  )
  # read PLINK result
  plink_assoc <- read.table(paste0(assoc_out_name_location, ".y.glm.linear"),    header = TRUE,    comment.char = "",    check.names = FALSE  )
  selected_snps <- plink_assoc$ID[plink_assoc$TEST == "ADD" & plink_assoc$P < pval_cutoff_adj]
  # length(selected_snps)  
  
  # dim(plink_assoc)
  # "chr22:10571559:G:A"
  #   
  #   keep <- which(colnames(G) %in% selected_snps)
  #   
  #   length(colnames(G))
  #   length(selected_snps)
  
  G_ids <- sub("_[ACGT]$", "", colnames(G))
  
  keep <- which(G_ids %in% selected_snps)  
  
  X_reduced <- as.matrix(G[, keep, drop = FALSE])
  
  variable_selection_elapsed_time <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs"))
  
  
  
  
  
  
  
  list(
    X_reduced = X_train[, keep, drop = FALSE],
    elapsed_time = variable_selection_elapsed_time, 
    selected_cols = keep,
    method = "TH",
    pval_cutoff = pval_cutoff,
    pvals = pvals)
  
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

