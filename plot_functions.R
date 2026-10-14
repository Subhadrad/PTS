library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(latex2exp)
library(ggbreak)
# install.packages("")



function_plots_xaxis_k=   function(data, set_n_train, set_LD_block_size, set_n_LD_blocks,
                                    set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR, 
                                    set_beta_type, set_subsample_percent_variable_selection,  variable = "r_sq_test_LASSO_Fullsample" , variable_like = "r_sq_test_LASSO_Subsample" ,
                                    summary_like = "r_sq_test_", values_to_var = "R_Squared", 
                                    algorithm_map = c(
                                      "r_sq_test_LASSO_Fullsample_mean" = "Full Sample",
                                      "r_sq_test_LASSO_Subsample_Uniform_mean" = "Uniform Subsample",
                                      "r_sq_test_LASSO_Subsample_TH_mean" = "TH IBOSS Subsample",
                                      "r_sq_test_LASSO_Subsample_CT_mean" = "CT IBOSS Subsample",
                                      "r_sq_test_LASSO_Subsample_SuSi_mean" = "SuSiE IBOSS Subsample",
                                      "r_sq_test_LASSO_Subsample_Perfect_mean" = "Oracle IBOSS Subsample"
                                    ), y_axis_label = "r_sq",
                                    x_axis_var = "k", file_name_start =   expression(r^2), legend_position = "bottom", 
                                    plot_title = set_beta_type , y_axis_break_range = NULL, y_axis_break_scales = 0.3, ymin = 0  )  {
  
  
  
  # All group varibles 
  # c("n_train", "n_test",
  #   "LD_block_size",  "n_LD_blocks", 
  #   "DataSet_Name", 
  #   "covariate_covariance_structure", 
  #   "noise_variance", "intercept",  "LD_rho", 
  #   "k", "subsample_percent",
  #   "sparsity", "SNR", "beta_type",                                                       
  #   "num_simulations",
  #   "subsample_percent_variable_selection"   
  # )
  
  #[1] "n_train", [2] "n_test",
  #[3] "LD_block_size", [4] "n_LD_blocks", 
  #[5] "DataSet_Name", # [6] "covariate_covariance_structure", [7] "noise_variance", [8] "intercept", [9] "LD_rho", 
  #[10] "k" [11] "subsample_percent"
  #[12] "sparsity", [13] "SNR", [14] "beta_type"                                                       
  #[15] "num_simulations",
  #[16] "subsample_percent_variable_selection"
  
  group_vars_plot_x_axis_k =    c("n_train", "n_test",
                                  "LD_block_size",  "n_LD_blocks",
                                  "DataSet_Name",
                                  "covariate_covariance_structure",
                                  "noise_variance", "intercept",  "LD_rho",
                                  "k", "subsample_percent",
                                  "sparsity_corrected", "SNR", "beta_type",
                                  "num_simulations",
                                  "subsample_percent_variable_selection"
  )
  
  selected_data <- data %>%
    filter( n_train %in% set_n_train,
            LD_block_size %in% set_LD_block_size, 
            n_LD_blocks %in% set_n_LD_blocks,
            DataSet_Name %in% set_DataSet_Name,
            LD_rho %in% set_LD_rho,
            sparsity_corrected %in% set_sparsity, 
            SNR %in% set_SNR, 
            beta_type %in% set_beta_type,
            subsample_percent_variable_selection %in% set_subsample_percent_variable_selection
    ) %>%
    dplyr::select( all_of(group_vars_plot_x_axis_k),  all_of(variable), starts_with(variable_like))
  
# sum( selected_data$r_sq_test_LASSO_Subsample_Uniform) 
  
  # selected_data$el  
  summary_data <- selected_data %>%
    group_by( 
      across(all_of(group_vars_plot_x_axis_k))  ) %>%
    summarise(   across(starts_with(summary_like), ~mean(.x, na.rm = TRUE) ,
                        .names = "{.col}_mean"),  .groups = "drop"   )
  
  long_data_summary <- pivot_longer(summary_data, 
                                    cols = starts_with(summary_like), 
                                    names_to = "Algorithm", 
                                    values_to = values_to_var)
  # names(long_data_summary)
  # dim(long_data_summary)
  # unique(long_data_summary$Algorithm) 
  long_data_summary <- long_data_summary %>%
    mutate(
      Algorithm = ifelse(
        Algorithm %in% names(algorithm_map),
        algorithm_map[Algorithm],
        Algorithm
      )
    )
  # unique(long_data_summary$Algorithm) 
  
  
  # Plot Parameters
  # Legend Parameters
  l1 = 25
  l1_1 = 10
  l1_2 = 20
  
  l2 = 0.6
  l3 = 1.8
  
  
  
  algorithm_style <- tibble::tribble(
    ~Algorithm,                         ~shape, ~linetype, ~alpha, ~linewidth, ~color,
    "Full Sample",                      NA,     "solid",   1.0,    1.5,        "#66C2A5",
    "Uniform Subsample",                2,     "dashed",   1.0,    1.0,        "red",
    "TH IBOSS Subsample",     4,      "dashed",  1.0,    1.0,        "#0072B2",
    "CT IBOSS Subsample",            8,      "dashed",  1.0,    1.0,        "#E69F00",
    "SuSiE IBOSS Subsample",            17,     "dashed",  1.0,    1.0,        "#9467BD",
    "Oracle IBOSS Subsample",      20,     "dashed",   0.7,    0.5,        "#013220"
  )
  
  shape_values <- setNames(algorithm_style$shape, algorithm_style$Algorithm)
  linetype_values <- setNames(algorithm_style$linetype, algorithm_style$Algorithm)
  alpha_values <- setNames(algorithm_style$alpha, algorithm_style$Algorithm)
  linewidth_values <- setNames(algorithm_style$linewidth, algorithm_style$Algorithm)
  color_values <- setNames(algorithm_style$color, algorithm_style$Algorithm)
  
  
  # shape_values <- c(
  #   "Full Sample"  = NA,
  #   "Uniform Subsample" = 2,
  #   "Thresholding IBOSS Subsample" = 4,
  #   "C + T IBOSS Subsample" = 8,
  #   "SuSi IBOSS Subsample" = 17,
  #   "Perfect(VS) IBOSS Subsample" = 20
  # )
  
  # linetype_values <- c(
  #   "Full Sample"  = "solid",
  #   "Uniform Subsample" = "dashed",
  #   "Thresholding IBOSS Subsample" = "dashed",
  #   "C + T IBOSS Subsample" = "dashed",
  #   "SuSi IBOSS Subsample" = "dashed",
  #   "Perfect(VS) IBOSS Subsample" = "dashed"
  # )
  
  # color_values <- c(
  #   "Full Sample"  ="#66C2A5",
  #   "Uniform Subsample" = "red",
  #   "Thresholding IBOSS Subsample" = "#0072B2",
  #   "C + T IBOSS Subsample" = "#E69F00",
  #   "SuSi IBOSS Subsample" = "#56B4E9",
  #   "Perfect(VS) IBOSS Subsample" =   "#013220")
  
  # y-axis upper limt mark  
  
  y_breaks <- sort(unique(c(
    pretty(long_data_summary[[values_to_var]]), 
    max(long_data_summary[[values_to_var]]), min(long_data_summary[[values_to_var]]) )))
  
  y_breaks <- pretty(long_data_summary[[values_to_var]],8) 
  # plot_name  
  # plot_title = paste0("p = ", unique(long_data_summary$LD_block_size) * unique(long_data_summary$n_LD_blocks) ) 
  
  
  # plot_file_name  
  plot_file_name <- paste(
    "n", format(unique(long_data_summary$n_train), scientific = FALSE, trim = TRUE),
    "p", unique(long_data_summary$LD_block_size) * unique(long_data_summary$n_LD_blocks),
    "bl", unique(long_data_summary$n_LD_blocks),
    "rho", gsub("\\.", "dt", unique(long_data_summary$LD_rho)) ,
    "s", gsub("\\.", "dt", unique(long_data_summary$sparsity_corrected)),
    "nu", gsub("\\.", "dt",unique(long_data_summary$SNR)),
    unique(long_data_summary$beta_type),
    "c", unique(long_data_summary$subsample_percent_variable_selection),
    substr(legend_position, 1, 1), 
    sep = "_"
  )
  
  # 
  
  legend_order <- c(
    "Full Sample",
    "TH IBOSS Subsample",
    "Oracle IBOSS Subsample",
    "CT IBOSS Subsample",
    "Uniform Subsample",
    "SuSiE IBOSS Subsample"
    )
  
  PlotCase <- 
    ggplot(long_data_summary,
           aes(x = .data[[x_axis_var]],
               y = .data[[values_to_var]] ,
               color   = Algorithm,    
               linetype = Algorithm,
               shape    = Algorithm,
               alpha    = Algorithm,    
               linewidth = Algorithm)) +
    geom_line(na.rm = TRUE) +
    geom_point(size = 3, na.rm = TRUE) +
    scale_alpha_manual( values = alpha_values, guide = "none" )  +
    scale_linewidth_manual(values = linewidth_values,
                           breaks = legend_order, guide = "none")+
    scale_shape_manual(
      values       = shape_values,
      breaks = legend_order,
      na.translate = FALSE
    ) +
    scale_linetype_manual(values = linetype_values,
                          breaks = legend_order) +
    scale_color_manual(                   
      values = color_values  ,
      breaks = legend_order                     
    ) +
    
    scale_x_continuous(breaks = sort(unique(long_data_summary[[x_axis_var]])) ) +
    labs(
      title    = plot_title,
      x        = x_axis_var,
      y        = y_axis_label,
      color    = "Algorithm",
      linetype = "Algorithm",
      shape    = "Algorithm"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "black", size = 0.6),
      axis.text         = element_text(size = l1),
      axis.title        = element_text(size = l1, face = "plain"),
      plot.title        = element_text(size = l1_2, face = "bold", hjust = 0.5),
      legend.title      = element_blank(),
      legend.text       = element_text(size = l1_1),
      legend.key        = element_rect(fill = "transparent"),
      legend.key.size   = unit(l2, 'cm'),
      legend.key.width  = unit(l3, 'cm'),
      legend.position   = legend_position ,
      legend.direction  = "horizontal",
      legend.box        = "horizontal",
      panel.border      = element_blank(),
    )
  
  
  # print(PlotCase)
  
  
  
  
  
  
  # y_axis_break_range = c(140, 1250)
  if (!is.null(y_axis_break_range)) {
    
      y_values <- long_data_summary[[values_to_var]]
      lower_values <- y_values[y_values <= y_axis_break_range[1]]
      upper_values <- y_values[y_values >= y_axis_break_range[2]]
      lower_breaks <- pretty( range(c(0, lower_values), na.rm = TRUE),  n = 5 )
      upper_breaks <- pretty(        range(upper_values, na.rm = TRUE),        n = 2      )
      y_breaks <- sort(unique(c(lower_breaks, upper_breaks)))
      
    
    
    PlotCase2 <- PlotCase    + 
      scale_y_continuous(breaks = y_breaks , limits = c(ymin, max(y_breaks)))  +
      ggbreak::scale_y_break(        y_axis_break_range,        scales = y_axis_break_scales ,
                                     symbol = "slash"      ) +
      theme(
        axis.text.y.right        = element_blank(),
        axis.ticks.y.right       = element_blank(),
        axis.line.y.right        = element_blank(),
        axis.title.y.right       = element_blank(),
        axis.ticks.length.y.right = unit(0, "pt")
      )
    print(PlotCase2)
    

  } else {
    PlotCase2 <- PlotCase +     scale_y_continuous(breaks = y_breaks) 
    print(PlotCase2)
  }
  
  
  ggsave(paste(getwd(), "/Plots/", file_name_start, "_",    plot_file_name, ".jpeg", sep= ""),
         plot = PlotCase2, width = 10, height = 6, units = "in")
  
  # legend_only <- cowplot::get_legend(PlotCase)  
  
  
  library(cowplot)
  library(grid)

  p <- PlotCase + theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

  legs <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
  leg <- legs[[which(!sapply(legs, inherits, "zeroGrob"))[1]]]

  png("by k legend_only.png", width = 4000, height = 350, res = 500)
  grid.newpage()
  grid.draw(leg)
  dev.off()

  
  
}


function_plots_xaxis_c =   function(data, set_n_train, set_k,  set_LD_block_size, set_n_LD_blocks,
                                    set_DataSet_Name, set_LD_rho,  set_sparsity, set_SNR, 
                                    set_beta_type, 
                                    variable = "r_sq_test_LASSO_Fullsample" , variable_like = "r_sq_test_LASSO_Subsample" ,
                                    summary_like = "r_sq_test_", values_to_var = "R_Squared", 
                                    algorithm_map = c(
                                      "r_sq_test_LASSO_Fullsample_mean" = "Full Sample",
                                      "r_sq_test_LASSO_Subsample_Uniform_mean" = "Uniform Subsample",
                                      "r_sq_test_LASSO_Subsample_TH_mean" = "Thresholding IBOSS Subsample",
                                      "r_sq_test_LASSO_Subsample_CT_mean" = "C + T IBOSS Subsample",
                                      "r_sq_test_LASSO_Subsample_SuSi_mean" = "SuSi IBOSS Subsample",
                                      "r_sq_test_LASSO_Subsample_Perfect_mean" = "Perfect(VS) IBOSS Subsample"
                                    ), y_axis_label = "r_sq",
                                    x_axis_var = "c", file_name_start =   expression(r^2),
                                    legend_position = "bottom",
                                     plot_title = set_beta_type)  {
  
  # All group varibles 
  # c("n_train", "n_test",
  #   "LD_block_size",  "n_LD_blocks", 
  #   "DataSet_Name", 
  #   "covariate_covariance_structure", 
  #   "noise_variance", "intercept",  "LD_rho", 
  #   "k", "subsample_percent",
  #   "sparsity", "SNR", "beta_type",                                                       
  #   "num_simulations",
  #   "subsample_percent_variable_selection"   
  # )
  
  #[1] "n_train", [2] "n_test",
  #[3] "LD_block_size", [4] "n_LD_blocks", 
  #[5] "DataSet_Name", # [6] "covariate_covariance_structure", [7] "noise_variance", [8] "intercept", [9] "LD_rho", 
  #[10] "k" [11] "subsample_percent"
  #[12] "sparsity", [13] "SNR", [14] "beta_type"                                                       
  #[15] "num_simulations",
  #[16] "subsample_percent_variable_selection"
  
  # variable = "r_sq_test_LASSO_Fullsample"
  # variable_like = "r_sq_test_LASSO_Subsample"
  # summary_like = "r_sq_test_"
  # values_to_var = "R_Squared"
  # algorithm_map = c(
  #   "r_sq_test_LASSO_Fullsample_mean" = "Full Sample",
  #   "r_sq_test_LASSO_Subsample_Uniform_mean" = "Uniform Subsample",
  #   "r_sq_test_LASSO_Subsample_TH_mean" = "Thresholding IBOSS Subsample",
  #   "r_sq_test_LASSO_Subsample_CT_mean" = "C + T IBOSS Subsample",
  #   "r_sq_test_LASSO_Subsample_SuSi_mean" = "SuSi IBOSS Subsample",
  #   "r_sq_test_LASSO_Subsample_Perfect_mean" = "Perfect(VS)- IBOSS Subsample"
  # )
  # y_axis_label = expression(r^2)
  # x_axis_var = "c"
  # file_name_start = "R_sq"

  group_vars_plot_x_axis_c =    c("n_train", "n_test",
                                  "LD_block_size",  "n_LD_blocks",
                                  "DataSet_Name",
                                  "covariate_covariance_structure",
                                  "noise_variance", "intercept",  "LD_rho",
                                  "k", "subsample_percent",
                                  "sparsity_corrected", "SNR", "beta_type",
                                  "num_simulations",
                                  "c")
  
  # table(data$beta_type)
  selected_data <- data %>%
    filter( n_train %in% set_n_train,
            k %in% set_k,
            LD_block_size %in% set_LD_block_size, 
            n_LD_blocks %in% set_n_LD_blocks,
            DataSet_Name %in% set_DataSet_Name,
            LD_rho %in% set_LD_rho,
            sparsity_corrected %in% set_sparsity, 
            SNR %in% set_SNR, 
            beta_type %in% set_beta_type
    ) %>%
    dplyr::select( all_of(group_vars_plot_x_axis_c),  all_of(variable), starts_with(variable_like))
  
  summary_data <- selected_data %>%
    group_by( 
      across(all_of(group_vars_plot_x_axis_c))  ) %>%
    summarise(   across(starts_with(summary_like), ~mean(.x) , .names = "{.col}_mean"),  .groups = "drop"   )
  
  long_data_summary <- pivot_longer(summary_data, 
                                    cols = starts_with(summary_like), 
                                    names_to = "Algorithm", 
                                    values_to = values_to_var)
  
  # unique(long_data_summary$Algorithm) 
  long_data_summary <- long_data_summary %>%
    mutate(
      Algorithm = ifelse(
        Algorithm %in% names(algorithm_map),
        algorithm_map[Algorithm],
        Algorithm
      )
    )
  # unique(long_data_summary$Algorithm) 
  
  
  # Plot Parameters
  # Legend Parameters
  l1 = 25
  l1_1 = 10
  l2 = 0.6
  l3 = 1.8

  
  
  algorithm_style <- tibble::tribble(
    ~Algorithm,                         ~shape, ~linetype, ~alpha, ~linewidth, ~color,
    "Full Sample",                      NA,     "solid",   1.0,    1.5,        "#66C2A5",
    "Uniform Subsample",                NA,     "solid",   1.0,    1.0,        "red",
    "TH IBOSS Subsample",     4,      "dashed",  1.0,    1.0,        "#0072B2",
    "CT IBOSS Subsample",            8,      "dashed",  1.0,    1.0,        "#E69F00",
    "SuSiE IBOSS Subsample",            17,     "dashed",  1.0,    1.0,        "#9467BD",
    "Oracle IBOSS Subsample",      NA,     "solid",   0.7,    0.5,        "#013220"
  )
  
  shape_values <- setNames(algorithm_style$shape, algorithm_style$Algorithm)
  linetype_values <- setNames(algorithm_style$linetype, algorithm_style$Algorithm)
  alpha_values <- setNames(algorithm_style$alpha, algorithm_style$Algorithm)
  linewidth_values <- setNames(algorithm_style$linewidth, algorithm_style$Algorithm)
  color_values <- setNames(algorithm_style$color, algorithm_style$Algorithm)
    
  # y-axis upper limt mark  
  y_breaks <- sort(unique(c(
    pretty(long_data_summary[[values_to_var]]), 
    max(long_data_summary[[values_to_var]]), min(long_data_summary[[values_to_var]]) )))
  
  y_breaks <- pretty(long_data_summary[[values_to_var]],8) 
  
  # plot_name  
  # plot_title = paste0("p = ", unique(long_data_summary$LD_block_size) * unique(long_data_summary$n_LD_blocks), " k = ", unique(long_data_summary$k) ) 
  # plot_title = set_beta_type
  
  # plot_file_name  
  plot_file_name <- paste(
    "n", format(unique(long_data_summary$n_train), scientific = FALSE, trim = TRUE),
    "p", unique(long_data_summary$LD_block_size) * unique(long_data_summary$n_LD_blocks),
    "bl", unique(long_data_summary$n_LD_blocks),
    "rho", gsub("\\.", "dt", unique(long_data_summary$LD_rho)) ,
    "s", gsub("\\.", "dt", unique(long_data_summary$sparsity_corrected)),
    "SNR", gsub("\\.", "dt",unique(long_data_summary$SNR)),
    unique(long_data_summary$beta_type),
    "k", unique(long_data_summary$k),  
    substr(legend_position, 1, 1),
    sep = "_"
  )
  
  # 
  
  legend_order <- c(
    "Full Sample",
    "Oracle IBOSS Subsample",
    "TH IBOSS Subsample",
    "CT IBOSS Subsample",
    "SuSiE IBOSS Subsample",
    "Uniform Subsample"  
  )
  
  PlotCase <- 
    ggplot(long_data_summary,
           aes(x = .data[[x_axis_var]],
               y = .data[[values_to_var]] ,
               color   = Algorithm,    
               linetype = Algorithm,
               shape    = Algorithm,
               alpha    = Algorithm,    
               linewidth = Algorithm)) +
    geom_line(na.rm = TRUE) +
    geom_point(size = 3, na.rm = TRUE) +
    scale_alpha_manual( values = alpha_values, guide = "none" )  +
    scale_linewidth_manual(values = linewidth_values,
                           breaks = legend_order, guide = "none")+
    scale_shape_manual(
      values       = shape_values,
      breaks = legend_order,
      na.translate = FALSE
    ) +
    scale_linetype_manual(values = linetype_values,
                          breaks = legend_order) +
    scale_color_manual(                   
      values = color_values ,
      breaks = legend_order                    
    ) +
    
    scale_x_continuous(trans = "reverse", breaks = sort(unique(long_data_summary[[x_axis_var]])) ) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      title    = plot_title,
      x        = x_axis_var,
      y        = y_axis_label,
      color    = "Algorithm",
      linetype = "Algorithm",
      shape    = "Algorithm"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "black", size = 0.6),
      axis.text         = element_text(size = l1),
      axis.title        = element_text(size = l1, face = "plain"),
      plot.title        = element_text(size = l1, face = "bold", hjust = 0.5),
      legend.title      = element_blank(),
      legend.text       = element_text(size = l1_1),
      legend.key        = element_rect(fill = "transparent"),
      legend.key.size   = unit(l2, 'cm'),
      legend.key.width  = unit(l3, 'cm'),
      legend.position   = legend_position,
      legend.direction  = "horizontal",
      legend.box        = "horizontal",
      panel.border      = element_blank(),
    )
  
  print(PlotCase)
  ggsave(paste(getwd(), "/PlotsSet2/", file_name_start, "_",    plot_file_name, ".jpeg", sep= ""),
         plot = PlotCase, width = 12, height = 7.2, units = "in", dpi = 300)
  
  

library(cowplot)
library(grid)

clean_legend_labels <- function(x) {
  x <- gsub("\\s*IBOSS Subsample\\s*$", "", x)
  x <- gsub("Perfect\\(VS\\)\\s*-\\s*$", "Perfect(VS)", x)
  x
}

# p <- PlotCase + theme(
#   legend.position = "bottom",
#   legend.title = element_blank(),
#   legend.text = element_text(size = 8)
# ) #+
#   # guides(
#   #   colour   = guide_legend(nrow = 2, byrow = TRUE),
#   #   linetype = guide_legend(nrow = 2, byrow = TRUE),
#   #   shape    = guide_legend(nrow = 2, byrow = TRUE)
#   # )

p <- PlotCase +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.spacing.x = unit(5.5, "cm"),
    legend.key.width = unit(3, "cm")
  ) +
  guides(
    color    = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1),
    shape    = guide_legend(nrow = 1)
  )

# Change labels without changing existing colors/shapes/linetypes
for (aes_name in c("colour", "color", "linetype", "shape")) {
  sc <- p$scales$get_scales(aes_name)
  if (!is.null(sc)) {
    sc$labels <- clean_legend_labels
  }
}




legs <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
leg <- legs[[which(!sapply(legs, inherits, "zeroGrob"))[1]]]

png("Initial prescreening legend_only.png", width = 4000, height = 200, res = 350)
grid.newpage()
grid.draw(leg)
dev.off()
  
}

library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(grid)
library(cowplot)
library(latex2exp)


library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(grid)
library(cowplot)
library(latex2exp)

function_plot_rsq_iboss_vs_uniform_by_k <- function(
    data, set_n_train, set_LD_block_size, set_n_LD_blocks,
    set_DataSet_Name, set_LD_rho, set_sparsity, set_SNR,
    set_beta_type, set_subsample_percent_variable_selection = 50,
    y_axis_label = TeX("${r}^{ 2}_{test}$"), x_axis_var = "k",
    file_name_start = "R_sq_IBOSS_vs_Uniform",
    legend_position = "bottom", plot_title = NULL, save_legend = FALSE) {
  
  group_vars_plot_x_axis_k <- c(
    "n_train", "n_test", "LD_block_size", "n_LD_blocks", "DataSet_Name",
    "covariate_covariance_structure", "noise_variance", "intercept", "LD_rho",
    "k", "subsample_percent", "sparsity_corrected", "SNR", "beta_type",
    "num_simulations", "subsample_percent_variable_selection"
  )
  
  variable <- c(
    "r_sq_test_LASSO_Subsample_Uniform",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_TH",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_CT",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_SuSi",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_Uniform_TH",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_Uniform_CT",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_Uniform_SuSi"
  )
  
  algorithm_map <- c(
    "r_sq_test_LASSO_ScreenedVariables_Subsample_TH_mean" = "TH IBOSS Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_Uniform_TH_mean" = "TH Uniform Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_CT_mean" = "CT IBOSS Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_Uniform_CT_mean" = "CT Uniform Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_SuSi_mean" = "SuSiE IBOSS Subsample",
    "r_sq_test_LASSO_ScreenedVariables_Subsample_Uniform_SuSi_mean" = "SuSiE Uniform Subsample",
    "r_sq_test_LASSO_Subsample_Uniform_mean" = "Uniform Subsample"
  )
  
  legend_order <- c(
    "TH IBOSS Subsample", "TH Uniform Subsample",
    "CT IBOSS Subsample", "CT Uniform Subsample",
    "SuSiE IBOSS Subsample", "SuSiE Uniform Subsample",
    "Uniform Subsample"
  )
  
  selected_data <- data %>%
    filter(
      n_train %in% set_n_train, LD_block_size %in% set_LD_block_size,
      n_LD_blocks %in% set_n_LD_blocks, DataSet_Name %in% set_DataSet_Name,
      LD_rho %in% set_LD_rho, sparsity_corrected %in% set_sparsity,
      SNR %in% set_SNR, beta_type %in% set_beta_type,
      subsample_percent_variable_selection %in% set_subsample_percent_variable_selection
    ) %>%
    dplyr::select(all_of(group_vars_plot_x_axis_k), all_of(variable))
  
  summary_data <- selected_data %>%
    group_by(across(all_of(group_vars_plot_x_axis_k))) %>%
    summarise(across(all_of(variable), ~mean(.x, na.rm = TRUE),
                     .names = "{.col}_mean"), .groups = "drop")
  
  long_data_summary <- summary_data %>%
    pivot_longer(cols = ends_with("_mean"),
                 names_to = "Algorithm", values_to = "R_Squared") %>%
    mutate(
      Algorithm = ifelse(Algorithm %in% names(algorithm_map),
                         algorithm_map[Algorithm], Algorithm),
      Algorithm = factor(Algorithm, levels = legend_order)
    )
  
  l1 = 25
  l1_1 = 10
  l2 = 0.6
  l3 = 2
  
  algorithm_style <- tibble::tribble(
    ~Algorithm,                 ~shape, ~linetype, ~alpha, ~linewidth, ~color,
    "TH IBOSS Subsample",                 4,      "dashed",  1.0,    1.0,        "#0072B2",
    "TH Uniform Subsample",     4,      "solid",   1.0,    1.0,        "#0072B2",
    "CT IBOSS Subsample",                 8,      "dashed",  1.0,    1.0,        "#E69F00",
    "CT Uniform Subsample",     8,      "solid",   1.0,    1.0,        "#E69F00",
    "SuSiE IBOSS Subsample",              17,     "dashed",  1.0,    1.0,        "#9467BD",
    "SuSiE Uniform Subsample",  17,     "solid",   1.0,    1.0,        "#9467BD",
    "Uniform Subsample",        2,      "solid",   1.0,    1.0,        "red"
  )
  
  shape_values     <- setNames(algorithm_style$shape, algorithm_style$Algorithm)
  linetype_values  <- setNames(algorithm_style$linetype, algorithm_style$Algorithm)
  alpha_values     <- setNames(algorithm_style$alpha, algorithm_style$Algorithm)
  linewidth_values <- setNames(algorithm_style$linewidth, algorithm_style$Algorithm)
  color_values     <- setNames(algorithm_style$color, algorithm_style$Algorithm)
  
  y_breaks <- pretty(long_data_summary$R_Squared, 8)
  
  plot_file_name <- paste(
    "n", format(unique(long_data_summary$n_train), scientific = FALSE, trim = TRUE),
    "p", unique(long_data_summary$LD_block_size) * unique(long_data_summary$n_LD_blocks),
    "bl", unique(long_data_summary$n_LD_blocks),
    "rho", gsub("\\.", "dt", unique(long_data_summary$LD_rho)),
    "s", gsub("\\.", "dt", unique(long_data_summary$sparsity_corrected)),
    "SNR", gsub("\\.", "dt", unique(long_data_summary$SNR)),
    unique(long_data_summary$beta_type),
    "c", unique(long_data_summary$subsample_percent_variable_selection),
    substr(legend_position, 1, 1), sep = "_"
  )
  
  PlotCase <- ggplot(
    long_data_summary,
    aes(x = .data[[x_axis_var]], y = R_Squared, color = Algorithm,
        linetype = Algorithm, shape = Algorithm, alpha = Algorithm,
        linewidth = Algorithm)
  ) +
    geom_line(na.rm = TRUE) +
    geom_point(size = 3, na.rm = TRUE) +
    scale_alpha_manual(values = alpha_values, guide = "none") +
    scale_linewidth_manual(values = linewidth_values,
                           breaks = legend_order, guide = "none") +
    scale_shape_manual(values = shape_values,
                       breaks = legend_order, na.translate = FALSE) +
    scale_linetype_manual(values = linetype_values, breaks = legend_order) +
    scale_color_manual(values = color_values, breaks = legend_order) +
    scale_x_continuous(breaks = sort(unique(long_data_summary[[x_axis_var]]))) +
    scale_y_continuous(breaks = y_breaks) +
    labs(title = plot_title, x = x_axis_var, y = y_axis_label) +
    guides(
      color = guide_legend(nrow = 2, byrow = FALSE),
      linetype = guide_legend(nrow = 2, byrow = FALSE),
      shape = guide_legend(nrow = 2, byrow = FALSE)
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.text = element_text(size = l1),
      axis.title = element_text(size = l1, face = "plain"),
      plot.title = element_text(size = l1, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = l1_1),
      legend.key = element_rect(fill = "transparent"),
      legend.key.size = unit(l2, "cm"),
      legend.key.width = unit(l3, "cm"),
      legend.position = legend_position,
      legend.direction = "horizontal",
      legend.box = "horizontal"
    )
  
  print(PlotCase)
  
  ggsave(
    paste0(getwd(), "/Sensitivity/", file_name_start, "_", plot_file_name, ".jpeg"),
    plot = PlotCase, width = 12, height = 7.2, units = "in", dpi = 300
  )
  
  # Legend-only plot
  if (save_legend) {
    p_legend <- PlotCase + theme(
      legend.position = "bottom", plot.title = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank()
    )
    
    legs <- cowplot::get_plot_component(p_legend, "guide-box", return_all = TRUE)
    leg <- legs[[which(!sapply(legs, inherits, "zeroGrob"))[1]]]
    
    legend_plot <- cowplot::ggdraw() +
      cowplot::draw_grob(leg, x = 0, y = 0, width = 1, height = 1)
    
    ggsave(
      paste0(getwd(), "/Sensitivity/IBOSS_vs_Uniform_legend_only.png"),
      plot = legend_plot, width = 11, height = 1.2, units = "in", dpi = 300
    )
  }
  
  invisible(PlotCase)
}








library(ggplot2)
library(dplyr)
library(grid)
library(latex2exp)
library(cowplot)

function_plot_TH_topq_vs_CV <- function(
    data, set_n_train, set_LD_block_size, set_n_LD_blocks,
    set_DataSet_Name, set_LD_rho, set_sparsity, set_SNR, set_beta_type,
    set_subsample_percent_variable_selection = 50,
    top_q_lower_cap = -Inf,
    top_q_upper_cap = Inf,
    y_axis_label = expression(r^2),
    x_axis_label = "q",
    file_name_start = "R_sq_TH_TopQ_vs_CV",
    legend_position = "bottom",
    plot_title = NULL,
    save_legend = FALSE) {
  
  # Select data
  selected_data <- data %>%
    filter(
      n_train %in% set_n_train,
      LD_block_size %in% set_LD_block_size,
      n_LD_blocks %in% set_n_LD_blocks,
      DataSet_Name %in% set_DataSet_Name,
      LD_rho %in% set_LD_rho,
      sparsity_corrected %in% set_sparsity,
      SNR %in% set_SNR,
      beta_type %in% set_beta_type,
      subsample_percent_variable_selection %in% set_subsample_percent_variable_selection
    )
  
  # TH with fixed number of retained variables
  TH_summary <- selected_data %>%
    filter(top_q >= top_q_lower_cap, top_q <= top_q_upper_cap) %>%
    group_by(top_q) %>%
    summarise(
      R_Squared = mean(r_sq_test_LASSO_ScreenedVariables_Subsample_TH, na.rm = TRUE),
      .groups = "drop"
    )
  
  # SuSiE benchmark: average number of selected variables
  CV_q <- mean(selected_data$Eval_SuSi_SignalsSelected, na.rm = TRUE)
  
  # SuSiE benchmark: average R-squared
  CV_rsq <- mean(
    selected_data$r_sq_test_LASSO_ScreenedVariables_Subsample_SuSi,
    na.rm = TRUE
  )
  
  # Y-axis breaks
  y_breaks <- pretty(c(TH_summary$R_Squared, CV_rsq), 8)
  
  # File name
  plot_file_name <- paste(
    "TopQ_n", format(unique(selected_data$n_train), scientific = FALSE, trim = TRUE),
    "p", unique(selected_data$LD_block_size) * unique(selected_data$n_LD_blocks),
    "bl", unique(selected_data$n_LD_blocks),
    "rho", gsub("\\.", "dt", unique(selected_data$LD_rho)),
    "s", gsub("\\.", "dt", unique(selected_data$sparsity_corrected)),
    "SNR", gsub("\\.", "dt", unique(selected_data$SNR)),
    unique(selected_data$beta_type),
    substr(legend_position, 1, 1),
    sep = "_"
  )
  
  # Plot
  PlotCase <- ggplot(TH_summary, aes(x = top_q, y = R_Squared)) +
    
    # TH curve
    geom_line(aes(color = "TH IBOSS Subsample", linetype = "TH IBOSS Subsample"),
              linewidth = 1) +
    
    # TH points
    geom_point(aes(color = "TH IBOSS Subsample"),
               shape = 4, size = 3, show.legend = FALSE) +
    
    # SuSiE horizontal benchmark
    geom_hline(aes(yintercept = CV_rsq,
                   color = "SuSiE IBOSS Subsample",
                   linetype = "SuSiE IBOSS Subsample"),
               linewidth = 1) +
    
    # Average number of variables retained by SuSiE
    geom_vline(xintercept = CV_q,
               color = "red", linetype = "solid", linewidth = 1,
               show.legend = FALSE) +
    
    # Colors
    scale_color_manual(
      values = c("TH IBOSS Subsample" = "#0072B2",
                 "SuSiE IBOSS Subsample" = "#9467BD"),
      breaks = c("TH IBOSS Subsample", "SuSiE IBOSS Subsample")
    ) +
    
    # Line types
    scale_linetype_manual(
      values = c("TH IBOSS Subsample" = "dashed",
                 "SuSiE IBOSS Subsample" = "solid"),
      breaks = c("TH IBOSS Subsample", "SuSiE IBOSS Subsample")
    ) +
    
    # X-axis: do not add CV_q as a labelled tick
    scale_x_continuous(
      breaks = sort(unique(TH_summary$top_q)),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    
    scale_y_continuous(breaks = y_breaks) +
    
    labs(
      title = plot_title,
      x = x_axis_label,
      y = y_axis_label,
      color = NULL,
      linetype = NULL
    ) +
    
    # Only lines in legend
    guides(
      color = guide_legend(override.aes = list(shape = NA)),
      linetype = guide_legend(override.aes = list(shape = NA))
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 25, face = "plain"),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.key = element_rect(fill = "transparent"),
      legend.key.width = unit(1.8, "cm"),
      legend.position = legend_position,
      legend.direction = "horizontal"
    )
  
  print(PlotCase)
  
  # Save plot
  ggsave(
    paste0(getwd(), "/Sensitivity/", file_name_start, "_", plot_file_name, ".jpeg"),
    plot = PlotCase,
    width = 12, height = 7.2, units = "in", dpi = 300
  )
  
  # Optional separate legend
  if (save_legend) {
    
    p_legend <- PlotCase + theme(legend.position = "bottom")
    
    legs <- cowplot::get_plot_component(
      p_legend, "guide-box", return_all = TRUE
    )
    
    leg <- legs[[which(!sapply(legs, inherits, "zeroGrob"))[1]]]
    
    legend_plot <- cowplot::ggdraw() +
      cowplot::draw_grob(leg, x = 0, y = 0, width = 1, height = 1)
    
    ggsave(
      paste0(getwd(), "/Sensitivity/TH_TopQ_vs_CV_legend_only.png"),
      plot = legend_plot,
      width = 6, height = 0.7, units = "in", dpi = 300
    )
  }
  
  invisible(PlotCase)
}