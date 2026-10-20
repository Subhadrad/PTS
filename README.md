PTS: Prescreen then Subsample

Usage: 
Download and extract this repository, keep all R scripts together, and set your R working directory to that folder.

Install the R packages loaded by the scripts before running them.
install.packages(c(  "susieR", "pracma", "MASS", "BEDMatrix", "data.table","glmnet", "doParallel", "doRNG", "foreach", "bigsnpr", "dplyr", "ggplot2", "tidyr", "cowplot", "latex2exp", "ggbreak", "tibble"))


A. Using PTS on Your Own Data

1. Open "PTS_Wrapper_and_Illustrations.R" for an example of running PTS and inspecting its output.

2. Load the PTS functions using source("Function_PTS.R"). The main function is function_select_PTS_subdata()
  
       Key Inputs:  X (covariate matrix), y (response vector), percentage of data to be used for prescreening, and subsample size.
                 i. TH: specify threshold (alpha; Bonferroni correction is applied), or top_q to select the top-ranked variables,
                ii. CT: additionally specify clumping_threshold, the absolute-correlation cutoff,
               iii. SuSiE: specify L (single effects per region), ld_regions, pip_cutoff (PIP threshold), and include_CS (TRUE adds the highest-PIP representative from each credible set).

                    See PTS_Illustrations_Examples.R for example usage for parameter settings.


B. Reproduce the simulation study results reported in Section 3 of our paper titled "Prescreen, Then Subsample: Efficient model training for large-scale biomedical data with correlated covariates" available on arXiv (arXiv:2610.21916).
 
The main script is "Main PTS Simulation Study.R".
source("Main PTS Simulation Study.R")

1. Select the desired setup at the beginning of "Main PTS Simulation Study.R", following Section 3 of the paper.
2. Run source("Main PTS Simulation Study.R") to generate and save the simulation results.
3. Set the matching results-file path and plotting settings in "Plot-Set 1.R", then run source("Plot-Set 1.R") to reproduce Figures 2 and 3.

