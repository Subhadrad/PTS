PTS Simulation Study

The main script, "Main PTS Simulation Study.R", can be used to reproduce the simulation study results reported in Section 3.

Usage

1. Download all R scripts in this repository and keep them in the same directory.
2. Select the desired simulation setup at the beginning of "Main PTS Simulation Study.R", using the settings specified in Section 3.
3. Run the main simulation script to generate the simulation data and save them in the data folder:

   source("Main PTS Simulation Study.R")

4. After the simulation has finished, run "Plot-Set 1.R" to use the generated data to reproduce Figures 2 and 3:

   source("Plot-Set 1.R")
