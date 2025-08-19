# rt_estimation
Rt estimation using local DSA

1. "rt_est_functions.R" contains 2 functions: 1 to simulate data given parameters and estimate Rt using DSA and Cori methods, and the other to plot the results against the true Rt value.
2. "data" folder contains all simulated data sets and results from the Rt estimation used in the manuscript (using rt_est function).
3. "rt_plots_file.R" contains all the code necessary to plot results from results in the data folder. Run this script to create all figures in the manuscript.
4. "run_all_rt.R" contains the code to recreate the files in the data folder. Please note that running this can take several days, which is why the csvs are provided in the data folder.
5. "localDSA_functions.R", "simulation.R", and "util2.R" are necessary only for rerunning the estimation procedure. "simulation.R" and "util2.R" are scripts from https://github.com/cobeylab/Rt_estimation/tree/master/code used to estimate Rt using the Corey method.
