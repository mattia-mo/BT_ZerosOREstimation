# Code for the simulation and analysis part of the Bachelor's Thesis
## Title: How the handling of zeros affects the performance of odds ratio estimation methods: A metascientific study

This repository is containing all needed files to conduct a simulation of 2x2 tables and their Odds Ratio estimations as well as the analysis of these estimations as it is described in the thesis. The structure is inspired by T.P. Morris based on the paper of "Using simulation studies to evaluate statistical methods"(2019). The code is implemented according to the model of Morris' code, which can be found [here](https://github.com/tpmorris/simtutorial/tree/master/R).

## Order of running the files
First the `DataGenerating.R` file needs to be runned where the tables and their estimations are computed. As a result of that file the file `estimates.rds` is created that contains all simulated data. That is number of the scenario, parameters of the scenario, all tables, the OR estimations of the methods and their respective confidence intervals. With the random number stored in `states.rds` a single table can be recreated.
Afterwards the analysis follows in the `EstimationAnalysis.R` file. This requires the `estimates.rds` file which needs to be loaded there.

In the plots folder all figures used in the thesis can be found.
