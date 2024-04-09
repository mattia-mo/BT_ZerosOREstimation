# Code for the simulation and analysis part of the Bachelor's Thesis
## Title: How the handling of zeros affects the performance of odds ratio estimation methods: A metascientific study

This repository is containing all needed files to conduct a simulation of 2x2 tables and their Odds Ratio estimations as well as the analysis of these estimations as it is described in the thesis. The structure is inspired by T.P. Morris based on the paper of "Using simulation studies to evaluate statistical methods"(2019). The code is implemented according to the model of Morris' code, which can be found [here](https://github.com/tpmorris/simtutorial/tree/master/R).

## Order of running the files
First the `DataGenerating.R` file needs to be runned where the tables and their estimations are computed. As a result of that file the file `estimates.rds` is created that contains all simulated data. That is number of the scenario, parameters of the scenario, all tables, the OR estimations of the methods and their respective confidence intervals. With the random numbers that are stored in the states folder for each scenario every single table can be recreated without having to run everything.

Note: For 10000 repetitions it took my computer about 6h to run all scenarios.

After the generation of the tables and the estimations, the analysis follows in the `EstimationAnalysis.R` file. This requires the `estimates.rds` file which needs to be loaded there. Here the data is prepared for the analysis, first insights are gained through an exporatory data analysis and the estimations of tables with zeros are handled. One time with the complete case analysis approach and one time with the available case analysis approach. Then the performance measurements are computed and visualized, for each scenario, method and handling approach individually.

In the plots folder all figures used in the thesis can be found that are created in the analysis step.
