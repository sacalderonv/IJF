# Replication repository for "Bayesian estimation of a multivariate TAR model when the noise process distribution belongs to the class of Gaussian variance mixtures"

** Vanegas, LH. & Calderón, SA & Rondón, LM.***
Mail: <name_surname@unal.edu.co>


Preprint available at
[Arxiv](https://www.arxiv.org/pdf/2503.04593). working paper conditionally accepted at [International Journal of Forecasting](httts://forecasters.org/ijf)



## Overview

This is a repository for  replications of some results for paper "Bayesian estimation of a multivariate TAR model when the noise process distribution belongs to the class of Gaussian variance mixtures"

The main contents of the repository are the following:

- `R Files for 2 reg estimation, missespecification and identify Distribution/`: Folder containing 3 R scripts and 1 rds file where were stored 1000 replications in order to get whole results of table 10 and student-t columns of tables 1,2 and 5 for M1 structure model.
- `R files for 3 reg/`: Folder containing 3 R scripts and 6 rds files where were stored 1000 replications for ecah distribution error in order to get whole results of table 3 and 6 for M2 structure model.
- `data/`: Folder containing the data and script file used in the empirical application.


## Usage and computational specifications 
All file paths are relative to the root of the replication repository.
Please set your working directory accordingly.

All the estimation and simulations and analysis is done in R. The main package used for
estimation and forecasting is `mtarm` (0.1.5) with major
dependency `GIGrvg` (),`Formula` () where the values in parenthesis indicate
the package versions we used. However, the simulations for M1 structure in the script called "IJFSimulChequeoDistribution2regbaseFinalParalelizar.R" requires `paralell`(), `forearch` () and `doParallel` () in other to parallelize the simulations to check the missespecification of the distribution error. In addition, the following add-on packages
are loaded and attached in the replication scripts (in alphabetical
order): `expm` (), `ltsa` (), `Rfast` () and `tsDyn` . The results are mainly printed in the console of R. Only results of the simulation are recorded in files *.rds.



### For the M2 structure, that is a 2-variate 3 regime MTAR model

#### Script for Resume simulations
The script named "Resumen_Replicas_3Reg_Final.R" is used to obtain for instance, the results for tables 3, 4 and 6. This script is divided in three parts. Only changes must be made in forst part. Second and three parts are only operatives and mainly compute and print the results in the R console.
These results are printed in the R console. The results are given for columns(i.e for each distribution of the errors). Therefore you have run this script  for each distribution only modifying the lines indicated in the script. For instance, if you want the results for the slash distribution you must:

- First, uncomment the line to load the results of the replications for this distribution, in this case, the line 25 load("replicas_slash_1000_3reg1.rds"). Next, comment the line 23(i.e, the line that there is without comment by defect.).
- Change the names in line 32 repl_estimation<-repl_gaussian_estimation. The list repl_slash_estimation contains the results for the estimation.
- Change the names in line 35 repl<-repl_gaussian. The list repl_hyperbolic contains the results for the forecasting.
- Change to TRUE the object para.extra in line 45, that is 45 para.extra=TRUE. This is because Slash distribution has extra parameter.
- Run all code in the script and wait for the results in the console.

#### <p style="color:red; font-size:18px;">Script for running simulations font</p> 

The script named "SimulayEstimaMtar_Replicas3Reg.R" is used to simmulating, get the estimation of the paramaters, obtain the h.ahead forecasting and store the results for the number of replications set. This script depends on function in "script summarymtar_simulation.R", therefore this must be ran before running this script.

What is stored? For each replication of the MTAR model simulated, it is stored mainly the estimation of the parameters including the credible intervals and the prediction for each forecast horizon with its respective prediction interval.

What can you modify in the code in order to explore all distributions for errors? Initially, you can modify the distribution error(for instace, you can use: "Gaussian", "Student-t", "Hyperbolic", "Laplace", "Slash", "Contaminated normal") and the extra parameter (check lines 79 and 81). The number of replications is controlled by object "n_rep"(see line 15) which is set in 1000. Length of the time series is establishe in object named "long"" in line 33, and the steps ahead for the forecasting is set in the object "h.ahead" in line 34. You can also change the names of the objects where it will be stored the outputs. These names suggested are in lines from 17 to 22 for estimation, and from 25 to 30 for prediction. In this script, the example is for the slash distribution. If you change thsese names you must ensure thta this names be changed in lines 226(store the estimation), 227(store the forecasting) and 237.




