# Replication repository for "Bayesian estimation of a multivariate TAR model when the noise process distribution belongs to the class of Gaussian variance mixtures"

** Vanegas, LH. & Calderón, SA & Rondón, LM.***
Mail: <name_surname@unal.edu.co>


Preprint available at
[Arxiv](https://www.arxiv.org/pdf/2503.04593). working paper conditionally accepted at [International Journal of Forecasting](httts://forecasters.org/ijf)



## Overview

This is a repository for  replications of some results for paper "Bayesian estimation of a multivariate TAR model when the noise process distribution belongs to the class of Gaussian variance mixtures"

The main contents of the repository are the following:

-`R Files for 2 reg estimation, missespecification and identify Distribution/`: Folder containing 3 R scripts and 1 rds file where were stored 1000 replications in order to get whole results of table 10 and student-t columns of tables 1,2 and 5 for M1 structure model.
- `R files for 3 reg/`: Folder containing 3 R scripts and 6 rds files where were stored 1000 replications for ecah distribution error in order to get whole results of table 3 and 6 for M2 structure model.
- `data/`: Folder containing the data used in the empirical application.


## Usage and computational specifications 
All file paths are relative to the root of the replication repository.
Please set your working directory accordingly.

All the estimation and simulations and analysis is done in R. The main package used for
estimation and forecasting is `mtarm` (0.1.5) with major
dependency `GIGrvg` (),`Formula` () where the values in parenthesis indicate
the package versions we used. In addition, the following add-on packages
are loaded and attached in the replication scripts (in alphabetical
order): `Rfast` (), `expm` (), `ltsa` (). The results are mainly printed in the console of R. Only results of the simulation are recorded in files *.rds.



### For the M2 structure, that is a bivariate 3 regime MTAR model

