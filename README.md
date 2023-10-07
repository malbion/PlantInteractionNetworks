Author: Malyon D Bimler

Includes all the data and code to run the models, analyses and figures for the paper and appendices, in a self-contained way. This folder should only contain 'necessary' code and figures.

Environmental categories are counted from '1' to '3' :
1 - 0 to 8% canopy cover
2 - 8 to 18% canopy cover
3 - 18 to 40% canopy cover


##Structure: 

**functions:** folder to store all functions called up by any of the scripts in parent dir. 

**analyses.R** - runs all results analyses 
**figures.R** - codes figures for manuscript
**joint_model.stan** - joint model framework code called by model.R
**model.R** - applies joint model framework, transforms interactions and saves output 
**run_me.sh** - bash script to run model.R for each environmental category, the analyses and the figures scripts

Scripts should be run in the following order: 
1 - model.R (for each environmental category)
2 - analyses.R 
3 - figures.R



The scripts provided also call to the following folders, not uploaded here, to retrieve data and store results - for best results, create these same folders 

**clean_data:** folder which should contain the datasets downloadable from Dryad (doi:10.5061/dryad.3tx95x6nf) 

The folders below should be empty previous to running the models and analyses:

**model:** folder to store files related to the joint model - output, validation tests, transformed parameters. Contains the following folders:
./output/1 /2 /3
raw output from the joint model including full model fit (model_fit.Rdata), posterior draws and raw parameter samples. Samples are taken from 80% of the posterior.
./transformed/1 /2 /3
tranformed parameters: lambda, unscaled interaction parameters, growth rates and scaled interactions.
./validation/1 /2 /3 
plots for joint model validation. Includes traceplots, posterior uncertainty intervals, posterior predictive check, interactions from the IFM vs REM and rstan diagnostic plots.

**analyses:** folder to store analysis of parameter outputs and results 

**figures:** folder to store figures included in the main text and appendices of the paper




