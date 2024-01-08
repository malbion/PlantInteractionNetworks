Author: Malyon D Bimler

Includes all the data and code to run the models, analyses and figures for the paper and appendices, in a self-contained way. This folder should only contain 'necessary' code and figures.

Environmental categories are counted from '1' to '3' and are identified with the 'comm' variable in the code:  
1 - 0 to 8% canopy cover  
2 - 8 to 18% canopy cover  
3 - 18 to 40% canopy cover  


## Structure: 

**clean_data:** folder containing the datasets, also downloadable from Dryad (doi:10.5061/dryad.3tx95x6nf)  
**functions:** folder to store all functions called up by any of the scripts in parent dir   
**uploaded_results:** folder containing intermediate data files and results to run analyses and figures without having to run the model

**analyses.R** - runs all results analyses  
**figures.R** - codes figures for manuscript  
**joint_model.stan** - joint model framework code called by model.R  
**model.R** - applies joint model framework, transforms interactions and saves output   
**run_me.sh** - bash script to run model.R for each environmental category, the analyses and the figures scripts  

Scripts should be run in the following order:  
1 - model.R (for each environmental category)  
2 - analyses.R  
3 - figures.R  

Note that analyses.R and figures.R depend on data files that are written and saved by the previous scripts, into the folders described below. Copies of those data files have been saved in the uploaded_results folder so that analyses and figures can still be conducted without having to run the model. Note that given every model run is unique, those data files will differ slightly from any data files output by running the model on your personal computer. 

The scripts above call to the following folders, not uploaded here, to store results. For best results, create these same folders. They should be empty previous to running the models and analyses:  

**model:** store files related to the joint model - output, validation tests, transformed parameters. Contains the following folders:  
./output/1 /2 /3  
raw output from the joint model including full model fit (model_fit.Rdata), posterior draws and raw parameter samples. Samples are taken from 80% of the posterior.  
./transformed/1 /2 /3  
tranformed parameters: lambda, unscaled interaction parameters, growth rates and scaled interactions.  
./validation/1 /2 /3   
plots for joint model validation. Includes traceplots, posterior uncertainty intervals, posterior predictive check, interactions from the IFM vs REM and rstan diagnostic plots.  

**analyses:** store analysis of parameter outputs and results   

**figures:** store figures included in the main text and appendices of the paper  





