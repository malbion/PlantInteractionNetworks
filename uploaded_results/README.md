# README - Uploaded results

This folder contains intermediate data and results files which are created by model.R and analyses.R and used for subsequent scripts. Copies of these files are stored here so that analyses and figures can be done without having to run the models. Running the models and analyses will store local versions of these files in separate folders so that the ones stored here are not at risk of being overwritten.   

### Structure   

**1, 2, 3:** folder containing the scaled interaction matrices for each environmental category.   
output from: model.R    
input into: analyses.R, figures.R  

**coms_metrics_mat.csv:** results file storing metrics calculated at the level of each environmental category   
output from: analyses.R   
input into: figures.R  

**species_metrics_mat.csv:** results file storing metrics calculated for each focal species in each environmental category 
output from: analyses.R   
input into: figures.R  






