#!/bin/bash


R CMD BATCH '--args 1' model.R model/1_model_script.Rout
R CMD BATCH '--args 2' model.R model/2_model_script.Rout
R CMD BATCH '--args 3' model.R model/3_model_script.Rout

R CMD BATCH analyses.R analyses.Rout
R CMD BATCH figures.R figures.Rout

