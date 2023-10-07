Author: Malyon D Bimler

Contains all R functions called up by the data, model, analyses and figures scripts.

##analyses

**calc_com_ntwk_metrics.R** - takes the output of the function below and returns community-level metrics

**calc_sp_ntwk_metrics.R** - takes a scaled alpha matrix (as samples) and returns species-level metrics - interaction strength, faciliation vs competition, interaction loops etc.  

##figures

**post_pred_checks.R** - runs posterior predictive checks on model outputs

##model

**dataprep.R** - preps a fecundity dataframe into the listed data format accepted by STAN

**scale_interactions.R** - scales the interaction parameters with demographic growth rates wo that they are in a per-capita unit

**stan_modelcheck_rem.R** - some functions for model validation and graphs
