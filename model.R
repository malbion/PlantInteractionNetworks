# Author: Malyon D. Bimler
# Objective: running a joint model combining to estimate all interactions at the same time. 
# 
# Run in Terminal: 
#   
#   R CMD BATCH '--args 0' model.R model/0_model_script.Rout
# 
# substituting '0' for the comm number (0, 1, 2, or 3) 
# 
# This calls up model.R, which runs joint_model.stan. 
# results got to model/output, /transformed and /validation

#-------------------------------------------------------------------------------------------
# PRELUDE 
#-------------------------------------------------------------------------------------------

Sys.time()
# Get arguments from bash script
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# take comm name as an argument from bash
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
if (length(args)>1) {
  stop("Model can only be run on 1 comm at a time.n", call.=FALSE)
}
comm <- args[1]

# set up R environment
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

setwd('~/Dropbox/Work/Projects/2018_Compnet/stormland/')

source('functions/model/data_prep.R')
source('functions/model/stan_modelcheck_rem.R')
source('functions/model/scale_interactions.R')

#-------------------------------------------------------------------------------------------
# PREPARE DATA
#-------------------------------------------------------------------------------------------

fecundities <- read.csv(paste0('clean_data/fecundities', comm, '.csv'), stringsAsFactors = F)

stan.data <- data_prep(perform = 'seeds', focal = 'focal', 
                       nonNcols = 4, df = fecundities)

# I need to keep note of how focals and neighbours are indexed
key_speciesID <- unlist(read.csv(paste0('clean_data/key_speciesID', comm, '.csv'), stringsAsFactors = F))
key_neighbourID <- unlist(read.csv(paste0('clean_data/key_neighbourID', comm, '.csv'), stringsAsFactors = F))

# ensure neighbours are linearly independent across the whole dataset
N_all <- apply(fecundities[ , 5:dim(fecundities)[2]], c(1,2), as.numeric)
X_all <- cbind(model.matrix(~as.factor(fecundities$focal)), N_all)
R_all <- pracma::rref(X_all)
Z_all <- t(R_all) %*% R_all
indep <- sapply(seq(1, dim(Z_all)[1], 1), function(k){ 
  ifelse(Z_all[k, k] == 1 & sum(Z_all[k, -k]) == 0, 1, 0)
}) #
all(indep == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep == 1)) message('WARNING neighbours are not linearly independent') 

message(paste0('Community selected: ', comm))
message(paste0('Fecundity data dimensions = ', dim(fecundities)[1], ', ', dim(fecundities)[2]))
message(paste0('Number of focals = ', length(key_speciesID)))
message(paste0('Number of neighbours = ', length(key_neighbourID)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data$Q)/(stan.data$S*stan.data$T)))

#-------------------------------------------------------------------------------------------
# RUN THE MODEL
#-------------------------------------------------------------------------------------------

# run model chains separately then combine then to avoid memory issues
stan.seed <- 52
fit <- list()
for (i in 1:4) {
  fit[i] <- stan(file = 'joint_model.stan',
                 data =  stan.data,               # named list of data
                 chains = 1,
                 warmup = 5000,          # number of warmup iterations per chain
                 iter = 10000,            # total number of iterations per chain
                 thin = 20,             # this ensures we get 1000 samples at the end
                 refresh = 100,         # show progress every 'refresh' iterations
                 control = list(max_treedepth = 20,
                                adapt_delta = 0.99),
                 seed = stan.seed,
                 chain_id = i
  )
}
fit <- sflist2stanfit(fit)  # combine chains

# select parameters of interest
param.vec <- fit@model_pars[!fit@model_pars %in% c('lp__')]

# Save the raw output:
save(fit, file = paste0('model/output/', comm, '/model_fit.Rdata')) # model fit
# draws from the posterior
joint.post.draws <- extract(fit)
save(joint.post.draws, file = paste0('model/output/', comm, '/post_draws.Rdata'))
# Save mean, 10% and 90% quantiles for each parameter, as well as n_eff and Rhat
fit_sum <- summary(fit, pars = param.vec, probs = c(0.1, 0.9))$summary
write.csv(fit_sum, file = paste0('model/output/', comm, '/summary_of_draws.csv'), row.names = T)
# Save the logarithm of the (unnormalized) posterior density (lp__)
log_post <- unlist(extract(fit, 'lp__'))
write.csv(log_post, file = paste0('model/output/', comm, '/log_post.csv'), row.names = F)

#-------------------------------------------------------------------------------------------
# MODEL VALIDATION AND DIAGNOSTICS
#-------------------------------------------------------------------------------------------

# Diagnostics: https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
stan_diagnostic(fit, paste0('model/validation/', comm, '/'), 
                params = c('gamma_i', 'ndd_betaij', 'weight', 'disp_dev')) # we can ignore response and effect
# Traceplots and posterior uncertainty intervals
stan_model_check(fit, paste0('model/validation/', comm, '/'), params = param.vec)
# Posterior predictive check
mu <- extract(fit, pars = 'mu')[[1]]
write.csv(mu, paste0('model/output/', comm, '/mu_samples.csv'), row.names = F)
disp_dev <- extract(fit, pars = 'disp_dev')[[1]]
write.csv(disp_dev, paste0('model/output/', comm, '/disp_dev_samples.csv'), row.names = F)
mu2 <- extract(fit, pars = 'mu2')[[1]]
write.csv(mu2, paste0('model/output/', comm, '/mu2_samples.csv'), row.names = F)

stan_post_pred_check(mu, disp_dev, 'mu', paste0('model/validation/', comm, '/'), stan.data)
stan_post_pred_check(mu2, disp_dev, 'mu2', paste0('model/validation/', comm, '/'), stan.data)

#-------------------------------------------------------------------------------------------
# SAVE AND VERIFY INTERACTION PARAMETER ESTIMATES
#-------------------------------------------------------------------------------------------

# JOINT interactions 
betaij <- joint.post.draws$ndd_betaij
betaijS <- as.data.frame(aperm(betaij, perm = c(1, 3, 2)))
colnames(betaijS) <- grep('ndd_betaij', rownames(fit_sum), value = T)
write.csv(betaijS, paste0('model/output/', comm, '/joint_betaij_samples.csv'), row.names = F)

# RIM interactions only
rim_betaij <- joint.post.draws$ri_betaij
rim_betaijS <- as.data.frame(aperm(rim_betaij, perm = c(1, 3, 2)))
colnames(rim_betaijS) <- grep('ri_betaij', rownames(fit_sum), value = T)
write.csv(rim_betaijS, paste0('model/output/', comm, '/RIM_betaij_samples.csv'), row.names = F)

# NDDM interactions only
Q <- as.vector(t(stan.data$Q)) # the t() is very important to fill Q by row! 
nddm_betaij <- sweep(betaijS, 2, Q, `*`) # replace uninferrable interactions with 0 to get the NDDM estimates only 
write.csv(nddm_betaij, paste0('model/output/', comm, '/NDDM_betaij_samples.csv'), row.names = F)

# VERIFY WITH SOME PLOTS
# interactions inferred by the NDDM only
nddm_betaij_inf <-as.matrix(nddm_betaij[  , which(colSums(nddm_betaij) != 0)]) # NDDM without non-inferrables
# RIM estimates for those inferrable interactions
rim_betaij_inf <- as.matrix(rim_betaijS[  , which(colSums(nddm_betaij) != 0)]) # NDDM without non-inferrables
# RIM estimates for non-inferrable interactions only
rim_betaij_noinf <- as.matrix(rim_betaijS[  , which(colSums(nddm_betaij) == 0)])

# Check estimates of inferrable interactions from both models 
png(paste0('model/validation/', comm, '/nddm_vs_rim_alphas.png'))
plot(nddm_betaij_inf, rim_betaij_inf, 
     xlab = 'NDDM interactions (inferrable only)', 
     ylab = 'RIM interactions (inferrable only)',
     xlim = c(min(nddm_betaij), max(nddm_betaij)),
     ylim = c(min(nddm_betaij), max(nddm_betaij)))
abline(0,1)
dev.off()

# check distribution of inferrable and non-inferrable interactions
png(paste0('model/validation/', comm, '/betaij_est_distr.png'))
par(mfrow=c(3,1))
hist(nddm_betaij_inf, xlab = "", breaks = 30,
     main = "Inferrable interactions (NDDM)", xlim = c(min(betaijS), max(betaijS)))
hist(rim_betaij_inf,  xlab = "", breaks = 30,
     main = 'Inferrable interactions (RIM)', xlim = c(min(betaijS), max(betaijS)))
hist(rim_betaij_noinf,  xlab = "", breaks = 30,
     main = 'Non-inferrable interactions (RIM)', xlim = c(min(betaijS), max(betaijS)))
dev.off()

#-------------------------------------------------------------------------------------------
# SCALE INTERACTIONS ACCORDING TO ANNUAL PLANT POP MODEL
#-------------------------------------------------------------------------------------------

# exponentiate gamma_i to get lambda_i
lambda_samples <- exp(joint.post.draws$gamma_i)
write.csv(lambda_samples, paste0('model/transformed/', comm, '/lambda_samples.csv'), row.names = F)

# Get seed rates (sruvival and germination)
seeds <- read.csv('clean_data/seed_rates.csv', row.names = 1) # average seed rates 
# scale interactions
scaled <- scale_interactions(betaij, lambda_samples, key_speciesID, seeds)
# save
write.csv(scaled$r_i, paste0('model/transformed/', comm, 'ri_demogrwth.csv'), row.names = T)
scaled.betas <- scaled$scaled.betas
save(scaled.betas, file = paste0('model/transformed/', comm, '/scaled_betaij_matrices.Rdata')) 
# the scaled betas are now the 80%CI 

#-------------------------------------------------------------------------------------------
# END 
gc(verbose = T)
Sys.time()

