# Scale the interaction parameters according to the annual plant pop model
# Returns a 3D array of dimensions [N_species, N_neighbs, N_samples]

#### This is a slightly modified version of the Natmeth code #### 

# NB: demographic scaling can only be done on those interactions where we have germination
# rates for both i and j aka interactions between focals only.

scale_interactions <- function(betas,     # joint interactions
                               lambdas,   # exp(gamma)i)
                               key_speciesID,
                               seeds,     # average seed rates
                               ...){ 
  
  N_samples <- dim(betas)[1]  # number of samples
  betamat <- aperm(betas, c(2, 3, 1))  # order as focals / neighbours / samples
 
  # SCALING: ----------------------------------------------
  # scaled beta = (g_j * model_param_ij) / ln(eta_i)      |
  # --- see Eq 15 Supps NATMETH ---------------------------
  
  # CALCULATE THE DIVISOR ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # r_i = ln(eta_i) = ln( (g_i * lambda_i ) / (1 - (1 - g_i) s_i) )
  # reminder: lambda_i = exp(gamma_i) 
    r_i <- sapply(c(1:length(key_speciesID)), function(x) {
    prod <- seeds[key_speciesID[x], 'germ'] / (1 - (1 - seeds[key_speciesID[x], 'germ'])*seeds[key_speciesID[x], 'surv'])
    if (prod < 0.1) {prod <- 0.1}
      log(lambdas[,x ] * prod)
    })
    
    dimnames(r_i) <- list('samples' = seq(1:N_samples), 'species' = key_speciesID)
    ri_means <- colMeans(r_i)[order(colMeans(r_i), decreasing = T)]
    message('Growth rates:')
    print(ri_means)
    
    # REMOVING GROWTH RATES < 1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(length(ri_means[ri_means < 1]) > 0){
      # If any growth rates are below 1 .... 
      to.rplc <- names(ri_means[ri_means < 1])
      # replace with the next lowest growth rate (that's over 1)
      rplc.ri <- min(ri_means[ri_means > 1])
      rplc.sp <- names(ri_means[ri_means == rplc.ri])
     
      message(paste0('Growth rates for ', to.rplc, ' replaced by ', rplc.sp, '\n'))
      for(i in to.rplc) {
        r_i[ , i] <- r_i[, rplc.sp]
      }
 }
    # SCALE THE INTERACTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # scaled beta = (g_j * model_param_ij) / r_i 
    ### 1. g_j * model_param_ij: 
    beta.gj <- sapply(c(1:length(key_speciesID)), function(j) {
      # dimensions: betamat[ i, j, samples]
      seeds[j, 'germ'] * betamat[ , j, ]
    }, simplify = 'array', USE.NAMES = T)
    # dimensions: beta.gj[i, samples, j]
    ### 2. beta.gj / r_i
    bij.scaled <- sapply(c(1:length(key_speciesID)), function(i) {
      beta.gj[i, , ] / r_i[ , i]
    }, simplify = 'array', USE.NAMES = T)
    # dimensions: bij.scaled[samples, j, i]
    
    # reorder into [i, j, samples]
    final.betas <- aperm(bij.scaled, c(3, 2, 1))
    dimnames(final.betas) <- list('species' = key_speciesID,
                                  'neighbour' = key_speciesID,
                                  'samples' = seq(1, 1000))
    
    attributes(dimnames(final.betas)[[1]]) <- NULL    
    attributes(dimnames(final.betas)[[2]]) <- NULL

    # END FUNCTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    return(list('scaled.betas' = final.betas,
                'r_i' = r_i))
  
}


