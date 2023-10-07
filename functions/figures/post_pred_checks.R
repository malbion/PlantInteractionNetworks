# Function for posterior predictive checks

post_pred_fig <- function(comm,    # community
                          mu_name,      # mu or mu2 ? 
                          ...) {
  
  # get real data
  fecundities <- read.csv(paste0('clean_data/fecundities', comm, '.csv'), stringsAsFactors = F)
  seeds <- fecundities$seeds
  focalobs <- as.numeric(as.factor(fecundities$focal))
  
  # extract mu and phi
  mu <- read.csv(paste0('model/output/', comm, '/', mu_name, '_samples.csv'))
  phi <- read.csv(paste0('model/output/', comm, '/disp_dev_samples.csv'))
  phi <- (phi^2)^(-1)
  
  # generating posterior predictions
  seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      seed_pred[i, j] <- rnbinom(1, mu = mu[i, j], size = phi[i, focalobs[j]])  
    }
  }
  # log transform seed predictions
  seed_pred <- log(seed_pred)
  # and using just the median of the parameters
  m.mu <- apply(mu, 2, median)
  m.phi <- apply(phi, 2, median)
  mean_seed_pred <- sapply(1:length(m.mu) , function(x) {
    rnbinom(1, mu = m.mu[x], size = m.phi[focalobs[x]])
  })
  mean_seed_pred <- log(mean_seed_pred)
  
  # get maximum density for plot limits
  max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                       max(density(seeds)$y)))
  
  
  png(paste0('figures/ppcheck_', comm, '_', mu_name, '.png'), width=500, height=500)
  # start a plot with the first draw 
  ppc.plot <- plot(density(seed_pred[1, ]), 
                   type = 'n',
                   bty = 'n',
                   ylim = c(0, max.density+0.03),  
                   xlim = c(-1, 10),
                   col = 'lightgrey',
                   ylab = 'Seed probability density',
                   xlab = 'Log seed production',
                   main = '',
                   # sub = '(grey = predicted, red = observed)'
  ) 
  for (i in 1:dim(seed_pred)[1]) {
    # add a line for each draw
    ppc.plot <- lines(density(seed_pred[i, ]), col = 'lightgrey')
  }
  # add the 'mean prediction'
  ppc.plot <- lines(density(mean_seed_pred), col = 'black', lwd = 1)  
  # add the actual data
  ppc.plot <- lines(density(log(seeds)), col = 'red', lwd = 1)  
  # add legend
  ppc.plot <- legend('topright', 
                     legend = c('predicted', 'observed'), 
                     col = c('lightgrey', 'red'), lwd = c(6, 1),
                     bty = 'n')
  print(ppc.plot)
  dev.off()
}
