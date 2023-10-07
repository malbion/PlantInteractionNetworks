# Malyon Bimler
# CompNet 2019 

# Community-level network analysis

# Given a dataframe of species network metrics, spit out the community-level metrics
# Takes the output of calc.sp.ntwk.metrics()

library(igraph)

calc.com.ntwk.metrics <- function(metrics,   # the df of metrics
                                  imat,   # the community matrices for RI 
                                  ...) {

  # total number of neighbours 
  Nb <- dim(metrics)[1]
  # total number of focals
  Sp <- Nb - sum(is.na(metrics[ , 1, 1]))

  # COMMUNITY MEDIANS
  #-----------------------------------------------------------------------------------------
  temp <- t(apply(metrics, c(2, 3), median, na.rm = T))
  temp <- as.data.frame(temp)
  
  # STANDARDISED INTERACTION STRENGTH
  #-----------------------------------------------------------------------------------------
  # net magnitude of interactions in each community = median strength for all species standardised by sp richness 
  temp[ , 'rC_sum_aij'] <- temp[ , 'C_sum_aij']/Nb
  temp[ , 'rF_sum_aij'] <- temp[ , 'F_sum_aij']/Nb
  
  # LINKAGE DENSITY 
  #-----------------------------------------------------------------------------------------
  temp[ , 'LD'] <- apply(metrics, 3, function(metsamp) {
      # LD 
      (sum(sapply(1:Sp, function(i){ 
        metsamp[i, 'sum_|aji|'] / sum(metsamp[ , 'sum_|aji|']) * 2^metsamp[i, 'H.out']
      })) + sum(sapply(1:Sp, function(i){ 
        metsamp[i, 'sum_|aij|'] / sum(metsamp[ , 'sum_|aji|']) * 2^metsamp[i, 'H.in']
      }))
      )* 0.5
  })

  # WEIGHTED CONNECTANCE 
  #-----------------------------------------------------------------------------------------
  temp[ , 'wc'] <- temp[ , 'LD']/Nb
  
  # RELATIVE INTRANSITIVITY
  #-----------------------------------------------------------------------------------------
  comp <- na.omit(metrics[ , 'comp.rank' , ]*Sp) # get the number of species each focal outcompetes
  temp[ , 'RI'] <- apply(comp, 2, function(cvec){
    # RI index = [var(obs) - var(min) / var(max) - var(min) ] 
    # 1 = very hierarchical (transitive)
    var(cvec) / var(seq(0, Sp-1))
  })
  
  # INTERACTION MEDIAN AND SD 
  #-----------------------------------------------------------------------------------------
  temp[ , 'mean_aij'] <- apply(imat, 3, function(x){diag(x) <- 0; mean(x)})
  temp[ , 'sd_aij'] <- apply(imat, 3, function(x){diag(x) <- 0; sd(x)})
  
  # MODULARITY
  #-----------------------------------------------------------------------------------------
  temp[ , 'mod'] <- apply(imat, 3, function(x){
    G <- graph.adjacency(-x, weighted = T, mode = 'directed')
    m <- igraph::cluster_spinglass(G, implementation = 'neg')
    return(m$modularity)
  })
  
  # # LIMITS
  # #-----------------------------------------------------------------------------------------
  # temp[ , 'per.capita.limit'] <- temp[ , 'aii']/temp[ , 'mean_aij']
  # temp[ , 'EF.limit'] <- temp[ , 'eff.comp.intra']/temp[ , 'eff.comp']

  return(temp)
  
}


