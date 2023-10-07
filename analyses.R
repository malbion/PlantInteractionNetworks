# Author: Malyon D. Bimler
# Objective: Time to analyse this! 

#-------------------------------------------------------------------------------------------
# PRELUDE
#-------------------------------------------------------------------------------------------

library(plyr)
library(dplyr)
library(magrittr)
library(qgraph)

source('functions/analyses/calc_sp_ntwk_metrics.R')
source('functions/analyses/calc_com_ntwk_metrics.R')

# Select communities
communities <- list.dirs('model/transformed/', full.names = F, recursive = F)
# communities <- c('1', '2', '3') 

#-------------------------------------------------------------------------------------------
# EXAMINE SCALED INTERACTION MATRICES - MEDIAN, % F vs C, OVERLAP w/ 0 
#-------------------------------------------------------------------------------------------
inter.mat <- list()
inter.mat <- sapply(communities, function(c){
  load(paste0('model/transformed/', c, '/scaled_betaij_matrices.Rdata'))
  inter.mat[[as.numeric(c)]] <- scaled.betas
})

message("Focal species in Open: ")
attributes(inter.mat$`1`)[[2]]$species
message("Focal species in Intermediate: ")
attributes(inter.mat$`2`)[[2]]$species
message("Focal species in Shady: ")
attributes(inter.mat$`3`)[[2]]$species


### Intraspecific Interactions
#-----------------------------
bii <- lapply(inter.mat, apply, 3, diag)
lapply(bii, range)
message('Lower (10%) and upper (90%) quantiles of intraspecific interactions')
lapply(bii, quantile, c(0.1, 0.9))
lapply(bii, apply, 2, median) -> bii_meds  # median
# median proportion of competitive interactions:
biiC <- lapply(bii, apply, 2, function(x) length(x[x>0])/length(x))
lapply(biiC, median)
# lower and upper quantiles:
lapply(biiC, quantile, c(0.1, 0.9))
# median proportion of facilitative interactions:
biiF <- lapply(bii, apply, 2, function(x) 1 - length(x[x>0])/length(x))
lapply(biiF, median)       
# lower and upper quantiles: 
lapply(biiF, quantile, c(0.1, 0.9))

lapply(bii, function(y) {
  apply(y, 1, function(x){
    CIx <- x[x > quantile(x, 0.1) & x < quantile(x, 0.9)]
    abs(sum(sign(CIx))) / length(CIx) * 100
  })
}) -> bii_overlap
bii_overlap
lapply(bii_overlap, function(z) {
  length(which(z == 100)) / length(z) * 100
})

### Interspecific Interactions
#-----------------------------
bij <- lapply(inter.mat, apply, 3, function(mat) {
  diag(mat) <- NA    # remove intraspecifics
  mat <- mat[!is.na(mat)==T]
  return(mat)
})
lapply(bij, range, na.rm = T) 
lapply(bij, quantile, c(0.1, 0.9), na.rm = T)   # quantiles
lapply(bij, apply, 2, median) -> bij_meds   # median

# median proportion of competitive interactions:
bijC <- lapply(bij, apply, 2, function(x) length(x[x>0])/length(x))
lapply(bijC, median)
# quantiles: 
lapply(bijC, quantile, c(0.1, 0.9))
# proportion of interactions that don't overlap 0:
overlaps <- lapply(bij, function(tens) {
  apply(tens, 1, function(vec){
    CIv <- vec[vec > quantile(vec, 0.1) & vec < quantile(vec, 0.9)]
    ifelse(sign(min(CIv)) == sign(max(CIv)), 1, 0) # return 1 if all samples are of the same sign
  })
})
lapply(overlaps, function(mat) { sum(mat) / length(mat) * 100})
# ratio of intra vs inter median:
ii.vs.ij <- mapply('/', bii_meds, bij_meds)  
apply(ii.vs.ij, 2, median)
apply(ii.vs.ij, 2, quantile, c(0.1, 0.9))

#-------------------------------------------------------------------------------------------
# SPECIES-LEVEL METRICS
#-------------------------------------------------------------------------------------------
# get plot-level abundances (not neighbourhoods!)
sp.abunds <- read.csv('clean_data/plot_species_abundances.csv', stringsAsFactors = F)
sp.abunds <- split(sp.abunds, as.factor(sp.abunds$Class.Canopy))
# get species-level metrics
sp.ntwk.metrics <- mapply(calc.sp.ntwk.metrics, communities, inter.mat, sp.abunds)
# flatten all into one big dataframe
smm <- lapply(sp.ntwk.metrics, adply, c(1,3))
smm <- do.call(rbind, smm)
smm$spcID <- paste(smm$species, smm$com, sep = '_')
# save
write.csv(smm, 'analyses/species_metrics_mat.csv', row.names = F)
# smm <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)

# there is no correlation between alpha_ii and abundance 
focals <- smm[is.na(smm$aii) == F, ]
cor(focals$aii, focals$rel.abund)
cor(focals$aii, focals$com.abund)

#-------------------------------------------------------------------------------------------
# COMMUNITY LEVEL METRICS & NETWORK PROPERTIES
#-------------------------------------------------------------------------------------------
# Get community-level network metrics
comm.ntwk.metrics <- mapply(calc.com.ntwk.metrics, sp.ntwk.metrics, inter.mat,
                            SIMPLIFY = F)
# flatten up and save
cmm <- do.call(rbind, comm.ntwk.metrics)
cmm$comm <- as.vector(sapply(as.numeric(communities), rep, nrow(cmm)/length(communities)))
write.csv(cmm, 'analyses/coms_metrics_mat.csv', row.names = F)

# Check pairs() and correlation of vars
cmm.vars <- cmm[ , c('RI', 'rF_sum_aij', 'aii', 'rC_sum_aij', 'wc', 'sd_aij', 'mod', 'comm')]
cmm.vars$rF_sum_aij <- abs(cmm.vars$rF_sum_aij)
png('analyses/comm_pairs.png', width = 1000, height = 1000)
pairs(cmm.vars[ ,1:7])
dev.off()
cor(cmm.vars[ , 1:7])
write.csv(cor(cmm.vars[ , 1:7]), 'analyses/comm_cor.csv')

cmm %>% group_by(com) %>% 
  summarise('medRI' = median(RI),
            'minRI' = quantile(RI, 0.1),
            'maxRI' = quantile(RI, 0.9),
            'medwc' = median(wc),
            'minwc' = quantile(wc, 0.1),
            'maxwc' = quantile(wc, 0.9),
            'medMod' = median(mod),
            'minMod' = quantile(mod, 0.1),
            'maxMod' = quantile(mod, 0.9),
            'medrF_sum_aij' = median(rF_sum_aij),
            'minrF_sum_aij' = quantile(rF_sum_aij, 0.1),
            'maxrF_sum_aij' = quantile(rF_sum_aij, 0.9),
            'medrC_sum_aij' = median(rC_sum_aij),
            'minrC_sum_aij' = quantile(rC_sum_aij, 0.1),
            'maxrC_sum_aij' = quantile(rC_sum_aij, 0.9),
            'medaii' = median(aii),
            'minaii' = quantile(aii, 0.1),
            'maxaii' = quantile(aii, 0.9),
            'medsd_aij' = median(sd_aij),
            'minsd_aij' = quantile(sd_aij, 0.1),
            'maxsd_aij' = quantile(sd_aij, 0.9)) %>% as.data.frame() -> cmm.sum
write.csv(cmm.sum, 'analyses/coms_metrics_summary.csv', row.names = F)
cmm.sum

cmm %>% group_by(com) %>% 
   summarise('medASYM' = median(perc.exploiter + perc.exploited),
             'minASYM' = quantile(perc.exploiter + perc.exploited, 0.1),
             'maxASYM' = quantile(perc.exploiter + perc.exploited, 0.9),
             'medCOOP' = median(perc.coop),
             'minCOOP' = quantile(perc.coop, 0.1),
             'maxCOOP' = quantile(perc.coop, 0.9), 
             'medCOMP' = median(perc.fullcomp),
             'minCOMP' = quantile(perc.fullcomp, 0.1),
             'maxCOMP' = quantile(perc.fullcomp, 0.9)) %>% as.data.frame() -> cmm.loops
cmm.loops


#-------------------------------------------------------------------------------------------
# SPECIES ACROSS COMMUNITIES
#-------------------------------------------------------------------------------------------
smm <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)
# select focals
smm <- smm[is.na(smm$aii) == F, ]
# Find species that occur across all three coms as focals 
sprws <- table(smm$species)
species_3 <- names(sprws[sprws > 2500])
species_3 <- species_3[species_3 != 'others']
# Select only those 9 species 
smm <- smm[smm$species %in% species_3, ]
smmlist <- split(smm, as.factor(smm$species))
smmlist <- lapply(smmlist, function(ss) {split(ss, as.factor(ss$com))})
lapply(smmlist, function(ss) {
  do.call(cbind, lapply(ss, function(sc) {
    rbind(mean(sc$perc.exploited) + mean(sc$perc.exploiter),
          mean(sc$perc.fullcomp),
          mean(sc$perc.coop))
  })) -> foo
  row.names(foo) <- c('ant', 'comp', 'coop')
  foo <- foo*100
  d12 <- foo[ , 2] - foo[ , 1]
  d13 <- foo[ , 3] - foo[ , 1]
  d23 <- foo[ , 3] - foo[ , 2]
  cumd <- abs(d12) + abs(d23) + abs(d13)
  foo <- cbind(foo, d12, d23, d13, cumd)
  return(foo)
}) -> strats

lapply(strats, function(foo) {sum(foo[ , 'cumd'])}) # to check which species moves the most or the least

#-------------------------------------------------------------------------------------------
# BETAS FROM COM TO COM
#-------------------------------------------------------------------------------------------

# load scaled interactions
inter.mat <- list()
inter.mat <- sapply(c(1,2,3), function(c){
  load(paste0('model/transformed/', c, '/scaled_betaij_matrices.Rdata'))
  inter.mat[[as.numeric(c)]] <- scaled.betas
})
# species which appear in all 3 communities 
species_3 <- c("ARCA", "HYGL", "PEAI", "PLDE", "POCA", "PTGA", "STPA", "VERO", "WAAC")
# select only interactions which appear in all 3 coms 
inter.3 <- lapply(inter.mat, function(betas) {
  betas <- betas[species_3, species_3, ]
  # select the median and upper / lower quantile limits
  med.b <- apply(betas, c(1,2), median)  
  min.b <- apply(betas, c(1,2), quantile, 0.1)
  max.b <- apply(betas, c(1,2), quantile, 0.9)
  return(list(med.b, min.b, max.b))
})
names(inter.3) <- ntw.name

delta12 <- inter.3$Open[[1]] - inter.3$Intermediate[[1]]
delta13 <- inter.3$Open[[1]] - inter.3$Shady[[1]]

signswitch <- vector()
for (i in species_3) {
  for (j in species_3) {
    signswitch <- c(signswitch, sum(sign(inter.3$Open[[1]][i,j]), 
        sign(inter.3$Intermediate[[1]][i,j]), 
        sign(inter.3$Shady[[1]][i,j])))
}}
table(abs(signswitch))
