# Author: Malyon D. Bimler
# Objective: code the figures that will appear in the manuscript 

###########################################################################################
################################### MAIN TEXT FIGURES #####################################
###########################################################################################

# Setting up stuff that gets used elsewhere
ntw.name <- c('Open', 'Intermediate', 'Shady')
# species that appear in all 3 networks:
species_3 <- c("ARCA", "HYGL", "PEAI", "PLDE", "POCA", "PTGA", "STPA", "VERO", "WAAC")

library(magrittr)
library(paletteer)
library(tidyverse)
library(ggjoy)
library(cowplot)
require(Ternary)
require(dplyr)
require(plotfunctions)
library(mixtools)
library(qgraph)
library(fmsb)

# defining path to file
library(here)
i_am('figures.R')

# if you haven't run the model or the analyses, you can still create the figures
# by accessing the appropriate files stored in the 'uploaded_results' folder


#-------------------------------------------------------------------------------------------
# FIGURE 1: COMPETITIVE VS FACILITATIVE OUTPUT
#-------------------------------------------------------------------------------------------

smets <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)
# replace with line below if not running your own model or analyses
# smets <- read.csv('uploaded_results/species_metrics_mat.csv', stringsAsFactors = F)
smm <- smets[is.na(smets$aii) == F, ]
smm <- smm[ , c('species', 'com', 'spcID', 'sum_aij', 
                'C_sum_aji', 'F_sum_aji', 'C_sum_aij', 'F_sum_aij')]
smm$F_sum_aji <- - smm$F_sum_aji
smm$F_sum_aij <- - smm$F_sum_aij


smm %>% split(., as.factor(smm$spcID)) %>%
  lapply(., function(spc) {c(apply(spc[ ,4:8], 2, median), 
                             unique(spc$com), unique(spc$species))}) %>%
  do.call(rbind, .) %>% as.data.frame(., stringsAsFactors = F) -> spcmediansCF
spcmediansCF[ , 1:5] <- apply(spcmediansCF[ , 1:5], 2, as.numeric)
colnames(spcmediansCF)[6:7] <- c('com', 'species')

# This sets individual colours for each species
ell.col <- c(paletteer_d("ggthemes::Classic_20"), '#B31166')
names(ell.col) <- sort(unique(smm$species))

png('figures/effects.png', width=2100, height=1400)
par(mfrow=c(2, 3), cex = 2, cex.lab=1.4, cex.main = 1.6,
    mai = c(.4,.4,.4,.4), mar = c(4.2,4.2,2,0.3), oma=c(0,0,0,0))
# upper panels: output
for(com in 1:3) {
  xlow <- min(smm[smm$com == com, c('C_sum_aji', 'F_sum_aji')])
  xmax <- max(smm[smm$com == com , c('C_sum_aji', 'F_sum_aji')])
  plot(spcmediansCF[spcmediansCF$com == com, 'F_sum_aji'] ~ 
         spcmediansCF[spcmediansCF$com == com, 'C_sum_aji'], 
       xlab = 'Competitive output',
       ylab = 'Facilitative output',
       xlim = c(xlow, xmax), 
       ylim = c(xlow, xmax),
       type='n', bty='n')
  points(smm[smm$com == com , 'F_sum_aji'] ~ smm[smm$com == com, 'C_sum_aji'],  
         pch=20, cex = 1.5, col = adjustcolor('grey90', alpha.f=0.5))
  abline(0, 1, lty = 'dashed')
  # Confidence intervals
  for(i in unique(spcmediansCF[spcmediansCF$com == com, ]$species)) {
    upper.x <- quantile(smm[smm$com == com & smm$species == i, ]$C_sum_aji, 0.9)
    lower.x <- quantile(smm[smm$com == com & smm$species == i, ]$C_sum_aji, 0.1)
    lines(x = c(lower.x, upper.x),
          y = c(spcmediansCF[spcmediansCF$com == com & spcmediansCF$species == i, ]$F_sum_aji,
                spcmediansCF[spcmediansCF$com == com & spcmediansCF$species == i, ]$F_sum_aji),
          lwd = 4, col = ell.col[i])
    upper.y <- quantile(smm[smm$com == com & smm$species == i, ]$F_sum_aji, 0.9)
    lower.y <- quantile(smm[smm$com == com & smm$species == i, ]$F_sum_aji, 0.1)
    lines(x = c(spcmediansCF[spcmediansCF$com == com & spcmediansCF$species == i, ]$C_sum_aji,
                spcmediansCF[spcmediansCF$com == com & spcmediansCF$species == i, ]$C_sum_aji),
          y = c(lower.y, upper.y),
          lwd = 4, col = ell.col[i])
  }
  # Species medians
  points(spcmediansCF[spcmediansCF$com == com, 'F_sum_aji'] ~ 
           spcmediansCF[spcmediansCF$com == com, 'C_sum_aji'], 
         pch=23, cex = 2, bg = ell.col[which(names(ell.col) %in% smm[smm$com == com, ]$species)])
  title(ntw.name[com])
  # correlation
  legend('topleft', bty='n', cex=1.2,
         legend = paste0('r = ', round(cor(smm[smm$com == com , 'C_sum_aji'], 
                                           smm[smm$com == com , 'F_sum_aji']), 2)))
}   
# lower panels: input
sapply(c(1,2,3), function(x){
  cg <- smm[smm$com == x, ]
  # Background
  plot(cg[ , 'F_sum_aij'] ~ cg[ , 'C_sum_aij'], 
       xlab = 'Competitive Input',
       ylab = 'Facilitative Input',
       xlim = c(min(cg[ , 'C_sum_aij']), max(cg[ , 'C_sum_aij'])), 
       ylim = c(min(cg[ , 'F_sum_aij']), max(cg[ , 'F_sum_aij'])),
       type='n', bty='n')
  #  title(ntw.name[x])
  points(cg[ , 'F_sum_aij'] ~ cg[ , 'C_sum_aij'], pch=20, cex = 1.5,  
         col = adjustcolor('grey90', alpha.f=0.5))
  abline(0, 1, lty = 'dashed')
  # Foreground
  sm <- spcmediansCF[spcmediansCF$com == x, ]
  # Confidence intervals
  for(i in unique(sm$species)) {
    upper.x <- quantile(smm[smm$com == x & smm$species == i, ]$C_sum_aij, 0.9)
    lower.x <- quantile(smm[smm$com == x & smm$species == i, ]$C_sum_aij, 0.1)
    lines(x = c(lower.x, upper.x),
          y = c(sm[sm$species == i, ]$F_sum_aij,
                sm[sm$species == i, ]$F_sum_aij),
          lwd = 4, col = ell.col[i])
    upper.y <- quantile(smm[smm$com == x & smm$species == i, ]$F_sum_aij, 0.9)
    lower.y <- quantile(smm[smm$com == x & smm$species == i, ]$F_sum_aij, 0.1)
    lines(x = c(sm[sm$species == i, ]$C_sum_aij,
                sm[sm$species == i, ]$C_sum_aij),
          y = c(lower.y, upper.y),
          lwd = 4, col = ell.col[i])
  }
  # Species medians
  points(sm[ , 'F_sum_aij'] ~ sm[ , 'C_sum_aij'], 
         pch=23, cex = 2, bg = ell.col[which(names(ell.col) %in% sm$species)])
  # correlation
  legend('topleft', legend = paste0('r = ', round(cor(cg$C_sum_aij, cg$F_sum_aij), 2)), 
         bty='n', cex=1.2)
})
dev.off()

#-------------------------------------------------------------------------------------------
# FIGURE 2: COMPETITIVE VS FACILITATIVE INPUT + Xii VS OUT-STRENGTH
#-------------------------------------------------------------------------------------------

smets <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)
# replace with line below if not running your own model or analyses
# smets <- read.csv('uploaded_results/species_metrics_mat.csv', stringsAsFactors = F)
smm <- smets[is.na(smets$aii) == F, ]
smm <- smm[ , c('species', 'com', 'spcID', 'aii', 'sum_aji')]

smm %>% split(., as.factor(smm$spcID)) %>%
  lapply(., function(spc) {c(apply(spc[ ,4:5], 2, median),
                             unique(spc$com), unique(spc$species))}) %>%
  do.call(rbind, .) %>% as.data.frame(., stringsAsFactors = F) -> spcmedians
spcmedians[ , 1:2] <- apply(spcmedians[ , 1:2], 2, as.numeric)
colnames(spcmedians)[3:4] <- c('com', 'species')

png('figures/intra_vs_out.png', width=2100, height=700)
par(mfrow=c(1, 3), cex = 2, cex.lab=1.4, cex.main = 1.6,
    mai = c(.4,.4,.4,.4), mar = c(4.2,4.2,2,0.3), oma=c(0,0,0,0))
# lower panels
sapply(c(1,2,3), function(x){
  cg <- smm[smm$com == x, ]
  # Background
  plot(cg[ , 'aii'] ~ cg[ , 'sum_aji'], 
       xlab = 'Out-strength',
       ylab = 'Intraspecific interaction strength',
       xlim = c(min(cg[ , 'sum_aji']), max(cg[ , 'sum_aji'])), 
       ylim = c(min(cg[ , 'aii'] -0.02), max(cg[ , 'aii'] +0.02)),
       type='n', bty='n')
  title(ntw.name[x])
  abline(h=0, v=0, lty = 2)
  points(cg[ , 'aii'] ~ cg[ , 'sum_aji'], pch=20, cex = 1.5,  
         col = adjustcolor('grey90', alpha.f=0.5))
  # Foreground
  sm <- spcmedians[spcmedians$com == x, ]
  # Confidence intervals
  for(i in unique(sm$species)) {
    upper.x <- quantile(smm[smm$com == x & smm$species == i, ]$sum_aji, 0.9)
    lower.x <- quantile(smm[smm$com == x & smm$species == i, ]$sum_aji, 0.1)
    lines(x = c(lower.x, upper.x),
          y = c(sm[sm$species == i, ]$aii,
                sm[sm$species == i, ]$aii),
          lwd = 4, col = ell.col[i])
    upper.y <- quantile(smm[smm$com == x & smm$species == i, ]$aii, 0.9)
    lower.y <- quantile(smm[smm$com == x & smm$species == i, ]$aii, 0.1)
    lines(x = c(sm[sm$species == i, ]$sum_aji,
                sm[sm$species == i, ]$sum_aji),
          y = c(lower.y, upper.y),
          lwd = 4, col = ell.col[i])
  }
  # Species medians
  points(sm[ , 'aii'] ~ sm[ , 'sum_aji'], 
         pch=23, cex = 2, bg = ell.col[which(names(ell.col) %in% sm$species)])
  # correlation
  legend('topleft', legend = paste0('r = ', round(cor(cg$aii, cg$sum_aji), 2)), 
         bty='n', cex=1.2)
})
dev.off()

#-------------------------------------------------------------------------------------------
# FIGURE 3: BETA DISTRIBUTIONS FOR A SAMPLE OF SPECIES 
#-------------------------------------------------------------------------------------------

# load scaled interactions
inter.mat <- list()
inter.mat <- sapply(c(1,2,3), function(c){
  load(paste0('model/transformed/', c, '/scaled_betaij_matrices.Rdata'))
  # replace with line below if not running your own models
  # load(paste0('uploaded_results/', c, '/scaled_betaij_matrices.Rdata'))
  inter.mat[[as.numeric(c)]] <- scaled.betas
})
# species which appear in all 3 communities 
species_3 <- c("ARCA", "HYGL", "PEAI", "PLDE", "POCA", "PTGA", "STPA", "VERO", "WAAC")

# select only interactions between 3 species only
sub3 <- c("HYGL", "PEAI", "POCA")
inter.bet <- lapply(inter.mat, function(betas) {
  betas <- betas[sub3, sub3, ]})

# set up plot! 
pl <- list()
for (i in 1:3) { # rows
  pl.j <-vector("list", length = 3)
  for (j in 1:3) { # columns
    
    s1 <- inter.bet[[1]][i, j, ]
    s1 <- s1[s1 > quantile(s1, 0.1) & s1 < quantile(s1, 0.9)]
    s2 <- inter.bet[[2]][i, j, ]
    s2 <- s2[s2 > quantile(s2, 0.1) & s2 < quantile(s2, 0.9)]
    s3 <- inter.bet[[3]][i, j, ]
    s3 <- s3[s3 > quantile(s3, 0.1) & s3 < quantile(s3, 0.9)]
    Samples <- c(s1, s2, s3)
    # Samples <- c(inter.bet[[1]][i, j, ], inter.bet[[2]][i, j, ], inter.bet[[3]][i, j, ])
    Comm <-  c(rep('1.Open', 800), rep('2.Inter', 800), rep('3.Shady', 800))
    df <- list(Comm, Samples)
    df <- as.data.frame(df, col.names = c('Community', 'Samples'))
    
    df2 <- group_by(df, Community) %>% summarise('Samples' = median(Samples))
    
    g <- ggplot(df, aes(x=Samples, y=Community, 
                        fill=as.factor(Community)))+
      scale_fill_manual(values=c('#fbc891', '#92dcaa', '#089941')) +
      geom_vline(xintercept=0, linetype = 'dashed') +
      geom_density_ridges(scale = 2, alpha=0.5) +
      {if(j==1) labs(y =sub3[i]) else labs(y = NULL)} +
      scale_y_discrete(#name = species_3[j], 
        expand=c(0.01, 0)) +
      scale_x_continuous(name = NULL, 
             #            limits = c(min(unlist(inter.bet)), max(unlist(inter.bet))),
                         expand=c(0.01, 0)) +
      geom_point(data = df2, col = 'red', cex = 4) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_text(size = 25),
            axis.text.x = element_text(size=13, colour = 'black'),
            legend.position='none',
            plot.title = element_text(hjust = 0.5, size = 25),
            plot.caption = element_text(hjust = 0, size = 25),
            plot.caption.position =  "panel") +
      {if(i==1) labs(title = sub3[j])}
    
    pl.j[[j]] <- g
  }
  pl <- c(pl, pl.j)
}

png('figures/betas_vary.png', width=650, height=700)
plot_grid(plotlist = pl, nrow = 3, align = 'hv')
dev.off()


#-------------------------------------------------------------------------------------------
# FIGURE 4: TERNARY PLOTS
#-------------------------------------------------------------------------------------------

grad <- paletteer_c("grDevices::Tropic", 30) 
grad <- grad[15:30]

smets <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)
# replace with line below if not running your own model or analyses
# smets <- read.csv('uploaded_results/species_metrics_mat.csv', stringsAsFactors = F)

smm <- smets[is.na(smets$aii) == F, ]
dff <- smm[ , c('spcID', 'species', 'com', 'perc.coop', 'perc.fullcomp')]
dff <- cbind(dff, smm$perc.exploited + smm$perc.exploiter)
colnames(dff) <- c('spcID', 'species', 'com', 'Cooperative', 'Competitive', 'Antagonistic')
dff <- dff[ , c(1:3, 6, 4:5)]  # rearrange the colums so that antagonistic comes first
dff %>% split(., as.factor(dff$spcID)) %>%
  lapply(., function(spc) {c(apply(spc[ ,4:6], 2, median), unique(spc$com), unique(spc$species))}) %>%
  do.call(rbind, .) %>% as.data.frame(., stringsAsFactors = F) -> spcmediansloops
spcmediansloops[ , 1:4] <- apply(spcmediansloops[ , 1:4], 2, as.numeric)
idsp <- names(table(spcmediansloops$V5)[table(spcmediansloops$V5)>2])
spcmediansloops <- spcmediansloops[spcmediansloops$V5 %in% idsp, ]
spcmediansloops %>% group_by(V5) %>% summarise('var' = var(Competitive)) -> x1
# create a null at 25 / 25/ 50 
null <- c(0.5, 0.25, 0.25)

png('figures/species_triangle_plots.png', width = 743, height = 457)
par(mfrow=c(1, 2),
    mar=c(0,0,3,4),
    oma=c(0,1.,0,1),
    cex = 1, cex.main=1.5, cex.lab = 3)
sp1 <- spcmediansloops[spcmediansloops$V5 == 'PEAI', ]
TernaryPlot(alab = 'Asymmetric    \u27F6', blab = 'Cooperative    \u27F6', 
            clab = '\u27F5   Competitive',
            grid.minor.lines=5, grid.col='white', lab.cex = 1.5, axis.cex = 1)
ColourTernary(TernaryDensity(dff[dff$species == 'PEAI', 4:6], resolution = 10L), 
              spectrum = grad)
AddToTernary(points, null, col = 'black', bg = 'darkgrey', pch = 21, cex = 2)
AddToTernary(points, sp1[ , 1:3], col = 'black',
             pch=c(8, 16, 3), cex=2, lwd = 3
         #    bg=rgb(.5,.5,.5, 180, maxColorValue=255)
)
title(main = expression(italic('Pentameris airoides ') (PEAI)))
sp1 <- spcmediansloops[spcmediansloops$V5 == 'POCA', ]
TernaryPlot(alab = 'Asymmetric    \u27F6', blab = 'Cooperative    \u27F6', 
            clab = '\u27F5   Competitive',
            grid.minor.lines=5, grid.col='white', lab.cex = 1.5, axis.cex = 1)
ColourTernary(TernaryDensity(dff[dff$species == 'POCA', 4:6], resolution = 10L), 
              spectrum = grad)
AddToTernary(points, null, col = 'black', bg = 'darkgrey', pch = 21, cex = 2)
AddToTernary(points, sp1[ , 1:3], #pch = 17, cex = 2.5
             pch=c(8, 16, 3), cex=2, col = 'black', lwd = 3
          #   bg=rgb(.5,.5,.5, 180, maxColorValue=255)
)
title(main = expression(italic('Podolepsis canescens ') (POCA)))
vrange <- range(TernaryDensity(dff[dff$species == 'POCA' | dff$species == 'PEAI', 4:6], resolution = 10L))
gradientLegend(valRange = round(vrange), color = grad)
dev.off()


#-------------------------------------------------------------------------------------------
# FIGURE 5: COMMUNITY NETWORKS & METRICS
#-------------------------------------------------------------------------------------------

# 1 - NETWORKS WITH QGRAPH
#-------------------------
inter.mat <- list()
inter.mat <- sapply(c(1,2,3), function(c){
  load(paste0('model/transformed/', c, '/scaled_betaij_matrices.Rdata'))
  # replace with line below if not running your own models
  # load(paste0('uploaded_results/', c, '/scaled_betaij_matrices.Rdata'))
  inter.mat[[as.numeric(c)]] <- scaled.betas
})
smets <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)
# replace with line below if not running your own model or analyses
# smets <- read.csv('uploaded_results/species_metrics_mat.csv', stringsAsFactors = F)
smm <- smets[is.na(smets$aii) == F, ]

netwks <- lapply(inter.mat, apply, c(1,2), median)
netwks <- lapply(netwks, function(ntw) {ntw <- ntw[ , 1:nrow(ntw)]}) # keep focals only
# this below is super ugly but will do for now
netwks <- lapply(netwks, function(ntw) {
  emptymat <- matrix(data = rep(0, 21*21), nrow = 21, ncol = 21, 
                     dimnames = list('species' = unique(smm$species)[order(unique(smm$species))], 
                                     'neighbour' = unique(smm$species)[order(unique(smm$species))]))
  for(i in row.names(ntw)) {
    for (j in colnames(ntw)) {
      emptymat[i, j] <- ntw[i, j]
    }
  }
  return(emptymat)
})

max.edge <- max(unlist(lapply(netwks, max)))  # have a meximum edge weight to comapre networks

png('figures/networks.png', 
    width = 600, height = 1200, units = 'px')
par(mfrow=c(3, 1))
lapply(netwks, function(ntw) {
  ntw <- t(ntw) # transpose so that arrows point towards i
  # set colours for species that repeat across networks
  all.sp <- rep('white', length(dimnames(ntw)$species))
  names(all.sp) <- dimnames(ntw)$species
  all.sp[species_3] <- 'grey60'
  #set white border for species absent from a network
  sp.in.ntw <- colSums(ntw)
  sp.in.ntw <- names(sp.in.ntw[sp.in.ntw > 0 | sp.in.ntw < 0])
  brd.col <- rep('white', length(dimnames(ntw)$species))
  names(brd.col) <- colnames(ntw)
  brd.col[which(names(brd.col) %in% sp.in.ntw)] <- 'black'
  # labels = dimnames(ntw)$species   # toggle on if wanting no labels for absent species
  # labels[!labels %in% sp.in.ntw] <- NA
  qgraph(ntw,  # plot interaction means
         #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
         layout = 'circle',
         edge.width = 2,
         negCol = '#78d0fd',   # facilitation = blue
         posCol = '#fda578',          
         vsize = 12,  # node size
         vTrans = 100,  # node transparency
         labels = dimnames(ntw)$species,
         label.scale = F,
         label.cex = 1.8,
         color = all.sp,
         border.color = brd.col,
         fade = T, directed = T,
         diag = T,
         maximum = max.edge,
         asize = 6   # arrow size
         #title = '', title.cex =5
  )
})
dev.off()

# 2 - NETWORK METRICS WITH SPIDER PLOTS
#--------------------------------------
require(fmsb)

cmm <- read.csv('analyses/coms_metrics_mat.csv', stringsAsFactors = F)
# replace with line below if not running your own model or analyses
# cmm <- read.csv('uploaded_results/coms_metrics_mat.csv', stringsAsFactors = F)

# Do one radar chart for each community, with sd of each var 

# modifying reported variables slightly
upres <- lapply(inter.mat, function(tens) {
  foo <- t(apply(tens, 3, function(net) {
    # % of intras that are F
    intras <- diag(net)
    piiF <- length(intras[intras<0])/length(intras)
    # % of inters that are C and F
    diag(net) <- 0
    pijF <- length(net[net<0])/(length(net)-nrow(net))
    pijC <- length(net[net>0])/(length(net)-nrow(net))
    return(c(piiF, pijF, pijC))
  }))
  colnames(foo) <- c('piiF', 'pijF', 'pijC')
  return(foo)
})
upres <- do.call(rbind, upres)
upres <- cbind(upres, as.integer(rownames(upres)), c(rep(1, 1000), rep(2, 1000), rep(3, 1000)))
colnames(upres)[4:5] <- c('samples', 'comm')
upres <- as.data.frame(upres)

cmm$samples <- rep(seq(1,1000), 3)

cmm <- merge(cmm, upres, by = c('comm', 'samples'))

cmm.vars <- cmm[ , c('RI', 'mod', 'piiF', 'pijF', 'pijC', 'wc', 'comm')]
colnames(cmm.vars) <- c('Hierarchy\n(RI index)', 'Modularity', '% Self-facilitation', 
                        '% Interspecific\nFacilitation', '% Interspecific\nCompetition',
                        'Weighted\nconnectance', 'comm')
max.min <- rbind(apply(cmm.vars, 2, max), 
                 apply(cmm.vars, 2, min))
max.min[1, ] <- 1
max.min[2, ] <- 0
cmm.vars <- split(cmm.vars, cmm.vars$comm)

q10 <- do.call(rbind, lapply(cmm.vars, apply, 2, quantile, .1))
q90 <- do.call(rbind, lapply(cmm.vars, apply, 2, quantile, .9))
medians <- do.call(rbind, lapply(cmm.vars, apply, 2, median))
radata2 <- as.data.frame(rbind(q10, medians, q90))
radata2 <- split(radata2, as.factor(radata2[ , 'comm']))
radata2 <- lapply(radata2, function(x) { 
  rbind(max.min, x)
})

# Plot! 
png('figures/spider_plot.png', height = 1800, width = 700)

par(mfrow=c(3, 1),
    cex = 1.2,
    mar=c(0,0,3,0),
    oma=c(0,1,0,1),
    mgp = c(2, 1, 0))

colours_in=c(rgb(0.984,0.784,0.569,0.5), rgb(0.573,0.863,0.667,0.5) , rgb(0.031,0.6,0.255,0.5) )
colours_border=c(rgb(0.965,0.4,0.086,1), rgb(0.227,0.706, 0.384,1) , rgb(0.016,0.314,0.133,1) )

lapply(radata2, function(rad){
  # get community number
  com <- as.numeric(names(table(rad$comm)[which.max(table(rad$comm))]))
  # and name 
  comname <- c("Open", "Intermediate", "Shady")
  
  radarchart(rad[ , 1:6],
             cglty = 1,  # axis line type
             cglcol = 'grey', # axis line colour
             cglwd=2,  # axis line width
             pty = c(NA, 16, NA),
             # custom polygon: 
             pcol=colours_border[com], # polygon line colours
             pfcol=c(NA, colours_in[com], NA), # only colour the mean polygon in
             plwd=c(3,4,3), # polygon line width
             plty= c(2,1,2), # polygon line type
             title = comname[com],
             vlcex = 1.3,  # axis label size
             cex.main=1.5
  )
})
dev.off()


# 3 - JOINT NETWORKS AND METRICS FIGURE
#--------------------------------------
png('figures/networks_and_metrics.png', 
    width = 1900, height = 2400, units = 'px')
par(mfrow=c(3, 2),
    mar=c(0,0,4,1))
comname <- c("Open", "Intermediate", "Shady")
#colours_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colours_in=c(rgb(0.984,0.784,0.569,0.5), rgb(0.573,0.863,0.667,0.5) , rgb(0.031,0.6,0.255,0.5) )
colours_border=c(rgb(0.965,0.4,0.086,1), rgb(0.227,0.706, 0.384,1) , rgb(0.016,0.314,0.133,1) )

# corresponds to c('#fbc891', '#92dcaa', '#089941')

lapply(c(1, 2, 3), function(com) {
  
  ntw <- netwks[[com]]
  ntw <- t(ntw) # transpose so that arrows point towards i
  # set colours for species that repeat across networks
  all.sp <- rep('white', length(dimnames(netwks[[com]])$species))
  names(all.sp) <- dimnames(netwks[[com]])$species
  all.sp[species_3] <- 'grey60'
  #set white border for species absent from a network
  sp.in.ntw <- colSums(ntw)
  sp.in.ntw <- names(sp.in.ntw[sp.in.ntw > 0 | sp.in.ntw < 0])
  brd.col <- rep('white', length(dimnames(ntw)$species))
  names(brd.col) <- colnames(ntw)
  brd.col[which(names(brd.col) %in% sp.in.ntw)] <- 'black'
  # and lower transparency
  transp <- brd.col
  transp[transp == 'white'] <- 100
  transp[transp == 'black'] <- 200
  transp <- as.numeric(transp)
  
  qgraph(ntw,  # plot interaction means
         #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
         edge.width = 2,
         layout = 'circle',
         negCol = '#78d0fd',   # facilitation = blue
         posCol = '#fda578',     # competition = orange     # competition = orange
         vsize = 12,  # node size
         vTrans = transp,  # node transparency
         labels = dimnames(ntw)$species,
         label.scale = F,
         label.cex = 2.7,
         color = all.sp,
         border.color = brd.col,
         fade = T, directed = T,
         diag = T,
         maximum = max.edge,
         asize = 6
         #   title = '', title.cex =5
  )
  
  radarchart(radata2[[com]][ , 1:6],
             cglty = 1,  # axis line type
             cglcol = 'grey', # axis line colour
             cglwd=2,  # axis line width
             pty = c(NA, 16, NA),
             # custom polygon: 
             pcol=colours_border[com], # polygon line colours
             pfcol=c(NA, colours_in[com], NA), # only colour the mean polygon in
             plwd=c(3,4,3), # polygon line width
             plty= c(2,1,2), # polygon line type
             #  title = comname[com],
             vlcex = 3.5,  # axis label size
             #  cex.main=4
  )
  title(list(comname[com], cex=6), adj = 0)
})
dev.off()

###########################################################################################
################################# SUPPLEMENTARY FIGURES ###################################
###########################################################################################

#-------------------------------------------------------------------------------------------
# SUPPS: POSTERIOR PREDICTIVE CHECK PLOTS
#-------------------------------------------------------------------------------------------

# Posterior Predictive checks
source('functions/figures/post_pred_checks.R')
# NB: This function will only work if you have run the stan models following model.R

post_pred_fig(1, 'mu2')
post_pred_fig(2, 'mu2')
post_pred_fig(3, 'mu2')

#-------------------------------------------------------------------------------------------
# SUPPS: COMPETITIVE VS FACIL OUTPUT SPECIES BY SPECIES 
#-------------------------------------------------------------------------------------------

smm <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)
# replace with line below if not running your own model or analyses
# smets <- read.csv('uploaded_results/species_metrics_mat.csv', stringsAsFactors = F)
smm <- smm[ , c('species', 'com', 'spcID', 'sum_aij', 'C_sum_aji', 'F_sum_aji')]
smm$F_sum_aji <- - smm$F_sum_aji

# This sets individual colours for each species
ell.col <- c(paletteer_d("ggthemes::Classic_20"), '#B31166')
names(ell.col) <- sort(unique(smm$species))

smm %>% split(., as.factor(smm$spcID)) %>%
  lapply(., function(spc) {c(apply(spc[ ,4:6], 2, median), 
                             unique(spc$com), unique(spc$species))}) %>%
  do.call(rbind, .) %>% as.data.frame(., stringsAsFactors = F) -> spcmediansCF
spcmediansCF[ , 1:3] <- apply(spcmediansCF[ , 1:3], 2, as.numeric)
colnames(spcmediansCF)[4:5] <- c('com', 'species')
xlow <- min(smm[ , c('C_sum_aji', 'F_sum_aji')])
xmax <- max(smm[ , c('C_sum_aji', 'F_sum_aji')])

png('figures/split_effects_out.png', width = 1240, height = 1754)
par(mfrow=c(7, 3),
    mai = c(.4,.4,.4,.4), mar = c(5,5,2,1), oma=c(1,1,2,1),
    cex.lab=2, cex.main=2.5)
spltmed <- split(spcmediansCF, as.factor(spcmediansCF$species))
spcid <- names(spltmed)
compt <- c(8, 16, 3)
# xlow <- min(smm[ , c('C_sum_aji', 'F_sum_aji')])
# xmax <- max(smm[ , c('C_sum_aji', 'F_sum_aji')])
for(i in 1:21) {
  spdf <- spltmed[[i]]
  xlow <- min(smm[smm$species == spcid[i] , c('C_sum_aji', 'F_sum_aji')])
  xmax <- max(smm[smm$species == spcid[i] , c('C_sum_aji', 'F_sum_aji')])
  plot(spdf[ , 'F_sum_aji'] ~ spdf[ , 'C_sum_aji'], 
       xlab = 'Competitive output',
       ylab = 'Facilitative output',
       xlim = c(xlow, xmax), 
       ylim = c(xlow, xmax),
       cex.axis = 1.5,
       type='n', bty='n')
  points(smm[smm$species == spcid[i] , 'F_sum_aji'] ~ smm[smm$species == spcid[i] , 'C_sum_aji'], 
         pch=20, cex = 1.5, col = adjustcolor(ell.col[i], alpha.f=0.5))
  abline(0, 1, lty = 2)
  cc <- spdf$com
  for (c in 1:length(cc)) {
    ellipse(mu = c(spdf[spdf$com == cc[c], 'C_sum_aji'], spdf[spdf$com == cc[c], 'F_sum_aji']), 
            sigma = cov(smm[smm$species == spcid[i] & smm$com == cc[c], c('C_sum_aji', 'F_sum_aji')]),
            alpha = 0.2, # prob to be excluded
            col = 'black')
    points(spdf[spdf$com == cc[c] , 'F_sum_aji'] ~ spdf[spdf$com == cc[c]  , 'C_sum_aji'], 
           pch=compt[as.numeric(cc[c])], cex = 2, lwd = 3)
  }
  if(spcid[i] %in% species_3) title(paste0(spcid[i], '*')) else title(spcid[i])
}

dev.off()

#-------------------------------------------------------------------------------------------
# SUPPS: COMPETITIVE VS FACIL INPUT SPECIES BY SPECIES 
#-------------------------------------------------------------------------------------------

smets <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)
# replace with line below if not running your own model or analyses
# smets <- read.csv('uploaded_results/species_metrics_mat.csv', stringsAsFactors = F)
smm <- smets[is.na(smets$aii) == F, ]
smm <- smm[ , c('species', 'com', 'spcID', 'aii', 'sum_aji', 'C_sum_aij', 'F_sum_aij')]
smm$F_sum_aij <- - smm$F_sum_aij

smm %>% split(., as.factor(smm$spcID)) %>%
  lapply(., function(spc) {c(apply(spc[ ,4:7], 2, median),
                             unique(spc$com), unique(spc$species))}) %>%
  do.call(rbind, .) %>% as.data.frame(., stringsAsFactors = F) -> spcmedians
spcmedians[ , 1:4] <- apply(spcmedians[ , 1:4], 2, as.numeric)
colnames(spcmedians)[5:6] <- c('com', 'species')

# This sets individual colours for each species
ell.col <- c(paletteer_d("ggthemes::Classic_20"), '#B31166')
names(ell.col) <- sort(unique(smm$species))

png('figures/split_effects_in.png', width = 1240, height = 1754)
par(mfrow=c(7, 3),
    mai = c(.4,.4,.4,.4), mar = c(5,5,2,1), oma=c(1,1,2,1),
    cex.lab=2, cex.main=2.5)
spltmed <- split(spcmedians, as.factor(spcmedians$species))
spcid <- names(spltmed)
compt <- c(8, 16, 3)
# xlow <- min(smm[ , c('C_sum_aij', 'F_sum_aij')])
# xmax <- max(smm[ , c('C_sum_aij', 'F_sum_aij')])
for(i in 1:21) {
  spdf <- spltmed[[i]]
  xlow <- min(smm[smm$species == spcid[i] , c('C_sum_aij', 'F_sum_aij')])
  xmax <- max(smm[smm$species == spcid[i] , c('C_sum_aij', 'F_sum_aij')])
  plot(spdf[ , 'F_sum_aij'] ~ spdf[ , 'C_sum_aij'], 
       xlab = 'Competitive input',
       ylab = 'Facilitative input',
       xlim = c(xlow, xmax), 
       ylim = c(xlow, xmax),
       cex.axis = 1.5,
       type='n', bty='n')
  points(smm[smm$species == spcid[i] , 'F_sum_aij'] ~ smm[smm$species == spcid[i] , 'C_sum_aij'], 
         pch=20, cex = 1.5, col = adjustcolor(ell.col[i], alpha.f=0.5))
  abline(0, 1, lty = 2)
  cc <- spdf$com
  for (c in 1:length(cc)) {
    ellipse(mu = c(spdf[spdf$com == cc[c], 'C_sum_aij'], spdf[spdf$com == cc[c], 'F_sum_aij']), 
            sigma = cov(smm[smm$species == spcid[i] & smm$com == cc[c], c('C_sum_aij', 'F_sum_aij')]),
            alpha = 0.2, # prob to be excluded
            col = 'black')
    points(spdf[spdf$com == cc[c] , 'F_sum_aij'] ~ spdf[spdf$com == cc[c]  , 'C_sum_aij'], 
           pch=compt[as.numeric(cc[c])], cex = 2, lwd = 3)
  }
  if(spcid[i] %in% species_3) title(paste0(spcid[i], '*')) else title(spcid[i])
}

dev.off()

#-------------------------------------------------------------------------------------------
# SUPPS: TERNARY PLOTS FOR ALL OTHER SPECIES
#-------------------------------------------------------------------------------------------

grad <- paletteer_c("grDevices::Tropic", 30) 
grad <- grad[15:30]
# triangle plots for all the other species 

smm <- read.csv('analyses/species_metrics_mat.csv', stringsAsFactors = F)
# replace with line below if not running your own model or analyses
# smets <- read.csv('uploaded_results/species_metrics_mat.csv', stringsAsFactors = F)
dff <- smm[ , c('spcID', 'species', 'com', 'perc.coop', 'perc.fullcomp')]
dff <- cbind(dff, smm$perc.exploited + smm$perc.exploiter)
colnames(dff) <- c('spcID', 'species', 'com', 'Cooperative', 'Competitive', 'Antagonistic')
dff <- dff[ , c(1:3, 6, 4:5)]  # rearrange the colums so that antagonistic comes first
dff %>% split(., as.factor(dff$spcID)) %>%
  lapply(., function(spc) {c(apply(spc[ ,4:6], 2, median), unique(spc$com), unique(spc$species))}) %>%
  do.call(rbind, .) %>% as.data.frame(., stringsAsFactors = F) -> spcmediansloops
spcmediansloops[ , 1:4] <- apply(spcmediansloops[ , 1:4], 2, as.numeric)
spcmediansloops %>% group_by(V5) %>% summarise('var' = var(Competitive)) -> x1
# Remove the two species that are plotted in the main text
spcmediansloops <- spcmediansloops[!spcmediansloops$V5 == 'PEAI', ]
spcmediansloops <- spcmediansloops[!spcmediansloops$V5 == 'POCA', ]
spctoplot <- unique(spcmediansloops$V5)
# create a null at 25 / 25/ 50 
null <- c(0.5, 0.25, 0.25)

png('figures/triplots_supps.png', width = 1240, height = 1754)
par(mfrow=c(5, 4),
    mar=c(0,0,3,4),
    oma=c(0,1,1,1),
    cex = 1, cex.main=2, cex.lab = 3)
for(i in 1:length(spctoplot)) {
  sp1 <- spcmediansloops[spcmediansloops$V5 == spctoplot[i], ]
  TernaryPlot(alab = 'Assymetric    \u27F6', blab = 'Cooperative    \u27F6', 
              clab = '\u27F5   Competitive',
              grid.minor.lines=5, grid.col='white', lab.cex = 1.5, axis.cex = 1)
  ColourTernary(TernaryDensity(dff[dff$species == spctoplot[i], 4:6], resolution = 10L), 
                spectrum = grad)
  AddToTernary(points, null, col = 'black', bg = 'darkgrey', pch = 21, cex = 2)
  AddToTernary(points, sp1[ , 1:3], #pch = 17, cex = 2.5
               col = 'black',
               pch=c(8, 16, 3), cex=2, lwd = 3
  )
  if(spctoplot[i] %in% species_3) title(paste0(spctoplot[i], '*')) else title(spctoplot[i])
  vrange <- range(TernaryDensity(dff[ , 4:6], resolution = 10L))
  gradientLegend(valRange = round(vrange), color = grad)
  
}

dev.off()

#-------------------------------------------------------------------------------------------
# SUPPS: SIGNIFICANT INTERACTIONS ONLY - NETWORK
#-------------------------------------------------------------------------------------------
inter.mat <- list()
inter.mat <- sapply(c(1,2,3), function(c){
  load(paste0('model/transformed/', c, '/scaled_betaij_matrices.Rdata'))
  # replace with line below if not running your own models
  # load(paste0('uploaded_results/', c, '/scaled_betaij_matrices.Rdata'))
  inter.mat[[as.numeric(c)]] <- scaled.betas
})
speciesnames <- read.csv('clean_data/key_speciesID0.csv', stringsAsFactors = F)$x
speciesnames <- speciesnames[-3]  # remove EROD

med.net <- lapply(inter.mat, apply, c(1,2), median)
# upper and lower quantiels of each interactions
q10 <- lapply(inter.mat, apply, c(1,2), quantile, 0.1)
q10 <- lapply(q10, sign)
q90 <- lapply(inter.mat, apply, c(1,2), quantile, 0.9)
q90 <- lapply(q90, sign)
# add them up, if 2 or -2 then upper and lower quantiles are of the same sign
sigs <- mapply('+', q10, q90)
sigs <- lapply(sigs, abs)
sigs <- lapply(sigs, function(x) x/2)
sig.ntwk <- mapply('*', med.net, sigs)

# add species that don't appear in that com
sig.ntwk <- lapply(sig.ntwk, function(ntw) {
  emptymat <- matrix(data = rep(0, 21*21), nrow = 21, ncol = 21, 
                     dimnames = list('species' = speciesnames, 
                                     'neighbour' = speciesnames))
  for(i in row.names(ntw)) {
    for (j in colnames(ntw)) {
      emptymat[i, j] <- ntw[i, j]
    }
  }
  return(emptymat)
})

max.edge <- max(unlist(lapply(sig.ntwk, max)))  # have a meximum edge weight to comapre networks
ntw.name <- c('Open', 'Intermediate', 'Shady')
species_3 <- c("ARCA", "HYGL", "PEAI", "PLDE", "POCA", "PTGA", "STPA", "VERO", "WAAC")

png('figures/networks_nonzero.png', 
    width = 600, height = 1200, units = 'px')
par(mfrow=c(3, 1))
lapply(sig.ntwk, function(ntw) {
  ntw <- t(ntw) # transpose so that arrows point towards i
  # set colours for species that repeat across networks
  all.sp <- rep('white', length(dimnames(ntw)$species))
  names(all.sp) <- dimnames(ntw)$species
  all.sp[species_3] <- 'grey60'
  #set white border for species absent from a network
  sp.in.ntw <- colSums(ntw)
  sp.in.ntw <- names(sp.in.ntw[sp.in.ntw > 0 | sp.in.ntw < 0])
  brd.col <- rep('white', length(dimnames(ntw)$species))
  names(brd.col) <- colnames(ntw)
  brd.col[which(names(brd.col) %in% sp.in.ntw)] <- 'black'
  qgraph(ntw,  # plot interaction medians
         #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
         layout = 'circle',         
         edge.width = 2,
         negCol = '#78d0fd',   # facilitation = blue
         posCol = '#fda578',  
         vsize = 12,  # node size
         vTrans = 100,  # node transparency
         labels = dimnames(ntw)$species,
         label.scale = F,
         label.cex = 1.8,
         color = all.sp,
         border.color = brd.col,
         fade = T, directed = T,
         diag = T,
         maximum = max.edge,
         asize = 6  # arrow size
         #title = '', title.cex =5
  )
})
dev.off()

#-------------------------------------------------------------------------------------------
# SUPPS: BETA DISTRIBUTIONS COM TO COM 
#-------------------------------------------------------------------------------------------

# load scaled interactions
inter.mat <- list()
inter.mat <- sapply(c(1,2,3), function(c){
  load(paste0('model/transformed/', c, '/scaled_betaij_matrices.Rdata'))
  # replace with line below if not running your own models
  # load(paste0('uploaded_results/', c, '/scaled_betaij_matrices.Rdata'))
  inter.mat[[as.numeric(c)]] <- scaled.betas
})
# species which appear in all 3 communities 
species_3 <- c("ARCA", "HYGL", "PEAI", "PLDE", "POCA", "PTGA", "STPA", "VERO", "WAAC")
# select only interactions which appear in all 3 coms 
inter.9 <- lapply(inter.mat, function(betas) {
  betas <- betas[species_3, species_3, ]})

# set up plot! 
pl <- list()

for (i in 1:9) { # rows
  pl.j <-vector("list", length = 9)
  for (j in 1:9) { # columns
    
    s1 <- inter.9[[1]][i, j, ]
    s1 <- s1[s1 > quantile(s1, 0.1) & s1 < quantile(s1, 0.9)]
    s2 <- inter.9[[2]][i, j, ]
    s2 <- s2[s2 > quantile(s2, 0.1) & s2 < quantile(s2, 0.9)]
    s3 <- inter.9[[3]][i, j, ]
    s3 <- s3[s3 > quantile(s3, 0.1) & s3 < quantile(s3, 0.9)]
    Samples <- c(s1, s2, s3)
    # Samples <- c(inter.9[[1]][i, j, ], inter.9[[2]][i, j, ], inter.9[[3]][i, j, ])
    Comm <-  c(rep('1.Open', 800), rep('2.Inter', 800), rep('3.Shady', 800))
    df <- list(Comm, Samples)
    df <- as.data.frame(df, col.names = c('Community', 'Samples'))
    
    df2 <- group_by(df, Community) %>% summarise('Samples' = median(Samples))
    
    g <- ggplot(df, aes(x=Samples, y=Community, 
                        fill=as.factor(Community)))+
      scale_fill_manual(values=c('#fbc891', '#92dcaa', '#089941')) +
      geom_vline(xintercept=0, linetype = 'dashed') +
      geom_density_ridges(scale = 1.5, alpha=0.5) +
      {if(j==1) labs(y = species_3[i]) else labs(y = NULL)} +
      scale_y_discrete(#name = species_3[j], 
                       expand=c(0.01, 0)) +
      scale_x_continuous(name = NULL, expand=c(0.01, 0)) +
      geom_point(data = df2, col = 'red', cex = 4) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_text(size = 40),
            legend.position='none',
            plot.title = element_text(hjust = 0.5, size = 40),
            plot.caption = element_text(hjust = 0, size = 40),
            plot.caption.position =  "panel") +
      {if(i==1) labs(title = species_3[j])}
    
    pl.j[[j]] <- g
  }
  pl <- c(pl, pl.j)
}

png('figures/betas_across_comms_density.png', width=1800, height=1800)
#par(oma=c(1,5,5,1))
plot_grid(plotlist = pl, nrow = 9, align = 'hv')
dev.off()

#-------------------------------------------------------------------------------------------
# SUPPS: CORRELATIONS BETWEEN NETWORK METRICS 
#-------------------------------------------------------------------------------------------
cmm <- read.csv('analyses/coms_metrics_mat.csv')
# replace with line below if not running your own model or analyses
# cmm <- read.csv('uploaded_results/coms_metrics_mat.csv', stringsAsFactors = F)

inter.mat <- list()
inter.mat <- sapply(c(1,2,3), function(c){
  load(paste0('model/transformed/', c, '/scaled_betaij_matrices.Rdata'))
  # replace with line below if not running your own models
  # load(paste0('uploaded_results/', c, '/scaled_betaij_matrices.Rdata'))
  inter.mat[[as.numeric(c)]] <- scaled.betas
})
upres <- lapply(inter.mat, function(tens) {
  foo <- t(apply(tens, 3, function(net) {
    # % of intras that are F
    intras <- diag(net)
    piiF <- length(intras[intras<0])/length(intras)
    # % of inters that are C and F
    diag(net) <- 0
    pijF <- length(net[net<0])/(length(net)-nrow(net))
    pijC <- length(net[net>0])/(length(net)-nrow(net))
    return(c(piiF, pijF, pijC))
  }))
  colnames(foo) <- c('piiF', 'pijF', 'pijC')
  return(foo)
})
upres <- do.call(rbind, upres)
upres <- cbind(upres, as.integer(rownames(upres)), c(rep(1, 1000), rep(2, 1000), rep(3, 1000)))
colnames(upres)[4:5] <- c('samples', 'comm')
upres <- as.data.frame(upres)

cmm$samples <- rep(seq(1,1000), 3)

cmm <- merge(cmm, upres, by = c('comm', 'samples'))
cmm.vars <- cmm[ , c('RI', 'mod', 'piiF', 'pijF', 'pijC', 'wc', 'comm')]
# get the medians
cmm.vars.med <- split(cmm.vars, cmm.vars$comm)
cmm.vars.med <- as.data.frame(do.call(rbind, lapply(cmm.vars.med, apply, 2, median)))
cmm.vars.med$points <- c(8, 16, 3) # different points for each networks median
cmm.vars.med$col <- rep('black', 3)
# add point type and colour
cmm.vars$points <- rep(20, dim(cmm.vars)[1])  
comcol <- c('#fbc891', '#92dcaa', '#089941')
cmm.vars$col <- comcol[cmm.vars$comm] 
colnames(cmm.vars) <- c('Hierarchy\n(RI index)', 'Modularity', '% Self-facilitation', 
                        '% Interspecific\nFacilitation', '% Interspecific\nCompetition',
                        'Weighted\nconnectance', 'comm', 'points', 'col')
# all together
colnames(cmm.vars.med) <- colnames(cmm.vars)
cmm.vars <- rbind(cmm.vars, cmm.vars.med)

png('figures/comm_pairs.png', width = 1200, height = 1200)
pairs(cmm.vars[ ,1:6], pch=cmm.vars$points, cex = 1.5, 
      col = adjustcolor(cmm.vars$col, alpha.f= 0.75))
dev.off()


