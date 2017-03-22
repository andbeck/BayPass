# setups
library(corrplot)
library(ape)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggthemes)
library(qqman)

source('baypass_utils.R') # this is in the repo

setwd('/Volumes/TTYLMF')

# Covariance matrix
omega<-as.matrix(read.table("anaprEnvfile__mat_omega.out"))
omega<-matrix(omega, 8,8, byrow = TRUE, dimnames=list(c(paste("pop:",1:8)), c(paste("pop:",1:8))))

omega_plot<-cov2cor(omega)
par(mfrow = c(2,1))
corrplot(omega_plot,method = "square", order = 'hclust',mar = c(2,1,1,0))
plot(hclust(dist(omega_plot), 'complete'))


# # This needs to be developed to establish the criteria for detecting XtX outlier
# # DETAIL provided by Gautier

# # omega read in above
# pi.params<-as.numeric(read.table("anaprEnvfile__summary_beta_params.out",h=T)[,2]) ##the parmeters of the Beta distribution
# poolsize<-as.numeric(read.table("SAMPLEFILE")) #the original pool sizes
# readcounts<-geno2YN("ALLELEFILE") # the orginal read counts (needed to sample from the observed coverages)

# simulate.baypass(omega,nsnp=250000,
	# beta.pi = pi.params,
	# sample.size = poolsize, 
	# coverage=readcounts$NN,
	# pi.maf = 0.01,
	# suffix="pods",
	# remove.fixed.loci = F)


# The real games begin here.

## =========================================================
## read in ALL baypass output and processing
## to be use in identification of outliers and plotting
## =========================================================

load('BayPassInterpret.RData') # this should eliminate data reading and prep for 

# # Step 1: use the pods analysis to identify the 99% ile of XtX
# XtX_Emp<-read.table('anaprEnvfile__summary_pi_xtx.out', header = TRUE)
# XtX_Sim<-read.table('anaPOD__summary_pi_xtx.out', header = TRUE)

# # calcuate the threshold XtX using quantile = 0.99, for the top 1%
# XtX_threshold<-quantile(XtX_Sim$M_XtX, 0.99)
# XtX_Sig<-XtX_Emp[XtX_Emp$M_XtX>=XtX_threshold,]

# # percent exceeding threshold
# dim(XtX_Sig)[1]/dim(XtX_Emp)[1]

# # Step 2: use the pods analysis to identify the 99% ile of betas
# # The tables are indexed by covariable 1-3
# # 1 = lat/pH
# # 2 = temp
# # 3 = predation
# beta_Emp<-read.table('anaprEnvfile__summary_betai_reg.out', header = TRUE)
# beta_Sim<-read.table('anaPOD__summary_betai_reg.out', header = TRUE)

# # the B
# one<-filter(beta_Emp, COVARIABLE == '1')
# two<-filter(beta_Emp, COVARIABLE == '2')
# three<-filter(beta_Emp, COVARIABLE == '3')

# betaThreshold<-beta_Sim %>%
	# group_by(COVARIABLE) %>%
	# summarise(
		# betaThreshold = quantile(BF.dB., 0.99))

# betaThreshold$COVARIABLE<-c("Lat/pH", "Temp", "Predation")
# betaThreshold

# # scaf positions and ID to add to dfs
# scaf_pos<-read.csv('scaf_pos.csv')

# # create df of bayes factors and XtX and scaf positions
# df_latpH<-data.frame(one, XtX = XtX_Emp$M_XtX, scaf_pos)
# df_temp<-data.frame(two, XtX = XtX_Emp$M_XtX, scaf_pos)
# df_predation<-data.frame(three, XtX = XtX_Emp$M_XtX, scaf_pos)

# these are the outliers based on XtX and bayes factor for each 
latpH_keyptsID<-filter(df_latpH, BF.dB. > as.numeric(betaThreshold[1,2]) & XtX > XtX_threshold) # 2890
temp_keyptsID<-filter(df_temp, BF.dB. > as.numeric(betaThreshold[2,2]) & XtX > XtX_threshold) # 2531
predation_keyptsID<-filter(df_predation, BF.dB. > as.numeric(betaThreshold[3,2]) & XtX > XtX_threshold) # 3094

# use VennDiagram and gplots
library(VennDiagram)
library(gplots)

# create list of MRKs
All<-list(latpH = latpH_keyptsID$MRK,temp = temp_keyptsID$MRK, predation = predation_keyptsID$MRK)
# build plot
vennPlot<-venn.diagram(All, NULL, fill = c('blue','orange','green'), alpha = c(0.5,0.5,0.5), cex = 3, cat.cex = 2)
# draw plot
grid.draw(vennPlot)

#  Hmmm - bit different numbers...
# # These are common to all three
# length(Reduce(intersect, list(latpH_keyptsID$MRK, temp_keyptsID$MRK, predation_keyptsID$MRK)))
# 
# # THese are common to pairs
# length(Reduce(intersect, list(latpH_keyptsID$MRK, predation_keyptsID$MRK)))
# length(Reduce(intersect, list(temp_keyptsID$MRK, predation_keyptsID$MRK)))


#save(latpH_keyptsID, temp_keyptsID, predation_keyptsID, file = 'outliers.RData')

# Produce outlier figures
par(mfrow = c(1,3))
#latitude/pH
plot(BF.dB. ~ XtX, data = df_latpH, pch = '.', col = '#00000033', ylim = c(-10,60))
points(BF.dB. ~ XtX, data = latpH_keyptsID, pch = 21, col = 'grey')
abline(v = XtX_threshold, lty = 3, lwd = 2, col = 'red')
abline(h = betaThreshold[1,2], lty = 3, lwd = 2, col = 'green')
legend(5,60, legend = c("XtX Threshold", "BayesFactor Threshold", "Outliers"), lty=c(3,3,NA), pch = c(NA, NA, 21), col = c('red', 'green', 'grey'))
title('Latitude/pH')

# temp
plot(BF.dB. ~ XtX, data = df_temp, pch = '.', col = '#00000033', ylim = c(-10,60))
points(BF.dB. ~ XtX, data = temp_keyptsID, pch = 21, col = 'grey')
abline(v = XtX_threshold, lty = 3, lwd = 2, col = 'red')
abline(h = betaThreshold[2,2], lty = 3, lwd = 2, col = 'green')
legend(5,60, legend = c("XtX Threshold", "BayesFactor Threshold", "Outliers"), lty=c(3,3,NA), pch = c(NA, NA, 21), col = c('red', 'green', 'grey'))
title('Temperature')

# predation
plot(BF.dB. ~ XtX, data = df_predation, pch = '.', col = '#00000033', ylim = c(-10,60))
points(BF.dB. ~ XtX, data = predation_keyptsID, pch = 21, col = 'grey')
abline(v = XtX_threshold, lty = 3, lwd = 2, col = 'red')
abline(h = betaThreshold[3,2], lty = 3, lwd = 2, col = 'green')
legend(5,60, legend = c("XtX Threshold", "BayesFactor Threshold", "Outliers"), lty=c(3,3,NA), pch = c(NA, NA, 21), col = c('red', 'green', 'grey'))
title('Predation')


