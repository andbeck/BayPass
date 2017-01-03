# setups
library(corrplot)
library(ape)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggthemes)
library(qqman)

source('baypass_utils.R') # this is in the repo

# Covariance matrices
omegaNoB1<-as.matrix(read.table("/Volumes/TTYLMF/ana.noB1__mat_omega.out"))
omegaNoBag<-as.matrix(read.table("/Volumes/TTYLMF/ana.noBag__mat_omega.out"))
omegaNoD8<-as.matrix(read.table("/Volumes/TTYLMF/ana.noD8__mat_omega.out"))
omegaNoD10<-as.matrix(read.table("/Volumes/TTYLMF/ana.noD10__mat_omega.out"))
omegaNoNM2<-as.matrix(read.table("/Volumes/TTYLMF/ana.noNM2__mat_omega.out"))
omegaNoStav<-as.matrix(read.table("/Volumes/TTYLMF/ana.noStav__mat_omega.out"))
omegaNoW1<-as.matrix(read.table("/Volumes/TTYLMF/ana.noW1__mat_omega.out"))
omegaNoW6<-as.matrix(read.table("/Volumes/TTYLMF/ana.noW6__mat_omega.out"))

pops<-c("B1","Bag","D8","D10","NM2","Stav","W1","W6")
labs<-combn(rev(pops),7)

omegas<-list(
	omegaNoB1<-matrix(omegaNoB1, 7,7, byrow = TRUE, dimnames=list(rev(labs[,1]),rev(labs[,1]))),
	omegaNoBag<-matrix(omegaNoBag, 7,7, byrow = TRUE, dimnames=list(rev(labs[,2]),rev(labs[,2]))),
	omegaNoD8<-matrix(omegaNoD8, 7,7, byrow = TRUE, dimnames=list(rev(labs[,3]),rev(labs[,3]))),
	omegaNoD10<-matrix(omegaNoD10, 7,7, byrow = TRUE, dimnames=list(rev(labs[,4]),rev(labs[,4]))),
	omegaNoNM2<-matrix(omegaNoNM2, 7,7, byrow = TRUE, dimnames=list(rev(labs[,5]),rev(labs[,5]))),
	omegaNoStav<-matrix(omegaNoStav, 7,7, byrow = TRUE, dimnames=list(rev(labs[,6]),rev(labs[,6]))),
	omegaNoW1<-matrix(omegaNoW1, 7,7, byrow = TRUE, dimnames=list(rev(labs[,7]),rev(labs[,7]))),
	omegaNoW6<-matrix(omegaNoW6, 7,7, byrow = TRUE, dimnames=list(rev(labs[,8]),rev(labs[,8])))
)

# # This needs to be developed to establish the criteria for detecting XtX outlier
# # DETAIL provided by Gautier

# # omega read in above

poolsize<-as.numeric(read.table("/Volumes/TTYLMF/SAMPLEFILE_Jack")) #the original pool sizes

##the parmeters of the Beta distribution

pi.params<-list(
	pi.params.NoB1<-as.numeric(read.table("/Volumes/TTYLMF/ana.noB1__summary_beta_params.out",h=T)[,2]),
	pi.params.NoBag<-as.numeric(read.table("/Volumes/TTYLMF/ana.noBag__summary_beta_params.out",h=T)[,2]), 
	pi.params.NoD8<-as.numeric(read.table("/Volumes/TTYLMF/ana.noD8__summary_beta_params.out",h=T)[,2]),
	pi.params.NoD10<-as.numeric(read.table("/Volumes/TTYLMF/ana.noD10__summary_beta_params.out",h=T)[,2]),
	pi.params.NoNM2<-as.numeric(read.table("/Volumes/TTYLMF/ana.noNM2__summary_beta_params.out",h=T)[,2]),
	pi.params.NoStav<-as.numeric(read.table("/Volumes/TTYLMF/ana.noStav__summary_beta_params.out",h=T)[,2]),
	pi.params.NoW1<-as.numeric(read.table("/Volumes/TTYLMF/ana.noW1__summary_beta_params.out",h=T)[,2]), 
	pi.params.NoW6<-as.numeric(read.table("/Volumes/TTYLMF/ana.noW6__summary_beta_params.out",h=T)[,2])
) 

# the orginal read counts (needed to sample from the observed coverages)
readcounts<-list(
	readcounts_NoB1<-geno2YN("/Volumes/TTYLMF/ALLELEFILE_noB1"),
	readcounts_NoBag<-geno2YN("/Volumes/TTYLMF/ALLELEFILE_noBag"),
	readcounts_NoD8<-geno2YN("/Volumes/TTYLMF/ALLELEFILE_noD8"),
	readcounts_NoD10<-geno2YN("/Volumes/TTYLMF/ALLELEFILE_noD10"),
	readcounts_NoNM2<-geno2YN("/Volumes/TTYLMF/ALLELEFILE_noNM2"),
	readcounts_NoStav<-geno2YN("/Volumes/TTYLMF/ALLELEFILE_noStav"),
	readcounts_NoW1<-geno2YN("/Volumes/TTYLMF/ALLELEFILE_noW1"),
	readcounts_NoW6<-geno2YN("/Volumes/TTYLMF/ALLELEFILE_noW6")
)

simDat2<-list()

for(i in 1:8){
	simDat2[[i]]<-
	simulate.baypass(omegas[[i]],nsnp=25000,
		beta.pi = pi.params[[i]],
		sample.size = poolsize, 
		coverage=readcounts[[i]]$NN,
		pi.maf = 0.01,
		suffix=paste("no_",pops[i],"_pods", sep=""),
		remove.fixed.loci = F)
}

### RUN BayPass on JackSimDat.RData files.
### Might need to export each independently from simDat list.

load('/Volumes/TTYLMF/JackSimDat.RData')

for(i in 1:8){
	
}

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
temp_keyptsID<-filter(df_temp, BF.dB. > as.numeric(betaThreshold[2,2]) & XtX > XtX_threshold) # 2890
predation_keyptsID<-filter(df_predation, BF.dB. > as.numeric(betaThreshold[3,2]) & XtX > XtX_threshold) # 2890

#save(latpH_keyptsID, temp_keyptsID, predation_keyptsID, file = 'outliers.RData')

# Produce outlier figures
par(mfrow = c(1,3))
#latitude/pH
plot(BF.dB. ~ XtX, data = df_latpH, pch = '.', col = '#00000033', ylim = c(-10,60))
points(BF.dB. ~ XtX, data = latpH_keypts, pch = 21, col = 'grey')
abline(v = XtX_threshold, lty = 3, lwd = 2, col = 'red')
abline(h = betaThreshold[1,2], lty = 3, lwd = 2, col = 'green')
legend(5,60, legend = c("XtX Threshold", "BayesFactor Threshold", "Outliers"), lty=c(3,3,NA), pch = c(NA, NA, 21), col = c('red', 'green', 'grey'))
title('Latitude/pH')

# temp
plot(BF.dB. ~ XtX, data = df_temp, pch = '.', col = '#00000033', ylim = c(-10,60))
points(BF.dB. ~ XtX, data = temp_keypts, pch = 21, col = 'grey')
abline(v = XtX_threshold, lty = 3, lwd = 2, col = 'red')
abline(h = betaThreshold[2,2], lty = 3, lwd = 2, col = 'green')
legend(5,60, legend = c("XtX Threshold", "BayesFactor Threshold", "Outliers"), lty=c(3,3,NA), pch = c(NA, NA, 21), col = c('red', 'green', 'grey'))
title('Temperature')

# predation
plot(BF.dB. ~ XtX, data = df_predation, pch = '.', col = '#00000033', ylim = c(-10,60))
points(BF.dB. ~ XtX, data = predation_keypts, pch = 21, col = 'grey')
abline(v = XtX_threshold, lty = 3, lwd = 2, col = 'red')
abline(h = betaThreshold[3,2], lty = 3, lwd = 2, col = 'green')
legend(5,60, legend = c("XtX Threshold", "BayesFactor Threshold", "Outliers"), lty=c(3,3,NA), pch = c(NA, NA, 21), col = c('red', 'green', 'grey'))
title('Predation')


