# APB - 30 Nov 2015 - Setup for BayPass - removing SoF - we don't know what it is.
# UPDATED 28 April 2016 with new ENVFILE with four covariates
# preping to use the AUX model
# Updated AGAIN - now working in GitHub
# ----------------------------------------------------------------------------------------

library(gdata)

# get original allele frequency data

# load('~/Dropbox/####WorkingProjects/ResequencingAnalysisTools/BerglandFiles/CrispOutput_Process/bed.obj_alan.RData')

# NEW DATA
load('~/alan/R/Robj/datUse.clean.Rdata')

# groups by fish/midge
ac<-bed.obj$allele.counts
rn<-rownames(ac)[c(1,2,3,4,5,7,8,9)]
rn2<-paste(rn,"_2",sep="")

# organise major allele count
alleleC1 <- rbind(
	bed.obj$allele.counts[1:2,],
	bed.obj$allele.counts[4,],
	bed.obj$allele.counts[3,],
	bed.obj$allele.counts[c(5,7:9),])
# Sample Size
ss <- rbind(
	bed.obj$sample.size[1:2,],
	bed.obj$sample.size[4,],
	bed.obj$sample.size[3,],
	bed.obj$sample.size[c(5,7:9),])
# minor allele count
alleleC2 <- ss	- alleleC1

# housekeeping
rownames(alleleC1)<-rn
rownames(alleleC2)<-rn2

# checking
alleleC1[1:8,1:2]
alleleC2[1:8,1:2]
ss[1:8,1:2]

# Do it with gdata function interleave
out.temp<-interleave(alleleC1,alleleC2) # keep it like this because writing to file transposes it, which is what we want

# make the covariate File
covs<-read.csv("~/baypass_2.1/DaphBayPass/pH_Temp_Lat_8_ponds.csv")

predR<-scale(covs$code)
temp<-scale(covs$Temp)
pH<-scale(covs$pH)
lat<-scale(covs$LatDec)

envfile<-t(cbind(predR, temp, pH, lat))
envfile

# make the samplesize file
ss<-rep(100,8)

# write
write(out.temp, "~/baypass_2.1/DaphBayPass/ALLELEFILE", ncolumns=16)
write(t(envfile), "~/baypass_2.1/DaphBayPass/ENVFILE", ncolumns = 8)
write(t(ss), "~/baypass_2.1/DaphBayPass/SAMPLEFILE", ncolumns = 8)