# APB - 30 Nov 2015 - Setup for BayPass - removing SoF - we don't know what it is.
# UPDATED 28 April 2016 with new ENVFILE with four covariates
# preping to use the AUX model
# Updated AGAIN - now working in GitHub
# ----------------------------------------------------------------------------------------

library(gdata)
library(dplyr)
library(data.table)
library(ggplot2)

# NEW DATA
load('~/alan/R/Robj/datUse.clean.Rdata')

wrk<-filter(datUse, pop!="SoF")

wrk2<-select(wrk, 1:10)

wrk %>% 
	group_by(pop) %>%
		


alleleC1<-data.frame(wrk$dp.down)
alleleC2<-data.frame(wrk$ad.down)

alleleCounts<-cbind(alleleC1, alleleC2)
alleleCounts[1:10,]








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