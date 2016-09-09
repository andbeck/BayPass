# APB - 30 Nov 2015 - Setup for BayPass - removing SoF - we don't know what it is.
# UPDATED 28 April 2016 with new ENVFILE with four covariates
# preping to use the AUX model
# Updated AGAIN - now working in GitHub
# ----------------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(magrittr)
library(data.table)
library(foreach)

# NEW DATA
load('~/Documents/Repos/BayPass/datUse.clean.Rdata')
#load('~/Repos/BayPass/datUse.clean.Rdata')

wrk<-filter(datUse, pop!="SoF")
names(wrk)

# METHOD FOR BayPass Data setup
 
 	methodBergland <- function(dat) {
 		setkey(dat, pop)
 		
 		pops <- unique(dat$pop)
 		
 		o <- foreach(i=pops)%do%{dat[J(i), c("dp.down", "ad.down"), with=F]}
 		bind_cols(o)
 	}

# test
wrkSmall <- wrk %>% group_by(pop) %>% slice(1:10)
outSmall<-methodBergland(wrkSmall)

# Massive
system.time(BayPassInput<-methodBergland(wrk))

# write the text file
write.table(BayPassInput, 'BayPassInputNew.txt', col.names = FALSE, row.names = FALSE)

# make the covariate File
covs<-read.csv("~/Documents/Repos/BayPass/pH_Temp_Lat_8_ponds.csv")

predR<-scale(covs$code)
temp<-scale(covs$Temp)
pH<-scale(covs$pH)
lat<-scale(covs$LatDec)

envfile<-t(cbind(predR, lat, temp, pH))
envfileB<-matrix(envfile[1,], nrow =1)
envfileBLat<-envfile[1:2,]
envfileBLatTemp<-envfile[1:3,]
envfileBLatTempPH<-envfile[1:4,]

# as per instruction from Gauthier, make a single envfile
# use PCA on env variables
# result is a phLat variable, a temp variable, and then our predation.

prenv<-t(prcomp(covs[,c("pH", "Temp", "LatDec")], center = TRUE, scale = TRUE)$x[,1:2])
prenv<-rbind(prenv, t(scale(covs[,'code'])))
rownames(prenv)<-c("pHLat", 'Temp','Predation')
prenv
write(t(prenv), "prENVFILE", ncolumns = 8)

# make the samplesize file
ss<-rep(100,8)

# write
write.table(BayPassInput, 'ALLELEFILE', col.names = FALSE, row.names = FALSE)
write(t(envfile), "~/Documents/Repos/BayPass/ENVFILE", ncolumns = 8)
write(t(envfileB), "~/Documents/Repos/BayPass/ENVFILE_B", ncolumns = 8)
write(t(envfileBLat), "~/Documents/Repos/BayPass/ENVFILE_BLat", ncolumns = 8)
write(t(envfileBLatTemp), "~/Documents/Repos/BayPass/ENVFILE_BLatTemp", ncolumns = 8)
write(t(envfileBLatTempPH), "~/Documents/Repos/BayPass/ENVFILE_BLatTempPH", ncolumns = 8)

write(t(ss), "~/Documents/Repos/BayPass/SAMPLEFILE", ncolumns = 8)