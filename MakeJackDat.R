# APB - Setup For Jacknife

# ----------------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(magrittr)
library(data.table)
library(foreach)

# NEW DATA
#load('~/Documents/Repos/BayPass/datUse.clean.Rdata')
load('~/GitHubRepos/BayPass/datUse.clean.Rdata')

# clean
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
wrkSmall <- data.table(wrk) %>% group_by(pop) %>% slice(1:10)
outSmall<-methodBergland(data.table(wrkSmall))

popsCode<-unique(wrk$pop)
ll<-length(popsCode)
JackList<-list()

for(i in 1:ll){
	tmp<-droplevels(filter(wrk, pop!=popsCode[i]))
	JackList[[i]]<-methodBergland(data.table(tmp))
	names(JackList)[i]<-paste("Allelfile_no",popsCode[i], sep="")
}

# write the text file
for (i in 1:ll){
	write.table(JackList[[i]], 
		file=paste("/Volumes/TTYLMF/",names(JackList)[i],sep=""),
		col.names = FALSE, row.names = FALSE)	
}

#### cov files
# make the covariate File
covs<-read.csv("~/GitHubRepos/BayPass/pH_Temp_Lat_8_ponds.csv")

predR<-scale(covs$code)
temp<-scale(covs$Temp)
pH<-scale(covs$pH)
lat<-scale(covs$LatDec)

# as per instruction from Gauthier, make a single envfile
# use PCA on env variables
# result is a phLat variable, a temp variable, and then our predation.

prenv<-t(prcomp(covs[,c("pH", "Temp", "LatDec")], center = TRUE, scale = TRUE)$x[,1:2])
prenv<-rbind(prenv, t(scale(covs[,'code'])))
rownames(prenv)<-c("pHLat", 'Temp','Predation')
prenv

JackEnvList<-list()
for(i in 1:ll){
	JackEnvList[[i]]<-prenv[,-c(i)]
}

for(i in 1:ll){
	write(t(JackEnvList[[i]]), 
	file = paste("/Volumes/TTYLMF/prENVFILE_no",popsCode[i], sep = ""), 
	ncolumns = 7)
}