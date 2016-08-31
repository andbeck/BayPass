library(dplyr)
library(ggplot2)
library(magrittr)
library(data.table)
library(foreach)


### functions
	makeData <- function(nSnps, nPops) {
		df <- data.table(ma=runif(nSnps*nPops))
		df[,mb:=1-ma]
		df[,code:=rep(LETTERS[1:nPops], each=nSnps)]
		
		return(df)
	}
	
	methodBeckermanOne <- function(dat) {
		i=1
		out<-bind_cols(split(dat, dat$code))
		return(out[1:(i+2)!=(i+2)])
	 }
  	
  	methodBeckermanTwo <- function(dat) {
  		out <- dat %>% 
		  group_by(code) %>% 
		  do(each_df = select(., -code)) %$%
		  bind_cols(each_df)
		  
		out
	}
 
 	methodBergland <- function(dat) {
 		setkey(dat, code)
 		
 		pops <- unique(dat$code)
 		
 		o <- foreach(i=pops)%do%{dat[J(i), c("ma", "mb"), with=F]}
 		bind_cols(o)
 	}
 	
 	
### testrun

	o <- foreach(i=10^c(1:6), .combine="rbind")%do%{
		foreach(j=1:10, .combine="rbind")%do%{
			print(paste(i, j, sep=" / "))
			testDat <- makeData(nSnps=i, nPops=j)
		
			beck1 <- system.time(methodBeckermanOne(testDat))
			beck2 <- system.time(methodBeckermanTwo(testDat))
			berg1 <- system.time(methodBergland(testDat))
			
			data.table(time=c(beck1[2], beck2[2], berg1[2]), 
						method=c("beck1", "beck2", "berg1"),
						nSNPs=i,
						nPops=j)
		
		}
	}
	
	
ggplot(data=o, aes(x=log10(nSNPs), y=time, group=nPops, color=nPops)) + 
	geom_line() + 
	facet_wrap(~method)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	













	
	AA<-data.frame(ma = runif(10))
	AA$mb<-1-AA$ma
	
	BB<-data.frame(ma = runif(10))
	BB$mb<-1-BB$ma
	
	AA
	BB
	
	df<-rbind(AA,BB)
	df$code<-rep(c("A","B"), each = 10)

# does work, but retains code columns
i=1
out<-bind_cols(split(df, df$code))
out[1:(i+2)!=(i+2)]

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # rbind(
	# c(AA[1,], BB[1,]),
	# c(AA[2,], BB[2,])
	# # etc
	# )

# YES!
bind_cols(AA,BB)

# Does not work
df %>%
	group_by(code) %>%
		do(vals = bind_cols(data.frame(.)))

