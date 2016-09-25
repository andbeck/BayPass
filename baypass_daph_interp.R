# setups
library(corrplot)
library(ape)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggthemes)

source('baypass_utils.R') # this is in the repo

# Covariance matrix
omega<-as.matrix(read.table("anaprEnvfile__mat_omega.out"))
omega<-matrix(omega, 8,8, byrow = TRUE, dimnames=list(c(paste("pop:",1:8)), c(paste("pop:",1:8))))

omega_plot<-cov2cor(omega)
par(mfrow = c(2,1))
corrplot(omega_plot,method = "square", order = 'hclust',mar = c(2,1,1,0))
plot(hclust(dist(omega_plot), 'complete'))

# This needs to be developed to establish the criteria for detecting XtX outlier
# see paper
#
simDat<-simulate.baypass(omega, sample.size = 10000)

beta<-read.table('anaprEnvfile__summary_betai_reg.out', header = TRUE)
XtX<-read.table('anaprEnvfile__summary_pi_xtx.out', header = TRUE)

beta_dB_XtX<-data.frame(beta = beta$Beta_is, dB = beta$BF.dB., XtX = XtX$M_XtX)
decisive_dB<-filter(beta_dB_XtX, dB>20 & XtX>quantile(beta_dB_XtX$XtX, 0.99))
strongAssoc<-filter(decisive_dB, beta>=0.2|beta<=c(-0.2))

# these take a long time - ggplot is slow and the criteria are wrong
p1<-ggplot(beta_dB_XtX, aes(x = XtX, y = dB))+
	geom_point(size = 0.5, alpha = 0.2)+
	geom_point(data = decisive_dB, aes(x = XtX, y = dB), size = 3, shape = 1, col = 'red')+
	geom_hline(yintercept = 20, col = "red", linetype = 'dashed')+
	geom_vline(xintercept = quantile(beta_dB_XtX$XtX, 0.99), col = "red", linetype = 'dashed')+
	theme_base()


p2<-ggplot(beta_dB_XtX, aes(x = XtX, y = beta))+
	geom_point(size = 2, alpha = 0.2)+
	geom_point(data = decisive_dB, aes(x = XtX, y = beta), size = 3, shape = 1,col = 'red')+
	geom_point(data = strongAssoc, aes(x = XtX, y = beta), size = 3, shape = 1,col = 'green')+
	geom_hline(yintercept = c(0.2,-0.2), col = "red", linetype = 'dashed')+
	geom_vline(xintercept = quantile(beta_dB_XtX$XtX, 0.99), col = "red", linetype = 'dashed')+
	theme_base()
	
grid.arrange(p1, p2, nrow = 2)