#### Step 0: setwd to ttylmf ####

setwd('/Volumes/TTYLMF')

#### Step 1: use the pods analysis to identify the 99% ile of XtX ####

# Jacked
XtX_NoB1<-read.table('ana.noB1__summary_pi_xtx.out', header = TRUE)
XtX_NoBag<-read.table('ana.noBag__summary_pi_xtx.out', header = TRUE)
XtX_NoD8<-read.table('ana.noD8__summary_pi_xtx.out', header = TRUE)
XtX_NoD10<-read.table('ana.noD10__summary_pi_xtx.out', header = TRUE)
XtX_NoNM2<-read.table('ana.noNM2__summary_pi_xtx.out', header = TRUE)
XtX_NoStav<-read.table('ana.noStav__summary_pi_xtx.out', header = TRUE)
XtX_NoW1<-read.table('ana.noW1__summary_pi_xtx.out', header = TRUE)
XtX_NoW6<-read.table('ana.noW6__summary_pi_xtx.out', header = TRUE)

# Simulated
XtX_POD_NoB1<-read.table('anaPOD.noB1__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoBag<-read.table('anaPOD.noBag__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoD8<-read.table('anaPOD.noD8__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoD10<-read.table('anaPOD.noD10__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoNM2<-read.table('anaPOD.noNM2__summary_pi_xtx.out', header = TRUE) # missing
XtX_POD_NoStav<-read.table('anaPOD.noStav__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoW1<-read.table('anaPOD.noW1__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoW6<-read.table('anaPOD.noW6__summary_pi_xtx.out', header = TRUE)

# calcuate the threshold XtX using quantile = 0.99, for the top 1%
XtX_NoB1_threshold<-quantile(XtX_POD_NoB1$M_XtX, 0.99)
XtX_NoB1_Sig<-XtX_NoB1[XtX_NoB1$M_XtX>=XtX_NoB1_threshold,]

# Temp Save Time
save.image(file = "JackOutliers.RData")

# percent exceeding threshold
dim(XtX_Sig)[1]/dim(XtX_Emp)[1]

# Step 2: use the pods analysis to identify the 99% ile of betas
# The tables are indexed by covariable 1-3
# 1 = lat/pH
# 2 = temp
# 3 = predation
beta_Emp<-read.table('anaprEnvfile__summary_betai_reg.out', header = TRUE)
beta_Sim<-read.table('anaPOD__summary_betai_reg.out', header = TRUE)

# the B
one<-filter(beta_Emp, COVARIABLE == '1')
two<-filter(beta_Emp, COVARIABLE == '2')
three<-filter(beta_Emp, COVARIABLE == '3')

betaThreshold<-beta_Sim %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold = quantile(BF.dB., 0.99))

betaThreshold$COVARIABLE<-c("Lat/pH", "Temp", "Predation")
betaThreshold

# scaf positions and ID to add to dfs
scaf_pos<-read.csv('scaf_pos.csv')

# create df of bayes factors and XtX and scaf positions
df_latpH<-data.frame(one, XtX = XtX_Emp$M_XtX, scaf_pos)
df_temp<-data.frame(two, XtX = XtX_Emp$M_XtX, scaf_pos)
df_predation<-data.frame(three, XtX = XtX_Emp$M_XtX, scaf_pos)
