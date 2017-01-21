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
XtX_POD_NoNM2<-read.table('anaPOD.noNM2__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoStav<-read.table('anaPOD.noStav__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoW1<-read.table('anaPOD.noW1__summary_pi_xtx.out', header = TRUE)
XtX_POD_NoW6<-read.table('anaPOD.noW6__summary_pi_xtx.out', header = TRUE)

# calcuate the threshold XtX using quantile = 0.99, for the top 1%
XtX_NoB1_threshold<-quantile(XtX_POD_NoB1$M_XtX, 0.99)
XtX_NoB1_Sig<-XtX_NoB1[XtX_NoB1$M_XtX>=XtX_NoB1_threshold,]

XtX_NoBag_threshold<-quantile(XtX_POD_NoBag$M_XtX, 0.99)
XtX_NoBag_Sig<-XtX_NoBag[XtX_NoBag$M_XtX>=XtX_NoBag_threshold,]

XtX_NoD8_threshold<-quantile(XtX_POD_NoD8$M_XtX, 0.99)
XtX_NoD8_Sig<-XtX_NoD8[XtX_NoD8$M_XtX>=XtX_NoD8_threshold,]

XtX_NoD10_threshold<-quantile(XtX_POD_NoD10$M_XtX, 0.99)
XtX_NoD10_Sig<-XtX_NoD10[XtX_NoD10$M_XtX>=XtX_NoD10_threshold,]

XtX_NoNM2_threshold<-quantile(XtX_POD_NoNM2$M_XtX, 0.99)
XtX_NoNM2_Sig<-XtX_NoNM2[XtX_NoNM2$M_XtX>=XtX_NoNM2_threshold,]

XtX_NoStav_threshold<-quantile(XtX_POD_NoStav$M_XtX, 0.99)
XtX_NoStav_Sig<-XtX_NoStav[XtX_NoStav$M_XtX>=XtX_NoStav_threshold,]

XtX_NoW1_threshold<-quantile(XtX_POD_NoW1$M_XtX, 0.99)
XtX_NoW1_Sig<-XtX_NoW1[XtX_NoW1$M_XtX>=XtX_NoW1_threshold,]

XtX_NoW6_threshold<-quantile(XtX_POD_NoW6$M_XtX, 0.99)
XtX_NoW6_Sig<-XtX_NoW6[XtX_NoW6$M_XtX>=XtX_NoW6_threshold,]

#### Step 2: use the pods analysis to identify the 99% ile of betais ####

# The tables are indexed by covariable 1-3
# 1 = lat/pH
# 2 = temp
# 3 = predation

# Empirical betai
beta_Emp_NoB1<-read.table('ana.noB1__summary_betai_reg.out', header = TRUE)
beta_Emp_NoBag<-read.table('ana.noBag__summary_betai_reg.out', header = TRUE)
beta_Emp_NoD8<-read.table('ana.noD8__summary_betai_reg.out', header = TRUE)
beta_Emp_NoD10<-read.table('ana.noD10__summary_betai_reg.out', header = TRUE)
beta_Emp_NoNM2<-read.table('ana.noNM2__summary_betai_reg.out', header = TRUE)
beta_Emp_NoStav<-read.table('ana.noStav__summary_betai_reg.out', header = TRUE)
beta_Emp_NoW1<-read.table('ana.noW1__summary_betai_reg.out', header = TRUE)
beta_Emp_NoW6<-read.table('ana.noW6__summary_betai_reg.out', header = TRUE)

# Simulate (randomise) betai
beta_Sim_NoB1<-read.table('anaPOD.noB1__summary_betai_reg.out', header = TRUE)
beta_Sim_NoBag<-read.table('anaPOD.noBag__summary_betai_reg.out', header = TRUE)
beta_Sim_NoD8<-read.table('anaPOD.noD8__summary_betai_reg.out', header = TRUE)
beta_Sim_NoD10<-read.table('anaPOD.noD10__summary_betai_reg.out', header = TRUE)
beta_Sim_NoNM2<-read.table('anaPOD.noNM2__summary_betai_reg.out', header = TRUE)
beta_Sim_NoStav<-read.table('anaPOD.noStav__summary_betai_reg.out', header = TRUE)
beta_Sim_NoW1<-read.table('anaPOD.noW1__summary_betai_reg.out', header = TRUE)
beta_Sim_NoW6<-read.table('anaPOD.noW6__summary_betai_reg.out', header = TRUE)

# empirical betai sorted to Covariates

# Lat/pH
one_NoB1<-filter(beta_Emp_NoB1, COVARIABLE == '1')
one_NoBag<-filter(beta_Emp_NoBag, COVARIABLE == '1')
one_NoD8<-filter(beta_Emp_NoD8, COVARIABLE == '1')
one_NoD10<-filter(beta_Emp_NoD10, COVARIABLE == '1')
one_NoNM2<-filter(beta_Emp_NoNM2, COVARIABLE == '1')
one_NoStav<-filter(beta_Emp_NoStav, COVARIABLE == '1')
one_NoW1<-filter(beta_Emp_NoW1, COVARIABLE == '1')
one_NoW6<-filter(beta_Emp_NoW6, COVARIABLE == '1')

# temp
two_NoB1<-filter(beta_Emp_NoB1, COVARIABLE == '2')
two_NoBag<-filter(beta_Emp_NoBag, COVARIABLE == '2')
two_NoD8<-filter(beta_Emp_NoD8, COVARIABLE == '2')
two_NoD10<-filter(beta_Emp_NoD10, COVARIABLE == '2')
two_NoNM2<-filter(beta_Emp_NoNM2, COVARIABLE == '2')
two_NoStav<-filter(beta_Emp_NoStav, COVARIABLE == '2')
two_NoW1<-filter(beta_Emp_NoW1, COVARIABLE == '2')
two_NoW6<-filter(beta_Emp_NoW6, COVARIABLE == '2')

# predation
three_NoB1<-filter(beta_Emp_NoB1, COVARIABLE == '3')
three_NoBag<-filter(beta_Emp_NoBag, COVARIABLE == '3')
three_NoD8<-filter(beta_Emp_NoD8, COVARIABLE == '3')
three_NoD10<-filter(beta_Emp_NoD10, COVARIABLE == '3')
three_NoNM2<-filter(beta_Emp_NoNM2, COVARIABLE == '3')
three_NoStav<-filter(beta_Emp_NoStav, COVARIABLE == '3')
three_NoW1<-filter(beta_Emp_NoW1, COVARIABLE == '3')
three_NoW6<-filter(beta_Emp_NoW6, COVARIABLE == '3')

#### Step 3: Note the scaf positions ####
# only need to do this once
scaf_pos<-read.csv('~/GithubRepos/BayPass/scaf_pos.csv')

#### Step 4: Process betai ####

# For each Jack - this creates the betai thresholds for each covariable
betaThreshold_NoB1<-beta_Sim_NoB1 %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold_NoB1 = quantile(BF.dB., 0.99))

betaThreshold_NoBag<-beta_Sim_NoBag %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold_NoBag = quantile(BF.dB., 0.99))

betaThreshold_NoD8<-beta_Sim_NoD8 %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold_NoD8 = quantile(BF.dB., 0.99))

betaThreshold_NoD10<-beta_Sim_NoD10 %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold_NoD10 = quantile(BF.dB., 0.99))

betaThreshold_NoNM2<-beta_Sim_NoNM2 %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold_NoNM2 = quantile(BF.dB., 0.99))

betaThreshold_NoStav<-beta_Sim_NoStav %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold_NoStav = quantile(BF.dB., 0.99))

betaThreshold_NoW1<-beta_Sim_NoW1 %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold_NoW1 = quantile(BF.dB., 0.99))

betaThreshold_NoW6<-beta_Sim_NoW6 %>%
  group_by(COVARIABLE) %>%
  summarise(
    betaThreshold_NoW6 = quantile(BF.dB., 0.99))

#### Step 6: Collect betai and XtX and scaf positions####

# create df of bayes factors and XtX and scaf positions
df_latpH_NoB1<-data.frame(one_NoB1, XtX_NoB1 = XtX_NoB1$M_XtX, scaf_pos)
df_latpH_NoBag<-data.frame(one_NoBag, XtX_NoBag = XtX_NoBag$M_XtX, scaf_pos)
df_latpH_NoD8<-data.frame(one_NoD8, XtX_NoD8 = XtX_NoD8$M_XtX, scaf_pos)
df_latpH_NoD10<-data.frame(one_NoD10, XtX_NoD10 = XtX_NoD10$M_XtX, scaf_pos)
df_latpH_NoNM2<-data.frame(one_NoNM2, XtX_NoNM2 = XtX_NoNM2$M_XtX, scaf_pos)
df_latpH_NoStav<-data.frame(one_NoStav, XtX_NoStav = XtX_NoStav$M_XtX, scaf_pos)
df_latpH_NoW1<-data.frame(one_NoW1, XtX_NoW1 = XtX_NoW6$M_XtX, scaf_pos)
df_latpH_NoW6<-data.frame(one_NoW6, XtX_NoW6 = XtX_NoW6$M_XtX, scaf_pos)

df_temp_NoB1<-data.frame(two_NoB1, XtX_NoB1 = XtX_NoB1$M_XtX, scaf_pos)
df_temp_NoBag<-data.frame(two_NoBag, XtX_NoBag = XtX_NoBag$M_XtX, scaf_pos)
df_temp_NoD8<-data.frame(two_NoD8, XtX_NoD8 = XtX_NoD8$M_XtX, scaf_pos)
df_temp_NoD10<-data.frame(two_NoD10, XtX_NoD10 = XtX_NoD10$M_XtX, scaf_pos)
df_temp_NoNM2<-data.frame(two_NoNM2, XtX_NoNM2 = XtX_NoNM2$M_XtX, scaf_pos)
df_temp_NoStav<-data.frame(two_NoStav, XtX_NoStav = XtX_NoStav$M_XtX, scaf_pos)
df_temp_NoW1<-data.frame(two_NoW1, XtX_NoW1 = XtX_NoW6$M_XtX, scaf_pos)
df_temp_NoW6<-data.frame(two_NoW6, XtX_NoW6 = XtX_NoW6$M_XtX, scaf_pos)

df_predation_NoB1<-data.frame(three_NoB1, XtX_NoB1 = XtX_NoB1$M_XtX, scaf_pos)
df_predation_NoBag<-data.frame(three_NoBag, XtX_NoBag = XtX_NoBag$M_XtX, scaf_pos)
df_predation_NoD8<-data.frame(three_NoD8, XtX_NoD8 = XtX_NoD8$M_XtX, scaf_pos)
df_predation_NoD10<-data.frame(three_NoD10, XtX_NoD10 = XtX_NoD10$M_XtX, scaf_pos)
df_predation_NoNM2<-data.frame(three_NoNM2, XtX_NoNM2 = XtX_NoNM2$M_XtX, scaf_pos)
df_predation_NoStav<-data.frame(three_NoStav, XtX_NoStav = XtX_NoStav$M_XtX, scaf_pos)
df_predation_NoW1<-data.frame(three_NoW1, XtX_NoW1 = XtX_NoW6$M_XtX, scaf_pos)
df_predation_NoW6<-data.frame(three_NoW6, XtX_NoW6 = XtX_NoW6$M_XtX, scaf_pos)

#### Step 5: ID the Ouliers ####

# these are the outliers based on XtX and bayes factor for each 
latpH_keyptsID_NoB1<-filter(df_latpH_NoB1, BF.dB. > as.numeric(betaThreshold_NoB1[1,2]) & XtX_NoB1 > XtX_NoB1_threshold) # 2890
latpH_keyptsID_NoBag<-filter(df_latpH_NoBag, BF.dB. > as.numeric(betaThreshold_NoBag[1,2]) & XtX_NoBag > XtX_NoBag_threshold) # 2890
latpH_keyptsID_NoD8<-filter(df_latpH_NoD8, BF.dB. > as.numeric(betaThreshold_NoD8[1,2]) & XtX_NoD8 > XtX_NoD8_threshold) # 2890
latpH_keyptsID_NoD10<-filter(df_latpH_NoD10, BF.dB. > as.numeric(betaThreshold_NoD10[1,2]) & XtX_NoD10 > XtX_NoD10_threshold) # 2890
latpH_keyptsID_NoNM2<-filter(df_latpH_NoNM2, BF.dB. > as.numeric(betaThreshold_NoNM2[1,2]) & XtX_NoNM2 > XtX_NoNM2_threshold) # 2890
latpH_keyptsID_NoStav<-filter(df_latpH_NoStav, BF.dB. > as.numeric(betaThreshold_NoStav[1,2]) & XtX_NoStav > XtX_NoStav_threshold) # 2890
latpH_keyptsID_NoW1<-filter(df_latpH_NoW1, BF.dB. > as.numeric(betaThreshold_NoW1[1,2]) & XtX_NoW1 > XtX_NoW1_threshold) # 2890
latpH_keyptsID_NoW6<-filter(df_latpH_NoW6, BF.dB. > as.numeric(betaThreshold_NoW6[1,2]) & XtX_NoW6 > XtX_NoW6_threshold) # 2890

temp_keyptsID_NoB1<-filter(df_temp_NoB1, BF.dB. > as.numeric(betaThreshold_NoB1[2,2]) & XtX_NoB1 > XtX_NoB1_threshold) # 2890
temp_keyptsID_NoBag<-filter(df_temp_NoBag, BF.dB. > as.numeric(betaThreshold_NoBag[2,2]) & XtX_NoBag > XtX_NoBag_threshold) # 2890
temp_keyptsID_NoD8<-filter(df_temp_NoD8, BF.dB. > as.numeric(betaThreshold_NoD8[2,2]) & XtX_NoD8 > XtX_NoD8_threshold) # 2890
temp_keyptsID_NoD10<-filter(df_temp_NoD10, BF.dB. > as.numeric(betaThreshold_NoD10[2,2]) & XtX_NoD10 > XtX_NoD10_threshold) # 2890
temp_keyptsID_NoNM2<-filter(df_temp_NoNM2, BF.dB. > as.numeric(betaThreshold_NoNM2[2,2]) & XtX_NoNM2 > XtX_NoNM2_threshold) # 2890
temp_keyptsID_NoStav<-filter(df_temp_NoStav, BF.dB. > as.numeric(betaThreshold_NoStav[2,2]) & XtX_NoStav > XtX_NoStav_threshold) # 2890
temp_keyptsID_NoW1<-filter(df_temp_NoW1, BF.dB. > as.numeric(betaThreshold_NoW1[2,2]) & XtX_NoW1 > XtX_NoW1_threshold) # 2890
temp_keyptsID_NoW6<-filter(df_temp_NoW6, BF.dB. > as.numeric(betaThreshold_NoW6[2,2]) & XtX_NoW6 > XtX_NoW6_threshold) # 2890

predation_keyptsID_NoB1<-filter(df_predation_NoB1, BF.dB. > as.numeric(betaThreshold_NoB1[3,2]) & XtX_NoB1 > XtX_NoB1_threshold) # 2890
predation_keyptsID_NoBag<-filter(df_predation_NoBag, BF.dB. > as.numeric(betaThreshold_NoBag[3,2]) & XtX_NoBag > XtX_NoBag_threshold) # 2890
predation_keyptsID_NoD8<-filter(df_predation_NoD8, BF.dB. > as.numeric(betaThreshold_NoD8[3,2]) & XtX_NoD8 > XtX_NoD8_threshold) # 2890
predation_keyptsID_NoD10<-filter(df_predation_NoD10, BF.dB. > as.numeric(betaThreshold_NoD10[3,2]) & XtX_NoD10 > XtX_NoD10_threshold) # 2890
predation_keyptsID_NoNM2<-filter(df_predation_NoNM2, BF.dB. > as.numeric(betaThreshold_NoNM2[3,2]) & XtX_NoNM2 > XtX_NoNM2_threshold) # 2890
predation_keyptsID_NoStav<-filter(df_predation_NoStav, BF.dB. > as.numeric(betaThreshold_NoStav[3,2]) & XtX_NoStav > XtX_NoStav_threshold) # 2890
predation_keyptsID_NoW1<-filter(df_predation_NoW1, BF.dB. > as.numeric(betaThreshold_NoW1[3,2]) & XtX_NoW1 > XtX_NoW1_threshold) # 2890
predation_keyptsID_NoW6<-filter(df_predation_NoW6, BF.dB. > as.numeric(betaThreshold_NoW6[3,2]) & XtX_NoW6 > XtX_NoW6_threshold) # 2890

# Save Time
save.image(file = "JackOutliers.RData")
setwd('/Volumes/TTYLMF/')
load("JackOutliers.RData")
