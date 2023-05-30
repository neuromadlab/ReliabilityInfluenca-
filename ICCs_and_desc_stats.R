library(ggplot2)
library(foreign)
library(MASS)
library(cowplot)
library(ggridges)
library(Hmisc)
library(car)
library(lme4)
library(viridis)
library(tidyverse)
library(broom.mixed)
library(tidybayes)

dRun <- read.csv('~Influenca_TUE002_Par_state_ID10.csv', header=T)
d10_lambda <- read.csv('~Results/TPlot10_2LR_add.csv', header=T)
fitVal10 <- read.csv('~/Results/Fitval10_2LR_add.csv', header=T)

fitVal10_original <- read.csv('~/Results/Fitval10_G.csv', header=T)


fitVal10$beta[fitVal10$beta > 10] <- 10
fitVal10_original$beta[fitVal10_original$beta > 10] <- 10

dRun$rHunger  <- dRun$Hunger / 100
dRun$rSatiety <- dRun$Satiety / 100
dRun$rThirst  <- dRun$Thirst / 100
dRun$rL_Meal  <- dRun$L_Meal / 6
dRun$rAlert   <- dRun$Alert / 100
dRun$rHappy   <- dRun$Happy / 100
dRun$rSad     <- dRun$Sad / 100
dRun$rStress  <- dRun$Stress / 100
dRun$rDist_Env  <- dRun$Dist_Env / 100
dRun$rDist_Tho  <- dRun$Dist_Tho / 100

fitVal10$logRun <- log(fitVal10$Run)
fitVal10_original$logRun <- log(fitVal10_original$Run)
dRun$logRun <- log(dRun$Run)


d10_lambda$logRT <- log(d10_lambda$RT)
d10_lambda$logRT[d10_lambda$RT == -9999] <- NA

d10_lambda$logRT[d10_lambda$RT >= 10000] <- NA
d10_lambda$logRT[d10$RT <50] <- NA
d10_lambda$logRun <- log(d10_lambda$Run)

library(lmerTest)
require(lattice)

# fit null models and conditional models to estimate ICCs
# model parameters
fm_n1.1 <- lmer(alpha_rew ~ (1|ID), fitVal10)
fm_n1.1_T <- lmer(alpha_rew ~ 1 + logRun + (1|ID), fitVal10)
fm_n1.1_b <- lmer(alpha_pun ~ (1|ID), fitVal10)
fm_n1.1_b_T <- lmer(alpha_pun ~ 1 + logRun + (1|ID), fitVal10)
#summary(fm_n1.1)
fm_n1.2 <- lmer(beta ~ (1|ID), fitVal10)
fm_n1.2_T <- lmer(beta ~ 1 + logRun + (1|ID), fitVal10)
#summary(fm_n1.2)
fm_n1.3 <- lmer(lambda ~ (1|ID), fitVal10)
fm_n1.3_T <- lmer(lambda ~ 1 + logRun + (1|ID), fitVal10)
#summary(fm_n1.3)
fm_n1.4 <- lmer(meanRew ~ (1 | ID), dRun)
fm_n1.4_T <- lmer(meanRew ~ 1 + logRun + (1 | ID), dRun)
fm_n1.5 <- lmer(fval ~ (1 | ID), fitVal10)
fm_n1.5_T <- lmer(fval ~ 1 + logRun + (1 | ID), fitVal10)
fm_RT <- lmer(logRT ~ (1|ID),d10_lambda)
fm_RT2 <- lmer(logRT ~ 1 + logRun + (1|ID),d10_lambda)

#same for original data
fm_a1.1 <- lmer(alpha ~ (1|ID), fitVal10_original)
fm_a1.1_T <- lmer(alpha ~ 1 + logRun + (1|ID), fitVal10_original)
fm_a1.2 <- lmer(beta ~ (1|ID), fitVal10_original)
fm_a1.2_T <- lmer(beta ~ 1 + logRun + (1|ID), fitVal10_original)
#summary(fm_n1.2)
fm_a1.3 <- lmer(gamma ~ (1|ID), fitVal10_original)
fm_a1.3_T <- lmer(gamma ~ 1 + logRun + (1|ID), fitVal10_original)
fm_a1.5 <- lmer(fval ~ (1 | ID), fitVal10_original)
fm_a1.5_T <- lmer(fval ~ 1 + logRun + (1 | ID), fitVal10_original)

VAR_mat_alpha_rew <- as.data.frame(VarCorr(fm_n1.1))
ICC_alpha_rew <- VAR_mat_alpha_rew$vcov[1] / sum(VAR_mat_alpha_rew$vcov)

VAR_mat_alphaT_rew <- as.data.frame(VarCorr(fm_n1.1_T))
ICC_alphaT_rew <- VAR_mat_alphaT_rew$vcov[1] / sum(VAR_mat_alphaT_rew$vcov)

VAR_mat_alpha_pun <- as.data.frame(VarCorr(fm_n1.1_b))
ICC_alpha_pun <- VAR_mat_alpha_pun$vcov[1] / sum(VAR_mat_alpha_pun$vcov)

VAR_mat_alphaT_pun <- as.data.frame(VarCorr(fm_n1.1_b_T))
ICC_alphaT_pun <- VAR_mat_alphaT_pun$vcov[1] / sum(VAR_mat_alphaT_pun$vcov)

VAR_mat_beta <- as.data.frame(VarCorr(fm_n1.2))
ICC_beta <- VAR_mat_beta$vcov[1] / sum(VAR_mat_beta$vcov)

VAR_mat_betaT <- as.data.frame(VarCorr(fm_n1.2_T))
ICC_betaT <- VAR_mat_betaT$vcov[1] / sum(VAR_mat_betaT$vcov)

VAR_mat_lambda <- as.data.frame(VarCorr(fm_n1.3))
ICC_lambda <- VAR_mat_lambda$vcov[1] / sum(VAR_mat_lambda$vcov)

VAR_mat_lambdaT <- as.data.frame(VarCorr(fm_n1.3_T))
ICC_lambdaT <- VAR_mat_lambdaT$vcov[1] / sum(VAR_mat_lambdaT$vcov)

VAR_mat_L <- as.data.frame(VarCorr(fm_n1.5))
ICC_L <- VAR_mat_L$vcov[1] / sum(VAR_mat_L$vcov)

VAR_mat_L_T <- as.data.frame(VarCorr(fm_n1.5_T))
ICC_L_T <- VAR_mat_L_T$vcov[1] / sum(VAR_mat_L_T$vcov)

VAR_mat_wins <- as.data.frame(VarCorr(fm_n1.4))
ICC_wins <- VAR_mat_wins$vcov[1] / sum(VAR_mat_wins$vcov)

VAR_mat_wins_T <- as.data.frame(VarCorr(fm_n1.4_T))
ICC_wins_T <- VAR_mat_wins_T$vcov[1] / sum(VAR_mat_wins_T$vcov)

VAR_RT <- as.data.frame(VarCorr(fm_RT))
ICC_RT <- VAR_RT$vcov[1] / sum(VAR_RT$vcov)

VAR_RT_T <- as.data.frame(VarCorr(fm_RT2))
ICC_RT_T <- VAR_RT_T$vcov[1] / sum(VAR_RT_T$vcov)


#original model
VAR_mat_alpha_orig <- as.data.frame(VarCorr(fm_a1.1))
ICC_alpha_orig <- VAR_mat_alpha_orig$vcov[1] / sum(VAR_mat_alpha_orig$vcov)

VAR_mat_alphaT_orig <- as.data.frame(VarCorr(fm_a1.1_T))
ICC_alphaT_orig <- VAR_mat_alphaT_orig$vcov[1] / sum(VAR_mat_alphaT_orig$vcov)

VAR_mat_beta_orig <- as.data.frame(VarCorr(fm_a1.2))
ICC_beta_orig <- VAR_mat_beta_orig$vcov[1] / sum(VAR_mat_beta_orig$vcov)

VAR_mat_betaT_orig <- as.data.frame(VarCorr(fm_a1.2_T))
ICC_betaT_orig <- VAR_mat_betaT_orig$vcov[1] / sum(VAR_mat_betaT_orig$vcov)

VAR_mat_gamma <- as.data.frame(VarCorr(fm_a1.3))
ICC_gamma <- VAR_mat_gamma$vcov[1] / sum(VAR_mat_gamma$vcov)

VAR_mat_gammaT <- as.data.frame(VarCorr(fm_a1.3_T))
ICC_gammaT <- VAR_mat_gammaT$vcov[1] / sum(VAR_mat_gammaT$vcov)

VAR_mat_L_orig <- as.data.frame(VarCorr(fm_a1.5))
ICC_L_orig <- VAR_mat_L_orig$vcov[1] / sum(VAR_mat_L_orig$vcov)

VAR_mat_L_T_orig <- as.data.frame(VarCorr(fm_a1.5_T))
ICC_L_T_orig <- VAR_mat_L_T_orig$vcov[1] / sum(VAR_mat_L_T_orig$vcov)


# state items (continuous)
fm_n2.1 <- lmer(Hunger ~ (1 | ID), dRun)
summary(fm_n2.1)
fm_n2.2 <- lmer(Satiety ~ (1 | ID), dRun)
summary(fm_n2.2)
fm_n2.3 <- lmer(Thirst ~ (1 | ID), dRun)
summary(fm_n2.3)
fm_n2.4 <- lmer(Alert ~ (1 | ID), dRun)
summary(fm_n2.4)
fm_n2.4a <- lmer(Alert ~ logRun + (1 +logRun | ID), dRun)
summary(fm_n2.4a)
fm_n2.5 <- lmer(Happy ~ (1 | ID), dRun)
summary(fm_n2.5)
fm_n2.6 <- lmer(Sad ~ (1 | ID), dRun)
summary(fm_n2.6)
fm_n2.7 <- lmer(Stress ~ (1 | ID), dRun)
summary(fm_n2.7)
fm_n2.8 <- lmer(Dist_Env ~ (1 | ID), dRun)
summary(fm_n2.8)
fm_n2.9 <- lmer(Dist_Tho ~ (1 | ID), dRun)
summary(fm_n2.9)


VAR_mat_Hunger <- as.data.frame(VarCorr(fm_n2.1))
ICC_Hunger <- VAR_mat_Hunger$vcov[1] / sum(VAR_mat_Hunger$vcov)

VAR_mat_Satiety <- as.data.frame(VarCorr(fm_n2.2))
ICC_Satiety <- VAR_mat_Satiety$vcov[1] / sum(VAR_mat_Satiety$vcov)


VAR_mat_Thirst <- as.data.frame(VarCorr(fm_n2.3))
ICC_Thirst <- VAR_mat_Thirst$vcov[1] / sum(VAR_mat_Thirst$vcov)

VAR_mat_Alert <- as.data.frame(VarCorr(fm_n2.4))
ICC_Alert <- VAR_mat_Alert$vcov[1] / sum(VAR_mat_Alert$vcov)

VAR_mat_Happy <- as.data.frame(VarCorr(fm_n2.5))
ICC_Happy <- VAR_mat_Happy$vcov[1] / sum(VAR_mat_Happy$vcov)

VAR_mat_Sad <- as.data.frame(VarCorr(fm_n2.6))
ICC_Sad <- VAR_mat_Sad$vcov[1] / sum(VAR_mat_Sad$vcov)

VAR_mat_Stress <- as.data.frame(VarCorr(fm_n2.7))
ICC_Stress <- VAR_mat_Stress$vcov[1] / sum(VAR_mat_Stress$vcov)

VAR_mat_Dist_Env <- as.data.frame(VarCorr(fm_n2.8))
ICC_Dist_Env <- VAR_mat_Dist_Env$vcov[1] / sum(VAR_mat_Dist_Env$vcov)

VAR_mat_Dist_Tho <- as.data.frame(VarCorr(fm_n2.9))
ICC_Dist_Tho <- VAR_mat_Dist_Tho$vcov[1] / sum(VAR_mat_Dist_Tho$vcov)

# ------------------- Descriptive Statistics
# Means
mean_RT <- mean(d10$logRT, na.rm=TRUE)
mean_win <- mean(dRun$meanRew)
mean_L <- mean(-1*fitVal10$fval)
mean_alpha_rew <- mean(fitVal10$alpha_rew)
mean_alpha_pun <- mean(fitVal10$alpha_pun)
mean_beta <- mean(fitVal10$beta)
mean_lambda <- mean(fitVal10$lambda)
mean_alert <- mean(dRun$Alert)
mean_happy <- mean(dRun$Happy)
mean_sad <- mean(dRun$Sad)
mean_stress <- mean(dRun$Stress)
mean_dEnv <- mean(dRun$Dist_Env)
mean_dTho <- mean(dRun$Dist_Tho)

mean_L_orig <- mean(-1*fitVal10_original$fval)
mean_alpha <- mean(fitVal10_original$alpha)
mean_beta_orig <- mean(fitVal10_original$beta)
mean_gamma <- mean(fitVal10_original$gamma)
# standard deviation
sd_RT <- sd(d10$logRT, na.rm=TRUE)
sd_win <- sd(dRun$meanRew)
sd_L <- sd(-1*fitVal10$fval)
sd_alpha_rew <- sd(fitVal10$alpha_rew)
sd_alpha_pun <- sd(fitVal10$alpha_pun)
sd_beta <- sd(fitVal10$beta)
sd_lambda <- sd(fitVal10$lambda)
sd_alert <- sd(dRun$Alert)
sd_happy <- sd(dRun$Happy)
sd_sad <- sd(dRun$Sad)
sd_stress <- sd(dRun$Stress)
sd_dEnv <- sd(dRun$Dist_Env)
sd_dTho <- sd(dRun$Dist_Tho)

sd_L_orig <- sd(-1*fitVal10_original$fval)
sd_alpha <- sd(fitVal10_original$alpha)
sd_beta_orig <- sd(fitVal10_original$beta)
sd_gamma <- sd(fitVal10_original$gamma)
# Median
med_RT <- median(d10$logRT, na.rm=TRUE)
med_win <- median(dRun$meanRew)
med_L <- median(-1*fitVal10$fval)
med_alpha_rew <- median(fitVal10$alpha_rew)
med_alpha_pun <- median(fitVal10$alpha_pun)
med_beta <- median(fitVal10$beta)
med_lambda <- median(fitVal10$lambda)
med_alert <- median(dRun$Alert)
med_happy <- median(dRun$Happy)
med_sad <- median(dRun$Sad)
med_stress <- median(dRun$Stress)
med_dEnv <- median(dRun$Dist_Env)
med_dTho <- median(dRun$Dist_Tho)

med_L_orig <- median(-1*fitVal10_original$fval)
med_alpha <- median(fitVal10_original$alpha)
med_beta_orig <- median(fitVal10_original$beta)
med_gamma <- median(fitVal10_original$gamma)
# 10th Percentile
p10_RT <- quantile(d10$logRT, c(.10,.90), na.rm=TRUE)
p10_win <- quantile(dRun$meanRew, c(.10,.90))
p10_L <- quantile(-1*fitVal10$fval, c(.10,.90))
p10_alpha_rew <- quantile(fitVal10$alpha_rew, c(.10,.90))
p10_alpha_pun <- quantile(fitVal10$alpha_pun, c(.10,.90))
p10_beta <- quantile(fitVal10$beta, c(.10,.90))
p10_lambda <- quantile(fitVal10$lambda, c(.10,.90))
p10_alert <- quantile(dRun$Alert, c(.10,.90))
p10_happy <- quantile(dRun$Happy, c(.10,.90))
p10_sad <- quantile(dRun$Sad, c(.10,.90))
p10_stress <- quantile(dRun$Stress, c(.10,.90))
p10_dEnv <- quantile(dRun$Dist_Env, c(.10,.90))
p10_dTho <- quantile(dRun$Dist_Tho, c(.10,.90))

p10_L_orig <- quantile(-1*fitVal10_original$fval, c(.10,.90))
p10_alpha <- quantile(fitVal10_original$alpha, c(.10,.90))
p10_beta_orig <- quantile(fitVal10_original$beta, c(.10,.90))
p10_gamma <- quantile(fitVal10_original$gamma, c(.10,.90))
# fit run effects
# model parameters
fm_n1.1 <- lmer(alpha_rew ~ 1 + logRun + (1 + logRun|ID), fitVal10)
summary(fm_n1.1)
fm_n1.1_b <- lmer(alpha_pun ~ 1 + logRun + (1 + logRun|ID), fitVal10)
summary(fm_n1.1_b)
fm_n1.2 <- lmer(beta ~ 1 + logRun + (1 + logRun|ID), fitVal10)
summary(fm_n1.2)
fm_n1.3 <- lmer(lambda ~1 + logRun + (1 + logRun|ID), fitVal10)
summary(fm_n1.3)
fm_n1.4 <- lmer(meanRew ~ 1 + logRun + (1 + logRun|ID), dRun)
summary(fm_n1.4)
fm_n1.5 <- lmer(fval ~ 1 + logRun + (1 + logRun|ID), fitVal10)
summary(fm_n1.5)
fm_RT <- lmer(logRT ~ 1 + logRun + (1 + logRun|ID),d10_lambda)
summary(fm_RT)