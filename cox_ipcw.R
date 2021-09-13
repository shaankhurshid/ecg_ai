# Script to perform traditional HF modeling using Ho models

# Dependencies
library(data.table)
library(plyr)
library(stringr)
library(rms)
library(survival)
library(ggplot2)
library(reshape2)
library(timeROC)
source(file='/functions/functions.R')

# Load wide file with clinical data
charge_set <- fread(file='charge_data.csv')

# Load ECG inferences from ECG-AI
charge_ecg_mgh <- fread(file='mgh_inferences.tsv')
charge_ecg_bwh <- fread(file='bwh_inferences.tsv')

# Logit transform
charge_ecg_mgh[,ecg_logit := log(survival_curve_af_prediction/(1-survival_curve_af_prediction))]
charge_ecg_bwh[,ecg_logit := log(survival_curve_af_prediction/(1-survival_curve_af_prediction))]

# Variable formatting
charge_ecg_mgh[,':='(start_fu_yrs = start_fu/365.25, start_fu_yrs5 = start_fu/365.25/5)]
charge_ecg_bwh[,':='(start_fu_yrs = start_fu/365.25, start_fu_yrs5 = start_fu/365.25/5)]

# Split out the train/validation and test splits
charge_ecg_mgh_train <- charge_ecg_mgh[split %in% c('train','validation')]
charge_ecg_mgh_test <- charge_ecg_mgh[split %in% c('test')]

################################# Fit models in MGH train/validation
age_sex <- coxph(Surv(af_5y.t,incd_af_5y) ~ start_fu_yrs + Dem.Gender.no_filter,data=charge_ecg_mgh_train)
ecg_only <- coxph(Surv(af_5y.t,incd_af_5y) ~ ecg_logit,data=charge_ecg_mgh_train)
charge_plus <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_startfu + ecg_logit,data=charge_ecg_mgh_train)
charge_plus_int <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_startfu + ecg_logit + ecg_logit:charge_startfu,data=charge_ecg_mgh_train)
charge_plus_int_age <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_startfu + ecg_logit + ecg_logit:start_fu_yrs5,data=charge_ecg_mgh_train)

# Standardized
charge_ecg_mgh_train[,':='(charge_std = (charge_startfu - mean(charge_startfu))/sd(charge_startfu),
                           ecg_logit_std = (ecg_logit - mean(ecg_logit))/sd(ecg_logit))]
charge_plus_int <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_std + ecg_logit_std,data=charge_ecg_mgh_train)
charge_plus_int <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_std + ecg_logit_std + charge_std:ecg_logit_std,data=charge_ecg_mgh_train)

# Export cox models for validation in other datasets as needed
# save(ecg_only,file='ecg_ai_cox.RData')
# save(charge_plus,file='charge_ecg_cox.RData')

# Define linear predictor of scores in MGH test
charge_ecg_mgh_test[,age_sex_lp := predict(age_sex,newdata=charge_ecg_mgh_test,type='lp')]
charge_ecg_mgh_test[,ecg_only_lp := predict(ecg_only,newdata=charge_ecg_mgh_test,type='lp')]
charge_ecg_mgh_test[,charge_ecg_lp := predict(charge_plus,newdata=charge_ecg_mgh_test,type='lp')]

# Define linear predictor of scores in MGH train
charge_ecg_mgh_train[,age_sex_lp := predict(age_sex,newdata=charge_ecg_mgh_train,type='lp')]
charge_ecg_mgh_train[,ecg_only_lp := predict(ecg_only,newdata=charge_ecg_mgh_train,type='lp')]
charge_ecg_mgh_train[,charge_ecg_lp := predict(charge_plus,newdata=charge_ecg_mgh_train,type='lp')]

## Concordance
cstat_age_sex <- timeROC(T=charge_ecg_mgh_test$af_5y.t, delta=charge_ecg_mgh_test$incd_af_5y,
                        marker=charge_ecg_mgh_test$age_sex_lp,cause=1,times=4.999)$AUC[2]
cstat_ecg_ai <- timeROC(T=charge_ecg_mgh_test$af_5y.t, delta=charge_ecg_mgh_test$incd_af_5y,
                  marker=charge_ecg_mgh_test$ecg_only_lp,cause=1,times=4.999)$AUC[2]
cstat_charge <- timeROC(T=charge_ecg_mgh_test$af_5y.t, delta=charge_ecg_mgh_test$incd_af_5y,
                  marker=charge_ecg_mgh_test$charge_startfu,cause=1,times=4.999)$AUC[2]
cstat_combo <- timeROC(T=charge_ecg_mgh_test$af_5y.t, delta=charge_ecg_mgh_test$incd_af_5y,
                 marker=charge_ecg_mgh_test$charge_ecg_lp,cause=1,times=4.999)$AUC[2]

# CH-AI vs. CHARGE AF
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='charge_startfu',
                         response2='charge_ecg_lp',eval.t=4.999,data=charge_ecg_mgh_test,runs=500)
z_cstat <- (cstat_combo[1] - cstat_charge[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. CHARGE-AF
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='charge_startfu',
                         response2='ecg_only_lp',eval.t=4.999,data=charge_ecg_mgh_test,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_charge[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. CH-AI
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='ecg_only_lp',
                         response2='charge_ecg_lp',eval.t=4.999,data=charge_ecg_mgh_test,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_combo[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='ecg_only_lp',
                           response2='age_sex_lp',eval.t=4.999,data=charge_ecg_mgh_test,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# CH-AI vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='charge_ecg_lp',
                           response2='age_sex_lp',eval.t=4.999,data=charge_ecg_mgh_test,runs=500)
z_cstat <- (cstat_combo[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# CHARGE vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='charge_startfu',
                           response2='age_sex_lp',eval.t=4.999,data=charge_ecg_mgh_test,runs=500)
z_cstat <- (cstat_charge[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# Create standardized variables
charge_ecg_mgh_test[,':='(age_sex_std = (age_sex_lp - mean(age_sex_lp))/sd(age_sex_lp),
                     charge_std = (charge_startfu - mean(charge_startfu))/sd(charge_startfu),
                     ecg_std = (ecg_only_lp - mean(ecg_only_lp))/sd(ecg_only_lp),
                     charge_ecg_std = (charge_ecg_lp - mean(charge_ecg_lp))/sd(charge_ecg_lp))]

# Standardized models
age_sex_std <- coxph(Surv(af_5y.t,incd_af_5y) ~ age_sex_std,data=charge_ecg_mgh_test)
charge_std <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_std,data=charge_ecg_mgh_test)
ecg_only_std <- coxph(Surv(af_5y.t,incd_af_5y) ~ ecg_std,data=charge_ecg_mgh_test)
charge_plus_std <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_ecg_std,data=charge_ecg_mgh_test)

################################# Carry to BWH
# Define linear predictor of scores in BWH
charge_ecg_bwh[,age_sex_lp := predict(age_sex,newdata=charge_ecg_bwh,type='lp')]
charge_ecg_bwh[,ecg_only_lp := predict(ecg_only,newdata=charge_ecg_bwh,type='lp')]
charge_ecg_bwh[,charge_ecg_lp := predict(charge_plus,newdata=charge_ecg_bwh,type='lp')]

# Create standardized variables in BWH
charge_ecg_bwh[,':='(age_sex_std = (age_sex_lp - mean(age_sex_lp))/sd(age_sex_lp),
                    charge_std = (charge_startfu - mean(charge_startfu))/sd(charge_startfu),
                    ecg_std = (ecg_only_lp - mean(ecg_only_lp))/sd(ecg_only_lp),
                    charge_ecg_std = (charge_ecg_lp - mean(charge_ecg_lp))/sd(charge_ecg_lp))]

## Concordance
cstat_age_sex <- timeROC(T=charge_ecg_bwh$af_5y.t, delta=charge_ecg_bwh$incd_af_5y,
                         marker=charge_ecg_bwh$age_sex_lp,cause=1,times=4.999)$AUC[2]
cstat_ecg_ai <- timeROC(T=charge_ecg_bwh$af_5y.t, delta=charge_ecg_bwh$incd_af_5y,
                        marker=charge_ecg_bwh$ecg_only_lp,cause=1,times=4.999)$AUC[2]
cstat_charge <- timeROC(T=charge_ecg_bwh$af_5y.t, delta=charge_ecg_bwh$incd_af_5y,
                        marker=charge_ecg_bwh$charge_startfu,cause=1,times=4.999)$AUC[2]
cstat_combo <- timeROC(T=charge_ecg_bwh$af_5y.t, delta=charge_ecg_bwh$incd_af_5y,
                       marker=charge_ecg_bwh$charge_ecg_lp,cause=1,times=4.999)$AUC[2]

# CH-AI vs. CHARGE AF
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='charge_startfu',
                           response2='charge_ecg_lp',eval.t=4.999,data=charge_ecg_bwh,runs=500)
z_cstat <- (cstat_combo[1] - cstat_charge[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. CHARGE-AF
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='charge_startfu',
                           response2='ecg_only_lp',eval.t=4.999,data=charge_ecg_bwh,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_charge[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. CH-AI
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='ecg_only_lp',
                           response2='charge_ecg_lp',eval.t=4.999,data=charge_ecg_bwh,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_combo[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='ecg_only_lp',
                           response2='age_sex_lp',eval.t=4.999,data=charge_ecg_bwh,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# CH-AI vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='charge_ecg_lp',
                           response2='age_sex_lp',eval.t=4.999,data=charge_ecg_bwh,runs=500)
z_cstat <- (cstat_combo[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# CHARGE vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_5y.t',status='incd_af_5y',response1='charge_startfu',
                           response2='age_sex_lp',eval.t=4.999,data=charge_ecg_bwh,runs=500)
z_cstat <- (cstat_charge[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# Standardized models
age_sex <- coxph(Surv(af_5y.t,incd_af_5y) ~ age_sex_std,data=charge_ecg_bwh)
charge <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_std,data=charge_ecg_bwh)
ecg_only <- coxph(Surv(af_5y.t,incd_af_5y) ~ ecg_std,data=charge_ecg_bwh)
charge_plus <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_ecg_std,data=charge_ecg_bwh)

################################# Density Plots
### RAW SCORES
# CHARGE-AF stratified by AF
x <- list(v1=charge_ecg_mgh_test[incd_af_5y==1]$charge_startfu,v2=charge_ecg_mgh_test[incd_af_5y==0]$charge_startfu)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(4,20,2),expand=c(0,0),limits=c(4,20)) +
  scale_y_continuous(breaks=seq(0,0.4,0.1),expand=c(0,0),limits=c(0,0.4)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='CHARGE-AF',y='Density') 
ggsave('density_charge_mgh_test.pdf',
       height=2,width=2.5,units='in',scale=4)

# ECG-AI stratified by AF
x <- list(v1=charge_ecg_mgh_test[incd_af_5y==1]$ecg_only_lp,v2=charge_ecg_mgh_test[incd_af_5y==0]$ecg_only_lp)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-3,5,1),expand=c(0,0),limits=c(-3,5)) +
  scale_y_continuous(breaks=seq(0,0.6,0.1),expand=c(0,0),limits=c(0,0.6)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='ECG-AI',y='Density') 
ggsave('density_ecg_mgh_test.pdf',
       height=2,width=2.5,units='in',scale=4)

# CH-AI stratified by AF
x <- list(v1=charge_ecg_mgh_test[incd_af_5y==1]$charge_ecg_lp,v2=charge_ecg_mgh_test[incd_af_5y==0]$charge_ecg_lp)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-454,1),expand=c(0,0),limits=c(-4,5)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1),expand=c(0,0),limits=c(0,0.5)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='CHARGE-AF',y='Density') 
ggsave('density_charge_ecg_mgh_test.pdf',
       height=2,width=2.5,units='in',scale=4)

# ECG pred stratified by AF
x <- list(v1=charge_ecg_bwh[incd_af_5y==1]$ecg_only_lp,
          v2=charge_ecg_bwh[incd_af_5y==0]$ecg_only_lp)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-5,8,1),expand=c(0,0),limits=c(-5,8)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1),expand=c(0,0),limits=c(0,0.5)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='ECG-based AF risk score',y='Density') 
ggsave('density_ecg_bwh.pdf',
       height=2,width=2.5,units='in',scale=4)

# CHARGE-AF stratified by AF
x <- list(v1=charge_ecg_bwh[incd_af_5y==1]$charge_ecg_lp,
          v2=charge_ecg_bwh[incd_af_5y==0]$charge_ecg_lp)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-3,5,1),expand=c(0,0),limits=c(-3,5)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1),expand=c(0,0),limits=c(0,0.5)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='CHARGE-ECG Score',y='Density') 
ggsave('density_charge_ecg_bwh.pdf',
       height=2,width=2.5,units='in',scale=4)

################################# ROC curves
## MGH TEST
age_sex <- timeROC(T=charge_ecg_mgh_test$af_5y.t, delta=charge_ecg_mgh_test$incd_af_5y,
                  marker=charge_ecg_mgh_test$age_sex_lp,cause=1,times=c(1,2,3,4,4.999))
ecg_ai <- timeROC(T=charge_ecg_mgh_test$af_5y.t, delta=charge_ecg_mgh_test$incd_af_5y,
                  marker=charge_ecg_mgh_test$ecg_only_lp,cause=1,times=c(1,2,3,4,4.999))
charge <- timeROC(T=charge_ecg_mgh_test$af_5y.t, delta=charge_ecg_mgh_test$incd_af_5y,
                  marker=charge_ecg_mgh_test$charge_startfu,cause=1,times=c(1,2,3,4,4.999))
combo <- timeROC(T=charge_ecg_mgh_test$af_5y.t, delta=charge_ecg_mgh_test$incd_af_5y,
                 marker=charge_ecg_mgh_test$charge_ecg_lp,cause=1,times=c(1,2,3,4,4.999))

pdf(file='roc_compare_timeroc_mgh.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))
plot.new() #?
plot(ecg_ai,4.999,add=T,col='#66c2a5',lwd=1.2)
par(new=TRUE)
plot(charge,4.999,add=T,col='#fc8d62',lwd=1.2)
par(new=TRUE)
plot(age_sex,4.999,add=T,col='darkgray',lwd=1.2)
par(new=TRUE)
plot(combo,4.999,add=T,col='#984ea3',lwd=1.2)
axis(1,at=seq(1,0,-0.2),labels=c('1.0','0.8','0.6','0.4','0.2','0.0'),cex.axis=1.6)
axis(2,at=seq(0,1,0.2),cex.axis=1.6,las=2)
title(xlab='1 - Specificity',line=2.5,cex.lab=1.8)
title(ylab='Sensitivity',line=3.2,cex.lab=1.8)
legend(0.5,0.25,legend=c('Age & Sex (0.768)','CHARGE-AF (0.802)','ECG-AI (0.823)','CH-AI (0.838)'),col=c('darkgray','#fc8d62','#66c2a5','#984ea3'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)
dev.off()

################################# ROC with CI over TIME
## Estimate standard errors
set.seed(1)
age_sex_ci <- timeroc_ci(data=charge_ecg_mgh_test,times=c(1,2,3,4,4.999),time='af_5y.t',status='incd_af_5y',marker='age_sex_lp',runs=500)
ecg_ai_ci <- timeroc_ci(data=charge_ecg_mgh_test,times=c(1,2,3,4,4.999),time='af_5y.t',status='incd_af_5y',marker='ecg_only_lp',runs=500)
charge_ci <- timeroc_ci(data=charge_ecg_mgh_test,times=c(1,2,3,4,4.999),time='af_5y.t',status='incd_af_5y',marker='charge_startfu',runs=500)
charge_ecg_ci <- timeroc_ci(data=charge_ecg_mgh_test,times=c(1,2,3,4,4.999),time='af_5y.t',status='incd_af_5y',marker='charge_ecg_lp',runs=500)

## Get point estimate with CI using standard errors
age_sex_est <- data.table(times=c(1,2,3,4,4.999),
                         AUC=age_sex$AUC,
                         lower=age_sex$AUC-1.96*apply(age_sex_ci,FUN=sd,MARGIN=2),
                         upper=age_sex$AUC+1.96*apply(age_sex_ci,FUN=sd,MARGIN=2))

ecg_ai_est <- data.table(times=c(1,2,3,4,4.999),
                         AUC=ecg_ai$AUC,
                         lower=ecg_ai$AUC-1.96*apply(ecg_ai_ci,FUN=sd,MARGIN=2),
                         upper=ecg_ai$AUC+1.96*apply(ecg_ai_ci,FUN=sd,MARGIN=2))

charge_est <- data.table(times=c(1,2,3,4,4.999),
                         AUC=charge$AUC,
                         lower=charge$AUC-1.96*apply(charge_ci,FUN=sd,MARGIN=2),
                         upper=charge$AUC+1.96*apply(charge_ci,FUN=sd,MARGIN=2))

charge_ecg_est <- data.table(times=c(1,2,3,4,4.999),
                             AUC=combo$AUC,
                             lower=combo$AUC-1.96*apply(charge_ecg_ci,FUN=sd,MARGIN=2),
                             upper=combo$AUC+1.96*apply(charge_ecg_ci,FUN=sd,MARGIN=2))

## The plot
pdf(file='auc_time_mgh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,4,1,1))
par(oma=c(1,1,1,1))

x1 <- 0.9:4.9; x2 <- 0.95:4.95; x3 <- 1.05:5.05; x4 <- 1.1:5.1

plot(x=x2,y=ecg_ai_est$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,5),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#66c2a5',cex=1.2)
par(new=TRUE)
plot(x=x3,y=charge_est$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,5),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#fc8d62',cex=1.2)
par(new=TRUE)
plot(x=x4,y=age_sex_est$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,5),ylim=c(0.6,0.9),bty='n',
     pch=19,col='darkgray',cex=1.2)
par(new=TRUE)
plot(x=x1,y=charge_ecg_est$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,5),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#984ea3',cex=1.2)
axis(2,las=2,cex.axis=1.5,pos=0.7)
axis(1,at=1:5,cex.axis=1.5)
mtext("AF Prediction Window (Years)",1,cex=1.8,line=3,at=3)
mtext("Time-dependent AUROC",2,cex=1.8,line=2)

segments(x2,ecg_ai_est$lower,x2,ecg_ai_est$upper,lwd=1,col='#66c2a5')
segments(x3,charge_est$lower,x3,charge_est$upper,lwd=1,col='#fc8d62')
segments(x4,age_sex_est$lower,x4,age_sex_est$upper,lwd=1,col='darkgray')
segments(x1,charge_ecg_est$lower,x1,charge_ecg_est$upper,lwd=1,col='#984ea3')

segments(0.7,0.5,5.5,0.5,lty=5,col='black')

legend(1,0.67,c('Age & Sex','CHARGE-AF','ECG-AI',"CH-AI"),bty='n',cex=1.5,pch=19,col=c('darkgray','#fc8d62','#66c2a5','#984ea3'))

dev.off()

## BWH TEST
age_sex_b <- timeROC(T=charge_ecg_bwh$af_5y.t, delta=charge_ecg_bwh$incd_af_5y,
                   marker=charge_ecg_bwh$age_sex_lp,cause=1,times=c(1,2,3,4,4.999))
ecg_ai_b <- timeROC(T=charge_ecg_bwh$af_5y.t, delta=charge_ecg_bwh$incd_af_5y,
                  marker=charge_ecg_bwh$ecg_only_lp,cause=1,times=c(1,2,3,4,4.999))
charge_b <- timeROC(T=charge_ecg_bwh$af_5y.t, delta=charge_ecg_bwh$incd_af_5y,
                  marker=charge_ecg_bwh$charge_startfu,cause=1,times=c(1,2,3,4,4.999))
combo_b <- timeROC(T=charge_ecg_bwh$af_5y.t, delta=charge_ecg_bwh$incd_af_5y,
                 marker=charge_ecg_bwh$charge_ecg_lp,cause=1,times=c(1,2,3,4,4.999))

pdf(file='roc_compare_timeroc_bwh.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))
plot.new() #?
plot(ecg_ai_b,4.999,add=T,col='#66c2a5',lwd=1.2)
par(new=TRUE)
plot(charge_b,4.999,add=T,col='#fc8d62',lwd=1.2)
par(new=TRUE)
plot(age_sex_b,4.999,add=T,col='darkgray',lwd=1.2)
par(new=TRUE)
plot(combo_b,4.999,add=T,col='#984ea3',lwd=1.2)
axis(1,at=seq(1,0,-0.2),labels=c('1.0','0.8','0.6','0.4','0.2','0.0'),cex.axis=1.6)
axis(2,at=seq(0,1,0.2),cex.axis=1.6,las=2)
title(xlab='1 - Specificity',line=2.5,cex.lab=1.8)
title(ylab='Sensitivity',line=3.2,cex.lab=1.8)
legend(0.5,0.25,legend=c('Age & Sex (0.730)','CHARGE-AF (0.752)','ECG-AI (0.748)','CH-AI (0.777)'),col=c('darkgray','#fc8d62','#66c2a5','#984ea3'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)
dev.off()

## Estimate standard errors
age_sex_ci_b <- timeroc_ci(data=charge_ecg_bwh,times=c(1,2,3,4,4.999),time='af_5y.t',status='incd_af_5y',marker='age_sex_lp',runs=500)
ecg_ai_ci_b <- timeroc_ci(data=charge_ecg_bwh,times=c(1,2,3,4,4.999),time='af_5y.t',status='incd_af_5y',marker='ecg_only_lp',runs=500)
charge_ci_b <- timeroc_ci(data=charge_ecg_bwh,times=c(1,2,3,4,4.999),time='af_5y.t',status='incd_af_5y',marker='charge_startfu',runs=500)
charge_ecg_ci_b <- timeroc_ci(data=charge_ecg_bwh,times=c(1,2,3,4,4.999),time='af_5y.t',status='incd_af_5y',marker='charge_ecg_lp',runs=500)

## Get point estimate with CI using standard errors
age_sex_est_b <- data.table(times=c(1,2,3,4,4.999),
                          AUC=age_sex_b$AUC,
                          lower=age_sex_b$AUC-1.96*apply(age_sex_ci_b,FUN=sd,MARGIN=2),
                          upper=age_sex_b$AUC+1.96*apply(age_sex_ci_b,FUN=sd,MARGIN=2))

ecg_ai_est_b <- data.table(times=c(1,2,3,4,4.999),
                         AUC=ecg_ai_b$AUC,
                         lower=ecg_ai_b$AUC-1.96*apply(ecg_ai_ci_b,FUN=sd,MARGIN=2),
                         upper=ecg_ai_b$AUC+1.96*apply(ecg_ai_ci_b,FUN=sd,MARGIN=2))

charge_est_b <- data.table(times=c(1,2,3,4,4.999),
                         AUC=charge_b$AUC,
                         lower=charge_b$AUC-1.96*apply(charge_ci_b,FUN=sd,MARGIN=2),
                         upper=charge_b$AUC+1.96*apply(charge_ci_b,FUN=sd,MARGIN=2))

charge_ecg_est_b <- data.table(times=c(1,2,3,4,4.999),
                             AUC=combo_b$AUC,
                             lower=combo_b$AUC-1.96*apply(charge_ecg_ci_b,FUN=sd,MARGIN=2),
                             upper=combo_b$AUC+1.96*apply(charge_ecg_ci_b,FUN=sd,MARGIN=2))

## The plot
pdf(file='auc_time_bwh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,4,1,1))
par(oma=c(1,1,1,1))

x1 <- 0.9:4.9; x2 <- 0.95:4.95; x3 <- 1.05:5.05; x4 <- 1.1:5.1

plot(x=x2,y=ecg_ai_est_b$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,5),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#66c2a5',cex=1.2)
par(new=TRUE)
plot(x=x3,y=charge_est_b$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,5),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#fc8d62',cex=1.2)
par(new=TRUE)
plot(x=x4,y=age_sex_est_b$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,5),ylim=c(0.6,0.9),bty='n',
     pch=19,col='darkgray',cex=1.2)
par(new=TRUE)
plot(x=x1,y=charge_ecg_est_b$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,5),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#984ea3',cex=1.2)
axis(2,las=2,cex.axis=1.5,pos=0.7,at=seq(0.6,0.9,0.1))
axis(1,at=1:5,cex.axis=1.5)
mtext("AF Prediction Window (Years)",1,cex=1.8,line=3,at=3)
mtext("Time-dependent AUROC",2,cex=1.8,line=2)

segments(x2,ecg_ai_est_b$lower,x2,ecg_ai_est_b$upper,lwd=1,col='#66c2a5')
segments(x3,charge_est_b$lower,x3,charge_est_b$upper,lwd=1,col='#fc8d62')
segments(x4,age_sex_est_b$lower,x4,age_sex_est_b$upper,lwd=1,col='darkgray')
segments(x1,charge_ecg_est_b$lower,x1,charge_ecg_est_b$upper,lwd=1,col='#984ea3')

legend(1,0.67,c('Age & Sex','CHARGE-AF','ECG-AI',"CH-AI"),bty='n',cex=1.5,pch=19,col=c('darkgray','#fc8d62','#66c2a5','#984ea3'))

dev.off()

################################# Assess calibration
############ PART 1: SLOPES
### MGH
af <- coxph(Surv(af_5y.t,incd_af_5y) ~ age_sex_lp,data=charge_ecg_mgh_test)
age_sex_mgh <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

af <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_startfu,data=charge_ecg_mgh_test)
charge_mgh <- c(af$coefficients[1],confint(af)[1],confint(af)[2])
  
af <- coxph(Surv(af_5y.t,incd_af_5y) ~ ecg_only_lp,data=charge_ecg_mgh_test)
ecg_mgh <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

af <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_ecg_lp,data=charge_ecg_mgh_test)
combo_mgh <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

### BWH
af <- coxph(Surv(af_5y.t,incd_af_5y) ~ age_sex_lp,data=charge_ecg_bwh)
age_sex_bwh <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

af <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_startfu,data=charge_ecg_bwh)
charge_bwh <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

af <- coxph(Surv(af_5y.t,incd_af_5y) ~ ecg_only_lp,data=charge_ecg_bwh)
ecg_bwh <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

af <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_ecg_lp,data=charge_ecg_bwh)
combo_bwh <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

###################################################################### FAMILY 1: MGH
## CHARGE
charge_ecg_mgh_test[,charge_pred5 := (1-0.9718412736^exp(charge_startfu-12.58156))*100]
## Age/sex
avg_beta <- mean(charge_ecg_mgh_train$age_sex_lp)
res <- coxph(Surv(af_5y.t,incd_af_5y) ~ age_sex_lp, data=charge_ecg_mgh_train)
km <- survfit(res, data=data.frame(x1=mean(age_sex_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(5))$surv
charge_ecg_mgh_test[,age_sex_pred5 := (1-(train_s0)^exp(age_sex_lp - avg_beta))*100]
## ECG
avg_beta <- mean(charge_ecg_mgh_train$ecg_only_lp)
res <- coxph(Surv(af_5y.t,incd_af_5y) ~ ecg_only_lp, data=charge_ecg_mgh_train)
km <- survfit(res, data=data.frame(x1=mean(ecg_only_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(5))$surv
charge_ecg_mgh_test[,ecg_pred5 := (1-(train_s0)^exp(ecg_only_lp - avg_beta))*100]
## COMBO
avg_beta <- mean(charge_ecg_mgh_train$charge_ecg_lp)
res <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_ecg_lp, data=charge_ecg_mgh_train)
km <- survfit(res, data=data.frame(x1=mean(charge_ecg_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(5))$surv
charge_ecg_mgh_test[,charge_ecg_pred5 := (1-(train_s0)^exp(charge_ecg_lp - avg_beta))*100]

# Use classifier to classify scores into quantiles (size 10)
charge_ecg_mgh_test$age_sex_decile <- classifier(risk=charge_ecg_mgh_test$age_sex_pred5,ncuts=10)
charge_ecg_mgh_test$charge_decile <- classifier(risk=charge_ecg_mgh_test$charge_pred5,ncuts=10)
charge_ecg_mgh_test$charge_cal_decile <- classifier(risk=charge_ecg_mgh_test$charge_pred5_cal,ncuts=10)
charge_ecg_mgh_test$ecg_decile <- classifier(risk=charge_ecg_mgh_test$ecg_pred5,ncuts=10)
charge_ecg_mgh_test$charge_ecg_decile <- classifier(risk=charge_ecg_mgh_test$charge_ecg_pred5,ncuts=10)

### CALCULATE OBSERVED RISK IN EACH QUANTILE
setDF(charge_ecg_mgh_test)
age_sex_obv <- survivor(data=charge_ecg_mgh_test,risk_data="age_sex_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
charge_obv <- survivor(data=charge_ecg_mgh_test,risk_data="charge_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
ecg_obv <- survivor(data=charge_ecg_mgh_test,risk_data="ecg_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
charge_ecg_obv <- survivor(data=charge_ecg_mgh_test,risk_data="charge_ecg_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
charge_cal_obv <- survivor(data=charge_ecg_mgh_test,risk_data="charge_cal_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
setDT(charge_ecg_mgh_test)

### CALCULATE AVERAGE PREDICTED RISK IN EACH QUANTILE
age_sex_pred <- charge_ecg_mgh_test[,mean(age_sex_pred5),by="age_sex_decile"][order(age_sex_decile)]
charge_pred <- charge_ecg_mgh_test[,mean(charge_pred5),by="charge_decile"][order(charge_decile)]
ecg_pred <- charge_ecg_mgh_test[,mean(ecg_pred5),by="ecg_decile"][order(ecg_decile)]
charge_ecg_pred <- charge_ecg_mgh_test[,mean(charge_ecg_pred5),by="charge_ecg_decile"][order(charge_ecg_decile)]
charge_cal_pred <- charge_ecg_mgh_test[,mean(charge_pred5_cal),by="charge_cal_decile"][order(charge_cal_decile)]

### Fit adaptive hazard model for age/sex
charge_ecg_mgh_test$cox.5yr.cll_age_sex <- log(-log(1-charge_ecg_mgh_test$age_sex_pred5/100))
calibrate.cox_age_sex <- hare(data=charge_ecg_mgh_test$af_5y.t,delta=charge_ecg_mgh_test$incd_af_5y,
                              cov=as.matrix(charge_ecg_mgh_test$cox.5yr.cll_age_sex))
predict.grid.cox_age_sex <- seq(quantile(charge_ecg_mgh_test$age_sex_pred5/100,probs=0.01),
                                quantile(charge_ecg_mgh_test$age_sex_pred5/100,probs=0.99),length=100)
predict.grid.cox.cll_age_sex <- log(-log(1-predict.grid.cox_age_sex))
predict.calibrate.cox_age_sex <- phare(5,predict.grid.cox.cll_age_sex,calibrate.cox_age_sex)

# Plots for visualization
pdf(file='cal_age_sex_mgh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- age_sex_pred$V1
y <- do.call(rbind,age_sex_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='darkgray',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_age_sex*100,predict.calibrate.cox_age_sex*100,type="l",lty=1,col="darkgray",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_age_sex <- phare(5,charge_ecg_mgh_test$cox.5yr.cll_age_sex,calibrate.cox_age_sex)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.5yr.cox <- mean(abs(charge_ecg_mgh_test$age_sex_pred5/100 - predict.calibrate.cox_age_sex))
E50.5yr.cox <- median(abs(charge_ecg_mgh_test$age_sex_pred5/100 - predict.calibrate.cox_age_sex))
E90.5yr.cox <- quantile(abs(charge_ecg_mgh_test$age_sex_pred5/100 - predict.calibrate.cox_age_sex),probs=0.9)

### Fit adaptive hazard model for CHARGE-AF
charge_ecg_mgh_test$cox.5yr.cll_charge <- log(-log(1-charge_ecg_mgh_test$charge_pred5/100))
calibrate.cox_charge <- hare(data=charge_ecg_mgh_test$af_5y.t,delta=charge_ecg_mgh_test$incd_af_5y,
                             cov=as.matrix(charge_ecg_mgh_test$cox.5yr.cll_charge))
predict.grid.cox_charge <- seq(quantile(charge_ecg_mgh_test$charge_pred5/100,probs=0.01),
                               quantile(charge_ecg_mgh_test$charge_pred5/100,probs=0.99),length=100)
predict.grid.cox.cll_charge <- log(-log(1-predict.grid.cox_charge))
predict.calibrate.cox_charge <- phare(5,predict.grid.cox.cll_charge,calibrate.cox_charge)

# Plots for visualization
pdf(file='cal_charge_mgh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- charge_pred$V1
y <- do.call(rbind,charge_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#fc8d62',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_charge*100,predict.calibrate.cox_charge*100,type="l",lty=1,col="#fc8d62",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_charge <- phare(5,charge_ecg_mgh_test$cox.5yr.cll_charge,calibrate.cox_charge)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.5yr.cox <- mean(abs(charge_ecg_mgh_test$charge_pred5/100 - predict.calibrate.cox_charge))
E50.5yr.cox <- median(abs(charge_ecg_mgh_test$charge_pred5/100 - predict.calibrate.cox_charge))
E90.5yr.cox <- quantile(abs(charge_ecg_mgh_test$charge_pred5/100 - predict.calibrate.cox_charge),probs=0.9)

### Fit adaptive hazard model for ECG score
charge_ecg_mgh_test$cox.5yr.cll_ecg <- log(-log(1-charge_ecg_mgh_test$ecg_pred5/100))
calibrate.cox_ecg <- hare(data=charge_ecg_mgh_test$af_5y.t,delta=charge_ecg_mgh_test$incd_af_5y,
                          cov=as.matrix(charge_ecg_mgh_test$cox.5yr.cll_ecg))
predict.grid.cox_ecg <- seq(quantile(charge_ecg_mgh_test$ecg_pred5/100,probs=0.01),
                            quantile(charge_ecg_mgh_test$ecg_pred5/100,probs=0.99),length=100)
predict.grid.cox.cll_ecg <- log(-log(1-predict.grid.cox_ecg))
predict.calibrate.cox_ecg <- phare(5,predict.grid.cox.cll_ecg,calibrate.cox_ecg)

# Plots for visualization
pdf(file='cal_ecg_mgh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- ecg_pred$V1
y <- do.call(rbind,ecg_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#66c2a5',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_ecg*100,predict.calibrate.cox_ecg*100,type="l",lty=1,col="#66c2a5",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_ecg <- phare(5,charge_ecg_mgh_test$cox.5yr.cll_ecg,calibrate.cox_ecg)
ICI.5yr.cox <- mean(abs(charge_ecg_mgh_test$ecg_pred5/100 - predict.calibrate.cox_ecg))
E50.5yr.cox <- median(abs(charge_ecg_mgh_test$ecg_pred5/100 - predict.calibrate.cox_ecg))
E90.5yr.cox <- quantile(abs(charge_ecg_mgh_test$ecg_pred5/100 - predict.calibrate.cox_ecg),probs=0.9)

### Fit adaptive hazard model for age/sex
charge_ecg_mgh_test$cox.5yr.cll_age_sex <- log(-log(1-charge_ecg_mgh_test$age_sex_pred5/100))
calibrate.cox_age_sex <- hare(data=charge_ecg_mgh_test$af_5y.t,delta=charge_ecg_mgh_test$incd_af_5y,
                              cov=as.matrix(charge_ecg_mgh_test$cox.5yr.cll_age_sex))
predict.grid.cox_age_sex <- seq(quantile(charge_ecg_mgh_test$age_sex_pred5/100,probs=0.01),
                                quantile(charge_ecg_mgh_test$age_sex_pred5/100,probs=0.99),length=100)
predict.grid.cox.cll_age_sex <- log(-log(1-predict.grid.cox_age_sex))
predict.calibrate.cox_age_sex <- phare(5,predict.grid.cox.cll_age_sex,calibrate.cox_age_sex)

# Plots for visualization
pdf(file='cal_age_sex_mgh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- age_sex_pred$V1
y <- do.call(rbind,age_sex_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='darkgray',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_age_sex*100,predict.calibrate.cox_age_sex*100,type="l",lty=1,col="darkgray",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_age_sex <- phare(5,charge_ecg_mgh_test$cox.5yr.cll_age_sex,calibrate.cox_age_sex)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.5yr.cox <- mean(abs(charge_ecg_mgh_test$age_sex_pred5/100 - predict.calibrate.cox_age_sex))
E50.5yr.cox <- median(abs(charge_ecg_mgh_test$age_sex_pred5/100 - predict.calibrate.cox_age_sex))
E90.5yr.cox <- quantile(abs(charge_ecg_mgh_test$age_sex_pred5/100 - predict.calibrate.cox_age_sex),probs=0.9)

###################################################################### FAMILY 2: BWH
## CHARGE
charge_ecg_bwh[,charge_pred5 := (1-0.9718412736^exp(charge_startfu-12.58156))*100]
## Age/sex
avg_beta <- mean(charge_ecg_mgh_train$age_sex_lp)
res <- coxph(Surv(af_5y.t,incd_af_5y) ~ age_sex_lp, data=charge_ecg_mgh_train)
km <- survfit(res, data=data.frame(x1=mean(age_sex_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(5))$surv
charge_ecg_bwh[,age_sex_pred5 := (1-(train_s0)^exp(age_sex_lp - avg_beta))*100]
## ECG
avg_beta <- mean(charge_ecg_mgh_train$ecg_only_lp)
res <- coxph(Surv(af_5y.t,incd_af_5y) ~ ecg_only_lp, data=charge_ecg_mgh_train)
km <- survfit(res, data=data.frame(x1=mean(ecg_only_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(5))$surv
charge_ecg_bwh[,ecg_pred5 := (1-(train_s0)^exp(ecg_only_lp - avg_beta))*100]
## COMBO
avg_beta <- mean(charge_ecg_mgh_train$charge_ecg_lp)
res <- coxph(Surv(af_5y.t,incd_af_5y) ~ charge_ecg_lp, data=charge_ecg_mgh_train)
km <- survfit(res, data=data.frame(x1=mean(charge_ecg_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(5))$surv
charge_ecg_bwh[,charge_ecg_pred5 := (1-(train_s0)^exp(charge_ecg_lp - avg_beta))*100]

# Use classifier to classify scores into quantiles (size 10)
charge_ecg_bwh$age_sex_decile <- classifier(risk=charge_ecg_bwh$age_sex_pred5,ncuts=10)
charge_ecg_bwh$charge_decile <- classifier(risk=charge_ecg_bwh$charge_pred5,ncuts=10)
charge_ecg_bwh$ecg_decile <- classifier(risk=charge_ecg_bwh$ecg_pred5,ncuts=10)
charge_ecg_bwh$charge_ecg_decile <- classifier(risk=charge_ecg_bwh$charge_ecg_pred5,ncuts=10)

### CALCULATE OBSERVED RISK IN EACH QUANTILE
setDF(charge_ecg_bwh)
age_sex_obv <- survivor(data=charge_ecg_bwh,risk_data="age_sex_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
charge_obv <- survivor(data=charge_ecg_bwh,risk_data="charge_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
ecg_obv <- survivor(data=charge_ecg_bwh,risk_data="ecg_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
charge_ecg_obv <- survivor(data=charge_ecg_bwh,risk_data="charge_ecg_decile",time='af_5y.t',status='incd_af_5y',eval.t=5)
setDT(charge_ecg_bwh)

### CALCULATE AVERAGE PREDICTED RISK IN EACH QUANTILE
age_sex_pred <- charge_ecg_bwh[,mean(age_sex_pred5),by="age_sex_decile"][order(age_sex_decile)]
charge_pred <- charge_ecg_bwh[,mean(charge_pred5),by="charge_decile"][order(charge_decile)]
ecg_pred <- charge_ecg_bwh[,mean(ecg_pred5),by="ecg_decile"][order(ecg_decile)]
charge_ecg_pred <- charge_ecg_bwh[,mean(charge_ecg_pred5),by="charge_ecg_decile"][order(charge_ecg_decile)]

### Fit adaptive hazard model for age/sex
charge_ecg_bwh$cox.5yr.cll_age_sex <- log(-log(1-charge_ecg_bwh$age_sex_pred5/100))
calibrate.cox_age_sex <- hare(data=charge_ecg_bwh$af_5y.t,delta=charge_ecg_bwh$incd_af_5y,
                              cov=as.matrix(charge_ecg_bwh$cox.5yr.cll_age_sex))
predict.grid.cox_age_sex <- seq(quantile(charge_ecg_bwh$age_sex_pred5/100,probs=0.01),
                                quantile(charge_ecg_bwh$age_sex_pred5/100,probs=0.99),length=100)
predict.grid.cox.cll_age_sex <- log(-log(1-predict.grid.cox_age_sex))
predict.calibrate.cox_age_sex <- phare(5,predict.grid.cox.cll_age_sex,calibrate.cox_age_sex)

# Plots for visualization
pdf(file='cal_age_sex_bwh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- age_sex_pred$V1
y <- do.call(rbind,age_sex_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='darkgray',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_age_sex*100,predict.calibrate.cox_age_sex*100,type="l",lty=1,col="darkgray",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_age_sex <- phare(5,charge_ecg_bwh$cox.5yr.cll_age_sex,calibrate.cox_age_sex)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.5yr.cox <- mean(abs(charge_ecg_bwh$age_sex_pred5/100 - predict.calibrate.cox_age_sex))
E50.5yr.cox <- median(abs(charge_ecg_bwh$age_sex_pred5/100 - predict.calibrate.cox_age_sex))
E90.5yr.cox <- quantile(abs(charge_ecg_bwh$age_sex_pred5/100 - predict.calibrate.cox_age_sex),probs=0.9)

### Fit adaptive hazard model for CHARGE-AF
charge_ecg_bwh$cox.5yr.cll_charge <- log(-log(1-charge_ecg_bwh$charge_pred5/100))
calibrate.cox_charge <- hare(data=charge_ecg_bwh$af_5y.t,delta=charge_ecg_bwh$incd_af_5y,
                             cov=as.matrix(charge_ecg_bwh$cox.5yr.cll_charge))
predict.grid.cox_charge <- seq(quantile(charge_ecg_bwh$charge_pred5/100,probs=0.01),
                               quantile(charge_ecg_bwh$charge_pred5/100,probs=0.99),length=100)
predict.grid.cox.cll_charge <- log(-log(1-predict.grid.cox_charge))
predict.calibrate.cox_charge <- phare(5,predict.grid.cox.cll_charge,calibrate.cox_charge)

# Plots for visualization
pdf(file='cal_charge_bwh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- charge_pred$V1
y <- do.call(rbind,charge_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#fc8d62',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_charge*100,predict.calibrate.cox_charge*100,type="l",lty=1,col="#fc8d62",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_charge <- phare(5,charge_ecg_bwh$cox.5yr.cll_charge,calibrate.cox_charge)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.5yr.cox <- mean(abs(charge_ecg_bwh$charge_pred5/100 - predict.calibrate.cox_charge))
E50.5yr.cox <- median(abs(charge_ecg_bwh$charge_pred5/100 - predict.calibrate.cox_charge))
E90.5yr.cox <- quantile(abs(charge_ecg_bwh$charge_pred5/100 - predict.calibrate.cox_charge),probs=0.9)

### Fit adaptive hazard model for ECG score
charge_ecg_bwh$cox.5yr.cll_ecg <- log(-log(1-charge_ecg_bwh$ecg_pred5/100))
calibrate.cox_ecg <- hare(data=charge_ecg_bwh$af_5y.t,delta=charge_ecg_bwh$incd_af_5y,
                          cov=as.matrix(charge_ecg_bwh$cox.5yr.cll_ecg))
predict.grid.cox_ecg <- seq(quantile(charge_ecg_bwh$ecg_pred5/100,probs=0.01),
                            quantile(charge_ecg_bwh$ecg_pred5/100,probs=0.99),length=100)
predict.grid.cox.cll_ecg <- log(-log(1-predict.grid.cox_ecg))
predict.calibrate.cox_ecg <- phare(5,predict.grid.cox.cll_ecg,calibrate.cox_ecg)

# Plots for visualization
pdf(file='cal_ecg_bwh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- ecg_pred$V1
y <- do.call(rbind,ecg_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#66c2a5',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_ecg*100,predict.calibrate.cox_ecg*100,type="l",lty=1,col="#66c2a5",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_ecg <- phare(5,charge_ecg_bwh$cox.5yr.cll_ecg,calibrate.cox_ecg)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.5yr.cox <- mean(abs(charge_ecg_bwh$ecg_pred5/100 - predict.calibrate.cox_ecg))
E50.5yr.cox <- median(abs(charge_ecg_bwh$ecg_pred5/100 - predict.calibrate.cox_ecg))
E90.5yr.cox <- quantile(abs(charge_ecg_bwh$ecg_pred5/100 - predict.calibrate.cox_ecg),probs=0.9)

### Fit adaptive hazard model for CHARGE ECG score
charge_ecg_bwh$cox.5yr.cll_charge_ecg <- log(-log(1-charge_ecg_bwh$charge_ecg_pred5/100))
calibrate.cox_charge_ecg <- hare(data=charge_ecg_bwh$af_5y.t,delta=charge_ecg_bwh$incd_af_5y,
                                 cov=as.matrix(charge_ecg_bwh$cox.5yr.cll_charge_ecg))
predict.grid.cox_charge_ecg <- seq(quantile(charge_ecg_bwh$charge_ecg_pred5/100,probs=0.01),
                                   quantile(charge_ecg_bwh$charge_ecg_pred5/100,probs=0.99),length=100)
predict.grid.cox.cll_charge_ecg <- log(-log(1-predict.grid.cox_charge_ecg))
predict.calibrate.cox_charge_ecg <- phare(5,predict.grid.cox.cll_charge_ecg,calibrate.cox_charge_ecg)

# Plots for visualization
pdf(file='cal_charge_ecg_bwh.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- charge_ecg_pred$V1
y <- do.call(rbind,charge_ecg_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#8da0cb',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_charge_ecg*100,predict.calibrate.cox_charge_ecg*100,type="l",lty=1,col="#8da0cb",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_charge_ecg <- phare(5,charge_ecg_bwh$cox.5yr.cll_charge_ecg,calibrate.cox_charge_ecg)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.5yr.cox <- mean(abs(charge_ecg_bwh$charge_ecg_pred5/100 - predict.calibrate.cox_charge_ecg))
E50.5yr.cox <- median(abs(charge_ecg_bwh$charge_ecg_pred5/100 - predict.calibrate.cox_charge_ecg))
E90.5yr.cox <- quantile(abs(charge_ecg_bwh$charge_ecg_pred5/100 - predict.calibrate.cox_charge_ecg),probs=0.9)

# Check correlation between linear predictors of CHARGE-AF and ECG-AI
cor.test(charge_ecg_mgh_test$ecg_logit,charge_ecg_mgh_test$charge_startfu)
cor.test(charge_ecg_bwh$ecg_logit,charge_ecg_bwh$charge_startfu)

# Save out processed datasets if desired
#write.csv(charge_ecg_mgh_test,file='charge_ecg_mgh_test_060821.csv',row.names=F)
#write.csv(charge_ecg_bwh,file='charge_ecg_bwh_060821.csv',row.names=F)

