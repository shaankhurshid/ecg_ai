# Depends
library(data.table)
library(survival)
library(polspline)
library(prodlim)
library(Cairo)
library(timeROC)
source(file='/functions/functions.R')

# Load ECG inferences from ECG-AI with UKBB clinical data
ecg <- fread(file='ukbb_inferences.tsv')

##### Step 2: Calculate CHARGE-AF score
# Age at ECG
ecg[,age_at_ecg := enroll_age + as.numeric(ecg_date - enroll_date)/365.25]

# Classify White race
ecg[,white := ifelse(self_race %in% c(1,1001,1002,1003),1,0)]

# Variable cleanup
## Remove people with missing continuous variables
ecg <- ecg[!is.na(ht) & !is.na(wt) & !is.na(sbp) & !is.na(dbp)] #N=74

## Assume zero for people without smoking data or BP med data
ecg[,prev_bpmed := ifelse(!is.na(prev_bpmed) & prev_bpmed==1,1,0)]
ecg[,smoker := ifelse(!is.na(smoker) & smoker==1,1,0)]

# Calculate necessary CHARGE-AF variables
ecg[,':='(age5 = (age_at_ecg/5),
          ht_10 = ht/10, wt_15 = wt/15,
          sbp20 = sbp/20, dbp10 = dbp/10)]

# Prevalent disease variables
ecg[,':='(prev_mi = ifelse(c(has_mi==1 & (mi_date <= ecg_date)),1,0),
          prev_dm = ifelse(c(has_dm==1 & (dm_date <= ecg_date)),1,0),
          prev_hf = ifelse(c(has_hf==1 & (hf_date <= ecg_date)),1,0))]

# AF variables
ecg[,':='(prev_af = ifelse(c(has_af==1 & (af_date <= ecg_date)),1,0),
          incd_af = ifelse(c(has_af==1 & (af_date > ecg_date)),1,0))]
ecg[,':='(incd_af_2y = ifelse(incd_af==1 & c((as.numeric(af_date) - as.numeric(ecg_date)) <= 365.25*2),1,0))]
ecg[,':='(af_2y.t = ifelse(incd_af_2y==1,(as.numeric(af_date) - as.numeric(ecg_date))/365.25,
                           pmin((as.numeric(af_date) - as.numeric(ecg_date))/365.25,2)))]
ecg[,':='(incd_af_3y = ifelse(incd_af==1 & c((as.numeric(af_date) - as.numeric(ecg_date)) <= 365.25*3),1,0))]
ecg[,':='(af_3y.t = ifelse(incd_af_3y==1,(as.numeric(af_date) - as.numeric(ecg_date))/365.25,
                           pmin((as.numeric(af_date) - as.numeric(ecg_date))/365.25,3)))]

# Remove prevalent AF
ecg <- ecg[prev_af==0] # 1293 
ecg <- ecg[af_2y.t > 0] # 43

# Consents
withdrawn <- fread(file='withdrawals.csv')
ecg <- ecg[!(sample_id %in% withdrawn$V1)] #0

# Calculate CHARGE-AF
ecg[,charge_startfu := age5*0.508 + (white==1)*0.465 + ht_10*0.248 + wt_15*0.115 + sbp20*0.197
    +dbp10*(-0.101) + (smoker==1)*0.359 + prev_bpmed*0.349
    +(prev_dm==1)*0.237 + (prev_hf==1)*0.701 + (prev_mi==1)*0.496]

# Load pre-fitted cox models from MGH
load(file='charge_ecg_cox.RData')
load(file='ecg_ai_cox.RData')
load(file='age_sex_cox.RData')

# Define variables for cox model
ecg[,ecg_logit := log(survival_curve_af_prediction/(1-survival_curve_af_prediction))]
ecg[,Dem.Gender.no_filter := ifelse(sex==1,'Male','Female')]
ecg[,start_fu_yrs := age_at_ecg]

# Calculate linear predictors in UKB set
ecg[,age_sex_lp := predict(age_sex,newdata=ecg,type='lp')]
ecg[,ecg_only_lp := predict(ecg_only,newdata=ecg,type='lp')]
ecg[,charge_ecg_lp := predict(charge_plus,newdata=ecg,type='lp')]

# Save out
#save(ecg,file='ecg_set.RData')

## Concordance
cstat_age_sex <- timeROC(T=ecg$af_2y.t, delta=ecg$incd_af_2y,
                        marker=ecg$age_sex_lp,cause=1,times=1.999)$AUC[2]
cstat_ecg_ai <- timeROC(T=ecg$af_2y.t, delta=ecg$incd_af_2y,
                        marker=ecg$ecg_only_lp,cause=1,times=1.999)$AUC[2]
cstat_charge <- timeROC(T=ecg$af_2y.t, delta=ecg$incd_af_2y,
                        marker=ecg$charge_startfu,cause=1,times=1.999)$AUC[2]
cstat_combo <- timeROC(T=ecg$af_2y.t, delta=ecg$incd_af_2y,
                       marker=ecg$charge_ecg_lp,cause=1,times=1.999)$AUC[2]

# CH-AI vs. CHARGE AF
cstat_boot <- timeroc_diff(time='af_2y.t',status='incd_af_2y',response1='charge_startfu',
                           response2='charge_ecg_lp',eval.t=1.999,data=ecg,runs=500)
z_cstat <- (cstat_combo[1] - cstat_charge[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. CHARGE-AF
cstat_boot <- timeroc_diff(time='af_2y.t',status='incd_af_2y',response1='charge_startfu',
                           response2='ecg_only_lp',eval.t=1.999,data=ecg,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_charge[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. CH-AI
cstat_boot <- timeroc_diff(time='af_2y.t',status='incd_af_2y',response1='ecg_only_lp',
                           response2='charge_ecg_lp',eval.t=1.999,data=ecg,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_combo[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# ECG only vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_2y.t',status='incd_af_2y',response1='ecg_only_lp',
                           response2='age_sex_lp',eval.t=1.999,data=ecg,runs=500)
z_cstat <- (cstat_ecg_ai[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# CH-AI vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_2y.t',status='incd_af_2y',response1='charge_ecg_lp',
                           response2='age_sex_lp',eval.t=1.999,data=ecg,runs=500)
z_cstat <- (cstat_combo[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# CHARGE vs. Age + Sex
cstat_boot <- timeroc_diff(time='af_2y.t',status='incd_af_2y',response1='charge_startfu',
                           response2='age_sex_lp',eval.t=1.999,data=ecg,runs=500)
z_cstat <- (cstat_charge[1] - cstat_age_sex[1])/sd(cstat_boot)
p_cstat <-  2*(1-pnorm(abs(z_cstat)))

# Create standardized variables
ecg[,':='(age_sex_std = (age_sex_lp - mean(age_sex_lp))/sd(age_sex_lp),
          charge_std = (charge_startfu - mean(charge_startfu))/sd(charge_startfu),
          ecg_std = (ecg_only_lp - mean(ecg_only_lp))/sd(ecg_only_lp),
          charge_ecg_std = (charge_ecg_lp - mean(charge_ecg_lp))/sd(charge_ecg_lp))]

# Standardized models
age_sex_std <- coxph(Surv(af_2y.t,incd_af_2y) ~ age_sex_std,data=ecg)
charge_std <- coxph(Surv(af_2y.t,incd_af_2y) ~ charge_std,data=ecg)
ecg_only_std <- coxph(Surv(af_2y.t,incd_af_2y) ~ ecg_std,data=ecg)
charge_plus_std <- coxph(Surv(af_2y.t,incd_af_2y) ~ charge_ecg_std,data=ecg)

################################# Density Plots
### RAW SCORES
# CHARGE-AF stratified by AF
x <- list(v1=ecg[incd_af_2y==1]$charge_startfu,v2=ecg[incd_af_2y==0]$charge_startfu)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(8,16,2),expand=c(0,0),limits=c(8,16)) +
  scale_y_continuous(breaks=seq(0,0.6,0.1),expand=c(0,0),limits=c(0,0.6)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='CHARGE-AF',y='Density') 
ggsave('density_charge_ukbb.pdf',
       height=2,width=2.5,units='in',scale=4)

# ECG-AI stratified by AF
x <- list(v1=ecg[incd_af_2y==1]$ecg_only_lp,v2=ecg[incd_af_2y==0]$ecg_only_lp)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-3,6,1),expand=c(0,0),limits=c(-3,6)) +
  scale_y_continuous(breaks=seq(0,0.6,0.1),expand=c(0,0),limits=c(0,0.6)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='ECG-AI',y='Density') 
ggsave('density_ecg_ukbb.pdf',
       height=2,width=2.5,units='in',scale=4)

# CH-AI stratified by AF
x <- list(v1=ecg[incd_af_2y==1]$charge_ecg_lp,v2=ecg[incd_af_2y==0]$charge_ecg_lp)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-4,6,1),expand=c(0,0),limits=c(-4,6)) +
  scale_y_continuous(breaks=seq(0,0.6,0.1),expand=c(0,0),limits=c(0,0.6)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='CHARGE-AF',y='Density') 
ggsave('density_charge_ecg_ukbb.pdf',
       height=2,width=2.5,units='in',scale=4)

######### KM
prod_af <- prodlim(Hist(af_2y.t,incd_af_2y)~1,data=ecg)

# Plot
CairoPDF(file='af_km_ukbb.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prod_af,"cuminc",ylim=c(0,0.03),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.03,0.01),axis2.las=2,axis2.cex.axis=2.5, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=2.5, #x-axis labeling parameters
     col=c("#FC4E2A"), # color of curves
     atrisk.col='black',
     confint=TRUE, # whether you want CI on the curves
     atrisk.title='',atrisk.pos=-0.4,atrisk.line=c(3.5), # position of the N at risk rows
     atrisk.cex=1.8,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,3,5), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk of AF (%)",side=2,line=-1.2,at=0.015,cex=2.5) # y-axis label
mtext("Years",side=1, line=1,cex=2.5) # x-axis label
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5) # descriptor for N at risk
dev.off()

################################# ROC curves
age_sex <- timeROC(T=ecg$af_2y.t, delta=ecg$incd_af_2y,
                  marker=ecg$age_sex_lp,cause=1,times=c(1,1.999))
ecg_ai <- timeROC(T=ecg$af_2y.t, delta=ecg$incd_af_2y,
                  marker=ecg$ecg_only_lp,cause=1,times=c(1,1.999))
charge <- timeROC(T=ecg$af_2y.t, delta=ecg$incd_af_2y,
                  marker=ecg$charge_startfu,cause=1,times=c(1,1.999))
combo <- timeROC(T=ecg$af_2y.t, delta=ecg$incd_af_2y,
                 marker=ecg$charge_ecg_lp,cause=1,times=c(1,1.999))

pdf(file='roc_compare_timeroc_ukbb.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))
plot.new() #?
plot(ecg_ai,1.999,add=T,col='#66c2a5',lwd=1.2)
par(new=TRUE)
plot(charge,1.999,add=T,col='#fc8d62',lwd=1.2)
par(new=TRUE)
plot(age_sex,1.999,add=T,col='darkgray',lwd=1.2)
par(new=TRUE)
plot(combo,1.999,add=T,col='#984ea3',lwd=1.2)
axis(1,at=seq(1,0,-0.2),labels=c('1.0','0.8','0.6','0.4','0.2','0.0'),cex.axis=1.6)
axis(2,at=seq(0,1,0.2),cex.axis=1.6,las=2)
title(xlab='1 - Specificity',line=2.5,cex.lab=1.8)
title(ylab='Sensitivity',line=3.2,cex.lab=1.8)
legend(0.5,0.2,legend=c('Age & Sex (0.728)','CHARGE-AF (0.732)','ECG-AI (0.705)','CH-AI (0.746)'),col=c('darkgray','#fc8d62','#66c2a5','#984ea3'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)
dev.off()

################################# ROC with CI over TIME
## Estimate standard errors
age_sex_ci <- timeroc_ci(data=ecg,times=c(1,1.999),time='af_2y.t',status='incd_af_2y',marker='age_sex_lp')
ecg_ai_ci <- timeroc_ci(data=ecg,times=c(1,1.999),time='af_2y.t',status='incd_af_2y',marker='ecg_only_lp')
charge_ci <- timeroc_ci(data=ecg,times=c(1,1.999),time='af_2y.t',status='incd_af_2y',marker='charge_startfu')
charge_ecg_ci <- timeroc_ci(data=ecg,times=c(1,1.999),time='af_2y.t',status='incd_af_2y',marker='charge_ecg_lp')

## Get point estimate with CI using standard errors
age_sex_est <- data.table(times=c(1,1.999),
                         AUC=age_sex$AUC,
                         lower=age_sex$AUC-1.96*apply(age_sex_ci,FUN=sd,MARGIN=2),
                         upper=age_sex$AUC+1.96*apply(age_sex_ci,FUN=sd,MARGIN=2))

ecg_ai_est <- data.table(times=c(1,1.999),
                         AUC=ecg_ai$AUC,
                         lower=ecg_ai$AUC-1.96*apply(ecg_ai_ci,FUN=sd,MARGIN=2),
                         upper=ecg_ai$AUC+1.96*apply(ecg_ai_ci,FUN=sd,MARGIN=2))

charge_est <- data.table(times=c(1,1.999),
                         AUC=charge$AUC,
                         lower=charge$AUC-1.96*apply(charge_ci,FUN=sd,MARGIN=2),
                         upper=charge$AUC+1.96*apply(charge_ci,FUN=sd,MARGIN=2))

charge_ecg_est <- data.table(times=c(1,1.999),
                             AUC=combo$AUC,
                             lower=combo$AUC-1.96*apply(charge_ecg_ci,FUN=sd,MARGIN=2),
                             upper=combo$AUC+1.96*apply(charge_ecg_ci,FUN=sd,MARGIN=2))

## The plot
pdf(file='auc_time_ukbb.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,4,1,1))
par(oma=c(1,1,1,1))

x1 <- 0.9:1.9; x2 <- 0.95:1.95; x3 <- 1.05:2.05; x4 <- 1.1:2.1

plot(x=x2,y=ecg_ai_est$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,2.2),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#66c2a5',cex=1.2)
# segments(0.7,0.724169137282253,2.2,0.724169137282253,lwd=1,col='#fc8d628C',lty=5)
# segments(0.7,0.7411131998548163,2.2,0.7411131998548163,lwd=1,col='#8da0cb8C',lty=5)
# segments(0.7,0.7008336048456216,2.2,0.7008336048456216,lwd=1,col='#66c2a58C',lty=5)
# segments(0.7,0.7201247097444149,2.2,0.7201247097444149,lwd=1,col='lightgray',lty=5)
par(new=TRUE)
plot(x=x3,y=charge_est$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,2.2),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#fc8d62',cex=1.2)
par(new=TRUE)
plot(x=x4,y=age_sex_est$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,2.2),ylim=c(0.6,0.9),bty='n',
     pch=19,col='darkgray',cex=1.2)
par(new=TRUE)
plot(x=x1,y=charge_ecg_est$AUC,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0.5,2.2),ylim=c(0.6,0.9),bty='n',
     pch=19,col='#984ea3',cex=1.2)
axis(2,las=2,cex.axis=1.5,pos=0.7,at=seq(0.6,0.9,0.1))
axis(1,at=1:2,cex.axis=1.5)
mtext("AF Prediction Window (Years)",1,cex=1.8,line=3,at=1.5)
mtext("Time-dependent AUROC",2,cex=1.8,line=1)

segments(x2,ecg_ai_est$lower,x2,ecg_ai_est$upper,lwd=1,col='#66c2a5')
segments(x3,charge_est$lower,x3,charge_est$upper,lwd=1,col='#fc8d62')
segments(x4,age_sex_est$lower,x4,age_sex_est$upper,lwd=1,col='darkgray')
segments(x1,charge_ecg_est$lower,x1,charge_ecg_est$upper,lwd=1,col='#984ea3')

legend(1.6,0.67,c('Age & Sex','CHARGE-AF','ECG-AI',"CH-AI"),bty='n',cex=1.5,pch=19,col=c('darkgray','#fc8d62','#66c2a5','#984ea3'))

dev.off()

################################# Assess calibration
############ PART 1: SLOPES
### MGH
af <- coxph(Surv(af_2y.t,incd_af_2y) ~ age_sex_lp,data=ecg)
age_sex_slope <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

af <- coxph(Surv(af_2y.t,incd_af_2y) ~ charge_startfu,data=ecg)
charge_slope <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

af <- coxph(Surv(af_2y.t,incd_af_2y) ~ ecg_only_lp,data=ecg)
ecg_slope <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

af <- coxph(Surv(af_2y.t,incd_af_2y) ~ charge_ecg_lp,data=ecg)
combo_slope <- c(af$coefficients[1],confint(af)[1],confint(af)[2])

############ PART 2: PLOTS and ICI
###################################################################### FAMILY 1: Calibrated to 2 year risk in UKBB
## CHARGE
avg_beta <- mean(ecg$charge_startfu)
res <- coxph(Surv(af_2y.t,incd_af_2y) ~ charge_startfu, data=ecg)
km <- survfit(res, data=data.frame(x1=mean(charge_startfu)),type="kaplan-meier")
train_s0 <- summary(km, times=c(2))$surv
ecg[,charge_pred2 := (1-(train_s0)^exp(charge_startfu - avg_beta))*100]
## Age/Sex
avg_beta <- mean(ecg$age_sex_lp)
res <- coxph(Surv(af_2y.t,incd_af_2y) ~ age_sex_lp, data=ecg)
km <- survfit(res, data=data.frame(x1=mean(age_sex_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(2))$surv
ecg[,age_sex_pred2 := (1-(train_s0)^exp(age_sex_lp - avg_beta))*100]
## ECG
avg_beta <- mean(ecg$ecg_only_lp)
res <- coxph(Surv(af_2y.t,incd_af_2y) ~ ecg_only_lp, data=ecg)
km <- survfit(res, data=data.frame(x1=mean(ecg_only_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(2))$surv
ecg[,ecg_pred2 := (1-(train_s0)^exp(ecg_only_lp - avg_beta))*100]
## COMBO
avg_beta <- mean(ecg$charge_ecg_lp)
res <- coxph(Surv(af_2y.t,incd_af_2y) ~ charge_ecg_lp, data=ecg)
km <- survfit(res, data=data.frame(x1=mean(charge_ecg_lp)),type="kaplan-meier")
train_s0 <- summary(km, times=c(2))$surv
ecg[,charge_ecg_pred2 := (1-(train_s0)^exp(charge_ecg_lp - avg_beta))*100]

# Use classifier to classify scores into quantiles (size 10)
ecg$age_sex_decile <- classifier(risk=ecg$age_sex_pred2,ncuts=10)
ecg$charge_decile <- classifier(risk=ecg$charge_pred2,ncuts=10)
ecg$ecg_decile <- classifier(risk=ecg$ecg_pred2,ncuts=10)
ecg$charge_ecg_decile <- classifier(risk=ecg$charge_ecg_pred2,ncuts=10)

### CALCULATE OBSERVED RISK IN EACH QUANTILE
setDF(ecg)
age_sex_obv <- survivor(data=ecg,risk_data="age_sex_decile",time='af_2y.t',status='incd_af_2y',eval.t=5)
charge_obv <- survivor(data=ecg,risk_data="charge_decile",time='af_2y.t',status='incd_af_2y',eval.t=5)
ecg_obv <- survivor(data=ecg,risk_data="ecg_decile",time='af_2y.t',status='incd_af_2y',eval.t=5)
charge_ecg_obv <- survivor(data=ecg,risk_data="charge_ecg_decile",time='af_2y.t',status='incd_af_2y',eval.t=5)
setDT(ecg)

### CALCULATE AVERAGE PREDICTED RISK IN EACH QUANTILE
age_sex_pred <- ecg[,mean(age_sex_pred2),by="age_sex_decile"][order(age_sex_decile)]
charge_pred <- ecg[,mean(charge_pred2),by="charge_decile"][order(charge_decile)]
ecg_pred <- ecg[,mean(ecg_pred2),by="ecg_decile"][order(ecg_decile)]
charge_ecg_pred <- ecg[,mean(charge_ecg_pred2),by="charge_ecg_decile"][order(charge_ecg_decile)]

### Fit adaptive hazard model for Age/Sex
ecg$cox.2yr.cll_age_sex <- log(-log(1-ecg$age_sex_pred2/100))
calibrate.cox_age_sex <- hare(data=ecg$af_2y.t,delta=ecg$incd_af_2y,
                             cov=as.matrix(ecg$cox.2yr.cll_age_sex))
predict.grid.cox_age_sex <- seq(quantile(ecg$age_sex_pred2/100,probs=0.01),
                               quantile(ecg$age_sex_pred2/100,probs=0.99),length=100)
predict.grid.cox.cll_age_sex <- log(-log(1-predict.grid.cox_age_sex))
predict.calibrate.cox_age_sex <- phare(2,predict.grid.cox.cll_age_sex,calibrate.cox_age_sex)

# Plots for visualization
pdf(file='cal_age_sex_ukbb.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- age_sex_pred$V1
y <- do.call(rbind,age_sex_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='darkgray',
     xaxt='n',xlim=c(0,5),ylim=c(0,5),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,5,1),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,5,1),las=1)

segments(-1,-1,6,6,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_age_sex*100,predict.calibrate.cox_age_sex*100,type="l",lty=1,col="darkgray",
     xlim=c(0,5),ylim=c(0,5),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 2 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 2 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_age_sex <- phare(2,ecg$cox.2yr.cll_age_sex,calibrate.cox_age_sex)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.2yr.cox <- mean(abs(ecg$age_sex_pred2/100 - predict.calibrate.cox_age_sex))
E50.2yr.cox <- median(abs(ecg$age_sex_pred2/100 - predict.calibrate.cox_age_sex))
E90.2yr.cox <- quantile(abs(ecg$age_sex_pred2/100 - predict.calibrate.cox_age_sex),probs=0.9)

### Fit adaptive hazard model for CHARGE-AF
ecg$cox.2yr.cll_charge <- log(-log(1-ecg$charge_pred2/100))
calibrate.cox_charge <- hare(data=ecg$af_2y.t,delta=ecg$incd_af_2y,
                             cov=as.matrix(ecg$cox.2yr.cll_charge))
predict.grid.cox_charge <- seq(quantile(ecg$charge_pred2/100,probs=0.01),
                               quantile(ecg$charge_pred2/100,probs=0.99),length=100)
predict.grid.cox.cll_charge <- log(-log(1-predict.grid.cox_charge))
predict.calibrate.cox_charge <- phare(2,predict.grid.cox.cll_charge,calibrate.cox_charge)

# Plots for visualization
pdf(file='cal_charge_ukbb.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- charge_pred$V1
y <- do.call(rbind,charge_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#fc8d62',
     xaxt='n',xlim=c(0,5),ylim=c(0,5),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,5,1),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,5,1),las=1)

segments(-1,-1,6,6,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_charge*100,predict.calibrate.cox_charge*100,type="l",lty=1,col="#fc8d62",
     xlim=c(0,5),ylim=c(0,5),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 2 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 2 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_charge <- phare(2,ecg$cox.2yr.cll_charge,calibrate.cox_charge)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.2yr.cox <- mean(abs(ecg$charge_pred2/100 - predict.calibrate.cox_charge))
E50.2yr.cox <- median(abs(ecg$charge_pred2/100 - predict.calibrate.cox_charge))
E90.2yr.cox <- quantile(abs(ecg$charge_pred2/100 - predict.calibrate.cox_charge),probs=0.9)

### Fit adaptive hazard model for ECG score
ecg$cox.2yr.cll_ecg <- log(-log(1-ecg$ecg_pred2/100))
calibrate.cox_ecg <- hare(data=ecg$af_2y.t,delta=ecg$incd_af_2y,
                          cov=as.matrix(ecg$cox.2yr.cll_ecg))
predict.grid.cox_ecg <- seq(quantile(ecg$ecg_pred2/100,probs=0.01),
                            quantile(ecg$ecg_pred2/100,probs=0.99),length=100)
predict.grid.cox.cll_ecg <- log(-log(1-predict.grid.cox_ecg))
predict.calibrate.cox_ecg <- phare(2,predict.grid.cox.cll_ecg,calibrate.cox_ecg)

# Plots for visualization
pdf(file='cal_ecg_ukbb.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- ecg_pred$V1
y <- do.call(rbind,ecg_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#66c2a5',
     xaxt='n',xlim=c(0,5),ylim=c(0,5),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,5,1),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,5,1),las=1)

segments(-1,-1,6,6,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_ecg*100,predict.calibrate.cox_ecg*100,type="l",lty=1,col="#66c2a5",
     xlim=c(0,5),ylim=c(0,5),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 2 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 2 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_ecg <- phare(2,ecg$cox.2yr.cll_ecg,calibrate.cox_ecg)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.2yr.cox <- mean(abs(ecg$ecg_pred2/100 - predict.calibrate.cox_ecg))
E50.2yr.cox <- median(abs(ecg$ecg_pred2/100 - predict.calibrate.cox_ecg))
E90.2yr.cox <- quantile(abs(ecg$ecg_pred2/100 - predict.calibrate.cox_ecg),probs=0.9)

### Fit adaptive hazard model for CHARGE ECG score
ecg$cox.2yr.cll_charge_ecg <- log(-log(1-ecg$charge_ecg_pred2/100))
calibrate.cox_charge_ecg <- hare(data=ecg$af_2y.t,delta=ecg$incd_af_2y,
                                 cov=as.matrix(ecg$cox.2yr.cll_charge_ecg))
predict.grid.cox_charge_ecg <- seq(quantile(ecg$charge_ecg_pred2/100,probs=0.01),
                                   quantile(ecg$charge_ecg_pred2/100,probs=0.99),length=100)
predict.grid.cox.cll_charge_ecg <- log(-log(1-predict.grid.cox_charge_ecg))
predict.calibrate.cox_charge_ecg <- phare(2,predict.grid.cox.cll_charge_ecg,calibrate.cox_charge_ecg)

# Plots for visualization
pdf(file='cal_charge_ecg_ukbb.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- charge_ecg_pred$V1
y <- do.call(rbind,charge_ecg_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#8da0cb',
     xaxt='n',xlim=c(0,5),ylim=c(0,5),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,5,1),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,5,1),las=1)

segments(-1,-1,6,6,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_charge_ecg*100,predict.calibrate.cox_charge_ecg*100,type="l",lty=1,col="#8da0cb",
     xlim=c(0,5),ylim=c(0,5),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 2 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 2 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_charge_ecg <- phare(2,ecg$cox.2yr.cll_charge_ecg,calibrate.cox_charge_ecg)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.2yr.cox <- mean(abs(ecg$charge_ecg_pred2/100 - predict.calibrate.cox_charge_ecg))
E50.2yr.cox <- median(abs(ecg$charge_ecg_pred2/100 - predict.calibrate.cox_charge_ecg))
E90.2yr.cox <- quantile(abs(ecg$charge_ecg_pred2/100 - predict.calibrate.cox_charge_ecg),probs=0.9)

###################################################################### FAMILY 2: Raw for CH-AI
## ECG
ecg[,ecg_pred2_raw := (1-(0.9831352)^exp(ecg_only_lp - (6.951522e-15)))*100]
## COMBO
ecg[,charge_ecg_pred2_raw := (1-(0.9851827)^exp(ecg_only_lp - (-8.472139e-16)))*100]

# Use classifier to classify scores into quantiles (size 10)
ecg$ecg_decile_raw <- classifier(risk=ecg$ecg_pred2_raw,ncuts=10)
ecg$charge_ecg_decile_raw <- classifier(risk=ecg$charge_ecg_pred2_raw,ncuts=10)

### CALCULATE OBSERVED RISK IN EACH QUANTILE
setDF(ecg)
ecg_obv <- survivor(data=ecg,risk_data="ecg_decile_raw",time='af_2y.t',status='incd_af_2y',eval.t=5)
charge_ecg_obv <- survivor(data=ecg,risk_data="charge_ecg_decile_raw",time='af_2y.t',status='incd_af_2y',eval.t=5)
setDT(ecg)

### CALCULATE AVERAGE PREDICTED RISK IN EACH QUANTILE
ecg_pred <- ecg[,mean(ecg_pred2_raw),by="ecg_decile_raw"][order(ecg_decile_raw)]
charge_ecg_pred <- ecg[,mean(charge_ecg_pred2_raw),by="charge_ecg_decile_raw"][order(charge_ecg_decile_raw)]

### Fit adaptive hazard model for ECG score
ecg$cox.2yr.cll_ecg <- log(-log(1-ecg$ecg_pred2_raw/100))
calibrate.cox_ecg <- hare(data=ecg$af_2y.t,delta=ecg$incd_af_2y,
                          cov=as.matrix(ecg$cox.2yr.cll_ecg))
predict.grid.cox_ecg <- seq(quantile(ecg$ecg_pred2_raw/100,probs=0.01),
                            quantile(ecg$ecg_pred2_raw/100,probs=0.99),length=100)
predict.grid.cox.cll_ecg <- log(-log(1-predict.grid.cox_ecg))
predict.calibrate.cox_ecg <- phare(2,predict.grid.cox.cll_ecg,calibrate.cox_ecg)

# Plots for visualization
pdf(file='cal_ecg_ukbb_raw.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- ecg_pred$V1
y <- do.call(rbind,ecg_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#66c2a5',
     xaxt='n',xlim=c(0,35),ylim=c(0,35),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,35,5),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,35,5),las=1)

segments(-1,-1,36,36,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_ecg*100,predict.calibrate.cox_ecg*100,type="l",lty=1,col="#66c2a5",
     xlim=c(0,35),ylim=c(0,35),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 2 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 2 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_ecg <- phare(2,ecg$cox.2yr.cll_ecg,calibrate.cox_ecg)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.2yr.cox <- mean(abs(ecg$ecg_pred2_raw/100 - predict.calibrate.cox_ecg))
E50.2yr.cox <- median(abs(ecg$ecg_pred2_raw/100 - predict.calibrate.cox_ecg))
E90.2yr.cox <- quantile(abs(ecg$ecg_pred2_raw/100 - predict.calibrate.cox_ecg),probs=0.9)

### Fit adaptive hazard model for CHARGE ECG score
ecg$cox.2yr.cll_charge_ecg <- log(-log(1-ecg$charge_ecg_pred2_raw/100))
calibrate.cox_charge_ecg <- hare(data=ecg$af_2y.t,delta=ecg$incd_af_2y,
                                 cov=as.matrix(ecg$cox.2yr.cll_charge_ecg))
predict.grid.cox_charge_ecg <- seq(quantile(ecg$charge_ecg_pred2_raw/100,probs=0.01),
                                   quantile(ecg$charge_ecg_pred2_raw/100,probs=0.99),length=100)
predict.grid.cox.cll_charge_ecg <- log(-log(1-predict.grid.cox_charge_ecg))
predict.calibrate.cox_charge_ecg <- phare(2,predict.grid.cox.cll_charge_ecg,calibrate.cox_charge_ecg)

# Plots for visualization
pdf(file='cal_charge_ecg_ukbb_raw.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))

x <- charge_ecg_pred$V1
y <- do.call(rbind,charge_ecg_obv)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#8da0cb',
     xaxt='n',xlim=c(0,35),ylim=c(0,35),pch=19,cex=1.5,bty='n')

axis(1,at=seq(0,35,5),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,35,5),las=1)

segments(-1,-1,36,36,lwd=1.2,lty=2)

par(new=TRUE)
plot(predict.grid.cox_charge_ecg*100,predict.calibrate.cox_charge_ecg*100,type="l",lty=1,col="#8da0cb",
     xlim=c(0,35),ylim=c(0,35),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

mtext("Predicted risk of AF at 2 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 2 years (%)",side=2,cex=1.8,line=4.5)

dev.off()

# Quant measures
predict.calibrate.cox_charge_ecg <- phare(2,ecg$cox.2yr.cll_charge_ecg,calibrate.cox_charge_ecg)

# Predicted probability of death within 1 year for all subjects in
# validation sample
ICI.2yr.cox <- mean(abs(ecg$charge_ecg_pred2_raw/100 - predict.calibrate.cox_charge_ecg))
E50.2yr.cox <- median(abs(ecg$charge_ecg_pred2_raw/100 - predict.calibrate.cox_charge_ecg))
E90.2yr.cox <- quantile(abs(ecg$charge_ecg_pred2_raw/100 - predict.calibrate.cox_charge_ecg),probs=0.9)

########## DENSITY PLOTS
### Age/sex
x <- list(v1=ecg[incd_af_2y==1]$age_sex_pred2_cal,v2=ecg[incd_af_2y==0]$age_sex_pred2_cal)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,3,1),expand=c(0,0),limits=c(0,3)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5),expand=c(0,0),limits=c(0,1.5)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted AF Risk (%)',y='Density') 
ggsave('pred_density_agesex_ukbb.pdf',
       height=2,width=2.5,units='in',scale=4)

### CHARGE-AF
x <- list(v1=ecg[incd_af_2y==1]$charge_pred2_cal,v2=ecg[incd_af_2y==0]$charge_pred2_cal)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,5,1),expand=c(0,0),limits=c(0,5)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5),expand=c(0,0),limits=c(0,1.5)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted AF Risk (%)',y='Density') 
ggsave('pred_density_charge_ukbb.pdf',
       height=2,width=2.5,units='in',scale=4)

### ECG-AI
x <- list(v1=ecg[incd_af_2y==1]$ecg_pred2_cal,v2=ecg[incd_af_2y==0]$ecg_pred2_cal)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,5,1),expand=c(0,0),limits=c(0,5)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5),expand=c(0,0),limits=c(0,1.5)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted AF Risk (%)',y='Density') 
ggsave('pred_density_ecg_ukbb.pdf',
       height=2,width=2.5,units='in',scale=4)

### CH-AI
x <- list(v1=ecg[incd_af_2y==1]$charge_ecg_pred2_cal,v2=ecg[incd_af_2y==0]$charge_ecg_pred2_cal)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,5,1),expand=c(0,0),limits=c(0,5)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5),expand=c(0,0),limits=c(0,1.5)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('AF','No AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted AF Risk (%)',y='Density') 
ggsave('~/Documents/MGH Research/af_prediction/pred_density_chai_ukbb.pdf',
       height=2,width=2.5,units='in',scale=4)

# Save processed output if desired
#save(ecg,file='ecg_set.RData')
