# Function to generate survival estimates per AF risk quantile
survivor <- function(data,risk_data,time,status,eval.t){
  est <- rep(NA,times=length(unique(data[,risk_data])))
  lower <- rep(NA,times=length(unique(data[,risk_data])))
  upper <- rep(NA,times=length(unique(data[,risk_data])))
  for (i in 1:length(unique(data[,risk_data]))){
    km <- survfit(Surv(data[data[,risk_data]==unique(data[,risk_data])[order(unique(data[,risk_data]))][i],time],
                       data[data[,risk_data]==unique(data[,risk_data])[order(unique(data[,risk_data]))][i],status]) ~ 1)
    est[i] <- 1-stepfun(km$time, c(1, km$surv))(eval.t)
    upper[i] <- 1-stepfun(km$time, c(1, km$lower))(eval.t)
    lower[i] <- 1-stepfun(km$time, c(1, km$upper))(eval.t)
  }
  return(data.frame(est=est,upper=upper,lower=lower))
}

# Cumulative risk of AF at 5 years
ir <- function(data,strata,time,status){
  out <- list()
  for (i in 1:length(unique(data[,get(strata)]))){
    subset <- data[get(strata)==unique(data[,get(strata)])[i]]
    events <- nrow(subset[get(status)==1])
    pt <- sum(subset[,get(time)])
    rate <- events/pt
    out[[i]] <- c(paste0(unique(data[,get(strata)])[i]),events,pt,rate,rate-1.96*(rate/sqrt(events)),rate+1.96*(rate/sqrt(events)))
  }
  return(do.call(rbind,out))
}

# Bootstrap function for difference in c-stats
cstat_diff <- function(time,status,response1,
                       response2,data,runs,size=nrow(data)){
  out <- list()
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    cstat1 <- summary(coxph(Surv(sample[,get(time)],sample[,get(status)]) ~ sample[,get(response1)],data=sample))$concordance[1]
    cstat2 <- summary(coxph(Surv(sample[,get(time)],sample[,get(status)]) ~ sample[,get(response2)],data=sample))$concordance[1]
    out[[i]] <- cstat1 - cstat2
    if (i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(unlist(out))
}

# Quantile sorter
classifier <- function(risk,ncuts){
  cuts <- quantile(risk,probs=seq(0,1,1/ncuts))
  index <- rep(NA,length(risk))
  for (i in 1:(length(cuts)-1)){
    for (j in 1:length(risk)){
      index[j] <- ifelse(risk[j] >= cuts[i],i,index[j])}}
  return(index)
}

# TimeROC CI
timeroc_ci <- function(data,times,time,status,marker,event_marker=1,runs=200){
  out <- list()
  for (i in 1:runs){
    n <- 1; est <- rep(NA,length(times))
    sample <- data[sample(1:nrow(data),size=nrow(data),replace=TRUE)]
    for (j in times){
      est[n] <- timeROC(T=sample[,get(time)], delta=sample[,get(status)],
                        marker=sample[,get(marker)],cause=event_marker,times=j)$AUC[2]
      n <- n+1
    }
    out[[i]] <- est
  }
  output <- do.call(rbind,out)
  return(output)
}

# Bootstrap function for difference in c-stats using timeROC
timeroc_diff <- function(time,status,response1,event_marker=1,
                         response2,data,runs,eval.t,size=nrow(data)){
  out <- list()
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    cstat1 <- timeROC(T=sample[,get(time)], delta=sample[,get(status)],marker=sample[,get(response1)],cause=event_marker,times=eval.t)$AUC[2]
    cstat2 <- timeROC(T=sample[,get(time)], delta=sample[,get(status)],marker=sample[,get(response2)],cause=event_marker,times=eval.t)$AUC[2]
    out[[i]] <- cstat1 - cstat2
    if (i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(unlist(out))
}

# Bootstrap function for difference in AUPRC
prc_diff <- function(time,status,response1,
                     response2,data,runs,eval.t,size=nrow(data)){
  out <- list()
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    prc1 <- APSurv(stime=sample[,get(time)],status=sample[,get(status)],marker=sample[,get(response1)],t0.list=eval.t,Plot=FALSE)
    prc_a <- prc1$ap_summary[6]
    prc2 <- APSurv(stime=sample[,get(time)],status=sample[,get(status)],marker=sample[,get(response2)],t0.list=eval.t,Plot=FALSE)
    prc_b <- prc2$ap_summary[6]
    out[[i]] <- prc_a - prc_b
    if (i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(unlist(out))
}

# Bootstrap function for difference in c-stats using timeROC
subgroup_discrim <- function(data,cut_column,cuts,time_var,eval.t,status,markers,runs){
  output <- list()
  for (i in 1:length(cuts)){
    
    if (i == 1){subset <- data[get(cut_column) < cuts[i]]}
    else {subset <- data[get(cut_column) < cuts[i] & get(cut_column) >= cuts[i-1]]}
    
    auc_out <- data.table(); prc_out <- data.table()
    for (j in markers){
      auc <- timeROC(T=subset[,get(time_var)], delta=subset[,get(status)],
                     marker=subset[,get(j)],cause=1,times=eval.t)$AUC[2]
      auc_se <- timeroc_ci(data=subset,time=time_var,times=eval.t,status=status,marker=j)
      auc_est <- data.table(model=paste(j),
                            n_total=nrow(subset),
                            n_events=sum(subset[,get(status)]),
                            times=eval.t,
                            AUC=auc,
                            lower=auc-1.96*apply(auc_se,FUN=sd,MARGIN=2),
                            upper=auc+1.96*apply(auc_se,FUN=sd,MARGIN=2))
      auc_out <- rbind(auc_out,auc_est)
      prc <- APSurv(stime=subset[,get(time_var)],status=subset[,get(status)],
                    marker=subset[,get(j)],t0.list=c(1,eval.t),method='bootstrap',B=runs,Plot=FALSE)
      prc_est <- data.table(model=paste(j),
                            times=eval.t,
                            n_total=nrow(subset),
                            n_events=sum(subset[,get(status)]),
                            no_skill=prc$ap_summary[4],
                            PRC=prc$ap_summary[6],
                            lower=prc$ap_summary[8],
                            upper=prc$ap_summary[10])
      prc_out <- rbind(prc_out,prc_est)
    }
    output[[paste(cuts[i])]] <- list(auc_out,prc_out)
  }
  return(output)
}

# Bootstrap function for difference in c-stats using timeROC
subgroup_compare <- function(data,cut_column,cuts,time_var,eval_t,status,marker1,marker2,B,all_pop=TRUE){
  output <- list()
  for (i in 1:length(cuts)){
    
    if (i == 1){subset <- data[get(cut_column) < cuts[i]]}
    else {subset <- data[get(cut_column) < cuts[i] & get(cut_column) >= cuts[i-1]]}

    auc1 <- timeROC(T=subset[,get(time_var)], delta=subset[,get(status)],
                   marker=subset[,get(marker1)],cause=1,times=eval_t)$AUC[2]
    auc2 <- timeROC(T=subset[,get(time_var)], delta=subset[,get(status)],
                    marker=subset[,get(marker2)],cause=1,times=eval_t)$AUC[2]
    auc_diff <- auc1 - auc2
    auc_se <- timeroc_diff(time=time_var,status=status,response1=marker1,
                                response2=marker2,eval.t=eval_t,data=subset,runs=B)
    auc_lower <- auc_diff - 1.96*sd(auc_se)
    auc_upper <- auc_diff + 1.96*sd(auc_se)
    auc_out <- data.table(stratum = paste0(cuts[i]), auc_diff = auc_diff, auc_se = sd(auc_se), auc_lower = auc_lower, auc_upper = auc_upper)

    prc1 <- APSurv(stime=subset[,get(time_var)],status=subset[,get(status)],
                   marker=subset[,get(marker1)],t0.list=c(1,eval_t),Plot=FALSE)
    prc2 <- APSurv(stime=subset[,get(time_var)],status=subset[,get(status)],
                   marker=subset[,get(marker2)],t0.list=c(1,eval_t),Plot=FALSE)
    
    prc_diff <- prc1$ap_summary[6] - prc2$ap_summary[6]
    prc_se <- prc_diff(data=subset,time=time_var,status=status,response1=marker1,
                       response2=marker2,runs=B,eval.t=c(1,eval_t))
    prc_lower <- prc_diff - 1.96*sd(prc_se)
    prc_upper <- prc_diff + 1.96*sd(prc_se)
    prc_out <- data.table(stratum = paste0(cuts[i]), prc_diff = prc_diff, prc_se = sd(prc_se), prc_lower = prc_lower, prc_upper = prc_upper)
    
    output[[paste(cuts[i])]] <- list(auc_out,prc_out)
  }
  if (all_pop == TRUE){
    auc1 <- timeROC(T=data[,get(time_var)], delta=data[,get(status)],
                    marker=data[,get(marker1)],cause=1,times=eval_t)$AUC[2]
    auc2 <- timeROC(T=data[,get(time_var)], delta=data[,get(status)],
                    marker=data[,get(marker2)],cause=1,times=eval_t)$AUC[2]
    auc_diff <- auc1 - auc2
    auc_se <- timeroc_diff(time=time_var,status=status,response1=marker1,
                           response2=marker2,eval.t=eval_t,data=data,runs=B)
    auc_lower <- auc_diff - 1.96*sd(auc_se)
    auc_upper <- auc_diff + 1.96*sd(auc_se)
    auc_out <- data.table(stratum = 'All', auc_diff = auc_diff, auc_se = sd(auc_se), auc_lower = auc_lower, auc_upper = auc_upper)
    
    prc1 <- APSurv(stime=data[,get(time_var)],status=data[,get(status)],
                   marker=data[,get(marker1)],t0.list=c(1,eval_t),Plot=FALSE)
    prc2 <- APSurv(stime=data[,get(time_var)],status=data[,get(status)],
                   marker=data[,get(marker2)],t0.list=c(1,eval_t),Plot=FALSE)
    
    prc_diff <- prc1$ap_summary[6] - prc2$ap_summary[6]
    prc_se <- prc_diff(data=data,time=time_var,status=status,response1=marker1,
                       response2=marker2,runs=B,eval.t=c(1,eval_t))
    prc_lower <- prc_diff - 1.96*sd(prc_se)
    prc_upper <- prc_diff + 1.96*sd(prc_se)
    prc_out <- data.table(stratum = 'All', prc_diff = prc_diff, prc_se = sd(prc_se), prc_lower = prc_lower, prc_upper = prc_upper)
    
    output[['All']] <- list(auc_out,prc_out)
  }
  return(output)
}

# Cumulative risk function
cuminc <- function(data,time,status,response){
  obj <- survfit(Surv(data[,get(time)],data[,get(status)]) ~ 1)
  ci <- c((1-obj$surv[length(obj$surv)])*100,
          (1-obj$upper[length(obj$upper)])*100,
          (1-obj$lower[length(obj$lower)])*100)
  n_event <- sum(data[,get(status)])
  n_total <- nrow(data)
  out <- data.table(n_total=n_total,n_event=n_event,ci=ci[1],ci_lower=ci[2],ci_upper=ci[3])
  return(out)
}

# Time dependent AUPRC curve
auprc <- function(data,time,status,marker,eval.t,tolerance=2){
cuts = unique(round(data[,get(marker)],tolerance))[order(unique(round(data[,get(marker)],tolerance)))]
points <- data.table(sens=NULL,ppv=NULL); n<-0
for (i in cuts){
  sens <- SeSpPPVNPV(cutpoint=i,T=data[,get(time)],delta=data[,get(status)],marker=data[,get(marker)],cause=1,times=eval.t)$TP[2]
  ppv <- SeSpPPVNPV(cutpoint=i,T=data[,get(time)],delta=data[,get(status)],marker=data[,get(marker)],cause=1,times=eval.t)$PPV[2]
  points <- rbind(points,data.table(sens=sens,ppv=ppv))
  if (n %% 50 == 0){print(paste0("I just finished ",n, "out of ",length(cuts)))}
  n <- n+1}
return(points)
}

# Linear predictor to probability
lp_to_prob <- function(data,lp,pred_var,lp_values){
out <- c()
  for (i in 1:length(lp)){
    out[i] <- mean(data[round(get(lp[i]),2) == round(lp_values[i],2),get(pred_var[i])])
  }
return(out)
}
