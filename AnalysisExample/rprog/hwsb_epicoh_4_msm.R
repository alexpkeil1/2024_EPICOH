setwd("~/Documents/Talks/Workshops/2024_EPICOH/AnalysisExample")
source("rprog/hwsb_epicoh_1_datamgmt.R")
library(dplyr) # data management functions from tidyverse
library(survival)

# Programming note: In many steps, I am writing code to a file and then running that file, 
#  this creates output for slides while also running the code without copy/paste errors. 
# It's not suggested as a general programming practice

########################################################################
# Finding 30 potential occupational limits for a marginal structural model
# strategy: pick limits that are well within the bounds of the data
########################################################################
# maximum annual exposures for each individual
maxexposures = as.numeric(tapply(sim_cohort$x, sim_cohort$id, max))

limits = quantile(maxexposures[maxexposures>0], seq(.20,0.90,length.out=30)) # 30 limits all within the range of observed lifetime max exposures

# do all analysis: note this is wrapped in a function to facilitate bootstrapping, if needed
doipw = function(sim_cohort, limits){
  
  #1b) create copy of the dataset for each limit
  clones <- do.call(rbind,
                    lapply(1:length(limits), function(x) {  
                      df = sim_cohort                  
                      df$limit = as.numeric(limits[x])
                      df$cloneid = paste0(x, "_", df$id)
                      df                    
                    }
                    ))
  #2) artificially censor data (keeping first observation that is censored)
  cens_data <- clones %>%
    group_by(cloneid) %>%
    mutate(
      cens = pmin(1,cumsum(x > limit)),
      drop = pmin(2, cumsum(cens))
    ) %>%
    group_by() %>%
    filter(
      drop < 2
    )
  
  # 3) create 1/0 weights to use for confounding/censoring during follow-up (alternative is to create new datasets)
  cens_data <- group_by(cens_data, cloneid) %>%
    mutate(
      one = 1,
      nobs = cumsum(one),
      fobs____ = as.numeric(nobs == 1)
    ) %>%
    mutate(
      conf_weight = as.numeric(nobs  == 1), # this duplicates fobs____ but is kept for clarity
      fu_weight = as.numeric((nobs  > 1) & (atwork == 1)),    
      dconf = 0,
      dcens = 0,
      intervened = 0,
    ) %>%
    select(-c(one, nobs)) %>%
    group_by()
  # check: which(tapply(sim_cohort$x, sim_cohort$id, max) < limit & tapply(sim_cohort$atwork, sim_cohort$id, sum)>3)[2] # index 427
  # check: names(tapply(sim_cohort$x, sim_cohort$id, max))[427]
  # check: print(select(cens_data, c(id, x, limit, conf_weight, fu_weight, cens, fobs____)) %>% filter(id==427), n=50)
  # check: print(select(filter(sim_cohort, id==427), c(id, x, atwork)))
  
  
  #4) fit censoring models (can include pre-baseline exposure in practice)

  # fit models if there is more than a small amount of censoring (seems to work either way, but this avoids convergence problems)
  for (l in limits){
    limidx = which(cens_data$limit == l)
    tempdat = cens_data[limidx,]
    tempdat$mxl = tempdat$cxl / (tempdat$cumatworkl + as.numeric(tempdat$time==1))
    agekn0 = attr(rcspline.eval(filter(tempdat, conf_weight==1)$age, nk = 4), "knots")
    yearkn0 = attr(rcspline.eval(filter(tempdat, conf_weight==1)$year, nk = 4), "knots")
    agekn = attr(rcspline.eval(filter(tempdat, fu_weight==1)$age, nk = 4), "knots")
    yearkn = attr(rcspline.eval(filter(tempdat, fu_weight==1)$year, nk = 4), "knots")
    #
    if(sum(tempdat[tempdat$conf_weight==1,"cens"]) > 10){
      summary(confnmod <- glm(cens ~ age + rcspline.eval(age, knots=agekn0) , data = tempdat, weight=conf_weight, family=binomial()))
      summary(confdmod <- glm(cens ~ age + rcspline.eval(age, knots=agekn0) + wagestatus + male + race, data = tempdat, weight=conf_weight, family=binomial()))
      # need to use base R operations here
      # only get non-zero predictions if "conf_weight" == 1 (eligible to be "censored" at baseline)
      cens_data[limidx,"nconf"] = tempdat$conf_weight*as.numeric(predict(confnmod, type="response"))
      cens_data[limidx,"dconf"] = tempdat$conf_weight*as.numeric(predict(confdmod, type="response"))
    } else{
      cens_data[limidx,"nconf"] = 0
      cens_data[limidx,"dconf"] = 0
    }
    if(sum(tempdat[tempdat$fu_weight==1,"cens"]) > 1){
      summary(censnmod <- glm(cens ~ age + rcspline.eval(age, knots=agekn), data = tempdat, weight=fu_weight, family=binomial()))
      summary(censdmod <- glm(cens ~  mxl + cumatworkl + age + rcspline.eval(age, knots=agekn) + wagestatus + male + race, data = tempdat, weight=fu_weight, family=binomial()))
      # only get non-zero predictions if "fu_weight" == 1 (eligible to be censored during follow-up)
      cens_data[limidx,"ncens"] = tempdat$fu_weight*as.numeric(predict(censnmod, type="response")) 
      cens_data[limidx,"dcens"] = tempdat$fu_weight*as.numeric(predict(censdmod, type="response")) 
    } else{
      cens_data[limidx,"ncens"] = 0
      cens_data[limidx,"dcens"] = 0
    }
  }
  
  # create final weights
  combined_wtd_data <- cens_data %>% 
    mutate(
      wtcontr = case_when(
        # weight: (unstabilized) inverse probability of NOT being censored
        # note regarding unstabilized weights: Cain et al 2010 show that weight stabilization is not guaranteed to reduce variance
        ((fobs____ == 1) & (atwork==1)) ~ (1-cens)*(1-nconf)/(1-dconf),
        ((fobs____ == 0) & (atwork==1)) ~ (1-cens)*(1-ncens)/(1-dcens),
        .default=1
      )
    ) %>%
    select(id,cloneid,agein,age,d1,d2,wtcontr,x,cumx,atwork,cens,limit, wagestatus, male, race, fobs____) %>% # optional step here to reduce comp. burden
    group_by(cloneid) %>%
    mutate(
      ipw = cumprod(wtcontr)
    ) %>%
    group_by() %>%
    mutate(
      # weight truncation at 10.0 following Cain et al 2010 (should be a user option)
      ipw_trunc = pmin(ipw, 10.0)
    )

  # check mean weights by intervention
  # note mean is taken across all possible observations, even those with weights = 0
  N = nrow(sim_cohort)
  Nid = length(unique(sim_cohort$id))
  wtdx = data.frame(
    limit = limits,
    mean_cumx = do.call(c,lapply(limits, function(x) mean(combined_wtd_data[combined_wtd_data$limit == x  & combined_wtd_data$cens == 0,]$cumx))),
    n_conf = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$fobs____==1,]$cens))),
    n_cens = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$fobs____==0,]$cens))),
    n_d1w = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw*combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$d1))),
    n_d2w = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw*combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$d2))),
    mean_ipw = do.call(c,lapply(limits, function(x) mean(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw)))
  )
  
  combined_wtd_data$event <- factor(combined_wtd_data$d2 + combined_wtd_data$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
  
  # marginal structural policy Cox model: effect of incremental change in the limit
  ttr1 = coxph(Surv(agein, age, d1)~limit, data=filter(combined_wtd_data, ipw>0), 
              id=cloneid, weight=ipw, cluster=id, x=FALSE, y=FALSE)
  ttr2 = coxph(Surv(agein, age, d2)~limit, data=filter(combined_wtd_data, ipw>0), 
              id=cloneid, weight=ipw, cluster=id, x=FALSE, y=FALSE)
  cat('
    ttr = coxph(Surv(agein, age, d1)~limit, data=filter(combined_wtd_data, ipw>0), 
              id=cloneid, weight=ipw, cluster=id)
      ', file = "routput/mscoxmodelcode.txt")
  #source("routput/mscoxmodelcode.txt") # does not work within a function
  list(msmr1=ttr1, msmr2=ttr2, dx=wtdx)
}

# 4) fit the MSM
# Note: this is a big model: on my computer I have to increase the "vector memory limit" for the robust variance estimator
# mem.maxVSize(vsize = 32768)
msm_estimates = doipw(sim_cohort, limits)
# check diagnostics
msm_estimates$dx

# ln-HR: effect of one unit change in the exposure limit (Marginal structural policy model)
ttr1 = msm_estimates$msmr1
summary(ttr1)

ttr2 = msm_estimates$msmr2
summary(ttr2)

# compare with standard approach
std_hr1 <- coxph(Surv(agein, age, d1)~cxl + male + race + wagestatus, data=sim_cohort)
summary(std_hr1)

std_hr2 <- coxph(Surv(agein, age, d2)~cxl + male + race + wagestatus, data=sim_cohort)
summary(std_hr2)

######### BOOTSTRAPPING #########
# We can optionally do bootstrapping for confidence intervals, which will generally be better
# than robust confidence intervals except in very large sample sizes
# Code to do this is not explicitly done here but it is a straightforward set of steps:
# 1. sample, with replacement, 10,000 individuals from the population (original size = 10,000)
# 2. Carry out clone-censor-weight + estimate effects on the sample from #1 (using doipw function above)
# 3. Repeat 1 and 2 200+ times, recording effect estimate for each sample
# 4. Bootstrap standard error is the standard deviation of the effect estimates across the 200+ bootstrap samples

