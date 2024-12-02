setwd("~/Documents/Talks/Workshops/2024_EPICOH/AnalysisExample")
source("rprog/hwsb_epicoh_1_datamgmt.R")
library(dplyr) # data management functions from tidyverse
library(survival)

# Programming note: In many steps, I am writing code to a file and then running that file, 
#  this creates output for slides while also running the code without copy/paste errors. 
# It's not suggested as a general programming practice


# clone the data for an intervention of an annual limit at 2
# create a "clone id" that is distinct to each clone
sink("routput/clonecode.txt", split=FALSE)
cat('
clones <- sim_cohort %>%
  mutate(
    limit = as.numeric(1.0),
    cloneid = paste0("clone2.0_", id)
  )')
sink()
source("routput/clonecode.txt")


#2) artificially censor data (keeping first observation that is censored)
# create output for slides
sink("routput/censcode.txt", split=FALSE)
cat('
cens_data <- clones %>%
  group_by(cloneid) %>%
  mutate(
    cens = pmin(1,cumsum(x > limit)),
    drop = pmin(2, cumsum(cens))
  ) %>%
  ungroup() %>%
  filter(
    drop < 2
  )
')
sink()
source("routput/censcode.txt")

# 3) create 1/0 weights to use for confounding/censoring during follow-up (alternative is to create new datasets with some observations dropped
#  note: these are not the inverse probability of censoring weights that will be used later
sink("routput/censmungecode.txt", split=FALSE)
cens_data <- group_by(cens_data, cloneid) %>%
  mutate(
    one = 1,
    nobs = cumsum(one),
    fobs____ = as.numeric(nobs == 1)
  ) %>%
  mutate(
    conf_weight = as.numeric(nobs  == 1), # all of these are at work
    fu_weight = as.numeric((nobs  > 1) & (atwork == 1)),    
    dconf = 0,
    dcens = 0,
    intervened = 0,
  ) %>%
  select(-c(one, nobs)) %>%
  ungroup()
sink()
source("routput/censmungecode.txt")



#4) fit censoring models (can include pre-baseline exposure in practice)
# use flexible functions of age and time based on restricted cubic splines
agekn0 = attr(rcspline.eval(filter(cens_data, time==1 & atwork == 1)$age, nk = 4), "knots")
yearkn0 = attr(rcspline.eval(filter(cens_data, time==1 & atwork == 1)$year, nk = 4), "knots")
agekn = attr(rcspline.eval(filter(cens_data, time>1 & atwork == 1)$age, nk = 4), "knots")
yearkn = attr(rcspline.eval(filter(cens_data, time>1 & atwork == 1)$year, nk = 4), "knots")
cxkn = attr(rcspline.eval(filter(cens_data, time>1 & atwork == 1)$cxl, nk = 4), "knots")
clkn = attr(rcspline.eval(filter(cens_data, time>1 & atwork == 1)$cumatworkl, nk = 4), "knots")

# fit models
sink("routput/censmodcode.txt", split=FALSE)
cat('# "censored at baseline" model 
confdmod <- glm(cens~ 
  year + rcspline.eval(year, knots=yearkn0) + 
  age + rcspline.eval(age, knots=agekn0) + 
  wagestatus + male + race, 
  data = cens_data, weight=conf_weight, family=binomial())

# "censored during follow-up" model 
censdmod <- glm(cens~ cxl + cumatworkl + 
  age + rcspline.eval(age, knots=agekn) + 
  wagestatus + male + race , data = cens_data, weight=fu_weight, 
  family=binomial())
')
sink()
# numerator for stabilizing weight
confnmod <- glm(cens~rcspline.eval(age, knots=agekn0) , data = cens_data, weight=conf_weight, family=binomial())
# numerator for stabilizing weight
censnmod <- glm(cens ~ rcspline.eval(age, knots=agekn), data = cens_data, weight=fu_weight, family=binomial())


sink("routput/censmodoutput.txt", split=TRUE)
source("routput/censmodcode.txt")
summary(confdmod) %>% print
summary(censdmod) %>% print
sink()

# 5) create weights
sink("routput/weightcode.txt", split=FALSE)
cat('cens_data$dconf = cens_data$conf_weight*as.numeric(predict(confdmod, type="response"))
cens_data$dcens = cens_data$fu_weight*as.numeric(predict(censdmod, type="response")) 
cens_data$nconf = cens_data$conf_weight*as.numeric(predict(confnmod, type="response"))
cens_data$ncens = cens_data$fu_weight*as.numeric(predict(censnmod, type="response"))   
cens_data <- cens_data %>% 
  mutate(
    wtcontr = case_when(
      ((fobs____ == 1) & (atwork==1)) ~ (1-cens)*(1-nconf)/(1-dconf),
      ((fobs____ == 0) & (atwork==1)) ~ (1-cens)*(1-ncens)/(1-dcens),
      .default=1
    )
  ) %>%
  group_by(cloneid) %>%
  mutate(
    ipw = cumprod(wtcontr),
    ipwt = pmin(10, cumprod(wtcontr)) # sometimes truncated weights are advocated
  )
')
sink()
source("routput/weightcode.txt")

# 5) examine weights
# crude histogram
hist(cens_data$ipw)
quantile(filter(cens_data, ipw>0)$ipw)
# top 5 weights
rnks = rank(-cens_data$ipw, ties.method="first")
highest = which(rnks< 5)
cens_data$ipw[highest]
cens_data[highest,c("id","x", "limit", "wagestatus", "male", "race", "age", "year", "cumatwork", "cumx", "wtcontr", "dconf", "dcens", "ipw", "ipwt")]
# id with max weight
whichclone = cens_data$cloneid[which(rnks== 1)]
allids = which(cens_data$cloneid == whichclone)
sink("routput/highestweight.txt", split=TRUE)
print(cens_data[allids,c("cloneid","x", "limit", "wagestatus", "male", "race", "age", "year", "cumatwork", "cumx", "wtcontr", "dconf", "dcens", "ipw", "ipwt")], n=100)
sink()
#6) estimate effect of introducing an exposure limit
sink("routput/weightedcoxcode.txt", split=FALSE)
cat('sim_cohort$exposed = 0
cens_data$exposed = 1
sim_cohort$ipw = 1
sim_cohort$cloneid = paste0("cloneobs_", sim_cohort$id)
combined_data <- bind_rows(sim_cohort, cens_data)

coxph(Surv(agein, age, d1)~exposed, 
     data=filter(combined_data, ipw != 0), weight=ipw, 
     id=cloneid, cluster=id) %>% 
      summary %>% print
')
sink()
source("routput/weightedcoxcode.txt")


# can also fit model for D2
coxph(Surv(agein, age, d2)~exposed, 
      data=filter(combined_data, ipw != 0), weight=ipw, 
      id=cloneid, cluster=id) %>% 
  summary %>% print

######### BOOTSTRAPPING #########
# We can optionally do bootstrapping for confidence intervals, which will generally be better
# than robust confidence intervals except in very large sample sizes
# Code to do this is not explicitly done here but it is a straightforward set of steps:
# 1. sample, with replacement, 10,000 individuals from the population (original size = 10,000)
# 2. Carry out clone-censor-weight + estimate effects on the sample from #1
# 3. Repeat 1 and 2 200+ times, recording effect estimate for each sample
# 4. Bootstrap standard error is the standard deviation of the effect estimates across the 200+ bootstrap samples


while(sink.number()>0){
  sink(NULL)
}
