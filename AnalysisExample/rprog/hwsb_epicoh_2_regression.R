setwd("~/Documents/Talks/Workshops/2024_EPICOH/AnalysisExample")
source("rprog/hwsb_epicoh_1_datamgmt.R")
library(dplyr) # data management functions from tidyverse
library(survival)

# Programming note: In many steps, I am writing code to a file and then running that file, 
#  this creates output for slides while also running the code without copy/paste errors. 
# It's not suggested as a general programming practice

#1) Assessing components
# is exposure associated with leaving work, given past exposure, confounders

# 1A) is employment associated with mortality, given exposure
sink("routput/mort1cumwork.txt", split=TRUE)
coxph(Surv(agein, age, d1)~ cumatwork + cxl + wagestatus + male + year, data = sim_cohort)
#The hazard of d1 is (1-exp(beta_1))% lower for each additional year of employment
sink()
sink("routput/mort2cumwork.txt", split=TRUE)
coxph(Surv(agein, age, d2)~ cumatwork + cxl + wagestatus + male + year, data = sim_cohort)
#The hazard of d2 is (1-exp(beta_1))% lower for each additional year of employment
sink()

# we could also just look at active vs. inactive work
sink("routput/mort1work.txt", split=TRUE)
coxph(Surv(agein, age, d1)~ atwork + cumatworkl + cxl + wagestatus + male + year, data = sim_cohort)
sink()
sink("routput/mort2work.txt", split=TRUE)
coxph(Surv(agein, age, d2)~ atwork + cumatworkl + cxl + wagestatus + male + year, data = sim_cohort)
sink()



# 1B) Is exposure associated with leaving work
# cumulative exposure
sink("routput/model_leavework_cum.txt", split=TRUE)
coxph(Surv(agein, age, leftwork)~ cxl + wagestatus + male + year, data = filter(sim_cohort, atwork == 1 | leftwork == 1))
sink()
# current exposure 
sink("routput/model_leavework_recent.txt", split=TRUE)
coxph(Surv(agein, age, leftwork)~ xl + cxl2 + wagestatus + male + year, data = filter(sim_cohort, atwork == 1 | leftwork == 1))
sink()


while(sink.number()>0){
  sink(NULL)
}

