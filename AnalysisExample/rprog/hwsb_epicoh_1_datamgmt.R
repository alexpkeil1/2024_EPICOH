setwd("~/Documents/Talks/Workshops/2024_EPICOH/AnalysisExample")

library(dplyr) # data management functions from tidyverse



#0) read in data)
sim_cohort = read.csv("data/sim_cohort.csv", stringsAsFactors = FALSE)

# 0b) first 10 observations in original data
sink("routput/first10obs.txt", split=FALSE)
print(as.data.frame(head(sim_cohort, n=15)), row.names = FALSE)
sink()

#1) data process (lagging variables)

# create output for slides
sink("routput/mungecode.txt", split=FALSE)
cat("sim_cohort <- sim_cohort %>%
  group_by(id) %>%
  mutate(
    one = 1,
    time = cumsum(one),
    timein = time-1,
    agein = age-1,
    xl = lag(x, default=0),
    cxl = lag(cumx, default=0),
    cumatwork = cumsum(atwork),
    cumatworkl = lag(cumatwork, default=0),
    atworkl = lag(atwork, default=0),  
    leftwork = as.numeric(atworkl==1 & atwork == 0),
    cxl2 = lag(cxl, default=0),  
  ) %>%
  select(-one) %>%
  ungroup()
sim_cohort <- sim_cohort %>%
  group_by(id) %>%
  mutate(
  male = as.numeric(gender == 'M')
  ) %>%
  ungroup()
")
sink()
source("routput/mungecode.txt")


# 1b) Table 1
bldat <- sim_cohort %>%
  filter(time == 1)
enddat <- sim_cohort %>%
  group_by(id) %>%
  filter(time == max(time)) %>%
  ungroup()

meansd <- function(x, digits=3) c(Mean=as.character(signif(mean(x),digits)), SD=as.character(signif(sd(x), digits)))
sumprop <- function(x, digits=2) c(N=as.character(sum(x)), Pct=as.character(signif(100*mean(x), digits)))
means = as.data.frame(with(bldat, rbind("Age"=meansd(agein), "Calendar year"=meansd(year-1, digits=4))))
means2 = as.data.frame(with(enddat, rbind(as.data.frame(with(sim_cohort, 
                                                             rbind(
                                                               "Years of follow-up"=meansd(time),
                                                               "Average X (at work)"=meansd(cumx/cumatwork), "Cumulative X"=meansd(cumx),
                                                               "Years worked"=meansd(cumatwork)
                                                             ))))))
props =  as.data.frame(with(bldat, rbind(
  "Waged"=sumprop(wagestatus==1), "Salaried"=sumprop(wagestatus==0),
  "Gender F"=sumprop(male==0), "Gender M"=sumprop(male==1), "White race"=sumprop(race=="W"), "Non-white race"=sumprop(race=="N")
)))
props2 =  as.data.frame(with(enddat, rbind("D1"=sumprop(d1), "D2"=sumprop(d2), "Survived"=sumprop(1-d1-d2))))
means$Variable = rownames(means)
means2$Variable = rownames(means2)
props$Variable = rownames(props)
props2$Variable = rownames(props2)

dat = merge(merge(merge(means, props, all=TRUE, on="Variable", sort=FALSE), means2, all=TRUE, on="Variable", sort=FALSE), props2, all=TRUE, on="Variable", sort=FALSE)
cf <- format(dat) ## use format to set other options like digits, justify , ...
cf[is.na(dat)] <- ""

sink("routput/table1.txt", split=FALSE)
print(cf, row.names = FALSE)
sink()

pdf("routput/ann_exposure_hist.pdf", width=3, height=4)
hist(filter(sim_cohort, atwork==1)$x, main=NULL, breaks=50, xlab="Annual exposure", ylab="Person years")
dev.off()


while(sink.number()>0){
  sink(NULL)
}
