## Analyses for non-adjacent dependency learning paper (Zettersten, Potter & Saffran)
library(tidyverse) # version 1.2.1
library(cowplot) # version 1.0.0
theme_set(theme_cowplot())
library(formattable) # version 0.2.0.1
library(lme4) # version 1.1-21
library(psych) # version 1.8.12
library(stats) # version 3.6.1
library(sciplot) #version 1.1-1
library(car) #version 3.0-3

#### read in complete data set ####
d <- read.table("APG_data.txt", header=T)

#add columns for hits and false alarms
#used for traditional signal detection analysis
d$hits = ifelse(d$correctResponse=="y",d$isRight,NA)
d$falseAlarms = ifelse(d$correctResponse=="n",1-d$isRight,NA)

######################################
#### demographics & overall means ####
######################################

#subject means
#computation of adjusted hit rates and false alarm rates uses Macmillan & Kaplan (1985) adjustment
subj_overall <- d %>%
  group_by(participantCode, exp, condition) %>%
  summarize(age=age[1],
            gender = gender[1],
            nativeLangEng = nativeLangEng[1],
            acc=mean(isRight),
            rt=mean(rt[rt<5000]),
            hitRate=mean(hits,na.rm=T),
            falseAlarmRate=mean(falseAlarms,na.rm=T)) %>%
  mutate(
    hitRateAdjusted= case_when(
      hitRate==1 ~ 1 - 1/(2*18),
      hitRate==0 ~ 1/(2*18),
      TRUE ~ hitRate
    ),
    falseAlarmRateAdjusted= case_when(
      falseAlarmRate==0 ~ 1/(2*18),
      falseAlarmRate==1 ~ 1 - 1/(2*18),
      TRUE ~ falseAlarmRate
    )
  ) %>%
  mutate(
    dprime=qnorm(hitRateAdjusted) - qnorm(falseAlarmRateAdjusted),
    c=-.5*(qnorm(hitRateAdjusted) + qnorm(falseAlarmRateAdjusted)))

#total demographics
demo <- subj_overall %>%
  group_by(exp) %>%
  summarize(
    N_learnable=sum(condition=="Learnable Pre-Exposure"),
    N_nonLearnable=sum(condition=="Non-Learnable Pre-Exposure"),
    N_noPre =sum(condition=="No Pre-Exposure"),
    N_unstructured = sum(condition=="Unstructured Pre-Exposure"),
    N_LearnablePreExpOnly = sum(condition=="Learnable Pre-Exposure Only"),
    N_female = sum(gender=="Female"),
    N_nativeLangEng=sum(nativeLangEng==1),
    mean_age = mean(age,na.rm=T),
    sd_age = sd(age))

#split accuracy by test type
subj_testType <- d %>%
  group_by(participantCode,exp,testType, condition) %>%
  summarize(
    acc=mean(isRight),
    rt=mean(rt[rt<5000]),
    hitRate=mean(hits,na.rm=T),
    falseAlarmRate=mean(falseAlarms,na.rm=T)) %>%
  mutate(
    hitRateAdjusted= case_when(
      hitRate==1 ~ 1 - 1/(2*9),
      hitRate==0 ~ 1/(2*9),
      TRUE ~ hitRate
    ),
    falseAlarmRateAdjusted= case_when(
      falseAlarmRate==0 ~ 1/(2*9),
      falseAlarmRate==1 ~ 1 - 1/(2*9),
      TRUE ~ falseAlarmRate
    )
  ) %>%
  mutate(
    dprime=qnorm(hitRateAdjusted) - qnorm(falseAlarmRateAdjusted),
    c=-.5*(qnorm(hitRateAdjusted) + qnorm(falseAlarmRateAdjusted)))

#create a wide format data set to correlate Familiar X and Novel X accuracy
subj_accuracy_wide <- spread(subset(subj_testType, select=-c(rt,dprime,c,hitRate,hitRateAdjusted,falseAlarmRate, falseAlarmRateAdjusted)), testType, acc)
  
#all experiment conditions - average accuracy, dprime, response bias
overall_exp <- subj_overall %>%
  group_by(exp,condition) %>%
  summarize(
    N=sum(!is.na(acc)),
    accuracy=mean(acc,na.rm=T),
    se=se(acc,na.rm=T),
    accuracy_ci=qt(0.975, N-1)*sd(acc,na.rm=T)/sqrt(N),
    accuracy_lower_ci=accuracy-accuracy_ci,
    accuracy_upper_ci=accuracy+accuracy_ci,
    d_prime=mean(dprime,na.rm=T),
    d_prime_ci=qt(0.975, N-1)*sd(dprime,na.rm=T)/sqrt(N),
    d_prime_lower_ci=d_prime-d_prime_ci,
    d_prime_upper_ci=d_prime+d_prime_ci,
    c_bias=mean(c,na.rm=T),
    c_bias_ci=qt(0.975, N-1)*sd(c,na.rm=T)/sqrt(N),
    c_bias_lower_ci=c_bias-c_bias_ci,
    c_bias_upper_ci=c_bias+c_bias_ci
  )

#split by test type - average accuracy, dprime, response bias
test_type_exp <- subj_testType %>%
  group_by(exp,condition,testType) %>%
  summarize(
    N=sum(!is.na(acc)),
    accuracy=mean(acc,na.rm=T),
    se=se(acc,na.rm=T),
    accuracy_ci=qt(0.975, N-1)*sd(acc,na.rm=T)/sqrt(N),
    accuracy_lower_ci=accuracy-accuracy_ci,
    accuracy_upper_ci=accuracy+accuracy_ci,
    d_prime=mean(dprime,na.rm=T),
    d_prime_ci=qt(0.975, N-1)*sd(dprime,na.rm=T)/sqrt(N),
    d_prime_lower_ci=d_prime-d_prime_ci,
    d_prime_upper_ci=d_prime+d_prime_ci,
    c_bias=mean(c,na.rm=T),
    c_bias_ci=qt(0.975, N-1)*sd(c,na.rm=T)/sqrt(N),
    c_bias_lower_ci=c_bias-c_bias_ci,
    c_bias_upper_ci=c_bias+c_bias_ci
  )

#collapsing across all experiments
#all conditions - average accuracy, dprime, response bias
overall <- subj_overall %>%
  group_by(condition) %>%
  summarize(
    N=sum(!is.na(acc)),
    se=se(acc,na.rm=T),
    accuracy=mean(acc,na.rm=T),
    accuracy_ci=qt(0.975, N-1)*sd(acc,na.rm=T)/sqrt(N),
    accuracy_lower_ci=accuracy-accuracy_ci,
    accuracy_upper_ci=accuracy+accuracy_ci,
    d_prime=mean(dprime,na.rm=T),
    d_prime_ci=qt(0.975, N-1)*sd(dprime,na.rm=T)/sqrt(N),
    d_prime_lower_ci=d_prime-d_prime_ci,
    d_prime_upper_ci=d_prime+d_prime_ci,
    c_bias=mean(c,na.rm=T),
    c_bias_ci=qt(0.975, N-1)*sd(c,na.rm=T)/sqrt(N),
    c_bias_lower_ci=c_bias-c_bias_ci,
    c_bias_upper_ci=c_bias+c_bias_ci
  )

#collapsing across all experiments
#split by test type - average accuracy, dprime, response bias
test_type <- subj_testType %>%
  group_by(condition,testType) %>%
  summarize(
    N=sum(!is.na(acc)),
    accuracy=mean(acc,na.rm=T),
    se=se(acc,na.rm=T),
    accuracy_ci=qt(0.975, N-1)*sd(acc,na.rm=T)/sqrt(N),
    accuracy_lower_ci=accuracy-accuracy_ci,
    accuracy_upper_ci=accuracy+accuracy_ci,
    d_prime=mean(dprime,na.rm=T),
    d_prime_ci=qt(0.975, N-1)*sd(dprime,na.rm=T)/sqrt(N),
    d_prime_lower_ci=d_prime-d_prime_ci,
    d_prime_upper_ci=d_prime+d_prime_ci,
    c_bias=mean(c,na.rm=T),
    c_bias_ci=qt(0.975, N-1)*sd(c,na.rm=T)/sqrt(N),
    c_bias_lower_ci=c_bias-c_bias_ci,
    c_bias_upper_ci=c_bias+c_bias_ci
  )

######################
#### Experiment 1 ####
######################

#### Experiment 1 - Logistic Mixed-Effects Analyses ####

#recode condition
d$conditionC <- ifelse(d$condition=="Learnable Pre-Exposure",0.5,
                    ifelse(d$condition=="Non-Learnable Pre-Exposure",-0.5,NA))
##all data
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp1"),family=binomial,glmerControl(optimizer="bobyqa"))
#removing the random intercept of stimulus which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+(1|participantCode)+(0+conditionC|stimulus),data=subset(d,exp=="exp1"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3:4,]
Anova(m, type="III") #Wald chi-squared test yields similar results to z-values across experiments

##by test type
#Familiar X Test
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp1"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa")) #convergence warning
#removing the random effect (by-item random intercept) with an estimated covariance of zero (identical results)
m <- glmer(isRight~conditionC+(1|participantCode)+(0+conditionC|stimulus),data=subset(d,exp=="exp1"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3:4,]
#Novel X Test
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp1"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
#removing the random intercept of stimulus which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+(1|participantCode)+(0+conditionC|stimulus),data=subset(d,exp=="exp1"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3:4,]

##comparison to chance
#Learnable Pre-Exposure
#Familiar X Test
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp1"&testType=="familiarX"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
#removing the random intercept of stimulus which leads to a singular model fit (identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp1"&testType=="familiarX"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#Novel X Test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp1"&testType=="novelX"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#Non-Learnable Pre-Exposure
#Familiar X Test
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp1"&testType=="familiarX"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
#removing the random intercept of stimulus which leads to a singular model fit (identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp1"&testType=="familiarX"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#Novel X Test
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp1"&testType=="novelX"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model (after removing random effects structure - nevertheless yields virtually identical results to maximal model)
m <- glm(isRight~1,data=subset(d,exp=="exp1"&testType=="novelX"&conditionC==-0.5),family=binomial)
summary(m)
confint(m, method="Wald")

## Interaction between test type and condition
#recode test type
d$testTypeC <- ifelse(d$testType=="novelX",0.5,
                      ifelse(d$testType=="familiarX",-0.5,NA))
#maximal model (does not converge)
m <- glmer(isRight~conditionC*testTypeC+(1+testTypeC|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp1"),family=binomial,glmerControl(optimizer="bobyqa"))
#final model with simplified random effects structure
m <- glmer(isRight~conditionC*testTypeC+(1|participantCode)+(0+conditionC|stimulus),data=subset(d,exp=="exp1"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3:6,]

#### Experiment 1 - Signal Detection Analyses ####

## Sensitivity (d-prime)
#recode condition
subj_overall$conditionC <- ifelse(subj_overall$condition=="Learnable Pre-Exposure",0.5,
                                  ifelse(subj_overall$condition=="Non-Learnable Pre-Exposure",-0.5,NA))
##all data
m <- lm(dprime~conditionC,data=subset(subj_overall,exp=="exp1"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp1"&conditionC==0.5))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp1"&conditionC==-0.5))
summary(m)

##Familiar X trials
#recode condition
subj_testType$conditionC <- ifelse(subj_testType$condition=="Learnable Pre-Exposure",0.5,
                                   ifelse(subj_testType$condition=="Non-Learnable Pre-Exposure",-0.5,NA))
m <- lm(dprime~conditionC,data=subset(subj_testType,exp=="exp1"&testType=="familiarX"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp1"&testType=="familiarX"&conditionC==0.5))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp1"&testType=="familiarX"&conditionC==-0.5))
summary(m)

##Novel X trials
m <- lm(dprime~conditionC,data=subset(subj_testType,exp=="exp1"&testType=="novelX"))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp1"&testType=="novelX"&conditionC==0.5))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp1"&testType=="novelX"&conditionC==-0.5))
summary(m)

# Interaction between Novel X and Familiar X trials
subj_testType$testTypeC <- ifelse(subj_testType$testType=="novelX",0.5,
                                  ifelse(subj_testType$testType=="familiarX",-0.5,NA))
m <- lm(dprime~conditionC*testTypeC,data=subset(subj_testType,exp=="exp1"))
summary(m)

## Response Bias (c)

##all data
m <- lm(c~conditionC,data=subset(subj_overall,exp=="exp1"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp1"&conditionC==0.5))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp1"&conditionC==-0.5))
summary(m)

##Familiar X trials
m <- lm(c~conditionC,data=subset(subj_testType,exp=="exp1"&testType=="familiarX"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp1"&testType=="familiarX"&conditionC==0.5))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp1"&testType=="familiarX"&conditionC==-0.5))
summary(m)

##Novel X trials
m <- lm(c~conditionC,data=subset(subj_testType,exp=="exp1"&testType=="novelX"))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp1"&testType=="novelX"&conditionC==0.5))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp1"&testType=="novelX"&conditionC==-0.5))
summary(m)

# Interaction between Novel X and Familiar X trials
m <- lm(c~conditionC*testTypeC,data=subset(subj_testType,exp=="exp1"))
summary(m)

#### Experiment 1 - Correlation Familiar X - Novel X ####

##Correlation
c <- corr.test(subset(subj_accuracy_wide, condition=="Learnable Pre-Exposure" & exp=="exp1")[,c("novelX","familiarX")])
c
c$p
c <- corr.test(subset(subj_accuracy_wide, condition=="Non-Learnable Pre-Exposure" & exp=="exp1")[,c("novelX","familiarX")])
c
c$p
#check outlier in Non-Learnable Pre-Exposure Condition
outlierTest(lm(novelX~familiarX,subset(subj_accuracy_wide,exp=="exp1"&condition=="Non-Learnable Pre-Exposure")))
#correlation without outlier
c <- corr.test(subset(subj_accuracy_wide, condition=="Non-Learnable Pre-Exposure" & exp=="exp1"&participantCode!="apg_exp1_p7")[,c("novelX","familiarX")])
c
c$p

#interaction
m <- lm(novelX~familiarX*condition, subset(subj_accuracy_wide,exp=="exp1"))
summary(m)

######################
#### Experiment 2 ####
######################

#### Experiment 2 - Logistic Mixed-Effects Analyses ####

#recode condition
d$conditionC <- ifelse(d$condition=="Learnable Pre-Exposure",0.5,
                    ifelse(d$condition=="Non-Learnable Pre-Exposure",-0.5,
                           ifelse(d$condition=="No Pre-Exposure",0,NA)))
d$conditionOrthContrast2 <- ifelse(d$condition=="Learnable Pre-Exposure",-1/3,
                                ifelse(d$condition=="Non-Learnable Pre-Exposure",-1/3,
                                       ifelse(d$condition=="No Pre-Exposure",2/3,NA)))

##Richter single-contrast approach
#overall analysis (Richter single-contrast apporach)
#model with maximal random effects structure (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - after removing the random effects of stimulus which lead to a singular model fit - virtually identical results)
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,exp=="exp2"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##check same analysis with Abelson & Prentice approach
#all data
#maximal  model (does not converge)
m=glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode)+(1+conditionC+conditionOrthContrast2|stimulus),data=subset(d,exp=="exp2"),family=binomial) #convergence warning
#final model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode),data=subset(d,exp=="exp2"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)

#compare conditions
#Learnable vs. Non-Learnable
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"&conditionC!=0),family=binomial,glmerControl(optimizer="bobyqa"))
#final model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,exp=="exp2"&conditionC!=0),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#storing p-values
p <- c(summary(m)$coefficients[2,4])

#Non-Learnable vs. No Pre-Exposure
d$conditionNonLearnableNoPreExp <- ifelse(d$conditionC==0.5,1.5,
                                    ifelse(d$conditionC==0,0.5,d$conditionC))
#maximal model (does not converge)
m <- glmer(isRight~conditionNonLearnableNoPreExp+(1|participantCode)+(1+conditionNonLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&conditionNonLearnableNoPreExp!=1.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionNonLearnableNoPreExp+(1|participantCode),data=subset(d,exp=="exp2"&conditionNonLearnableNoPreExp!=1.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#storing p-values
p <- c(p,summary(m)$coefficients[2,4])

#No Pre-Exposure vs. Learnable
d$conditionLearnableNoPreExp <- ifelse(d$conditionC==-0.5,-1.5,
                                    ifelse(d$conditionC==0,-0.5,d$conditionC))
#maximal model (does not converge)
m <- glmer(isRight~conditionLearnableNoPreExp+(1|participantCode)+(1+conditionLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&conditionLearnableNoPreExp!=-1.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final model - removing the by-stimulus random intercept which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionLearnableNoPreExp+(1|participantCode)+(0+conditionLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&conditionLearnableNoPreExp!=-1.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#storing p-values
p <- c(p,summary(m)$coefficients[2,4])
#if using Holm-Bonferroni correction:
p.adjust(p)

#Compare each condition to chance performance (collapsing across test trials)
#Learnable
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#No Pre-Exposure
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&conditionC==0),family=binomial,glmerControl(optimizer="bobyqa"))
#final model - removing the by-stimulus random intercept which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp2"&conditionC==0),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#non-Learnable
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final model - removing the by-stimulus random intercept which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp2"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##### Familiar X Trials

## Condition Effect
#Richer single-contrast approach
#Familiar X Test
#model with maximal random effects structure (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - after removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,exp=="exp2"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##check same analysis with Abelson & Prentice approach
#Familiar X Test
#maximal  model (does not converge)
m <- glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode)+(1+conditionC+conditionOrthContrast2|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
#final model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode),data=subset(d,exp=="exp2"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)

##Pairwise comparisons
#Learnable vs. Non-Learnable Pre-Exposure
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"&conditionC!=0),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,exp=="exp2"&testType=="familiarX"&conditionC!=0),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[2:3,]

#storing p-values
p <- c(summary(m)$coefficients[2,4])

#No Pre-Exposure vs. Non-Learnable Pre-Exposure
d$conditionNonLearnableNoPreExp <- ifelse(d$conditionC==0.5,1.5,
                                          ifelse(d$conditionC==0,0.5,d$conditionC))
#maximal model (does not converge)
m <- glmer(isRight~conditionNonLearnableNoPreExp+(1|participantCode)+(1+conditionNonLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"&conditionNonLearnableNoPreExp!=1.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionNonLearnableNoPreExp+(1|participantCode),data=subset(d,exp=="exp2"&testType=="familiarX"&conditionNonLearnableNoPreExp!=1.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[2:3,]

#storing p-values
p <- c(p,summary(m)$coefficients[2,4])

#Learnable vs. No Pre-Exposure
d$conditionLearnableNoPreExp <- ifelse(d$conditionC==-0.5,-1.5,
                                       ifelse(d$conditionC==0,-0.5,d$conditionC))
#maximal model
m <- glmer(isRight~conditionLearnableNoPreExp+(1|participantCode)+(1+conditionLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"&conditionLearnableNoPreExp!=-1.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing the by-stimulus random intercept which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionLearnableNoPreExp+(1|participantCode)+(0+conditionLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"&conditionLearnableNoPreExp!=-1.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3:4,]

#storing p-values
p <- c(p,summary(m)$coefficients[2,4])
#if using Holm-Bonferroni correction:
p.adjust(p)

#Testing against chance: Learnable Pre-Exposure
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"&condition=="Learnable Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3,]

#Testing against chance: No Pre-Exposure
#maximal model
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"&condition=="No Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing the by-stimulus random intercept which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp2"&testType=="familiarX"&condition=="No Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[2,]

#Testing against chance: Non-Learnable Pre-Exposure
#maximal model
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&testType=="familiarX"&condition=="Non-Learnable Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing random effects (virtually identical results)
m <- glm(isRight~1,data=subset(d,exp=="exp2"&testType=="familiarX"&condition=="Non-Learnable Pre-Exposure"),family=binomial)
summary(m)
confint(m, method="Wald")

##### Novel X Trials

##Condition Effect
#Richer single-contrast approach
#Novel X Test
#model with maximal random effects structure (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
#final model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,exp=="exp2"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##check same analysis with Abelson & Prentice approach
#Novel X Test
#maximal  model (does not converge)
m <- glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode)+(1+conditionC+conditionOrthContrast2|stimulus),data=subset(d,exp=="exp2"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
#final model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode),data=subset(d,exp=="exp2"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)

##Pairwise comparisons
#Learnable vs. Non-Learnable Pre-Exposure
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"&testType=="novelX"&conditionC!=0),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,exp=="exp2"&testType=="novelX"&conditionC!=0),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[2:3,]

#storing p-values
p <- c(summary(m)$coefficients[2,4])

#No Pre-Exposure vs. Non-Learnable Pre-Exposure
d$conditionNonLearnableNoPreExp <- ifelse(d$conditionC==0.5,1.5,
                                          ifelse(d$conditionC==0,0.5,d$conditionC))
#maximal model (does not converge)
m <- glmer(isRight~conditionNonLearnableNoPreExp+(1|participantCode)+(1+conditionNonLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&testType=="novelX"&conditionNonLearnableNoPreExp!=1.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing the random effects of stimulus which lead to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionNonLearnableNoPreExp+(1|participantCode),data=subset(d,exp=="exp2"&testType=="novelX"&conditionNonLearnableNoPreExp!=1.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[2:3,]

#storing p-values
p <- c(p,summary(m)$coefficients[2,4])

#Learnable vs. No Pre-Exposure
d$conditionLearnableNoPreExp <- ifelse(d$conditionC==-0.5,-1.5,
                                       ifelse(d$conditionC==0,-0.5,d$conditionC))
#maximal model
m <- glmer(isRight~conditionLearnableNoPreExp+(1|participantCode)+(1+conditionLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&testType=="novelX"&conditionLearnableNoPreExp!=-1.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing the by-stimulus random intercept which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~conditionLearnableNoPreExp+(1|participantCode),data=subset(d,exp=="exp2"&testType=="novelX"&conditionLearnableNoPreExp!=-1.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[2:3,]

#storing p-values
p <- c(p,summary(m)$coefficients[2,4])
#if using Holm-Bonferroni correction:
p.adjust(p)

#Testing against chance: Learnable Pre-Exposure
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&testType=="novelX"&condition=="Learnable Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3,]

#Testing against chance: No Pre-Exposure
#maximal model
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&testType=="novelX"&condition=="No Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing the by-stimulus random intercept which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp2"&testType=="novelX"&condition=="No Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[2,]

#Testing against chance: Non-Learnable Pre-Exposure
#maximal model
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&testType=="novelX"&condition=="Non-Learnable Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model - removing random effect for stimulus (virtually identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp2"&testType=="novelX"&condition=="Non-Learnable Pre-Exposure"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#### Interaction between test type and condition
#recode test type
d$testTypeC <- ifelse(d$testType=="novelX",0.5,
                      ifelse(d$testType=="familiarX",-0.5,NA))
#maximal model (does not converge)
m <- glmer(isRight~conditionC*testTypeC+(1+testTypeC|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"),family=binomial,glmerControl(optimizer="bobyqa"))
##final converging model
m <- glmer(isRight~conditionC*testTypeC+(1|participantCode),data=subset(d,exp=="exp2"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#### Experiment 2 - Signal Detection Analysis ####

#### Sensitivity

#### All Test Trials

## Condition Effect

#recode condition
subj_overall$conditionC <- ifelse(subj_overall$condition=="Learnable Pre-Exposure",0.5,
                                  ifelse(subj_overall$condition=="Non-Learnable Pre-Exposure",-0.5,
                                         ifelse(subj_overall$condition=="No Pre-Exposure",0,NA)))
subj_overall$conditionOrthContrast2 <- ifelse(subj_overall$condition=="Learnable Pre-Exposure",-1/3,
                                              ifelse(subj_overall$condition=="Non-Learnable Pre-Exposure",-1/3,
                                                     ifelse(subj_overall$condition=="No Pre-Exposure",2/3,NA)))
#overall analysis (Richter single-contrast apporach)
m <- lm(dprime~conditionC,data=subset(subj_overall,exp=="exp2"))
summary(m)

##check same analysis with Abelson & Prentice approach
m <- lm(dprime~conditionC+conditionOrthContrast2,data=subset(subj_overall,exp=="exp2"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp2"&conditionC==0.5))
summary(m)

#Testing against chance: No Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp2"&conditionC==0))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp2"&conditionC==-0.5))
summary(m)

## Pairwise condition comparisons
#Learnable vs. Non-Learnable
m <- lm(dprime~conditionC,data=subset(subj_overall,exp=="exp2"&conditionC!=0))
summary(m)

#Non-Learnable vs. No Pre-Exposure
subj_overall$conditionNonLearnableNoPreExp <- ifelse(subj_overall$conditionC==0.5,1.5,
                                                     ifelse(subj_overall$conditionC==0,0.5,subj_overall$conditionC))
m <- lm(dprime~conditionNonLearnableNoPreExp,data=subset(subj_overall,exp=="exp2"&conditionNonLearnableNoPreExp!=1.5))
summary(m)

#No Pre-Exposure vs. Learnable
subj_overall$conditionLearnableNoPreExp <- ifelse(subj_overall$conditionC==-0.5,-1.5,
                                                  ifelse(subj_overall$conditionC==0,-0.5,subj_overall$conditionC))
m <- lm(dprime~conditionLearnableNoPreExp,data=subset(subj_overall,exp=="exp2"&conditionLearnableNoPreExp!=-1.5))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp2"&conditionC==0.5))
summary(m)

#Testing against chance: No Pre-Exposure
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp2"&conditionC==0))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp2"&conditionC==-0.5))
summary(m)

##### Familiar X Trials

##Condition Effect

subj_testType$conditionC <- ifelse(subj_testType$condition=="Learnable Pre-Exposure",0.5,
                                   ifelse(subj_testType$condition=="Non-Learnable Pre-Exposure",-0.5,
                                          ifelse(subj_testType$condition=="No Pre-Exposure",0,NA)))
subj_testType$conditionOrthContrast2 <- ifelse(subj_testType$condition=="Learnable Pre-Exposure",-1/3,
                                               ifelse(subj_testType$condition=="Non-Learnable Pre-Exposure",-1/3,
                                                      ifelse(subj_testType$condition=="No Pre-Exposure",2/3,NA)))

#overall analysis (Richter single-contrast apporach)
m <- lm(dprime~conditionC,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"))
summary(m)

##check same analysis with Abelson & Prentice approach
m <- lm(dprime~conditionC+conditionOrthContrast2,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"))
summary(m)

##Pairwise comparisons
#compare conditions
#Learnable vs. Non-Learnable
m <- lm(dprime~conditionC,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionC!=0))
summary(m)

#Non-Learnable vs. No Pre-Exposure
subj_testType$conditionNonLearnableNoPreExp <- ifelse(subj_testType$conditionC==0.5,1.5,
                                                      ifelse(subj_testType$conditionC==0,0.5,subj_testType$conditionC))
m <- lm(dprime~conditionNonLearnableNoPreExp,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionNonLearnableNoPreExp!=1.5))
summary(m)

#No Pre-Exposure vs. Learnable
subj_testType$conditionLearnableNoPreExp <- ifelse(subj_testType$conditionC==-0.5,-1.5,
                                                   ifelse(subj_testType$conditionC==0,-0.5,subj_testType$conditionC))
m <- lm(dprime~conditionLearnableNoPreExp,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionLearnableNoPreExp!=-1.5))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionC==0.5))
summary(m)

#Testing against chance: No Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionC==0))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionC==-0.5))
summary(m)

##### Novel X Trials
#Condition Effect
subj_testType$conditionC <- ifelse(subj_testType$condition=="Learnable Pre-Exposure",0.5,
                                   ifelse(subj_testType$condition=="Non-Learnable Pre-Exposure",-0.5,
                                          ifelse(subj_testType$condition=="No Pre-Exposure",0,NA)))
subj_testType$conditionOrthContrast2 <- ifelse(subj_testType$condition=="Learnable Pre-Exposure",-1/3,
                                               ifelse(subj_testType$condition=="Non-Learnable Pre-Exposure",-1/3,
                                                      ifelse(subj_testType$condition=="No Pre-Exposure",2/3,NA)))

#Novel X test
#overall analysis (Richter single-contrast apporach)
m <- lm(dprime~conditionC,data=subset(subj_testType,exp=="exp2"&testType=="novelX"))
summary(m)

##check same analysis with Abelson & Prentice approach
m <- lm(dprime~conditionC+conditionOrthContrast2,data=subset(subj_testType,exp=="exp2"&testType=="novelX"))
summary(m)

##Pairwise comparisons
#compare conditions
#Learnable vs. Non-Learnable
m <- lm(dprime~conditionC,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionC!=0))
summary(m)

#Non-Learnable vs. No Pre-Exposure
subj_testType$conditionNonLearnableNoPreExp <- ifelse(subj_testType$conditionC==0.5,1.5,
                                                      ifelse(subj_testType$conditionC==0,0.5,subj_testType$conditionC))
m <- lm(dprime~conditionNonLearnableNoPreExp,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionNonLearnableNoPreExp!=1.5))
summary(m)

#No Pre-Exposure vs. Learnable
subj_testType$conditionLearnableNoPreExp <- ifelse(subj_testType$conditionC==-0.5,-1.5,
                                                   ifelse(subj_testType$conditionC==0,-0.5,subj_testType$conditionC))
m <- lm(dprime~conditionLearnableNoPreExp,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionLearnableNoPreExp!=-1.5))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionC==0.5))
summary(m)

#Testing against chance: No Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionC==0))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionC==-0.5))
summary(m)

## Interaction between Novel X and Familiar X trials
#recode test type
subj_testType$testTypeC <- ifelse(subj_testType$testType=="novelX",0.5,
                                  ifelse(subj_testType$testType=="familiarX",-0.5,NA))
m <- lm(dprime~conditionC*testTypeC,data=subset(subj_testType,exp=="exp2"))
summary(m)

#### Response Bias

##### All Test Trials
#Unlike with the sensitivty/ accuracy measures, we had no specific hypothesis (linear or otherwise) with respect to response bias. 
#Hence, we test for an overall effect of condition (multi-df, omnibus effect) with condition dummy-coded.
subj_testType$conditionC <- ifelse(subj_testType$condition=="Learnable Pre-Exposure",0.5,
                                   ifelse(subj_testType$condition=="Non-Learnable Pre-Exposure",-0.5,
                                          ifelse(subj_testType$condition=="No Pre-Exposure",0,NA)))

#overall analysis (condition dummy-coded, since no directional hypothesis for response bias)
m <- lm(c~condition,data=subset(subj_overall,exp=="exp2"))
summary(m)
Anova(m, type="III")
#There is no significant overall effect of condition on response bias _c_ - indicating that response bias does not differ overall across conditions.

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&conditionC==0.5))
summary(m)

#Testing against chance: No Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&conditionC==0))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&conditionC==-0.5))
summary(m)

##### Pairwise condition comparisons
#Learnable vs. Non-Learnable
m <- lm(c~conditionC,data=subset(subj_overall,exp=="exp2"&conditionC!=0))
summary(m)

#Non-Learnable vs. No Pre-Exposure
subj_overall$conditionNonLearnableNoPreExp <- ifelse(subj_overall$conditionC==0.5,1.5,
                                                     ifelse(subj_overall$conditionC==0,0.5,subj_overall$conditionC))
m <- lm(c~conditionNonLearnableNoPreExp,data=subset(subj_overall,exp=="exp2"&conditionNonLearnableNoPreExp!=1.5))
summary(m)

#No Pre-Exposure vs. Learnable
subj_overall$conditionLearnableNoPreExp <- ifelse(subj_overall$conditionC==-0.5,-1.5,
                                                  ifelse(subj_overall$conditionC==0,-0.5,subj_overall$conditionC))
m <- lm(c~conditionLearnableNoPreExp,data=subset(subj_overall,exp=="exp2"&conditionLearnableNoPreExp!=-1.5))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_overall,exp=="exp2"&conditionC==0.5))
summary(m)

#Testing against chance: No Pre-Exposure
m <- lm(c~1,data=subset(subj_overall,exp=="exp2"&conditionC==0))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_overall,exp=="exp2"&conditionC==-0.5))
summary(m)

##### Familiar X Trials
##Condition Effect
#overall analysis (condition dummy-coded, since no directional hypothesis for response bias)
m <- lm(c~condition,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"))
summary(m)
Anova(m, type="III")

##Pairwise comparisons
#Learnable vs. Non-Learnable
m <- lm(c~conditionC,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionC!=0))
summary(m)

#Non-Learnable vs. No Pre-Exposure
subj_testType$conditionNonLearnableNoPreExp <- ifelse(subj_testType$conditionC==0.5,1.5,
                                                      ifelse(subj_testType$conditionC==0,0.5,subj_testType$conditionC))
m <- lm(c~conditionNonLearnableNoPreExp,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionNonLearnableNoPreExp!=1.5))
summary(m)

#No Pre-Exposure vs. Learnable
subj_testType$conditionLearnableNoPreExp <- ifelse(subj_testType$conditionC==-0.5,-1.5,
                                                   ifelse(subj_testType$conditionC==0,-0.5,subj_testType$conditionC))
m <- lm(c~conditionLearnableNoPreExp,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionLearnableNoPreExp!=-1.5))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionC==0.5))
summary(m)
#Testing against chance: No Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionC==0))
summary(m)
#Testing against chance: Non-Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&testType=="familiarX"&conditionC==-0.5))
summary(m)

####Novel X Trials
##Condition Effect
#overall analysis (condition dummy-coded, since no directional hypothesis for response bias)
m <- lm(c~condition,data=subset(subj_testType,exp=="exp2"&testType=="novelX"))
summary(m)
Anova(m, type="III")

##Pairwise comparisons
#Learnable vs. Non-Learnable
m <- lm(c~conditionC,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionC!=0))
summary(m)

#Non-Learnable vs. No Pre-Exposure
subj_testType$conditionNonLearnableNoPreExp <- ifelse(subj_testType$conditionC==0.5,1.5,
                                                      ifelse(subj_testType$conditionC==0,0.5,subj_testType$conditionC))
m <- lm(c~conditionNonLearnableNoPreExp,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionNonLearnableNoPreExp!=1.5))
summary(m)

#No Pre-Exposure vs. Learnable
subj_testType$conditionLearnableNoPreExp <- ifelse(subj_testType$conditionC==-0.5,-1.5,
                                                   ifelse(subj_testType$conditionC==0,-0.5,subj_testType$conditionC))
m <- lm(c~conditionLearnableNoPreExp,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionLearnableNoPreExp!=-1.5))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionC==0.5))
summary(m)

#Testing against chance: No Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionC==0))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp2"&testType=="novelX"&conditionC==-0.5))
summary(m)

##Interaction between Novel X and Familiar X trials
m <- lm(c~condition*testTypeC,data=subset(subj_testType,exp=="exp2"))
summary(m)
Anova(m, type="III")


#### Experiment 2 - Correlation Familiar  X - Novel X ####

##Correlations between Familiar X and Novel X Test
c <- corr.test(subset(subj_accuracy_wide, condition=="Learnable Pre-Exposure" & exp=="exp2")[,c("novelX","familiarX")])
c
c$p
c <- corr.test(subset(subj_accuracy_wide, condition=="No Pre-Exposure" & exp=="exp2")[,c("novelX","familiarX")])
c
c$p
c <- corr.test(subset(subj_accuracy_wide, condition=="Non-Learnable Pre-Exposure" & exp=="exp2")[,c("novelX","familiarX")])
c
c$p

#interaction
subj_accuracy_wide$conditionC=ifelse(subj_accuracy_wide$condition=="Learnable Pre-Exposure",0.5,
                           ifelse(subj_accuracy_wide$condition=="Non-Learnable Pre-Exposure",-0.5,
                                  ifelse(subj_accuracy_wide$condition=="No Pre-Exposure",0,NA)))
m <- lm(novelX~familiarX*conditionC, subset(subj_accuracy_wide,exp=="exp2"))
summary(m)

####Hits vs. False Alarm Rates Experiment 2
#split by test type and correct response type
subjType <- d %>%
  group_by(participantCode,exp,testType, correctResponse,condition) %>%
  summarize(acc=mean(isRight))

#test type
response_type_exp <- subjType %>%
  group_by(exp,condition,testType,correctResponse) %>%
  summarize(
    N=sum(!is.na(acc)),
    accuracy=mean(acc),
    accuracy_ci=qt(0.975, N-1)*sd(acc,na.rm=T)/sqrt(N),
    accuracy_lower_ci=accuracy-accuracy_ci,
    accuracy_upper_ci=accuracy+accuracy_ci,
  )

t.test(subjType$acc[subjType$correctResponse=="y"&subjType$exp=="exp2"&subjType$condition=="No Pre-Exposure"&subjType$testType=="familiarX"],subjType$acc[subjType$correctResponse=="y"&subjType$exp=="exp2"&(subjType$condition=="Learnable Pre-Exposure"|subjType$condition=="Non-Learnable Pre-Exposure")&subjType$testType=="familiarX"],var.equal = T)
t.test(subjType$acc[subjType$correctResponse=="n"&subjType$exp=="exp2"&subjType$condition=="Learnable Pre-Exposure"&subjType$testType=="familiarX"],subjType$acc[subjType$correctResponse=="n"&subjType$exp=="exp2"&(subjType$condition=="No Pre-Exposure"|subjType$condition=="Non-Learnable Pre-Exposure")&subjType$testType=="familiarX"],var.equal = T)

######################
#### Experiment 3 ####
######################

#### Experiment 3 - Logistic Mixed-Effects Analyses ####
#recode condition
d$conditionC <- ifelse(d$condition=="Learnable Pre-Exposure",0.5,
                       ifelse(d$condition=="Unstructured Pre-Exposure",-0.5,NA))
##all data
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp3"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~conditionC+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3:4,]

#Familiar X Trials
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp3"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,exp=="exp3"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[2:3,]
#Novel X Test
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp3"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~conditionC+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3:4,]

#Testing Conditions against chance
#Learnable Condition - Familiar X Test
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="familiarX"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp3"&testType=="familiarX"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#Learnable Condition - Novel X Test
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="novelX"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model 
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp3"&testType=="novelX"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#Unstructured Pre-Exposure - Familiar X Test
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="familiarX"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="exp3"&testType=="familiarX"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#Unstructured Pre-Exposure - Novel X Test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="novelX"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3,]

##Interaction between test type and condition
#recode test type
d$testTypeC <- ifelse(d$testType=="novelX",0.5,
                      ifelse(d$testType=="familiarX",-0.5,NA))
#maximal model (does not converge)
m <- glmer(isRight~conditionC*testTypeC+(1+testTypeC|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp3"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~conditionC*testTypeC+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m,method="Wald")

#### Experiment 3 - Signal Detection Analyses ####

## Sensitivity (d-prime)
#recode condition
subj_overall$conditionC <- ifelse(subj_overall$condition=="Learnable Pre-Exposure",0.5,
                                  ifelse(subj_overall$condition=="Unstructured Pre-Exposure",-0.5,NA))
##all data
m <- lm(dprime~conditionC,data=subset(subj_overall,exp=="exp3"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp3"&conditionC==0.5))
summary(m)

#Testing against chance: Non-Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_overall,exp=="exp3"&conditionC==-0.5))
summary(m)

##Familiar X trials
#recode condition
subj_testType$conditionC <- ifelse(subj_testType$condition=="Learnable Pre-Exposure",0.5,
                                   ifelse(subj_testType$condition=="Unstructured Pre-Exposure",-0.5,NA))
m <- lm(dprime~conditionC,data=subset(subj_testType,exp=="exp3"&testType=="familiarX"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp3"&testType=="familiarX"&conditionC==0.5))
summary(m)

#Testing against chance: Unstructured Pre-Exposure Condition
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp3"&testType=="familiarX"&conditionC==-0.5))
summary(m)

##Novel X trials
m <- lm(dprime~conditionC,data=subset(subj_testType,exp=="exp3"&testType=="novelX"))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp3"&testType=="novelX"&conditionC==0.5))
summary(m)

#Testing against chance: Unstructured Pre-Exposure
m <- lm(dprime~1,data=subset(subj_testType,exp=="exp3"&testType=="novelX"&conditionC==-0.5))
summary(m)

# Interaction between Novel X and Familiar X trials
subj_testType$testTypeC <- ifelse(subj_testType$testType=="novelX",0.5,
                                  ifelse(subj_testType$testType=="familiarX",-0.5,NA))
m <- lm(dprime~conditionC*testTypeC,data=subset(subj_testType,exp=="exp3"))
summary(m)

## Response Bias (c)

##all data
m <- lm(c~conditionC,data=subset(subj_overall,exp=="exp3"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp3"&conditionC==0.5))
summary(m)

#Testing against chance: Unstructured Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp3"&conditionC==-0.5))
summary(m)

##Familiar X trials
m <- lm(c~conditionC,data=subset(subj_testType,exp=="exp3"&testType=="familiarX"))
summary(m)

#Testing against chance: Learnable Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp3"&testType=="familiarX"&conditionC==0.5))
summary(m)

#Testing against chance: Unstructured Pre-Exposure Condition
m <- lm(c~1,data=subset(subj_testType,exp=="exp3"&testType=="familiarX"&conditionC==-0.5))
summary(m)

##Novel X trials
m <- lm(c~conditionC,data=subset(subj_testType,exp=="exp3"&testType=="novelX"))
summary(m)

#Testing against chance: Learnable Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp3"&testType=="novelX"&conditionC==0.5))
summary(m)

#Testing against chance: Unstructured Pre-Exposure
m <- lm(c~1,data=subset(subj_testType,exp=="exp3"&testType=="novelX"&conditionC==-0.5))
summary(m)

# Interaction between Novel X and Familiar X trials
m <- lm(c~conditionC*testTypeC,data=subset(subj_testType,exp=="exp3"))
summary(m)

#### Experiment 3 - Correlation Familiar X - Novel X ####

##Correlation
#Learnable Pre-Exposure
c <- corr.test(subset(subj_accuracy_wide, condition=="Learnable Pre-Exposure" & exp=="exp3")[,c("novelX","familiarX")])
c
c$p
#Unstructured Pre-Exposure
c <- corr.test(subset(subj_accuracy_wide, condition=="Unstructured Pre-Exposure" & exp=="exp3")[,c("novelX","familiarX")])
c
c$p

#interaction
m <- lm(novelX~familiarX*condition, subset(subj_accuracy_wide,exp=="exp3"))
summary(m)

#### OVERALL ANALYSIS EXPS 1 - 3 ####

### Overall effect of condition across all experiments

##All test trials
#Reference level is Learnable Pre-Exposure condition - the estimates of pairwise comparisons to this condition can be read in the summary() command.
#maximal model (does not converge)
m <- glmer(isRight~condition+(1|participantCode)+(1+condition|stimulus),data=filter(d, condition!="Learnable Pre-Exposure Only"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~condition+(1|participantCode),data=filter(d, condition!="Learnable Pre-Exposure Only"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
Anova(m, type="III")

##Familiar X test trials
#Reference level is Learnable Pre-Exposure condition - the estimates of pairwise comparisons to this condition can be read in the summary() command.
#maximal model (does not converge)
m <- glmer(isRight~condition+(1|participantCode)+(1+condition|stimulus),data=filter(d, testType=="familiarX"&condition!="Learnable Pre-Exposure Only"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~condition+(1|participantCode)+(1|stimulus),data=filter(d, testType=="familiarX"&condition!="Learnable Pre-Exposure Only"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
Anova(m, type="III")

##Novel X test trials
#Reference level is Learnable Pre-Exposure condition - the estimates of pairwise comparisons to this condition can be read in the summary() command.
#maximal model (does not converge)
m <- glmer(isRight~condition+(1|participantCode)+(1+condition|stimulus),data=filter(d, testType=="novelX"&condition!="Learnable Pre-Exposure Only"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~condition+(1|participantCode),data=filter(d, testType=="novelX"&condition!="Learnable Pre-Exposure Only"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
Anova(m, type="III")

#### Signal Detection Analysis

####Sensitivity

##All Test trials
m <- lm(dprime~condition,data=filter(subj_overall, condition!="Learnable Pre-Exposure Only"))
summary(m)
Anova(m, type="III")

##Familiar X
m <- lm(dprime~condition,data=filter(subj_testType, testType=="familiarX"&condition!="Learnable Pre-Exposure Only"))
summary(m)
Anova(m, type="III")

##Novel X
m <- lm(dprime~condition,data=filter(subj_testType, testType=="novelX"&condition!="Learnable Pre-Exposure Only"))
summary(m)
Anova(m, type="III")

#######################
#### Experiment S1 ####
#######################

#### Experiment S1 - Logistic Mixed-Effects Analyses ####
##all data
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="s1"),family=binomial,glmerControl(optimizer="bobyqa"))
#removing the random effect of stimulus which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="s1"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##by test type
#Familiar X Test
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="s1"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa")) #convergence warning
#removing the random effect of stimulus which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="s1"&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#Novel X Test
#maximal model (does not converge)
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="s1"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
#removing the random effect of stimulus which leads to a singular model fit (virtually identical results)
m <- glmer(isRight~1+(1|participantCode),data=subset(d,exp=="s1"&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##Comparison to No Pre-Exposure condition from Exp 2
#All trials
#recode condition
d$conditionC <- ifelse(d$condition=="Learnable Pre-Exposure Only",0.5,
                       ifelse(d$condition=="No Pre-Exposure",-0.5,NA))
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure")),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure")),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)

#Familiar X
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,(exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure"))&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,(exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure"))&testType=="familiarX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)

#Novel X
#maximal model (does not converge)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,(exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure"))&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
#final converging model
m <- glmer(isRight~conditionC+(1|participantCode),data=subset(d,(exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure"))&testType=="novelX"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)

#### Experiment S1 - Signal Detection Analyses ####

#### Sensitivity (d-prime)

#Overall- all data
m <- lm(dprime~1,data=subset(subj_overall,exp=="s1"))
summary(m)

#Familiar X
m <- lm(dprime~1,data=subset(subj_testType,exp=="s1"&testType=="familiarX"))
summary(m)

#Novel X
m <- lm(dprime~1,data=subset(subj_testType,exp=="s1"&testType=="novelX"))
summary(m)

##Comparison to No Pre-Exposure condition from Exp 2
#Overall
#recode condition
subj_overall$conditionC <- ifelse(subj_overall$condition=="Learnable Pre-Exposure Only",0.5,
                                  ifelse(subj_overall$condition=="No Pre-Exposure",-0.5,NA))
m <- lm(dprime~conditionC,data=subset(subj_overall,exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure")))
summary(m)

#Familiar X
#recode condition
subj_testType$conditionC <- ifelse(subj_testType$condition=="Learnable Pre-Exposure Only",0.5,
                                   ifelse(subj_testType$condition=="No Pre-Exposure",-0.5,NA))
m <- lm(dprime~conditionC,data=subset(subj_testType,(exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure"))&testType=="familiarX"))
summary(m)

#Novel X
m <- lm(dprime~conditionC,data=subset(subj_testType,(exp=="s1"|(exp=="exp2"&condition=="No Pre-Exposure"))&testType=="novelX"))
summary(m)

#### Experiment S1 - Correlation Familiar X - Novel X ####
##Correlation
c <- corr.test(subset(subj_accuracy_wide, exp=="s1")[,c("novelX","familiarX")])
c
c$p
#check for outliers (none found)
outlierTest(lm(novelX~familiarX,subset(subj_accuracy_wide,exp=="s1")))


################
#### Plots #####
################

#Experiment 1
#Familiar X Test
p1 <- ggplot(subset(test_type_exp,testType=="familiarX"&exp=="exp1"),aes(condition,accuracy,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="familiarX"&exp=="exp1"), aes(y=acc), width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#4DAF4A","#E41A1C"))+
  scale_color_manual(values=c("#4DAF4A","#E41A1C"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Familiar X Test")+
  xlab("Condition")

#Novel X Test
p2 <- ggplot(subset(test_type_exp,testType=="novelX"&exp=="exp1"),aes(condition,accuracy,fill=condition,color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="novelX"&exp=="exp1"),aes(y=acc), width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#4DAF4A","#E41A1C"))+
  scale_color_manual(values=c("#4DAF4A","#E41A1C"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Novel X Test")+
  xlab("Condition")
plot_grid(p1,p2, labels=c("A","B"),label_size=20)
ggsave("figures/exp1_test.tiff",width=11,height=6,dpi=300)

#Experiment 2
#Familiar X Test
p1 <- ggplot(subset(test_type_exp,testType=="familiarX"&exp=="exp2"),aes(condition,accuracy,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="familiarX"&exp=="exp2"), aes(y=acc),width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","No\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#E41A1C","#377EB8","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Familiar X Test")+
  xlab("Condition")

#Novel X Test
p2 <- ggplot(subset(test_type_exp,testType=="novelX"&exp=="exp2"),aes(condition,accuracy,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="novelX"&exp=="exp2"),aes(y=acc), width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","No\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#E41A1C","#377EB8","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Novel X Test")+
  xlab("Condition")
plot_grid(p1,p2, labels=c("A","B"),label_size=18)
ggsave("figures/exp2_test.tiff",width=13,height=7,dpi=300)

#Experiment 3
#Familiar X Test
p1 <- ggplot(subset(test_type_exp,testType=="familiarX"&exp=="exp3"),aes(condition,accuracy,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="familiarX"&exp=="exp3"),aes(y=acc), width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Unstructured Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Unstructured\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#4DAF4A","#54278f"))+
  scale_color_manual(values=c("#4DAF4A","#54278f"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Familiar X Test")+
  xlab("Condition")

#Novel X Test
p2 <- ggplot(subset(test_type_exp,testType=="novelX"&exp=="exp3"),aes(condition,accuracy,fill=condition,color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="novelX"&exp=="exp3"), aes(y=acc),width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Unstructured Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Unstructured\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#4DAF4A","#54278f"))+
  scale_color_manual(values=c("#4DAF4A","#54278f"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Novel X Test")+
  xlab("Condition")
plot_grid(p1,p2, labels=c("A","B"),label_size=20)
ggsave("figures/exp3_test.tiff",width=11,height=6,dpi=300)

#s1
#Familiar X Test
p1 <- ggplot(subset(test_type_exp,testType=="familiarX"&exp=="s1"),aes(condition,accuracy,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="familiarX"&exp=="s1"),aes(y=acc), width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Learnable Pre-Exposure Only"),
                   labels=c("Learnable\nPre-Exposure\nOnly"))+
  scale_fill_manual(values=c("#265725"))+
  scale_color_manual(values=c("#265725"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Familiar X Test")+
  xlab("Condition")

#Novel X Test
p2 <- ggplot(subset(test_type_exp,testType=="novelX"&exp=="s1"),aes(condition,accuracy,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="novelX"&exp=="s1"),aes(y=acc), width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Learnable Pre-Exposure Only"),
                   labels=c("Learnable\nPre-Exposure\nOnly"))+
  scale_fill_manual(values=c("#265725"))+
  scale_color_manual(values=c("#265725"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Novel X Test")+
  xlab("Condition")
plot_grid(p1,p2, labels=c("A","B"),label_size=20)
ggsave("figures/exps1_test.tiff",width=11,height=6,dpi=300)

#S4: Overall accuracy plot
#Familiar X test
p1 <- ggplot(subset(test_type,testType=="familiarX"&condition!="Learnable Pre-Exposure Only"),aes(condition,accuracy,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="familiarX"&condition!="Learnable Pre-Exposure Only"),aes(y=acc), width = 0.05,height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Unstructured Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","No\n Pre-Exposure","Unstructured\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#E41A1C","#377EB8","#54278f","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_color_manual(values=c("#E41A1C","#377EB8","#54278f","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Familiar X Test")+
  xlab("Condition")

#Novel X
p2 <- ggplot(subset(test_type,testType=="novelX"&condition!="Learnable Pre-Exposure Only"),aes(condition,accuracy,fill=condition,color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="novelX"&condition!="Learnable Pre-Exposure Only"),aes(y=acc), width = 0.05, height=0.01, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=20)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Unstructured Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","No\n Pre-Exposure","Unstructured\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#E41A1C","#377EB8","#54278f","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_color_manual(values=c("#E41A1C","#377EB8","#54278f","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(legend.position="none")+
  ylab("Accuracy - Novel X Test")+
  xlab("Condition")
plot_grid(p1,p2, labels=c("A","B"),label_size=20)
ggsave("figures/all_exps_test.tiff",width=14.8,height=8,dpi=300)


#####correlation plots
#Experiment 1
ggplot(subset(subj_accuracy_wide,exp=="exp1"),aes(familiarX,novelX, color=condition,linetype=condition,shape=condition))+
  geom_jitter(width=0.02,height=0.02,size=3)+
  geom_smooth(method=lm,se =F,size=1.5)+
  scale_color_manual(name="Condition",values=c("#4DAF4A","#E41A1C"))+
  scale_linetype_discrete(name="Condition")+
  scale_shape_discrete(name="Condition")+
  scale_x_continuous(breaks=seq(0,1,0.1))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  ylab("Accuracy - Novel X")+
  xlab("Accuracy - Familiar X")+
  theme_classic(base_size=18)+
  theme(legend.position=c(0.25,0.9))
ggsave("figures/exp1_correlation.tiff",width=9,height=7,dpi=300)


#Experiment 2
ggplot(subset(subj_accuracy_wide,exp=="exp2"),aes(familiarX,novelX, color=condition,linetype=condition,shape=condition))+
  geom_jitter(width=0.02,height=0.02,size=3)+
  geom_smooth(method=lm,se =F,size=1.5)+
  scale_color_manual(name="Condition",values=c("#4DAF4A","#377EB8","#E41A1C"))+
  scale_linetype_manual(name="Condition",values=c(1,4,2))+
  scale_shape_manual(name="Condition",values=c(16,15,17))+
  scale_x_continuous(breaks=seq(0,1,0.1),limits=c(0.1,1.08))+
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0.1,1.08))+
  ylab("Accuracy - Novel X")+
  xlab("Accuracy - Familiar X")+
  theme_classic(base_size=18)+
  theme(legend.position=c(0.25,0.9))
ggsave("figures/exp2_correlation.tiff",width=9,height=7,dpi=300)

#Experiment 3
ggplot(subset(subj_accuracy_wide,exp=="exp3"),aes(familiarX,novelX, color=condition,linetype=condition,shape=condition))+
  geom_jitter(width=0.02,height=0.02,size=3)+
  geom_smooth(method=lm,se =F,size=1.5)+
  scale_color_manual(name="Condition",values=c("#4DAF4A","#54278f"))+
  scale_linetype_manual(name="Condition",values=c(1,5))+
  scale_shape_manual(name="Condition",values=c(16,18))+
  scale_x_continuous(breaks=seq(0,1,0.1),limits=c(0.15,1.06))+
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0.15,1.06))+
  ylab("Accuracy - Novel X")+
  xlab("Accuracy - Familiar X")+
  theme_classic(base_size=18)+
  theme(legend.position=c(0.25,0.9))
ggsave("figures/exp3_correlation.tiff",width=9,height=7,dpi=300)


#S1
ggplot(subset(subj_accuracy_wide,exp=="s1"),aes(familiarX,novelX, color=condition,linetype=condition,shape=condition))+
  geom_jitter(width=0.02,height=0.02,size=3)+
  geom_smooth(method=lm,se =F,size=1.5)+
  scale_color_manual(name="Condition",values=c("#265725"))+
  scale_linetype_manual(name="Condition",values=c(1))+
  scale_shape_manual(name="Condition",values=c(16))+
  scale_x_continuous(breaks=seq(0,1,0.1),limits=c(0.2,1.06))+
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0.2,1.06))+
  ylab("Accuracy - Novel X")+
  xlab("Accuracy - Familiar X")+
  theme_classic(base_size=18)+
  theme(legend.position=c(0.25,0.9))
ggsave("figures/exps1_correlation.tiff",width=9,height=7,dpi=300)


