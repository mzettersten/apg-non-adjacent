## Analyses for non-adjacent dependency learning paper (Zettersten, Potter & Saffran)
library(car) #version 3.0-2
library(cowplot) # version 0.9.3
library(formattable) # version 0.2.0.1
library(lme4) # version 1.1-18-1
library(plyr) # version 1.8.4
library(psych) # version 1.8.4
library(stats) # version 3.5.2
library(tidyverse) # version 1.2.1
source("summarizeData.R") #set of helper functions for summarizing data

#### read in complete data set ####
d <- read.table("APG_data.txt", header=T)

#### demographics & overall means ####
#subject means
subj_overall <- d %>%
  group_by(participantCode, exp, condition) %>%
  summarize(age=age[1],
            gender = gender[1],
            nativeLangEng = nativeLangEng[1],
            acc=mean(isRight),
            rt=mean(rt[rt<5000]))

#total demographics
demo <- subj_overall %>%
  group_by(exp) %>%
  summarize(
    N_learnable=sum(condition=="Learnable Pre-Exposure"),
    N_nonLearnable=sum(condition=="Non-Learnable Pre-Exposure"),
    N_noPre =sum(condition=="No Pre-Exposure"),
    N_unstructured = sum(condition=="Unstructured Pre-Exposure"),
    N_female = sum(gender=="Female"),
    N_nativeLangEng=sum(nativeLangEng==1),
    mean_age = mean(age,na.rm=T),
    sd_age = sd(age))

#split accuracy by test type
subj_testType <- d %>%
  group_by(participantCode,exp,testType, condition) %>%
  summarize(acc=mean(isRight),rt=mean(rt[rt<5000]))


#all experiment conditions
overall_exp <-  summarySE(subj_overall, "acc", groupvars=c("exp","condition"))
overall_exp$lowerCI <-  overall_exp$acc - overall_exp$ci
overall_exp$upperCI <-  overall_exp$acc + overall_exp$ci

#test type
test_type_exp <-  summarySE(subj_testType, "acc", groupvars=c("exp","condition","testType"))
test_type_exp$lowerCI <-  test_type_exp$acc - test_type_exp$ci
test_type_exp$upperCI <-  test_type_exp$acc + test_type_exp$ci

#create a wide format data set to correlate recognition and generalization accuracy.
subj_accuracy_wide <- spread(subset(subj_testType, select=-c(rt)), testType, acc)

#### Experiment 1 ####
#recode condition
d$conditionC <- ifelse(d$condition=="Learnable Pre-Exposure",0.5,
                    ifelse(d$condition=="Non-Learnable Pre-Exposure",-0.5,NA))
##all data
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp1"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
Anova(m, type="III") #Wald chi-squared test yields similar results to z-values across experiments

##by test type
#recognition test
#maximal model
#m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp1"&testType=="recognition"),family=binomial,glmerControl(optimizer="bobyqa")) #convergence warning
#removing the random effect with an estimated covariance of zero (identical results)
m <- glmer(isRight~conditionC+(1|participantCode)+(0+conditionC|stimulus),data=subset(d,exp=="exp1"&testType=="recognition"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#generalization test
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp1"&testType=="generalization"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##comparison to chance
#Learnable Pre-Exposure
#recognition test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp1"&testType=="recognition"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#generalization test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp1"&testType=="generalization"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#Non-Learnable Pre-Exposure
#recognition test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp1"&testType=="recognition"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#generalization test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp1"&testType=="generalization"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##Correlation
c <- corr.test(subset(subj_accuracy_wide, condition=="Learnable Pre-Exposure" & exp=="exp1")[,c("generalization","recognition")])
c
c <- corr.test(subset(subj_accuracy_wide, condition=="Non-Learnable Pre-Exposure" & exp=="exp1")[,c("generalization","recognition")])
c
c$p
#check outlier in Non-Learnable Pre-Exposure Condition
outlierTest(lm(generalization~recognition,subset(subj_accuracy_wide,exp=="exp1"&condition=="Non-Learnable Pre-Exposure")))
#correlation without outlier
c <- corr.test(subset(subj_accuracy_wide, condition=="Non-Learnable Pre-Exposure" & exp=="exp1"&participantCode!="apg_exp1_p7")[,c("generalization","recognition")])
c
c$p

#interaction
m <- lm(generalization~recognition*condition, subset(subj_accuracy_wide,exp=="exp1"))
summary(m)

#### Experiment 2 ####
#recode condition
d$conditionC <- ifelse(d$condition=="Learnable Pre-Exposure",0.5,
                    ifelse(d$condition=="Non-Learnable Pre-Exposure",-0.5,
                           ifelse(d$condition=="No Pre-Exposure",0,NA)))
d$conditionOrthContrast2 <- ifelse(d$condition=="Learnable Pre-Exposure",-1/3,
                                ifelse(d$condition=="Non-Learnable Pre-Exposure",-1/3,
                                       ifelse(d$condition=="No Pre-Exposure",2/3,NA)))

#overall analysis (Richter single-contrast apporach)
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#recognition test
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"&testType=="recognition"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#generalization test
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"&testType=="generalization"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")


##check same analysis with Abelson & Prentice approach
#m=glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode)+(1+conditionC+conditionOrthContrast2|stimulus),data=subset(d,exp=="exp2"),family=binomial) #convergence warning
m <- glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode)+(0+conditionC+conditionOrthContrast2|stimulus),data=subset(d,exp=="exp2"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
#recognition test
m <- glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode)+(0+conditionC+conditionOrthContrast2|stimulus),data=subset(d,exp=="exp2"&testType=="recognition"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
#generalization test
m <- glmer(isRight~conditionC+conditionOrthContrast2+(1|participantCode)+(0+conditionC+conditionOrthContrast2|stimulus),data=subset(d,exp=="exp2"&testType=="generalization"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)

#compare conditions
#Learnable vs. Non-Learnable
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp2"&conditionC!=0),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#Non-Learnable vs. No Pre-Exposure
d$conditionNonLearnableNoPreExp <- ifelse(d$conditionC==0.5,1.5,
                                    ifelse(d$conditionC==0,0.5,d$conditionC))
m <- glmer(isRight~conditionNonLearnableNoPreExp+(1+conditionNonLearnableNoPreExp|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&conditionNonLearnableNoPreExp!=1.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#No Pre-Exposure vs. Learnable
d$conditionLearnableNoPreExp <- ifelse(d$conditionC==-0.5,-1.5,
                                    ifelse(d$conditionC==0,-0.5,d$conditionC))
m <- glmer(isRight~conditionLearnableNoPreExp+(1|participantCode)+(1+conditionLearnableNoPreExp|stimulus),data=subset(d,exp=="exp2"&conditionLearnableNoPreExp!=-1.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

#Compare each condition to chance performance (collapsing across test trials)
#Learnable
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#No Pre-Exposure
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&conditionC==0),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")
#non-Learnable
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp2"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")

##Correlations between Recognition and Generalization Test
c <- corr.test(subset(subj_accuracy_wide, condition=="Learnable Pre-Exposure" & exp=="exp2")[,c("generalization","recognition")])
c
c$p
c <- corr.test(subset(subj_accuracy_wide, condition=="No Pre-Exposure" & exp=="exp2")[,c("generalization","recognition")])
c
c$p
c <- corr.test(subset(subj_accuracy_wide, condition=="Non-Learnable Pre-Exposure" & exp=="exp2")[,c("generalization","recognition")])
c
c$p

#interaction
subj_accuracy_wide$conditionC=ifelse(subj_accuracy_wide$condition=="Learnable Pre-Exposure",0.5,
                           ifelse(subj_accuracy_wide$condition=="Non-Learnable Pre-Exposure",-0.5,
                                  ifelse(subj_accuracy_wide$condition=="No Pre-Exposure",0,NA)))
m <- lm(generalization~recognition*conditionC, subset(subj_accuracy_wide,exp=="exp2"))
summary(m)

####Hits vs. False Alarm Rates Experiment 2
#split by test type and correct response type
subjType <- ddply(d,.(participantCode,exp,testType, correctResponse,condition),summarize,
              acc=mean(isRight))

#test type
response_type_exp <-  summarySE(subjType, "acc", groupvars=c("exp","condition","testType","correctResponse"))
response_type_exp$lowerCI <-  response_type_exp$acc - response_type_exp$ci
response_type_exp$upperCI <-  response_type_exp$acc + response_type_exp$ci

t.test(subjType$acc[subjType$correctResponse=="y"&subjType$exp=="exp2"&subjType$condition=="No Pre-Exposure"&subjType$testType=="recognition"],subjType$acc[subjType$correctResponse=="y"&subjType$exp=="exp2"&(subjType$condition=="Learnable Pre-Exposure"|subjType$condition=="Non-Learnable Pre-Exposure")&subjType$testType=="recognition"],var.equal = T)
t.test(subjType$acc[subjType$correctResponse=="n"&subjType$exp=="exp2"&subjType$condition=="Learnable Pre-Exposure"&subjType$testType=="recognition"],subjType$acc[subjType$correctResponse=="n"&subjType$exp=="exp2"&(subjType$condition=="No Pre-Exposure"|subjType$condition=="Non-Learnable Pre-Exposure")&subjType$testType=="recognition"],var.equal = T)

####Experiment 3####
#recode condition
d$conditionC <- ifelse(d$condition=="Learnable Pre-Exposure",0.5,
                       ifelse(d$condition=="Unstructured Pre-Exposure",-0.5,NA))
##all data
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp3"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[5:6,]

#recognition test
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp3"&testType=="recognition"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[5:6,]
#generalization test
m <- glmer(isRight~conditionC+(1|participantCode)+(1+conditionC|stimulus),data=subset(d,exp=="exp3"&testType=="generalization"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[5:6,]

#Testing Conditions against chance
#Learnable Condition - recognition test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="recognition"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3,]
#Learnable Condition - generalization test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="generalization"&conditionC==0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3,]
#Unstructured Pre-Exposure - recognition test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="recognition"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3,]
#Unstructured Pre-Exposure - generalization test
m <- glmer(isRight~1+(1|participantCode)+(1|stimulus),data=subset(d,exp=="exp3"&testType=="generalization"&conditionC==-0.5),family=binomial,glmerControl(optimizer="bobyqa"))
summary(m)
confint(m, method="Wald")[3,]

##Correlation
#Learnable Pre-Exposure
c <- corr.test(subset(subj_accuracy_wide, condition=="Learnable Pre-Exposure" & exp=="exp3")[,c("generalization","recognition")])
c
c$p
#Unstructured Pre-Exposure
c <- corr.test(subset(subj_accuracy_wide, condition=="Unstructured Pre-Exposure" & exp=="exp3")[,c("generalization","recognition")])
c
c$p

#interaction
m <- lm(generalization~recognition*condition, subset(subj_accuracy_wide,exp=="exp3"))
summary(m)

####accuracy plots####

#Experiment 1
#Recognition Test
p1 <- ggplot(subset(test_type_exp,testType=="recognition"&exp=="exp1"),aes(condition,acc,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="recognition"&exp=="exp1"), width = 0.05, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=acc-se,ymax=acc+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=18)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#4DAF4A","#E41A1C"))+
  scale_color_manual(values=c("#4DAF4A","#E41A1C"))+
  theme(legend.position="none")+
  ylab("Accuracy - Recognition Test")+
  xlab("Condition")

#Generalization Test
p2 <- ggplot(subset(test_type_exp,testType=="generalization"&exp=="exp1"),aes(condition,acc,fill=condition,color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="generalization"&exp=="exp1"), width = 0.05, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=acc-se,ymax=acc+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=18)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#4DAF4A","#E41A1C"))+
  scale_color_manual(values=c("#4DAF4A","#E41A1C"))+
  theme(legend.position="none")+
  ylab("Accuracy - Generalization Test")+
  xlab("Condition")
plot_grid(p1,p2, labels=c("A","B"),label_size=20)

#Experiment 2
#Recognition Test
p1 <- ggplot(subset(test_type_exp,testType=="recognition"&exp=="exp2"),aes(condition,acc,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="recognition"&exp=="exp2"), width = 0.05, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=acc-se,ymax=acc+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=18)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","No\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#E41A1C","#377EB8","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  theme(legend.position="none")+
  ylab("Accuracy - Recognition Test")+
  xlab("Condition")

#Generalization Test
p2 <- ggplot(subset(test_type_exp,testType=="generalization"&exp=="exp2"),aes(condition,acc,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="generalization"&exp=="exp2"), width = 0.05, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=acc-se,ymax=acc+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=18)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Non-Learnable\nPre-Exposure","No\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#E41A1C","#377EB8","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A"),limits=c("Non-Learnable Pre-Exposure","No Pre-Exposure","Learnable Pre-Exposure"))+
  theme(legend.position="none")+
  ylab("Accuracy - Generalization Test")+
  xlab("Condition")
plot_grid(p1,p2, labels=c("A","B"),label_size=18)

#experiment 3
#Recognition Test
p1 <- ggplot(subset(test_type_exp,testType=="recognition"&exp=="exp3"),aes(condition,acc,fill=condition, color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="recognition"&exp=="exp3"), width = 0.05, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=acc-se,ymax=acc+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=18)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Unstructured Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Unstructured\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#4DAF4A","#54278f"))+
  scale_color_manual(values=c("#4DAF4A","#54278f"))+
  theme(legend.position="none")+
  ylab("Accuracy - Recognition Test")+
  xlab("Condition")

#Generalization Test
p2 <- ggplot(subset(test_type_exp,testType=="generalization"&exp=="exp3"),aes(condition,acc,fill=condition,color=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=1.2,alpha=0.3, width=0.7)+
  geom_jitter(data=subset(subj_testType,testType=="generalization"&exp=="exp3"), width = 0.05, alpha=0.6,shape=21)+
  geom_errorbar(aes(ymin=acc-se,ymax=acc+se),color="black",position=position_dodge(.9),width=0.05, size=0.8)+
  theme_classic(base_size=18)+
  geom_hline(yintercept=0.5, linetype="dotted")+
  scale_x_discrete(name="Condition",
                   breaks=c("Unstructured Pre-Exposure","Learnable Pre-Exposure"),
                   labels=c("Unstructured\nPre-Exposure","Learnable\nPre-Exposure"),
                   limits=c("Unstructured Pre-Exposure","Learnable Pre-Exposure"))+
  scale_fill_manual(values=c("#4DAF4A","#54278f"))+
  scale_color_manual(values=c("#4DAF4A","#54278f"))+
  theme(legend.position="none")+
  ylab("Accuracy - Generalization Test")+
  xlab("Condition")
plot_grid(p1,p2, labels=c("A","B"),label_size=20)

#####correlation plots####
#Experiment 1
ggplot(subset(subj_accuracy_wide,exp=="exp1"),aes(recognition,generalization, color=condition,linetype=condition,shape=condition))+
  geom_jitter(width=0.02,height=0.02,size=3)+
  geom_smooth(method=lm,se =F,size=1.5)+
  scale_color_manual(name="Condition",values=c("#4DAF4A","#E41A1C"))+
  scale_linetype_discrete(name="Condition")+
  scale_shape_discrete(name="Condition")+
  ylab("Accuracy - Generalization")+
  xlab("Accuracy - Recognition")+
  theme_classic(base_size=18)+
  theme(legend.position=c(0.25,0.9))

#Experiment 2
ggplot(subset(subj_accuracy_wide,exp=="exp2"),aes(recognition,generalization, color=condition,linetype=condition,shape=condition))+
  geom_jitter(width=0.02,height=0.02,size=3)+
  geom_smooth(method=lm,se =F,size=1.5)+
  scale_color_manual(name="Condition",values=c("#4DAF4A","#377EB8","#E41A1C"))+
  scale_linetype_manual(name="Condition",values=c(1,4,2))+
  scale_shape_manual(name="Condition",values=c(16,15,17))+
  ylab("Accuracy - Generalization")+
  xlab("Accuracy - Recognition")+
  theme_classic(base_size=18)+
  theme(legend.position=c(0.25,0.9))

#Experiment 3
ggplot(subset(subj_accuracy_wide,exp=="exp3"),aes(recognition,generalization, color=condition,linetype=condition,shape=condition))+
  geom_jitter(width=0.02,height=0.02,size=3)+
  geom_smooth(method=lm,se =F,size=1.5)+
  scale_color_manual(name="Condition",values=c("#4DAF4A","#54278f"))+
  scale_linetype_manual(name="Condition",values=c(1,5))+
  scale_shape_manual(name="Condition",values=c(16,18))+
  ylab("Accuracy - Generalization")+
  xlab("Accuracy - Recognition")+
  theme_classic(base_size=18)+
  theme(legend.position=c(0.25,0.9))
