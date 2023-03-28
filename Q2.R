data<-read.csv("Data/Compiled2.csv")

library(ape)
library(dplyr)
library(lme4)
library(phytools)
library(ggplot2)
library(effects)
library(lmtest)

data

head(data)
str(data)

Mod1<- lmer(Heritability ~ Vert + Thermy + D.W + L.F + Behavior + (1 | Study.ID), data = data) 
summary(Mod1)

#Code fits a mixed effects model with a fixed effect for each of the variables Vert, Thermy, D.W., L.F., and Behavior, and a random intercept for each Study.ID. 
#The response variable is Heritability. The lmer() function estimates the coefficients of the fixed and random effects, as well as the residual variance. 
#The summary() function provides the estimated coefficients and their standard errors, as well as the t-value and p-value for each coefficient. 
#Additionally, the summary includes information about the estimated variance components and the residual standard deviation.


Mod2<- lmer(Heritability~Behavior + (1 | Study.ID), data=data)

#Mod 2 tests the relationship between Heritability and Behavior while accounting for the clustering of the data by Study.ID

Mod3<- lmer(Zr~Behavior + (1 | Study.ID), data=data)
summary(Mod3)

#Model looking at the effect of different behavior on the transformed value of heritability. From the summary none are statistically signficant
#However "Social" Behavior had the lowest t-value (0.023)

Mod4<- lmer(Zr~Behavior + Phylum*Class + (1 | Study.ID), data=data)
summary(Mod4)

#Model looking at transformed heritability and behavior with phylum as a fixed effect and study ID as random effect.

#Study ID needs to be kept as a random effect because this was a meta analysis and different studies were used. 

Mod5<- lmer(Zr~Behavior + Class*Family + (1 | Study.ID), data=data)
summary(Mod5)

#Test fitted Model for (Near) Singularity - > variances of one or more linear combinations of effects are zero.

#More Models
Mod6 <-lmer(Zr ~ Behavior + Genus + (1 | Study.ID), data=data)
summary(Mod6)

Mod7 <- lmer(Zr~Behavior + (1 | Study.ID)+(1 | Family)+(1 | Genus)+(1 | Phylum)+(1 | Class)+(1 | Order), data=data)
summary(Mod7)

Mod8 <- lmer(Zr~Behavior +species.epithet +(1 | Study.ID)+(1 | Family)+(1 | Genus)+(1 | Phylum)+(1 | Class)+(1 | Order), data=data)
summary(Mod8)

Mod9 <-lmer(Zr ~ Behavior + Class + (1 | Study.ID), data=data)
summary(Mod9)

Mod10 <-lmer(Zr ~ Behavior + Order + (1 | Study.ID), data=data)
summary(Mod10)

Mod11 <-lmer(Zr ~ Behavior + species.epithet + (1 | Study.ID), data=data)
summary(Mod11)

#Comparing models using AIC criteria
AIC1=AIC(Mod1)
AIC2=AIC(Mod2)
AIC3=AIC(Mod3)
AIC4=AIC(Mod4)
AIC5=AIC(Mod5)
AIC6=AIC(Mod6)
AIC7=AIC(Mod7)
AIC8=AIC(Mod8)
AIC9=AIC(Mod9)
AIC10=AIC(Mod10)
AIC11=AIC(Mod11)

#Comparing likelihood ratios
LR<- anova(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9, Mod10, Mod11)

#Likelyhood ratio testing
LRT=lrtest(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9, Mod10, Mod11)

#Model 3 is the best model using these analyses

#Pedict heritibility using model 3 (best model)
data$Predicted <- predict(Mod3, data)

#Plot data
ggplot(data, aes(x=Behavior, y=Predicted)) +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+ylab("Predicted Heritability")+xlab("Behaviour")

#Save workspace
save.image(file = "Q2.RData")
