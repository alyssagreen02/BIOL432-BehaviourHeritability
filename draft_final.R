
getwd()
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

# Code fits a mixed effects model with a fixed effect for each of the variables Vert, Thermy, D.W., L.F., and Behavior, and a random intercept for each Study.ID. 
#The response variable is Heritability. The lmer() function estimates the coefficients of the fixed and random effects, as well as the residual variance. 
#The summary() function provides the estimated coefficients and their standard errors, as well as the t-value and p-value for each coefficient. 
#Additionally, the summary includes information about the estimated variance components and the residual standard deviation.


Mod2<- lmer(Heritability~Behavior + (1 | Study.ID), data=data)

# Mod 2 tests the relationship between Heritability and Behavior while accounting for the clustering of the data by Study.ID

Mod3<- lmer(Zr~Behavior + (1 | Study.ID), data=data)
summary(Mod3)

#Model looking at the effect of different behavior on the transformed value of heritability. From the summary none are statistically signficant
# However "Social" Behavior had the lowest t-value (0.023)


# Not able to run model with Behavior as Response variable because it is not numeric

Mod4<- lmer(Zr~Behavior + Phylum*Class + (1 | Study.ID), data=data)
summary(Mod4)

#Model looking at transformed heritability and behavior with phylum as a fixed effect and study ID as random effect.

# Study ID needs to be kept as a random effect because this was a meta analysis and different studies were used. 

Mod5<- lmer(Zr~Behavior + Class*Family + (1 | Study.ID), data=data)
summary(Mod5)

# Test fitted Model for (Near) Singularity - > variances of one or more linear combinations of effects are zero.


Mod6 <-lmer(Zr ~ Behavior + Genus + (1 | Study.ID), data=data)
summary(Mod6)

Mod7 <-lmer(Zr ~ Behavior + Family + (1 | Study.ID), data=data)
Mod8 <-lmer(Zr ~ Behavior + Class + (1 | Study.ID), data=data)
Mod9 <-lmer(Zr ~ Behavior + Order + (1 | Study.ID), data=data)
Mod10 <-lmer(Zr ~ Behavior + species.epithet + (1 | Study.ID), data=data)

summary(Mod10)

#Models with phylogenetic information as random variables
Mod11=lmer(Zr~Behavior + (1 | Study.ID)+(1 | Family)+(1 | Genus)+(1 | Phylum)+(1 | Class)+(1 | Order), data=data)
Mod12=lmer(Zr~Behavior +species.epithet +(1 | Study.ID)+(1 | Family)+(1 | Genus)+(1 | Phylum)+(1 | Class)+(1 | Order), data=data)

#Likelihood ratio testing
lrtest(Mod9, Mod10)

#Likelyhood ratio test for all models
lrtest(Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9, Mod10, Mod11, Mod12)
#According to this LRT 

# AIC Testing
AIC(Mod3) #Lowest AIC value
AIC(Mod4)
AIC(Mod5)
AIC(Mod6)
AIC(Mod7)
AIC(Mod8)
AIC(Mod9)
AIC(Mod10)
AIC(Mod11)
AIC(Mod12)

AIC(Mod10) - AIC(Mod9)

exp(0.5*AIC(Mod9)-0.5*AIC(Mod10))


# rough visualization of data and models 
plot(Mod6, type="residuals")

ggplot(data, aes(x = Heritability, y = predict(Mod1))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Observed Heritability") +
  ylab("Predicted Heritability")

ggplot(data, aes(x = Heritability, y = predict(Mod3))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Observed Heritability") +
  ylab("Predicted Heritability")

# Model not fitting the data

#. Trying a marginal effects model
plot(effect("Behavior", Mod6))
plot(effect("Behavior", Mod3))

# Looking at phylogenetic tree
Tree<-read.newick("Data/Dochter_etal(TREE).txt")

