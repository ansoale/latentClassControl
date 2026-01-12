library(rstudioapi)
library(MASS)
library(plotly)
library(sjmisc)
library(sjPlot)
library(GGally)
library(ggplot2)
library(lmtest)
library(sandwich)
library(lmerTest)
library(MAVE)
library(dplyr)
library(nlme)
library(ggpubr)
library(stargazer)

## Set working directory
setwd(dirname(getActiveDocumentContext()$path))

######################## Insurance data  ###################
#
############################################################
rm(list=ls())
source("functions.R")

insured0 <- read.csv("insurance.csv", header=T)
insured0 <- insured0[,c(5,1:4,6,7)]
dim(insured0)
head(insured0)

n <- dim(insured0)[1]
insured <- insured0

## charges by smoking habit
insured %>% ggplot(aes(x=charges,fill=smoker)) + 
            geom_histogram(alpha=0.5, position='identity')

## charges by region habit
insured %>% ggplot(aes(x=charges,fill=sex)) + 
  geom_histogram(alpha=0.5, position='identity')

g1 <- insured %>% ggplot(aes(x=age, y=charges, colour=smoker)) + 
            geom_point() # + theme(panel.background = element_blank())

g2 <- insured %>% ggplot(aes(x=bmi, y=charges, colour=smoker)) + 
        geom_point() +  xlab("BMI") # + theme(panel.background = element_blank()) +
       
ggarrange(g1, g2, common.legend = TRUE) 

## extract residuals and model matrix
X <- model.matrix(charges ~.-smoker, data=insured)[,-1]
Z <- ifelse(insured$smoker =="yes", 1, 0)
Y <- insured$charges

## Naive model
fit.naive <- lm(charges ~., data=insured)
summary(fit.naive)

par(mfrow=c(1,2))
plot(fit.naive,1,main="Residuals vs Fitted",font.main=1,caption="")
plot(density(fit.naive$residuals),main="Residual Density",font.main=1)

## sufficient latent confounder
SLC <- est_slc(Z,X,Y, d=2)
round(t(SLC$Gamma), 3)

par(mfrow=c(1,2))
plot(SLC$U ~ SLC$Xw[,1])
plot(SLC$U ~ SLC$Xw[,2])

## Latent class recovery
class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="ward.D", k=2)
W.hat <- class_est$`2_class`

par(mfrow=c(1,2))
plot(SLC$U ~ SLC$Xw[,1], col=W.hat, main="ESSP", ylab=expression(hat(U)),
     xlab="dir1", font.main=1)
abline(h=0, lty=2, col="red3")
plot(SLC$U ~ SLC$Xw[,2], col=W.hat, main="ESSP", ylab=expression(hat(U)),
     xlab="dir2", font.main=1)
abline(h=0, lty=2, col="red3")

### final model
fit.LCC <- est_LCC(Z,X,Y,W.hat,covType = "HC3")
fit.LCC$summary



stargazer::stargazer(fit.naive, fit.LCC$summary,  type="text", 
                     single.row = FALSE, column.labels=c("Naive", "LCC (K=2)"))

## residual diagnostics after correction
par(mfrow=c(2,2))
plot(fit.naive,1, col=W.hat, caption = "Naive model")
plot(density(fit.naive$resid), main="(b) Naive residual density", font.main=1)
plot(fit.LCC$resid ~ fit.LCC$fitted, col=W.hat, main="LCC model",
     ylab="Residual",xlab="Fitted values", font.main=1)
abline(h=0, lty=2)
plot(density(fit.LCC$resid), main="LCC residual density", font.main=1)