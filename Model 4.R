### Load libraries
library(rstudioapi)
library(MASS)
library(MAVE)
library(LaplacesDemon)
library(MatchIt)
library(cobalt)
library(sjPlot)
library(dplyr)
library(lme4)
library(lmtest)
library(sandwich)
library(ggplot2)
library(dbscan)
library(grid)
library(stargazer)

## Set working directory
setwd(dirname(getActiveDocumentContext()$path))



##--------------- Model 4: ill-defined categorical -----##

## load functions
rm(list=ls())
source("functions.R")
ncores <- parallel::detectCores()

## Begin simulation
set.seed(1)
n <- 200

## Control variables
p <- 5
rho <- 0.5
SigX <- rho^abs(outer(1:p,1:p,"-"))
X <- LaplacesDemon::rmvn(n,mu=rep(0,p), Sigma=SigX)

## ill-defined categorical variable
V <- rbinom(n, 1, 0.5)  ## ill-defined classes
V0 <- numeric(n)    
V0[V==0] <- sample(c(1,2), sum(V==0), replace=T, prob = c(0.35,0.65))
V0[V==1] <- sample(c(3,4), sum(V==1), replace=T, prob = c(0.7,0.3))
V0 <- relevel(factor(V0), ref=2) ## true classes
V1 <- model.matrix(~V0)[,-1] ## true dummies

## treatment variables
Z <- rbinom(n, 1, plogis(0.5*X[,1] + 2*V1[,1] +  5*V1[,2]))

## check correlations
cor(Z,V1)
cor(Z,V)

## true model
B <- c(1,0,1, -0.5,-1)
eps <- rnorm(n,0,1) 
Y <- 5*Z + X%*%B - 2*V1[,2] + 15*V1[,3] + eps

## Naive model
fit.naive <- lm(Y ~ Z + X)
summary(fit.naive)
plot(fit.naive,1)

## sufficient latent confounder
XV <- cbind(X,V)
colnames(XV) <- c(paste0("X",1:p),"V")

SLC <- est_slc(Z,XV,Y, d=2)
round(SLC$Gamma, 4)


## Latent class recovery
class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="centroid", k=2)
W.hat <- class_est$`2_class`

par(mfrow=c(1,2))
plot(SLC$U ~ SLC$Xw[,1], col=W.hat, main="ESSP", ylab=expression(hat(U)),
     xlab="dir1", font.main=1)
plot(SLC$U ~ SLC$Xw[,2], col=W.hat, main="ESSP", ylab=expression(hat(U)),
     xlab="dir2", font.main=1)


### final model
fit.LCC <- est_LCC(Z,XV,Y,W.hat)
fit.LCC$summary

## residual diagnostics after correction
par(mfrow=c(1,2))
plot(fit.LCC$resid ~ fit.LCC$fitted, main="(a) Residual vs Fitted",font.main=1)
abline(h=0, lty=2)
plot(density(fit.LCC$resid), main="(b) Residual density", font.main=1)


stargazer::stargazer(fit.naive, fit.LCC$summary,  type="text", 
                     single.row = FALSE, column.labels=c("Naive", "LCC (K=2)"))


## residual diagnostics before and after correction
par(mfrow=c(2,2))
plot(fit.naive,1, col=W.hat, caption = "Naive model")
plot(density(fit.naive$resid), main="(b) Naive residual density", font.main=1)
plot(fit.LCC$resid ~ fit.LCC$fitted, col=W.hat, main="LCC model",
     ylab="Residual",xlab="Fitted values", font.main=1)
abline(h=0, lty=2)
plot(density(fit.LCC$resid), main="LCC residual density", font.main=1)
