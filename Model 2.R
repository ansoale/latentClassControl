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



##--------------- Model 2: backdoor association -----##

## load functions
rm(list=ls())
source("functions.R")
ncores <- parallel::detectCores()

## Begin simulation
set.seed(100)
n <- 200

## Control variables
p <- 5
rho <- 0.5
SigX <- rho^abs(outer(1:p,1:p,"-"))
X <- LaplacesDemon::rmvn(n,mu=rep(0,p), Sigma=SigX)

## Latent Class
W <- rbinom(n, 1, plogis(-0.2*X[,1] - 1.2*X[,2]))

## treatment variables
Z <- rbinom(n, 1, plogis(0.05 + 1.5*X[,2]))

## check correlations
cor(Z,W)

## true model
B <- c(1,0,1, -0.5,-1)
eps <- rnorm(n,0,1) 
Y <- 5*Z + X%*%B + 10*W  + eps


## Naive model
fit.naive <- lm(Y ~ Z + X)
summary(fit.naive)
plot(fit.naive,1)

## sufficient latent confounder
SLC <- est_slc(Z,X,Y, d=1)
round(SLC$Gamma, 4)

## Latent class recovery
class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="average", k=2)
W.hat <- class_est$`2_class`

par(mfrow=c(1,2))
plot(SLC$U ~ SLC$Xw, col=W.hat, main="ESSP", ylab=expression(hat(U)),
     xlab="dir1", font.main=1)
abline(h=0, lty=2, col="red3")

### final model
fit.LCC <- est_LCC(Z,X,Y,W.hat)
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

