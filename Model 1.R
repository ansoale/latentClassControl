### Load libraries
library(rstudioapi)
library(MASS)
library(mvtnorm)
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
library(energy)
library(mgcv)
library(pracma)
library(dr)
library(energy)
library(lsa)
library(cdcsis)
library(FNN)
library(igraph)

## Set working directory
setwd(dirname(getActiveDocumentContext()$path))



##-------------   Model I  ------------ ##

## load functions
rm(list=ls())
source("functions.R")

## Begin simulation
set.seed(1)
n <- 200

## Control variables
p <- 5
rho <- 0.5
SigX <- rho^abs(outer(1:p,1:p,"-"))
X <- LaplacesDemon::rmvn(n,mu=rep(0,p), Sigma=SigX)

## treatment variables
Z <- rbinom(n, 1, plogis(-0.5 + 0.1*X[,1]+0.2*X[,2]))

# latent class indicator
W <- numeric(n)
W[Z==0] <- rbinom(sum(Z==0), 1, prob = 0.3)
W[Z==1] <- rbinom(sum(Z==1), 1, prob = 0.7)


## true model
B <- c(1, 0, 1, -0.5, -1)
eps <- rnorm(n,0,1)
Y <- as.vector(5*Z + X%*%B + 8*W + eps)

## Correlations
cor(Z,W)
cor(Y,W)

## Naive model
fit.naive <- lm(Y ~ Z + X)
summary(fit.naive)
plot(fit.naive,1)

## sufficient latent confounder
SLC <- est_slc(Z,X,Y, d=1)
round(SLC$Gamma, 4)

par(mfrow=c(1,1))
plot(SLC$U ~ SLC$Xw)

## Latent class recovery
hc <- rob_slink(cbind(SLC$U, SLC$Xw),k=20, alpha=0.9)
plot(hc)
W.hat <- cutree(hc, k=2, h=1.2)

class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="ward.D", k=2)
W.hat <- class_est$`2_class`

par(mfrow=c(1,2))
plot(SLC$U ~ SLC$Xw, col=W.hat, main="ESSP", ylab=expression(hat(U)),
     xlab="dir1", font.main=1)
abline(h=0, lty=2, col="red3")

### final model
fit.LCC <- est_LCC(Z,X,Y,W.hat)
fit.LCC$summary

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
