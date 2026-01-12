### Load libraries
library(rstudioapi)
library(MASS)
library(lme4)
library(lmerTest)
library(lmtest)
library(sandwich)
library(dr)
library(MAVE)
library(parallel)
library(LaplacesDemon)
library(MatchIt)
library(cobalt)
library(foreach)
library(doParallel)
library(GGally)
library(reshape2)
library(ggpubr)

setwd(dirname(getActiveDocumentContext()$path))

rm(list=ls())
source("functions.R")
ncores <- parallel::detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)



####  Model 1 #############
#
# --- Define Ranges for Power Curves --
run_reg_sim1 <- function(n, z_null) {
  
  ## fixed parameters
  p <- 5
  rho <- 0.5
  SigX <- rho^abs(outer(1:p,1:p,"-"))
  
  ## observed covariates
  X <- LaplacesDemon::rmvn(n, mu=rep(0,p), Sigma=SigX)
  
  ## treatment variables
  Z <- rbinom(n, 1, plogis(-0.5 + 0.1*X[,1]+0.2*X[,2]))
  
  # latent class indicator
  W <- numeric(n)
  W[Z==0] <- rbinom(sum(Z==0), 1, prob = 0.3)
  W[Z==1] <- rbinom(sum(Z==1), 1, prob = 0.7)
  
  ## true model
  B <- c(1, 0, 1, -0.5, -1)
  eps <- rnorm(n,0,1)
  Y <- 5*Z + X%*%B + 8*W + eps
  
  #### Naive model
  fit.naive <- summary(lm(Y ~ Z + X))
  t.naive <- (fit.naive$coefficients[2,1] - z_null)/fit.naive$coefficients[2,2]
  pval.naive <- 2*pt(abs(t.naive), df=fit.naive$df[2], lower.tail=FALSE)
  # dec.naive <- 1*(pval.naive < alpha)

  ## Latent class control
  
  ## sufficient latent confounder and latent class recovery
  SLC <- est_slc(Z,X,Y, d=1)
  class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="ward.D", k=2:4)
  hat.W2 <- class_est$`2_class`
  hat.W3 <- class_est$`3_class`
  hat.W4 <- class_est$`4_class`
  
  ## K=2
  LCC2 <- est_LCC(Z,X,Y,hat.W2,covType = "HC3")
  t.LCC2 <- (LCC2$summary[2,1] - z_null)/LCC2$summary[2,2]
  pval.LCC2 <- 2*pt(abs(t.LCC2), df=fit.naive$df[2]-2, lower.tail=FALSE)
  # dec.LCC2 <- 1*(pval.LCC2 < alpha)
  
  ## K=3
  LCC3 <- est_LCC(Z,X,Y,hat.W3,covType = "HC3")
  t.LCC3 <- (LCC3$summary[2,1] - z_null)/LCC3$summary[2,2]
  pval.LCC3 <- 2*pt(abs(t.LCC3), df=fit.naive$df[2]-3, lower.tail=FALSE)
  # dec.LCC3 <- 1*(pval.LCC3 < alpha)
  
  ## K=4
  LCC4 <- est_LCC(Z,X,Y,hat.W4,covType = "HC3")
  t.LCC4 <- (LCC4$summary[2,1] - z_null)/LCC4$summary[2,2]
  pval.LCC4 <- 2*pt(abs(t.LCC4), df=fit.naive$df[2]-4, lower.tail=FALSE)
  # dec.LCC4 <- 1*(pval.LCC4 < alpha)
  
  return(c(pval.naive, pval.LCC2, pval.LCC3, pval.LCC4))
}

# --- Define Ranges for Power Curves ---
sample_sizes <- seq(200, 1000, by=200) # Range of sample sizes
z_null0 <- 5
nsim <- 500

power_results <- list()
for(j in 1:length(sample_sizes)){
    
pvals <- foreach(r=1:nsim, .packages = c("LaplacesDemon","MatchIt",
                "cobalt","MAVE","pracma","sandwich","lmtest")) %dopar% {
                  set.seed(r)
          run_reg_sim1(n=sample_sizes[j], z_null=z_null0)
        }
power_results[[j]]  <- do.call(rbind, pvals)  
}

### alpha = 0.01
alpha0 <- 0.01
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y <= alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

### alpha = 0.05
alpha0 <- 0.05
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y <= alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

### alpha = 0.10
alpha0 <- 0.10
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y <= alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate


# df_power <- data.frame(do.call(rbind, power_results) )             
# names(df_power) <- c("n","Naive","LCC_2","LCC_3","LCC_4")
# powdat1 <- melt(df_power, id="n") 
# 
# G1 <- ggplot(data=powdat1, aes(x=n, y=value, colour=variable)) + geom_line() + geom_point() + 
#   geom_hline(yintercept=alpha0,linetype="dashed",color="black") + annotate("text", x = max(sample_sizes),
#   y=alpha0, label=paste0("Target size (",alpha0,")"),vjust =-0.5, hjust=1, color="black") + xlab("Sample size (n)") + ylab("Power") +
#   ggtitle("Model I") + theme_minimal() + scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
#   theme(legend.title = element_blank()) + scale_colour_discrete(labels = c("Naive", expression("LCC("~hat(K)~"=2)"), 
#   expression("LCC("~hat(K)~"=3)"),expression("LCC("~hat(K)~"=4)")) )
#   




############# Model 2  ##################
#
######################################
run_reg_sim2 <- function(n, z_null){
  
  ## fixed parameters
  p <- 5
  rho <- 0.5
  SigX <- rho^abs(outer(1:p,1:p,"-"))
  
  ## observed covariates
  X <- rmvn(n,mu=rep(0,p), Sigma=SigX)
  
  ## Latent Class
  W <- rbinom(n, 1, plogis(-0.2*X[,1] - 1.2*X[,2]))
  
  ## treatment variables
  Z <- rbinom(n, 1, plogis(0.05 + 1.5*X[,2]))
  
  ## true model
  B <- c(1,0,1, -0.5,-1)
  eps <- rnorm(n,0,1) 
  Y <- 5*Z + X%*%B  + 10*W  + eps
  
  #### Naive model
  fit.naive <- summary(lm(Y ~ Z + X))
  t.naive <- (fit.naive$coefficients[2,1] - z_null)/fit.naive$coefficients[2,2]
  pval.naive <- 2*pt(abs(t.naive), df=fit.naive$df[2], lower.tail=FALSE)
  
  ## Latent class control
  
  ## sufficient latent confounder and latent class recovery
  SLC <- est_slc(Z,X,Y, d=1)
  class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="ward.D", k=2:4)
  hat.W2 <- class_est$`2_class`
  hat.W3 <- class_est$`3_class`
  hat.W4 <- class_est$`4_class`
  
  ## K=2
  LCC2 <- est_LCC(Z,X,Y,hat.W2)
  t.LCC2 <- (LCC2$summary[2,1] - z_null)/LCC2$summary[2,2]
  pval.LCC2 <- 2*pt(abs(t.LCC2), df=fit.naive$df[2], lower.tail=FALSE)
  # dec.LCC2 <- 1*(pval.LCC2 > alpha)
  
  ## K=3
  LCC3 <- est_LCC(Z,X,Y,hat.W3)
  t.LCC3 <- (LCC3$summary[2,1] - z_null)/LCC3$summary[2,2]
  pval.LCC3 <- 2*pt(abs(t.LCC3), df=fit.naive$df[2], lower.tail=FALSE)
  # dec.LCC3 <- 1*(pval.LCC3 > alpha)
  
  ## K=4
  LCC4 <- est_LCC(Z,X,Y,hat.W4)
  t.LCC4 <- (LCC4$summary[2,1] - z_null)/LCC4$summary[2,2]
  pval.LCC4 <- 2*pt(abs(t.LCC4), df=fit.naive$df[2], lower.tail=FALSE)
  # dec.LCC4 <- 1*(pval.LCC4 > alpha)
  
  return(c(pval.naive, pval.LCC2, pval.LCC3, pval.LCC4))
}


# --- Define Ranges for Power Curves ---
sample_sizes <- seq(200, 1000, by=200) # Range of sample sizes
z_null0 <- 5
nsim <- 500

power_results <- list()
for(j in 1:length(sample_sizes)){
  
  pvals <- foreach(r=1:nsim, .packages = c("LaplacesDemon","MatchIt",
                             "cobalt","MAVE","pracma","sandwich","lmtest")) %dopar% {
                               set.seed(r)
                               run_reg_sim2(n=sample_sizes[j], z_null=z_null0)
                                           }
  power_results[[j]]  <- do.call(rbind, pvals)  
}

### alpha = 0.01
alpha0 <- 0.01
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y <= alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

### alpha = 0.05
alpha0 <- 0.05
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y <= alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

### alpha = 0.10
alpha0 <- 0.10
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y <= alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

# df_power <- data.frame(do.call(rbind, power_results) )             
# names(df_power) <- c("n","Naive","LCC_2","LCC_3","LCC_4")
# powdat2 <- melt(df_power, id="n") 
# 
# G2 <- ggplot(data=powdat2, aes(x=n, y=value, colour=variable)) + geom_line() + geom_point() + 
#   geom_hline(yintercept=0.8,linetype="dashed",color="black") + annotate("text", x = max(sample_sizes),
#   y=0.8, label="Target Power (0.8)",vjust =-0.5, hjust=1, color="black") + xlab("Sample size (n)") + ylab("Power") +
#   ggtitle("Model II") + theme_minimal() + scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
#   theme(legend.title = element_blank()) + scale_colour_discrete(labels = c("Naive", expression("LCC("~hat(K)~"=2)"), 
#   expression("LCC("~hat(K)~"=3)"),expression("LCC("~hat(K)~"=4)")) )






############# Model 3 ##################
#
######################################
run_reg_sim3 <- function(n, z_null) {
  
  ## fixed parameters
  p <- 5
  rho <- 0.5
  SigX <- rho^abs(outer(1:p,1:p,"-"))
  
  ## observed covariates
  X <- rmvn(n,mu=rep(0,p), Sigma=SigX)
  
  ## latent variables
  Z <- rbinom(n, 1, plogis(-0.5 + 0.1*X[,1]+0.2*X[,2]))
  
  ## treatment variable
  W <- rbinom(n, 1, plogis(0.05 - 1.2*Z - 0.2*X[,1]))
  
  ## true model
  B <- c(1,0,1, -0.5,-1) 
  eps <- rnorm(n,0,1)
  Y <- 5*Z + X%*%B + 15*W + eps
  
  #### Naive model
  fit.naive <- summary(lm(Y ~ Z + X))
  t.naive <- (fit.naive$coefficients[2,1] - z_null)/fit.naive$coefficients[2,2]
  pval.naive <- 2*pt(abs(t.naive), df=fit.naive$df[2], lower.tail=FALSE)
  
  ## Latent class control
  
  ## sufficient latent confounder and latent class recovery
  SLC <- est_slc(Z,X,Y, d=1)
  class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="ward.D", k=2:4)
  hat.W2 <- class_est$`2_class`
  hat.W3 <- class_est$`3_class`
  hat.W4 <- class_est$`4_class`
  
  ## K=2
  LCC2 <- est_LCC(Z,X,Y,hat.W2)
  t.LCC2 <- (LCC2$summary[2,1] - z_null)/LCC2$summary[2,2]
  pval.LCC2 <- 2*pt(abs(t.LCC2), df=fit.naive$df[2], lower.tail=FALSE)
  
  ## K=3
  LCC3 <- est_LCC(Z,X,Y,hat.W3)
  t.LCC3 <- (LCC3$summary[2,1] - z_null)/LCC3$summary[2,2]
  pval.LCC3 <- 2*pt(abs(t.LCC3), df=fit.naive$df[2], lower.tail=FALSE)
  
  ## K=4
  LCC4 <- est_LCC(Z,X,Y,hat.W4)
  t.LCC4 <- (LCC4$summary[2,1] - z_null)/LCC4$summary[2,2]
  pval.LCC4 <- 2*pt(abs(t.LCC4), df=fit.naive$df[2], lower.tail=FALSE)
  
  return(c(pval.naive, pval.LCC2, pval.LCC3, pval.LCC4))
}


# --- Define Ranges for Power Curves ---
sample_sizes <- seq(200, 1000, by=200) # Range of sample sizes
z_null0 <- 5
nsim <- 500

power_results <- list()
for(j in 1:length(sample_sizes)){
  
  pvals <- foreach(r=1:nsim, .packages = c("LaplacesDemon","MatchIt",
                                 "cobalt","MAVE","pracma","sandwich","lmtest")) %dopar% {
                                   set.seed(r)
                                   run_reg_sim3(n=sample_sizes[j], z_null=z_null0)
                              }
  power_results[[j]]  <- do.call(rbind, pvals)  
}

### alpha = 0.01
alpha0 <- 0.01
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y < alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

### alpha = 0.05
alpha0 <- 0.05
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y < alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

### alpha = 0.10
alpha0 <- 0.10
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y < alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

# df_power <- data.frame(do.call(rbind, power_results) )             
# names(df_power) <- c("n","Naive","LCC_2","LCC_3","LCC_4")
# powdat3 <- melt(df_power, id="n") 
# 
# G3 <- ggplot(data=powdat3, aes(x=n, y=value, colour=variable)) + geom_line() + geom_point() + 
#   geom_hline(yintercept=0.8,linetype="dashed",color="black") + annotate("text", x = max(sample_sizes),
#   y=0.8, label="Target Power (0.8)",vjust =-0.5, hjust=1, color="black") + xlab("Sample size (n)") + ylab("Power") +
#   ggtitle("Model III") + theme_minimal() + scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
#   theme(legend.title = element_blank()) + scale_colour_discrete(labels = c("Naive", expression("LCC("~hat(K)~"=2)"), 
#   expression("LCC("~hat(K)~"=3)"),expression("LCC("~hat(K)~"=4)")) )





############# Model 4  ##################
#
######################################
run_reg_sim4 <- function(n, z_null) {

  ## fixed parameters
  p <- 5
  rho <- 0.5
  SigX <- rho^abs(outer(1:p,1:p,"-"))
  
  ## observed covariates
  X <- LaplacesDemon::rmvn(n,mu=rep(0,p), Sigma=SigX)
  
  ## Ill-defined categorical variable
  V <- rbinom(n, 1, 0.5)  ## ill-defined classes
  V0 <- numeric(n)    
  V0[V==0] <- sample(c(1,2), sum(V==0), replace=T, prob = c(0.35,0.65))
  V0[V==1] <- sample(c(3,4), sum(V==1), replace=T, prob = c(0.7,0.3))
  V0 <- relevel(factor(V0), ref=2) ## true classes
  V1 <- model.matrix(~V0)[,-1] ## true dummies
  
  ## treatment variables
  Z <- rbinom(n, 1, plogis(0.5*X[,1] + 2*V1[,1] + 5*V1[,2]))
  
  ## true model
  B <- c(1,0,1, -0.5,-1)
  eps <- rnorm(n,0,1) 
  Y <- 5*Z + X%*%B - 2*V1[,2] + 15*V1[,3] + eps
  
  XV <- cbind(X,V)
  colnames(XV) <- c(paste0("X",1:p),"V")
  
  #### Naive model
  fit.naive <- summary(lm(Y ~ Z + XV))
  t.naive <- (fit.naive$coefficients[2,1] - z_null)/fit.naive$coefficients[2,2]
  pval.naive <- 2*pt(abs(t.naive), df=fit.naive$df[2], lower.tail=FALSE)
  
  ## Latent class control
  
  ## sufficient latent confounder and latent class recovery
  colnames(XV) <- c(paste0("X",1:p),"V")
  SLC <- est_slc(Z,XV,Y, d=2)
  class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="single", k=2:4)
  hat.W2 <- class_est$`2_class`
  hat.W3 <- class_est$`3_class`
  hat.W4 <- class_est$`4_class`
  dat2 <- data.frame(Z,XV,Y,hat.W2,hat.W3,hat.W4)
  
  ## K=2
  LCC2 <- est_LCC(Z,XV,Y,hat.W2, covType = "HC1")
  t.LCC2 <- (LCC2$summary[2,1] - z_null)/LCC2$summary[2,2]
  pval.LCC2 <- 2*pt(abs(t.LCC2), df=fit.naive$df[2], lower.tail=FALSE)

  # K=3
  LCC3 <- est_LCC(Z,XV,Y,hat.W3, covType = "HC1")
  t.LCC3 <- (LCC3$summary[2,1] - z_null)/LCC3$summary[2,2]
  pval.LCC3 <- 2*pt(abs(t.LCC3), df=fit.naive$df[2], lower.tail=FALSE)

  # K=4
  LCC4 <- est_LCC(Z,XV,Y,hat.W4, covType = "HC1")
  t.LCC4 <- (LCC4$summary[2,1] - z_null)/LCC4$summary[2,2]
  pval.LCC4 <- 2*pt(abs(t.LCC4), df=fit.naive$df[2], lower.tail=FALSE)
 
  return(c(pval.naive, pval.LCC2, pval.LCC3, pval.LCC4))
}

# --- Define Ranges for Power Curves ---
sample_sizes <- seq(200, 1000, by=200) # Range of sample sizes
z_null0 <- 5
nsim <- 500

power_results <- list()
for(j in 1:length(sample_sizes)){
  
  pvals <- foreach(r=1:nsim, .packages = c("LaplacesDemon","MatchIt",
                               "cobalt","MAVE","pracma","sandwich","lmtest")) %dopar% {
                                 set.seed(r)
                                 run_reg_sim4(n=sample_sizes[j], z_null=z_null0)
                             }
  power_results[[j]]  <- do.call(rbind, pvals)  
}

### alpha = 0.01
alpha0 <- 0.01
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y < alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

### alpha = 0.05
alpha0 <- 0.05
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y < alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate

### alpha = 0.10
alpha0 <- 0.10
rej.rate <- do.call(rbind, lapply(power_results, function(x) apply(x,2, function(y) mean(y < alpha0))))
colnames(rej.rate) <- c("Naive","LCC(K=2)","LCC(K=3)","LCC(K=4)")
rownames(rej.rate) <- sample_sizes
rej.rate


# df_power <- data.frame(do.call(rbind, power_results) )             
# names(df_power) <- c("n","Naive","LCC_2","LCC_3","LCC_4")
# powdat4 <- melt(df_power, id="n") 
# 
# G4 <- ggplot(data=powdat4, aes(x=n, y=value, colour=variable)) + geom_line() + geom_point() + 
#   geom_hline(yintercept=0.8,linetype="dashed",color="black") + annotate("text", x = max(sample_sizes),
#   y=0.8, label="Target Power (0.8)",vjust =-0.5, hjust=1, color="black") + xlab("Sample size (n)") + ylab("Power") +
#   ggtitle("Model IV") + theme_minimal() + scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
#   theme(legend.title = element_blank()) + scale_colour_discrete(labels = c("Naive", expression("LCC("~hat(K)~"=2)"), 
#   expression("LCC("~hat(K)~"=3)"),expression("LCC("~hat(K)~"=4)")) )



## Combine plots
ggarrange(G1, G2, common.legend = TRUE)

ggarrange(G3, G4, common.legend = TRUE)

# Stop the cluster
stopCluster(cl)
