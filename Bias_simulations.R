### Load libraries
library(rstudioapi)
library(MASS)
library(lmerTest)
library(lmtest)
library(sandwich)
library(dr)
library(MAVE)
library(parallel)
library(foreach)
library(doParallel)
library(LaplacesDemon)
library(MatchIt)
library(cobalt)
library(GGally)
library(reshape2)
library(ggpubr)

## Set working directory
setwd(dirname(getActiveDocumentContext()$path))

## load functions
rm(list=ls())
source("functions.R")
ncores <- parallel::detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)



######### Model I #######

## Fixed parameters
n <- 200
p <- 5
rho <- 0.5
SigX <- rho^abs(outer(1:p,1:p,"-"))

nsim <- 500

effect_sim <- foreach(j=1:nsim, .packages = c("LaplacesDemon","MatchIt",
                    "cobalt","MAVE","pracma","sandwich","lmtest")) %dopar% {
                              
    set.seed(j)
  
  # tryCatch({ 
    
    ## observed covariates
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
    Y <- 5*Z + X%*%B + 8*W + eps
    
    #### Naive model
    fit.naive <- lm(Y ~ Z + X)
    tau.naive <- coef(fit.naive)[2]

    ## Latent class control
    
    ## sufficient latent confounder and latent class recovery
    SLC <- est_slc(Z,X,Y, d=1)
    class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="ward.D", k=2:4)
    hat.W2 <- class_est$`2_class`
    hat.W3 <- class_est$`3_class`
    hat.W4 <- class_est$`4_class`
    dat2 <- data.frame(Z,X,Y,hat.W2,hat.W3,hat.W4)
    
    ## K=2
    LCC2 <- est_LCC(Z,X,Y,hat.W2)
    tau.lcc2 <- LCC2$summary[2,1]
    
    ## K=3
    LCC3 <- est_LCC(Z,X,Y,hat.W3)
    tau.lcc3 <- LCC3$summary[2,1]
    
    ## K=4
    LCC4 <- est_LCC(Z,X,Y,hat.W4)
    tau.lcc4 <- LCC4$summary[2,1]
   
    return(c(tau.naive,tau.lcc2,tau.lcc3,tau.lcc4) )
  
      # }, error=function(e) NULL)
}

effect_size <- do.call(rbind, effect_sim)
colnames(effect_size) <- c("Naive","LCC2","LCC3","LCC4")

## bias
bias <- effect_size - 5
longdat1 <- melt(bias)
G1 <- ggplot(longdat1) + geom_boxplot(aes(x=Var2, y=value, color=Var2)) + 
      geom_hline(yintercept = 0, linetype="dashed") + 
      stat_summary(aes(x=Var2, y=value),fun=mean, geom="point", shape=18, size=3, color="red") +
      labs(title = "Model I", y=expression(hat(beta)[z]-beta[z]), x="") + 
      scale_x_discrete(labels = c("Naive", expression("LCC("~hat(K)~"=2)"), 
          expression("LCC("~hat(K)~"=3)"),expression("LCC("~hat(K)~"=4)") )) + 
          theme(legend.title = element_blank(), legend.position = "none")







########## Model II #########

## Fixed parameters
n <- 200
p <- 5
rho <- 0.5
SigX <- rho^abs(outer(1:p,1:p,"-"))

nsim <- 500

effect_sim <- foreach(j=1:nsim, .packages = c("LaplacesDemon","MatchIt",
                "cobalt","MAVE","pracma","sandwich","lmtest")) %dopar% {
                                                
  set.seed(j)
  
  tryCatch({ 
    
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
    fit.naive <- lm(Y ~ Z + X)
    tau.naive <- coef(fit.naive)[2]
    
    ## Latent class control
    
    ## sufficient latent confounder and latent class recovery
    SLC <- est_slc(Z,X,Y, d=1)
    class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="ward.D", k=2:4)
    hat.W2 <- class_est$`2_class`
    hat.W3 <- class_est$`3_class`
    hat.W4 <- class_est$`4_class`
    dat2 <- data.frame(Z,X,Y,hat.W2,hat.W3,hat.W4)
    
    ## K=2
    LCC2 <- est_LCC(Z,X,Y,hat.W2)
    tau.lcc2 <- LCC2$summary[2,1]
    
    ## K=3
    LCC3 <- est_LCC(Z,X,Y,hat.W3)
    tau.lcc3 <- LCC3$summary[2,1]
    
    ## K=4
    LCC4 <- est_LCC(Z,X,Y,hat.W4)
    tau.lcc4 <- LCC4$summary[2,1]
    
    
    return(c(tau.naive,tau.lcc2,tau.lcc3,tau.lcc4) )
    
  }, error=function(e) NULL)
  
}

effect_size <- do.call(rbind, effect_sim)
colnames(effect_size) <- c("Naive","LCC2","LCC3","LCC4")

## bias
bias <- effect_size - 5
longdat2 <- melt(bias)
G2 <- ggplot(longdat2) +  geom_boxplot(aes(x=Var2, y=value, color=Var2)) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  stat_summary(aes(x=Var2, y=value),fun=mean, geom="point", shape=18, size=3, color="red") +
  labs(title = "Model II", y=expression(hat(beta)[z]-beta[z]), x="") + 
  scale_x_discrete(labels = c("Naive", expression("LCC("~hat(K)~"=2)"), 
                              expression("LCC("~hat(K)~"=3)"),expression("LCC("~hat(K)~"=4)") )) + 
  theme(legend.title = element_blank(), legend.position = "none")













########## Model III #########

## Fixed parameters
n <- 200
p <- 5
rho <- 0.5
SigX <- rho^abs(outer(1:p,1:p,"-"))

nsim <- 500

effect_sim <- foreach(j=1:nsim, .packages = c("LaplacesDemon","MatchIt",
               "cobalt","MAVE","pracma","sandwich","lmtest")) %dopar% {
                                                
  set.seed(j)
  
  tryCatch({ 
    
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
    fit.naive <- lm(Y ~ Z + X)
    tau.naive <- coef(fit.naive)[2]
    
    ## Latent class control
    
    ## sufficient latent confounder and latent class recovery
    SLC <- est_slc(Z,X,Y, d=1)
    class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="ward.D", k=2:4)
    hat.W2 <- class_est$`2_class`
    hat.W3 <- class_est$`3_class`
    hat.W4 <- class_est$`4_class`
    dat2 <- data.frame(Z,X,Y,hat.W2,hat.W3,hat.W4)
    
    ## K=2
    LCC2 <- est_LCC(Z,X,Y,hat.W2)
    tau.lcc2 <- LCC2$summary[2,1]
    
    ## K=3
    LCC3 <- est_LCC(Z,X,Y,hat.W3)
    tau.lcc3 <- LCC3$summary[2,1]
    
    ## K=4
    LCC4 <- est_LCC(Z,X,Y,hat.W4)
    tau.lcc4 <- LCC4$summary[2,1]
    
    
    return(c(tau.naive,tau.lcc2,tau.lcc3,tau.lcc4) )
    
  }, error=function(e) NULL)
  
}

effect_size <- do.call(rbind, effect_sim)
colnames(effect_size) <- c("Naive","LCC2","LCC3","LCC4")

## bias
bias <- effect_size - 5
longdat3 <- melt(bias)

G3 <- ggplot(longdat3) + geom_boxplot(aes(x=Var2, y=value, color=Var2)) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  stat_summary(aes(x=Var2, y=value),fun=mean, geom="point", shape=18, size=3, color="red") +
  labs(title = "Model III", y=expression(hat(beta)[z]-beta[z]), x="") + 
  scale_x_discrete(labels = c("Naive", expression("LCC("~hat(K)~"=2)"), 
  expression("LCC("~hat(K)~"=3)"),expression("LCC("~hat(K)~"=4)") )) + 
  theme(legend.title = element_blank(), legend.position = "none")





############### Model IV ###########

## Fixed parameters
n <- 200
p <- 5
rho <- 0.5
SigX <- rho^abs(outer(1:p,1:p,"-"))

nsim <- 500

effect_sim <- foreach(j=1:nsim, .packages = c("LaplacesDemon","MatchIt",
                "cobalt","MAVE","pracma","sandwich","lmtest")) %dopar% {
                                              
  set.seed(j)
  
  tryCatch({ 
    
    ## observed covariates
    X <- rmvn(n,mu=rep(0,p), Sigma=SigX)
    
    ## Ill-defined categorical variable
    V <- rbinom(n, 1, 0.5)  ## ill-defined classes
    V0 <- numeric(n)    
    V0[V==0] <- sample(c(1,2), sum(V==0), replace=T, prob = c(0.35,0.65))
    V0[V==1] <- sample(c(3,4), sum(V==1), replace=T, prob = c(0.7,0.3))
    V0 <- relevel(factor(V0), ref=2) ## true classes
    V1 <- model.matrix(~V0)[,-1] ## true dummies
    
    ## treatment variables
    Z <- rbinom(n, 1, plogis(0.5*X[,1] + 2*V1[,1] +  5*V1[,2]))
    
    ## true model
    B <- c(1,0,1, -0.5,-1)
    eps <- rnorm(n,0,1) 
    Y <- 5*Z + X%*%B - 2*V1[,2] + 15*V1[,3] + eps

    XV <- cbind(X,V)
    colnames(XV) <- c(paste0("X",1:p),"V")
    
    #### Naive model
    fit.naive <- lm(Y ~ Z + XV)
    tau.naive <- coef(fit.naive)[2]
    
    ## Latent class control
    
    ## sufficient latent confounder and latent class recovery
    SLC <- est_slc(Z,XV,Y, d=2)
    class_est <- est_latClass(cbind(SLC$U, SLC$Xw), linkage="single", k=2:4)
    hat.W2 <- class_est$`2_class`
    hat.W3 <- class_est$`3_class`
    hat.W4 <- class_est$`4_class`
    dat2 <- data.frame(Z,XV,Y,hat.W2,hat.W3,hat.W4)
    
    ## K=2
    LCC2 <- est_LCC(Z,XV,Y,hat.W2)
    tau.lcc2 <- LCC2$summary[2,1]
    
    ## K=3
    LCC3 <- est_LCC(Z,XV,Y,hat.W3)
    tau.lcc3 <- LCC3$summary[2,1]
    
    ## K=4
    LCC4 <- est_LCC(Z,XV,Y,hat.W4)
    tau.lcc4 <- LCC4$summary[2,1]
    
    return(c(tau.naive,tau.lcc2,tau.lcc3,tau.lcc4) )
    
  }, error=function(e) NULL)
  
}

effect_size <- do.call(rbind, effect_sim)
colnames(effect_size) <- c("Naive","LCC2","LCC3","LCC4")

## bias
bias <- effect_size - 5
longdat4 <- melt(bias)

G4 <- ggplot(longdat4) + geom_boxplot(aes(x=Var2, y=value, color=Var2)) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  stat_summary(aes(x=Var2, y=value),fun=mean, geom="point", shape=18, size=3, color="red") +
  labs(title = "Model IV", y=expression(hat(beta)[z]-beta[z]), x="") + 
  scale_x_discrete(labels = c("Naive", expression("LCC("~hat(K)~"=2)"), 
 expression("LCC("~hat(K)~"=3)"),expression("LCC("~hat(K)~"=4)") )) + 
  theme(legend.title = element_blank(), legend.position = "none")

## Combine plots
ggarrange(G1, G2, legend=NULL)
ggarrange(G3, G4, legend=NULL)
