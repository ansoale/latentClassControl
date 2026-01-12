proj <- function(x) {
  return(x%*%solve(crossprod(x))%*%t(x))
}

## estimate sufficient predictor for latent confounder
## require MAVE package
est_slc <- function(Z,X,Y, d=1){
  n <- dim(X)[1]
  ZX <- scale(cbind(Z,X), center = T, scale = F)
  PX <- ZX%*%solve(crossprod(ZX))%*%t(ZX)
  Yc <- scale(Y, center = T, scale = F)
  U <- (diag(n)-PX)%*%Yc
  sdr <- mave(U ~ ZX, method="CSOPG")
  G.hat <- sdr$dir[[d]]
  Xw <- sdr$x%*%G.hat
  return(list(Gamma = G.hat, U=as.vector(U), Xw= Xw))
}

## estimate the latent classes
est_latClass <- function(dat, linkage="ward.D", k=c(2,3,4)){
  
  hc <- hclust(dist(dat), method=linkage)
  plot(hc)
  
  if(length(k)==1){
    W.hat <- factor(cutree(hc,k))
  } else {
    W.hat <- factor(cutree(hc,k[1])) 
    for(i in 2:length(k))
    W.hat <- cbind(W.hat, factor(cutree(hc,k[i])))  
  }
  W <- as.data.frame(W.hat)
  colnames(W) <- paste0(k,"_class")
  return(W)
}
  

est_LCC <- function(Z,X,Y,W,covType="HC3"){
  
  ## Y on W
  fit_Y <- lm(Y ~ W)
  resid_Y <- fit_Y$residuals
  
  ## Z on W
  fit_Z <- lm(Z ~ W)
  resid_Z <- fit_Z$residuals
  
  ## last fit
  dat <- data.frame(Z=resid_Z, X, Y=resid_Y)
  fit_obs <- lm(Y ~., data=dat)
  
  ## robust est.
  robust_results <- coeftest(fit_obs, vcov=vcovCL, type=covType, cluster=~W)
  return(list(summary = robust_results, fitted=fit_obs$fitted.values, resid=fit_obs$residuals))
}  
  

standvec <- function(x) return((x - mean(x))/sd(x))

ar1cor <- function(p,rho){
  if(rho==0) sigma = diag(p)
  else sigma = rho^abs(outer(1:p,1:p,'-'))
  return(sigma)
}

matpower <- function(A,n){
  A = round((A+t(A))/2,7)
  tmp = eigen(A)
  return(tmp$vectors%*%diag((tmp$values)^n)%*%t(tmp$vectors))
}


### Robust single linkage for hierarchical clustering
rob_slink <- function(X, k=10, alpha=0.5) {
  n <- nrow(X)
  d <- ncol(X)
  
  # Step 1: k-NN and r_k
  knn <- get.knn(X, k = k)
  r_k <- knn$nn.dist[, k]
  
  # Step 2: Mutual k-NN graph
  edges <- c()
  for (i in 1:n) {
    for (j in knn$nn.index[i, ]) {
      if (i %in% knn$nn.index[j, ]) {
        edges <- c(edges, i, j)
      }
    }
  }
  
  g <- graph(edges = edges, directed = FALSE)
  g <- simplify(g)
  
  # Step 3: Edge pruning based on density
  edge_list <- as_edgelist(g)
  keep <- logical(nrow(edge_list))
  
  for (e in seq_len(nrow(edge_list))) {
    i <- edge_list[e, 1]
    j <- edge_list[e, 2]
    dist_ij <- sqrt(sum((X[i, ] - X[j, ])^2))
    
    if (max(r_k[i], r_k[j]) <= alpha * dist_ij) {
      keep[e] <- TRUE
    }
  }
  
  pruned_edges <- edge_list[keep, , drop = FALSE]
  
  g_pruned <- graph_from_edgelist(pruned_edges, directed = FALSE)
  
  # Step 4: Single linkage on pruned graph
  dist_matrix <- as.matrix(dist(X))
  dist_matrix[!as_adj(g_pruned, sparse = FALSE)] <- Inf
  
  hc <- hclust(as.dist(dist_matrix), method = "single")
  return(hc)
}


