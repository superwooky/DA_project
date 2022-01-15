setwd("/Users/jongwookkim/Box/IU/DA")
source("functions.R")

sim <- function(a1=0.8, a2=0.1, a3=1, P=200, t=20){
  
  #a1 <- 0.8 #should be a1<1 #a parameter for spatial term
  #a2 <- 0.1 #a parameter for time
  #a3 <- 1 #scale parameter
  
  #P <- 200 #location data size
  #t <- 20 #time data size
  
  n <- P*t #total data size
  
  lat <- runif(P, min=-pi/2, max=pi/2) #latitude
  long <- runif(P, min=0, max=2*pi) #longitude
  time <- rep(1:t, each=P)
  
  dat <- data.frame(lat,long,time)
  
  R <- Cov_mat(dat, a1, a2, a3) #covariance matrix
  
  # Use svd for R^(1/2)
  sv <- svd(R)
  ud <- sv$u %*% diag(sqrt(sv$d))
  
  # Sample Z(x)
  sigma <- 1 #scale parameter
  Z <- sigma*ud %*% matrix(rnorm(n,0,1),n,1)
  
  ######################Compute MOM estimator#########################

  H <- 10         # Number of distances at which to estimate
  Hvec <- as.vector(seq(0, 1, by = 1/H))[1:H]
  eps <- 0.10     # use to bin distances in estimation
  
  #H <- 20         # Number of distances at which to estimate
  #Hvec <- as.vector(seq(0, 1.5, by = 1/H))[1:H]
  #eps <- 0.05     # use to bin distances in estimation
  
  sDvec <- Dfun(Lmat=dat[,1:2], P, t) #spherical distances
  
  tDmat <- as.matrix(dist(dat$time, method="manhattan")) #time distances
  tDvec  <-tDmat[upper.tri(tDmat,diag = TRUE)]
  
  mom_est <- G_Hat(Z, sDvec, Hvec, tDvec ,eps, P, t) #MOM estimate
  R_mom <- mom_est$G
  W <- 1/mom_est$N #weight
  
  #theoretical covariance matrix
  R_mat <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))

  for(i in 1:length(Hvec)){
    for(j in 0:(t/2-1)){
      R_mat[i,j+1] <- Rvec(dist=Hvec[i], time=j, a1, a2, a3)
    }
  }
  
  
  #Loss function (object function) for optimization
  obj <- function(a){
    Rvec <- function(dist, time, a1=a[1], a2=a[2], a3=a[3]){
      Ct <- a[1]*exp(-a[2]*abs(time))
      R <- a[3]*(1 - Ct^2)/(1-2*cos(dist)*Ct+Ct^2)^(3/2)
      return(R)
    }
    
    #theoretical covariance matrix
    R_mat <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
    for(i in 1:length(Hvec)){
      for(j in 0:(t/2-1)){
        R_mat[i,j+1] <- Rvec(dist=Hvec[i], time=j, a1=a[1], a2=a[2], a3=a[3])
      }
    }
    sum(W*(R_mom- R_mat)^2)
  }
  
  fitted_para <- nlminb(c(0,0,0),obj)$par #fitted parameters by nlminb
  #fitted_para <- nlminb(c(0,0,0),obj,lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters by nlminb
  
  fitted_para2 <- optim(c(0,0,0), obj, method = "BFGS")$par #fitted parameters
  #fitted_para <- optim(c(0,0,0), obj, method = "L-BFGS-B", lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters
  return(list(nlminb=fitted_para, optim=fitted_para2))
}

sim_size <- 1000
para <- data.frame(matrix(ncol=3, nrow=sim_size))
colnames(para) <- c("a1","a2","a3")
para_list <- list(nlminb=para, optim=para)

for(i in 1:sim_size){
  #para[i,] <- sim(a1=0.8, a2=0.1, a3=1, P=10, t=20)
  #cat("\n Iteration : ", i, "\n")
  #print(para[[1]][i,])
  
  para1 <- sim(a1=0.8, a2=0.1, a3=1, P=200, t=20)
  para_list[[1]][i,] <- para1[[1]]
  para_list[[2]][i,] <- para1[[2]]
  cat("\n Iteration : ", i, "\n")
  print(para1)
}

