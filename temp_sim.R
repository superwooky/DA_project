####################################################################
#########################Generate data##############################
####################################################################

#setwd("/Users/jongwookkim/Box/IU/DA")
setwd("/Users/jongwookkim/OneDrive - Indiana University/IU/DA")
source("functions.R")

sim2 <- function(year, sample_size=700 ,Dvec=Dvec){
  if(year == 1980){
    temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1980.Rds")
    temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1981.Rds")
  }else if(year == 1985){
    temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1985.Rds")
    temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1986.Rds")
  }else if(year == 1990){
    temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1990.Rds")
    temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1991.Rds")
  }else if(year == 1995){
    temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1995.Rds")
    temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1996.Rds")
  }else if(year == 2000){
    temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2000.Rds")
    temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2001.Rds")
  }else if(year == 2005){
    temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2005.Rds")
    temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2006.Rds")
  }else if(year == 2010){
    temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2010.Rds")
    temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2011.Rds")
  }
  
  #Dvec <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda02V2-Dvec.Rds") #Spherical distance vector

  #index <- sample(1:nrow(temp1), size=700, replace=FALSE)
  index <- sample(1:nrow(temp1), sample_size, replace=FALSE)
  temp1 <- temp1[index,]
  temp2 <- temp2[index,]
  
  #remove the redundancy of the column names.
  temp <- cbind(temp1, temp2[,3:14])
  colnames(temp)[15:26] <- c(13:24)
  
  
  t <- seq(1, ncol(temp)-2) #time term
  dat_lat <- rep(temp[,1], length(t))
  dat_long <- rep(temp[,2], length(t))
  dat_time <- rep(t, each=nrow(temp))
  Z <- as.vector(temp[,c(-1,-2)])
  
  dat <- data.frame(lat=dat_lat, long=dat_long, time=dat_time, Z=Z) #Create Data
  
  ####################################################################
  ######################Compute MOM estimator#########################
  ####################################################################
  
  P <- nrow(temp) #location data size
  t <- ncol(temp)-2 #time data size
  n <- P*t #total data size
  
  H <- 10         # Number of distances at which to estimate
  Hvec <- as.vector(seq(0, 1, by = 1/H))[1:H]
  eps <- 0.10     # use to bin distances in estimation
  
  #H <- 20         # Number of distances at which to estimate
  #Hvec <- as.vector(seq(0, 1.5, by = 1/H))[1:H]
  #eps <- 0.05     # use to bin distances in estimation
  
  sDvec <- Dfun(Lmat=dat[,1:2], P, t) #spherical distances
  
  tDmat <- as.matrix(dist(dat$time, method="manhattan")) #time distances
  tDvec  <-tDmat[upper.tri(tDmat,diag = TRUE)]
  
  mom_est <- G_Hat(dat$Z, sDvec, Hvec, tDvec ,eps, P, t) #MOM estimate
  R_mom <- mom_est$G
  W <- mom_est$N #weight
  
  
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
    sum((R_mom- R_mat)^2) #without weight
    #sum(W*(R_mom- R_mat)^2)
  }
  
  par1 <- nlminb(c(0,0,0),obj) #fitted parameters by nlminb
  #par1 <- nlminb(c(0,0,0),obj,lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters by nlminb
  
  par2 <- optim(c(0,0,0), obj, method = "BFGS") #fitted parameters
  #par2 <- optim(c(0,0,0), obj, method = "L-BFGS-B", lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters
  
  return(list(nlminb=par1$par, optim=par2$par, convergence1=par1$convergence, convergence2=par2$convergence))
}


Dvec <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda02V2-Dvec.Rds") #Spherical distance vector

sim_size <- 50


para <- data.frame(matrix(ncol=3, nrow=sim_size))
colnames(para) <- c("a1","a2","a3")
para_list2 <- list(nlminb=para, optim=para, convergence1=1, convergence2=1)
#para_list2$convergence1 <- para_list2$convergence2 <- 0

for(i in 1:sim_size){

  while(any(para_list2$nlminb[i,] < 0) | para_list2$convergence1==1 | is.na(para_list2[[1]][i,1])){
    para1 <- sim2(year=1980, sample_size=700)
  
    para_list2[[1]][i,] <- para1[[1]]
    para_list2[[2]][i,] <- para1[[2]]
    para_list2[[3]] <- para1[[3]]
    para_list2[[4]] <- para1[[4]]
    
    cat("\n Iteration : ", i, "\n")
    print(para1)
    print(colMeans(para_list2[[1]][1:i,]))
    #print(colMeans(para_list2[[2]][1:i,]))
    
    print(apply(para_list2[[1]][1:i,],2,sd))
    
    para1980 <- para_list2
    save(para1980,file="para1980.Rda")
  }
}


para1980
para_list2 <- para1980
#para_list2[[2]][12,] <- NA


colMeans(para1980[[1]])
apply(para1980$nlminb,2,sd)

colMeans(para1990[[1]])
apply(para1990$nlminb,2,sd)

colMeans(para2000[[1]])
apply(para2000$nlminb,2,sd)

colMeans(para2010[[1]])
apply(para2010$nlminb,2,sd)

summary(para1980[[2]])
summary(para1990[[2]])
summary(para2000[[2]])
summary(para2010[[2]])


par(mfrow=c(4,3))
#boxplot(para1980[[1]][,1])
#boxplot(para1980[[1]][,2])
#boxplot(para1980[[1]][,3])

boxplot(para1980[[2]][,1])
boxplot(para1980[[2]][,2])
boxplot(para1980[[2]][,3])

#############################

#boxplot(para1990[[1]][,1])
#boxplot(para1990[[1]][,2])
#boxplot(para1990[[1]][,3])

boxplot(para1990[[2]][,1])
boxplot(para1990[[2]][,2])
boxplot(para1990[[2]][,3])

#############################

#boxplot(para2000[[1]][,1])
#boxplot(para2000[[1]][,2])
#boxplot(para2000[[1]][,3])

boxplot(para2000[[2]][,1])
boxplot(para2000[[2]][,2])
boxplot(para2000[[2]][,3])

#############################

#boxplot(para2010[[1]][,1])
#boxplot(para2010[[1]][,2])
#boxplot(para2010[[1]][,3])

boxplot(para2010[[2]][,1])
boxplot(para2010[[2]][,2])
boxplot(para2010[[2]][,3])
par(mfrow=c(1,1))

