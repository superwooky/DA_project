#setwd("/Users/jongwookkim/Box/IU/DA")
setwd("/Users/jongwookkim/OneDrive - Indiana University/IU/DA")
source("functions.R")

#####################################################################
####################check proper parameter values####################
#####################################################################

a1 <- 0.3 #should be a1<1 #a parameter for spatial term
a2 <- 0.02 #a parameter for time
a3 <- 1 #scale parameter

a1 <- 0.8 #should be a1<1 #a parameter for spatial term
a2 <- 0.1 #a parameter for time
a3 <- 1 #scale parameter

distance <- 0
time <- 0
Ct <- a1*exp(-a2*time)
variance <- (1 - Ct^2)/(1-2*cos(distance)*Ct+Ct^2)^(3/2)

#check the spatial term
par(mfrow=c(2,2))

distance <- seq(0,pi,length=100)
time <- 0
Ct <- a1*exp(-a2*time)
R1 <- a3*(1 - Ct^2)/(1-2*cos(distance)*Ct+Ct^2)^(3/2)/variance
plot(distance,R1,type="l", ylab="Correlation", xlab="Spherical Distance", main=bquote(p[1]==.(a1) ~~~ p[2]==.(a2) ~~~ p[3]==.(a3)))

#check the time term
distance <- 1
time <- seq(0,50, length=50)
Ct <- a1*exp(-a2*time)
R1 <- a3*(1 - Ct^2)/(1-2*distance*Ct+Ct^2)^(3/2)/variance
plot(time,R1,type="l", ylab="Correlation", xlab="Time Lag", main=bquote(p[1]==.(a1) ~~~ p[2]==.(a2) ~~~ p[3]==.(a3))) 

par(mfrow=c(1,1))

####################################################################
#########################Generate data##############################
####################################################################

#a1 <- 0.75 #should be a1<1 #a parameter for spatial term
#a2 <- 0.05 #a parameter for time
#a3 <- 1 #scale parameter


P <- 200 #location data size
t <- 20 #time data size
n <- P*t #total data size

#P <- 10 #location data size
#t <- 5 #time data size
#n <- P*t #total data size

set.seed(117)
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


####################################################################
######################Compute MOM estimator#########################
####################################################################


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
W <- mom_est$N #weight
#W <- 1

#theoretical covariance matrix
R_mat <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
for(i in 1:length(Hvec)){
  for(j in 0:(t/2-1)){
    R_mat[i,j+1] <- Rvec(dist=Hvec[i], time=j, a1, a2, a3)
  }
}



fitted_para <- nlminb(c(0,0,0),obj)$par #fitted parameters by nlminb
#fitted_para <- nlminb(c(0,0,0),obj,lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters by nlminb

#fitted_para <- optim(c(0,0,0), obj, method = "BFGS")$par #fitted parameters
#fitted_para <- optim(c(0,0,0), obj, method = "L-BFGS-B", lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters

nlminb(c(0,0,0),obj)$par
#nlminb(c(0,0,0),obj,lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par
optim(c(0,0,0), obj, method = "BFGS")$par
#optim(c(0,0,0), obj, method = "L-BFGS-B", lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters


#fitted covariance matrix
R_fitted <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
for(i in 1:length(Hvec)){
  for(j in 0:(t/2-1)){
    R_fitted[i,j+1] <- Rvec(dist=Hvec[i], time=j, fitted_para[1], fitted_para[2], fitted_para[3])
  }
}

#sum(W*(R_mom- R_mat)^2) #weighted sum of least square
sum((R_mom- R_mat)^2) #weighted sum of least square

nlminb(c(0,0,0),obj)$par #fitted values of the parameters
a1;a2;a3 #true values of the parameters

optim(c(0,0,0), obj, method = "BFGS")


#R_fitted/variance  #fitted value
#R_mom/variance     #mom estimate
#R_mat/variance     #true value


#plot for spherical distances
s_plot <- function(j){
  min_r <- min(R_mat[,j], R_fitted[,j], R_mom[,j])
  max_r <- max(R_mat[,j], R_fitted[,j], R_mom[,j])
  plot1 <- plot(Hvec, R_mat[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("h=",(j-1)), ylim=c(0, 50)); 
  lines(Hvec, R_fitted[,j], col="red", lty=2); lines(Hvec, R_mom[,j], col="blue", lty=4)
  if(j==1){
    legend("topright", legend=c("True", "Fitted", "MoM"), lty=c(1,2,4), col=c("black","red","blue"))
  }
}

#plot for time
t_plot <- function(i){
  time <- 0:(ncol(R_mat)-1)
  min_r <- min(R_mat[i,], R_fitted[i,], R_mom[i,])
  max_r <- max(R_mat[i,], R_fitted[i,], R_mom[i,])
  plot1 <- plot(time, R_mat[i,], type="l", lty=1, xlab="Time Lag", ylab="Covariance", main=bquote(psi== .(Hvec[i])), ylim=c(0,50)); 
  lines(time, R_fitted[i,], col="red", lty=2); lines(time, R_mom[i,], col="blue", lty=4)
  if(i==1){
    legend("topright", legend=c("True", "Fitted", "MoM"), lty=c(1,2,4), col=c("black","red","blue"))
  }
}

par(mfrow=c(2,3))

s_plot(1)
s_plot(2)
s_plot(3)
#s_plot(4)
#s_plot(5)
s_plot(6)
#s_plot(7)
s_plot(8)
#s_plot(9)
s_plot(10)


t_plot(1)
t_plot(2)
t_plot(3)
t_plot(4) 
t_plot(5)
#t_plot(6)
t_plot(7)
#t_plot(8)
#t_plot(9)
#t_plot(10)

par(mfrow=c(1,1))
