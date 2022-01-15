####################################################################
#########################Generate data##############################
####################################################################
setwd("/Users/jongwookkim/Downloads/git2/DA_project")
source("functions.R")

temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1980.Rds")
temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1981.Rds")

temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1990.Rds")
temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-1991.Rds")

temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2000.Rds")
temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2001.Rds")

temp1 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2010.Rds")
temp2 <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-2011.Rds")

#Dvec <- readRDS("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda02V2-Dvec.Rds")                  #Spherical distance vector
#Lmat <- readRDS("~/Box/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-Lmat.Rds")       #Location matrix
load("/Users/jongwookkim/Downloads/git2/DA_project/Dmat.Rda") #Spherical distance matrix
sDmat <- Dmat

#sample data
#set.seed(115)

index <- sample(1:nrow(temp1), size=5, replace=FALSE)
temp1 <- temp1[index,]
temp2 <- temp2[index,]

nrow(temp1)


#remove the redundancy of the column names.
temp <- cbind(temp1, temp2[,3:14])
colnames(temp)[15:26] <- c(13:24)


#create a spherical distance maxtrix
P <- nrow(temp)
sDmat <- matrix(NA,P,P)
for (i in 1:P) {
  for (j in i:P) {
    sDmat[i,j] <- sDmat[j,i] <- SphDist(temp[i,1],temp[i,2],temp[j,1],temp[j,2])
  }
}



t <- seq(1, ncol(temp)-2) #time term
#t <- 1:5
dat_lat <- rep(temp[,1], length(t))
dat_long <- rep(temp[,2], length(t))
dat_time <- rep(t, each=nrow(temp))
Z <- as.vector(temp[,c(-1,-2)])

dat <- data.frame(lat=dat_lat, long=dat_long, time=dat_time, Z=Z) #Create Data
dat$index <- rep(1:nrow(temp), length(t)) #Create location index


#nrow(temp)
dim(dat)
#dat[1:10,]
#as.vector(temp[,c(-1,-2)])[1:10]
#temp[1:10,c(-1,-2)]

choose(nrow(dat),2) + nrow(dat)
####################################################################
######################Compute MOM estimator#########################
####################################################################


P <- nrow(temp) #location data size
t <- ncol(temp)-2 #time data size
n <- P*t #total data size

H <- 10         # Number of distances at which to estimate
Hvec <- as.vector(seq(0, 1, by = 1/H))[1:H]
eps <- 0.10     # use to bin distances in estimation

#Use a C function

dyn.load("C_code.so")
is.loaded('G_hat3')
dyn.unload("C_code.so")



mom_est3 <- G_hat3(dat=dat, sDmat, Hvec=Hvec, eps=eps, P=P, t=t/2)
R_mom <- matrix(mom_est3[[1]]/mom_est3[[2]], nrow=length(Hvec)) #MOM estimate
W <- matrix(mom_est3[[2]], nrow=length(Hvec)) #weight





par1 <- nlminb(c(0,0,0),obj)$par #fitted parameters by nlminb
#fitted_para <- nlminb(c(0,0,0),obj,lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters by nlminb

par2 <- optim(c(0,0,0), obj, method = "BFGS")$par #fitted parameters
#fitted_para <- optim(c(0,0,0), obj, method = "L-BFGS-B", lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters

#fitted covariance matrix

fitted <- function(fitted_para){
  R_fitted <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
  for(i in 1:length(Hvec)){
    for(j in 0:(t/2-1)){
      R_fitted[i,j+1] <- Rvec(dist=Hvec[i], time=j, fitted_para[1], fitted_para[2], fitted_para[3])
    }
  }
  return(R_fitted)
}


R_fitted1 <- fitted(par1)
R_fitted2 <- fitted(par2)

#####
#par1980_1 <- par1
#R_fitted1980_1 <- R_fitted1
#R_mom1980_1 <-R_mom

#par1980_2 <- par2
#R_fitted1980_2 <- R_fitted2
#R_mom1980_2 <-R_mom

#####
#par1990_1 <- par1
#R_fitted1990_1 <- R_fitted1
#R_mom1990_1 <-R_mom

#par1990_2 <- par2
#R_fitted1990_2 <- R_fitted2
#R_mom1990_2 <-R_mom

#####
#par2000_1 <- par1
#R_fitted2000_1 <- R_fitted1
#R_mom2000_1 <-R_mom

#par2000_2 <- par2
#R_fitted2000_2 <- R_fitted2
#R_mom2000_2 <-R_mom

####
par2010_1 <- par1
R_fitted2010_1 <- R_fitted1
R_mom2010_1 <-R_mom

par2010_2 <- par2
R_fitted2010_2 <- R_fitted2
R_mom2010_2 <-R_mom

optim_par <- list(par1980_1, par1990_1, par2000_1, par2010_1)
nlminb_par <- list(par1980_2, par1990_2, par2000_2, par2010_2)
mom_list <- list(R_mom1980_1, R_mom1990_1, R_mom2000_1, R_mom2010_1)
optim_fitted <- list(R_fitted1980_1, R_fitted1990_1, R_fitted2000_1, R_fitted2010_1)
nlminb_fitted <- list(R_fitted1980_2, R_fitted1990_2, R_fitted2000_2, R_fitted2010_2)
names(mom_list) <- c(1980, 1990, 2000, 2010)
names(optim_fitted) <- c(1980, 1990, 2000, 2010)
names(nlminb_fitted) <- c(1980, 1990, 2000, 2010)

nlminb(c(0,0,0),obj)$par #fitted values of the parameters
nlminb(c(0,0,0),obj,lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par
optim(c(0,0,0), obj, method = "BFGS")$par
optim(c(0,0,0), obj, method = "L-BFGS-B", lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par


#par2 <- fitted_para
#R_fitted2 <- R_fitted
#R_mom2 <-R_mom

#R_fitted/variance  #fitted value
#R_mom/variance     #mom estimate
#R_mat/variance     #true value


#plot for spherical distances
s_plot <- function(j, R_fitted){
  min_r <- min(R_fitted[,j], R_mom[,j])
  max_r <- max(R_fitted[,j], R_mom[,j])
  #plot1 <- plot(Hvec, R_fitted[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("t=",(j-1)), ylim=c(0, 1.5));
  plot1 <- plot(Hvec, R_fitted[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("t=",(j-1)), ylim=c(0, 2)); 
  lines(Hvec, R_mom[,j], col="blue", lty=2);
  legend("topright", legend=c("Fitted", "MOM"), lty=c(1,2), col=c("black","blue"))
}

#plot for time
t_plot <- function(i, R_fitted){
  time <- 0:(ncol(R_fitted)-1)
  min_r <- min(R_fitted[i,], R_mom[i,])
  max_r <- max(R_fitted[i,], R_mom[i,])
  #plot1 <- plot(time, R_fitted[i,], type="l", lty=1, xlab="Time leg", ylab="Covariance", main=paste(expression(theta),"=",(i-1)), ylim=c(0, 1.5));
  plot1 <- plot(time, R_fitted[i,], type="l", lty=1, xlab="Time leg", ylab="Covariance", main=paste(expression(theta),"=",(i-1)), ylim=c(0, 2));
  lines(time, R_mom[i,], col="blue", lty=2);
  legend("topright", legend=c("Fitted", "MOM"), lty=c(1,2), col=c("black","blue"))
}

s_plot(1, R_fitted1)
s_plot(2, R_fitted1)
s_plot(3, R_fitted1)
s_plot(4, R_fitted1)
s_plot(5, R_fitted1)
s_plot(6, R_fitted1)
s_plot(7, R_fitted1)
s_plot(8, R_fitted1)
s_plot(9, R_fitted1)
s_plot(10, R_fitted1)


t_plot(1, R_fitted1)
t_plot(2, R_fitted1)
t_plot(3, R_fitted1)
t_plot(4, R_fitted1)
t_plot(5, R_fitted1)
t_plot(6, R_fitted1)
t_plot(7, R_fitted1)
t_plot(8, R_fitted1)
t_plot(9, R_fitted1)
t_plot(10, R_fitted1)

#plot for spherical distances
s_plot2 <- function(j, R_mom1 ,R_fitted1, R_fitted2){
  min_r <- min(R_fitted1[,j], R_mom1[,j])
  max_r <- max(R_fitted1[,j], R_mom1[,j])
  #plot1 <- plot(Hvec, R_fitted[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("t=",(j-1)), ylim=c(0, 1.5));
  plot1 <- plot(Hvec, R_mom1[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("t=",(j-1)), ylim=c(0, 2)); 
  lines(Hvec, R_fitted1[,j], col="blue" ,lty=2);
  lines(Hvec, R_fitted2[,j], col="red", lty=3);
  #lines(Hvec, R_mom2[,j], col="red", lty=2);
  #legend("topright", legend=c("Fitted2009", "MOM2009","Fitted1980", "MOM1980"), lty=c(1,2,1,2), col=c("black","black","red","red"))
  legend("topright", legend=c("MOM","optim", "nlminb"), lty=1:3, col=c("black","blue","red"))
}

#plot for time
t_plot2 <- function(i, R_mom1 ,R_fitted1, R_fitted2){
  time <- 0:(ncol(R_fitted1)-1)
  min_r <- min(R_fitted1[i,], R_mom1[i,])
  max_r <- max(R_fitted1[i,], R_mom1[i,])
  #plot1 <- plot(time, R_fitted[i,], type="l", lty=1, xlab="Time leg", ylab="Covariance", main=paste(expression(theta),"=",(i-1)), ylim=c(0, 1.5));
  plot1 <- plot(time, R_mom1[i,], type="l", lty=1, xlab="Time leg", ylab="Covariance", main=paste(expression(theta),"=",(i-1)), ylim=c(0, 2));
  lines(time, R_fitted1[i,], col="blue", lty=2);
  lines(time, R_fitted2[i,], col="red", lty=3);
  #lines(time, R_mom2[i,], col="red", lty=2);
  #legend("topright", legend=c("Fitted2009", "MOM2009","Fitted1980", "MOM1980"), lty=c(1,2,1,2), col=c("black","black","red","red"))
  legend("topright", legend=c("MOM","optim", "nlminb"), lty=1:3, col=c("black","blue","red"))
}


n=4
s_plot2(1, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(2, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(3, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(4, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(5, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(6, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(7, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(8, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(9, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
s_plot2(10, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])


t_plot2(1, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(2, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(3, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(4, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(5, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(6, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(7, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(8, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(9, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])
t_plot2(10, mom_list[[n]], optim_fitted[[n]], nlminb_fitted[[n]])



#plot for spherical distances
s_plot3 <- function(j, R_fitted){
  #plot1 <- plot(Hvec, R_fitted[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("t=",(j-1)), ylim=c(0, 1.5));
  plot1 <- plot(Hvec, R_fitted[[1]][,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("t=",(j-1)), ylim=c(0, 2)); 
  lines(Hvec, R_fitted[[2]][,j], col="blue" ,lty=2);
  lines(Hvec, R_fitted[[3]][,j], col="red", lty=3);
  lines(Hvec, R_fitted[[4]][,j], col="green", lty=4);
  legend("topright", legend=c("1980","1990", "2000","2010"), lty=1:4, col=c("black","blue","red","green"))
}

#plot for time
t_plot3 <- function(i, R_fitted){
  time <- 0:(ncol(R_fitted[[1]])-1)
  #plot1 <- plot(time, R_fitted[i,], type="l", lty=1, xlab="Time leg", ylab="Covariance", main=paste(expression(theta),"=",(i-1)), ylim=c(0, 1.5));
  plot1 <- plot(time, R_fitted[[1]][i,], type="l", lty=1, xlab="Time leg", ylab="Covariance", main=paste(expression(theta),"=",(i-1)), ylim=c(0, 2));
  lines(time, R_fitted[[2]][i,], col="blue", lty=2);
  lines(time, R_fitted[[3]][i,], col="red", lty=3);
  lines(time, R_fitted[[4]][i,], col="green", lty=4);
  legend("topright", legend=c("1980","1990", "2000","2010"), lty=1:4, col=c("black","blue","red","green"))
}


s_plot3(1, mom_list)
s_plot3(2, mom_list)
s_plot3(3, mom_list)
s_plot3(4, mom_list)
s_plot3(5, mom_list)
s_plot3(6, mom_list)
s_plot3(7, mom_list)
s_plot3(8, mom_list)
s_plot3(9, mom_list)
s_plot3(10, mom_list)

t_plot3(1, mom_list)
t_plot3(2, mom_list)
t_plot3(3, mom_list)
t_plot3(4, mom_list)
t_plot3(5, mom_list)
t_plot3(6, mom_list)
t_plot3(7, mom_list)
t_plot3(8, mom_list)
t_plot3(9, mom_list)
t_plot3(10, mom_list)


s_plot3(1, nlminb_fitted)
s_plot3(2, nlminb_fitted)
s_plot3(3, nlminb_fitted)
s_plot3(4, nlminb_fitted)
s_plot3(5, nlminb_fitted)
s_plot3(6, nlminb_fitted)
s_plot3(7, nlminb_fitted)
s_plot3(8, nlminb_fitted)
s_plot3(9, nlminb_fitted)
s_plot3(10, nlminb_fitted)

t_plot3(1, nlminb_fitted)
t_plot3(2, nlminb_fitted)
t_plot3(3, nlminb_fitted)
t_plot3(4, nlminb_fitted)
t_plot3(5, nlminb_fitted)
t_plot3(6, nlminb_fitted)
t_plot3(7, nlminb_fitted)
t_plot3(8, nlminb_fitted)
t_plot3(9, nlminb_fitted)
t_plot3(10, nlminb_fitted)


s_plot3(1, optim_fitted)
s_plot3(2, optim_fitted)
s_plot3(3, optim_fitted)
s_plot3(4, optim_fitted)
s_plot3(5, optim_fitted)
s_plot3(6, optim_fitted)
s_plot3(7, optim_fitted)
s_plot3(8, optim_fitted)
s_plot3(9, optim_fitted)
s_plot3(10, optim_fitted)

t_plot3(1, optim_fitted)
t_plot3(2, optim_fitted)
t_plot3(3, optim_fitted)
t_plot3(4, optim_fitted)
t_plot3(5, optim_fitted)
t_plot3(6, optim_fitted)
t_plot3(7, optim_fitted)
t_plot3(8, optim_fitted)
t_plot3(9, optim_fitted)
t_plot3(10, optim_fitted)

par1980_1;par1980_2
par1990_1;par1990_2
par2000_1;par2000_2
par2010_2;par2010_2

#save(Dmat,file="Dmat.Rda")