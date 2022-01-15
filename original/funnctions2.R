#####################################################
##############Homogeneous and Stationary#############
#####################################################

#Spherical distance
SphDist <- function(lat1, long1, lat2, long2){
  # Great Circle Formula
  cosd <- sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(abs(long2 - long1))
  suppressWarnings(d <- acos(cosd)) 
  if(is.nan(d)) d <- 0 #acos has an issue that it returns NaN when d=1 although it actually should return 0.
  return(d)
}

#Covariance Function
Cov_fn <- function(dat1, dat2, a1, a2, a3=1){
  lat1 <- dat1[[1]]
  lat2 <- dat2[[1]]
  
  long1 <- dat1[[2]]
  long2 <- dat2[[2]]
  
  time1 <- dat1[[3]]
  time2 <- dat2[[3]]
  
  d <- SphDist(lat1, long1, lat2, long2) #spherical distance
  #Ct <- exp(-a2*abs(time1-time2))
  #R <- a1/(4*pi)*(1-Ct^2)/(1-2*d*Ct+Ct^2)^(3/2)
  
  Ct <- a1*exp(-a2*abs(time1-time2))
  R <- a3*(1 - Ct^2)/(1-2*cos(d)*Ct+Ct^2)^(3/2)
  return(R)
}

#Function to create covariance matrix
Cov_mat <- function(dat, a1, a2, a3){
  P <- nrow(dat)
  R <- matrix(nrow=P,ncol=P)
  for(i in 1:P){
    for(j in i:P){
      R[i,j] <- R[j,i] <- Cov_fn(dat[i,], dat[j,], a1, a2, a3)
    }
  }
  return(R)
}

#log-likelihood function
loglike <- function(dat, Z, a1, a2, a3, mu.vec=rep(0,nrow(dat))){
  Sigma <- R <- Cov_mat(dat, a1, a2, a3)
  d <- nrow(dat)
  
  sum_i <- sum(apply(t(Z), 1, function(x) t(x-mu.vec)%*%solve(Sigma)%*%(x-mu.vec)))
  #sum_i <- t(matrix(Z))%*%solve(Sigma)%*%matrix(Z)
  logl <- as.numeric(-.5*d*log(2*pi) - .5*log(det(Sigma)) - .5*sum_i)
  
  return(logl)
}


#Function for matrix of spherical distance
Dfun <- function(Lmat, P, t){
  # returns vector of all distance pairs
  
  Dmat <- matrix(NA,P,P)
  
  for (i in 1:P) {
    for (j in i:P) {
      Dmat[i,j] <- Dmat[j,i] <- SphDist(Lmat[i,1],Lmat[i,2],Lmat[j,1],Lmat[j,2])
    }
  }
  
  Dmat <- do.call("cbind", replicate(t, Dmat, simplify = FALSE)) #all the repeated column values of spherical distances (by time)
  Dmat <- do.call("rbind", replicate(t, Dmat, simplify = FALSE)) #all the repeated row values of spherical distances (by time)
  
  
  #Dmat2 <- Dmat
  #for(k in 1:(t-1)){
  #  Dmat2 <- cbind(Dmat, Dmat2) #all the repeated column values of spherical distances (by time)
  #}
  #Dmat3 <- Dmat2
  #for(l in 1:(t-1)){
  #  Dmat3 <- rbind(Dmat2, Dmat3) #all the repeated row values of spherical distances (by time)
  #}
  #Dmat <- Dmat3
  
  Dvec  <- Dmat[upper.tri(Dmat,diag = TRUE)] # length(which(Dvec == 0))
  
  return(Dvec)
}

#Function for MOM estimator
G_Hat <- function(Z, sDvec, Hvec, tDvec ,eps, P, t){
  
  # Create Second Moment Matrix and Vector
  
  ZZmat <- matrix(0,P*t,P*t)
  
  ZZmat <- Z%*%t(Z)
  #  for (i in 1:P) {
  #    for (j in i:P) {
  #      ZZmat[i,j] <- Z[i]*Z[j]
  #    }
  #  } # length(which(ZZmat == 0)) # (P - 1)*P/2
  
  ZZvec <- ZZmat[upper.tri(ZZmat,diag = TRUE)] # (P + 1)*P/2
  
  # Calculate Estimates using MOM Estimator
  G <- N <- matrix(nrow=length(Hvec), ncol=t/2)
  
  for(i in 0:(t/2-1)){
    s <- which(sDvec <= eps/2 & tDvec == i)
    #s <- which(sDvec <=  Hvec[2] & tDvec == i)
    N[1,i+1] <- length(s)
    G[1,i+1] <- mean(ZZvec[s])# - mean(Z)^2
    
    for (j in 2:(length(Hvec))){
      s <- which(((Hvec[j] - eps/2) <  sDvec) & (sDvec <  (Hvec[j] + eps/2)) & tDvec == i)
      #s <- which((sDvec > Hvec[j]) & (sDvec <=  Hvec[j]+eps) & tDvec == i)
      N[j,i+1] <- length(s)
      G[j,i+1] <- mean(ZZvec[s])# - mean(Z)^2
    }
  }
  
  return(list(G = G, N = N))
}

#Function to compute theoretical covariance values
Rvec <- function(dist, time, a1, a2, a3=1){
  Ct <- a1*exp(-a2*abs(time))
  R <- a3*(1 - Ct^2)/(1-2*cos(dist)*Ct+Ct^2)^(3/2)
  return(R)
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
  sum((R_mom- R_mat)^2) #without the weights
  #sum(W*(R_mom- R_mat)^2)
}


#plot for spherical distances
s_plot <- function(j){
  min_r <- min(R_mat[,j], R_fitted[,j], R_mom[,j])
  max_r <- max(R_mat[,j], R_fitted[,j], R_mom[,j])
  plot1 <- plot(Hvec, R_mat[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("t=",(j-1)), ylim=c(min_r, max_r)); 
  lines(Hvec, R_fitted[,j], col="red", lty=2); lines(Hvec, R_mom[,j], col="blue", lty=4);
  legend("topright", legend=c("True", "Fitted", "MoM"), lty=c(1,2,4), col=c("black","red","blue"))
}

#plot for time
t_plot <- function(i){
  time <- 0:(ncol(R_mat)-1)
  min_r <- min(R_mat[i,], R_fitted[i,], R_mom[i,])
  max_r <- max(R_mat[i,], R_fitted[i,], R_mom[i,])
  plot1 <- plot(time, R_mat[i,], type="l", lty=1, xlab="Time lag", ylab="Covariance", main=paste(expression(theta),"=",(i-1)), ylim=c(min_r, max_r)); 
  lines(time, R_fitted[i,], col="red", lty=2); lines(time, R_mom[i,], col="blue", lty=4);
  legend("topright", legend=c("True", "Fitted", "MoM"), lty=c(1,2,4), col=c("black","red","blue"))
}