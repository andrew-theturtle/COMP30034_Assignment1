# import required matrix
library(plot.matrix)
library(abind)

# import dataset
X = as.matrix(read.csv("../data/data.csv"))
TC = as.matrix(read.csv("../data/TC.csv"))
SM = as.matrix(read.csv("../data/SM.csv"))

# estimate A(retrieval of SM) and D(TC) using least square solution
A_LSR = solve(t(TC)%*%TC)%*%t(TC)%*%X
D_LSR = X%*%t(A_LSR)

# plot six retrieved spatial and temporal sources
par(mfrow=c(3,2), mar=c(2,3,2,3))
for (i in 1:6) {
  # apply abs on retrieved coefficient values A, making it easer to visualize
  spatial = abs(A_LSR[i,])
  spatial = matrix(spatial, 21, 21)
  plot(cor(spatial), main=paste("Retrieved SM ", i), border=NA, ylab='', mar=c(2,4,2,4))
  temporal = D_LSR[,i]
  plot(temporal, type='l', ylab='', xlab='', main=paste('Retrieved TM ', i), cex=0.5)
}

par(mfrow=c(2,1), mar=c(4,4,2,2))
# scatter plot between 3rd column of D and 30th column of X
plot(x=D_LSR[,3], y=X[,30], xlab='3rd column of D', ylab='30th column of X')
# scatter plot between 3rd column of D and 30th column of X
plot(x=D_LSR[,4], y=X[,30], xlab='4th column of D', ylab='30th column of X')
# SM[3,30] = 1, SM[4,30] = 0, only TC3 contribute in the creation of 30th data element
# and it also explain the perfect linear relationship between these two columns.

# ridge regression
lambda = 0.7 * 441
A_RR = solve(t(TC)%*%TC + lambda*diag(6))%*%t(TC)%*%X
D_RR = X%*%t(A_RR)

# correlation vectors retaining only maximum absolute
c_TLSR = abs(diag(cor(TC, D_LSR)))
c_TRR = abs(diag(cor(TC, D_RR)))
sum(c_TRR)>sum(c_TLSR)

# ridge regression with lambda=1000
lambda = 1000 * 441
A_RR_1000 = solve(t(TC)%*%TC + lambda*diag(6))%*%t(TC)%*%X
cbind("a_RR" = A_RR_1000[,1], "a_LSR" = A_LSR[,1])
# We indeed observe all value in a_RR are being pushed close to zero, comparing
# with a_LSR

# Code for LR given in assignment specification
lasso_regression <- function(TC, rho, N, nsrcs, x1, x2) {
  step = 1/(norm(TC %*% t(TC)) * 1.1)
  thr = rho*N*step
  Ao = matrix(0, nsrcs, 1)
  A = matrix(0, nsrcs, 1)
  Alr = matrix(0, nsrcs, x1*x2)
  
  for (k in 1:(x1*x2)) {
    A = Ao+step*(t(TC) %*% (X[,k]-(TC%*%Ao)))
    A = (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
    
    for (i in 1:10) {
      Ao = A
      A = Ao+step * (t(TC)%*%(X[,k]-(TC%*%Ao)))
      A = (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
    }
    Alr[,k] = A
  }
  return(Alr)
}

# calculate MSEs for lasso regression with range of different
# rho values: 21 values between 0 and 1 with interval of 0.05.
calculate_MSEs <- function(X, TC) {
  rhos = seq(0,1,0.05)
  N = 240; V=441; nsrcs=6; x1=21; x2=21
  MSEs = c()
  for (i in 1:length(rhos)) {
    A_LR = lasso_regression(TC, rhos[i], N=N, nsrcs=6, x1=x1, x2=x2)
    D_LR = X%*%t(A_LR)
    MSE = sum(sum((X-D_LR%*%A_LR)^2))/(N*V)
    MSEs = c(MSEs, MSE)
  } 
  return(MSEs)
}

# generate new standardized dataset X 
generate_X <- function(TC, SM, seed) {
  set.seed(seed)
  # generate white Gaussian noise
  temporal_noise = rnorm(240*6, mean=0, sd=sqrt(0.25))
  temporal_noise = matrix(temporal_noise, nrow=240, ncol=6)
  
  spatial_noise = rnorm(6*441, mean=0, sd=sqrt(0.015))
  spatial_noise = matrix(spatial_noise, nrow=6, ncol=441)
  
  X = (TC + temporal_noise) %*% (SM + spatial_noise)
  # standardized dataset
  X = scale(X)
  return(X)
}

# seed=156 is the seed used in question 1 to generate X
SEED = c(156,157,158,159,160,161,162,163,164,165)
MSE_all = c()
# perform 10 realizations
for (i in 1:10) {
  X = generate_X(TC, SM, SEED[i])
  MSEs = calculate_MSEs(X, TC)
  MSE_all = abind(MSE_all, MSEs, along=2)
}

# average of MSE over 10 realizations
MSE_avg = rowMeans(MSE_all)
par(mfrow=c(1,1))
plot(x=seq(0,1,0.05), y=MSE_avg, col='red', xlab=expression(paste("",rho)))
# at rho=0.3, MSE value start to increase again. 
# at rho=0.25, MSE value is at its minimum

# using generated data from question 1
X = as.matrix(read.csv("../data/data.csv"))
# estimate LR parameters for rho=0.25
A_LR = lasso_regression(TC, 0.25, N=240, nsrcs=6, x1=21, x2=21)
D_LR = X%*%t(A_LR)

# correlation vectors retaining only maximum absolute value
c_TRR = abs(diag(cor(TC, D_RR)))
c_SRR = abs(diag(cor(t(SM), t(A_RR))))
c_TLR = abs(diag(cor(TC, D_LR)))
c_SLR = abs(diag(cor(t(SM), t(A_LR))))
# verify that our choice of regularization term rho is correct
sum(c_TLR)>sum(c_TRR)
sum(c_SLR)>sum(c_SRR)

# plot D and A side by side
par(mfrow=c(3,4), mar=c(2,3,3,3))
for (i in 1:6) {
  plot(abs(matrix(A_RR[i,],21,21)), border=NA, main=paste("A_RR ", i))
  plot(D_RR[,i], type='l', main=paste("D_RR ", i))
  plot(abs(matrix(A_LR[i,],21,21)), border=NA, main=paste("A_LR ", i))
  plot(D_LR[,i], type='l', main=paste("D_LR ", i))
}

# principal component regression
s = svd(TC)
PC = s$u
eigenvalues = s$d

# plot eigenvalues
par(mfrow=c(1,1), mar=c(4,4,3,3))
plot(eigenvalues, type='l', xlab="PC", ylab='Eigenvalue') 
# The eigenvalue of the 6 PC is the lowest

# plot the regressors in Z and source TC
par(mfrow=c(3,2), mar=c(3,3,3,3))
for (i in 1:3) {
  plot(TC[,i], type='l', main=paste("TC ", i))
  plot(PC[,i], type='l', main=paste("PC ", i))
}

# estimate PCR parameters for rho=0.001
A_PCR = lasso_regression(PC, 0.001, N=240, nsrcs=6, x1=21, x2=21)
D_PCR = X%*%t(A_PCR)
# plot D_PCR and A_PCR side by side
par(mfrow=c(3,2), mar=c(2,3,3,3))
for (i in 1:6) {
  plot(abs(matrix(A_PCR[i,],21,21)), border=NA, main=paste("A_PCR ", i))
  plot(D_PCR[,i], type='l', main=paste("D_PCR ", i))
  #plot(abs(matrix(A_RR[i,],21,21)), border=NA, main=paste("A_LR ", i))
}
