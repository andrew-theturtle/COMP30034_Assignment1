#import required library
library(plot.matrix)
library(abind)
library(reshape2)
library(ggplot2)

AV = c(0,20,0,0,0,0)
VI = c(30,45,60,40,40,40)
Duration = c(15,20,25,15,20,25)
# create empty matrix to hold the data
TC = matrix(nrow=240)
# generate one temporal source given Arrival, Increment, and Duration
generate_TC <- function(Arrival, Increment, Duration) {
  x = rep(0,240)
  i = Arrival + 1
  while (i <= 240) {
    j = i + Duration - 1
    x[i:j] = 1
    i = i + Increment
  } 
  return(x)
}
# construct a matrix TC of six temporal sources
for (i in 1:length(AV)) {
  single_TC = generate_TC(AV[i], VI[i], Duration[i])
  TC = cbind(TC, single_TC)
}
TC = TC[,-1]
colnames(TC) = NULL
# standardize each temporal source to have mean 0 and std 1
TC_scaled = scale(TC)
# save TC to file
write.csv(TC_scaled, file="../data/TC.csv", row.names = FALSE)

# plot all TCs as six subplots
par(mfrow=c(2,3))
for (i in 1:ncol(TC_scaled)) {
  data = TC_scaled[,i]
  plot(data, type='l', xlab='', ylab='', main=sprintf("TC %s", i))
}

# normalization (divide by l2-norm)
l2_norm <- function(x) sqrt(sum(x^2))
TC_norm_value = apply(TC, 2, l2_norm)
TC_norm = sweep(TC, 2, TC_norm_value, '/')

# plot all normed TCs
par(mfrow=c(2,3))
for (i in 1:ncol(TC_norm)) {
  data = TC_norm[,i]
  plot(data, type='l', xlab='', ylab='', main=sprintf("TC %s", i))
}

par(mfrow=c(1,3))
plot(TC[,1], type='l', xlab='', ylab='', main=paste("original TC ", 1))
plot(TC_scaled[,1], type='l', xlab='', ylab='', main=paste("standardized TC ", 1))
plot(TC_norm[,1], type='l', xlab='', ylab='', main=paste("normalized TC ", 1))

# generate correlation matrix
CM = cor(TC)
par(mfrow=c(1,1))
plot(CM, border=NA)

# generate SMs
tmpSM = c()
vertical = c(c(2,6), c(2,6), c(8,13), c(8,13), c(15,19), c(15,19))
horizontal = c(c(2,6), c(15,19), c(2,6), c(15,19), c(2,6), c(15,19))
i = 1
while (i <= 12) {
  tmp = matrix(nrow=21, ncol=21, 0)
  tmp[vertical[i]:vertical[i+1],horizontal[i]:horizontal[i+1]] = 1
  tmpSM = abind(tmpSM, tmp, along=3)
  i = i + 2
}

# plot these SMs in six subplot
par(mfrow=c(2,3), mar=c(2,4,2,4))
for (i in 1:6) {
  data = tmpSM[,,i]
  plot(data, border=NA, xlab='', ylab='', main=paste("SM ", i))
}

# reshape array into 2D matrix with shape 6*441
SM = t(matrix(tmpSM, 21*21, 6))
# save to file
write.csv(SM, file='../data/SM.csv', row.names = FALSE)

# generate white Gaussian noise
set.seed(156)
temporal_noise = rnorm(240*6, mean=0, sd=sqrt(0.25))
temporal_noise = matrix(temporal_noise, nrow=240, ncol=6)

spatial_noise = rnorm(6*441, mean=0, sd=sqrt(0.015))
spatial_noise = matrix(spatial_noise, nrow=6, ncol=441)

# correlation matrix for each noise type
par(mfrow=c(1,2))
plot(cor(temporal_noise), border=NA, xlab="", ylab="", main="Temporal noise")
plot(cor(t(spatial_noise)), border=NA, xlab="", ylab="", main="Spatial noise")

#histogram of both noise sources
par(mfrow=c(1,2), mar=c(2,2,2,2))
hist(as.vector(temporal_noise), main="Temporal noise", freq=FALSE)
lines(density(as.vector(temporal_noise)),col='red')
hist(as.vector(spatial_noise), main="Spatial noise", freq=FALSE)
lines(density(as.vector(spatial_noise)),col='red')

# p-value = 0.3973 > 0.05, can not reject the null hypothesis, 
# implying temporal noise follow normal distribution
shapiro.test(as.vector(temporal_noise))
# p-value = 0.9198 > 0.05, can not reject the null hypothesis, 
# implying spatial noise also follow normal distribution
shapiro.test(as.vector(spatial_noise))

# product of temporal and spatial source
temporal_spatial_noise = temporal_noise %*% spatial_noise
par(mfrow=c(1,1))
plot(cor(temporal_spatial_noise), border=NA)

# generate a synthetic dataset
X = (TC + temporal_noise) %*% (SM + spatial_noise)
# plot 100 randomly selected time-series from X
X_df = data.frame(X)
sample_col_idx = sample(ncol(X_df), 100)
X_df = X_df[,sample_col_idx]
X_df['index']=1:240

# convert data into long format to use ggplot2
X_df_long = melt(X_df, id.vars = c('index'))

par(mfrow=c(1,1), mar=c(2,2,2,2))
ggplot(X_df_long, aes(x=index, y=value, group=variable, color=variable)) + 
    geom_line()

# plot variance of all 441 variable
par(mfrow=c(1,1), mar=c(4,4,2,2))
plot(diag(var(X)), main='Variance of all 441 variables', xlab="Variable", ylab="Variance")

# standardized dataset
X = scale(X)
# save the dataset
write.csv(X, file="../data/data.csv", row.names=FALSE)
