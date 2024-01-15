############################################################################################################################
# With no measurement error data
############################################################################################################################
rm(list=ls())
set.seed(613)
library(sp);library(fda);library(refund);library(RandomFields)
library(nimble, warn.conflicts = FALSE);library(mnormt);library(matrixcalc)
library(spatial);library(mvtnorm);library(fields);library(rdist)
library(classInt); library(INLA)
library(gstat) ;library(nlme) ; library(fda.usc)
library(fdagstat)
library(geoFKF)
library(fda)

#No measurement Error Data
load("simulation_data_no_measurement_error.RData")
xft<-Xfts%*%bmat
cov <- cbind(as.matrix(ncoord), as.matrix(xft))
cov <- data.frame(cov)

g <- fstat(NULL,  vName = "Y", data.frame(cov[ind,]), data.frame(Y_train))
g <- estimateDrift("~.", g, Intercept = TRUE)
g <- fvariogram("~X1+X2", g, Nlags = 100, LagMax = 1, ArgStep = 1, comments=FALSE)
g <- fitVariograms(g, model=vgm("Gau"))
g <- addCovariance(g, 'omni')

pt <- proc.time()
forecasts <- predictFstat(g, .newCoordinates = data.frame(cov[-ind,]), .what = "Y", .type = "UK")
time <- proc.time()-pt
print(time)

#90*300 test data
par(mfrow=c(4,4))
for (i in 1:16){
  matplot(forecasts$Forecast[,i], type="l", col="blue", xlab="Time(days)", ylab="Value", main="Simulation data-fdagstat", ylim=c(-3, 4))
  matplot(Y_test[,i], type="l", col="red", add = TRUE, ylim = c(-3,4))
}

sqrt(mean((forecasts$Forecast-Y_test)^2, na.rm=T))
############################################################################################################################
# With measurement error data
############################################################################################################################
rm(list=ls())
set.seed(613)

library(sp);library(fda);library(refund);library(RandomFields)
library(nimble, warn.conflicts = FALSE);library(mnormt);library(matrixcalc)
library(spatial);library(mvtnorm);library(fields);library(rdist)
library(classInt); library(INLA)
library(gstat) ;library(nlme) ; library(fda.usc)
library(fdagstat)
library(geoFKF)
library(fda)

set.seed(613)
load("simulation_data_tau1.RData")
xft<-Xfts%*%bmat
cov <- cbind(as.matrix(ncoord), as.matrix(xft))
cov <- data.frame(cov)

g <- fstat(NULL,  vName = "Y", data.frame(cov[ind,]), data.frame(Y_train))
g <- estimateDrift("~.", g, Intercept = TRUE)
g <- fvariogram("~X1+X2", g, Nlags = 100, LagMax = 1, ArgStep = 1, comments=FALSE)
g <- fitVariograms(g, model=vgm("Gau"))
g <- addCovariance(g, 'omni')

pt <- proc.time()
forecasts <- predictFstat(g, .newCoordinates = data.frame(cov[-ind,]), .what = "Y", .type = "UK")
time <- proc.time()-pt
print(time)

par(mfrow=c(4,4))
for (i in 1:16){
  matplot(forecasts$Forecast[,i], type="l", col="blue", xlab="Time(days)", ylab="Value", main="Simulation data-fdagstat", ylim=c(-3, 4))
  matplot(Y_test[,i], type="l", col="red", add = TRUE, ylim = c(-3,4))
}

sqrt(mean((forecasts$Forecast-Y_test)^2, na.rm=T))