rm(list=ls())
library(fdagstat)
library(geoFKF)
library(fda)

library(sp);library(fda);library(refund);library(RandomFields)
library(nimble, warn.conflicts = FALSE);library(mnormt)
library(spatial);library(mvtnorm);library(fields);library(rdist)
library(classInt); library(INLA)
library(gstat) ;library(nlme) ; library(fda.usc)

load("japan_pm25_dat_align.Rdata")
load("japan_coord_dat_align.Rdata")
load("japan_no2_dat_align.RData")
load("japan_align_train_index.RData")

dat<-as.matrix(log(jap_pm + 0.1))
cov <- cbind(m_coord,no2)

#Define fstat object for fdagstat
g <- fstat(g=NULL,  vName = "PM25", Coordinates=data.frame(cov[ind,]),Functions=data.frame(dat[,ind]))
g <- estimateDrift("~.",g, Intercept = TRUE)
g <- fvariogram("~lat+lon", g, Nlags = 100, LagMax = 1, ArgStep = 1, comments=FALSE)
plotVariogram(g)
g <- fitVariograms(g, model=vgm("Gau"))
g <- addCovariance(g)

forecasts <- predictFstat(g, .newCoordinates = data.frame(cov[-ind,]), .what = "PM25", .type = "UK")

time <- proc.time()-pt
print(time)

#Draw Kriging of fdagstat
par(mfrow=c(4,4))
for (i in 1:16){
  matplot(forecasts$Forecast[,i], type="l", col="blue", xlab="Time(days)", ylab="Value", main=paste(i,"th place Japan, PM2.5 Fdagstat"))
  matplot(dat[,-ind][,i], type="l", col="red", add = TRUE)
}

#MSE
sqrt(mean((forecasts$Forecast-dat[,-ind])^2, na.rm=T))

