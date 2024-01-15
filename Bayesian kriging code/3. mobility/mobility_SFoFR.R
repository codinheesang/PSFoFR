rm(list=ls())
## install and load packages
library(sp);library(fda);library(refund);library(RandomFields)
library(nimble, warn.conflicts = FALSE);library(mnormt)
library(spatial);library(mvtnorm);library(fields);library(rdist)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(classInt); library(INLA)
library(gstat) ;library(nlme) ; library(fda.usc)
library(ICglm) ;library(plot3D); library(ggplot2); library(mcmcse)

### need to load
#########################################################################################################
#----------------------------------------------------------------------------
# Exponential Covariance Function
expCov<-function(distMat,phi){
  exp(-distMat/phi)
}

sqeCov<-function(distMat,phi){
  exp(-0.5*(distMat/phi)^2)
}

matCov<-function(distMat,phi){
  (1+(sqrt(5)*(distMat/phi))+((5*distMat^2)/(3*(phi^2))))*exp(-(sqrt(5)*(distMat/phi)))
}

# NIMBLE FUNCTIONS

expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-dists[i,j]/phi)
      }
    }
    
    return(result)
  })

matcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- (1+(sqrt(5)*(dists[i,j]/phi))+((5*dists[i,j]^2)/(3*(phi^2))))*exp(-(sqrt(5)*(dists[i,j]/phi)))
      }
    }
    
    return(result)
  })


## make orthonormalized B-spline basis function ###########
getSplineInfo_d = function(tau, m_eff, orthonormalize = TRUE){
  
  # Just in case, reform as matrix
  tau = as.matrix(tau)
  
  # Number of observation points
  m = nrow(tau)
  
  # Dimension:
  d = ncol(tau)
  
  # Order of derivative in penalty:
  m_deriv = 2
  
  # This is the linear component
  X = cbind(1, tau)
  
  # Number of effective observation points:
  #if(is.null(m_eff)) m_eff = m
  
  # Number of knots: if more than 25 effective observation points, likely can use fewer knots
  #if(m_eff > 25){
  # Guaranteed to be between 20 and 150 knots (but adjust as needed)
  # num_knots = max(20, min(ceiling(m_eff/4), 150))
  #} else num_knots = max(3, m_eff)
  
  num_knots=m_eff
  
  # Locations of knots:
  if(num_knots < m){
    # Usual case: fewer knots than TOTAL observation points
    if(d == 1){
      # d = 1, just use quantiles of the observed data points:
      knots = as.matrix(quantile(unique(tau), seq(0,1,length=(num_knots+2))[-c(1,(num_knots+2))]))
    } else {
      # d > 1, use space-filling algorithm:
      knots = cover.design(tau, num_knots)$design
    }
  } else knots = tau
  
  # For the penalty matrix, need to compute distances between obs. points and knots
  dist.mat <- matrix(0, num_knots, num_knots); dist.mat[lower.tri(dist.mat)] <- dist(knots); dist.mat <- dist.mat + t(dist.mat)
  if(d%%2 == 0){
    # Even dim:
    Omega = dist.mat^(2*m_deriv - d)*log(dist.mat)
  } else {
    # Odd dim:
    Omega = dist.mat^(2*m_deriv - d)
  }
  # For numerical stability:
  diag(Omega) = 0
  
  # Compute the "random effects" matrix
  Zk = matrix(0, nrow=m, ncol=num_knots)
  for (k in 1:num_knots){
    di = sqrt(rowSums((tau - matrix(rep(knots[k,], each = m), nrow=m))^2)) # di = 0; for(j in 1:d) di = di + (tau[,j] - knots[k,j])^2; di = sqrt(di)
    if(d%%2 == 0){# Even dim:
      Zk[,k] = di^(2*m_deriv - d)*log(di)
    } else { # Odd dim:
      Zk[,k] = di^(2*m_deriv - d)
    }
  }
  Zk[is.nan(Zk)] = 0
  
  # Natural constraints, if necessary:
  if(num_knots > m - 1){Q2 = qr.Q(qr(X), complete=TRUE)[,-(1:2)]; Zk = Zk%*%Q2; Omega = crossprod(Q2, Omega)%*%Q2}
  
  # SVD of penalty matrix
  # So that the "random effects" have diagonal prior variance
  svd.Omega = svd(Omega)
  sqrt.Omega = t(svd.Omega$v %*%(t(svd.Omega$u)*sqrt(svd.Omega$d)))
  Z = t(solve(sqrt.Omega,t(Zk)))
  
  # Now combine the linear and nonlinear pieces to obtain the matrix of basis functions evaluated at the obs. points
  Bmat = cbind(X, Z);
  
  # The penalty matrix:
  Omega = diag(c(rep(0, ncol(X)), rep(1, ncol(Z))))
  
  if(orthonormalize){
    # QR decomposition:
    qrd = qr(Bmat, complete = TRUE);  R.t = t(qr.R(qrd));
    # Update hte basis and the penalty matrix:
    Bmat = qr.Q(qrd); Omega = forwardsolve(R.t, t(forwardsolve(R.t, Omega, upper.tri = FALSE)), upper.tri = FALSE)
    
    BtB = diag(1, ncol(Bmat))
  } else BtB = crossprod(Bmat)
  
  # Return the matrix, the penalty, and the cross product (of the basis)
  list(Bmat = Bmat, Omega = Omega, BtB = BtB)
}

#----------------------------------------------------------------------------------------------------------

## get best number of basis function with GCV
get_best_fdobj <- function(x, y, basistype = "fourier", nbasis, lambda) {
  result <- expand.grid(nb = nbasis, l = lambda)
  result$gcv <- NA
  
  for (i in 1:nrow(result)) {
    nb <- result$nb[i]
    l <- result$l[i]
    
    if (basistype == "fourier") {
      data_basis <- create.fourier.basis(c(min(x), max(x)), nbasis=nb)
    } else if (basistype == "bspline") {
      data_basis <- create.bspline.basis(c(min(x), max(x)), nbasis=nb)
    } else if (basistype == "pc") {
      data_basis <- create.pc.basis(y, 1:nb, class.out="fdata")$data
    }
    
    data_fdPar <- fdPar(data_basis, 2, l)
    data_fdobj <- smooth.basis(x, y, data_fdPar)
    result$gcv[i] <- mean(data_fdobj$gcv)
  }
  
  nb <- result$nb[which.min(result$gcv)]
  l <- result$l[which.min(result$gcv)]
  
  if (basistype == "fourier") {
    data_basis <- create.fourier.basis(c(min(x), max(x)),nbasis=nb)
  } else if (basistype == "bspline") {
    data_basis <- create.bspline.basis(c(min(x), max(x)), nbasis=nb)
  } else if (basistype == "pc") {
    data_basis <- create.pc.basis(y, 1:nb, class.out="fdata")$data
  }
  
  data_fdPar <- fdPar(data_basis, 2, l)
  data_fdobj <- smooth.basis(x, y, data_fdPar)
  
  return(list(summary = result, fdobj = data_fdobj))
}

## standardized function
std.scale = function(x){
  return(abs(x-mean(x))/sd(x))
}


##########################################################################################################
## load data

load("ad_mat_korea_246site.RData") # adjacency matrix
load("age_korea_mob_246site.RData") # covariate : age curve
load("weekly_mob_inout_korea_246site.RData") # response : mobility
load("coord_korea_246site.RData") # coordinate

## response - (time grid point) x (location) 
## covarite - (location) x (time grid point)

set.seed(613)
m_data<-weekly_mob_inout_korea
N<-dim(m_data)[2] #location

m=dim(m_data)[1]
pts=c(1:m)
mm=dim(age_korea)[2]
pts2=c(1:mm)

dat<-m_data
dat<-as.matrix(dat)
ind<-sample(c(1:N), round(0.7*N,0))
N<-round(0.7*N,0) 

dat<-(dat-min(dat))/(max(dat)-min(dat)) # min-max scaling for response
age_korea<-(age_korea-min(age_korea))/(max(age_korea)-min(age_korea)) # min-max scaling for covariate

#################################################################################################################
################## get best number of basis functions with GCV for different basis type #########################
#################################################################################################################
##1. B-spline
### for functional response
mobout.best <- get_best_fdobj(1:182, dat, basistype="bspline", nbasis = c(21:40), lambda = c(0))
mobout.best$summary[which.min(mobout.best$summary$gcv),]
### for functional covariate
age.best <- get_best_fdobj(1:16, t(age_korea),basistype="bspline", nbasis = c(4:10), lambda = c(0))
age.best$summary[which.min(age.best$summary$gcv),]

K<- 19 ## basis number for functional response (best summary result - 2)
KK<-10 ## basis number for functional covariate (best summary result - 2)

# construct orthonormal basis function
spinfo<-getSplineInfo_d(pts, K, orthonormalize = T) 
fmat<-spinfo$Bmat
K<-dim(fmat)[2]

spinfo2<-getSplineInfo_d(pts2, KK, orthonormalize = T) 
ffmat<-spinfo2$Bmat
KK<-dim(ffmat)[2]


## 2. Fourier basis
### Fourier basis - nbasis : only odd number
### functional response
#mobout.best <- get_best_fdobj(1:182, dat, basistype="fourier", nbasis = c(27,29,31,33,35,37,39,41,43,45,47,49,51,53), lambda = c(0))
#mobout.best$summary[which.min(mobout.best$summary$gcv),]
### function covariate
#age.best <- get_best_fdobj(1:16, t(age_korea),basistype="fourier", nbasis = c(7,9,11), lambda = c(0))
#age.best$summary[which.min(age.best$summary$gcv),]

#K<- 28 ## basis number for functional response (best summary result)
#KK<-30 ## basis number for functional response (best summary result)

#tdat<-t(dat)
#dat_fd<-fdata(tdat)
#age_korea_fd<-fdata(age_korea)
#fmat<-create.fdata.basis(dat_fd, 1:K, type.basis="fourier", class.out = "fdata")$data
#ffmat<-create.fdata.basis(age_korea_fd, 1:KK, type.basis="fourier", class.out = "fdata")$data

#fmat<-t(as.matrix(fmat))
#ffmat<-t(as.matrix(ffmat))


## 3. FPC basis

#tdat<-t(dat)
#dat_fd<-fdata(tdat)
#age_korea_fd<-fdata(age_korea)

#summary(create.pc.basis(dat_fd, 1:2))
#summary(create.pc.basis(age_korea_fd, 1:3))

#K<- 28 ## basis number for functional response 
#KK<-30 ## basis number for functional response 

#fmat<-create.pc.basis(dat_fd, 1:K)$basis$data
#ffmat<-create.pc.basis(age_korea_fd, 1:KK)$basis$data

#################################################################################################################
#################################################################################################################
#################################################################################################################
Y_train<-dat[,ind]
Y_test<-dat[,-ind]

#### set iteration, burn-in, and thining ####
niter = 200000
nburnin = 100000
thin = 100

Xft<-as.matrix(age_korea)[ind,]
Xft_test<-as.matrix(age_korea)[-ind,]
Xft<-as.matrix(Xft)
Xft_test<-as.matrix(Xft_test)

## ICAR
Qmat<-(diag(rowSums(ad_mat_korea, length(ad_mat_korea))) - ad_mat_korea) 
Qmat_train<-Qmat[ind,ind]
Qmat_test<-Qmat[-ind,-ind]
xft<-Xft%*%ffmat

model_SFoFR <- nimbleCode({
  
  #likelihood part, Data Model
  for (j in 1:k){
    for (i in 1:n){
      XB[i,j]<-XP[i,j] + zero_beta[j]
      lambda[i,j]<-XB[i,j]+W[i,j]
      Y[i,j] ~dnorm(lambda[i,j],tau2) 
    }
  }
  precMat[1:n,1:n]<-nu*covMat[1:n,1:n]
  for (j in 1:k){
    W[1:n,j] ~ dmnorm(mean=ma[1:n], cov=precMat[1:n,1:n])
  }
  XP[1:n,1:k]<-Xf[1:n,1:t]%*%psi[1:t,1:k]
  #parameter model
  for (j in 1:k){
    psi[1:t,j] ~ dmnorm(mean=mn[1:t], cov=covmat[1:t,1:t])
  }
  zero_beta[1:k] ~ dmnorm(mean=mn2[1:k], cov=covmat2[1:k,1:k])
  nu~dgamma(0.5,1/2000)
  tau2~dinvgamma(2,0.1)
})

#constant and input partshttp://127.0.0.1:20307/graphics/f456fec6-f2d8-4a2a-8177-06f4f0bb2430.png
consts<- list(n=dim(Xft)[1], k=dim(fmat)[2], t=dim(ffmat)[2],mn=rep(0,dim(ffmat)[2]),mn2=rep(0,dim(fmat)[2]), 
              covmat=diag(sqrt(100),dim(ffmat)[2]),covmat2=diag(sqrt(100),dim(fmat)[2]),Xf=xft,ma=rep(0,N), covMat=Qmat_train)

inits_ <- list(psi=matrix(rep(0.1,dim(fmat)[2]*dim(ffmat)[2]), nrow=dim(ffmat)[2]), nu=0.001, tau2=1, zero_beta=rep(0,dim(fmat)[2]), 
               W=matrix(rnorm(N*K),nrow=N,ncol=K))
dat_ <- list(Y=t(Y_train)%*%fmat)

pt <- proc.time()
samples_SFoFR <- nimbleMCMC(model_SFoFR, data=dat_, inits=inits_,
                                   constants = consts,
                                   monitors = c("psi","nu","tau2","zero_beta", "lambda"),
                                   samplesAsCodaMCMC = TRUE, WAIC=FALSE, summary=FALSE,
                                   niter = niter, nburnin = nburnin, thin = thin, nchains = 1)


ptFinal_glm <- proc.time()-pt
ptFinal_glm




psi.sample<-samples_SFoFR[,which(substr(colnames(samples_SFoFR),1,3)=="psi")]
lambda.sample<-samples_SFoFR[,which(substr(colnames(samples_SFoFR),1,6)=="lambda")]
beta0.sample<-samples_SFoFR[,which(substr(colnames(samples_SFoFR),1,4)=="zero")]
nu.sample<-samples_SFoFR[,which(substr(colnames(samples_SFoFR),1,2)=="nu")]
tau2.sample<-samples_SFoFR[,which(colnames(samples_SFoFR)=="tau2")]


B=(niter-nburnin)/thin #number of iteration - number of burn-in
BB=B #number of MCMC samples using in kriging(backward). If BB=B, we use all samples.


mean_psi<-colMeans(psi.sample)
mean_psi<-t(as.matrix(mean_psi))

psi_ma<-NULL

for(i in 1:K){
  psi_mat<-mean_psi[,(KK*(i-1)+1):(KK*i)]
  psi_ma<-rbind(psi_ma,psi_mat)
}

psi_ori<-ffmat%*%t(psi_ma)%*%t(fmat)

### 3D plot of regression function
persp3D(pts2,pts,psi_ori)

psii_ma<-NULL
for (i in 1:BB){
  for (j in 1:K){
    psii_mat<-psi.sample[B-(BB-i),(KK*(j-1)+1):(KK*j)]
    psii_ma<-rbind(psii_ma,psii_mat)
  }
}

psii_maa<-NULL
for(j in 1:BB){
  psii_mata<-ffmat%*%t(psii_ma[(K*(j-1)+1):(K*j),])%*%t(fmat)
  psii_maa<-rbind(psii_maa,psii_mata)
}

psi_orisam<-matrix(rep(0,BB*mm*m), nrow=BB)

for(j in 1:m){
  for(k in 1:mm){
    for(i in 1:BB){
      psi_orisam[i,k+mm*(j-1)] <- psii_maa[k+mm*(i-1),j]
    }
  }
}


arr<-t(psi_orisam)
arr_M_a<-t(apply(arr,1,std.scale))
M_a<-t(apply(arr_M_a,2,max))
M_a_q=quantile(M_a, probs=1-0.05)


low=as.matrix(rowMeans(arr)-M_a_q*apply(arr,1,sd))
up=as.matrix(rowMeans(arr)+M_a_q*apply(arr,1,sd))

sig<-rep(0.5,(mm*m))

alpha<-c(0.5,0.4,0.3,0.2,0.1,0.075,0.05,0.025,0.01,0.005,0.001)
for(i in 1:length(alpha)){
  M_a_quan=quantile(M_a, probs=1-alpha[i])
  lows=as.matrix(rowMeans(arr)-M_a_quan*apply(arr,1,sd))
  ups=as.matrix(rowMeans(arr)+M_a_quan*apply(arr,1,sd))
  
  for(j in 1:(mm*m)){
    if(lows[j]>0){
      sig[j]<-alpha[i]
    }
    
    if(ups[j]<0){
      sig[j]<-alpha[i]
    }
  }  
}
sig<-as.matrix(sig)

sig_ma<-NULL
for(j in 1:m){
  sig_mat<-sig[(mm*(j-1)+1):(mm*j),]
  sig_ma<-rbind(sig_ma,sig_mat)
}

### draw significant band score ###
x <- pts2
y <- pts
data <- expand.grid(X=x, Y=y)
data$Z <- sig
graph<-ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()

print(graph)

### draw 95% significant credible interval ###
data$Z<-(sig<=0.05)*sign((rowMeans(arr)))
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()


arr<-fmat%*%t(beta0.sample)
arr_M_a<-t(apply(arr,1,std.scale))
M_a<-t(apply(arr_M_a,2,max))
M_a_q=quantile(M_a, probs=1-0.05)
rm(arr_M_a)
rm(M_a)

low=rowMeans(arr)-M_a_q*apply(arr,1,sd)
up=rowMeans(arr)+M_a_q*apply(arr,1,sd)


## plot estimated beta function with credible interval
par(mfrow=c(1,1))
pos_mean<-colMeans(beta0.sample)
pos_mean<-as.matrix(pos_mean)
meanpos<-fmat%*%pos_mean
#matplot(pts,arr,type="l")


matplot(pts,meanpos, type="l", lty=3, lwd=2,col="red", ylim=c(min(low), max(up))) # posterior mean
lines(pts,low, lty=1, col="red")  # lower bound of credible interval
lines(pts,up,lty=1, col="red")    # upper bound of credible interval


### draw trace plots of posterior samples ###
plot(as.vector(tau2.sample), type="l", main="tau2")
plot(as.vector(nu.sample), type="l", main="nu")

par(mfrow=c(2,2))
plot(as.vector(psi.sample[,1]), type="l", main="psi1_1")
plot(as.vector(psi.sample[,2]), type="l", main="psi2_1")
plot(as.vector(psi.sample[,3]), type="l", main="psi3_1")
plot(as.vector(psi.sample[,4]), type="l", main="psi4_1")


plot(as.vector(beta0.sample[,1]), type="l", main="b0_1")
plot(as.vector(beta0.sample[,2]), type="l", main="b0_2")
plot(as.vector(beta0.sample[,3]), type="l", main="b0_3")
plot(as.vector(beta0.sample[,4]), type="l", main="b0_4")

### draw all trace plots of MCMC samples ###
#par(mfrow=c(4,4))
#for(num in 1:dim(beta0.sample)[2]){
#  plot(as.vector(beta0.sample[,num]), type="l", main=paste("b0_", num))
#}
#par(mfrow=c(4,4))
#for(num in 1:dim(psi.sample)[2]){
#  plot(as.vector(psi.sample[,num]), type="l", main=paste("psi_", num))
#}

para_mcmcse<-mcse.mat(x = samples_SFoFR, method = "bm", g = NULL)
para_mcmcse<-data.frame(para_mcmcse)
mean(para_mcmcse$se)

save.image("mobility_car_inout_all_gamma_200000.RData")
save(Y_test, psi_ori,psi_orisam, file="final_mobility_car_inout_gamma_200000.RData")



