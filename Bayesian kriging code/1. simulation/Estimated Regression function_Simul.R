
#####################################################################################################################################

# 3D plot, Regression function Estimate, Significant area Plot For Simulation Data

#####################################################################################################################################


#####################################################################################################################################
#WITHOUT MEASUREMENT ERROR
#bspline-regression function
#load("final_simu_final_tauno_sigma_0.5_nu0.5_rank5per_gammaprior.RData") #PSFoFR
#load("final_simu_matern_tauno_nu0.5_kriging_mean_gamma.RData") #SFoFR

#fourier-regression function
#load("final_simu_final_tauno_sigma_0.5_nu0.5_rank5per_gammaprior_fourier.RData")
#pc-regression function
#load("final_simu_final_tauno_sigma_0.5_nu0.5_rank5per_gammaprior_pc.RData")
#####################################################################################################################################

#####################################################################################################################################
#WITH MEASUREMENT ERROR
#bspline-regression function
#load("final_simu_final_tau1_sigma_0.5_nu0.5_rank5per_gammaprior.RData") #PSFoFR
#load("final_simu_matern_tau1_nu0.5_kriging_mean_gamma.RData") #SFoFR

#fourier-regression function
#load("final_simu_final_tau1_sigma_0.5_nu0.5_rank5per_gammaprior.RData")
#pc-regression function
#load("final_simu_final_tau1_sigma_0.5_nu0.5_rank5per_gammaprior_pc.RData")
#####################################################################################################################################

m=225
mm=225
pts=c(1:m)
pts2=c(1:mm)

psi_m = matrix(rep(0,m*mm), ncol=m)

for(j in 1:m){
  for(i in 1:mm){
    psi_m[i,j]=(7/500)*(1/sqrt(2*pi*0.003))*exp(-(1/(2*0.003))*(i/225-j/225)^2)
  }
}

std.scale = function(x){
  return(abs(x-mean(x))/sd(x))
}

arr<-t(psi_orisam)
arr_M_a<-t(apply(arr,1,std.scale))
M_a<-t(apply(arr_M_a,2,max))
M_a_q=quantile(M_a, probs=1-0.05)
#rm(arr_M_a)
#rm(M_a)


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


library(ggplot2)
x <- pts2
y <- pts
data <- expand.grid(X=x, Y=y)
data$Z <- sig
graph<-ggplot(data, aes(X, Y, fill= Z)) + geom_tile()
print(graph)

data$Z<-(sig<=0.05)*sign((rowMeans(arr)))

par(mfrow=c(1,2))
library(plot3D)
persp3D(pts2,pts,psi_m, zlim=c(0,0.1), clim=c(-0.02,0.1), xlab="t", ylab="r", zlab="value") #TRUE
persp3D(pts2,pts,psi_ori, zlim=c(0,0.1), clim=c(-0.02,0.1), xlab="t", ylab="r", zlab="value")


library(cowplot)

#check color scale
data$Z<-as.vector(psi_m)
true<-ggplot(data, aes(X, Y, fill= Z)) +scale_fill_continuous(limits = c(-0.02,0.16))+
  geom_tile() +
  labs(x="t", y="r", colour="value")+
  theme(
    legend.title=element_blank(),
    axis.text=element_text(size=13),
    axis.title=element_text(size=14,face="bold"),
    legend.text = element_text(size=12)
  )


data$Z<-as.vector(psi_ori)
print(range(data$Z))
esti<-ggplot(data, aes(X, Y, fill= Z)) +scale_fill_continuous(limits = c(-0.02,0.16))+geom_tile() +
  labs(x="t", y="r", colour="value")+
  theme(
    legend.title=element_blank(),
    axis.text=element_text(size=13),
    axis.title=element_text(size=14,face="bold"),
    legend.text = element_text(size=12)
  )

data$Z<-(sig<=0.05)*sign((rowMeans(arr)))
data$Z<-as.factor(data$Z)
band<-ggplot(data, aes(X, Y, fill= Z)) + 
  scale_fill_manual(values=c("#56B1F7", "#FFFFFF", "#132B43"), 
                    breaks=c("1","0","-1"),
                    labels=c("+", "no", "-"))+
  guides(fill=guide_legend(title=NULL))+
  geom_tile() +
  labs(x="t", y="r")+
  theme(
    axis.text=element_text(size=13),
    axis.title=element_text(size=14,face="bold"),
    legend.text = element_text(size=12)
  )
plot_grid(true, esti, band, ncol=3)
