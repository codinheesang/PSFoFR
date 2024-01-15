#####################################################################################################################################

# 3D plot, Regression function Estimate, Significant area Plot For Mobility Data

#####################################################################################################################################

#Data load
load("final_mobility_inout_rank10per_kriging_mean_gamma.RData")

m=182
mm=16
pts=c(1:m)
pts2=c(1:mm)


std.scale = function(x){
  return(abs(x-mean(x))/sd(x))
}

arr<-t(psi_orisam)
arr_M_a<-t(apply(arr,1,std.scale))
M_a<-t(apply(arr_M_a,2,max))
M_a_q=quantile(M_a, probs=1-0.05)
# rm(arr_M_a)
# rm(M_a)


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
x <- (pts2)
y <- (pts)
data <- expand.grid(X=x, Y=y)
data$Z <- sig
graph<-ggplot(data, aes(X, Y, fill= Z)) + 
  scale_y_continuous(breaks=c(0,50,100,150),labels=c('19/04', '20/03', '21/02', '22/02'))+
  scale_x_continuous(breaks=c(0,5,10,15),labels=c('0','30-34', '55-59', '80-84'))+
  geom_tile()

print(graph)

data$Z<-(sig<=0.05)*sign((rowMeans(arr)))



library(plot3D)
persp3D(pts2,pts,psi_ori, xlab="t", ylab="r", zlab="value")


library(cowplot)

data$Z<-as.vector(psi_ori)
esti<-ggplot(data, aes(X, Y, fill= Z))+
  labs(x="t", y="r")+
  geom_tile() +
  theme(
    legend.title=element_blank(),
    axis.text=element_text(size=13),
    axis.title=element_text(size=14,face="bold"),
    legend.text = element_text(size=12)
  )+
  scale_y_continuous(breaks=c(0,50,100,150),labels=c('19/04', '20/03', '21/02', '22/02'))+
  scale_x_continuous(breaks=c(0,5,10,15),labels=c('0','30-34', '55-59', '80-84'))+
  geom_tile()

data$Z<-(sig<=0.05)*sign((rowMeans(arr)))
data$Z<-as.factor(data$Z)
band<-ggplot(data, aes(X, Y, fill= Z)) + 
  scale_fill_manual(values=c("#56B1F7", "#FFFFFF", "#132B43"), 
                    breaks=c("1","0","-1"),
                    labels=c("+", "no", "-"))+
  geom_tile() +
  guides(fill=guide_legend(title=NULL))+
  labs(x="t", y="r")+
  theme(
    axis.text=element_text(size=13),
    axis.title=element_text(size=14,face="bold"),
    legend.text = element_text(size=12)
  )+scale_y_continuous(breaks=c(0,50,100,150),labels=c('19/04', '20/03', '21/02', '22/02'))+
  scale_x_continuous(breaks=c(0,5,10,15),labels=c('0','30-34', '55-59', '80-84'))

plot_grid(esti, band, ncol=2)

