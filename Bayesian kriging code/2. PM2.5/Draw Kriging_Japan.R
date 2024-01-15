################################
# Draw Kriging Map for Japan   #
################################

load("JapanFdagstat_Align_final.RData")
load("final_japan_matern_nu0.5_with_lowup_align.RData")
load("final_japan_matern_nu0.5_align.RData")  
load("final_pm25_final_rank5per_kriging_mean_gammaprior_align_bspline.RData")
load("final_pm25_final_rank5per_kriging_mean_gammaprior_lowup_align_bspline.RData")
library(ggplot2)
m=52
pt<-c(1:m)

lambdaa <- lambdaa/2
PSFoFR<-etaa[,11]
SFoFR<-lambdaa[,11]
fdagstat<-k2[,11]
true<-Y_test[,11]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))

p1 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_x_continuous(labels=c('22/07', '22/09', '22/11', '23/01', "23/04", "23/06"))+
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")



PSFoFR<-etaa[,29]
SFoFR<-lambdaa[,29]
fdagstat<-k2[,29]
true<-Y_test[,29]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p2 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_x_continuous(labels=c('22/07', '22/09', '22/11', '23/01', "23/04", "23/06"))+
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")

PSFoFR<-etaa[,41]
SFoFR<-lambdaa[,41]
fdagstat<-k2[,41]
true<-Y_test[,41]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p3 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_x_continuous(labels=c('22/07', '22/09', '22/11', '23/01', "23/04", "23/06"))+
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")


PSFoFR<-etaa[,51]
SFoFR<-lambdaa[,51]
fdagstat<-k2[,51]
true<-Y_test[,51]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p4 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_x_continuous(labels=c('22/07', '22/09', '22/11', '23/01', "23/04", "23/06"))+
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")

library(cowplot)
plot_grid(p1 , p2, p3, p4,  ncol=2)

