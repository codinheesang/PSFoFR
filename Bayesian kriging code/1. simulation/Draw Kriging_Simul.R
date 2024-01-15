########################################################################################################################
##kriging plot : WITH TAU
##simulation
##bspline
load("final_simu_matern_tau1_nu0.5_kriging_mean_gamma_lowup.RData")
load("final_simu_final_tau1_sigma_0.5_nu0.5_rank5per_gammaprior_lowup.RData")
load("Fdagstat_WithErrorSimul.RData")
########################################################################################################################

library(ggplot2)
m=225
pt<-c(1:m)

PSFoFR<-etaa[,1]
SFoFR<-lambdaa[,1]
fdagstat<-k2[,1]
true<-Y_test[,1]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))

p1 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")



PSFoFR<-etaa[,2]
SFoFR<-lambdaa[,2]
fdagstat<-k2[,2]
true<-Y_test[,2]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p2 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")

PSFoFR<-etaa[,3]
SFoFR<-lambdaa[,3]
fdagstat<-k2[,3]
true<-Y_test[,3]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p3 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")


PSFoFR<-etaa[,4]
SFoFR<-lambdaa[,4]
fdagstat<-k2[,4]
true<-Y_test[,4]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p4 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
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
#legend <- get_legend(p1)
plot_grid(p1 , p2, p3, p4,  ncol=2)

########################################################################################################################
##kriging plot : WITHOUT TAU
##simulation
##bspline
load("final_simu_matern_tauno_nu0.5_kriging_mean_gamma_lowup.RData")
load("final_simu_final_tauno_sigma_0.5_nu0.5_rank5per_gammaprior_lowup.RData")
load("Fdagstat_WithoutErrorSimul.RData")
########################################################################################################################

library(ggplot2)
m=225
pt<-c(1:m)

PSFoFR<-etaa[,1]
SFoFR<-lambdaa[,1]
fdagstat<-k2[,1]
true<-Y_test[,1]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))

p1 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")



PSFoFR<-etaa[,2]
SFoFR<-lambdaa[,2]
fdagstat<-k2[,2]
true<-Y_test[,2]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p2 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")

PSFoFR<-etaa[,3]
SFoFR<-lambdaa[,3]
fdagstat<-k2[,3]
true<-Y_test[,3]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p3 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
  scale_linetype_manual("",values = c(rep("solid",4))) +
  scale_colour_manual("", breaks = c("true", "PSFoFR", "SFoFR", "fdagstat"), values = c("black","red","blue", "green")) +
  xlab('time') + ylab('values')+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")


PSFoFR<-etaa[,4]
SFoFR<-lambdaa[,4]
fdagstat<-k2[,4]
true<-Y_test[,4]


df<-as.data.frame(cbind(PSFoFR,SFoFR,fdagstat,pt,true))


p4 = ggplot() +
  geom_line(data = df, aes(x = pt, y = true, colour="true"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = PSFoFR, colour="PSFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = SFoFR, colour="SFoFR"), linewidth=0.5) +
  geom_line(data = df, aes(x = pt, y = fdagstat, colour="fdagstat"), linewidth=0.5) +
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
#legend <- get_legend(p1)
plot_grid(p1 , p2, p3, p4,  ncol=2)

