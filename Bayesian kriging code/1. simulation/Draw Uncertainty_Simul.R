#############################################
# Draw Uncertainty map for simulation data  #
#############################################

library(sp);library(fda);library(refund)
library(spatial);library(mvtnorm);library(fields);library(rdist)
library(classInt)
library(gstat) ;library(nlme) ; library(fda.usc)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(ggmap)
library(raster)

##############################
# Without Measurement Error  #
##############################

load("final_simu_final_tauno_sigma_0.5_nu0.5_rank5per_gammaprior_lowup.RData")
load("simu_all_coord.RData")
load("simu_train_index.RData")

par(mfrow=c(1,1))

ncoord <- ncoord[-ind,]
Longitude<-as.matrix(ncoord[,1])
Latitude<-as.matrix(ncoord[,2])
Uncertainty<-as.matrix(colMeans(up-low))
#
#
myPalette <- colorRampPalette(rev(brewer.pal(9, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,0.65))
#
#
df<-cbind(Longitude, Latitude, Uncertainty)
df<-data.frame(df)
ggplot(df, aes(x=Longitude, y=Latitude, color=Uncertainty)) +
geom_point(size=3)+
   sc+
   theme(
     legend.title=element_blank(),
     axis.text=element_text(size=13),
     axis.title=element_text(size=14,face="bold"),
     legend.text = element_text(size=12)
   ) + theme_bw() + theme(axis.line = element_line(color='black'),
                          plot.background = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank())


###########################
# With Measurement Error  #
###########################

load("final_simu_final_tau1_sigma_0.5_nu0.5_rank5per_gammaprior_lowup.RData")
load("simu_all_coord.RData")
load("simu_train_index.RData")


par(mfrow=c(1,1))

ncoord <- ncoord[-ind,]

Longitude<-as.matrix(ncoord[,1])
Latitude<-as.matrix(ncoord[,2])
Uncertainty<-as.matrix(colMeans(up-low))
#
#
myPalette <- colorRampPalette(rev(brewer.pal(9, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,0.66))
#
#
df<-cbind(Longitude, Latitude, Uncertainty)
df<-data.frame(df)
ggplot(df, aes(x=Longitude, y=Latitude, color=Uncertainty)) +
  geom_point(size=3)+
  sc+
  theme(
    legend.title=element_blank(),
    axis.text=element_text(size=13),
    axis.title=element_text(size=14,face="bold"),
    legend.text = element_text(size=12)
  ) + theme_bw() + theme(axis.line = element_line(color='black'),
                         plot.background = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank())

