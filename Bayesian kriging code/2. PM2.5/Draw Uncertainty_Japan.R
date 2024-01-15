#############################################
# Draw Uncertainty map for Japan data  #
#############################################

library(sp);library(fda);library(refund)
library(spatial);library(mvtnorm);library(fields);library(rdist)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(classInt)
library(gstat) ;library(nlme) ; library(fda.usc)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(ggmap)
library(raster)
library(ggmap)
library(ggplot2)

load("final_pm25_final_rank5per_kriging_mean_gammaprior_lowup_align_bspline.RData")
load("JapanFdagstat_Align_final.RData")
load("japan_coord_dat_align.RData")
load("japan_align_train_index.RData")
load("japan_pm25_dat_align.RData")

map<-shapefile("JapanMap/JPN_adm1.shp")
map <- spTransform(map, CRSobj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
new_map <- fortify(map)
df_for<-new_map
#ne_map = ne_countries(scale = "small", country = "Japan")
#df_for = fortify(ne_map)

ggplot(data = df_for) + 
  geom_polygon(aes(x = long, y = lat,
                   group = group),
               fill = "#FFFFFF",
               color = "#000000")


Latitude<-as.matrix(m_coord[-ind,1])
Longitude<-as.matrix(m_coord[-ind,2])
Uncertainty<-as.matrix(colMeans(up-low))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(Uncertainty), max(Uncertainty)))

df<-cbind(Longitude, Latitude, Uncertainty)
df<-data.frame(df)
ggplot(df, aes(x=Longitude, y=Latitude, color=Uncertainty)) +
  geom_polygon(data = df_for, aes(x = long, y = lat,
                                  group = group),
               fill = "#FFFFFF",
               color = "#000000")+
  geom_point(size=2)+
  #scale_fill_gradient(low = "#ffe5e5", high = "#ff3232", space = "Lab", guide = "colourbar")+
  coord_cartesian(xlim = c(127,148), ylim = c(30,46)) +
  sc+
  theme_bw()+
  theme(axis.line=element_line(color="black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(
    legend.title=element_blank(),
    axis.text=element_text(size=13),
    axis.title=element_text(size=14,face="bold"),
    legend.text = element_text(size=12)
  )

