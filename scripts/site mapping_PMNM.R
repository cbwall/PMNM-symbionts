rm(list=ls())
ls()

library(RgoogleMaps)
library(SDMTools)
library(devtools)

data<-read.csv("data/PMNM_latlong.csv", header=TRUE, na.string=NA)
data<-subset(data, !(is.na(Site.ID)))

par(mfrow=c(2,1))

#RGoogleMaps
# all samples
PMNM <- c(24.407320, -168.089023)
PMNM.map <- GetMap(center = PMNM, zoom = 5, maptype = "satellite", SCALE = 2)
PlotOnStaticMap(PMNM.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km

#zoomed in to just NWHI
PMNM2 <- c(24.907320, -171.9023)
PMNM2.map <- GetMap(center = PMNM2, zoom = 6, maptype = "satellite", SCALE = 2)
PMNM.image<-PlotOnStaticMap(PMNM2.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km

dev.copy(pdf, "PMNM.image.pdf")
dev.off()

par(mfrow=c(2,2))

# FFS
FFS <- c(23.75,	-166.206)
FFS.map <- GetMap(center = FFS, zoom = 11, maptype = "satellite", SCALE = 2)
FFS.image<-PlotOnStaticMap(FFS.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km

# Maro
MAR <- c(25.40866667,	-170.5874)
MAR.map <- GetMap(center = MAR, zoom = 12, maptype = "satellite", SCALE = 2)
MAR.image<-PlotOnStaticMap(MAR.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km

# Laysan
LAY <- c(25.7695,	-171.7382)
LAY.map <- GetMap(center = LAY, zoom = 14, maptype = "satellite", SCALE = 2)
LAY.image<-PlotOnStaticMap(LAY.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km


# Lisianski
LIS <- c(26.0038,	-173.9506)
LIS.map <- GetMap(center = LIS, zoom = 12, maptype = "satellite", SCALE = 2)
LIS.image<-PlotOnStaticMap(LIS.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km

dev.copy(pdf, "PMNM_South_Central.pdf")
dev.off()

par(mfrow=c(2,2))
# Pearl and Hermes
PHA <- c(27.81994,	-175.89474)
PHA.map <- GetMap(center = PHA, zoom = 11, maptype = "satellite", SCALE = 2)
PHA.image<-PlotOnStaticMap(PHA.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km


# Midway
MID <- c(28.22913,	-177.38635)
MID.map <- GetMap(center = MID, zoom = 12, maptype = "satellite", SCALE = 2)
MID.image<-PlotOnStaticMap(MID.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km


# Kure
KUR <- c(28.42348333,	-178.3286)
KUR.map <- GetMap(center = KUR, zoom = 13, maptype = "satellite", SCALE = 2)
KUR.image<-PlotOnStaticMap(KUR.map, data$latitude, data$longitude, col="orange", bg="red", pch=21)
Scalebar(x=-300,y=-275,distance=200) #show values in km
save.image("KUR.image")


dev.copy(pdf, "PMNM_North.pdf")
dev.off()
