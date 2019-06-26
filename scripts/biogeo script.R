########
rm(list=ls())
ls()

library(effects)
library(gplots)
library(plotrix)
library(psych)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(MASS)
library(dplyr)

data<-read.csv("data/POM/chla_PMNM.csv")

names(data)
data<-data[, c(1,3:5,7)]
colnames(data)<-c("Island", "Site", "depth.ft", "Loc", "chla")
data$Island<-factor(data$Island, levels(data$Island)[c(2,6,7,4,3,5,1)])

ggplot(data, aes(x= Island, y=chla, colour=factor(Island))) +
  ylab(expression(paste("chlorophyll", ~italic("a"), ~(mu*g~ L^-1), sep=""))) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  facet_grid(. ~Island, scales ="free") +
  theme(strip.text.x = element_text(size = 5)) +
  geom_point(size=2)
dev.copy(pdf, "figures/POM/chla.site.pdf", height=5, width=7)
dev.off()


#### Site level means
chla.site.mean<-aggregate(chla~Island+Loc, data, FUN=mean)
chla.site.SE<-aggregate(chla~Island+Loc, data, FUN=std.error)
chla.site.count<-aggregate(chla~Island+Loc, data, FUN=length)

#write.csv(chla.site.mean, "data/chla.mean.csv")

ggplot(chla.site.mean, aes(x= Island, y=chla, colour=factor(Island))) +
  ylab(expression(paste("chlorophyll", ~italic("a"), ~(mu*g~ L^-1), sep=""))) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  facet_grid(. ~Island, scales ="free") +
  theme(strip.text.x = element_text(size = 5)) +
  geom_point(size=2)
dev.copy(pdf, "figures/POM/chla.sitemean.pdf", height=5, width=7)
dev.off()


############
iso.data<-read.csv("data/POM/SPM.PMNM.isotopes.csv")
names(iso.data)
colnames(iso.data)<-c("Sample.ID", "Island", "Site", "depth.ft", "filter.number", "acidified", "N.ug", "d15N", "C.ug", "d13C")
iso.data$Island<-factor(iso.data$Island, levels(iso.data$Island)[c(2,6,7,4,3,5,1)])

acid.iso<-iso.data[(iso.data$acidified=="Y"),] ## this data is for the d13C, acidified 
#write.csv(acid.iso, "d13C.csv")

noacid.iso<-iso.data[(iso.data$acidified=="N"),] ## this data is for d15N, non-acidified
#write.csv(noacid.iso, "d15N.csv")

### non acidified samples d15N
ggplot(noacid.iso, aes(x= Island, y=d15N, colour=factor(Island))) +
  ylab("d15N permil") +
  ylim(2, 8) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  facet_grid(. ~Island, scales ="free") +
  theme(strip.text.x = element_text(size = 5)) +
  geom_point(size=2)
dev.copy(pdf, "figures/POM/d15N.SPM.site.pdf", height=5, width=7)
dev.off()


## acidifed d13C 
ggplot(acid.iso, aes(x= Island, y=d13C, colour=factor(Island))) +
  ylab("d13C permil") +
  ylim(-27, -17) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  facet_grid(. ~Island, scales ="free") +
  theme(strip.text.x = element_text(size = 5)) +
  geom_point(size=2)
dev.copy(pdf, "figures/POM/d13C.SPM.site.pdf", height=5, width=7)
dev.off()


## looking at chla and d15N, d13C
##########
dfCNchl<-read.csv("data/POM/iso.chla.csv")
dfCNchl$Island<-factor(dfCNchl$Island, levels(dfCNchl$Island)[c(2,6,7,4,3,5,1)])
names(dfCNchl)

dfCNchl$depth..m<- (dfCNchl$depth..ft*0.3048)

summary(lm(d13C~chla, data=dfCNchl)) # no effect, but a negative relationship, more chla and more depleted d13C from photosynthesis removing light carbon, remineralization of detritus at sediment surface, sources?

par(mfrow=c(3,1), mar=c(3,4,1,1), mgp=c(2,0.5,0))
plot(d13C~chla, data=dfCNchl, pch=16, col="coral",
     xlab=expression(paste("chl", ~italic("a"), ~(mu*g~L^-1), sep="")),
     ylab=expression(paste(delta^{13}, C[POM], " (\u2030)")))
abline(lm(d13C~chla, data=dfCNchl), col="red")

plot(d15N~chla, data=dfCNchl, pch=16, col="dodgerblue",
     xlab=expression(paste("chl", ~italic("a"), ~(mu*g~L^-1), sep="")),
     ylab=expression(paste(delta^{15}, N[POM], " (\u2030)")))
abline(lm(d15N~chla, data=dfCNchl), col="dodgerblue3")
summary(lm(d15N~chla, data=dfCNchl))

plot(d15N~d13C, data=dfCNchl, col="darkseagreen", pch=16,
     ylab=expression(paste(delta^{15}, N[POM], " (\u2030)")),
     xlab=expression(paste(delta^{13}, C[POM], " (\u2030)")))
abline(lm(d15N~d13C, data=dfCNchl), col="forestgreen")
summary(lm(d15N~d13C, data=dfCNchl))

dev.copy(pdf, "figures/POM/chla plus ios.pdf", width=4, height=6,  encod="MacRoman")
dev.off()


## looking at DEPTH and d15N, d13C
##########

par(mfrow=c(3,1), mar=c(3,4,1,1), mgp=c(2,0.5,0))
plot(d13C~depth..m, data=dfCNchl, pch=16, col="coral",
     xlab="depth (m)",
     ylab=expression(paste(delta^{13}, C[POM], " (\u2030)")))
abline(lm(d13C~depth..m, data=dfCNchl), col="red")
summary(lm(d13C~depth..m, data=dfCNchl)) # depth is significant here

plot(d15N~depth..m, data=dfCNchl, pch=16, col="dodgerblue", 
     xlab="depth (m)",
     ylab=expression(paste(delta^{15}, N[POM], " (\u2030)")))
abline(lm(d15N~depth..m, data=dfCNchl), col="dodgerblue3")
summary(lm(d15N~depth..m, data=dfCNchl))


## scatter by site
p1<-ggplot(dfCNchl, aes(x=depth..m, y=d13C, colour=factor(Island))) +
  theme(legend.position = "none", axis.text.x = element_blank()) + 
  theme(strip.text.x = element_text(size = 5)) +
  ylab(expression(paste(delta^{13}, C[POM], " (\u2030)"))) +
  xlab("") +
  geom_point(size=2) + facet_grid(. ~Island)

p2<-ggplot(dfCNchl, aes(x=depth..m, y=d15N, colour=factor(Island))) +
  theme(legend.position = "none") + theme(strip.text.x = element_text(size = 5)) +
  ylab(expression(paste(delta^{15}, N[POM], " (\u2030)"))) +
  xlab("depth (m)")+
  geom_point(size=2) + facet_grid(. ~Island)

grid.arrange(p1, p2, ncol = 1)

dev.copy(pdf, "figures/POM/depth.iso.pdf", width=4, height=6,  encod="MacRoman")
dev.off()
