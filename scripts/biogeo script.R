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

data<-read.csv("data/chla_PMNM.csv")

names(data)
data<-data[, c(1,3,4,6)]
colnames(data)<-c("Island", "Site", "Loc", "chla")
data$Island<-factor(data$Island, levels(data$Island)[c(2,6,7,4,3,5,1)])

ggplot(data, aes(x= Island, y=chla, colour=factor(Island))) +
  ylab(expression(paste("chlorophyll", ~italic("a"), ~(mu*g~ L^-1), sep=""))) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  facet_grid(. ~Island, scales ="free") +
  theme(strip.text.x = element_text(size = 5)) +
  geom_point(size=2)
dev.copy(pdf, "figures/chla.site.pdf", height=5, width=7)
dev.off()


#### Site level means
chla.site.mean<-aggregate(chla~Island+Loc, data, FUN=mean)
chla.site.SE<-aggregate(chla~Island+Loc, data, FUN=std.error)
chla.site.count<-aggregate(chla~Island+Loc, data, FUN=length)

as.data.frame(chla.site.mean)
chla.site.mean$Island<-factor(chla.site.mean$Island, levels(chla.site.mean$Island)[c(2,6,7,4,3,5,1)])

write.csv(chla.site.mean, "data/chla.mean.csv")


############
iso.data<-read.csv("data/SPM.PMNM.isotopes.csv")
names(iso.data)
colnames(iso.data)<-c("Sample.ID", "Island", "Site", "filter.number", "acidified", "N.ug", "d15N", "C.ug", "d13C")
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
dev.copy(pdf, "figures/d15N.SPM.site.pdf", height=5, width=7)
dev.off()


## acidifed d13C 
ggplot(acid.iso, aes(x= Island, y=d13C, colour=factor(Island))) +
  ylab("d13C permil") +
  ylim(-27, -17) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  facet_grid(. ~Island, scales ="free") +
  theme(strip.text.x = element_text(size = 5)) +
  geom_point(size=2)
dev.copy(pdf, "figures/d13C.SPM.site.pdf", height=5, width=7)
dev.off()


## looking at chla and d15N, d13C
##########
dfCNchl<-read.csv("data/iso.chla.csv")
names(dfCNchl)
summary(lm(d13C~chla, data=dfCNchl)) # no effect, but a negative relationship, more chla and more depleted d13C from photosynthesis removing light carbon, remineralization of detritus at sediment surface, sources?

par(mfrow=c(3,1), mar=c(3,4,1,1), mgp=c(2,0.5,0))
plot(d13C~chla, data=dfCNchl, pch=16, col="coral",
     xlab=expression(paste("chl", ~italic("a"), ~(mu*g~L^-1), sep="")),
     ylab=expression(paste(delta^{13}, C[POM], " (\u2030)")))
abline(lm(d13C~chla, data=dfCNchl), col="red")

plot(d15N~chla, data=dfCNchl, pch=16, col="dodgerblue",
     xlab=expression(paste("chl", ~italic("a"), ~(mu*g~L^-1), sep="")),
     ylab=expression(paste(delta^{13}, C[POM], " (\u2030)")))
abline(lm(d15N~chla, data=dfCNchl), col="dodgerblue3")
summary(lm(d15N~chla, data=dfCNchl))

plot(d15N~d13C, data=dfCNchl, col="darkseagreen", pch=16,
     ylab=expression(paste(delta^{15}, N[POM], " (\u2030)")),
     xlab=expression(paste(delta^{13}, C[POM], " (\u2030)")))
abline(lm(d15N~d13C, data=dfCNchl), col="forestgreen")
summary(lm(d15N~d13C, data=dfCNchl))

dev.copy(pdf, "figures/chla plus ios.pdf", width=4, height=6,  encod="MacRoman")
dev.off()



