### Andreas PMNM data
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

# use pacman to load all the packages you are missing!
pacman::p_load("ggplot2", "dplyr", "vegan")

# put data in data folder, in the working directory, you should use R projects 
# check working directory "wd"
getwd()
dat<-read.csv("data/Kealoha collab/Coral_isotopes_andrea.csv")

# check structure
str(dat)

# make these factors
make.fac<-c("Island.atoll", "Species", "Site.name", "Type")
dat[make.fac] <- lapply(dat[make.fac], factor) # make all these factors

# make these numeric
make.numeric<-c("NO3", "SiO2", "PO4", "NH3")
dat[make.numeric] <- lapply(dat[make.numeric], as.numeric) # make all these factors

# make difference colum
dat$d13C.diff<- dat$d13C.host-dat$d13C.alg
dat["d13C.diff"][dat["d13C.diff"] == -14.7] <- NA # remove an outlier
dat["d15N.alg"][dat["d15N.alg"] == 15.8] <- NA # remove an outlier

dat$d15N.diff<- dat$d15N.host-dat$d15N.alg

######### ######### ######### 
######### plots
######### ######### ######### 

### difference of host and symbiont d15N
d13C.diff.species<-ggplot(dat, aes(x=Type, y=d13C.diff, fill =Species)) + 
  geom_boxplot() +
  geom_point(pch = 21, position = position_jitterdodge(), alpha=0.6) +
  ylab("d13C host-symb") + theme_bw()

d13C.diff.species
ggsave("figures/d13C.diff.species.pdf", encod="MacRoman", height=5, width=8)


### difference of host and symbiont d15N
d15N.diff.species<-ggplot(dat, aes(x=Type, y=d15N.diff, fill =Species)) + 
  geom_boxplot() +
  geom_point(pch = 21, position = position_jitterdodge(), alpha=0.6) +
  ylab("d15N host-symb") + theme_bw()

d15N.diff.species
ggsave("figures/d15N.diff.species.pdf", encod="MacRoman", height=5, width=8)


## biplot of color d13C and d15N
host.biplot<-ggplot(dat, aes(x=d13C.host, y=d15N.host, color=Species)) +
  geom_point() +
  stat_ellipse()+
  theme_bw()

host.biplot

## biplot of symbiont d13C and d15N
alga.biplot<-ggplot(dat, aes(x=d13C.alg, y=d15N.alg, color=Species)) +
  geom_point() +
  stat_ellipse()+
  theme_bw()

alga.biplot

########### NMDS, select columns we want
nmds.df.full<- dat %>% 
  select(Island.atoll, Species, Site.name, Type, Lat, Lon, Depth..m, d15N.alg, d15N.host,
         d15N.diff, d13C.alg, d13C.host, d13C.diff, C.N.alga, C.N.host,
         d13C.POM, d15N.POM, Temperature, Salinity, SiO2, PO4, NO3, NH3, Chla..ugL, 
         TA..umolkg, DIC..umolkg)

#replace NA with 0.01 since below detection
nmds.df.full$NO3[is.na(nmds.df.full$NO3)] <- 0.01
nmds.df.full$SiO2[is.na(nmds.df.full$SiO2)] <- 0.01
nmds.df.full$PO4[is.na(nmds.df.full$PO4)] <- 0.01
nmds.df.full$NH3[is.na(nmds.df.full$NH3)] <- 0.01


# below has the most data possible, with replacement values
nmds.df<- nmds.df.full %>%
  select(Island.atoll, Species, Site.name, Type, Lat, Lon, Depth..m, d13C.alg, d13C.host, 
         d13C.diff, C.N.alga, C.N.host, d13C.POM, d15N.POM, Temperature, Salinity,
         Chla..ugL, TA..umolkg, DIC..umolkg, NO3, SiO2, PO4, NH3)

nmds.df.noNA<-na.omit(nmds.df)

# make POM +
nmds.df.noNA$d13C.POM<-abs(nmds.df.noNA$d13C.POM)

# NMDS dataframe
nmd<- nmds.df.noNA %>%
  select(d13C.alg, d13C.host, d13C.diff, C.N.alga, C.N.host)

# make + for d13C and add a number (7) to the host-symbiont diff
nmd$d13C.alg<-abs(nmd$d13C.alg)
nmd$d13C.host<-abs(nmd$d13C.host)
nmd$d13C.diff<-7 + nmd$d13C.diff

# continuous and categorical environmental data
env<- nmds.df.noNA %>%
  select(Species, Type, Lat, Lon, Depth..m, d13C.POM, d15N.POM,
         Temperature, Chla..ugL, NO3, SiO2, PO4, NH3, TA..umolkg, DIC..umolkg, Salinity)

# factors -- categorical
facs<- nmds.df.noNA %>%
  select(Island.atoll, Species, Site.name, Type)
  
# run NMDS, euclidean distance because different scales
NMDS<-metaMDS(nmd, distance='euclidean', k=2, trymax=100) 

# see stress plot, looks good
stressplot(NMDS)

# run envfit -- lat, long, depth, NO3, SiO2, TA, Salinity all significant
vars<-envfit(NMDS, env, permu=999)


###### build the dataframe
# gets NMDS 1 and 2 into a dataframe
data.scores = as.data.frame(scores(NMDS)$sites)

# extract vectors and factors from envfit
en_coord_cont = as.data.frame(scores(vars, "vectors")) * ordiArrowMul(vars)
en_coord_cat = as.data.frame(scores(vars, "factors")) * ordiArrowMul(vars)

# add back in site variables with NMDS1 and NMDS2
data.scores.df<-cbind(facs, data.scores)


gg = ggplot(data = data.scores.df, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores.df, aes(colour = Species), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c("orange", "steelblue", "coral","mediumseagreen")) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(en_coord_cat), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Species")

gg



# plot
par(mfrow=c(1,1))
plot(NMDS)
plot(vars, p.max=0.05, cex=0.8, lwd=1)

