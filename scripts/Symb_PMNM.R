rm(list=ls())
ls()

library(devtools)
library(plyr)
library(reshape2)
library(lme4)
library(steponeR)   # devtools::install_github("jrcunning/steponeR")


#### qPCR data ####

################################
## C and D mutiplexed plates ##
################################
CD.plates <- list.files(path="data/clade C D", pattern = "txt$", full.names = T); CD.plates
CD.qPCR <- steponeR(files=CD.plates, delim="\t",
                 target.ratios=c("C.D"),
                 fluor.norm=list(C=2.26827, D=0),
                 copy.number=list(C=33, D=3),
                 ploidy=list(C=1, D=1), 
                 extract=list(C=0.813, D=0.813))

CD.qPCR <- CD.qPCR$result
colnames(CD.qPCR)

# removes clade A on sample PMNM-333 (single sample run)
CD.qPCR<-CD.qPCR[!(CD.qPCR$Sample.Name=="A-86"),] # drop A sample
drop.A<-c("A.CT.mean","A.CT.sd","A.reps") # all A columns
CD.qPCR<-CD.qPCR[, (!names(CD.qPCR) %in% drop.A)] # drop columns by name
names(CD.qPCR) # check names for C and D only

head(CD.qPCR)

# to remove any early-amplification CT noise
CD.qPCR$C.CT.mean[which(CD.qPCR$C.CT.mean < 15)] <- 0
CD.qPCR$D.CT.mean[which(CD.qPCR$D.CT.mean < 15)] <- 0


# If C or D only detected in one technical replicate, set C:D ratio to...
CD.qPCR$C.D[which(CD.qPCR$D.reps==1)] <- 1 # where 1 = 100% C 
CD.qPCR$C.D[which(CD.qPCR$D.reps==0)] <- 1 # where 0 = 100% C 
CD.qPCR$C.D[which(CD.qPCR$C.reps==1)] <- 0 # where 1 = 100% D
CD.qPCR$C.D[which(CD.qPCR$C.reps==0)] <- 0 # where 0 = 100% D

#Remove failed samples, i.e., those where either C or D were NOT found in both reps
CD.qPCR$fail <- ifelse(CD.qPCR$C.reps < 2 & CD.qPCR$D.reps < 2, TRUE, FALSE)
fails <- CD.qPCR[CD.qPCR$fail==TRUE, ]
CD.qPCR <- CD.qPCR[which(CD.qPCR$fail==FALSE),] #dataframe now only "successes" of 2 reps

# replace CT means with 'NA' as zero
CD.qPCR$C.CT.mean[is.na(CD.qPCR$C.CT.mean)] <-0
CD.qPCR$D.CT.mean[is.na(CD.qPCR$D.CT.mean)] <-0

CD.qPCR$C.D[is.na(CD.qPCR$C.D)] <- 1 # sets all infinity (= 100% C) to 1.0

#---------------------------------------
#---------------------------------------
# calculate proportion C and proprtion D where C and D are both present
CD.qPCR$propC<- CD.qPCR$C.D / (CD.qPCR$C.D + 1)
CD.qPCR$propD<- 1 / (CD.qPCR$C.D + 1)

# where C and D are not cooccuring...
# if C.D = 1 = 100% C, make 'PropC' = 1 and 'PropD' = 0
# if C.D = 0 = 100% D, make 'PropD' = 1 and 'PropC' = 0

CD.qPCR$propC[which(CD.qPCR$C.D==1)] <- 1 #C.D = 1 is 100% and propC = 1
CD.qPCR$propD[which(CD.qPCR$propC==1)] <- 0 #C.D = 1, then D=0%, and propD = 0
CD.qPCR$propD[which(CD.qPCR$C.D==0)] <- 1 # C.D = 0, then 100% D, and propD = 1

#---------------------------------------
#---------------------------------------

# calculate FOUR COMMUNITY categories: C, C>D, D>C, D
CD.qPCR$syms <- factor(ifelse(CD.qPCR$propC > CD.qPCR$propD, ifelse(CD.qPCR$propD!= 0, "CD", "C"), ifelse(CD.qPCR$propD > CD.qPCR$propC, ifelse(CD.qPCR$propC!=0, "DC", "D"), NA)), levels=c("C", "CD", "DC", "D"))

# Identify SINGLE dominant symbiont clade: C or D
CD.qPCR$dom <- factor(substr(as.character(CD.qPCR$syms), 1, 1))

# Set zeros to NA to facilitate log transformation
CD.qPCR$propC[which(CD.qPCR$propC==0)] <- NA
CD.qPCR$propD[which(CD.qPCR$propD==0)] <- NA

CD.qPCR <- CD.qPCR[grep("cntl_1", CD.qPCR$Sample.Name, fixed=T, invert = T), ]
CD.qPCR <- CD.qPCR[grep("ddH20", CD.qPCR$Sample.Name, fixed=T, invert = T), ]
CD.qPCR <- CD.qPCR[grep("+C", CD.qPCR$Sample.Name, fixed=T, invert = T), ]
CD.qPCR <- CD.qPCR[grep("PMNM_+C", CD.qPCR$Sample.Name, fixed=T, invert = T), ]         
CD.qPCR <- CD.qPCR[grep("blank", CD.qPCR$Sample.Name, fixed=T, invert = T), ]  
CD.qPCR <- CD.qPCR[grep("cntl", CD.qPCR$Sample.Name, fixed=T, invert = T), ]  
CD.qPCR <- CD.qPCR[grep("+C_52", CD.qPCR$Sample.Name, fixed=T, invert = T), ]
CD.qPCR <- CD.qPCR[grep("H2O", CD.qPCR$Sample.Name, fixed=T, invert = T), ]

## see structure and break Sample.Name into Project = PMNM, and Sample.ID = number
df<-cbind(CD.qPCR, colsplit(CD.qPCR$Sample.Name, pattern= "_", c("Project", "Sample.ID")))
colnames(df)

# order columns 
CD.df<-df[, c(15:16, 1:14)]
colnames(CD.df)

# order all samples by Sample ID
CD.df<-CD.df[order(CD.df$Sample.ID), ]

#write.csv(CD.df, "output/precursor/PMNM CD qPCR_cleaned.csv")
head(CD.qPCR)


####### remove replicate samples from re-run plates
# look for duplicates in dataset
CD.df[duplicated(CD.df$Sample.Name), ]

# two part conditional statement to remove duplicates by ID and plate
CD.df<-CD.df[!(CD.df$Sample.ID=="4" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="24" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="28" & CD.df$File.Name=="Kelly_PMNM_plate11.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="28" & CD.df$File.Name=="Kelly_PMNM_plate7.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="32" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="63" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="72" & CD.df$File.Name=="Kelly_PMNM_plate4.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="72" & CD.df$File.Name=="Kelly_PMNM_plate11.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="73" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="82" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="82" & CD.df$File.Name=="Kelly_PMNM_plate4.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="104" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="107" & CD.df$File.Name=="Kelly_PMNM_plate6.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="107" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="123" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="141" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="303" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="311" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="314" & CD.df$File.Name=="Kelly_PMNM_plate4.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="314" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="317" & CD.df$File.Name=="Kelly_PMNM_plate6.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="317" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="319" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="319" & CD.df$File.Name=="Kelly_PMNM_plate8.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="331" & CD.df$File.Name=="Kelly_PMNM_plate2.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="343" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="343" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="353" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="356" & CD.df$File.Name=="Kelly_PMNM_plate6.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="361" & CD.df$File.Name=="Kelly_PMNM_plate11.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="361" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="362" & CD.df$File.Name=="Kelly_PMNM_plate11.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="362" & CD.df$File.Name=="Kelly_PMNM_plate2.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="362" & CD.df$File.Name=="Kelly_PMNM_plate10.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="365" & CD.df$File.Name=="Kelly_PMNM_plate5.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="385" & CD.df$File.Name=="Kelly_PMNM_plate5.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="30" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="35" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="97" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="106" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="172" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="175" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="195" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="220" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="228" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="240" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="251" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="271" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="333" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="337" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="358" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="360" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="363" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="366" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="372" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="386" & CD.df$File.Name=="Wall_PMNM_CD_plate14.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="28" & CD.df$File.Name=="Wall_PMNM_CD_plate15.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="362" & CD.df$File.Name=="Wall_PMNM_CD_plate15.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="15" & CD.df$File.Name=="Kelly_PMNM_plate4.txt"),]

#below are not duplicates, but amplification errors
CD.df<-CD.df[!(CD.df$Sample.ID=="392" & CD.df$File.Name=="Kelly_PMNM_plate3.txt"),]
CD.df<-CD.df[!(CD.df$Sample.ID=="396" & CD.df$File.Name=="Kelly_PMNM_plate6.txt"),]

CD.df[duplicated(CD.df$Sample.Name), ] # check that all duplicates removed.
#look for a sample of interest?
#CD.df[CD.df$Sample.Name=="PMNM_15",]

# which samples had D at all?
D.samples<-CD.df[CD.df$sym=="DC",] # 3 out of 4 from Laysan, 1 from FFS

# write.csv(CD.df, "output/precursor/PMNM CD qPCR_duplremov.csv")

##########################
# symb community table
CD.symb=table(CD.df$syms); CD.symb
CD.dom=table(CD.df$dom); CD.dom

# order all samples by dominant type
CD.qPCR<-CD.qPCR[order(CD.qPCR$dom), ]
# write.csv(CD.qPCR, "PMNM qPCR_dom.types.csv")


####################################
######## Clade A samples ###########
####################################

A.plates <- list.files(path="data/clade A", pattern = "txt$", full.names = T); A.plates
A.qPCR <- steponeR(files=A.plates, delim="\t", target.ratios = NULL, 
                   fluor.norm = NULL, copy.number = NULL, ploidy = NULL, extract = NULL)

# Clade A results
A.qPCR <- A.qPCR$result; head(A.qPCR)

# remove all controls, ddh20 samples
A.qPCR <- A.qPCR[grep("A-86", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("+C_86", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("+C_52", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("ddH2O", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("A-33", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("A-43", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("A-52", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("A-38", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("A-33", A.qPCR$Sample.Name, fixed=T, invert = T), ]
A.qPCR <- A.qPCR[grep("A-43", A.qPCR$Sample.Name, fixed=T, invert = T), ]

# If amplification NOT in BOTH replicate wells, give a "TRUE"
A.qPCR$fail <- ifelse(A.qPCR$A.reps < 2, TRUE, FALSE) # FALSE shows a "non-fail"

# to remove any early-amplification CT noise
A.qPCR$A.CT.mean[which(A.qPCR$A.CT.mean < 19)] <- NA

A.df<-cbind(A.qPCR, colsplit(A.qPCR$Sample.Name, pattern= "_", c("Project", "Sample.ID")))
A.df<-A.df[, c(7:8, 1:6)]

# order all samples by Sample ID
A.df<-A.df[order(A.df$Sample.ID), ]
head(A.df)

# write.csv(A.df, "output/precursor/PMNM A qPCR_cleaned.csv") # final output data with a FAIL=TRUE where no amplification

####### remove replicate samples from re-run plates
# look for duplicates in dataset
A.df[duplicated(A.df$Sample.Name), ]

# remove duplicates, or amplification errors
A.df<-A.df[!(A.df$Sample.ID=="38" & A.df$File.Name=="Kelly_PMNM_plateA8.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="122" & A.df$File.Name=="Kelly_PMNM_plateA6.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="199" & A.df$File.Name=="Kelly_PMNM_plateA9.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="223" & A.df$File.Name=="Kelly_PMNM_plateA2.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="285" & A.df$File.Name=="Kelly_PMNM_plateA7.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="285" & A.df$File.Name=="Kelly_PMNM_plateA9.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="306" & A.df$File.Name=="Kelly_PMNM_plateA3.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="309" & A.df$File.Name=="Kelly_PMNM_plateA9.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="357" & A.df$File.Name=="Kelly_PMNM_plateA7.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="375" & A.df$File.Name=="Kelly_PMNM_plateA3.txt"),]

## not duplicates, but errors
A.df<-A.df[!(A.df$Sample.ID=="391" & A.df$File.Name=="Kelly_PMNM_plateA1.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="395" & A.df$File.Name=="Kelly_PMNM_plateA3.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="399" & A.df$File.Name=="Kelly_PMNM_plateA5.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="196" & A.df$File.Name=="Kelly_PMNM_plateA8.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="308" & A.df$File.Name=="Kelly_PMNM_plateA4.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="369" & A.df$File.Name=="Kelly_PMNM_plateA5.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="377" & A.df$File.Name=="Kelly_PMNM_plateA4.txt"),]
A.df<-A.df[!(A.df$Sample.ID=="387" & A.df$File.Name=="Kelly_PMNM_plateA4.txt"),]

A.df[A.df$Sample.Name=="PMNM_370",]

# write.csv(A.df, "output/precursor/PMNM A qPCR_dupl remove.csv") # final output data with a FAIL=TRUE where no amplification


# if want to remove all samples where A was NOT found in both reps....
fails <- A.df[A.df$fail==TRUE, ]
A.df.both.reps <- A.df[which(A.df$fail==FALSE),]
# data file now is ONLY samples where A was found, and where A occurred in both technical reps
# "A.df.both.reps" needs QA/QC




################################
################################
# compare D samples and A samples
D.samples # recall, 1 FFS and 3 LAY
A.df.both.reps # MAR (2), KUR (1), LIS (12), LAY (5)

# import A dataframe, labeled
df.A<-read.csv("PMNM A qPCR_labeled.csv")
A.reps <- df.A[which(df.A$fail==FALSE),]
df.table.A.allcorals=table(df.A$Island) # number of corals sampled
df.table.Apres=table(df.A$Island, df.A$fail) # number of corals with A (fail) and without
A.prop.table<-round(prop.table(df.table.Apres,1), 2) # proportion table for corals A/site. The second argument specifies the total for each should sum to "1" and round to 2 digits

# import C/D dataframe labeled
df.CD<-read.csv("PMNM CD qPCR_labeled.csv")
df.table.CD.allcorals=table(df.CD$Island) # number of corals sampled
df.table.CDpres=table(df.CD$Island, df.CD$dom) # corals with CD (fail) and without
CD.prop.table<-round(prop.table(df.table.CDpres,1), 2)


##### combine dataframes #####
# note that C and D on same, multiplexed run. A alone. At this point the # of reps is not equal due to a need for 100% of frags to be assessed for qPCR. Still gives good approximate.

df.ACD<-as.data.frame.matrix(A.prop.table)
df.ACD$C<-CD.prop.table[,1]
df.ACD$D.C<-CD.prop.table[,2]
colnames(df.ACD)<-c("C.w.A", "C.no.A", "C.no.D", "D.C")
# C.w.A =  corals with A, but all showed C, so this is C and A together from different runs
# C.no.A = corals from A runs that did not show A--but all had  C from multiplexed run
# C.no.D = corals with C alone, no D, but some may have been C+A
# D.C = corals where D was observed, always in presence of C, but D in high abundance

# order prop.table and make df
df.ACD.ordered<-as.data.frame.matrix(df.ACD[c(2,6,7,4,3,5,1),])
df.ACD.ordered$Island<-factor(cbind("KUR", "MID", "PHA", "LIS", "LAY", "MAR", "FFS")); 
df.ACD.ordered$C<-1-(df.ACD.ordered$C.w.A+df.ACD.ordered$D.C)

############
m<-layout(matrix(1:2), 1, 1, respect = FALSE)
layout.show(m)
Sites.ordered<-c("KUR", "MID", "PHA", "LIS", "LAY", "MAR", "FFS")

########## separate plots
####### Figure of corals with C dominant
par(mfrow=c(2,1), mar=c(4,4,0.5,4), xpd=TRUE)
plot(df.ACD.ordered$C, xaxt="n", ylab="proportion of corals", xlab="", type="b", pch=20, lwd=3, ylim=c(0.6, 1.0), col="darkseagreen3", lty="dotted", cex=2)
axis(1, at=1:7, labels=Sites.ordered)
legend("bottomleft", horiz=FALSE, title="symbiont clade", col=c("darkseagreen3", "coral", "cadetblue3"), legend=c("C only","D>C", "C+A"), pch=c(19, 1, 19), lty="dotted", lwd=3, bty="n", inset=c(0,0), x.intersp=0.9, y.intersp=0.5)

plot(df.ACD.ordered$C.w.A, xaxt="n", ylab="proportion of corals", xlab="Island or Atoll--North (L) to South (R)", type="b", pch=20, lwd=3, ylim=c(0, 0.3), col="cadetblue3", lty="dotted", cex=2)
axis(1, at=1:7, labels=Sites.ordered)
with(points(df.ACD.ordered$D.C, xaxt="n", col="coral", type="b", lwd=3, pch=1, lty="dotted", cex=1.5))

############# combined plot
####### Figure of corals with C dominant
par(mfrow=c(1,1), mar=c(4,4,0.5,4), xpd=TRUE)
plot(df.ACD.ordered$C, xaxt="n", ylab="proportion of corals", xlab="Island or Atoll--North (L) to South (R)", type="b", pch=20, lwd=3, ylim=c(0.0, 1.0), col="darkseagreen3", lty="dotted", cex=3, cex.axis=1.3, cex.lab=1.3)
axis(1, at=1:7, labels=Sites.ordered, cex.axis=1.3, cex=1.2)
legend("left", horiz=FALSE, title="symbionts", col=c("darkseagreen3", "coral", "cadetblue3"), legend=c("C only","D+C", "C+A"), pch=c(19, 1, 19), cex=1.5, lty="dotted", lwd=3, bty="n", inset=c(0,0), x.intersp=0.9, y.intersp=0.5)
with(points(df.ACD.ordered$C.w.A, xaxt="n", type="b", pch=20, lwd=3, col="cadetblue3", lty="dotted", cex=3))
#axis(1, at=1:7, labels=Sites.ordered)
with(points(df.ACD.ordered$D.C, xaxt="n", col="coral", type="b", lwd=3, pch=1, lty="dotted", cex=2))

dev.copy(pdf, "CDA plot2.pdf", width=6, height=5)
dev.off()


##########
df.CD<-read.csv("PMNM CD qPCR_labeled.csv") # labeled dataframe
df.CD$Depth..m<-(df.CD$Depth..ft*0.348) # converts ft. to m.

#########
par(mfrow=(c(1,1)), mar=c(5,5,2,1))
#########
# depths for each site--MID and MAR are noticably shallow relative to oother sites

depth.df<-aggregate(Depth..m~Site+Island+Region, data=df.CD, FUN=mean) # dataframe for PMNM
levels(depth.df$Island)
depth.df$Island<-factor(depth.df$Island,levels(depth.df$Island)[c(2,6,7,4,3,5,1)])
plot(depth.df$Depth..m~depth.df$Island, main="", xlab="Island or Atoll--North (L) to South (R)", ylab="Depth (m)", col="lightblue", cex=1, cex.axis=1.5, cex.lab=1.5) # depths across PMNM

dev.copy(pdf, "Depthplot.pdf", width=8, height=5)
dev.off()


#########
#########

df.CD.depth<-df.CD[(df.CD$dom=="D"),]
df.CD.depth$Island<-factor(df.CD.depth$Island,levels(df.CD.depth$Island)[c(2,6,7,4,3,5,1)])
plot(df.CD.depth$Depth..m~df.CD.depth$Island, ylim=c(0,17), ylab="Depth (m)", xlab="Island or Atoll--North (L) to South (R)")

levels(df.CD.depth$Region)
df.CD.depth$Region<-factor(df.CD.depth$Region,levels(df.CD.depth$Region)[c(3,1,2)])
plot(df.CD.depth$Depth..m~df.CD.depth$Region, ylim=c(0,17),  ylab="Depth (m)", xlab="Geographic Region") # depth from 3 - 50 ft

###########
df.A<-read.csv("PMNM A qPCR_labeled.csv")
df.A # labeled dataframe
df.A$Depth..m<-(df.A$Depth..ft*0.348) # converts ft. to m.
df.A$Island<-factor(df.A$Island,levels(df.A$Island)[c(2,6,7,4,3,5,1)])
levels(df.A$Island)

df.A.depth<-df.A[(df.A$fail=="FALSE"),]
levels(df.A.depth$Island)
plot(df.A.depth$Depth..m~df.A.depth$Island, ylim=c(0,17), ylab="Depth (m)", xlab="Island or Atoll--North (L) to South (R)")

df.A.depth$Region<-factor(df.A.depth$Region,levels(df.A.depth$Region)[c(3,1,2)])
plot(df.A.depth$Depth..m~df.A.depth$Region, ylim=c(0,17),  ylab="Depth (m)", xlab="Geographic Region") # depth from 3 - 50 ft


Islands=c(1:7)
clades<-c("D>C", "C+A")
colors=c("coral", "cadetblue3")
####### Figure of corals with A and D across depths and atolls 
par(mfrow=c(2,1), mar=c(4,5,0.5,5), xpd=TRUE)
m<-layout(matrix(1:2), 1, 1, respect = FALSE)
plot(df.CD.depth$Depth..m~df.CD.depth$Island, ylim=c(0,17), xaxt="n", ylab="Depth (m)", xlab="", col="coral")
axis(side=1, at=Islands, cex.axis=0.1)
legend("topleft", legend=clades, col=colors, pch=15, pt.cex=2.5, cex=1, bty="n", x.intersp=0.4, y.intersp=0.6, title="symbiont clades")
plot(df.A.depth$Depth..m~df.A.depth$Island, ylim=c(0,17), ylab="Depth (m)", xlab="Island or Atoll", col="cadetblue3")


##############
df.ACD.ordered

par(mfrow=c(1,1))
slices<-c(361,4)
lab<-c("C", "D+C")
pie(slices, labels=lab, main="Symbiont community", col=c("cadetblue2", "coral"))

pie.df<-df.ACD.ordered[,c(1:4)] # proportions
isl<-df.ACD.ordered[,5] # names
labs<- c("C+A", "C", "D+C")

#Kure
slices<-c(3, 97, 100, 0)
pie(slices, labels=labs, main="Symbiont community", col=c("cadetblue3", "darkseagreen3", "coral"), main="Kure"))


##### chi-square tests

# df.ACD.ordered combined df
df.table.Apres
df.table.CDpres
chisq.test(df.ACD.ordered$C.w.A, df.ACD.ordered$Region)

# Region
df.ACD.ordered
df.ACD.ordered$Region<-c("North", "North", "North", "Central", "Central", "Central", "South")

chisq.test(df.ACD.ordered$C.w.A, df.ACD.ordered$Region)
chisq.test(df.ACD.ordered$C.w.A, df.ACD.ordered$Island)


# Site
table.CD=table(CD.data$Island, CD.data$dom) # corals with CD (fail) and without
CD.prop.table<-round(prop.table(table.CD,1), 2)
chisq.test(CD.prop.table)


##### chi-square tests
# Region
table.A=table(A.data$Region, A.data$fail) # corals with CD (fail) and without
df<-as.data.frame.matrix(table.A)
df$Region<-c("North", "Central", "South")
colnames(df)<-c("A", "w.o A", "Region")
A.prop.table<-round(prop.table(table.A,1), 2)


# Site
table.A=table(A.data$Island, A.data$fail) # corals with CD (fail) and without
A.prop.table<-round(prop.table(table.A,1), 2)


chisq.test(df.ACD.ordered$C.w.A, df.ACD.ordered$Region)
chisq.test(df.ACD.ordered$C.w.A, df.ACD.ordered$Island)

