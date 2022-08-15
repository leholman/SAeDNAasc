#############################################
####== South African Ascidian Analysis ==####
####==== Luke E. Holman====08.02.2021====####
#############################################

###Script 1 - Analysis###

####====0.0 Packages====######

library("reshape")
library("seqinr")
library("adegenet")
library("haplotypes")
library("ape")
library("RColorBrewer")
library("pegas")
library("hierfstat")
library("poppr")
library("maps")
library("prettymapr")
library("sp")

#read in distance
dist <- read.csv("data/distance.csv")

####====0.1 Map====######
pdf("figures/figure1/map.pdf",width=8,height=5)
m<- map("world", xlim=c(16,33), ylim=c(-35,-28), col="gray", fill=TRUE,mar=c(4, 4, 4, 4))
addnortharrow("topleft",scale=0.75)
addscalebar(pos="topright")


xat <- pretty(m$range[1:2],n = 10)
xlab <- parse(text=degreeLabelsEW(xat))

yat <- pretty(m$range[3:4])
ylab <- parse(text=degreeLabelsNS(yat))

box()
axis(2, las=TRUE, at=yat, labels=ylab)
axis(3, at=xat, labels=xlab)
#axis(4, las=TRUE, at=yat, labels=ylab)

points(dist$Long,dist$Lat,pch=16,cex=1.7,col="black")
points(dist$Long,dist$Lat,pch=16,cex=1.5,col="white")
points(dist$Long,dist$Lat,pch=16,cex=1,col="darkblue")


dev.off()


####====1.0 Ascidian Distribution====######

##Here we first subset the data to include only the four main species we are interested in 
##1.1 Rapid Assessment 
focus.species <- c("Ciona robusta","Clavelina lepadiformis","Microcosmus squamiger","Styela plicata")
RA.s.09 <- read.csv(file="data/distribution/09.RA.csv",header=TRUE,stringsAsFactors = T)
RA.s.09 <- RA.s.09[RA.s.09$Species %in% focus.species,]
RA.s.09$Species <- droplevels(RA.s.09$Species)
RA.s.09 <- melt(RA.s.09, id=c("Species"))
RA.s.17 <- read.csv("data/distribution/17.RA.csv",header=TRUE,stringsAsFactors = T)
RA.s.17 <- RA.s.17[RA.s.17$Species %in% focus.species,]
RA.s.17$Species <- droplevels(RA.s.17$Species)
RA.s.17 <- melt(RA.s.17, id=c("Species"))


##RABS sitewise 09 & 17

pdf("figures/RA.09-17.s.pdf",width = 8,height = 2.5)
par(mar=c(1.1,11.1,3.1,2.1))
plot(as.numeric(RA.s.09$variable),as.numeric(RA.s.09$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n",ylim=c(0.5,4.5))
for (line in 1:12){
  segments(line,0,line,10,lty=5,col="grey")
}
points(as.numeric(RA.s.09$variable)-0.2,as.numeric(factor(RA.s.09$Species, levels=rev(levels(RA.s.09$Species)))),cex=(as.numeric(RA.s.09$value))*1.2,pch=19,col="darkred")
points(as.numeric(RA.s.09$variable[RA.s.09$value=="NS"])-0.2,as.numeric(factor(RA.s.09$Species[RA.s.09$value=="NS"], levels=rev(levels(RA.s.09$Species[RA.s.09$value=="NS"])))),cex=2*1.2,pch=1,col="darkred")
points(as.numeric(RA.s.17$variable)+0.2,as.numeric(factor(RA.s.17$Species, levels=rev(levels(RA.s.17$Species)))),cex=(as.numeric(RA.s.17$value))*1.2,pch=19,col="darkblue")
axis(3,at=1:12,labels=levels(RA.s.09$variable))
axis(2,at=1:4,labels=rev(levels(RA.s.09$Species)),las=1,font.axis=3)

dev.off()



#RABS by distance 09 & 17
distances <- dist$DistFromWest[match(as.character(RA.s.09$variable),dist$sitecode)]
pdf("figures/figure1/RA.09-17.s.range.pdf",width = 8,height = 2.5)
par(mar=c(3.1,11.1,2.1,2.1))
plot(distances,as.numeric(RA.s.09$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n",ylim=c(0.5,4.5))


counter <- 0
for (line in levels(RA.s.09$variable)){
  segments(counter*183.3333333333,4.5,counter*183.3333333333,5,lty=5,col="grey")
  segments(dist$DistFromWest[dist$sitecode==line],4.2,counter*183.3333333333,4.5,lty=5,col="grey")
  segments(dist$DistFromWest[dist$sitecode==line],0,dist$DistFromWest[dist$sitecode==line],4.2,lty=5,col="grey")
  counter <- counter+1
}


for (species in levels(RA.s.09$Species)){
  min09 <- min(distances[RA.s.09$value>0 & RA.s.09$value!="NS" & RA.s.09$Species==species ])
  max09 <- max(distances[RA.s.09$value>0 & RA.s.09$value!="NS" & RA.s.09$Species==species ])
  # min09ns <- min(distances[RA.s.09$value>0  & RA.s.09$Species==species ])
  # max09ns <- max(distances[RA.s.09$value>0  & RA.s.09$Species==species ])
  min17 <- min(distances[RA.s.17$value>0 & RA.s.17$value!="NS" & RA.s.17$Species==species])
  max17 <- max(distances[RA.s.17$value>0 & RA.s.17$value!="NS" & RA.s.17$Species==species])
  speciesindex <- match(species,rev(levels(RA.s.09$Species)))
  #  segments(min09ns,speciesindex+0.1,max09ns,speciesindex+0.1,col="darkred",lwd=3,lty=5)
  segments(min09,speciesindex+0.1,max09,speciesindex+0.1,col="darkred",lwd=5)
  segments(min17,speciesindex-0.1,max17,speciesindex-0.1,col="darkblue",lwd=5)
  #add a dotted line for species found in DU in 09 to indicate non-survey
  #if(max(distances[RA.s.09$value>0 & RA.s.09$value!="NS" & RA.s.09$Species==species ])>1800){
  #  segments(1879.18,speciesindex+0.1,2047.53,speciesindex+0.1,col="darkred",lwd=5,lty=3)
  #}
}

axis(3,at=(0:11)*183.3333333333,labels=levels(RA.s.09$variable))
axis(1,at=seq(0,2200,200),labels=seq(0,2200,200),cex.axis=0.8)
axis(2,at=1:4,labels=rev(levels(RA.s.09$Species)),las=1,font.axis=3)
dev.off()

#1.2 What about eDNA surveys? (dropping site SY as not covered in other surveys)
#COI
eDNA.coi <- read.csv("data/eDNA/cleaned/COI.incidence.csv")
eDNA.coi <- eDNA.coi[-c(1,15)]
eDNA.coi <- eDNA.coi[eDNA.coi$Assignment %in% focus.species,]
eDNA.coi <- melt(eDNA.coi, id=c("Assignment"))
eDNA.coi <- eDNA.coi[eDNA.coi$variable!="SY",]

eDNA.coi.u <- read.csv("data/eDNA/cleaned/COI.u.incidence.csv")
eDNA.coi.u <- eDNA.coi.u[-c(1,15)]
eDNA.coi.u <- eDNA.coi.u[eDNA.coi.u$Assignment %in% focus.species,]
eDNA.coi.u <- melt(eDNA.coi.u, id=c("Assignment"))
eDNA.coi.u <- eDNA.coi.u[eDNA.coi.u$variable!="SY",]

#Calculate mean of DADA2 and UNOISE3 detection because there are some small differences! 
eDNA.coi$value <- (eDNA.coi$value + eDNA.coi.u$value[match(paste(eDNA.coi.u$Assignment,eDNA.coi.u$variable,sep="."),paste(eDNA.coi$Assignment,eDNA.coi$variable,sep="."))])/2


#18S
eDNA.18S <- read.csv("data/eDNA/cleaned/18S.incidence.csv")
eDNA.18S <- eDNA.18S[-c(1,15)]
eDNA.18S <- eDNA.18S[eDNA.18S$Assignment %in% focus.species,]
eDNA.18S <- melt(eDNA.18S, id=c("Assignment"))
eDNA.18S <- eDNA.18S[eDNA.18S$variable!="SY",]

#Orders for plot
speciesOrder <- rev(levels(RA.s.09$Species))
siteOrder <- levels(RA.s.09$variable)
binnedReads <- cut(eDNA.coi$value,breaks = c(0,0.000001,0.00001,0.0001,0.001,0.01),right=F,labels=c("0","1","2","3","4"))
binnedReads18S <- cut(eDNA.18S$value,breaks = c(0,0.00001,0.0001,0.001,0.01,1),right=F,labels=c("0","1","2","3","4"))

#plot (just eDNA)
pdf("figures/eDNA.17.BS.pdf",width = 8,height = 2.5)
par(mar=c(1.1,11.1,3.1,2.1))
plot(as.numeric(RA.s.09$variable),as.numeric(RA.s.09$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n")
for (line in 1:12){
  segments(line,0,line,10,lty=5,col="grey")
}
points(match(eDNA.coi$variable,siteOrder),match(eDNA.coi$Assignment,speciesOrder),cex= as.numeric(as.character(binnedReads)),pch=0,col="darkgreen")
points(match(eDNA.18S$variable,siteOrder),match(eDNA.18S$Assignment,speciesOrder),cex= as.numeric(as.character(binnedReads18S)),pch=5,col="darkviolet")

axis(3,at=1:12,labels=levels(RA.s.09$variable))
axis(2,at=1:4,labels=rev(levels(RA.s.09$Species)),las=1)

legend("topright",c("COI","18S"),pch=c(0,5),col=c("darkgreen","darkviolet"),box.lty = 0)
dev.off()


#Combined plot including only the eDNA and 2017 data 


pdf("figures/figure1/RA.17.eDNA.s.pdf",width = 8,height = 2.5)
par(mar=c(1.1,11.1,3.1,2.1))
plot(as.numeric(RA.s.09$variable),as.numeric(RA.s.09$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n",ylim=c(0.5,4.5))
for (line in 1:12){
  segments(line,0,line,10,lty=5,col="grey")
}
points(match(eDNA.coi$variable,siteOrder)-0.2,match(eDNA.coi$Assignment,speciesOrder),cex= as.numeric(as.character(binnedReads)),pch=0,col="darkgreen")
points(match(eDNA.18S$variable,siteOrder)-0.2,match(eDNA.18S$Assignment,speciesOrder),cex= as.numeric(as.character(binnedReads18S)),pch=5,col="darkviolet")
points(as.numeric(RA.s.17$variable)+0.2,as.numeric(factor(RA.s.17$Species, levels=rev(levels(RA.s.17$Species)))),cex=(as.numeric(RA.s.17$value))*1.2,pch=19,col="darkblue")
axis(1,at=1:12,)
axis(2,at=1:4,labels=rev(levels(RA.s.09$Species)),las=1,font.axis=3)
#legend("topright",c("Rapid Assessment","eDNA COI","eDNA 18S"),pch=c(16,0,5),col=c("darkblue","darkgreen","darkviolet"),box.lty = 0,cex=0.6)
dev.off()



#Legend for bubbles
pdf("figures/figure1/legend.pdf",width = 1.5,height = 2)
par(mar=c(1.1,1.1,1.1,1.1))
plot(NA,xlim=c(1,3),ylim=c(1,4),xaxt="n",yaxt="n",bty="n")
points(rep(1.5,3),seq(3.5,3.75,0.09),cex=c(1.2,2.4,3.6),pch=1,col="darkblue")
text(rep(3.2,3),seq(3.5,3.9,length.out=3),labels=c("< 10%","10-50%","> 50%"),cex=0.6,pos=2,col="darkblue")
points(rep(1.5,4),seq(2.5,2.75,0.065),cex=1:4,pch=0,col="darkgreen")
text(rep(3.2,4),seq(2.4,2.9,length.out=4),labels=c("> 0.000001","> 0.00001","> 0.0001","> 0.001"),cex=0.6,pos=2,col="darkgreen")
points(rep(1.5,4),seq(1.5,1.75,0.065),cex=1:4,pch=5,col="darkviolet")
text(rep(3.2,4),seq(1.4,1.9,length.out=4),labels=c("> 0.00001","> 0.0001","> 0.001","> 0.01"),cex=0.6,pos=2,col="darkviolet")

dev.off()


#Legend for dates
pdf("figures/figure1/legend.dates.pdf",width = 1.5,height = 2)
par(mar=c(1.1,1.1,1.1,1.1))
plot(NA,xlim=c(1,3),ylim=c(1,4),xaxt="n",yaxt="n",bty="n")
text(2.5,3.5,pos=2,"2017",col="darkblue")
segments(1.5,3.2,2.3,3.2,lwd=5,col="darkblue")
text(2.5,2.5,pos=2,"2009",col="darkred")
segments(1.5,2.2,2.3,2.2,lwd=5,col="darkred")

dev.off()



##1.1Rapid Assessment Results including additional species 
RA.09 <- read.csv(file="data/distribution/09.RA.csv",header=TRUE,stringsAsFactors = T)
RA.09 <- melt(RA.09, id=c("Species"))
RA.17 <- read.csv("data/distribution/17.RA.csv",header=TRUE,stringsAsFactors = T)
RA.17 <- melt(RA.17, id=c("Species"))


#RABS
pdf("figures/RA.09-17.pdf",width = 8,height = 5)
par(mar=c(1.1,11.1,3.1,2.1))
plot(as.numeric(RA.09$variable),as.numeric(RA.09$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n")
for (line in 1:12){
  segments(line,0,line,10,lty=5,col="grey")
}
points(as.numeric(RA.09$variable)-0.2,as.numeric(factor(RA.09$Species, levels=rev(levels(RA.09$Species)))),cex=(as.numeric(RA.09$value))*1.2,pch=19,col="darkred")
points(as.numeric(RA.09$variable[RA.09$value=="NS"])-0.2,as.numeric(factor(RA.09$Species[RA.09$value=="NS"], levels=rev(levels(RA.09$Species[RA.09$value=="NS"])))),cex=2*1.2,pch=1,col="darkred")
points(as.numeric(RA.17$variable)+0.2,as.numeric(factor(RA.17$Species, levels=rev(levels(RA.17$Species)))),cex=(as.numeric(RA.17$value))*1.2,pch=19,col="darkblue")
axis(3,at=1:12,labels=levels(RA.09$variable))
axis(2,at=1:9,labels=rev(levels(RA.09$Species)),las=1)

dev.off()

#RABS by distance
distances <- dist$DistFromWest[match(as.character(RA.09$variable),dist$sitecode)]
pdf("figures/RA.09-17.range.pdf",width = 8,height = 3)
par(mar=c(3.1,11.1,1.1,2.1))
plot(distances,as.numeric(RA.09$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n")
for (line in seq(0,2200,200)){
  segments(line,0,line,10,lty=5,col="grey")
}
for (species in levels(RA.09$Species)){
  min09 <- min(distances[RA.09$value>0 & RA.09$value!="NS" & RA.09$Species==species ])
  max09 <- max(distances[RA.09$value>0 & RA.09$value!="NS" & RA.09$Species==species ])
 # min09ns <- min(distances[RA.09$value>0  & RA.09$Species==species ])
 # max09ns <- max(distances[RA.09$value>0  & RA.09$Species==species ])
  min17 <- min(distances[RA.17$value>0 & RA.17$value!="NS" & RA.17$Species==species])
  max17 <- max(distances[RA.17$value>0 & RA.17$value!="NS" & RA.17$Species==species])
  speciesindex <- match(species,rev(levels(RA.09$Species)))
#  segments(min09ns,speciesindex+0.1,max09ns,speciesindex+0.1,col="darkred",lwd=3,lty=5)
  segments(min09,speciesindex+0.1,max09,speciesindex+0.1,col="darkred",lwd=3)
  segments(min17,speciesindex-0.1,max17,speciesindex-0.1,col="darkblue",lwd=3)
  #add a dotted line for species found in DU in 09 to indicate non-survey
  if(max(distances[RA.09$value>0 & RA.09$value!="NS" & RA.09$Species==species ])>1800){
    segments(1879.18,speciesindex+0.1,2047.53,speciesindex+0.1,col="darkred",lwd=3,lty=3)
  }
  
}
axis(1,at=seq(0,2200,200),labels=seq(0,2200,200),cex.axis=0.8)
axis(2,at=1:9,labels=rev(levels(RA.09$Species)),las=1)
dev.off()

##Historical data 

historical <- read.csv("data/historical/RawData.csv")

historical <- read.csv("data/historical/rawdata.update2022.csv")



pdf("figures/HistroicalRangeNew.pdf",width = 4,height = 7)
palette(c("#d55e00","#cc79a7","#0072b2","#009e73"))
par(mfrow=c(4,1),mar=c(2.1, 2.1, 1, 1))
counter <-  1
for (spp in unique(historical$Species)){
  plot(0,xlim=c(1935,2020),ylim=c(0,2000),cex.axis=0.92)
  polygon(c(historical$Date[historical$Species==spp],rev(historical$Date[historical$Species==spp])),
          c(historical$Range[historical$Species==spp],rep(0,length(historical$Range[historical$Species==spp]))),
          col=counter,border = NA)
  #text(1945,1500,bquote(italic(.(spp))),pos=3)
  counter <- counter + 1
}
dev.off()


#1.2 What about eDNA surveys? (dropping site SY as not covered in other surveys)
#COI
eDNA.coi <- read.csv("data/eDNA/cleaned/COI.incidence.csv")
eDNA.coi <- eDNA.coi[-c(1,15)]
eDNA.coi <- eDNA.coi[eDNA.coi$Assignment %in% levels(RA.09$Species),]
eDNA.coi <- melt(eDNA.coi, id=c("Assignment"))
eDNA.coi <- eDNA.coi[eDNA.coi$variable!="SY",]

#18S
eDNA.18S <- read.csv("data/eDNA/cleaned/18S.incidence.csv")
eDNA.18S <- eDNA.18S[-c(1,15)]
eDNA.18S <- eDNA.18S[eDNA.18S$Assignment %in% levels(RA.09$Species),]
eDNA.18S <- melt(eDNA.18S, id=c("Assignment"))
eDNA.18S <- eDNA.18S[eDNA.18S$variable!="SY",]

#Orders for plot
speciesOrder <- rev(levels(RA.09$Species))
siteOrder <- levels(RA.09$variable)
binnedReads <- cut(eDNA.coi$value,breaks = c(0,0.000001,0.00001,0.0001,0.001,0.01),right=F,labels=c("0","1","2","3","4"))
binnedReads18S <- cut(eDNA.18S$value,breaks = c(0,0.00001,0.0001,0.001,0.01,1),right=F,labels=c("0","1","2","3","4"))

#plot
pdf("figures/eDNA.17.BS.pdf",width = 8,height = 5)
par(mar=c(1.1,11.1,3.1,2.1))
plot(as.numeric(RA.09$variable),as.numeric(RA.09$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n")
for (line in 1:12){
  segments(line,0,line,10,lty=5,col="grey")
}
points(match(eDNA.coi$variable,siteOrder),match(eDNA.coi$Assignment,speciesOrder),cex= as.numeric(as.character(binnedReads)),pch=0,col="darkgreen")
points(match(eDNA.18S$variable,siteOrder),match(eDNA.18S$Assignment,speciesOrder),cex= as.numeric(as.character(binnedReads18S)),pch=5,col="darkviolet")

axis(3,at=1:12,labels=levels(RA.09$variable))
axis(2,at=1:9,labels=rev(levels(RA.09$Species)),las=1)

legend("topright",c("COI","18S"),pch=c(0,5),col=c("darkgreen","darkviolet"),box.lty = 0)
dev.off()

##Add in pyura for comparison
eDNA.TEMP <- read.csv("data/eDNA/cleaned/18S.incidence.csv")
eDNA.TEMP <- eDNA.TEMP[-c(1,15)]
eDNA.py <- eDNA.TEMP[grep("Pyura",eDNA.TEMP$Assignment),]
eDNA.py <- eDNA.py[,-12]
eDNA.py <- melt(eDNA.py, id=c("Assignment"))
binnedReads.py <- cut(eDNA.py$value,breaks = c(0,0.00001,0.0001,0.001,0.01,1),right=F,labels=c("0","1","2","3","4"))


pdf("figures/eDNA.Py.pdf",width = 8,height = 2)
par(mar=c(1.1,9.1,3.1,2.1))
plot(as.numeric(eDNA.py$variable),as.numeric(as.factor(eDNA.py$Assignment)),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n",ylim=c(0,2.75))
for (line in 1:12){
  segments(line,-5,line,10,lty=5,col="grey")
}
points(match(eDNA.py$variable,siteOrder),as.factor(eDNA.py$Assignment),cex= as.numeric(as.character(binnedReads.py)),pch=16,col="darkgreen")

axis(3,at=1:12,labels=siteOrder)
axis(2,at=1:2,labels=levels(as.factor(eDNA.py$Assignment)),las=1)

dev.off()

##Compare abundance measures between eDNA and RA
#18S
Data <- paste(RA.17$Species,RA.17$variable)
eDNA0 <- eDNA.18S$value>0
RA0 <- RA.17$value[match(paste(eDNA.18S$Assignment,eDNA.18S$variable),Data)]>0
RAdat <- RA.17$value[match(paste(eDNA.18S$Assignment,eDNA.18S$variable),Data)]
RAdat18 <- RAdat
NonZero <- eDNA0 & RA0
NonZero18 <-NonZero

#plot(RAdat[NonZero],log10(eDNA.18S$value[NonZero]))
summary(lm(log10(eDNA.18S$value[NonZero])~RAdat[NonZero]))


#COI
eDNA0 <- eDNA.coi$value>0
RA0 <- RA.17$value[match(paste(eDNA.coi$Assignment,eDNA.coi$variable),Data)]>0
RAdat <- RA.17$value[match(paste(eDNA.coi$Assignment,eDNA.coi$variable),Data)]
NonZero <- eDNA0 & RA0


#plot(RAdat[NonZero],log10(eDNA.coi$value[NonZero]))
summary(lm(log10(eDNA.coi$value[NonZero])~RAdat[NonZero]))


#Plots
pdf("figures/eDNA.RA.Quant.pdf",width = 8,height = 4)
par(mar=c(5.1, 4.1, 4.1, 2.1),mfrow=c(1,2))
plot(RAdat18[NonZero18],log10(eDNA.18S$value[NonZero18]),ylab="log10 eDNA Read Proportion",xlab="Rapid Assessment Abundance",pch=16,main="18S")
abun <- RAdat18[NonZero18]
reads <- log10(eDNA.18S$value[NonZero18])
#predictions <- predict( lm(reads~abun),data.frame(abun=seq(0.5,3.5,0.01)),interval = "confidence")
#lines(seq(0.5,3.5,0.01),predictions[,1],col="red")
#lines(seq(0.5,3.5,0.01),predictions[,2],col="red",lty=3)
#lines(seq(0.5,3.5,0.01),predictions[,3],col="red",lty=3)
plot(RAdat[NonZero],log10(eDNA.coi$value[NonZero]),ylab="log10 eDNA Read Proportion",xlab="Rapid Assessment Abundance",pch=16,main="COI")
dev.off()

####====2.0 Pop Gen w/ Sanger====######

##How much data do we have from our sanger seqs
sanger.dat <- read.csv("data/SangerData.csv")
table(paste0(as.character(sanger.dat$Sample.Species),".",as.character(sanger.dat$Sample.Site)))

##2.1 Summary


#Make dataframe for summary stats

for (species in c("Ciona.fa","Clav.fa","Micro.fa",'Styela.fa')){
  

critter <- gsub(".fa","",species)


HaploSummary <- matrix(data=NA,nrow=length(dist$sitecode),ncol=6)
row.names(HaploSummary) <- dist$sitecode
colnames(HaploSummary) <- c("NumH.09","NumH.17","DivH.09","DivH.17","DivN.09","DivN.17")

#Calculate Haplotype richness

HaploSummary[,1] <- colSums( read.csv(paste0("data/sanger/HaploDistribution/",critter,"09.csv"),row.names = 1)[,1:12]!= 0)
HaploSummary[,2] <- colSums( read.csv(paste0("data/sanger/HaploDistribution/",critter,"17.csv"),row.names = 1)[,1:12]!= 0)


#Calculate haplotype / nucleotide diversity 

LoopDat <- read.dna(paste0('data/sanger/Alignments/',critter,".fa"), "fasta")
location <- sapply(strsplit(labels(LoopDat),"_"),"[[",2)
year <- sapply(strsplit(labels(LoopDat),"_"),"[[",1)

unique(paste(location,year,sep="_"))
match(paste(location,year,sep="_"),unique(paste(location,year,sep="_")))

for (item in unique(paste(location,year,sep="_"))){
  loopingsite <- unlist(strsplit(item,"_"))[1]
  loopingyear <- unlist(strsplit(item,"_"))[2]
  if (loopingyear=="09"){
    HaploSummary[loopingsite,3] <- paste(round(hap.div(LoopDat[which(paste(location,year,sep="_") == item),],variance=TRUE)[1],3),
                                         round(hap.div(LoopDat[which(paste(location,year,sep="_") == item),], variance=TRUE)[2],3),sep = "±") 
    HaploSummary[loopingsite,5] <- paste(round(nuc.div(LoopDat[which(paste(location,year,sep="_") == item),],variance=TRUE)[1],5),
                                         round(nuc.div(LoopDat[which(paste(location,year,sep="_") == item),], variance=TRUE)[2],5),sep = "±") 
  }
  if (loopingyear=="17"){
    HaploSummary[loopingsite,4] <- paste(round(hap.div(LoopDat[which(paste(location,year,sep="_") == item),],variance=TRUE)[1],3),
                                         round(hap.div(LoopDat[which(paste(location,year,sep="_") == item),], variance=TRUE)[2],3),sep = "±") 
    HaploSummary[loopingsite,6] <- paste(round(nuc.div(LoopDat[which(paste(location,year,sep="_") == item),],variance=TRUE)[1],5),
                                         round(nuc.div(LoopDat[which(paste(location,year,sep="_") == item),], variance=TRUE)[2],5),sep = "±")
  }
}

write.csv(HaploSummary,paste0("data/sanger/HaploStats/",critter,".stats.csv"))

#calculate pairwise Fst
LoopHierfstat <- genind2hierfstat(DNAbin2genind(LoopDat),pop=paste(location,year,sep="_"))
Fst <- genet.dist(LoopHierfstat,diploid=F,method="Nei87")
Fst[Fst<0] <- 0
Fst[Fst=="NaN"] <- 0
write.csv(as.matrix(Fst),paste0("data/sanger/Fst/",critter,"Fst.csv"))


}

#now lets plot haplotype / nuc diversity over time 

#we are going to dynamically grow a data which is bad scripting
HaploSumData  <- c()


input.files <- list.files("data/sanger/HaploStats/")

for (spp in input.files){
  critter <- gsub(".stats.csv","",spp)
  loopfile <- read.csv(paste0("data/sanger/HaploStats/",spp))
  loopfile <- loopfile[!is.na(loopfile$DivH.09) & !is.na(loopfile$DivH.17),]
  
  HaploSumData <- rbind(HaploSumData,cbind(c(as.numeric(unlist(lapply(strsplit(loopfile$DivH.09,"±"),"[[",1))),
          as.numeric(unlist(lapply(strsplit(loopfile$DivH.17,"±"),"[[",1)))),
        c(rep("S.09",length(loopfile$DivH.09)),rep("S.17",length(loopfile$DivH.17))),
        rep("DivH",length(loopfile$DivH.09)),
        rep(critter,length(loopfile$DivH.09))))
  HaploSumData <- rbind(HaploSumData,cbind(c(as.numeric(unlist(lapply(strsplit(loopfile$DivN.09,"±"),"[[",1))),
                             as.numeric(unlist(lapply(strsplit(loopfile$DivN.17,"±"),"[[",1)))),
                           c(rep("S.09",length(loopfile$DivN.09)),rep("S.17",length(loopfile$DivN.17))),
                           rep("DivN",length(loopfile$DivN.09)),
                           rep(critter,length(loopfile$DivN.09))))

  
}

#HaploSummary Charts and Stats
HaploSumData <- as.data.frame(HaploSumData)
names(HaploSumData) <- c("value","year","stat","critter")
HaploSumData$value <- as.numeric(HaploSumData$value)

sink("data/sanger/HaploStatT/t.test.txt")
pdf("figures/HaploDiv.pdf",heigh=5.5,width = 5.5)
par(mfrow=c(2,2),mar=c(4.1, 4.1, 2.1, 1.1))
count <- 1
letters <- c("a","b","c","d")
for (spp in list.files("data/sanger/HaploStats/")){
critter <- gsub(".stats.csv","",spp)
print(critter)

test <- HaploSumData[HaploSumData$stat=="DivH" & HaploSumData$critter==critter,]
t.out <- t.test(test$value[test$year=="S.09"],test$value[test$year=="S.17"],paired = T)
boxplot(test$value~test$year,border = "white",xaxt='n',xlab = "Survey Year",ylab="Haplotype Diversity")
axis(1,at=1:2,labels=c("2009","2017"))
points(as.factor(test$year),test$value,pch=16)
segments(rep(1,length(test$year[test$year=="S.09"])),
          test$value[test$year=="S.09"],
          rep(2,length(test$year[test$year=="S.09"])),
          test$value[test$year=="S.17"],col="red",lty=2)
mtext(letters[count],side=3,line=1,at=0,font=2,cex=1.2)
print(t.out)
count <- count+1
}
dev.off()
count <- 1
pdf("figures/NuclDiv.pdf",heigh=5.5,width = 5.5)
par(mfrow=c(2,2),mar=c(4.1, 4.1, 2.1, 1.1))
for (spp in list.files("data/sanger/HaploStats/")){
  critter <- gsub(".stats.csv","",spp)
  
  test <- HaploSumData[HaploSumData$stat=="DivN" & HaploSumData$critter==critter,]
  t.out <- t.test(test$value[test$year=="S.09"],test$value[test$year=="S.17"],paired = T)
  boxplot(test$value~test$year,border = "white",xaxt='n',xlab = "Survey Year",ylab="Nucleotide Diversity")
  axis(1,at=1:2,labels=c("2009","2017"))
  points(as.factor(test$year),test$value,pch=16)
  segments(rep(1,length(test$year[test$year=="S.09"])),
           test$value[test$year=="S.09"],
           rep(2,length(test$year[test$year=="S.09"])),
           test$value[test$year=="S.17"],col="red",lty=2)
  mtext(letters[count],side=3,line=1,at=0,font=2,cex=1.2)
  print(t.out)
  count <- count+1
}
dev.off()
sink()


#2.2 Haplotype Maps

#Let's plot a haplotype map (with POPART)


#We'll count haplotypes in the next loop so lets make a landing pad for the output data
files <- c("Ciona.fa","Clav.fa","Styela.fa","Micro.fa")

HaploCounts <- as.data.frame(matrix(0,ncol = 3,nrow = length(files)))
colnames(HaploCounts) <- c("Shared","Unique17","Unique09")
rownames(HaploCounts) <- gsub(".fa","",files)

for (loopfile in files){

critter <- gsub(".fa","",loopfile)
file <- paste0("data/sanger/Alignments/",loopfile)

#Read in a FASTA and detect haplotypes
h<-haplotypes::haplotype(read.fas(file),indels="s")

#get metadata (works only if you format your fasta with underscores)
location <- sapply(strsplit(names(read.fas(file)),"_"),"[[",2)
year <- sapply(strsplit(names(read.fas(file)),"_"),"[[",1)

#pull in seqs
seqs <- read.fasta(file)

#outputHaplotypes
haplo <- seqs[unname(h@uniquehapind)]
names(haplo) <- paste0("Haplo_",1:length(haplo))
write.fasta(haplo,names(haplo),paste0("data/sanger/Alignments/",basename(file),"_haplo.fasta")) 

#now we output the nexus file 
write.nexus.data(haplo,paste0("data/sanger/Alignments/NexusFiles/",basename(file),"_nexus.nex"))

#and a trait file for colours
traitMatrix <- matrix(0,nrow = length(haplo),ncol = length(haplo))
colnames(traitMatrix) <- names(haplo)
rownames(traitMatrix) <- names(haplo)
for (num in 1:length(haplo)){
  traitMatrix[num,num] <- 10
}
write.csv(traitMatrix,paste0("data/sanger/Alignments/NexusFiles/",basename(file),"_table.csv"),quote = FALSE)

##Now lets create a table with all our haplotypes in 

haploTab <- matrix(0,nrow=length(haplo),ncol=length(levels(RA.09$variable)))
colnames(haploTab) <- levels(RA.09$variable)
rownames(haploTab) <- names(haplo)

#09 data

haploTab09 <- haploTab
for (hap in 1:length(h@hapind)){
  looplocation <- sapply(strsplit(names(unlist(h@hapind[hap])),"_"),"[[",2)
  loopyear <- sapply(strsplit(sapply(strsplit(names(unlist(h@hapind[hap])),"_"),"[[",1),"\\."),"[[",2)
  looplocation <- looplocation[loopyear=="09"]
  for (item in names(table(looplocation))){
    haploTab09[hap,item] <- unname(table(looplocation)[item])
  }
}
#17 data
haploTab17 <- haploTab
for (hap in 1:length(h@hapind)){
  looplocation <- sapply(strsplit(names(unlist(h@hapind[hap])),"_"),"[[",2)
  loopyear <- sapply(strsplit(sapply(strsplit(names(unlist(h@hapind[hap])),"_"),"[[",1),"\\."),"[[",2)
  looplocation <- looplocation[loopyear=="17"]
  for (item in names(table(looplocation))){
    haploTab17[hap,item] <- unname(table(looplocation)[item])
  }
}
#Save Data
write.csv(cbind(as.data.frame(haploTab09),"Seqs"=unlist(getSequence(haplo,as.string = T))),paste0("data/sanger/HaploDistribution/",critter,"09.csv"))
write.csv(cbind(as.data.frame(haploTab17),"Seqs"=unlist(getSequence(haplo,as.string = T))),paste0("data/sanger/HaploDistribution/",critter,"17.csv"))

#Count number of unique haplotypes per year 

YearComp <- as.data.frame(cbind("Y09"=rowSums(haploTab09),"Y17"=rowSums(haploTab17)))
YearComp$result <- ifelse(YearComp$Y09>0&YearComp$Y17>0,"Shared",
       ifelse(YearComp$Y09==0&YearComp$Y17>0,"Unique17",
              ifelse(YearComp$Y09>0&YearComp$Y17==0,"Unique09","None")))

HaploCounts[critter,names(table(YearComp$result))] <- table(YearComp$result)


#get many colours in the palette
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
palette(col_vector[20:(20+length(h@hapind))])


##Here we plot the haplotypes over time 


pdf(paste0("figures/figure2/",critter,"HapDist.pdf"),height = 5.5,width = 9)
par(mfrow=c(2,1),mar=c(3.1, 3.5, 1.1, 1.1))
barplot(prop.table(haploTab09,2),axisnames=F,col=1:length(h@hapind),cex.axis=2,las=1)
barplot(prop.table(haploTab17,2),col=1:length(h@hapind),cex.names=1.5,cex.axis=2,las=1)
dev.off()

pdf(paste0("data/sanger/HaploTypeColours/",critter,"_cols.pdf"),height = 5.5,width = 9)
par(mfrow=c(1,1),mar=c(3.1, 2.1, 1.1, 1.1))
barplot(rep(10,length(h@hapind)),col=1:length(h@hapind),names.arg = 1:length(h@hapind))
dev.off()

}


##AMOVA
files <- c("Ciona.fa","Clav.fa","Styela.fa","Micro.fa")
sink("data/sanger/AMOVAoutput.txt")
for (spp in files){ 

seqs <- read.FASTA(paste0("data/sanger/Alignments/",spp))

looplocation <- as.factor(sapply(strsplit(names(seqs),"_"),"[[",2))
loopyear <- as.factor(sapply(strsplit(names(seqs),"_"),"[[",1))

seqs2 <- DNAbin2genind(seqs)

strata(seqs2) <- data.frame(cbind(loopyear,looplocation))
nameStrata(seqs2) <- ~Loopyear/Looplocation


#Seq.dist <- dist.dna(seqs)
looplocation <- as.factor(sapply(strsplit(names(seqs),"_"),"[[",2))
loopyear <- as.factor(sapply(strsplit(names(seqs),"_"),"[[",1))

#test AMOVA
AMOVA <- poppr.amova(seqs2,~Loopyear/Looplocation)

amova.result <- randtest(AMOVA,nrepet = 1000)
plot(amova.result)
print(spp)
print(AMOVA)
print(amova.result)

}

sink()
####====3.0 eDNA Genotypes====######

#Now lets analyse the eDNA haplotypes to see if they reflect the collected data 

eDNA.haplo.dada <- read.csv("data/eDNA/cleaned/COI.metaphy.csv",row.names = 1)
eDNA.haplo.unoise <- read.csv("data/eDNA/cleaned/COI.metaphy.u.csv",row.names = 1)

files <- c("Ciona.fasta","Clavelina.fasta","Micro.fasta","Styela.fasta")
for (loopfile in files){
  
  
  critter <- gsub(".fasta","",loopfile)
  file <- paste0("data/eDNA/HaplotypeAlignments/",loopfile)
  
  #Read in a FASTA and detect haplotypes
  h<-haplotypes::haplotype(read.fas(file),indels="s")
  
  #get metadata (works only if you format your fasta with underscores)
  #location <- sapply(strsplit(names(read.fas(file)),"_"),"[[",2)
  #year <- sapply(strsplit(names(read.fas(file)),"_"),"[[",1)
  
  #pull in seqs
  seqs <- read.fasta(file)
  
  #outputHaplotypes
  haplo <- seqs[unname(h@uniquehapind)]
  names(haplo) <- paste0("Haplo_",1:length(haplo))
  write.fasta(haplo,names(haplo),paste0("data/eDNA/HaplotypeAlignments/HaplotypeFasta/",basename(file),"_haplo.fasta")) 
  
  #now we output the nexus file 
  write.nexus.data(haplo,paste0("data/eDNA/HaplotypeAlignments/Nexus/",basename(file),"_nexus.nex"))
  
  #and a trait file for colours
  traitMatrix <- matrix(0,nrow = length(haplo),ncol = length(haplo))
  colnames(traitMatrix) <- names(haplo)
  rownames(traitMatrix) <- names(haplo)
  for (num in 1:length(haplo)){
    traitMatrix[num,num] <- 10
  }
  write.csv(traitMatrix,paste0("data/eDNA/HaplotypeAlignments/Nexus/",basename(file),"_table.csv"),quote = FALSE)
  
  ##Now lets create a table with all our haplotypes in 
  
  haploTab <- matrix(0,nrow=length(haplo),ncol=length(levels(RA.09$variable)))
  colnames(haploTab) <- levels(RA.09$variable)
  rownames(haploTab) <- names(haplo)
  
  #17 data
  haploTab17 <- haploTab
  haplodata <- names(unlist(h@hapind))
  loopyear <- sapply(strsplit(sapply(strsplit(names(unlist(h@hapind)),"_"),"[[",1),"\\."),"[[",2)
  haplodata <- haplodata[loopyear=="17"]
  loopsamples <- sapply(strsplit(haplodata,"\\."),"[[",2)
  haploloop <- sapply(strsplit(haplodata,"\\."),"[[",1)
  rownames(haploTab17) <- paste0("haplotype",1:h@nhap)
  
  for (hap in unique(haploloop)){
    for (item in names(table(sapply(strsplit(loopsamples,"_"),"[[",2)))){
      haploTab17[hap,item] <- unname(table(sapply(strsplit(loopsamples[haploloop==hap],"_"),"[[",2))[item])
    }
  }
  
  haploTab17[is.na(haploTab17)] <- 0
  

  #eDNA data dada2
  haploTab.eDNA.dada <- haploTab
  
  #loop over haplotypes in big file 
  
  for (hap in 1:length(haplo)){
    loopseq <- unlist(getSequence(haplo[hap],as.string = T))
    OTUindex <- grep(loopseq,eDNA.haplo.dada$seqs)
    loopdata <- eDNA.haplo.dada[OTUindex,1:13]
    if(length(loopdata[,1])==0){next}else{
      haploTab.eDNA.dada[hap,] <-unlist(colSums(loopdata[colnames(haploTab)]))
  }}
  
  
  
  #eDNA data unoise3
  haploTab.eDNA.u3 <- haploTab
  
  #loop over haplotypes in big file 
  
  for (hap in 1:length(haplo)){
    loopseq <- unlist(getSequence(haplo[hap],as.string = T))
    OTUindex <- grep(loopseq,eDNA.haplo.unoise$seqs)
    loopdata <- eDNA.haplo.unoise[OTUindex,1:13]
    if(length(loopdata[,1])==0){next}else{
      haploTab.eDNA.u3[hap,] <-unlist(colSums(loopdata[colnames(haploTab)]))
    }}
  
  

  
  
  #Save Data
  write.csv(cbind(as.data.frame(haploTab17),"Seqs"=unlist(getSequence(haplo,as.string = T))),paste0("data/eDNA/HaplotypeAlignments/HaplotypeDistribution/",critter,"17.csv"))
  write.csv(cbind(as.data.frame(haploTab.eDNA.dada),"Seqs"=unlist(getSequence(haplo,as.string = T))),paste0("data/eDNA/HaplotypeAlignments/HaplotypeDistribution/",critter,".eDNA.dada.csv"))
  write.csv(cbind(as.data.frame(haploTab.eDNA.u3),"Seqs"=unlist(getSequence(haplo,as.string = T))),paste0("data/eDNA/HaplotypeAlignments/HaplotypeDistribution/",critter,".eDNA.u3.csv"))
  
  
  #get many colours in the palette
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  palette(col_vector[20:(20+length(h@hapind))])
  
  pdf(paste0("figures/eDNAcomparison/",critter,".eDNA.HapDist.pdf"),height = 5.5,width = 9)
  par(mfrow=c(3,1),mar=c(3.1, 4.1, 1.1, 1.1))
  barplot(prop.table(haploTab.eDNA.dada,2),axisnames=F,col=1:length(h@hapind),cex.axis=2,las=1)
  mtext("DADA2",4,cex=1.5)
  barplot(prop.table(haploTab.eDNA.u3,2),axisnames=F,col=1:length(h@hapind),cex.axis=2,las=1)
  mtext("UNOISE3",4,cex=1.5)
  barplot(prop.table(haploTab17,2),col=1:length(h@hapind),cex.axis=2, cex.names=2,las=1)
  mtext("Sanger",4,cex=1.5)
  dev.off()
  
  pdf(paste0("data/eDNA/HaplotypeAlignments/colours/",critter,"_cols.pdf"),height = 5.5,width = 9)
  par(mfrow=c(1,1),mar=c(3.1, 2.1, 1.1, 1.1))
  barplot(rep(10,length(h@hapind)),col=1:length(h@hapind),names.arg = 1:length(h@hapind))
  dev.off()
  
}  



