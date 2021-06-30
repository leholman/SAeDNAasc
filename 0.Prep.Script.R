#############################################
####== South African Ascidian Analysis ==####
####==== Luke E. Holman====01.10.2020====####
#############################################

###Script 0 - QC & Prep Script###

####====0.0 Packages====######

library(metabarTOAD)
library("reshape")
library("seqinr")


####====0.1 eDNA Data Prep====######


metadat <- read.csv("data/eDNA/metadata/locations.csv")
COI <- read.csv("data/eDNA/COI.dada2.tab.csv")
z18S <- read.csv("data/eDNA/18S.dada2.tab.csv") 
COI.unoise3 <- read.table("data/eDNA/COI.unoise3.a5.relaxed.tab.tsv",sep="\t",header=1,row.names = 1)

#Minimum number of reads
minreads <- 2


#Lets prep the taxonomy for the ascidian data

#COIassignments <- ParseTaxonomy(blastoutput = "data/eDNA/taxonomy/COI.results.txt",lineages = "data/eDNA/taxonomy/lineages-2019-08-14.csv.gz",
#                                pctThreshold = 97,
#                                lwrpctThreshold = 95,
#                               lwrcovpct = 90)
#write.csv(COIassignments,file="data/eDNA/taxonomy/COI.assignments.csv")

COIassignments <- read.csv("data/eDNA/taxonomy/COI.assignments.csv",row.names = 1,stringsAsFactors = F)


#z18Sassignments <- ParseTaxonomy(blastoutput = "data/eDNA/taxonomy/18S.results.txt",lineages = "data/eDNA/taxonomy/lineages-2019-08-14.csv.gz",
#                                pctThreshold = 99,
#                                lwrpctThreshold = 97,
#                                lwrcovpct = 90)
#write.csv(z18Sassignments,file="data/eDNA/taxonomy/18S.assignments.csv"

z18Sassignments <- read.csv("data/eDNA/taxonomy/18S.assignments.csv",row.names = 1,stringsAsFactors = F)

#we pull in the unoise3 assignments (there are so many ParseTaxonomy has to be done on HPC)

COI.u3.assignments <- read.csv("data/eDNA/taxonomy/COI.unoise3.results.csv",row.names = 1,stringsAsFactors = F)



##Now we can pull in the data

aCOI <- COI[match(na.omit(as.character(COIassignments$OTU[COIassignments$class=="Ascidiacea"])),row.names(COI)),]
aCOIassign <- COIassignments[match(na.omit(as.character(COIassignments$OTU[COIassignments$class=="Ascidiacea"])),as.character(COIassignments$OTU)),]

a18S <- z18S[match(na.omit(as.character(z18Sassignments$OTU[z18Sassignments$class=="Ascidiacea"])),row.names(z18S)),]
a18Sassign <- z18Sassignments[match(na.omit(as.character(z18Sassignments$OTU[z18Sassignments$class=="Ascidiacea"])),as.character(z18Sassignments$OTU)),]

aCOI.u3 <- COI.unoise3[match(na.omit(as.character(COI.u3.assignments$OTU[COI.u3.assignments$class=="Ascidiacea"])),row.names(COI.unoise3)),]
aCOI.u3.assign <- COI.u3.assignments[match(na.omit(as.character(COI.u3.assignments$OTU[COI.u3.assignments$class=="Ascidiacea"])),as.character(COI.u3.assignments$OTU)),]

#Filtering

rawdat <-aCOI


#Seperate controls and samples
samples <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="sample","COI"]]
controls <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="control","COI"]]

#Filter 1 - minimum number of reads for any ID
samples[samples< minreads] <- 0
samples <- samples[rowSums(samples) > 0,]

#Filter 2 - within samples OTU must appear in more than one sample (this works becuase there are lots of reps per site and sample)
filtersam <- samples
filtersam[filtersam>0 ] <- 1
filtersam <-filtersam[rowSums(filtersam) > 1,]
samples <- samples[rownames(samples) %in% rownames(filtersam),]

#Filter 3 -Make the maximum number of reads for each OTU in the contam the zero value in the main data
controlsCONTAM <- controls[rowSums(controls) > 0,]
for (contamOTU in 1:length(controlsCONTAM[,1])){
  loopOTU <- row.names(controlsCONTAM[contamOTU,])
  loopMax <- max(as.numeric(controlsCONTAM[contamOTU,]))
  if (any(is.na(samples[loopOTU,]))){next}
  samples[loopOTU,samples[loopOTU,]<loopMax] <- 0
  print(paste("Cleaning contaminants",contamOTU))
}


#get rid of natural sites
samples <- samples[colnames(samples) %in% metadat[metadat$sitetype=="m","COI"]]


#Let's think in relative proportion of reads
##This calculates the num of reads per sample 
readcount <-colSums(COI[match(colnames(samples),colnames(COI))])
##this divides all rows of the samples by the vector created above
samples <-sweep(samples, 2,readcount , "/")


###Add taxonomy and sequence
OTUsCOI <- read.fasta("data/eDNA/COI.dada2.OTUs.fasta")

seqs <- unlist(getSequence(OTUsCOI[match(rownames(samples),names(OTUsCOI))]
                                  ,as.string = TRUE))

test <- OTUsCOI[match(rownames(samples),names(OTUsCOI))]
write.fasta(test,names=names(test),file.out = "data/eDNA/taxonomy/BLASTAscidianCheck/COI_OTUS.fa")


#read in checked taxonomy
tax.checked.COI <- read.csv("data/eDNA/taxonomy/BLASTAscidianCheck/COI.check.csv")

#output data before being averaged
outputCOI <- cbind(samples,seqs,aCOIassign[match(rownames(samples),aCOIassign$OTU),])

##Save this data for metaphylogeogrphy 
write.csv(outputCOI,"data/eDNA/cleaned/COI.metaphy.notpooled.csv")

#Calculate average as below 
#Finally we collapse the technical replicates (by averaging them) into a matrix 

Csamples <- matrix(ncol=length(unique(substr(colnames(samples),3,4))),nrow = length(samples[,1]))
colnames(Csamples) <- unique(substr(colnames(samples),3,4))
rownames(Csamples) <- rownames(samples)
for (site in unique(substr(colnames(samples),3,4))){
  Csamples[,site] <- rowMeans(samples[,substr(colnames(samples),3,4) == site])
}
COI.asc.site <- as.data.frame(Csamples)
COI.asc.siteCollapsed <- COI.asc.site

#lets output the data for eDNA haplotypes
seqs.2 <- unlist(getSequence(OTUsCOI[match(rownames(COI.asc.siteCollapsed ),names(OTUsCOI))]
                           ,as.string = TRUE))

outputCOI <- cbind(COI.asc.siteCollapsed,seqs,"Asignment"=tax.checked.COI$CheckedTax[match(rownames(COI.asc.siteCollapsed),tax.checked.COI$OTU)])
##Save this data for metaphylogeogrphy 
write.csv(outputCOI,"data/eDNA/cleaned/COI.metaphy.csv")

##Collapse by taxonomy 

#Get all spp
species <- tax.checked.COI$CheckedTax


for (name in names(table(species)[table(species)>1])){
  collapseOTUs <- as.character(rownames(COI.asc.site[species==name,])) 
  MotherOTU <- names(sort(rowSums(COI.asc.site[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  COI.asc.siteCollapsed [MotherOTU,] <- COI.asc.siteCollapsed [MotherOTU,] + colSums(COI.asc.siteCollapsed [collapseOTUs,])
  COI.asc.siteCollapsed  <- COI.asc.siteCollapsed [-match(collapseOTUs,rownames(COI.asc.siteCollapsed )),]
}

seqs <- unlist(getSequence(OTUsCOI[match(rownames(COI.asc.siteCollapsed),names(OTUsCOI))]
                           ,as.string = TRUE))

outputCOIcoll <- cbind(COI.asc.siteCollapsed,seqs,species[match(rownames(COI.asc.siteCollapsed),tax.checked.COI$OTU)])
names(outputCOIcoll)[15] <- "Assignment"

## Save the data for incidence
write.csv(outputCOIcoll,"data/eDNA/cleaned/COI.incidence.csv")



#Save this data for abundance / incidence


#18S

#####RERUN ANALYSIS WITH 18S Data###
#Filtering

rawdat <-a18S


#Seperate controls and samples
samples <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="sample","ZHAN"]]
controls <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="control","ZHAN"]]

#Filter 1 - minimum number of reads for any ID
samples[samples< minreads] <- 0
samples <- samples[rowSums(samples) > 0,]

#Filter 2 - within samples OTU must appear in more than one sample (this works becuase there are lots of reps per site and sample)
filtersam <- samples
filtersam[filtersam>0 ] <- 1
filtersam <-filtersam[rowSums(filtersam) > 1,]
samples <- samples[rownames(samples) %in% rownames(filtersam),]

#Filter 3 -Make the maximum number of reads for each OTU in the contam the zero value in the main data
controlsCONTAM <- controls[rowSums(controls) > 0,]
for (contamOTU in 1:length(controlsCONTAM[,1])){
  loopOTU <- row.names(controlsCONTAM[contamOTU,])
  loopMax <- max(as.numeric(controlsCONTAM[contamOTU,]))
  if (any(is.na(samples[loopOTU,]))){next}
  samples[loopOTU,samples[loopOTU,]<loopMax] <- 0
  print(paste("Cleaning contaminants",contamOTU))
}


#get rid of natural sites
samples <- samples[colnames(samples) %in% metadat[metadat$sitetype=="m","ZHAN"]]


#Let's think in relative proportion of reads
##This calculates the num of reads per sample 
readcount <-colSums(z18S[match(colnames(samples),colnames(z18S))])
##this divides all rows of the samples by the vector created above
samples <-sweep(samples, 2,readcount , "/")


###Add taxonomy and sequence
OTUs18S <- read.fasta("data/eDNA/18S.dada2.OTUs.fasta")

seqs <- unlist(getSequence(OTUs18S[match(rownames(samples),names(OTUs18S))]
                           ,as.string = TRUE))
test <- OTUs18S[match(rownames(samples),names(OTUs18S))]
write.fasta(test,names=names(test),file.out = "data/eDNA/taxonomy/BLASTAscidianCheck/18S_OTUS.fa")

output18S <- cbind(samples,seqs,a18Sassign[match(rownames(samples),a18Sassign$OTU),])

##Save this data for metaphylogeogrphy 
write.csv(output18S,"data/eDNA/cleaned/18S.metaphy.csv")

#Pull in taxa check data
tax.checked.18S <- read.csv("data/eDNA/taxonomy/BLASTAscidianCheck/18S.check.csv")


#Calculate average as below 
#Finally we collapse the technical replicates (by averaging them) into a matrix 

Csamples <- matrix(ncol=length(unique(substr(colnames(samples),3,4))),nrow = length(samples[,1]))
colnames(Csamples) <- unique(substr(colnames(samples),3,4))
rownames(Csamples) <- rownames(samples)
for (site in unique(substr(colnames(samples),3,4))){
  Csamples[,site] <- rowMeans(samples[,substr(colnames(samples),3,4) == site])
}
z18S.asc.site <- as.data.frame(Csamples)
z18S.asc.siteCollapsed <- z18S.asc.site
##Collapse by taxonomy 

#NOTE - we mark ciona as ciona spp. as NA's mess with the below. 
species <- tax.checked.18S$CheckedTax
species[is.na(species)] <- "NoAssign"

for (name in names(table(species)[table(species)>1])){
  if (name=="NoAssign"){next}
  collapseOTUs <- as.character(rownames(z18S.asc.site[species==name,])) 
  MotherOTU <- names(sort(rowSums(z18S.asc.site[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  z18S.asc.siteCollapsed [MotherOTU,] <- z18S.asc.siteCollapsed [MotherOTU,] + colSums(z18S.asc.siteCollapsed [collapseOTUs,])
  z18S.asc.siteCollapsed  <- z18S.asc.siteCollapsed [-match(collapseOTUs,rownames(z18S.asc.siteCollapsed )),]
}

seqs <- unlist(getSequence(OTUs18S[match(rownames(z18S.asc.siteCollapsed),names(OTUs18S))]
                           ,as.string = TRUE))

output18Scoll <- cbind(z18S.asc.siteCollapsed,seqs,species[match(rownames(z18S.asc.siteCollapsed),tax.checked.18S$OTU)])
names(output18Scoll)[15] <- "Assignment"
output18Scoll <- output18Scoll[output18Scoll$Assignment!="NoAssign",]


## Save the data for incidence
write.csv(output18Scoll,"data/eDNA/cleaned/18S.incidence.csv")



###Finally UNOISE3 data


rawdat <-aCOI.u3


#Seperate controls and samples
samples <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="sample","COI"]]
controls <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="control","COI"]]

#Filter 1 - minimum number of reads for any ID
samples[samples< minreads] <- 0
samples <- samples[rowSums(samples) > 0,]

#Filter 2 - within samples OTU must appear in more than one sample (this works becuase there are lots of reps per site and sample)
filtersam <- samples
filtersam[filtersam>0 ] <- 1
filtersam <-filtersam[rowSums(filtersam) > 1,]
samples <- samples[rownames(samples) %in% rownames(filtersam),]

#Filter 3 -Make the maximum number of reads for each OTU in the contam the zero value in the main data
controlsCONTAM <- controls[rowSums(controls) > 0,]
for (contamOTU in 1:length(controlsCONTAM[,1])){
  loopOTU <- row.names(controlsCONTAM[contamOTU,])
  loopMax <- max(as.numeric(controlsCONTAM[contamOTU,]))
  if (any(is.na(samples[loopOTU,]))){next}
  samples[loopOTU,samples[loopOTU,]<loopMax] <- 0
  print(paste("Cleaning contaminants",contamOTU))
}


#get rid of natural sites
samples <- samples[colnames(samples) %in% metadat[metadat$sitetype=="m","COI"]]


#Let's think in relative proportion of reads
##This calculates the num of reads per sample 
readcount <-colSums(COI[match(colnames(samples),colnames(COI))])
##this divides all rows of the samples by the vector created above
samples <-sweep(samples, 2,readcount , "/")


###Add taxonomy and sequence
OTUsCOI <- read.fasta("data/eDNA/COI.unoise3.OTU.fasta")

seqs <- unlist(getSequence(OTUsCOI[match(rownames(samples),names(OTUsCOI))]
                           ,as.string = TRUE))

test <- OTUsCOI[match(rownames(samples),names(OTUsCOI))]
write.fasta(test,names=names(test),file.out = "data/eDNA/taxonomy/BLASTAscidianCheck/COI_U3_OTUS.fa")


#read in checked taxonomy
tax.checked.COI.u <- read.csv("data/eDNA/taxonomy/BLASTAscidianCheck/COI.U.check.csv")

#output data before being averaged
outputCOI.u <- cbind(samples,seqs,aCOI.u3.assign[match(rownames(samples),aCOI.u3.assign$OTU),])

##Save this data for metaphylogeogrphy 
write.csv(outputCOI,"data/eDNA/cleaned/COI.u.metaphy.notpooled.csv")

#Calculate average as below 
#Finally we collapse the technical replicates (by averaging them) into a matrix 

Csamples <- matrix(ncol=length(unique(substr(colnames(samples),3,4))),nrow = length(samples[,1]))
colnames(Csamples) <- unique(substr(colnames(samples),3,4))
rownames(Csamples) <- rownames(samples)
for (site in unique(substr(colnames(samples),3,4))){
  Csamples[,site] <- rowMeans(samples[,substr(colnames(samples),3,4) == site])
}
COI.asc.u.site <- as.data.frame(Csamples)
COI.asc.u.siteCollapsed <- COI.asc.u.site

#lets output the data for eDNA haplotypes
seqs.2 <- unlist(getSequence(OTUsCOI[match(rownames(COI.asc.u.siteCollapsed ),names(OTUsCOI))]
                             ,as.string = TRUE))

outputCOI.u <- cbind(COI.asc.u.siteCollapsed,seqs,"Assignment"=tax.checked.COI.u$CheckedTax[match(rownames(COI.asc.u.siteCollapsed),tax.checked.COI.u$OTU)])
##Save this data for metaphylogeogrphy 
write.csv(outputCOI.u,"data/eDNA/cleaned/COI.metaphy.u.csv")



##Collapse by taxonomy 

#Get all spp
species <- tax.checked.COI.u$CheckedTax


for (name in names(table(species)[table(species)>1])){
  collapseOTUs <- as.character(rownames(COI.asc.u.site[species==name,])) 
  MotherOTU <- names(sort(rowSums(COI.asc.u.site[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  COI.asc.u.siteCollapsed [MotherOTU,] <- COI.asc.u.siteCollapsed [MotherOTU,] + colSums(COI.asc.u.siteCollapsed [collapseOTUs,])
  COI.asc.u.siteCollapsed  <- COI.asc.u.siteCollapsed [-match(collapseOTUs,rownames(COI.asc.u.siteCollapsed )),]
}

seqs <- unlist(getSequence(OTUsCOI[match(rownames(COI.asc.u.siteCollapsed),names(OTUsCOI))]
                           ,as.string = TRUE))

outputCOIcoll <- cbind(COI.asc.u.siteCollapsed,seqs,species[match(rownames(COI.asc.u.siteCollapsed),tax.checked.COI.u$OTU)])
names(outputCOIcoll)[15] <- "Assignment"

## Save the data for incidence
write.csv(outputCOIcoll,"data/eDNA/cleaned/COI.u.incidence.csv")







 