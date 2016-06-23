## load libraries
library(Biobase)
library(RColorBrewer)
library(vegan)
library(WGCNA)
library(vegetarian)
library(gplots)


source("heatmap.microbe.R")

heatcol <- rainbow(256,start=0.2,end=1)
## SET up data
## Read in, remove unnecessary columns, etc.

source_level <- "species"
abund <- read.table(file.path("Data_Jan25_2013",paste("Adolescence",source_level,"abundance","txt",sep=".")),sep="\t",header=T,row.names=1,as.is=T)

da <- dim(abund)
abund <- abund[-da[1],]
abund <- abund[,-da[2]]
level <- abund$Level
abund$Level <- NULL


################################################## Read Counts
readData <- read.table(file.path("Data_Jan25_2013","Adolescence.readcounts.txt"),sep="\t",header=TRUE)

################################################## Parse names and generatate metaData object
getchars = function(str,n){ 
  split <- strsplit(str,split="") 
  substr <- lapply(split,"[",n)
  sapply(substr,paste,collapse="")	
}

Sample_ID <- sub("^X","",colnames(abund))

State <- sapply(getchars(Sample_ID,1),switch,"1"="girl","2"="mom",NA)
Sample <- sapply(getchars(Sample_ID,c(6,7,8)),switch,"Vag"="Vagina","Vul"="Vulva",NA)
Subject <- as.numeric(getchars(Sample_ID,c(1,2,3)))
Visit <- as.numeric(getchars(Sample_ID,c(4,5)))
Replicate<- as.numeric(getchars(Sample_ID,c(9)))

metaData <- data.frame(Sample_ID,State,Sample,Subject,Visit,Replicate)
metaData$fail[match(levels(factor(readData$Sample_ID)),metaData$Sample_ID)] <- tapply(readData$Fail,readData$Sample_ID,sum)
metaData$pass[match(levels(factor(readData$Sample_ID)),metaData$Sample_ID)] <- tapply(readData$Pass,readData$Sample_ID,sum)

#### check
abund["Lactobacillus gasseri;Lactobacillus johnsonii","X10201Vag1"]/sum(abund[,"X10201Vag1"]) ## sould equal 0.1079849
################################################## make and check sample order
ord <- order(metaData$Subject,metaData$State,metaData$Sample,metaData$Visit,metaData$Replicate)
metaData <- metaData[ord,]
abund <- abund[,ord]
colnames(abund) <- metaData$Sample_ID

#### check
abund["Lactobacillus gasseri;Lactobacillus johnsonii","10201Vag1"]/sum(abund[,"10201Vag1"]) ## should equal 0.1079849

rem <- which(rowSums(abund)==0)
if (length(rem) > 0){
  abund <- abund[-rem,]
  level <- level[-rem]
}

## SAVE ORINGINAL DATA
abund.orig <- abund
level.orig <- level
metaData.orig <- metaData
### filter Samples (NOT YET)
remSample <- colSums(abund) >= 1000 & !is.na(colSums(abund))
table(remSample)
#remSample
#FALSE  TRUE 
#    7   451 
#abund <- abund[,remSample]
#metaData <- metaData[remSample,]

## Histogram of Read Counts per Sample
pdf("ReadCountsHist.pdf")
hist(colSums(abund),breaks=50,xlab="Number of Reads",main="Histgram of Sample Read Counts")
dev.off()

### combine species level data
abund.red <- abund
level.red <- level
freqAbund.red <- sweep(abund.red,2,colSums(abund.red),"/")

## check
freqAbund.red["Lactobacillus gasseri;Lactobacillus johnsonii","10201Vag1"] ## should equal 0.1079849

id.nm <- rownames(abund.red)
id.nm[level=="species"] <- sapply(strsplit(rownames(abund.red)[level=="species"],split=" "),"[[",1L)
species2red <- unique(id.nm[level=="species"])

for(tax in species2red){
  ind <- id.nm %in% tax
  keepRows <- ((apply(freqAbund.red[ind,] >= 0.01,1,sum,na.rm=TRUE) > 1) | (apply(freqAbund.red[ind,] >= 0.05, 1, sum,na.rm=TRUE) > 0)) ## remove taxa with < 0.05 and only found in 1
  # Combine non species
  Other <- colSums(abund.red[ind,][!keepRows,])

  abund.tmp <- abund.red[ind,][keepRows,]
  level.tmp <- level.red[ind][keepRows]
  nms.tmp <- c(rownames(abund.tmp),paste(tax,"Other",sep=" "))
  abund.tmp <- rbind(abund.tmp,Other)
  level.tmp <- c(level.tmp,"species")
  rownames(abund.tmp) <- nms.tmp
  
  abund.red <- abund.red[!ind,]
  abund.red <- rbind(abund.red,abund.tmp)
  level.red <- level.red[!ind]
  level.red <- c(level.red,level.tmp)  

  freqAbund.red <- sweep(abund.red,2,colSums(abund.red),"/")
  id.nm <- rownames(abund.red)
  id.nm[level.red=="species"] <- sapply(strsplit(rownames(abund.red)[level.red=="species"],split=" "),"[[",1L)  
}

## check
freqAbund.red["Lactobacillus gasseri;Lactobacillus johnsonii","10201Vag1"] ## should equal 0.1079849

ord <- order(rownames(abund.red))
abund.red <- abund.red[ord,]
level.red <- level.red[ord]
freqAbund.red <- freqAbund.red[ord,]

## filter out genus for table
keepRows <- ((apply(freqAbund.red >= 0.01,1,sum,na.rm=TRUE) > 1) | (apply(freqAbund.red >= 0.05, 1, sum,na.rm=TRUE) > 0)) | level.red == "species"## remove taxa with < 0.05 and only found in 1

keepRows["Bacteria"] <- FALSE
# Combine non species
Other <- colSums(abund.red[!keepRows,])
abund.red <- abund.red[keepRows,]
level.red <- level.red[keepRows]

abund.red <- rbind(abund.red,"Other Bacteria"=Other)
level.red <- c(level.red,"Bacteria")
freqAbund.red <- sweep(abund.red,2,colSums(abund.red),"/")
## check
freqAbund.red["Lactobacillus gasseri;Lactobacillus johnsonii","10201Vag1"] ## should equal 0.1079849

#########################################################################################################################
#### WRITE OUT NEW DATA
### compute frequencies and write out
metaData$read_count <- colSums(abund.red)
freqAbund.red <- sweep(abund.red,2,colSums(abund.red),"/")
write.table(data.frame(Names=rownames(freqAbund.red),Level=level.red,freqAbund.red),"adolescence.analysis.relative.abundance",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(data.frame(Names=rownames(abund.red),Level=level.red,abund.red),"adolescence.analysis.species.abundance",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(metaData,file="adolescence.anno",sep="\t",row.names=FALSE,col.names=TRUE)
#########################################################################################################################

### GRAB GARDNERELLA
Gard <- grep("Gardnerella",rownames(abund.red))
abund.gard <- abund.red[,colSums(abund.red[Gard,])/colSums(abund.red) > 0.3]
freqAbund.gard <- sweep(abund.gard,2,colSums(abund.gard),"/")

write.table(metaData[colSums(freqAbund.red[Gard,]) > 0.3,],file="adolescence.Gardnerella.metadata",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(data.frame(Names=rownames(abund.gard),Level=level.red,abund.gard),"adolescence.Gardnerella.abundance",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(data.frame(Names=rownames(freqAbund.gard),Level=level.red,freqAbund.gard),"adolescence.Gardnerella.relative.abundance",sep="\t",row.names=FALSE,col.names=TRUE)


##############


girlMom <- c("pink","red")[as.numeric(as.factor(metaData$State))]
vagVul <- c("salmon","orange")[as.numeric(as.factor(metaData$Sample))]
others <- labels2colors(apply(metaData[,c("Subject","Visit")],2,function(x)as.numeric(as.factor(x))))
labelColors <- cbind(girlMom,vagVul,others)
colnames(labelColors) <-  c("mom/girl","vagina/vulva","subject","visit")

### done preparing data
##################################################################################################


##################################################################################################
### Paper 1 based Figures visit 1 to 3 only
##################################################################################################
require(vegan)
library(plotrix)
library(gplots)
mdata <- metaData
splitAbund <- abund.red
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors
subjectstring <- paste(mdata$Subject,mdata$Visit,substr(mdata$Sample,1,3),sep="-")


abund.standardized <- decostand(t(splitAbund), method = "hellinger")
abdist <- dist(abund.standardized, method="euc")

maintext <- "Cluster Analysis, all samples"

hcl <- hclust(abdist,method="ward")

pdf(file="cluster-analysis-all.pdf",width=18,height=5.5,pointsize=6)
plotDendroAndColors(hcl, 
                    colors = lColors,groupLabels=names(labelColors), rowText=mdata[,c("State","Sample","Subject","Visit")],
                    dendroLabels = subjectstring, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,cex.dendroLabels=0.4,
                    main = paste("Hellinger standardization with euclidean distances and average linkage clustering","[",maintext,"]"))
dev.off()

## Girls Only
### first reduce to vagina only
mdata <- metaData[metaData$State == "girl",]
splitAbund <- abund.red[,metaData$State == "girl"]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[metaData$State == "girl",]
subjectstring <- paste(mdata$Subject,mdata$Visit,substr(mdata$Sample,1,3),sep="-")


abund.standardized <- decostand(t(splitAbund), method = "hellinger")
abdist <- dist(abund.standardized, method="euc")

maintext <- "Cluster Analysis, girls only"

hcl <- hclust(abdist,method="ward")

pdf(file="cluster-analysis-girls.pdf",width=18,height=5.5,pointsize=6)
plotDendroAndColors(hcl, 
                   colors = lColors[,2:4],groupLabels=names(labelColors)[2:4], rowText=mdata[,c("Sample","Subject","Visit")],
                    dendroLabels = subjectstring, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,cex.dendroLabels=0.4,
                    main = paste("Hellinger standardization with euclidean distances and average linkage clustering","[",maintext,"]"))
dev.off()


mdata <- metaData[metaData$State == "girl" & metaData$Sample=="Vagina",]
splitAbund <- abund.red[,metaData$State == "girl" & metaData$Sample=="Vagina"]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[metaData$State == "girl" & metaData$Sample=="Vagina",]
subjectstring <- paste(mdata$Subject,mdata$Visit,substr(mdata$Sample,1,3),sep="-")


abund.standardized <- decostand(t(splitAbund), method = "hellinger")
abdist <- dist(abund.standardized, method="euc")

maintext <- "Cluster Analysis, girls (vagina) only"

hcl <- hclust(abdist,method="ward")

pdf(file="cluster-analysis-girls-vagina.pdf",width=18,height=5.5,pointsize=6)
plotDendroAndColors(hcl, 
                    colors = lColors[,3:4],groupLabels=names(labelColors)[3:4], rowText=mdata[,c("Subject","Visit")],
                    dendroLabels = subjectstring, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,cex.dendroLabels=0.4,
                    main = paste("Hellinger standardization with euclidean distances and average linkage clustering","[",maintext,"]"))
dev.off()

mdata <- metaData[metaData$State == "girl" & metaData$Sample=="Vulva",]
splitAbund <- abund.red[,metaData$State == "girl" & metaData$Sample=="Vulva"]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[metaData$State == "girl" & metaData$Sample=="Vulva",]
subjectstring <- paste(mdata$Subject,mdata$Visit,substr(mdata$Sample,1,3),sep="-")


abund.standardized <- decostand(t(splitAbund), method = "hellinger")
abdist <- dist(abund.standardized, method="euc")

maintext <- "girls (vulva) only"

hcl <- hclust(abdist,method="ward")

pdf(file="cluster-analysis-girls-vulva.pdf",width=18,height=5.5,pointsize=6)
plotDendroAndColors(hcl, 
                    colors = lColors[,3:4],groupLabels=names(labelColors)[3:4], rowText=mdata[,c("Subject","Visit")],
                    dendroLabels = subjectstring, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,cex.dendroLabels=0.4,
                    main = paste("Hellinger standardization with euclidean distances and average linkage clustering","[",maintext,"]"))
dev.off()

renamed <- rownames(abund.red)
renamed[c(50,91,92)] <- c(
"Lactobacillus crispatus",
"Streptococcus anginosus;etc.",
"Streptococcus australis;etc."
)
########### Per Sample Plot
load("phColorTbl_ct.2k.RData")
names(phColorTbl) <- sub("^L_","Lactobacillus ",names(phColorTbl))
names(phColorTbl) <- sub("_[0-9]+$","",names(phColorTbl))
phColorTbl <- phColorTbl[!duplicated(names(phColorTbl))]
my_colors <- unique(sub("[0-9]+$","",colors()))
my_colors <- setdiff(my_colors,unique(sub("[0-9]+$","",phColorTbl)))

cmatch <- unlist(sapply(names(phColorTbl),grep,rownames(abund.red)),use.names=FALSE)
taxa_colors <- my_colors[sample(1:length(my_colors),length(rownames(abund.red)))]
taxa_colors[cmatch] <- phColorTbl[match( names(sapply(names(phColorTbl),grep,rownames(abund.red))),names(phColorTbl))]

taxa_colors[is.na(taxa_colors)] <- setdiff(my_colors,taxa_colors)[sample(1:length(setdiff(my_colors,taxa_colors)),sum(is.na(taxa_colors)))]
names(taxa_colors) <- rownames(abund.red)

pdf(file=paste("sampleplots.new.pdf"),width=10.5,height=8,pointsize=12)
for (i in sort(unique(substr(metaData$Subject,2,3)))){
  #i = "01"
  pick <- which(substr(metaData$Subject,2,3) == i)
  if (length(pick) > 1){
  sample <- abund.red[,pick]
  freq <- sweep(sample,2,colSums(sample),"/")
  meta <- metaData[pick,]
  ord <- order(meta$Sample,meta$Visit,meta$State)
  freq <- freq[,ord]
  meta <- meta[ord,]
  
  layout(matrix(c(1,1,1,1,3,3,
                  1,1,1,1,3,3,
                  1,1,1,1,3,3,
                  1,1,1,1,3,3,
                  2,2,2,2,3,3,
                  2,2,2,2,3,3,
                  2,2,2,2,3,3,
                  2,2,2,2,3,3,
                  0,0,0,0,0,0),ncol=6,byrow=T))
  
  par(mar=c(4, 4, 6, 2) + 0.1)
  Top20 <- order(rowMeans(as.matrix(sample[,meta$Sample=="Vagina"])),decreasing=T)[1:10]
  Top20 <- order(rowSums(freq > 0.05),decreasing=T)[1:30]
  fullTaxa <- Top20
  stackpoly(t(freq[Top20,meta$Sample=="Vagina"]),stack=T,col=taxa_colors[Top20], xat=NA,xaxlab=rep("",length(meta$Visit[meta$Sample=="Vagina"])),ylim=c(0,1),main="Vagina",ylab="Proportion" )
  axis(paste("V",meta$Visit[meta$Sample=="Vagina"],sep=""),at=1:length(meta$Visit[meta$Sample=="Vagina"]),side=1,line=0,tick=F,las=3)  
  axis(meta$State[meta$Sample=="Vagina"],at=1:length(meta$State[meta$Sample=="Vagina"]),side=1,line=2,tick=F,las=3)

#  Top20 <- order(rowMeans(as.matrix(sample[,meta$Sample=="Vulva"])),decreasing=T)[1:10]
#  fullTaxa <- union(fullTaxa,Top20)
  stackpoly(t(freq[Top20,meta$Sample=="Vulva"]),stack=T,col=taxa_colors[Top20], xat=NA,xaxlab=rep("",length(meta$Visit[meta$Sample=="Vulva"])),ylim=c(0,1),main="Vulva",ylab="Proportion" )
  axis(paste("V",meta$Visit[meta$Sample=="Vulva"],sep=""),at=1:length(meta$Visit[meta$Sample=="Vulva"]),side=1,line=0,tick=F,las=3)  
  axis(meta$State[meta$Sample=="Vulva"],at=1:length(meta$State[meta$Sample=="Vulva"]),side=1,line=2,tick=F,las=3)
  
#  mtext("Visit",side=1,line=4,cex=0.8)
  mtext(paste("Subject",i),side=1,line=8,cex=3.0)
  par(mar = c(5,2,4,0))
  
  ord <- order(rownames(freq)[fullTaxa])
  plot(seq(1,10), seq(1,10), type = "n", axes = F,xlab="",ylab="")
  legend("bottomleft", renamed[fullTaxa][ord],bty="n",bg="#ffffff55",fill=taxa_colors[fullTaxa][ord],inset=0)
  
  #  QFM <- meta[1,c(2,,8,9,10,11)]
  #  textplot("Subject 2",show.rownames=FALSE,show.colnames = FALSE,valign="top")
  }  
}
dev.off()


freqa <- sweep(abund.red,2,colSums(abund.red),"/")
pdf("TaxaRepresentation.pdf",width=7,height=10,pointsize=8)
par(mar=c(3,15,4,2))
for (i in unique(metaData$Sample)){
  redmat <- freqa[,metaData$Sample == i]
  redmeta <- metaData[metaData$Sample == i,]
  maintext <- paste("Taxa represention across sample (",redmeta[1,"Sample"],")")
  boxplot(t(redmat[order(rowSums(redmat),decreasing=F),]),names=renamed[order(rowSums(redmat),decreasing=F)],horizontal=T,las=1,cex=0.3,cex.axis=0.75,
          main=maintext)
}

dev.off()

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################


### first reduce to vagina only
mdata <- metaData[metaData$Sample == "vagina",]
splitAbund <- abund[,metaData$Sample == "vagina"]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[metaData$Sample == "vagina",]

subjectstring <- paste(mdata$Subject,mdata$Visit,mdata$State,sep="-")
abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)

freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
rsum <- rowSums(freqAbund)
pick <- order(rsum,decreasing=TRUE)[1:25]
fcol <- lColors[,"mom/girl"]

pdf(file="MicrobeHeatTop25-DynamicV1-V3.pdf",bg="white",width=9,height=6.5,pointsize=10)

#tiff(file="MicrobeHeatTop25-DynamicV1-V3.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=10)
par(oma=c(2,2,2,10))
heatmap.microbe(as.matrix(freqAbund)[pick,],
                col=heatcol,
                distfun=vegdist,
                labCol=subjectstring,
                ColSideColors=fcol,
                ColSideLabels="Girl/Mother",
                colsep=which(!duplicated(mdata$Subject)[-1]),
                cexRow=1.3,
                cexCol=0.8,
                mar=c(7,5),
                dendrogram="none",
                Rowv=F,
                Colv=F,
                sepcolor="white",
                keysize=0.5,
                trace="none",
                key=T,
                density.info="none",family="sans")
#                main="Heatmap including only top 25 organisms (V1 to V3)")
dev.off()


### first reduce to vagina only
mdata <- metaData[metaData$Sample == "vagina" & metaData$State=="girl",]
splitAbund <- abund[,metaData$Sample == "vagina" & metaData$State=="girl"]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[metaData$Sample == "vagina" & metaData$State=="girl",]

subjectstring <- paste(mdata$Subject,mdata$Visit,sep="-")
abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)

freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
rsum <- rowSums(freqAbund)
pick <- order(rsum,decreasing=TRUE)[1:25]
#fcol <- lColors[,"mom/girl"]

pdf(file="VegenHeatTop25-DynamicV1-V3.pdf",bg="white",width=9,height=6.5,pointsize=8)
#tiff(file="VegenHeatTop25-DynamicV1-V3.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=8)
par(oma=c(2,2,5,10))
heatmap.2.2(as.matrix(freqAbund)[pick,],
col=heatcol,
distfun=vegdist,
labCol=subjectstring,
#ColSideColors=fcol,
#          ColSideLabels="Girl/Mom",
cexRow=1.3,
cexCol=1.3,

mar=c(3,5),
dendrogram="column",
Rowv=F,
Colv=T,
keysize=0.5,
trace="none",
key=F,
          density.info="none")#,
#          main="Heatmap including only top 25 organisms (V1 to V3)")
dev.off()
















subjectstring <- sub("V","",paste(mdata$V12,sep="."))

png(file="VegenHeatTop25-DynamicV1-V3.png",bg="white",width=6.9,height=5,units="in",res=600,pointsize=8)

par(oma=c(0,1,2,10))
heatmap.microbe(as.matrix(freqAbund)[order(rsum,decreasing=TRUE)[1:25],order(mdata$V11,mdata$V12)],col=col,
                distfun=vegdist,
                #labCol=NA,
                labCol=subjectstring[order(mdata$V11,mdata$V9,mdata$V12)],
                cexRow=1.4,cexCol=1.00, mar=c(5,5),
                dendrogram="none",Rowv=F,Colv=F,ColSideColors=fcol[order(mdata$V11,mdata$V9,mdata$V12)],
                ColSideNote=matrix(sub("V","",mdata$V12[order(mdata$V11,mdata$V12)]),nrow=1),
                ColSideLabels="Girl/Mom",
                colsep=(which(!duplicated(mdata$V11[ order(mdata$V11,mdata$V9,mdata$V12)]))-1)[-1],sepcolor="white",
                keysize=0.5,trace="none",key=T,density.info="none",colRot=TRUE)


dev.off()

png(file="VegenHeatTop25-DynamicV1-V3wCluster.png",bg="white",width=6.9,height=5,units="in",res=600,pointsize=8)

par(oma=c(5,1,2,10))
heatmap.2(as.matrix(freqAbund)[order(rsum,decreasing=TRUE)[1:25],],col=col,
          distfun=vegdist,
          #labCol=NA,
          labCol=paste(mdata$V11,mdata$V9,mdata$V10,mdata$V12),
          cexRow=1.4,cexCol=1.00,
          dendrogram="column",Rowv=F,Colv=T,ColSideColors=fcol,
          keysize=0.5,trace="none",key=T,density.info="none")


dev.off()


##################################################################################################
### Paper 1 based Figures visit 1 to 3 only
##################################################################################################
sm <- which(((metaData$V9 == "girl" & metaData$V12 %in% c("V1","V2","V3")) | (metaData$V9 == "mom" & metaData$V12=="V1")))

mdata <- metaData[sm,]
splitAbund <- abund[,sm]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[sm,]

ord <- order(mdata$V11,mdata$V9,mdata$V12,mdata$V10)
mdata <- mdata[ord,]
splitAbund <- splitAbund[,ord]
lColors <-  lColors[ord,]

## have moms and girls
#sm <-  which(mdata$V11 %in% names(tapply(mdata$V9,mdata$V11,function(x) length(table(x))))[tapply(mdata$V9,mdata$V11,function(x)length(table(x))) == 2])
## Reduce down to only those paired vagina/vulva
sm <- paste(mdata$V9,mdata$V11,mdata$V12) %in% names(table(paste(mdata$V9,mdata$V11,mdata$V12)))[table(paste(mdata$V9,mdata$V11,mdata$V12)) == 2]
#sm[which(mdata$V9 == "mom")] <- TRUE

mdata <- mdata[sm,]
splitAbund <- splitAbund[,sm]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- lColors[sm,]


subjectstring <- paste(mdata$V11,mdata$V12,mdata$V10,sep=".")
abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)

freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
rsum <- rowSums(freqAbund)
pick <- order(rsum,decreasing=TRUE)[1:25]


############## END CLUSTER

sm <- which(paste(mdata$V9,mdata$V11,mdata$V12) %in% names(table(paste(mdata$V9,mdata$V11,mdata$V12)))[table(paste(mdata$V9,mdata$V11,mdata$V12)) == 2])
#sm <- which(mdata$V11 %in% names(table(mdata$V11[mdata$V9 == "girl"]))[table(mdata$V11[mdata$V9 == "girl"]) > 1])

#sm <- which((mdata$V9 == "girl" & mdata$V12 %in% c("V1","V2","V3")))
#sm <- which(mdata$V11 %in% names(table(mdata$V11))[table(mdata$V11) > 3])
mdata <- mdata[sm,]
splitAbund <- splitAbund[,sm]



pdf(file="Cluster-DynamicV1-V1.pdf",bg="white",width=10,height=7,pointsize=8)

labelColors <- apply((mdata[,c(12,10)]),2,function(x)as.numeric(as.factor(x)))
labelColors <- cbind(labels2colors(labelColors[,1]),labels2colors(labelColors[,2]))

plotDendroAndColors(
  hclust(abdist), 
  colors=labelColors, 
  groupLabels = c("mom/girl","vagina/vulva"), 
  colorText=NULL,
  setLayout = TRUE, 
  autoColorHeight = TRUE, 
  colorHeight = 0.2, 
  dendroLabels = subjectstring, 
  addGuide = FALSE, guideAll = FALSE, 
  guideCount = 50, guideHang = 0.2, 
  cex.colorLabels = 0.8, cex.dendroLabels = 0.6, 
  marAll = c(1, 5, 3, 1), saveMar = TRUE, 
  abHeight = NULL, abCol = "red",main="Bray's clustering with metaData")
dev.off()


cmdfit <- cmdscale(abdist,eig=TRUE, k=2) # k is the number of dim

fcol <- c("pink","red")[as.numeric(as.factor(paste(mdata$V9)))]

# plot solution
pdf(file="mds-Dynamic.pdf",bg="white",width=10,height=7,pointsize=8)
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS, Bray's Distance", type="p",cex=2,pch=19,col=fcol)
#text(x, y, labels = subjectstring,col= fcol)
dev.off()

































# Generate A cluster (vegen distance) with metadata below based on all the available data
subjectstring <- paste(metaData$V11,metaData$V9,metaData$V10,metaData$V12,sep=".")
abdist <- vegdist(t(abund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)


pdf(file="Cluster-AllData.pdf",bg="white",width=10,height=7,pointsize=8)
plotDendroAndColors(
  hclust(abdist), 
  labelColors, 
  colorText=NULL,
  setLayout = TRUE, 
  autoColorHeight = TRUE, 
  colorHeight = 0.2, 
  dendroLabels = subjectstring, 
  addGuide = FALSE, guideAll = FALSE, 
  guideCount = 50, guideHang = 0.2, 
  cex.colorLabels = 0.8, cex.dendroLabels = 0.6, 
  marAll = c(1, 5, 3, 1), saveMar = TRUE, 
  abHeight = NULL, abCol = "red",main="Bray's clustering with metaData (All Data)")
dev.off()

## classical multidimentional scaling
cmdfit <- cmdscale(abdist,eig=TRUE, k=2) # k is the number of dim

fcol <- c("pink4","pink1","red4","red1")[as.numeric(as.factor(paste(metaData$V9,metaData$V10)))]

# plot solution
pdf(file="mds-AllData.pdf",bg="white",width=10,height=7,pointsize=4)
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS, Bray's Distance (All Data)", type="n")
text(x, y, labels = subjectstring,col= fcol)
dev.off()


### plot ordination, using meta MDS
ord <- metaMDS(t(abund))

pdf(file="metricmdsLARGE-AllData.pdf",bg="white",width=30,height=30,pointsize=8)
plot(ord,type="n")
points(ord,display="sites",cex=2,pch=21,col=fcol,bg=fcol)
text(ord,display="spec",cex=1,col="black")
dev.off()

pdf(file="metricmdsSmall-AllData.pdf",bg="white",width=10,height=7,pointsize=12)
plot(ord)
dev.off()

#### Results employing 'vegetarian' distance rather than vegan
vegetarianDist <- as.dist(1-sim.table(t(abund)))

pdf(file="VegetarianCluster-AllData.pdf",bg="white",width=10,height=7,pointsize=8)
plotDendroAndColors(
  hclust(vegetarianDist), 
  colors=labelColors, 
  colorText=NULL,
  setLayout = TRUE, 
  autoColorHeight = TRUE, 
  colorHeight = 0.2, 
  dendroLabels = subjectstring, 
  addGuide = FALSE, guideAll = FALSE, 
  guideCount = 50, guideHang = 0.2, 
  cex.colorLabels = 0.8, cex.dendroLabels = 0.6, 
  marAll = c(1, 5, 3, 1), saveMar = TRUE, 
  abHeight = NULL, abCol = "red",main="Vegetarian clustering with metaData (All Data)")
dev.off()


### compute frequencies and write out
freqAbund <- sweep(abund,2,colSums(abund),"/")
write.table(data.frame(Names=rownames(freqAbund),Level=level,freqAbund),"adolescence.species.relative.abundance",sep="\t",row.names=FALSE,col.names=TRUE)

## genrate full headmap of vegen results
pdf(file="VegenHeat-AllData.pdf",bg="white",width=7,height=10,pointsize=10)
heatmap.2(as.matrix(freqAbund),
          distfun=vegdist,
          labCol=metaData$ID,
          cexRow=0.4,
          ColSideColors=fcol,
          trace="none",
          col=cm.colors(100))
dev.off()

## only top 25 now
rsum <- rowSums(freqAbund)
pick <- order(rsum,decreasing=TRUE)[1:25]

pdf(file="VegenHeatTop25noRowD-AllData.pdf",bg="white",width=9,height=6.5,pointsize=10)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[pick,],
          col=heatcol,
          distfun=vegdist,
          labCol=subjectstring,
          cexRow=0.8,
          ColSideColors=fcol,
          dendrogram="column",
          Rowv=F,
          keysize=1.2,
          trace="none",
          key=T,
          density.info="none")
dev.off()

##################################################################################################
### END ALL
##################################################################################################



##################################################################################################
### Reduce to Mom and vagina only
##################################################################################################

mdata <- metaData[metaData$V10 == "vagina" & metaData$V9 == "mom",]
splitAbund <- abund[,metaData$V10 == "vagina"& metaData$V9 == "mom"]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[metaData$V10 == "vagina" & metaData$V9 == "mom",]

ord <- order(mdata$V11,mdata$V12)
mdata <- mdata[ord,]
splitAbund <- splitAbund[,ord]
lColors <-  lColors[ord,]

subjectstring <- paste(mdata$V11,mdata$V12,sep=".")
abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)

pdf(file="Cluster-mom-vagina.pdf",bg="white",width=10,height=7,pointsize=8)
plotDendroAndColors(
  hclust(abdist), 
  colors=lColors[, c("subject","visit","run")],
  colorText=NULL,
  setLayout = TRUE, 
  autoColorHeight = TRUE, 
  colorHeight = 0.2, 
  dendroLabels = subjectstring, 
  addGuide = FALSE, guideAll = FALSE, 
  guideCount = 50, guideHang = 0.2, 
  cex.colorLabels = 0.8, cex.dendroLabels = 0.6, 
  marAll = c(1, 5, 3, 1), saveMar = TRUE, 
  abHeight = NULL, abCol = "red",main="Bray's clustering (Mothers,Vaginal) with metaData")
dev.off()

freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
rsum <- rowSums(freqAbund)
pick <- order(rsum,decreasing=TRUE)[1:25]
fcol <- lColors[,"subject"]

pdf(file="VegenHeatTop25noRowD-mom-vagina.pdf",bg="white",width=9,height=6.5,pointsize=10)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[pick,],
          col=heatcol,
          distfun=vegdist,
          labCol=subjectstring,
          ColSideColors=fcol,
          cexRow=0.8,
          cexCol=0.8,
          mar=c(7,5),
          dendrogram="column",
          Rowv=F,
          Colv=T, 
          keysize=0.8,
          trace="none",
          key=T,
          density.info="none",
          main="Heatmap including only top 25 organisms (Mother,Vaginal)")
dev.off()


##################################################################################################
### Reduce to vagina only
##################################################################################################

mdata <- metaData[metaData$V10 == "vagina",]
splitAbund <- abund[,metaData$V10 == "vagina"]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[metaData$V10 == "vagina",]

ord <- order(mdata$V11,mdata$V12)
mdata <- mdata[ord,]
splitAbund <- splitAbund[,ord]
lColors <-  lColors[ord,]

subjectstring <- paste(mdata$V11,mdata$V12,sep=".")
abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)

freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
rsum <- rowSums(freqAbund)
pick <- order(rsum,decreasing=TRUE)[1:25]
fcol <- lColors[,"subject"]

pdf(file="VegenHeatTop25-vagina.pdf",bg="white",width=9,height=6.5,pointsize=10)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[pick,],
          col=heatcol,
          distfun=vegdist,
          labCol=subjectstring,
          ColSideColors=fcol,
          cexRow=0.8,
          cexCol=0.8,
          mar=c(7,5),
          dendrogram="column",
          Rowv=F,
          Colv=T, 
          keysize=0.8,
          trace="none",
          key=T,
          density.info="none",
          main="Heatmap including only top 25 organisms (Vaginal)")
dev.off()

tiff(file="VegenHeatTop25-vagina.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=8)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[pick,],
          col=heatcol,
          distfun=vegdist,
          labCol=subjectstring,
          ColSideColors=fcol,
          cexRow=0.8,
          cexCol=0.8,
          mar=c(7,5),
          dendrogram="column",
          Rowv=F,
          Colv=T, 
          keysize=0.8,
          trace="none",
          key=T,
          density.info="none",
          main="Heatmap including only top 25 organisms (Vaginal)")
dev.off()


##################################################################################################
### Paper 1 based Figures visit 1 to 3 only
##################################################################################################

### first reduce to vagina only
mdata <- metaData[metaData$V10 == "vagina",]
splitAbund <- abund[,metaData$V10 == "vagina"]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[metaData$V10 == "vagina",]

### then reduce to only those with girl V1 to V3 and Mom V1, remove c("S104","S117","S118","S121") due to small number of samples
sm <- which(((mdata$V9 == "girl" & mdata$V12 %in% c("V1","V2","V3")) | (mdata$V9 == "mom" & mdata$V12=="V1")) & !(mdata$V11 %in% c("S104","S117","S118","S121")))
mdata <- mdata[sm,]
splitAbund <- splitAbund[,sm]
lColors <- lColors[sm,]

ord <- order(mdata$V11,mdata$V9,mdata$V12)
mdata <- mdata[ord,]
splitAbund <- splitAbund[,ord]
lColors <-  lColors[ord,]

subjectstring <- paste(mdata$V11,mdata$V9,mdata$V12,sep=".")
abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)

freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
rsum <- rowSums(freqAbund)
pick <- order(rsum,decreasing=TRUE)[1:25]
fcol <- lColors[,"subject"]

tiff(file="VegenHeatTop25-DynamicV1-V3.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=8)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[pick,],
          col=heatcol,
          distfun=vegdist,
          labCol=subjectstring,
          ColSideColors=fcol,
          cexRow=0.8,
          cexCol=0.8,
          mar=c(7,5),
          dendrogram="column",
          Rowv=F,
          Colv=T, 
          keysize=0.8,
          trace="none",
          key=T,
          density.info="none",
          main="Heatmap including only top 25 organisms (V1 to V3)")
dev.off()

tiff(file="MicrobeHeatTop25-DynamicV1-V3.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=8)
par(oma=c(2,2,2,5))
heatmap.microbe(as.matrix(freqAbund)[pick,],
                col=heatcol,
                distfun=vegdist,
                labCol=subjectstring,
                ColSideColors=fcol,
                colsep=(which(!duplicated(mdata$V11)-1)[-1]),
                cexRow=0.8,
                cexCol=0.8,
                mar=c(7,5),
                dendrogram="none",
                Rowv=F,
                Colv=F,
                sepcolor="white",
                keysize=0.8,
                trace="none",
                key=T,
                density.info="none",
                main="Heatmap including only top 25 organisms (V1 to V3)")
dev.off()

subjectstring <- sub("V","",paste(mdata$V12,sep="."))

png(file="VegenHeatTop25-DynamicV1-V3.png",bg="white",width=6.9,height=5,units="in",res=600,pointsize=8)

par(oma=c(0,1,2,10))
heatmap.microbe(as.matrix(freqAbund)[order(rsum,decreasing=TRUE)[1:25],order(mdata$V11,mdata$V12)],col=col,
  distfun=vegdist,
  #labCol=NA,
  labCol=subjectstring[order(mdata$V11,mdata$V9,mdata$V12)],
  cexRow=1.4,cexCol=1.00, mar=c(5,5),
  dendrogram="none",Rowv=F,Colv=F,ColSideColors=fcol[order(mdata$V11,mdata$V9,mdata$V12)],
  ColSideNote=matrix(sub("V","",mdata$V12[order(mdata$V11,mdata$V12)]),nrow=1),
  ColSideLabels="Girl/Mom",
  colsep=(which(!duplicated(mdata$V11[ order(mdata$V11,mdata$V9,mdata$V12)]))-1)[-1],sepcolor="white",
  keysize=0.5,trace="none",key=T,density.info="none",colRot=TRUE)


dev.off()

png(file="VegenHeatTop25-DynamicV1-V3wCluster.png",bg="white",width=6.9,height=5,units="in",res=600,pointsize=8)

par(oma=c(5,1,2,10))
heatmap.2(as.matrix(freqAbund)[order(rsum,decreasing=TRUE)[1:25],],col=col,
  distfun=vegdist,
  #labCol=NA,
  labCol=paste(mdata$V11,mdata$V9,mdata$V10,mdata$V12),
  cexRow=1.4,cexCol=1.00,
  dendrogram="column",Rowv=F,Colv=T,ColSideColors=fcol,
  keysize=0.5,trace="none",key=T,density.info="none")


dev.off()



##################################################################################################
### Paper 1 based Figures visit 1 to 3 only
##################################################################################################
sm <- which(((metaData$V9 == "girl" & metaData$V12 %in% c("V1","V2","V3")) | (metaData$V9 == "mom" & metaData$V12=="V1")))

mdata <- metaData[sm,]
splitAbund <- abund[,sm]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- labelColors[sm,]

ord <- order(mdata$V11,mdata$V9,mdata$V12,mdata$V10)
mdata <- mdata[ord,]
splitAbund <- splitAbund[,ord]
lColors <-  lColors[ord,]

## have moms and girls
#sm <-  which(mdata$V11 %in% names(tapply(mdata$V9,mdata$V11,function(x) length(table(x))))[tapply(mdata$V9,mdata$V11,function(x)length(table(x))) == 2])
## Reduce down to only those paired vagina/vulva
sm <- paste(mdata$V9,mdata$V11,mdata$V12) %in% names(table(paste(mdata$V9,mdata$V11,mdata$V12)))[table(paste(mdata$V9,mdata$V11,mdata$V12)) == 2]
#sm[which(mdata$V9 == "mom")] <- TRUE

mdata <- mdata[sm,]
splitAbund <- splitAbund[,sm]
splitAbund <- splitAbund[rowSums(splitAbund) >0,]
lColors <- lColors[sm,]


subjectstring <- paste(mdata$V11,mdata$V12,mdata$V10,sep=".")
abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)

freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
rsum <- rowSums(freqAbund)
pick <- order(rsum,decreasing=TRUE)[1:25]

pdf(file="Cluster-DynamicV1-V3.pdf",bg="white",width=10,height=7,pointsize=8)
plotDendroAndColors(
  hclust(abdist,method="ward"), 
  colors=lColors[, c("vagina/vulva","subject","visit")],
  colorText=NULL,
  setLayout = TRUE, 
  autoColorHeight = TRUE, 
  colorHeight = 0.2, 
  dendroLabels = subjectstring, 
  addGuide = FALSE, guideAll = FALSE, 
  guideCount = 50, guideHang = 0.2, 
  cex.colorLabels = 0.8, cex.dendroLabels = 0.6, 
  marAll = c(1, 5, 3, 1), saveMar = TRUE, 
  abHeight = NULL, abCol = "red",main="Bray's clustering (V1 to V3, paired vulva/vagina) with metaData")
dev.off()

############## END CLUSTER

sm <- which(paste(mdata$V9,mdata$V11,mdata$V12) %in% names(table(paste(mdata$V9,mdata$V11,mdata$V12)))[table(paste(mdata$V9,mdata$V11,mdata$V12)) == 2])
#sm <- which(mdata$V11 %in% names(table(mdata$V11[mdata$V9 == "girl"]))[table(mdata$V11[mdata$V9 == "girl"]) > 1])

#sm <- which((mdata$V9 == "girl" & mdata$V12 %in% c("V1","V2","V3")))
#sm <- which(mdata$V11 %in% names(table(mdata$V11))[table(mdata$V11) > 3])
mdata <- mdata[sm,]
splitAbund <- splitAbund[,sm]



pdf(file="Cluster-DynamicV1-V1.pdf",bg="white",width=10,height=7,pointsize=8)

labelColors <- apply((mdata[,c(12,10)]),2,function(x)as.numeric(as.factor(x)))
labelColors <- cbind(labels2colors(labelColors[,1]),labels2colors(labelColors[,2]))

plotDendroAndColors(
  hclust(abdist), 
  colors=labelColors, 
  groupLabels = c("mom/girl","vagina/vulva"), 
  colorText=NULL,
  setLayout = TRUE, 
  autoColorHeight = TRUE, 
  colorHeight = 0.2, 
  dendroLabels = subjectstring, 
  addGuide = FALSE, guideAll = FALSE, 
  guideCount = 50, guideHang = 0.2, 
  cex.colorLabels = 0.8, cex.dendroLabels = 0.6, 
  marAll = c(1, 5, 3, 1), saveMar = TRUE, 
  abHeight = NULL, abCol = "red",main="Bray's clustering with metaData")
dev.off()


cmdfit <- cmdscale(abdist,eig=TRUE, k=2) # k is the number of dim

fcol <- c("pink","red")[as.numeric(as.factor(paste(mdata$V9)))]

# plot solution
pdf(file="mds-Dynamic.pdf",bg="white",width=10,height=7,pointsize=8)
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric MDS, Bray's Distance", type="p",cex=2,pch=19,col=fcol)
#text(x, y, labels = subjectstring,col= fcol)
dev.off()


########################################################################################################################################
## BY VISIT ???
############## MOM VAGINA

mdata <- metaData[metaData$V10 == "vagina" & metaData$V9 == "mom",]
splitAbund <- abund[,metaData$V10 == "vagina"& metaData$V9 == "mom"]


subjectstring <- paste(mdata$V11,mdata$V9,mdata$V10,mdata$V12,sep=".")

abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
fcol <- brewer.pal(3,"Paired")[as.numeric(as.factor(paste(mdata$V9,mdata$V10)))]

tiff(file="VegenHeatTop25-mom-vagina-byVisit.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=8)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[order(rsum,decreasing=TRUE)[1:25],order(mdata$V11,mdata$V12)],col=col,
distfun=vegdist,
labCol=subjectstring[order(mdata$V11,mdata$V12)],
cexRow=0.8,cexCol=0.8,  mar=c(7,5),
dendrogram="column",Rowv=F,Colv=T, 
#labRow=paste(names(ap),": ",ap,"%",sep=""),
keysize=0.8,trace="none",key=T,density.info="none")
dev.off()



############## MOM VAGINA

mdata <- metaData[metaData$V10 == "vagina" & metaData$V9 == "mom",]
splitAbund <- abund[,metaData$V10 == "vagina"& metaData$V9 == "mom"]


subjectstring <- paste(mdata$V11,mdata$V9,mdata$V10,mdata$V12,sep=".")

abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
fcol <- brewer.pal(3,"Paired")[as.numeric(as.factor(paste(mdata$V9,mdata$V10)))]

tiff(file="VegenHeatTop25noRowD-mom-vagina-byVisit.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=8)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[order(rsum,decreasing=TRUE)[1:25],order(mdata$V11,mdata$V12)],col=col,  mar=c(7,5),
distfun=vegdist,labCol=subjectstring[order(mdata$V11,mdata$V12)],cexRow=0.8,cexCol=0.8,
dendrogram="none",Rowv=F,Colv=F,
#labRow=paste(names(ap),": ",ap,"%",sep=""),
keysize=0.8,trace="none",key=T,density.info="none")
dev.off()


############## GIRL VAGINA

mdata <- metaData[metaData$V10 == "vagina" & metaData$V9 == "girl",]
splitAbund <- abund[,metaData$V10 == "vagina"& metaData$V9 == "girl"]


sm <- which(mdata$V11 %in% names(table(mdata$V11))[table(mdata$V11) > 1])
mdata <- mdata[sm,]
splitAbund <- splitAbund[,sm]

subjectstring <- paste(mdata$V11,mdata$V9,mdata$V10,mdata$V12,sep=".")

abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
fcol <- brewer.pal(3,"Paired")[as.numeric(as.factor(paste(mdata$V9,mdata$V10)))]

tiff(file="VegenHeatTop25noRowD-girl-vagina-byVisitWHITE.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=8)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[order(rsum,decreasing=TRUE)[1:25],order(mdata$V11,mdata$V12)],col=col,
distfun=vegdist,labCol=subjectstring[order(mdata$V11,mdata$V12)],cexRow=0.8,cexCol=0.8,  mar=c(7,5),
dendrogram="none",Rowv=F,Colv=F,
colsep=(which(!duplicated(mdata$V11[ order(mdata$V11,mdata$V12)]))-1)[-1],sepcolor="white",
keysize=0.8,trace="none",key=T,density.info="none")
dev.off()


############## GIRL VULVA

mdata <- metaData[metaData$V10 == "vulva" & metaData$V9 == "girl",]
splitAbund <- abund[,metaData$V10 == "vulva"& metaData$V9 == "girl"]

sm <- which(mdata$V11 %in% names(table(mdata$V11))[table(mdata$V11) > 1])
mdata <- mdata[sm,]
splitAbund <- splitAbund[,sm]


subjectstring <- paste(mdata$V11,mdata$V9,mdata$V10,mdata$V12,sep=".")

abdist <- vegdist(t(splitAbund),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
freqAbund <- sweep(splitAbund,2,colSums(splitAbund),"/")
fcol <- brewer.pal(3,"Paired")[as.numeric(as.factor(paste(mdata$V9,mdata$V10)))]

tiff(file="VegenHeatTop25noRowD-girl-vulva-byVisitGREY.tiff",bg="white",width=9,height=6.5,units="in",res=300,pointsize=8)
par(oma=c(2,2,2,5))
heatmap.2(as.matrix(freqAbund)[order(rsum,decreasing=TRUE)[1:25],order(mdata$V11,mdata$V12)],col=col,
distfun=vegdist,labCol=subjectstring[order(mdata$V11,mdata$V12)],cexRow=0.8,cexCol=0.8,  mar=c(7,5),
dendrogram="none",Rowv=F,Colv=F,
colsep=(which(!duplicated(mdata$V11[ order(mdata$V11,mdata$V12)]))-1)[-1],sepcolor="grey",
keysize=0.8,trace="none",key=T,density.info="none")
dev.off()



