# source("http://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# install.packages("WGCNA")

library(gplots)
library(vegan)
library(RColorBrewer)
library(colorRamps)

## Read data and set it up for plotting
pA = read.table("./well12A2016.abundance.txt", header=T, as.is=T, sep='\t', comment.char='')
Taxon_Name = pA$Taxon_Name
Level = pA$Level
pA = pA[,-c(1,2)]
pAmat = as.matrix(pA)
rownames(pAmat) = Taxon_Name
pAprop = t(t(pAmat)/colSums(pAmat))


# build a filter which removes bacteria which are at < 1% in 1% or less of samples:
rare.idx = rowSums(pAprop > .01) > (ncol(pAprop) * .01)

# Figure 1 - Reads per sample, with samples < 10000 reads higlighted
pdf(file="Data Summaries 2016.pdf", width=8 + .7*ncol(pA), height=11)
    smallSums <- colSums(pA) >= 10000
    cols <- c("yellow", "forestgreen")
    b = barplot(colSums(pA), xaxt='n', main="Total reads per sample", col=cols[smallSums+1])
    #axis(side=1, at=b, labels=F)
    text(b, par("usr")[3], labels = colnames(pA), srt = 45, xpd = TRUE, cex=.6, adj=c(1,0))
    legend('topleft', fill="yellow", col="yellow", legend = "Samples < 10,000 reads")

# Figure 2 - What level were reads classified to?
levels_colors=rainbow(8)
names(levels_colors) = c('domain', 'kingdom','phylum','class','order','family','genus','species')
levels_order = c('domain', 'kingdom','phylum','class','order','family','genus','species')
# Count reads by classification level:
cLevelTable = apply(pA, 2, function(x){tapply(x, FUN=sum, I=Level)})
# Order table rows by classification level:
cLevelTable = cLevelTable[order(match(rownames(cLevelTable), levels_order), decreasing=T),]
# pdf(file="Reads_per_classification_level.pdf", width=8 + .1*ncol(pA), height=8)
    par(mfrow=c(1, 1), mar=c(5, 5, 4, 3))
    c = barplot(cLevelTable, col=levels_colors[rownames(cLevelTable)], xaxt='n', main='Reads per classification level', ylim = c(0,3000000))
    #c = barplot(try, col=levels_colors[rownames(try)], xaxt='n', main='Reads per classification level')
    text(c, par("usr")[3], labels = colnames(try), srt = 45, xpd = TRUE, cex=.6, adj=1)
    legend('topright', col=levels_colors[rev(rownames(cLevelTable))], legend=rev(rownames(cLevelTable)), pch=20, cex=1.2, ncol=1, bty="n")
    
taxa_info <- read.delim("./well12A2016.taxa_info.txt")
Counts <- taxa_info$Total
family <- as.data.frame(cbind(lapply(strsplit(as.character(taxa_info$Taxon_Name), ";"), "[",5), Counts))
genus <- as.data.frame(cbind(lapply(strsplit(as.character(taxa_info$Taxon_Name), ";"), "[",6), Counts))
family1 <- data.frame(matrix(unlist(family), nrow=nrow(family)),stringsAsFactors=FALSE)
genus1 <- data.frame(matrix(unlist(genus), nrow=nrow(genus)),stringsAsFactors=FALSE)
colnames(family1) <- c("Names", "Reads")
colnames(genus1) <- c("Names", "Reads")
family2 <- family1[complete.cases(family1),]
genus2 <- genus1[complete.cases(genus1),]

# Convert df to characters so that aggregate will work
family2[] <- lapply(family2, function(x) type.convert(as.character(x)))
genus2[] <- lapply(genus2, function(x) type.convert(as.character(x)))
family.agg <- aggregate(. ~ Names, family2, sum)
genus.agg <- aggregate(. ~ Names, genus2, sum)

# Order by decreasing abundance, take reads > 1000, and write to pdf
family.agg2 <- family.agg[order(-family.agg$Reads),]
family.topRows <- family.agg2[family.agg2$Reads>5000,]
genus.agg2 <- genus.agg[order(-genus.agg$Reads),]
genus.topRows <- genus.agg2[genus.agg2$Reads>10000,]

# Figure 3 and 4 - Reads per taxonomic level, top genera > 1000 reads
# pdf(file="Reads_per_FamilyLevel.pdf", width= 18, height=8)
par(mar=c(11,3,4,0))
f <- barplot(family.topRows$Reads, main = "Reads for each family > 5000, all samples", ylab="Number of Reads", col="firebrick", width=8 + .1*ncol(family.topRows))
axis(side=1, at=f, labels=NA, lwd=0, lwd.ticks=1, pos=1, tck=-.009)
text(f, par("usr")[3], labels = family.topRows$Names, srt = 45, xpd = TRUE, cex=.7, adj=1)
genus.tops <- barplot(genus.topRows$Reads, main = "Top genera with reads > 10000, all samples", ylab="Number of Reads", col="darkslategray2")
axis(side=1, at=genus.tops, labels=NA, lwd=0, lwd.ticks=1, pos=1, tck=-.009)
text(genus.tops, par("usr")[3], labels = genus.topRows$Names, srt = 45, xpd = TRUE, cex=.7, adj=1)


# Figure 5 and 5.1)
#     Stacked barplot for community composition with colors for species and shannon-weaver index for each sample

#read and cleanup abundance table:
header <- read.table("./well12A2016.abundance.txt", header=F, as.is=T, sep='\t', comment.char='', nrow=1)
header= unlist(header[1,])
header= c(header[1:2], paste("s", header[3:length(header)], sep='_'))
div <- read.table("./well12A2016.abundance.txt", header=T, as.is=T, sep='\t', comment.char='')
colnames(div) = header

div.genus <- subset(div, div$Level == "genus")

# build a filter which removes bacteria which are at < 1% in 1% or less of samples:
pAg = div.genus[,-c(1,2)]
pAmat = as.matrix(pAg)
rownames(pAmat) = div.genus$Taxon_Name
pAprop = t(t(pAmat)/colSums(pAmat))

# build a filter which keeps bacteria which are at > 1% in 1% or more of samples:
common.idx = rowSums(pAprop > .01) > (ncol(pAprop) * .01)

div.genus.common = div.genus[common.idx,]
div.genus.common = div.genus.common[order(div.genus.common$Taxon_Name, decreasing=F), ]
genus.colors = primary.colors(nrow(div.genus.common))
names(genus.colors) = div.genus.common$Taxon_Name
dgc.mat = as.matrix(div.genus.common[,-c(1,2)])

layout(mat=matrix(nrow=1, c(1,2)), widths=c(8,2))
par(las=1)
par(mar=c(8,4.5,1,1))
m <- barplot(dgc.mat, col=genus.colors, xaxt='n', main='Genera in each sample', width=8 + .1*ncol(dgc.mat))
text(m, par("usr")[3], labels = colnames(dgc.mat), srt = 45, xpd = TRUE, cex=.7, adj=1)
par(mar=c(.1,.1,.1,.1))
plot(NA, xaxt='n', yaxt='n', ylim=c(-1,1), xlim=c(-1,1), xlab='', ylab='')
legend('topleft', legend=names(genus.colors), col=genus.colors, pch=20, cex=.65, pt.cex = 1)

dev.off()

### Write
genus <- div.genus[,-c(1,2)]
genust <- t(genus)
shann <- diversity(genust)
write.table(shann, file= "shannon.diversity2016.txt", quote=F, sep='\t', row.names=T, col.names=F)
