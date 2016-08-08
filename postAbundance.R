# source("http://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# install.packages("WGCNA")

library(gplots)
# library(WGCNA)
library(RColorBrewer)

## Read data and set it up for plotting
pA = read.table("./testData/postAbundance.abundance.txt", header=T, as.is=T, sep='\t', comment.char='')
Taxon_Name = pA$Taxon_Name
Level = pA$Level
pA = pA[,-c(1,2)]
pAmat = as.matrix(pA)
rownames(pAmat) = Taxon_Name
pAprop = t(t(pAmat)/colSums(pAmat))


# build a filter which removes bacteria which are at < 1% in 1% or less of samples:
rare.idx = rowSums(pAprop > .01) > (ncol(pAprop) * .01)

#---Come up with unique colors per species/classification level.

# Figure 1 - Reads per sample, with samples < 10000 reads higlighted
pdf(file="Data Summaries", width=8 + .1*ncol(pA), height=11)
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
    c = barplot(cLevelTable, col=levels_colors[rownames(cLevelTable)], xaxt='n', main='Reads per classification level')
    text(c, par("usr")[3], labels = colnames(pA), srt = 45, xpd = TRUE, cex=.6, adj=1)
    legend('topleft', col=levels_colors[rev(rownames(cLevelTable))], legend=rev(rownames(cLevelTable)), pch=20, cex=1)

taxa_info <- read.delim("./testData/postAbundance.taxa_info.txt")
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
family.topRows <- family.agg2[family.agg2$Reads>100,]
genus.agg2 <- genus.agg[order(-genus.agg$Reads),]
genus.topRows <- genus.agg2[genus.agg2$Reads>1000,]

# Figure 3 and 4 - Reads per taxonomic level, top genera > 1000 reads
# pdf(file="Reads_per_FamilyLevel.pdf", width= 18, height=8)
f <- barplot(family.topRows$Reads, main = "Reads for each family > 100", ylab="Number of Reads", col="firebrick")
text(f, par("usr")[3], labels = family.topRows$Names, srt = 45, xpd = TRUE, cex=.8, adj=1)
genus.tops <- barplot(genus.topRows$Reads, main = "Top genera with reads > 1000", ylab="Number of Reads", col="darkslategray2")
text(genus.tops, par("usr")[3], labels = genus.topRows$Names, srt = 45, xpd = TRUE, cex=.75, adj=1)
 

# Figure 5 and 5.1)
#     Stacked barplot for community composition with colors for species and shannon-weaver index for each sample

div <- read.table("./testData/postAbundance.abundance.txt", header=T, as.is=T, sep='\t', comment.char='')
div.genus <- subset(div, div$Level == "genus")
names <- div.genus$Taxon_Name

genus.colors <- rainbow(nrow(div.genus))
names(genus.colors) <- div.genus$Taxon_Name
levels_order <- div.genus$Taxon_Name

legcols <- cbind(div.genus$Taxon_Name, rainbow(nrow(div.genus)))
suum <- rowSums(div.genus[,-c(1,2)])
allcols <- as.data.frame(cbind(legcols, suum))
allcols1 <- allcols[order(-suum),]
leg <- head(allcols1, n=10)
cols <-as.character(leg$V2)
nams <-as.character(leg$V1)

genes <- apply(div.genus[,-c(1,2)], 2, function(x){tapply(x, FUN=sum, I=names)})
genes2 <- genes[order(match(rownames(genes), levels_order), decreasing=F),]

m <- barplot(genes2, col=genus.colors[rownames(genes2)], xaxt='n', main='Genera in each sample', ylim=c(0,100000))
text(m, par("usr")[3], labels = colnames(div.genus), srt = 45, xpd = TRUE, cex=.7, adj=1)
legend('topleft', legend=nams, col=cols, pch=20, cex=1, pt.cex = 1.2)

dev.off()


