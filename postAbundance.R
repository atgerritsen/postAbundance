# source("http://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# install.packages("WGCNA")

library(gplots)
library(WGCNA)
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

## TODO: READ METADATA

# temporary: generate random metadata for testing

metadata = data.frame(samples=colnames(pAmat),
                      group1=sample(c('a','b','c'), size=ncol(pAmat), replace=T),
                      group2=sample(c(1,2,3), size=ncol(pAmat), replace=T),
                      group3=sample(c('tissue1','tissue2','tissue3','tissue4','tissue5'), size=ncol(pAmat), replace=T)
                      )



#---Come up with unique colors per species/classification level.

#-Create subsets for level of classification.


#Figure 1)
pdf(file="Reads_per_sample.pdf", width=8 + .1*ncol(pA), height=8)
    b = barplot(colSums(pA), xaxt='n', main="Reads per sample")
    #axis(side=1, at=b, labels=F)
    text(b, par("usr")[3], labels = colnames(pA), srt = 45, xpd = TRUE, cex=.6, adj=1)
dev.off()


#Figure 2)
# what level were reads classified to?
levels_colors=rainbow(8)
names(levels_colors) = c('domain', 'kingdom','phylum','class','order','family','genus','species')
levels_order = c('domain', 'kingdom','phylum','class','order','family','genus','species')
#count reads by classification level:
cLevelTable = apply(pA, 2, function(x){tapply(x, FUN=sum, I=Level)})
#order table rows by classification level:
cLevelTable = cLevelTable[order(match(rownames(cLevelTable), levels_order), decreasing=T),]
pdf(file="Reads_per_classification_level.pdf", width=8 + .1*ncol(pA), height=8)
    b = barplot(cLevelTable, col=levels_colors[rownames(cLevelTable)], xaxt='n', main='Reads per classification level')
    text(b, par("usr")[3], labels = colnames(pA), srt = 45, xpd = TRUE, cex=.6, adj=1)
    legend('topleft', col=levels_colors[rev(rownames(cLevelTable))], legend=rev(rownames(cLevelTable)), pch=20, cex=1)
dev.off()


Figure 2)
    Stacked barplot for community composition
    -using colors for species

Figure 3)
    Stacked barplot for subset of common species

# Figure N)
# hclust for community membership based on percentages
mdcolors = data.frame(colorRampPalette(c('blue','red'))(length(unique(metadata$group1)))[metadata$group1],
                 colorRampPalette(c('green','yellow'))(length(unique(metadata$group1)))[metadata$group2],
                 colorRampPalette(c('purple','pink'))(length(unique(metadata$group1)))[metadata$group3])

colnames(mdcolors) = colnames(metadata)[2:4]


pdf(file="Hierarchical_clustering_proportions_wards_euclidean.pdf", width=8 + .1*ncol(pA), height=8)
    par(mar=c(1,6, 6, .5))
    layout(matrix(c(1,2), nrow=2), heights=c(20,6))
    hc = hclust(dist(t(pAprop)), method='ward.D2')
    plot(hc, main=paste("Hierarchical clustering using Ward's agglomeration and euclidean distances with proportional values"),
        xlab='', ylab='')
    par(mar=c(1,6,0, .5))
    plotColorUnderTree(hc, colors=mdcolors)
dev.off()

Figure 4)
    MDS-scaling 2d plot of clusters colored by Grouping Factor 1, with pcr controlled by grouping factor 2


Figure 5)
    -Hclust with up to 3 factors as bars underneath

Figure 6-
Metrics:
    -Shannon diversity  -- maybe a barplot with bars colored by grouping factor, ordered by diversity?
    -etc

Protocol for comparing composition of two groups?


Figure N)
    -Time course type plot by subject.
