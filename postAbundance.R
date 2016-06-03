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
                      group3=sample(c('tissue1','tissue2','tissue3'), size=ncol(pAmat), replace=T)
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
cLevelTable = apply(pA, 2, function(x){tapply(x, FUN=sum, I=Level)})
b = barplot(cLevelTable, col=levels_colors[rownames(cLevelTable)], xaxt='n', main='Reads per classification level')
text(b, par("usr")[3], labels = colnames(pA), srt = 45, xpd = TRUE, cex=.6, adj=1)
legend('topleft', col=levels_colors[rownames(cLevelTable)], legend=rownames(cLevelTable), pch=20, cex=1)

heatmap.2(pAprop[rare.idx,], trace='none', margins=c(8,8))

Figure 2)
    Stacked barplot for community composition
    -using colors for species

Figure 3)
    Stacked barplot for subset of common species

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
