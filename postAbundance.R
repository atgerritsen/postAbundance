
































                       width = 2000, reuse = TRUE)
                 colorRampPalette(c('green','yellow'))(length(unique(metadata$group1)))[metadata$group2],
                 colorRampPalette(c('purple','pink'))(length(unique(metadata$group1)))[metadata$group3])
        xlab='', ylab='')
    #axis(side=1, at=b, labels=F)
    -etc
    -Hclust with up to 3 factors as bars underneath
    -Shannon diversity  -- maybe a barplot with bars colored by grouping factor, ordered by diversity?
    -Time course type plot by subject.
    -using colors for species
    b = barplot(cLevelTable, col=levels_colors[rownames(cLevelTable)], xaxt='n', main='Reads per classification level')
    b = barplot(colSums(pA), xaxt='n', main="Reads per sample")
    hc = hclust(dist(t(pAprop)), method='ward.D2')
    layout(matrix(c(1,2), nrow=2), heights=c(20,6))
    legend('topleft', col=levels_colors[rev(rownames(cLevelTable))], legend=rev(rownames(cLevelTable)), pch=20, cex=1)
    par(mar=c(1,6, 6, .5))
    par(mar=c(1,6,0, .5))
    plot(hc, main=paste("Hierarchical clustering using Ward's agglomeration and euclidean distances with proportional values"),
    plotColorUnderTree(hc, colors=mdcolors)
    Stacked barplot for community composition
    Stacked barplot for subset of common species
    text(b, par("usr")[3], labels = colnames(pA), srt = 45, xpd = TRUE, cex=.6, adj=1)
    text(b, par("usr")[3], labels = colnames(pA), srt = 45, xpd = TRUE, cex=.6, adj=1)
    xlim=range(x) * 1.2, ylim=range(y) * 1.01)
  legend('topright', col=unique(mdcolors[,1]), legend=unique(metadata$group1), pch=20)
  plot(x,y, type='p', xlab='', ylab='', asp=1, axes=F, main="Metric Multidimensional Scaling plot", pch=20, col=mdcolors[,1],
  text(x,y, rownames(loc), cex=0.4, adj=1.1)
#                       group1=sample(c('a','b','c'), size=ncol(pAmat), replace=T),
#                       group2=sample(c(1,2,3), size=ncol(pAmat), replace=T),
#                      )
#                      group3=sample(c('tissue1','tissue2','tissue3','tissue4','tissue5'), size=ncol(pAmat), replace=T)
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# build a filter which removes bacteria which are at < 1% in 1% or less of samples:
# Figure N)
# Figure N+1)  MDS-scaling 2d plot of clusters colored by Grouping Factor 1, with pcr controlled by grouping factor 2
# hclust for community membership based on unfiltered percentages
# install.packages("WGCNA")
# metadata = data.frame(samples=colnames(pAmat),
# read in metadata and format sample names column to match figures
# source("http://bioconductor.org/biocLite.R")
# what level were reads classified to?
## Make an interactive 3d model in a web page:
## Read data and set it up for plotting
######
#---Come up with unique colors per species/classification level.
#-Create subsets for level of classification.
#count reads by classification level:
#Figure 1)
#Figure 2)
#order table rows by classification level:
cLevelTable = apply(pA, 2, function(x){tapply(x, FUN=sum, I=Level)})
cLevelTable = cLevelTable[order(match(rownames(cLevelTable), levels_order), decreasing=T),]
colnames(mdcolors) = colnames(metadata)[2:4]
dev.off()
dev.off()
dev.off()
dev.off()
FecalMeta.2014 <- read.csv("./testData/postAbundance/FecalMeta.2014.csv", header=T)
FecalMeta.2014$Taxon_Name <- sub("^", "X", FecalMeta.2014$Taxon_Name)
Figure 2)
Figure 3)
Figure 5)
Figure 6-
Figure N)
filename <- writeWebGL(dir = file.path('./', "webGL"), 
Level = pA$Level
levels_colors=rainbow(8)
levels_order = c('domain', 'kingdom','phylum','class','order','family','genus','species')
library(gplots)
library(RColorBrewer)
library(rgl)
library(WGCNA)
loc = cmdscale(dist(t(pAprop)), k=3)
mdcolors = data.frame(colorRampPalette(c('blue','red'))(length(unique(metadata$group1)))[metadata$group1],
Metrics:
names(levels_colors) = c('domain', 'kingdom','phylum','class','order','family','genus','species')
pA = pA[,-c(1,2)]
pA = read.table("./testData/postAbundance.abundance.txt", header=T, as.is=T, sep='\t', comment.char='')
pAmat = as.matrix(pA)
pAprop = t(t(pAmat)/colSums(pAmat))
pdf(file="Hierarchical_clustering_proportions_wards_euclidean.pdf", width=8 + .1*ncol(pA), height=8)
pdf(file="Metric_mds_scaling_2d_euclidean.pdf", width=8, height=8)
pdf(file="Reads_per_classification_level.pdf", width=8 + .1*ncol(pA), height=8)
pdf(file="Reads_per_sample.pdf", width=8 + .1*ncol(pA), height=8)
plot3d(x, y, z, col=mdcolors[,1], radius=.01, type='s')
Protocol for comparing composition of two groups?
rare.idx = rowSums(pAprop > .01) > (ncol(pAprop) * .01)
rownames(pAmat) = Taxon_Name
Taxon_Name = pA$Taxon_Name
treatment <- FecalMeta.2014$Treatment
x = loc[,1]
y = loc[,2]
z = loc[,3]