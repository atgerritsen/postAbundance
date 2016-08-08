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

# read in metadata groupings
FecalMeta.2014 <- read.csv("./testData/FecalMeta.2014.csv", header=T, sep='\t')
FecalMeta.2014$Taxon_Name <- sub("^", "X", FecalMeta.2014$Taxon_Name)
foodTreatment <- FecalMeta.2014$Treatment
tissueOrigin <- FecalMeta.2014$Tissue

# metadata = data.frame(samples=colnames(pAmat),
#                       group1=sample(c('a','b','c'), size=ncol(pAmat), replace=T),
#                       group2=sample(c(1,2,3), size=ncol(pAmat), replace=T),
#                       group3=sample(c('tissue1','tissue2','tissue3','tissue4','tissue5'), size=ncol(pAmat), replace=T)
#                       )


#---Come up with unique colors per species/classification level.

#-Create subsets for level of classification.

# Might be more useful to 
# Figure 2)
#     Stacked barplot for community composition
#     -using colors for species

# Figure 3)
#     Stacked barplot for subset of common species

# Figure N)
# hclust for community membership based on unfiltered percentages
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

# Figure N+1)  MDS-scaling 2d plot of clusters colored by Grouping Factor 1, with pcr controlled by grouping factor 2
loc = cmdscale(dist(t(pAprop)), k=3)
x = loc[,1]
y = loc[,2]
z = loc[,3]
pdf(file="Metric_mds_scaling_2d_euclidean.pdf", width=8, height=8)
  plot(x,y, type='p', xlab='', ylab='', asp=1, axes=F, main="Metric Multidimensional Scaling plot", pch=20, col=mdcolors[,1],
    xlim=range(x) * 1.2, ylim=range(y) * 1.01)
  text(x,y, rownames(loc), cex=0.4, adj=1.1)
  legend('topright', col=unique(mdcolors[,1]), legend=unique(metadata$group1), pch=20)
dev.off()

## Make an interactive 3d model in a web page:
library(rgl)
plot3d(x, y, z, col=mdcolors[,1], radius=.01, type='s')

filename <- writeWebGL(dir = file.path('./', "webGL"), 
                       width = 2000, reuse = TRUE)
######


# Figure 5)
#     -Hclust with up to 3 factors as bars underneath

# Figure 6-
# Metrics:
#     -Shannon diversity  -- maybe a barplot with bars colored by grouping factor, ordered by diversity?
#     -etc

# Protocol for comparing composition of two groups?


# Figure N)
#     -Time course type plot by subject.