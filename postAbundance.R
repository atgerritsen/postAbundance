library(ggplots)

## Read data
pA = read.table("./testData/postAbundance.abundance.txt", header=T, as.is=T, sep='\t', comment.char='')
Taxon_Name = pA$Taxon_Name
Level = pA$Level
pA = pA[,-c(1,2)]
pA2 = as.matrix(pA)
heatmap(pA2)




---Come up with unique colors per species/classification level.

-Create subsets for level of classification.


Figure 1)
    Reads per sample

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
