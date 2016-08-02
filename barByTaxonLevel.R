#!/usr/bin/env Rscript

library("reshape2")

taxa_info <- read.delim("./testData/postAbundance.taxa_info.txt")
Counts <- taxa_info$Total

# Separate each taxon list into desired component parts and combine with read counts; in this case, Family

family <- as.data.frame(cbind(lapply(strsplit(as.character(taxa_info$Taxon_Name), ";"), "[",5), Counts))
# genus <- as.data.frame(cbind(lapply(strsplit(as.character(taxa_info$Taxon_Name), ";"), "[",6), Counts))
family1 <- data.frame(matrix(unlist(family1), nrow=nrow(family1)),stringsAsFactors=FALSE)
# genus <- as.data.frame(cbind(lapply(strsplit(as.character(taxa_info$Taxon_Name), ";"), "[",6), Counts))
colnames(family1) <- c("Names", "Reads")
# colnames(genus1) <- c("Names", "Reads")
family2 <- family1[complete.cases(family1),]
# genus2 <- genus1[complete.cases(genus1),]

# Convert df to characters so that aggregate will work

family2[] <- lapply(family2, function(x) type.convert(as.character(x)))
# genus2[] <- lapply(genus2, function(x) type.convert(as.character(x)))
agg <- aggregate(. ~ Names, family2, sum)
# genus.agg <- aggregate(. ~ Names, genus2, sum)

# head(agg)
#                   Names Reads
# 1   f__Acetobacteraceae   162
# 2 f__Acidaminococcaceae     2
# 3  f__Acidimicrobiaceae    68
# 4   f__Actinomycetaceae   427
# 5      f__Aerococcaceae   321
# 6     f__Aeromonadaceae    68

# Order by decreasing abundance, take reads > 1000, and write to pdf

agg2 <- agg[order(-agg$Reads),]
topRows <- agg2[agg2$Reads>1000,]

pdf(file="Reads_per_FamilyLevel.pdf", width= 18, height=8)
b <- barplot(agg2$Reads, main = "Reads for each family", ylab="Number of Reads")
text(b, par("usr")[3], labels = agg2$Names, srt = 45, xpd = TRUE, cex=.6, adj=1)

tops <- barplot(topRows$Reads, main = "Top Families", ylab="Number of Reads")
text(tops, par("usr")[3], labels = topRows$Names, srt = 45, xpd = TRUE, cex=.8, adj=1)


dev.off()

