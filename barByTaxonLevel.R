#!/usr/bin/env Rscript

library("reshape2")

taxa_info <- read.delim("./testData/postAbundance.taxa_info.txt")
Counts <- taxa_info$Total

# Separate each taxon list into desired component parts and combine with read counts; in this case, Family

family <- as.data.frame(cbind(lapply(strsplit(as.character(taxa_info$Taxon_Name), ";"), "[",5), Counts))
family1 <- data.frame(matrix(unlist(family1), nrow=nrow(family1)),stringsAsFactors=FALSE)
colnames(family1) <- c("Names", "Reads")
family2 <- family1[complete.cases(family1),]

# Convert df to characters so that aggregate will work

family2[] <- lapply(family2, function(x) type.convert(as.character(x)))
agg <- aggregate(. ~ Names, family2, sum)

# head(agg)
#                   Names Reads
# 1   f__Acetobacteraceae   162
# 2 f__Acidaminococcaceae     2
# 3  f__Acidimicrobiaceae    68
# 4   f__Actinomycetaceae   427
# 5      f__Aerococcaceae   321
# 6     f__Aeromonadaceae    68

# Order by decreasing abundance and write to file

agg2 <- agg[order(-agg$Reads),]
pdf(file="Reads_per_FamilyLevel.pdf", width= 17, height=8)
b <- barplot(agg2$Reads)
text(b, par("usr")[3], labels = agg2$Names, srt = 45, xpd = TRUE, cex=.6, adj=1)
dev.off()
