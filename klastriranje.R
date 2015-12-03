library(xlsx)
library(ggplot2)
library(gridExtra)

# uvozimo podatke
spe <- read.xlsx("./data/data_20142005.xlsx", sheetName = "morphospace1")

# izberemo ven samo telesno dolžino in znake gpI*
sub.spe <- spe[, grepl("gpI|log\\.body", names(spe))]

# naredimo clustering samo s telesno dolžino in gnatopodi
spe.hc <- hclust(dist(sub.spe), method = "ward.D2")

plot(spe.hc, labels = spe$Species)

# naredimo clustering z vsemi znaki
plot(hclust(dist(spe[, c(-1, -19)]), method = "ward.D2"), labels = spe$Species)

