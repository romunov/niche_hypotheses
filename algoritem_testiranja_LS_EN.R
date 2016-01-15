# iščemo funkcionalno redundanco v ekosistemu (functional redundancy) -

library(xlsx)
library(ggplot2)
library(gridExtra)
source("functions.R")

# uvozimo podatke
spe <- read.xlsx("./data/data_20142005.xlsx", sheetName = "morphospace1")
coms <- read.xlsx("./data/data_20142005.xlsx", sheetName = "communities", encoding = "UTF-8")

# naredimo clustering samo s telesno dolžino in gnatopodi
spe.hc <- hclust(dist(spe[, grepl("gpI|log\\.body", names(spe))]), method = "ward.D2")
plot(spe.hc, labels = spe$Species)

# 1. najprej naredimo clustering in režemo tako, da dobimo dve skupini
spe$group <- cutree(spe.hc, k = 2)

# subset of only key variables
spe <- spe[, grepl("log\\.bodyLenght|gp|Species|ecomorph|group", names(spe))]

# 2. naredimo nulti model razlike med klastri (guildi) iz vseh vrst tako, da iz naključnih parov izračunamo
#    evklidske razdalje -> histogram

# find how many unique pairs can be generated
# all possible combinations
# all.combos <- combn(x = 1:nrow(spe), m = 2)
# contrained between cluster pairs
all.combos <- t(expand.grid(predator = 1:10, microf = 11:18))

# for each pair, calculate distance
null.model.diff.clust <- apply(all.combos, MARGIN = 2, function(x, d) {
    euclidDistance(
        d[sample(1:nrow(d), size = 2, replace = FALSE), ]
    )
}, d = spe[, !grepl("Species|ecomorph", names(spe))])

# draw histogram of counts/density
ggplot(data.frame(x = null.model.diff.clust), aes(x = x)) +
    theme_bw() +
    geom_histogram()

# 3. za vsako združbo (8) izračunamo razliko med centroidoma klastrov (zelen x),
#    ta delta X narišeš (izračunaš p) na histogram. pričakujemo, da bo na zgornjem koncu histograma
#    to ponovimo za vse združbe (8x)
diff.btw.clust <- sapply(levels(coms$Locality), FUN = actualDiffClusters, coms = coms)


diff.btw.clust <- data.frame(distance = diff.btw.clust, community = names(diff.btw.clust))
rownames(diff.btw.clust) <- NULL

#  calculate p-value
diff.btw.clust$p.val <- sapply(diff.btw.clust$distance, FUN = function(x) prop(table(x < null.model.diff.clust)))
diff.btw.clust$label <- sprintf("%s: %f", diff.btw.clust$community, round(diff.btw.clust$p.val, 3))
diff.btw.clust

last_plot() +
    geom_vline(data = diff.btw.clust, aes(xintercept = distance)) +
    geom_text(data = diff.btw.clust, aes(x = distance, y = c(12, 12, 10, 10, 12, 10, 7.5, 12) - 5, label = label),
              angle = 90, vjust = 1, size = 3)

ggsave("./figure/pairs_from_constrained_clusters.pdf", width = 5, height = 5)

t.test(x = null.model.diff.clust, y = diff.btw.clust$distance)
var.test(x = null.model.diff.clust, y = diff.btw.clust$distance)

sd(diff.btw.clust$distance)
sd(null.model.diff.clust)

# 4. za vsak klastr:
#    preštejemo število vrst iz klastru v združbi in glede na število (2, 3) naredimo nulti model - izbiramo 2 ali 3 vrste
#    iz vseh vrst glavnega drevesa
null.model.2 <- apply(X = combn(x = 1:nrow(spe), m = 2), MARGIN = 2, FUN = nullCluster2, spe = spe)

ggplot(data.frame(x = null.model.2), aes(x = x)) +
    theme_bw() +
    geom_histogram()

# tole nariši na nulti model skupaj s p-ji
list.o.specs2 <- list(c("N_grandii", "N_microcerberus"),
     c("N_longidactylus", "N_pupetta"),
     c("Carinurella", "N_transitivus"),
     c("N_kochianus", "N_longidactylus"),
     c("N_multipennatus", "N_serbicus"),
     c("N_aquilex", "N_fontanus"))

dist.2spec <- sapply(X = list.o.specs2, FUN = function(x, spe) {
         euclidDistance(spe[spe$Species %in% x, !grepl("Species|ecomorph|group", x = names(spe))])
     }, spe = spe)

dist.2specs <- data.frame(dists = dist.2spec, specs = sapply(list.o.specs2, paste, collapse = "-"))
dist.2specs$pval <- sapply(dist.2specs$dists, FUN = function(x) prop(table(x < null.model.2)))
dist.2specs$label <- sprintf("%s: %.3f", dist.2specs$specs, dist.2specs$pval)

last_plot() +
    geom_vline(data = dist.2specs, aes(xintercept = dists)) +
    geom_text(data = dist.2specs, aes(x = dists, y = 7, label = label),
              angle = 90, vjust = 1, size = 3)

ggsave("./figure/unconstrained_2_species.pdf", width = 5, height = 5)

t.test(x = null.model.2, y = c(ng.nm, nl.np, c.nt, nk.nl, nm.ns, na.nf))

# 5. izračunamo minimum spanning tree (MST) iz evklidskih razdalj oz. evklidsko
#    razdaljo za dve vrsti, in to narišemo na histogram oz. izračunamo p (10 grafov)

null.model.3 <- apply(X = combn(x = 1:nrow(spe), m = 3), MARGIN = 2, FUN = nullCluster3, spe = spe)

ggplot(data.frame(x = null.model.3), aes(x = x)) +
    theme_bw() +
    geom_histogram()

list.o.specs3 <- list(c("N_kenki_5", "N_multipennatus", "N_pectinicauda"),
                      c("N_pupetta", "N_labacensis", "N_longidactylus"),
                      c("N_grandii", "N_microcerberus", "N_dolenianesis"))

dist.3spec <- sapply(X = list.o.specs3,
       FUN = function(x, spe) mst(spe[spe$Species %in% x, !grepl("Species|ecomorph|group", x = names(spe))]),
       spe = spe)

dist.3specs <- data.frame(dists = dist.3spec, specs = sapply(list.o.specs3, paste, collapse = "-"))
dist.3specs$pval <- sapply(dist.3specs$dists, FUN = function(x) prop(table(x < null.model.3)))
dist.3specs$label <- sprintf("%s: %.3f", dist.3specs$specs, dist.3specs$pval)

last_plot() +
    geom_vline(data = dist.3specs, aes(xintercept = dists)) +
    geom_text(data = dist.3specs, aes(x = dists, y = 40, label = label),
              angle = 90, vjust = 1.2, size = 2.5)

ggsave("./figure/unconstrained_3_species.pdf", width = 5, height = 5)

t.test(x = null.model.3, y = c(kmp, pll, gmd))
