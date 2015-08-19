# Test of competition in the intersticial
# Cene Fišer & Roman Luštrik, 20.5.2014

setwd("./niche_within_niche")

source("functions_nisa_v_nisi_po_zdruzbah.R")
library(ggplot2)
library(dplyr)

data <- read.table("./data/data_20142005.csv", header = TRUE, sep = ",", check.names = FALSE)
comms <- read.table("./data/data_20142005_communities.csv", header = TRUE, sep = ",")
data$ecomorph <- NULL

# used characteristics
tipi.body <- "log(bodyLenght)"
tipi.shape <- c("aI/body_length", "aII/body_length", "coxaII-shape", "coxaIII-shape", "ppV/body_length",
   "ppVI/body_length", "ppVII/body_length", "width/length_ppV/2", "width/length_ppVI/2", 
   "width/length_ppVII/2")
tipi.gnath <- c("gpI-wrist", "gpI-size", "gpI_angle", "gpII-wrist", "gpII-size", "gpII_angle")
tipi.all <- c(tipi.body, tipi.shape, tipi.gnath)

# constants
nsim <- 1000
comm.split <- split(comms, comms$Locality)

out <- vector("list", length(comm.split))
names(out) <- names(comm.split)

for (community in names(comm.split)) {
   data.comm <- data[data$Species %in% comm.split[[community]]$Species, ]
   data.comm$Species <- NULL
   
   # body length
   bl <- calcAxis(data.comm = data.comm, tipi = tipi.body, 
      community = community, os = "body_length", nsim = nsim)
   
   # body shape
   bs <- calcAxis(data.comm = data.comm, tipi = tipi.shape, 
      community = community, os = "body_shape", nsim = nsim)
   
   # gnathopods
   gn <- calcAxis(data.comm = data.comm, tipi = tipi.gnath, 
      community = community, os = "gnathopods", nsim = nsim)
   
   # all
   tot <- calcAxis(data.comm = data.comm, tipi = tipi.all, 
      community = community, os = "all", nsim = nsim)
   
   meta.bl <- data.frame(TEI = bl$TEI, TCI = bl$TCI, os = bl$os, community = bl$community)
   meta.bs <- data.frame(TEI = bs$TEI, TCI = bs$TCI, os = bs$os, community = bs$community)
   meta.gn <- data.frame(TEI = gn$TEI, TCI = gn$TCI, os = gn$os, community = gn$community)
   meta.tot <- data.frame(TEI = tot$TEI, TCI = tot$TCI, os = tot$os, community = tot$community)
   meta <- rbind(meta.bl, meta.bs, meta.gn, meta.tot)
   
   out[[community]] <- list(body_length = bl, body_shape = bs, gnathopods = gn, all = tot, meta = meta)
   
}

final <- vector("list", length(out))
names(final) <- names(out)

for (i in names(out)) {
   x <- out[[i]]
   rslt <- lapply(x[!names(x) == "meta"], FUN = function(y) {
         data.frame(value = y$distribution, real.sd = y$real.sd$sd, axis = y$os, 
            community = y$community, pval = unname(y$pval))
      })
   rslt <- do.call("rbind", rslt)
   rslt$community <- i
   final[[i]] <- rslt
}
final <- do.call("rbind", final)

put.vals <- final %.% group_by(community, axis) %.% summarise(pval = unique(pval), real.sd = unique(real.sd))

ggplot(final, aes(x = value)) +
   theme_bw() +
   geom_density() +
   geom_vline(data = put.vals, aes(xintercept = real.sd)) +
   geom_text(data = put.vals, aes(x = real.sd, y = 20, label = pval), hjust = -0.25, size = 3) +
   facet_grid(community ~ axis)
ggsave("competition_of_comms.pdf", width = 10, height = 10)

# Extract results for meta analysis adapted from Ingram and Shurin 2009.
meta <- sapply(out, "[[", "meta", simplify = FALSE)
meta <- do.call("rbind", meta)
rownames(meta) <- NULL

library(gridExtra)
pdf("meta_table.pdf", width = 7, height = 10)
grid.table(meta)
dev.off()