# Simulations done in the accompanying paper.
# Cene Fišer & Roman Luštrik, 20.5.2014

# Calculate Euclidean distance pair-wise, sort by size and take shortest num.spec-1.
calcEuclid <- function(cmbn, data) {
   num.spec <- nrow(data)
   
   euc <- apply(cmbn, MARGIN = 2, FUN = function(x, data) {
         sqrt(sum(apply(data[x, , drop = FALSE], MARGIN = 2, FUN = diff)^2))
      }, data = data)
   euc <- euc[order(euc)]
   sd.euc <- head(euc, num.spec-1)
   list(sd = sd(sd.euc), R = sum(sd.euc))
}

# get distribution of phantom species
getDistro <- function(phant.spec, num.spec, nsim, cmbn) {
   distro <- replicate(nsim, {
         # For NULL model:
         # For each sampled "community" calculate standard deviation and range of traits, which
         # is the sum of all standard deviations.
         sub.spec <- phant.spec[sample(1:nrow(phant.spec), size = num.spec), , drop = FALSE]
         sub.sd <- calcEuclid(cmbn = cmbn, data = sub.spec)
         
         out <- list(sd = sub.sd$sd, R = sub.sd$R)
      }, simplify = FALSE)
   
   sdNR <- sapply(distro, "[[", "sd")
   R <- sapply(distro, "[[", "R")
   sdNR_byR <- mean(sdNR/R)
   
   list(sd = sdNR, mean_sdNR_byR = sdNR_byR, sd_sdNR_byR_null = sd(sdNR), R = R)
}

calcAxis <- function(data.comm, tipi, community, os, nsim) {
   data <- data.comm[, tipi, drop = FALSE]
   num.spec <- nrow(data)
   
   bs.range <- sapply(data, range)
   phant.spec <- apply(bs.range, MARGIN = 2, FUN = function(x) {
         runif(100, min = x[1], max = x[2]) # generate number of phantom species
      })
   cmbn <- combn(1:num.spec, m = 2)
   
# Calculate Euclidean distance between species. Number of species
# is equal to the number in the community.
   bs.distro <- getDistro(phant.spec = phant.spec, num.spec = num.spec, nsim = nsim, cmbn = cmbn)
   bs.sd <- calcEuclid(cmbn = cmbn, data = data)
   pval <- prop.table(table(bs.distro$sd < bs.sd$sd))["TRUE"]
   
   # TEI and TCI are calculate as follows (see page 2449 for formulas (3) and (4)).
   # We calculate euclidean distance between all pairs of species and select the best
   # n-1 pairs (where n is number of species) for real and every simulated community.
   # These are the values used to calculate mean and sd of observations and null model.
   # Range is calculated by summing these 'best' euclidean distances.
   # Indices are calculated by the following formulae.
   # TEI = -((sdND/R - mean(sdND/R)_null)/sd(sdND/R)_null)
   # TCI = -((range_obs - mean(range_null))/sd(range)_null)
   
   # Calculate TEI as in Ingram and Shurin 2009 (see page 2449).
   R.obs <- bs.sd$sd / bs.sd$R 
   R.null <- bs.distro$mean_sdNR_byR
   SD.null <- bs.distro$sd_sdNR_byR_null
   TEI <- -((R.obs - R.null)/SD.null)

   # Calculate TCI as in Ingram and Shurin 2009 (see page 2449).
   tci.robs <- bs.sd$R
   tci.range.null <- mean(bs.distro$R)
   tci.sd.range.null <- sd(bs.distro$R)
   TCI <- -((tci.robs - tci.range.null)/tci.sd.range.null)
   
   result <- list(distribution = bs.distro$sd, real.sd = bs.sd, pval = pval, 
      TEI = TEI, TCI = TCI, os = os, community = community)
}

