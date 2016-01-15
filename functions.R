## ---- functions ----

# Function to calculate ratio between TRUE/FALSE
# x = boolean vector
prop <- function(x) {
   dolz <- sum(x)
   x["TRUE"]/dolz
}

# Calculate euclidean distance between two species.
euclidDistance <- function(x) {
   sqrt(sum(sapply(x, diff)^2))
}

# Simulation function where we compare two species using euclidean distance against
# a 17-dimensional hypervolume of virtual species. Data for generating these species
# (that constitute a null model) is constructed using a uniform distribution from
# parameters of actual species.

calculateEuclideanForTwoSpecies <- function(input, MLP = 2, desc, N = 10000) {
   # input = raw data, where one column is species name and the rest are morphological
   #         characters (already transformed, e.g. log of body size)
   # MLP = to get 10000 simulations, multiply by this constant
   # desc = description of species, one column is location, one is species name and
   # other columns are ignored
   # N = number of simulated (pair-wise) values

   # calculate range of variable (character)
   character.range <- sapply(input[, -which(names(input) == "species")], range)

   # simulate N *MLP virtual species with 17 characters each
   vir.spec.two <- t(replicate(N * MLP, apply(character.range, MARGIN = 2,
            FUN = function(x) runif(n = 1, min = x[1], max = x[2]))))

   # construct a vector that will enable us to compare two and two species
   split.along <- rep(1:((N * MLP)/2), each = 2)

   sim.euc <- data.frame(vals = unclass(by(vir.spec.two, INDICES = split.along,
            FUN = function(x) sqrt(sum(sapply(x, diff)^2)))))

   # Find which localities have two species in a community.
   twos <- rle(as.character(desc$locality))$lengths == 2
   twos <- droplevels(unique(desc$locality)[twos])
   twos <- droplevels(desc[desc$locality %in% twos, ])
   twos <- split(twos, f = twos$locality)

   # Calculate Euclidean distance of real pairs.
   real.ed <- unlist(lapply(twos, FUN = function(x, data) {
            cte <- data[data$species %in% x$species, ]
            ed <- sqrt(sum(sapply(cte[, -1], diff)^2))
            names(ed) <- paste(cte$species, collapse = " - ")
            ed
         }, data = input))

   # Find what proportion of simulated values are bigger than values from real pairs.
   sim.bigger.real <- sapply(real.ed, FUN = function(x, sim.vals) {
         rslt <- prop(table(sim.vals > x))
         names(rslt) <- names(x)
         rslt
      }, sim.vals = unname(sim.euc$vals))

   # We round the numbers
   real.ed <- data.frame(vals = real.ed)

   lbls <- sapply(twos, FUN = function(x) sapply(x$species, FUN = gsub, pattern = "_", replacement = ". "))
   lbls <- apply(lbls, MARGIN = 2, FUN = function(x) paste(x, collapse = " - "))

   output <- list(sim.euc = sim.euc,
      real.euc = real.ed,
      pval = sim.bigger.real,
      labels = lbls)
   output
}

# Simulation function to compare three species using standard deviation. Null model
# is constructed as in calculateEuclideanForTwoSpecies().

calculateSDforThreeSpecies <- function(input, MLP, desc, N = 10000) {
   # input = raw data, where one column is species name and the rest are morphological
   # characters (already transformed, e.g. log of body size)
   # MLP = to get 10000 simulations, multiply by this constant
   # desc = description of species, one column is location, one is species name and
   # other columns are ignored

   # calculate range of variable (character)
   character.range <- sapply(input[, -which(names(input) == "species")], range)

   vir.spec.three <- t(replicate(N * MLP, apply(character.range, MARGIN = 2,
            FUN = function(x) runif(n = 1, min = x[1], max = x[2]))))
   split.along <- rep(1:((N * MLP)/3), each = 3)

   # Find which localities have three species in a community.
   threes <- rle(as.character(desc$locality))$lengths == 3
   threes <- droplevels(unique(desc$locality)[threes])
   threes <- droplevels(desc[desc$locality %in% threes, ])
   threes <- split(threes, f = threes$locality)

   # Function that calculates standard deviation of a triplet.
   calculateSD <- function(x) {
      combos <- matrix(c(1,2, 1,3, 2, 3), byrow = TRUE, ncol = 2)
      # for each combination, calculate euclidean distance
      get.sd <- apply(combos, MARGIN = 1, FUN = function(m, z) euclidDistance(z[m, ]), z = x)
      top3 <- which(get.sd == max(get.sd))
      sd(get.sd[-top3]) # calculate SD based on two shortest euclidean distances
   }

   # Calculate standard deviation from all combinations of euclidean distances
   # of three species in a community.
   sim.sd <- data.frame(vals = unclass(by(vir.spec.three, INDICES = split.along, FUN = calculateSD)))

   # Calculate sd of euclidean distances for real communities.
   real.sd <- data.frame(vals = unlist(lapply(threes, FUN = function(m, data) {
               to.clclt <- data[data$species %in% m$species, ]
               calculateSD(to.clclt[, -which(names(to.clclt) == "species")])
            }, data = input)))

   # Find what proportion of simulated values are bigger compared to real.
   sim.bigger.real <- apply(real.sd, MARGIN = 1, FUN = function(x, sim.vals) {
         rslt <- prop(table(sim.vals > x))
         names(rslt) <- names(x)
         rslt
      }, sim.vals = sim.sd$vals)

   lbls <- sapply(threes, FUN = function(x) sapply(x$species, FUN = gsub, pattern = "_", replacement = ". "))
   lbls <- apply(lbls, MARGIN = 2, FUN = function(x) paste(x, collapse = " - "))

   output <- list(sim.euc = sim.sd,
      real.euc = real.sd,
      pval = sim.bigger.real,
      labels = lbls)
}

# Draw the whole shebang.
drawFigure <- function(x, xlab = c("sd", "ed"), label = TRUE) {
   # x = object from a simulation function (calculateEuclideanForTwoSpecies
   # or calculateSDforThreeSpecies)
   # xlab = what to put as label for x axis

   if (xlab == "sd") {
      xlab <- "St. dev. of euclidean distances between virtual species"
   } else {
      xlab <- "Euclidean distances between virtual species"
   }

   gg.object <- ggplot(x[["sim.euc"]], aes(x = vals)) +
      theme_bw() +
      xlab(xlab) +
      theme(axis.title.x = element_text(size = 10)) +
      ylab("Frequency") +
      geom_histogram(fill = "light gray") +
      geom_vline(xintercept = x[["real.euc"]]$vals, size = 0.8)

   if (label) {
      # Calulate relative position of labels, thanks to Paul
      # http://stackoverflow.com/questions/7705345/how-can-i-extract-plot-axes-ranges-for-a-ggplot2-object
      get.counts <- ggplot_build(gg.object)
      position.l <- max(get.counts$panel$ranges[[1]]$y.range)
      position.l <- position.l - 0.5 * position.l
      # and construct labels
      lbls <- paste(x[["pval"]], x[["labels"]], sep = "    ")

      gg.object +
         geom_text(data = x[["real.euc"]], aes(x = vals), y = position.l, label = lbls, angle = -90,
            vjust = -0.7, size = 3)
   }
   gg.object

   # preberi koliko ima vrst združba
   # naredi test
}

#
actualDiffClusters <- function(loc, coms) {
    one.com <- droplevels(coms[coms$Locality == loc, "Species"])
    one.com <- spe[spe$Species %in% one.com, ]

    if (length(unique(one.com$group)) == 1) {
        message("Only one group present at locality")
        return(NA)
    }

    # for each group, find cluster mean...
    clust.means <- aggregate(. ~ group,
                             FUN = mean, na.rm = TRUE,
                             data = one.com[, !grepl("Species|ecomorph", names(one.com))])

    # ... and find euclidean difference between them
    clust.group.diff <- euclidDistance(clust.means[, !grepl("group", names(clust.means))])
    return(clust.group.diff)
}
# Calculate null distribution of a cluster given the number of species.
nullCluster2 <- function(x, spe) {
    # browser()
    ss <- spe[x, ]
    euclidDistance(ss[, !grepl("group|ecomorph|Species", names(ss))])
}

nullCluster3 <- function(x, spe) {
    mst(spe[x, !grepl("Species|ecomorph|gro", x = names(spe))])
}

mst <- function(x) {
    combos <- combn(3, 2)
    euc <- apply(combos, MARGIN = 2, FUN = function(m, x) {
        euclidDistance(x[m, ])
    }, x = x)
    sum(sort(euc)[1:2])
    }
