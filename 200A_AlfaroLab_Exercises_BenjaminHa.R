##### EE BIOL 200A - Alfaro Lab Exercises #####
### Benjamin Ha ~ October 5, 2017 ###

setwd("~/Box Sync/200A/Alfaro Labs/Exercises")
getwd()

# Three separate packages required before using phytools library
install.packages("ape")
install.packages("maps")
install.packages("phytools")
library(phytools)
install.packages("geiger")
library(geiger)

### Exercise 1
# Calculate the gamma statistic for this phylogeny of homalopsid snakes from Alfaro et al., 2008: snake.tre.
snake.tree <- read.tree("homalops.phy")
obj <- ltt(snake.tree,
           log.lineages = FALSE)
snake.gamma <- obj$gamma # save only the gamma value into an object
snake.gamma

### Exercise 2
# Given this gamma value, what would you conclude about the tempo of speciation in this clade?
# Since gamma statistic = -3.2411, I would conclude the phylogenetic tree may have nodes closer to the root than the tips. This means the speciation rate was initally high (hence, closer to the root), but slowed down over time (hence, further away from the tips).


### Exercise 3
# Given that the crown age of the snake radiation is 22 MY and the total richness of the clade is 34 species, determine whether the observed gamma could be due to the amount of incomplete sampling in the empirical tree. On the basis of the MCCR test what can you conclude about the tempo of homalopsid snake diversification?

# input data based on exercise as the starting point. This is observed data
age <- 22
richness <- 34
snake.birth = (log(richness) - log(2)) / age 
snake.birth
# birth rate of new species over time as opposed to an individual species

# now simulate gamma values when trees are undersampled
richness <- 34
missing <- 12
num_simulations <- 200

gamma_null <- numeric(num_simulations) #gamma_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(snake.birth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->gamma_null[i]
} # long loop, so let it run on its own

# create a histogram of the simulated birth rate
par(mfrow=c(1,1), mar=c(4, 4, 3, 2))
hist(gamma_null) 

# now have arrove to show where the gamma value lies on the simulated tree
arrows(snake.gamma, 40, snake.gamma, 0, col="red", lwd=2)
snake.gamma

# Conclusion: The gamma value is -3.241081. Due to incomplete sampling, the snake.gamma value suggests the birth rate is skewed in relation to the simulated data. Again, this suggests the speciation rate was initially high, and then slowed down over time.

### Exercise 4
# What is the birth rate and death rate of the homalopsid tree
library(phytools)
fitbd <- birthdeath(snake.tree)
fitbd
bd(fitbd)
# Birth rate is 0.06839495. Death rate is zero. The death rate is definitely skeptical, but I assume that since the snake.tree has incomplete sampling that it would skew the birth-death rates since birth-death rates are conditioned to full sampling of a tree.


### Exercise 5
# Find a time-calibrated phylogeny for a group that interests you (ideally with more than 30 tips and fewer than 200). Do the following: 

# (1) Describe the clade (including a description of the number of tips in the tree and the total number of species in the clades) and provide a reference or citation to the source.
fish.tree <- read.tree("etheostoma_percina_chrono.tre")
# This phylogenetic tree shows the relationships of darters, a species-rich clade of North American freshwater fishes. The data is a near-complete taxon sampling of 245 out of 248 species.
# Source: Near et al. (2011) Phylogeny and temporal diversification of darters (Percidae: Etheostomatinae). Systematic Biology.

# (2) Fit a birthdeath model to this tree and report b and d. 
fitbd.fish <- birthdeath(fish.tree) # code originally not working because is.binary = FALSE
is.binary(fish.tree)
fish.tree <- multi2di(fish.tree)
fitbd
bd(fitbd) # death rate is also zero???

# (3) Perform an MCCR test and describe whether the gamma value is extreme or not given the level of sampling in the tree.
obj <- ltt(fish.tree,
           log.lineages = FALSE)
obj
fish.gamma <- obj$gamma # save only the gamma value into an object
fish.gamma

age <- 20 # arbitrary number that you selected
richness <- 248 # total number of species
fish.birth = (log(richness) - log(2)) / age 

missing.fish <- 3 # taken from fact that 245 out of 248 species were sequenced
num_simulations <- 200

gamma.null.fish <- numeric(num_simulations)
for(i in 1:num_simulations) {
  sim.bdtree(fish.birth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->gamma.null.fish[i]
}

# create a histogram of the simulated birth rate
hist(gamma.null.fish) 
arrows(fish.gamma, 40, fish.gamma, 0, col="forest green", lwd=3)

# Conclusion: The gamma value is 0.2007 and when comparing the gamma value to the simulated data, the gamma value is not too extreme. It is still somewhat skewed (positive), which suggests the speciation rate was slow initially and then increased  later; however, since the gamma value is relatively close to zero, I would anticipate the speciation rate did not increase significantly over time.