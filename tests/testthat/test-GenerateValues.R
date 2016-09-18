context("generateValues")

test_that("generateValues: use", {
  skip("Work in Progress")
  #generateValues <- function(params, species, abundances, community_number, n_traits)
  params <- c(0,1,0)
  set.seed(42)
  n_traits <- 2
  x <- generate.Artificial.Data(numSpecies = 10, numTraits = n_traits,
                                numCommunities = 10,
                                occurence_distribution = 0.5,
                                average_richness = 1,
                                SD_richness = 1,
                                random.Mechanism = FALSE)

  data_species <- x$traits
  data_species$trait1 <- c(1,1.1,3,4,5,7,8,9,10,11)
  data_species$trait2 <- c(1,1.1,3,4,5,7,8,9,10,11)

  data_abundances <- x$abundances
  community_number <- 1
  scaled_species <- scaleSpeciesvalues(data_species,n_traits)
  scaled_species <- as.data.frame(cbind(scaled_species[, c(1:(n_traits + 1))], data_frequencies)) ;
  traitnames <- names(data_species)[-1]
  names(scaled_species) <- c("sp", traitnames[1:n_traits], "freq")
  row.names(scaled_species) <- c(1:taxa)

  v <- generateValues(params, scaled_species, data_abundances,
                      community_number, n_traits)


  abundances <- data_abundances
  species <-data_species

})
