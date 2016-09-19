context("ordinationAxes")

test_that("ordinationAxes: use", {
  set.seed(42)
  n_traits <- 3
  n_plots <- 10
  num_species <- 10;
  x <- generate.Artificial.Data(n_species = num_species, n_traits = n_traits,
                                n_communities = n_plots,
                                occurence_distribution = 0.5,
                                average_richness = 10,
                                sd_richness = 1,
                                mechanism_random = TRUE)

  data_species <- x$traits
 # data_species$trait1 <- 1:10
#  data_species$trait2 <- 1:10
#  data_species$trait3 <- 1:10

  data_abundances <- x$abundances

  species  <- scaleSpeciesvalues(data_species,n_traits)
  abundances <- data_abundances


  row.names(abundances) <- c(1:n_plots)
  abundances2 <- as.data.frame(abundances)
  species2 <- species[,c(2:(n_traits + 1))] ;
  #species2 <- cbind(names(abundances2),species2)
  species2 <- as.matrix(species2)
  row.names(species2) <- names(abundances2)

  # calculate observed FD values
  Ord <- ordinationAxes(x = species2, stand.x = FALSE)
  Ord <- ordinationAxes(x = species2, stand.x = TRUE)

  v2 <- species2[1,]
  v2[1:3] <- runif(3,0,1)
  Ord <- ordinationAxes(x = v2, stand.x = FALSE)

  v2[1:3] <- as.factor(c(1,1,2))
  v2 <- as.dist(v2)
  Ord <- ordinationAxes(x = v2, stand.x = FALSE)
})
