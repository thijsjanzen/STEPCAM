context("calcSD")

test_that("calcSD: use", {
  set.seed(42)
  n_traits <- 3
  n_plots <- 10
  num_species <- 10
  x <- generate.Artificial.Data(numSpecies = num_species, numTraits = n_traits,
                                numCommunities = n_plots,
                                occurence_distribution = 0.5,
                                average_richness = 10,
                                SD_richness = 1,
                                random.Mechanism = FALSE)

  data_species <- x$traits
  data_species$trait1 <- 1:10
  data_species$trait2 <- 1:10
  data_species$trait3 <- 1:10

  scaled_species <- scaleSpeciesvalues(data_species, n_traits)

  x <- calcSD(scaled_species,data_abundances, n_plots, n_traits)

  expect_equal(
    x[[1]] ,
    0
  )

  expect_equal(
    x[[2]],
    0
  )

  expect_equal(
    x[[4]],
    0
  )
})
