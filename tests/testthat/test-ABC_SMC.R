context("ABC_SMC")

test_that("ABC_SMC: use", {
  skip("Work in Progress")

  #source("/Users/janzen/GitHub/STEPCAM/R/Preparation.R")

  set.seed(42)
  n_traits <- 2
  n_plots <- 10
  x <- generate.Artificial.Data(n_species = 10, n_traits = n_traits,
                                n_communities = ncomm,
                                occurence_distribution = 0.5,
                                average_richness = 1,
                                sd_richness = 1,
                                mechanism_random = FALSE)

  data_species <- x$traits;
  data_species$trait1 <- c(1,1.1,1.2,1.3,5,7.6,7.7,7.8,7.9,8)
  data_species$trait2 <- c(1,1.1,1.2,1.3,5,7.6,7.7,7.8,7.9,8)
  taxa <- nrow(data_species);

  data_abundances <- x$abundances;

  numParticles <- 100
  plot_number <- 1
  stopRate <- 0.1


  params <- c(0,1,0)
  scaled_species  <- scaleSpeciesvalues(data_species,n_traits)
  abundances <- data_abundances
  row.names(abundances) <- c(1:n_plots)
  data_frequencies <- generateFrequencies(data_abundances)

  scaled_species <- as.data.frame(cbind(scaled_species[, c(1:(n_traits + 1))], data_frequencies)) ;
  traitnames <- names(data_species)[-1]
  names(scaled_species) <- c("sp", traitnames[1:n_traits], "freq")
  row.names(scaled_species) <- c(1:taxa)

  summary_stats <- generateValues(params, scaled_species, data_abundances,
                                   plot_number, n_traits)

  esppres <- which(data_abundances[plot_number, ] > 0) ;
  S <- length(esppres);
  species_fallout <- taxa - S
  Ord <- ordinationAxes(x = scaled_species[,-1], stand.x = FALSE)
  sd_vals <- calcSD(scaled_species, data_abundances, n_plots, n_traits);

  species <- scaled_species
  abundances <- data_abundances






  scaled_species <- as.data.frame(cbind(scaled_species[, c(1:(n_traits + 1))], data_frequencies)) ;
  traitnames <- names(data_species)[-1]
  names(scaled_species) <- c("sp", traitnames[1:n_traits], "freq")
  row.names(scaled_species) <- c(1:taxa)
  stopRate <- 0.1

  output <- ABC_SMC(numParticles, species_fallout, taxa, esppres, n_traits,
                                sd_vals, summary_stats, community_number, scaled_species,
                                data_abundances, data_frequencies, stopRate, Ord)


})


