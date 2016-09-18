context("TernPlot")

test_that("TernPlot: use", {
  set.seed(42)
  Artificial.Data <- generate.Artificial.Data(n_species = 40, n_traits = 3,
                                              n_communities = 5, occurence_distribution = 0.2,
                                              average_richness = 0.5, sd_richness = 0.2,
                                              mechanism_random = TRUE)
  O <- STEPCAM_ABC(Artificial.Data$abundances, Artificial.Data$traits,
                   numParticles = 100, n_traits = 3, plot_number = 1, stopRate = 0.008) ;
  TernPlot(O);
})
