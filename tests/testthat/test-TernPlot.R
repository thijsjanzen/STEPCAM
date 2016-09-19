context("TernPlot")

test_that("TernPlot: use", {
  set.seed(42)
  Artificial.Data <- generate.Artificial.Data(n_species = 10, n_traits = 3,
                                              n_communities = 5, occurence_distribution = 0.2,
                                              average_richness = 0.5, sd_richness = 0.2,
                                              mechanism_random = TRUE)
  O <- STEPCAM_ABC(Artificial.Data$abundances, Artificial.Data$traits,
                   numParticles = 100, n_traits = 3, plot_number = 1, stopRate = 0.008,
                   stop_at_iteration = 7, continue_from_file = TRUE) ;
  TernPlot(O);
  d <- cbind(O$DA, O$HF, O$LS);
  ternaryplot2(d, scale=1, col="black", grid=T, cex = 0.5, labels = c("inside"),
               dimnames = c("DA", "HF", "LS"), main="", coordinates = TRUE)
  ternaryplot2(d, scale=1, col="black", grid=T, cex = 0.5, labels = c("outside"),
               dimnames = c("DA", "HF", "LS"), main="", coordinates = TRUE,
               dimnames_position = "edge")
  ternaryplot2(d, scale=1, col="black", grid=T, cex = 0.5, labels = c("outside"),
               main="", coordinates = TRUE,
               dimnames_position = "edge")
})

test_that("TernPlot: abuse:", {
  Artificial.Data <- generate.Artificial.Data(n_species = 10, n_traits = 3,
                                              n_communities = 5, occurence_distribution = 0.2,
                                              average_richness = 0.5, sd_richness = 0.2,
                                              mechanism_random = TRUE)
  O <- STEPCAM_ABC(Artificial.Data$abundances, Artificial.Data$traits,
                   numParticles = 100, n_traits = 3, plot_number = 1, stopRate = 0.008,
                   stop_at_iteration = 7, continue_from_file = TRUE)

  d <- cbind(O$DA, O$HF);

  expect_error(
  ternaryplot2(d, scale=1, col="black", grid=T, cex = 0.5, labels = c("inside"),
               dimnames = c("DA", "HF", "LS"), main="", coordinates = TRUE),
  "Need a matrix with 3 columns")

  d <- cbind(O$DA, O$HF, O$LS);
  d[1,1] <- -1
  expect_error(
    ternaryplot2(d, scale=1, col="black", grid=T, cex = 0.5, labels = c("inside"),
                 dimnames = c("DA", "HF", "LS"), main="", coordinates = TRUE),
    "X must be non-negative")
  d[1,2] <- -1
  d[1,3] <- 1
})
