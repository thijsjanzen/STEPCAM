# this function calculates the sd values for the different traits, 
# in order to normalize the traits
calcSD <- function(species, abundances, n_plots, n_traits)  {

  a <- rowSums(abundances)
  if(length(which(a == 0)) > 0) {
    stop("calcSD: ",
         "One of your communities doesn't have species in it")
  }
  a <- rowSums(abundances)
  if(length(which(a < 3)) > 0) {
    stop("calcSD: ",
         "One of your communities doesn't have at least three species",
         " in it in order to calculate summary statistics")
  }
  optimum <- matrix(nrow = n_plots, ncol = n_traits)

  for (n in seq_len(n_plots)) {
    # select species present in focal community: abundance above 0
    esppres <- which(abundances[n, ] > 0)
    # traits and abundances of species present in community
    tr <- species[esppres,]
    # Community Trait Means (CTMs). 
    # These will be used as the 'optimal' trait value in a
    # community for the filtering model. 
    # Species very dissimilar from this value will be filtered out
    for (i in seq_len(n_traits)){
      optimum[n, i] <- mean(tr[, (i + 1)])
    }
  }

  row.names(abundances) <- c(1:n_plots)
  abundances2 <- as.data.frame(abundances)
  species2 <- species[,c(2:(n_traits + 1))]
  #species2 <- cbind(names(abundances2),species2)
  species2 <- as.matrix(species2)
  row.names(species2) <- names(abundances2)

  # calculate observed FD values
  Ord <- ordinationAxes(x = species2, stand.x = FALSE)
  res <- detMnbsp(Ord, abundances)
  FD_output <- strippedDbFd(Ord, abundances2, res[[1]], res[[2]])

  # calculate average CTM values across plots
  average_optimums <- c()
  for (i in seq_len(n_traits)){
    average_optimums[i] <- mean(optimum[, i])
  }
  optimums_plus_average <- rbind(optimum, average_optimums)

  # calculate trait distances from average optima
  mean_multi_trait_difference <- mean( as.matrix( 
                                 dist( optimums_plus_average))
                                 [n_plots + 1, c(1:n_plots)])

  # calculate SD of FD values
  sdFRic <- sd(FD_output$FRic)
  sdFEve <- sd(FD_output$FEve)
  sdFDiv <- sd(FD_output$FDiv)

  # standard deviations of observed FD and trait mean values in plots
  sd_values <- c(sdFRic, sdFEve, sdFDiv, mean_multi_trait_difference)

  return(sd_values)
}

# calculate frequencies (number of plots in which it occurs) of all species
generateFrequencies <- function(abundances) {
  if(length(which(is.na(abundances) == TRUE)) > 0) {
    stop("generateFrequencies: ",
         "One or more entries in the abundance matrix is NA")
  }
  if(min(abundances) < 0) {
    stop("generateFrequencies: ",
         "One or more entries in the abundance matrix is negative")
  }

  # number of plots
  samples <- length(abundances[, 1])

  # total number of species in species pool
  taxa <- length(abundances[1, ])

  # new (empty) presence matrix: species x plot matrix 
  # with only 1's (present) and 0's (absent)
  presences <- matrix(nrow = samples, ncol = taxa)

  # vector with frequency of each species
  frequencies <- c();

  for(i in seq_len(samples)){
    for(j in seq_len(taxa)){
      ifelse(abundances[i, j] > 0 , 
             presences[i, j] <- 1 , 
             presences[i, j] <- 0)
    }
  }
  for (i in seq_len(taxa)){
    frequencies[i] <- sum(presences[, i])
  }
  return(frequencies)
}

# function to transform species trait variables to a standard 
# normal distribution (mean = 0, sd = 1)
scaleSpeciesvalues <- function(species, n_traits) {
  #check dimensions
  if (length(species[1,]) != (n_traits+1)) {
    stop("scaleSpeciesvalues: ",
         "incorrect trait matrix dimensions ", "\n",
         " did you perhaps remove the species names?")
  }

  # scale the trait values
  standard_deviations <- c();
  means <- c()
  for (i in 2:(n_traits + 1)) {
    standard_deviations[i] <- sd(species[, i])
    means[i] <- mean(species[, i])
  }

  if (length(which(standard_deviations == 0)) > 0) {
    stop("scaleSpeciesvalues: ",
         "one of your traits has no variation ", "\n",
         "  most likely trait(s): ", -1 + which(standard_deviations==0))
  }

  for (i in 2:(n_traits + 1)) {
    species[, i] <- (species[, i] - means[i]) / standard_deviations[i]
  }
  return(species)
}

# function to generate simulated FD / trait mean values given 
# a certain species/trait pool and certain community assembly 
# parameter settings
generateValues <- function(params, species, 
                           abundances, community_number, 
                           n_traits) {
  
  
  if (n_traits == 1) {
    stop("generateValues: ",
         "need more than 1 trait")
  }
  
  # calculate for each species in how many plots it occurs
  samples <- length(abundances[, 1])

  # total number of species in species pool
  taxa <- length(abundances[1, ])

  # put all info about species in one data frame
  traitnames <- names(species[-1])
  species <- as.data.frame(cbind(species))
  names(species) <- c("sp", traitnames)
  row.names(species) <- c(1:taxa)
  species[,1] <- c(1:taxa)

  # select species present in focal community: abundance above 0
  esppres <- which(abundances[community_number, ] > 0)

  #  number of species in the community
  S <- length(esppres)

  # total species falling out:
  species_fallout <- taxa - S

  # species that fall out through stochasticity
  dispersal_fallout <- round(params[1] * species_fallout, 0)

  # species that fall out through filtering
  filtering_fallout <- round((params[2] / (params[2] + params[3])) *
                               (species_fallout - dispersal_fallout),0)

  # species that fall out through limiting similarity
  competition_fallout <- species_fallout - 
                         dispersal_fallout - 
                         filtering_fallout

  params2 <- c(dispersal_fallout, filtering_fallout, competition_fallout)

  # the Kraft.generator function: function that runs hybrid Kraft models
  allfinaloutput <- STEPCAM(params2, species, abundances, taxa, esppres,
                            community_number, n_traits, species_fallout)

  traits <- as.data.frame(species[, c(2:(n_traits + 1))],
                          row.names = seq_along(allfinaloutput))
  presences <- as.data.frame(t(allfinaloutput), seq_along(allfinaloutput))
  names(presences) <- seq_len(allfinaloutput)

  esppres <- which(presences > 0)
  traits <- traits[esppres, ]
  presences <- presences[ ,esppres]

  # now make a new species x trait matrix, 
  # with only an n_traits (see settings) 
  # amount of PCA traits
  "+" <- function(...) UseMethod("+")
  "+.default" <- .Primitive("+")
  "+.character" <- function(...) paste(..., sep = "")
  traitnames <- c()
  for (i in 1:n_traits) {
    traitnames[i] <- "PC" + i
  }
  names(traits) <- traitnames

  Ord <- ordinationAxes(x = traits, stand.x = FALSE)
  res <- detMnbsp(Ord, presences)
  FD_output <- strippedDbFd(Ord, presences, res[[1]], res[[2]])

  trait_means <- c()
  for (i in 1:n_traits){
    trait_means[i] <- mean(traits[, i])
  }

  summary_stats <- cbind(FD_output$FRic, 
                         FD_output$FEve, 
                         FD_output$FDiv, 
                         t(trait_means))

  return(summary_stats)
}


