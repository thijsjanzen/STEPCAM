# function to generate a random combination of dispersal, filtering and competition parameter settings (uninformed prior assumed)
getRandomVals <- function(max_val)  {
  x <- stats::runif(3, min = 0, max = 1)
  x <- x / sum(x); #normalize to 1
  x <- x * max_val; #translate to integers
  x2 <- floor(x); #round to integers
  while (sum(x2) != max_val) {
    a <- sample(1:3, size = 1, prob = x - floor(x) )
    x2[a] <- x2[a] + 1;
  }
  return(x2)
}

# function to randomly draw a particle depending on it's weight
getFromPrevious <- function(inds, ws, disps, filts, comps)  {
  index <- sample(x = inds, size = 1, replace = TRUE, prob = ws)
  output <- c(disps[index], filts[index], comps[index])
  return(output);
}

# function to calculate the weight of a particle
calculateWeight <- function(params, target, sigma,
                            disp_vals, filt_vals, comp_vals,
                            weights)  {
  sum <- 0.0
  vals <- c()

  for (i in seq_along(disp_vals)) {
    previous_particles <- c(disp_vals[i], filt_vals[i], comp_vals[i])
    diff <- params[target] - previous_particles[target]
    vals[i] <- weights[i] * dnorm(diff, mean = 0, sd = sigma)
  }
  return(1 / sum(vals));
}

# normalize all the weights of all particles such that they sum to 1
normalizeWeights <- function(x) {
  sum_x <- sum(x)
  x <- x / sum_x
  return(x)
}

# function to randomly change the contribution of one of the processes:
perturb <- function(p, sigma)  {
  params <- p;
  max_number <- sum(p)
  numbers <- 1:3

  x <- sample(numbers, 3, replace = FALSE)

  oldval <- params[x[1]]

  params[x[1]] <- round(params [x[1]] + rnorm(1, mean = 0, sd = sigma), 0)
  params[x[1]] <- max(0, params[x[1]])
  params[x[1]] <- min(max_number, params[x[1]]);

  diff <- params [x[1]] - oldval

  params[x[2]] <- params[x[2]] - diff
  params[x[2]] <- max(0, params[x[2]])
  params[x[2]] <- min(max_number, params[x[2]])

  params[x[3]] <- max_number - (params[x[1]] + params[x[2]])
  params[x[3]] <- max(0, params[x[3]])
  params[x[3]] <- min(max_number, params[x[3]])

  return(c(params, x[1]))
}

# function to calculate the fit
calculateDistance <- function(rich, even, div, opt_diff, obs, sd_vals)  {
  fit_rich <- (abs( (rich - obs[, 1]) ) / sd_vals[1]) ^ 2;
  fit_even <- (abs( (even - obs[, 2]) ) / sd_vals[2]) ^ 2;
  fit_div  <- (abs( (div -  obs[, 3]) ) / sd_vals[3]) ^ 2;
  fit_optima <- (opt_diff / sd_vals[4])^2;

  full_fit <- fit_rich + fit_even + fit_div + fit_optima;

  return(full_fit)
}


ABC_SMC <- function(numParticles, species_fallout, taxa, esppres, n_traits,
                    sd_vals, summary_stats, community_number, species,
                    abundances, frequencies, stopRate, Ord, continue_from_file = TRUE,
                    stop_at_iteration = 50)  {

  for(i in seq_along(sd_vals)) {
    if(sd_vals[[i]] == 0.000) {
      stop("ABC_SMC: ",
           "one of the community summary statistics shows no variation in your dataset")
    }
  }

  optimum <- summary_stats[, 4:(3 + n_traits)]

  disp_vals <- 1:numParticles
  filt_vals <- 1:numParticles
  comp_vals <- 1:numParticles

  fits <- 1:numParticles
  rich_vec <- 1:numParticles
  eve_vec <- 1:numParticles
  div_vec <- 1:numParticles
  opt_vec <- 1:numParticles

  next_disp <- disp_vals
  next_filt <- filt_vals
  next_comp <- comp_vals

  weights <- rep(1, numParticles)
  next_weights <- rep(1, numParticles)
  indices <- 1:numParticles

  sigma <- 1
  t <- 1

  f <- list.files(pattern = "particles_t=")
  if (length(f) > 0 && continue_from_file == TRUE) {
    cat("Found previous output, continuing from that output\n"); flush.console();
    f <- gtools::mixedsort(f)
    t1 <- 1 + length(f)
    d <- read.table(f[length(f)], header = F)
    if(d[numParticles,1] == numParticles) {
      d <- read.table(f[length(f) - 1], header = F)
      t1 <- t1 - 1
    }

    disp_vals <- d[, 1]
    filt_vals <- d[, 2]
    comp_vals <- d[, 3]
    fits <-     d[, 8]
    weights <-  d[, 9]

    t <- t1
  }

  # continuously sampling
  while (t < 50)  {
    cat("\nGenerating Particles for iteration\t", t, "\n")
    cat("0--------25--------50--------75--------100\n")
    cat("*"); flush.console()
    PRINT_FREQ <- 20

    numberAccepted <- 0
    if (t != 1) weights <- normalizeWeights(weights)

    threshold <- 200 * exp(-0.5 * t)

    stop_iteration <- 0
    changed <- 1
    tried <- 1

    while (numberAccepted < numParticles) {
      params <- c(species_fallout, 0, 0)
      # get a parameter combination
      if (t == 1)  {
        params <- getRandomVals(species_fallout)
      } else {
        params <- getFromPrevious(indices, weights,
                                  disp_vals, filt_vals, comp_vals)
        params <- perturb(params, sigma)

        # we need to know which parameter was perturbed,
        # to be able to calculate its weight later
        changed <- params[4]
        params <- params[1:3]

        if ( sum(params) > species_fallout) {
          stop("ABC_SMC: ",
               "too much params after perturb!");
        }
      }

      # total number of species in species pool
      taxa <- length(abundances[1, ])
      allcommunities <- matrix(nrow = taxa)
      allcommunities <- STEPCAM(params, species, abundances, taxa,
      esppres, community_number, n_traits, species_fallout)

      # make traits (total species richness x number of traits)
      # and community ((66 * permutations * communities) x
      # total species richness) matrices
      traits <- as.data.frame(species[, c(2:(n_traits + 1))],
                              row.names = c(1:taxa))
      communities <- as.data.frame(t(allcommunities), names = (1:taxa))
      names(communities) <- c(1:taxa)
      present_species <- as.vector(which(colSums(communities) > 0))

      # calculate several measures of FD of modeled community
      FD_output <- strippedDbFd(Ord, communities[, present_species])
      FRic <- FD_output$FRic # FRic = functional richness (Villeger et al, 2008, Ecology)
      FEve <- FD_output$FEve # FEve = functional evenness (Villeger et al, 2008, Ecology)
      FDiv <- FD_output$FDiv # FDiv = functional diversity (Villeger et al, 2008, Ecology)

      trait_means <- c()
      for (i in 1:n_traits) {
        # trait means of simulated community
        trait_means[i] <- mean(traits[present_species, i])
      }
      optimum_plus_trait_means <- rbind(optimum, trait_means)

      # calculate distance of trait mean simulated community from that of observed
      mean_optimum <- dist(optimum_plus_trait_means)

      # (inverse) fit of model: euclidian distance of FD and trait mean values of
      # observed community from that of simulated
      fit <- calculateDistance(FRic[[1]], FEve[[1]], FDiv[[1]],
                               mean_optimum[1], summary_stats, sd_vals)

      # function to accept / reject models based on the fit
      if (fit < threshold) {
        numberAccepted <- numberAccepted + 1
        next_disp[numberAccepted] <- params[1]
        next_filt[numberAccepted] <- params[2]
        next_comp[numberAccepted] <- params[3]

        fits[numberAccepted] <- fit
        rich_vec[numberAccepted] <- FRic[[1]]
        eve_vec[numberAccepted] <- FEve[[1]]
        div_vec[numberAccepted] <- FDiv[[1]]
        opt_vec[numberAccepted] <- mean_optimum;

        if (t == 1) {
          next_weights[numberAccepted] = 1
        } else {
          next_weights[numberAccepted] =
            calculateWeight(params, changed, sigma, disp_vals, filt_vals,
                            comp_vals, weights)
        }

        if ((numberAccepted) %% (numParticles / PRINT_FREQ) == 0) {
          cat("**") ; flush.console()
        }
      }

      tried <- tried + 1
      if (tried > (1/stopRate) && tried > 50)  {
        # do not check every particle if the acceptance rate is OK
        if (numberAccepted / tried < stopRate) {
          stop_iteration <- 1; break;
        }
      }

      if (t >= stop_at_iteration) {
        stop_iteration <- 1; break;
      }
    }

    # replace values
    disp_vals <- next_disp
    filt_vals <- next_filt
    comp_vals <- next_comp
    weights <- next_weights

    output <- cbind(disp_vals, filt_vals, comp_vals, rich_vec,
                    eve_vec, div_vec, opt_vec, fits, weights)
    file_name <- paste("particles_t=", t, ".txt", sep="", collapse = NULL)
    write.table(output, file_name, row.names = F, col.names = F)

    cat(" ", mean(disp_vals), mean(filt_vals), mean(comp_vals), "\t", "accept rate = ",
        numberAccepted / (tried-1), "\n")

    # and reset
    next_weights <- rep(1,numParticles)
    next_disp <- 1:numParticles
    next_filt <- 1:numParticles
    next_comp <- 1:numParticles

    if (stop_iteration == 1){
      break
    }
    t <- t + 1
  }


  if (t >= 1) {
    d <- read.table(paste("particles_t=", t - 1, ".txt", sep="",
                          collapse = NULL), header = F)
  } else {
      d <- read.table(paste("particles_t=", t, ".txt", sep="",
	                          collapse = NULL), header = F)
  }
  output <- list( DA = d[, 1], HF = d[, 2], LS = d[, 3])
  return(output)
}
