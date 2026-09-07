# function to generate a random combination of dispersal, 
# filtering and competition parameter settings (uninformed prior assumed)
getRandomVals <- function(max_val)  {
  x <- stats::runif(3, min = 0, max = 1)
  x <- x / sum(x) #normalize to 1
  x <- x * max_val #translate to integers
  x2 <- floor(x) #round to integers
  while (sum(x2) != max_val) {
    a <- sample(1:3, size = 1, prob = x - floor(x) )
    x2[a] <- x2[a] + 1
  }
  x2 <- c(x2, 1)
  return(x2)
}

# function to randomly draw a particle depending on it's weight
getFromPrevious <- function(ws, disps, filts, comps, orders)  {
  index <- sample(x = seq_along(ws), size = 1, replace = TRUE, prob = ws)
  output <- c(disps[index], filts[index], comps[index], orders[index])
  return(output)
}

# function to calculate the weight of a particle
calculateWeight <- function(params, target, sigma,
                            disp_vals, filt_vals, comp_vals, order_vals,
                            weights)  {
  diff <- c() 
  if (target == 1) diff <- params[target] - disp_vals
  if (target == 2) diff <- params[target] - filt_vals
  if (target == 3) diff <- params[target] - comp_vals
  
  diff_prob <- dnorm(diff, mean = 0, sd = sigma)
  
  # we have to multiply with the ordering as well
  diff_order <- 1 - (params[4] != order_vals)
  # 90 % prob of remaining the same
  diff_order <- diff_order * 0.9
  diff_order[diff_order == 0] <- 0.1

  #vals <- cbind(log(weights), diff_prob, log(diff_order))
  #vals <- rowSums(vals)
  #vals <- exp(vals)
  
  a1 <- length(weights)
  a2 <- length(diff_prob)
  a3 <- length(diff_order)
  if (a1 != a2 || a1 != a3 || a2 != a3) {
    cat("oh oh")
    b <- 78
  }
  
  
  vals <- weights * diff_prob * diff_order 
  
  return( 1 / sum(vals) ) # prior density is 1.
}

# normalize all the weights of all particles such that they sum to 1
normalizeWeights <- function(x) {
  sum_x <- sum(x)
  x <- x / sum_x
  return(x)
}

# function to randomly change the contribution of one of the processes:
perturb <- function(p, sigma, fit_order)  {
  params <- p
  max_number <- sum(p[1:3])
  numbers <- 1:3

  x <- sample(numbers, 3, replace = FALSE)

  oldval <- params[x[1]]

  params[x[1]] <- round(params[x[1]] + rnorm(1, mean = 0, sd = sigma), 0)
  params[x[1]] <- max(0, params[x[1]])
  params[x[1]] <- min(max_number, params[x[1]])

  diff <- params[x[1]] - oldval

  params[x[2]] <- params[x[2]] - diff
  params[x[2]] <- max(0, params[x[2]])
  params[x[2]] <- min(max_number, params[x[2]])

  params[x[3]] <- max_number - (params[x[1]] + params[x[2]])
  params[x[3]] <- max(0, params[x[3]])
  params[x[3]] <- min(max_number, params[x[3]])

  if (fit_order == -1) {
    if (stats::runif(1, 0, 1) < 0.1) {
      new_order <- 1:6
      new_order <- new_order[-params[4]]
      params[4] <- sample(new_order, 1)
    }
  }
  
  return(c(params, x[1]))
}

# function to calculate the fit
calculateDistance <- function(rich, even, div, opt_diff, obs, sd_vals)  {
  fit_rich <- (abs( (rich - obs[, 1]) ) / sd_vals[1]) ^ 2
  fit_even <- (abs( (even - obs[, 2]) ) / sd_vals[2]) ^ 2
  fit_div  <- (abs( (div -  obs[, 3]) ) / sd_vals[3]) ^ 2
  fit_optima <- (opt_diff / sd_vals[4])^2

  full_fit <- fit_rich + fit_even + fit_div + fit_optima

  return(full_fit)
}

get_fit <- function(params, species, abundances, taxa,
                    esppres, community_number, n_traits,
                    species_fallout, fit_order, Ord, res, optimum,
                    summary_stats, sd_vals) {
  taxa <- length(abundances[1, ])
  
  allcommunities <- STEPCAM(params, species, abundances, taxa,
                            esppres, community_number, n_traits,
                            species_fallout,
                            fit_order)
  traits <- as.data.frame(species[, c(2:(n_traits + 1))],
                          row.names = c(1:taxa))
  
  communities <- as.data.frame(t(allcommunities))
  present_species <- which(communities > 0)
  
  FD_output <- strippedDbFd(Ord, communities, 
                            m = res[[1]], nb.sp = res[[2]]) 
  
  # FRic = functional richness (Villeger et al, 2008, Ecology)
  FRic <- FD_output$FRic 
  # FEve = functional evenness (Villeger et al, 2008, Ecology)
  FEve <- FD_output$FEve 
  # FDiv = functional diversity (Villeger et al, 2008, Ecology)
  FDiv <- FD_output$FDiv 
  
  trait_means <- vector("numeric", n_traits)
  for (i in seq_len(n_traits)) {
    # trait means of simulated community
    trait_means[i] <- mean(traits[present_species, i])
  }
  optimum_plus_trait_means <- rbind(optimum, trait_means)
  
  # calculate distance of trait mean between simulated community 
  # and observed community
  mean_optimum <- dist(optimum_plus_trait_means)
  
  # (inverse) fit of model: euclidian distance of FD and trait mean 
  # values of observed community from that of simulated
  fit <- calculateDistance(FRic[[1]], FEve[[1]], FDiv[[1]],
                           mean_optimum[1], summary_stats, sd_vals)
  
  return(list(fit = fit,
              FRic = FRic,
              FEve = FEve,
              FDiv = FDiv,
              mean_optimum = mean_optimum,
              params = params))
}


ABC_SMC <- function(numParticles, species_fallout, taxa, esppres, n_traits,
                    sd_vals, summary_stats, community_number, species,
                    abundances, frequencies, stopRate, Ord, 
                    continue_from_file = TRUE, stop_at_iteration = 50,
                    fit_order = FALSE,
                    num_threads = 1)  {

  for (i in seq_along(sd_vals)) {
    if (sd_vals[[i]] == 0.000) {
      stop("ABC_SMC: ",
           "one of the community summary statistics 
            shows no variation in your dataset")
    }
  }
  res <- detMnbsp(Ord, abundances)
  optimum <- summary_stats[, 4:(3 + n_traits)]

  disp_vals  <- c()
  filt_vals  <- c()
  comp_vals  <- c()
  order_vals <- c()

  fits <- c()
  rich_vec <- c()
  eve_vec <- c()
  div_vec <- c()
  opt_vec <- c()

  next_disp <- disp_vals
  next_filt <- filt_vals
  next_comp <- comp_vals
  next_order <- order_vals

  weights <- c()
  next_weights <- c()

  sigma <- 1
  t <- 1

  f <- list.files(pattern = "particles_t=")
  if (length(f) > 0 && continue_from_file == TRUE) {
    cat("Found previous output, continuing from that output\n")
    flush.console()
    f <- gtools::mixedsort(f)
    t1 <- 1 + length(f)
    d <- read.table(f[length(f)], header = FALSE)
    if (d[numParticles,1] == numParticles) {
      d <- read.table(f[length(f) - 1], header = FALSE)
      t1 <- t1 - 1
    }

    disp_vals <- d[, 1]
    filt_vals <- d[, 2]
    comp_vals <- d[, 3]
    fits <-     d[, 8]
    weights <-  d[, 9]
    order_vals <- d[, 10]

    t <- t1
  }

  # continuously sampling
  while (t < 50)  {
    cat("\nGenerating Particles for iteration\t", t, "\n")
    cat("0--------25--------50--------75--------100\n")
    cat("*")
    flush.console()
    PRINT_FREQ <- 20

    numberAccepted <- 0
    if (t != 1) weights <- normalizeWeights(weights)

    threshold <- 200 * exp(-0.5 * t)

    stop_iteration <- 0
    changed <- 1
    tried <- 0

    while (numberAccepted <= (numParticles)) {
      
      remaining <- numParticles - numberAccepted
      if (remaining < 1) break
      
      block_size <- numParticles - numberAccepted
      
      if (tried > 0 && numberAccepted > 0)
        block_size <- block_size * tried / numberAccepted # 1 / (number_accepted / tried)
      
      block_size <- floor(block_size)
      block_size <- min(block_size, 10000)
      
      cat(numberAccepted, block_size, "\n")
      
      param_matrix <- list()
      
      for (i in 1:block_size) {
        params <- c(species_fallout, 0, 0, 1)
        # get a parameter combination
        if (t == 1)  {
          params <- getRandomVals(species_fallout)
          if (fit_order == -1) {
            params[4] <- sample(1:6, 1)
          } else {
            params[4] <- fit_order
          }
        } else {
          params <- getFromPrevious(weights,
                                    disp_vals, filt_vals, comp_vals, order_vals)
          params <- perturb(params, sigma, fit_order)
  
          # we need to know which parameter was perturbed,
          # to be able to calculate its weight later
          changed <- params[5]
          params <- params[1:4]
        }
        
        param_matrix[[i]] <- params
      }

      # now we simulate them all
      process_particle <- function(local_params) {
        return(get_fit(local_params, species, abundances, taxa,
                       esppres, community_number, n_traits,
                       species_fallout, fit_order, Ord, res, optimum,
                       summary_stats, sd_vals))
      }
      
      local_results <- parallel::mclapply(param_matrix, process_particle,
                                          mc.cores = num_threads)
      
      #local_results <- list()
      #for (i in 1:length(param_matrix)) {
      #  local_results[[i]] <- process_particle(param_matrix[[i]])
      #}
      
      
      
      for (i in 1:length(local_results)) {
        local_res <- local_results[[i]]
        if (is.na(local_res$fit) || is.nan(local_res$fit)) next
        
        if (local_res$fit < threshold) {
          numberAccepted <- numberAccepted + 1
          
          next_disp [numberAccepted]  <- local_res$params[1]
          next_filt [numberAccepted]  <- local_res$params[2]
          next_comp [numberAccepted]  <- local_res$params[3]
          next_order[numberAccepted] <- local_res$params[4]
          
          fits[numberAccepted] <- local_res$fit
          rich_vec[numberAccepted] <- local_res$FRic[[1]]
          eve_vec[numberAccepted] <- local_res$FEve[[1]]
          div_vec[numberAccepted] <- local_res$FDiv[[1]]
          opt_vec[numberAccepted] <- local_res$mean_optimum
          
          if (t == 1) {
            next_weights[numberAccepted] <- 1
          } else {
            next_weights[numberAccepted] <-
              calculateWeight(params = local_res$params, 
                              target = changed, 
                              sigma = sigma, 
                              disp_vals = disp_vals, 
                              filt_vals = filt_vals, 
                              comp_vals = comp_vals,
                              order_vals = order_vals,
                              weights = weights)
          }
          
          if ((numberAccepted) %% (numParticles / PRINT_FREQ) == 0) {
            cat("**")
            flush.console()
          }
        }
      }
      
       # function to accept / reject models based on the fit
      tried <- tried + block_size
      if (tried > (1 / stopRate) && tried > 50)  {
        # do not check every particle if the acceptance rate is OK
        if (numberAccepted / tried < stopRate) {
          stop_iteration <- 1
          break
        }
      }

      if (t >= stop_at_iteration) {
        stop_iteration <- 1
        break
      }
    }

    # replace values
    disp_vals  <- next_disp[1:numParticles]
    filt_vals  <- next_filt[1:numParticles]
    comp_vals  <- next_comp[1:numParticles]
    order_vals <- next_order[1:numParticles]
    weights    <- next_weights[1:numParticles]
    
    rich_vec <- rich_vec[1:numParticles]
    eve_vec  <- eve_vec[1:numParticles]
    div_vec  <- div_vec[1:numParticles]
    opt_vec  <- opt_vec[1:numParticles]
    fits     <- fits[1:numParticles]
    

    output <- cbind(disp_vals, filt_vals, comp_vals, rich_vec,
                    eve_vec, div_vec, opt_vec, fits, weights, order_vals)
    file_name <- paste("particles_t=", t, ".txt", sep = "", collapse = NULL)
    write.table(output, file_name, row.names = FALSE, col.names = FALSE)

    if (fit_order) {
      cat(" ", mean(disp_vals), mean(filt_vals),
          mean(comp_vals), mean(order_vals), 
          "\t", "accept rate = ", numberAccepted / (tried - 1), "\n")
    } else {
    cat(" ", mean(disp_vals), mean(filt_vals), mean(comp_vals), 
        "\t", "accept rate = ", numberAccepted / (tried - 1), "\n")
    }
    # and reset
    next_weights <- rep(1,numParticles)
    next_disp <- 1:numParticles
    next_filt <- 1:numParticles
    next_comp <- 1:numParticles
    next_order <- 1:numParticles

    if (stop_iteration == 1) {
      break
    }
    t <- t + 1
  }


  if (t >= 2) {
    d <- read.table(paste("particles_t=", t - 1, ".txt", sep = "",
                          collapse = NULL), header = FALSE)
  } else {
      stop("ABC_SMC: ",
           "Can't stop at iteration 1 - 
           please set stop_at_iteration to 2 if you only 
           want to generate from the prior")
  }
  output <- list( DA = d[, 1], HF = d[, 2], LS = d[, 3])
  if (fit_order) {
    output <- list( DA = d[, 1], HF = d[, 2], LS = d[, 3],
                    OR = d[, 10])
  }
  return(output)
}
