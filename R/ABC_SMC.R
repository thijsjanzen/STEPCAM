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
#
# This has to evaluate the density of the kernel perturb() ACTUALLY uses.
# perturb() draws an ordered pair (a, b) from the three processes, moves
# component a by a rounded normal step with sd sigma[a], and absorbs the
# same amount out of component b; the third component is untouched. So a
# move is a 1-D step along one of the three edge directions (1,-1,0),
# (1,0,-1), (0,1,-1) - NOT an isotropic step in the 2 free dimensions of the
# simplex. The kernel is therefore degenerate, supported on a union of three
# lines through the current particle, and has no 2-D density: scoring it with
# a product of two normals is as wrong as scoring a single component.
#
# For a candidate reached from a previous particle by delta, the probability
# is a mixture over the six equally likely ordered pairs. Each unordered pair
# {a, b} contributes (1/6) * (P(step = delta[a] | sigma[a])
#                             + P(step = -delta[a] | sigma[b])),
# and the step is rounded to an integer, so the probability is the normal
# mass over the rounding interval rather than a density. A delta of exactly
# zero is reachable through all three pairs and the triple contribution is
# the correct total. Clamping at 0 and max_number is ignored; it moves a
# little extra mass onto the boundary.
calculateWeight <- function(params, target, sigma,
                            disp_vals, filt_vals, comp_vals, order_vals,
                            weights)  {
  sigma <- expandSigma(sigma)
  prev <- cbind(disp_vals, filt_vals, comp_vals)

  delta <- sweep(-prev, 2, params[1:3], "+")   # params - prev, row-wise

  step_mass <- function(eps, sd_val) {
    stats::pnorm(eps + 0.5, 0, sd_val) - stats::pnorm(eps - 0.5, 0, sd_val)
  }

  diff_prob <- numeric(nrow(prev))
  for (pr in list(c(1, 2), c(1, 3), c(2, 3))) {
    a <- pr[1]
    b <- pr[2]
    untouched <- setdiff(1:3, pr)
    eps <- delta[, a]
    # the move is only reachable along this edge if the third component did
    # not change and the other two are equal and opposite
    reachable <- (delta[, untouched] == 0) & (delta[, b] == -eps)
    diff_prob <- diff_prob +
      reachable * (1 / 6) * (step_mass(eps, sigma[a]) +
                             step_mass(-eps, sigma[b]))
  }

  # we have to multiply with the ordering as well
  diff_order <- 1 - (params[4] != order_vals)
  # 90 % prob of remaining the same
  diff_order <- diff_order * 0.9
  diff_order[diff_order == 0] <- 0.1

  vals <- weights * diff_prob * diff_order

  denom <- sum(vals)
  # a particle no previous particle could plausibly have produced would
  # otherwise get an infinite weight and wipe out the whole population
  if (!is.finite(denom) || denom <= 0) return(0)

  return( 1 / denom ) # prior density is 1.
}

# normalize all the weights of all particles such that they sum to 1
normalizeWeights <- function(x) {
  sum_x <- sum(x)
  x <- x / sum_x
  return(x)
}

# sigma may be a single width (the historical form, and what the tests pass)
# or one width per process. Expand to length 3 so that callers can always
# index it by component. Both rnorm() and dnorm() silently recycle a vector
# rather than erroring, so indexing has to be explicit everywhere sigma is
# used; this helper is the single place that decides what sigma means.
expandSigma <- function(sigma) {
  if (length(sigma) == 1) return(rep(sigma, 3))
  if (length(sigma) != 3) {
    stop("expandSigma: ",
         "sigma must have length 1 or 3, got ", length(sigma))
  }
  return(sigma)
}

# function to randomly change the contribution of one of the processes:
perturb <- function(p, sigma, fit_order)  {
  sigma <- expandSigma(sigma)
  params <- p
  max_number <- sum(p[1:3])
  numbers <- 1:3

  x <- sample(numbers, 3, replace = FALSE)

  oldval <- params[x[1]]

  # sd must be the width for THIS component: passing the whole vector would
  # silently use sigma[1] for every process
  params[x[1]] <- round(params[x[1]] +
                          rnorm(1, mean = 0, sd = sigma[x[1]]), 0)
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

# Run the model once and return the four summary statistics plus the
# resulting distance to the observed community.
get_fit_single <- function(params, species, abundances, taxa,
                           esppres, community_number, n_traits,
                           species_fallout, Ord, res, optimum,
                           summary_stats, sd_vals) {
  allcommunities <- STEPCAM(params, species, abundances, taxa,
                            esppres, community_number, n_traits,
                            species_fallout)

  traits <- as.data.frame(species[, c(2:(n_traits + 1))],
                          row.names = c(1:taxa))

  communities <- as.data.frame(t(allcommunities))
  present_species <- which(communities > 0)

  # nb.sp must describe THIS simulated community, not plot 1 of the
  # observed data set. The simulated community always holds
  # taxa - species_fallout species.
  nb_sp_sim <- length(present_species)

  FD_output <- strippedDbFd(Ord, communities,
                            m = res[[1]], nb.sp = nb_sp_sim)

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
              FRic = FRic[[1]],
              FEve = FEve[[1]],
              FDiv = FDiv[[1]],
              mean_optimum = mean_optimum[1]))
}

# Dispersal is stochastic but filtering and limiting similarity are
# deterministic, so the variance of `fit` grows with the number of dispersal
# events. Under a hard acceptance threshold that hands high-dispersal
# particles a free advantage: they get repeated draws at the lower tail,
# while a particle with zero dispersal events produces one single fit value
# and is either always or never accepted. Averaging over n_replicates
# simulations makes acceptance depend on the expected fit instead, which
# removes that bias. n_replicates = 1 reproduces the old behaviour.
get_fit <- function(params, species, abundances, taxa,
                    esppres, community_number, n_traits,
                    species_fallout, Ord, res, optimum,
                    summary_stats, sd_vals, n_replicates = 1) {
  taxa <- length(abundances[1, ])

  reps <- vector("list", n_replicates)
  for (r in seq_len(n_replicates)) {
    reps[[r]] <- get_fit_single(params, species, abundances, taxa,
                                esppres, community_number, n_traits,
                                species_fallout, Ord, res, optimum,
                                summary_stats, sd_vals)
  }

  pull <- function(nm) {
    v <- vapply(reps, function(z) as.numeric(z[[nm]]), numeric(1))
    if (all(is.na(v))) return(NA_real_)
    mean(v, na.rm = TRUE)
  }

  return(list(fit = pull("fit"),
              FRic = pull("FRic"),
              FEve = pull("FEve"),
              FDiv = pull("FDiv"),
              mean_optimum = pull("mean_optimum"),
              params = params))
}


ABC_SMC <- function(numParticles, species_fallout, taxa, esppres, n_traits,
                    sd_vals, summary_stats, community_number, species,
                    abundances, frequencies, stopRate, Ord, 
                    continue_from_file = TRUE, stop_at_iteration = 50,
                    fit_order = FALSE,
                    num_threads = 1,
                    n_replicates = 1,
                    threshold_quantile = NULL,
                    adapt_sigma = !is.null(threshold_quantile))  {

  # threshold_quantile is a PROBABILITY, not a switch. Passing TRUE reads
  # naturally as "turn this on", but TRUE coerces to 1, and quantile(fits,
  # probs = 1) is the maximum of the previous distances - a threshold almost
  # nothing fails. The run then accepts nearly every particle, applies no
  # selection at all, and (with adapt_sigma on) the kernel width runs away.
  # That failure is silent, so reject it here rather than let it through.
  if (!is.null(threshold_quantile)) {
    if (is.logical(threshold_quantile)) {
      stop("ABC_SMC: ",
           "threshold_quantile must be a probability, not TRUE/FALSE. ",
           "TRUE would be read as quantile 1 (the worst previous distance), ",
           "which accepts almost everything. Use a value such as 0.5, ",
           "or NULL to keep the fixed threshold schedule.")
    }
    if (!is.numeric(threshold_quantile) ||
        length(threshold_quantile) != 1 ||
        !is.finite(threshold_quantile) ||
        threshold_quantile <= 0 || threshold_quantile >= 1) {
      stop("ABC_SMC: ",
           "threshold_quantile must be a single number strictly between ",
           "0 and 1 (0.5 is a reasonable default), or NULL. Got: ",
           paste(utils::capture.output(str(threshold_quantile)),
                 collapse = " "))
    }
  }

  if (!is.logical(adapt_sigma) || length(adapt_sigma) != 1 ||
      is.na(adapt_sigma)) {
    stop("ABC_SMC: ",
         "adapt_sigma must be TRUE or FALSE")
  }

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

  # one perturbation width per process (dispersal, filtering, competition)
  sigma <- rep(1, 3)
  t <- 1

  f <- list.files(pattern = "particles_t=")
  if (length(f) > 0 && continue_from_file == TRUE) {
    cat("Found previous output, continuing from that output\n")
    flush.console()

    d <- NULL
    for (t1 in stop_at_iteration:1) {
      file_name <- paste0("particles_t=", t1, ".txt")
      if (!file.exists(file_name)) next

      candidate <- try(utils::read.table(file_name, header = FALSE),
                       silent = TRUE)
      # an unreadable or truncated file is not a usable restart point
      if (inherits(candidate, "try-error")) next
      if (nrow(candidate) < numParticles || ncol(candidate) < 10) next
      # a partially written iteration contains NA rows - skip it
      if (anyNA(candidate[seq_len(numParticles), 1:10])) next
      # isTRUE() guards against NA propagating into the if()
      if (isTRUE(candidate[numParticles, 1] == numParticles)) next

      d <- candidate
      t <- t1 + 1
      break
    }

    # no usable file: fall through and start from the prior rather than
    # calling break outside a loop (which is a runtime error in R)
    if (is.null(d)) {
      cat("No usable previous output found, starting from the prior\n")
      flush.console()
      t <- 1
    } else {
      disp_vals  <- d[, 1]
      filt_vals  <- d[, 2]
      comp_vals  <- d[, 3]
      fits       <- d[, 8]
      weights    <- d[, 9]
      order_vals <- d[, 10]
    }
  }

  # continuously sampling
  while (t < stop_at_iteration)  {
    cat("\nGenerating Particles for iteration\t", t, "\n")
    cat("0--------25--------50--------75--------100\n")
    cat("*")
    flush.console()
    PRINT_FREQ <- 20

    numberAccepted <- 0
    if (t != 1) weights <- normalizeWeights(weights)

    # The hard-coded schedule 200 * exp(-0.5 * t) takes no account of what
    # distance is actually achievable for this data set: if the best possible
    # fit is large the chain dies after one or two iterations and returns the
    # prior, and if it is tiny the chain over-converges. With
    # threshold_quantile set, each threshold is instead a quantile of the
    # distances accepted in the previous iteration, which is the usual
    # adaptive ABC-SMC choice. Iteration 1 still needs a fixed starting value.
    if (is.null(threshold_quantile) || t == 1 || length(fits) == 0) {
      threshold <- 200 * exp(-0.5 * t)
    } else {
      prev_fits <- fits[is.finite(fits)]
      if (length(prev_fits) == 0) {
        threshold <- 200 * exp(-0.5 * t)
      } else {
        threshold <- stats::quantile(prev_fits, probs = threshold_quantile,
                                     names = FALSE)
        # min(threshold, max(prev_fits)) was a no-op, since a quantile can
        # never exceed the maximum. Compare against the PREVIOUS threshold
        # instead, so the sequence can genuinely never grow or stall.
        if (exists("prev_threshold", inherits = FALSE)) {
          threshold <- min(threshold, prev_threshold)
        }
      }
    }
    prev_threshold <- threshold
    cat("threshold =", signif(threshold, 4), "\n")
    flush.console()

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
        changed <- NA
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
        
        param_matrix[[i]] <- list(params = params, changed = changed)
      }

      # now we simulate them all
      # NB: the result is named fit_res, not res - naming it `res` shadows the
      # `res` holding detMnbsp() output that is passed straight back into
      # get_fit() on the same line.
      process_particle <- function(local_params) {
        fit_res <- get_fit(local_params$params,
                           species, abundances, taxa,
                           esppres, community_number, n_traits,
                           species_fallout, Ord, res, optimum,
                           summary_stats, sd_vals,
                           n_replicates = n_replicates)
        return(list(res = fit_res, changed = local_params$changed))
      }
      
      local_results <- list()
      if (.Platform$OS.type == "unix" && num_threads > 1) {
        local_results <- parallel::mclapply(param_matrix, process_particle,
                                            mc.cores = num_threads)
      } else {
        for (i in 1:length(param_matrix)) {
          local_results[[i]] <- process_particle(param_matrix[[i]])
        }
      }
      
      
      
      
      
      for (i in 1:length(local_results)) {
        local_res <- local_results[[i]]$res
        if (is.na(local_res$fit) || is.nan(local_res$fit)) next
        
        if (local_res$fit < threshold) {
          numberAccepted <- numberAccepted + 1
          
          next_disp[numberAccepted]  <- local_res$params[1]
          next_filt[numberAccepted]  <- local_res$params[2]
          next_comp[numberAccepted]  <- local_res$params[3]
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
                              target = local_results[[i]]$changed, 
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

    # An iteration that stopped early has fewer than numParticles accepted
    # particles. Writing it would produce a file padded with NA, which both
    # corrupts the posterior and breaks continue_from_file on the next run.
    if (numberAccepted < numParticles) {
      cat("\nIteration", t, "stopped early with only", numberAccepted,
          "of", numParticles, "particles accepted - not written to file.\n")
      flush.console()
      stop_iteration <- 1
      break
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

    # Adapt the perturbation kernel to the spread of the current population.
    # A kernel pinned at sigma = 1 moves one or two species per step, so with
    # a large species pool the chain barely explores the simplex within the
    # number of iterations actually run. The usual ABC-SMC choice is twice the
    # weighted variance of the previous population (Beaumont et al. 2009).
    #
    # This is deliberately tied to threshold_quantile by default. The kernel
    # width is derived from the accepted population and the threshold decides
    # what gets accepted, so the two are coupled: adapting sigma underneath
    # the fixed 200 * exp(-0.5 * t) schedule lets the threshold outrun the
    # kernel, acceptance falls away iteration on iteration, and sigma then
    # collapses with the degenerate population that survives. Adapting both
    # together keeps the acceptance rate roughly stable instead.
    # One width per process rather than one pooled width: the posterior is
    # routinely tighter on one process than another, and a single pooled
    # sigma then over-perturbs the tight one and under-perturbs the wide one.
    if (adapt_sigma) {
      w_norm <- normalizeWeights(weights)
      wvar <- function(x) {
        mu <- sum(w_norm * x)
        sum(w_norm * (x - mu) ^ 2)
      }
      sigma <- sqrt(2 * c(wvar(disp_vals), wvar(filt_vals), wvar(comp_vals)))
      # keep each kernel wide enough to still move at least one species
      sigma[!is.finite(sigma) | sigma < 1] <- 1
      cat("sigma (D,F,C) =", paste(signif(sigma, 4), collapse = ", "), "\n")
      flush.console()
    }

    output <- cbind(disp_vals, filt_vals, comp_vals, rich_vec,
                    eve_vec, div_vec, opt_vec, fits, weights, order_vals)
    file_name <- paste("particles_t=", t, ".txt", sep = "", collapse = NULL)
    utils::write.table(output, file_name,
                       row.names = FALSE, col.names = FALSE)

    accept_rate <- numberAccepted / max(tried, 1)
    if (fit_order) {
      cat(" ", mean(disp_vals), mean(filt_vals),
          mean(comp_vals), mean(order_vals),
          "\t", "accept rate = ", accept_rate, "\n")
    } else {
      cat(" ", mean(disp_vals), mean(filt_vals), mean(comp_vals),
          "\t", "accept rate = ", accept_rate, "\n")
    }
    # and reset
    next_weights <- rep(NA, numParticles)
    next_disp    <- rep(NA, numParticles)
    next_filt    <- rep(NA, numParticles)
    next_comp    <- rep(NA, numParticles)
    next_order   <- rep(NA, numParticles)

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
