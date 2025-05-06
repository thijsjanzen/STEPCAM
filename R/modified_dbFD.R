# the ordinationAxes function was adapted from the supplementary
# material from Hauffe et al. 2016, which be accessed here:
# doi:10.5194/bg-13-2901-2016-supplement.
# The original paper:
# Hauffe, Torsten, Christian Albrecht, and Thomas Wilke. "Assembly processes
# of gastropod community change with horizontal and vertical zonation
# in ancient Lake Ohrid: a metacommunity speciation perspective."
# Biogeosciences 13.10 (2016): 2901-2911.
# Which can be accessed here:
# http://www.biogeosciences.net/13/2901/2016/bg-13-2901-2016.pdf
ordinationAxes <- function(x, corr = c("sqrt", "cailliez", "lingoes", "none"),
                           ord = c("podani", "metric"), w,
                           asym.bin = NULL, messages = FALSE, stand.x = TRUE) {
  dist.bin <- 2
  tol <- .Machine$double.eps
  corr <- match.arg(corr)
  ord <- match.arg(ord)
  if (is.matrix(x) | is.data.frame(x)) {
    is.dist.x <- FALSE
    s.x <- nrow(x)
    t.x <- ncol(x)
    if (is.null(row.names(x)))
      stop("'x' must have row names.", "\n")
    else x.rn <- row.names(x)
  }
  if (is.vector(x) | is.factor(x)) {
    is.dist.x <- FALSE
    s.x <- length(x)
    t.x <- 1
    if (is.null(names(x)))
      stop("'x' must have names.", "\n")
    else x.rn <- names(x)
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    is.dist.x <- TRUE
    s.x <- attr(x, "Size")
    t.x <- 1
    if (is.null(attr(x, "Labels")))
      stop("'x' must have labels.", "\n")
    else x.rn <- attr(x, "Labels")
  }
  if (missing(w))
    w <- rep(1, t.x)/sum(rep(1, t.x))
  if (is.matrix(x) | is.data.frame(x)) {
    x <- data.frame(x)
    if (t.x >= 2) {
      x.class <- sapply(x, data.class)
      if (any(x.class == "character")) {
        index <- which(x.class == "character")
        for (i in index) {
          x[, i] <- as.factor(x[, i])
        }
      } else {
        x <- x
      }
      if (all(x.class == "numeric") & all(!is.na(x))) {
        if (length(unique(w)) == 1) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        else {
          x.dist <- FD::gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
        }
      }
      else {
        x.dist <- FD::gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
      }
    }
    if (t.x == 1) {
      if (is.numeric(x[, 1])) {
        if (all(!is.na(x))) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          if (messages)
            cat("Warning:",
                "Species with missing trait values have been excluded.",
                "\n")
        }
      }
      if (is.factor(x[, 1]) | is.character(x[, 1])) {
        if (is.ordered(x[, 1]))
          x <- x
        else x[, 1] <- as.factor(x[, 1])
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          x.rn <- x.rn[-pos.NA]
          if (messages)
            cat("Warning:",
                "Species with missing trait values have been excluded.",
                "\n")
        }
        if (is.ordered(x[, 1])) {
          x.s <- data.frame(rank(x[, 1]))
          names(x.s) <- x.rn
          x.dist <- dist(x.s)
        }
        else {
          x.f <- as.factor(x[, 1])
          x.dummy <- diag(nlevels(x.f))[x.f, ]
          x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
          sequence <- 1:10
          if (all(dist.bin != sequence[any(sequence)]))
            stop("'dist.bin' must be an integer between 1 and 10.",
                 "\n")
          x.dist <- ade4::dist.binary(x.dummy.df, method = dist.bin)
        }
      }
    }
  }
  if (is.vector(x) & is.numeric(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      x.rn <- x.rn[-pos.NA]
      if (messages)
        cat("Warning: Species with missing trait values have been excluded.",
            "\n")
    }
    else x <- x
    x.s <- scale(x, center = TRUE, scale = stand.x)
    x.dist <- dist(x.s)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (is.vector(x) & is.character(x)) {
    x <- as.factor(x)
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      x.rn <- x.rn[-pos.NA]
      if (messages)
        cat("Warning: Species with missing trait values have been excluded.",
            "\n")
    }
    else x <- x
    #dimnames(x) <- list(x.rn, "Trait")
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)]))
      stop("'dist.bin' must be an integer between 1 and 10.",
           "\n")
    x <- data.frame(x)
    x.dist <- ade4::dist.binary(x.dummy.df, method = dist.bin)
  }
  if (is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      x.rn <- x.rn[-pos.NA]
      cat("Warning: Species with missing trait values have been excluded.",
          "\n")
    }
    else x <- x
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
    x.dist <- FD::gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
  }
  if (is.factor(x) & !is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      x.rn <- x.rn[-pos.NA]
      if (messages)
        cat("Warning: Species with missing trait values have been excluded.",
            "\n")
    }
    else x <- x
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)]))
      stop("'dist.bin' must be an integer between 1 and 10.",
           "\n")
    x.dist <- ade4::dist.binary(x.dummy.df, method = dist.bin)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    if (any(is.na(x)))
      stop("When 'x' is a distance matrix,",
           "it cannot have missing values (NA).",
           "\n")
    x.dist <- x
  }
  if (any(is.na(x.dist)))
    stop("NA's in the distance matrix.", "\n")
  if (!is.dist.x) {
    no.traits <- apply(x, 1, function(v) length(v[!is.na(v)]))
    if (any(no.traits == 0))
      stop("At least one species has no trait data.", "\n")
  }
  attr(x.dist, "Labels") <- x.rn
  if (ade4::is.euclid(x.dist))
    x.dist2 <- x.dist
  if (!ade4::is.euclid(x.dist)) {
    if (corr == "lingoes") {
      x.dist2 <- ade4::lingoes(x.dist)
      if (messages)
        cat("Species x species distance matrix was not Euclidean.",
            "Lingoes correction was applied.",
            "\n")
    }
    if (corr == "cailliez") {
      x.dist2 <- ade4::cailliez(x.dist)
      if (messages)
        cat("Species x species distance matrix was not Euclidean.",
            "Cailliez correction was applied.",
            "\n")
    }
    if (corr == "sqrt") {
      x.dist2 <- sqrt(x.dist)
      if (!ade4::is.euclid(x.dist2))
        stop("Species x species distance matrix was still not Euclidean
             after 'sqrt' correction. Use another correction method.",
             "\n")
      if (ade4::is.euclid(x.dist2))
        if (messages)
          cat("Species x species distance matrix was not Euclidean.",
              "'sqrt' correction was applied.",
              "\n")
    }
    if (corr == "none") {
      x.dist2 <- ade4::quasieuclid(x.dist)
      if (messages)
        cat("Species x species distance was not Euclidean,",
            "but no correction was applied. Only the PCoA axes",
            "with positive eigenvalues were kept.",
            "\n")
    }
  }
  x.pco <- ade4::dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
  # Best is to return the whole object, because dbFD needs sometimes
  # more than the axes alone?!
  return(x.pco)
}

# Function to dermine the number of PCoA axes used and number of species.
# (Both constant in STEPCAM and therefor only camculated once.)
detMnbsp <- function(x.pco, a){
  tol <- .Machine$double.eps
  a <- as.matrix(a)
  c <- nrow(a) # Number of communities
  traits <- x.pco$li
  # Number of species per community where traits are present
  # do I need it later? Yes!
  nb.sp <- numeric(c)
  for (i in 1:c) {
    sp.pres <- which(a[i, ] > 0)
    traits.sp.pres <- traits[sp.pres, , drop = FALSE]
    traits.sp.pres[traits.sp.pres != 0 & abs(traits.sp.pres) < tol] <- 0
    nb.sp[i] <- nrow(unique(traits.sp.pres))
  }
  min.nb.sp <- min(nb.sp)
  # Minimum number of species in one of the communities - 1
  m.max <- min.nb.sp - 1
  #Warning <- FALSE
  if (min.nb.sp < 3) {
    nb.sp2 <- nb.sp[nb.sp > 2]
    m.max <- min(nb.sp2) - 1
  } else {
    m.max <- m.max
  }
  m <- m.max
  Res <- list()
  Res[[1]] <- m
  Res[[2]] <- nb.sp
  return(Res)
}

# the strippeddbFD function was adapted from the supplementary
# material from Hauffe et al. 2016, which be accessed here:
# doi:10.5194/bg-13-2901-2016-supplement.
# The original paper:
# Hauffe, Torsten, Christian Albrecht, and Thomas Wilke. "Assembly processes
# of gastropod community change with horizontal and vertical zonation
# in ancient Lake Ohrid: a metacommunity speciation perspective."
# Biogeosciences 13.10 (2016): 2901-2911.
# Which can be accessed here:
# http://www.biogeosciences.net/13/2901/2016/bg-13-2901-2016.pdf



# Strip away radically everything except FRic, FDiv, and FEve
# Because this function needs to be much much faster!
# All communities are run in a for loop which cannot be vectorized
strippedDbFd <- function(x.pco, a, m, nb.sp){
  #tol <- .Machine$double.eps
  a <- as.matrix(a)
  c <- nrow(a) # Number of communities
  traits <- x.pco$li
  Warning <- FALSE
  # If there is a community with less species than ordination axes
  if (m < x.pco$nf){
    traits.FRic <- x.pco$li[, 1:m] # Use only the first m axes
  }
  if (m >= x.pco$nf){
    traits.FRic <- x.pco$li
  }
  FRic <- rep(NA, c)
  names(FRic) <- row.names(a)
  FEve <- FRic
  FDiv <- FRic
  #AbundRel <- a/rowSums(a)


  for (i in 1:c) { # For each community
    sppres <- which(a[i, ] > 0) # Present species
    S <- length(sppres)
    tr <- data.frame(traits[sppres, ]) # Axes coordinates of the present species
    tr.FRic <- data.frame(traits.FRic[sppres, ])
    # Will I need relative abundances of the species?
    ab <- as.matrix(a[i, sppres])
    abundrel <- ab/sum(ab)
    abund2 <- sapply( c(abundrel), function(x) x + abundrel)
    abund2vect <- as.dist(abund2)

    # New part: check range of axes values, because if range is very small,
    # convhulln will fail
    #apply(tr.FRic, 2, range)

    # If there are more than 3 species present
    if (ncol(tr.FRic) > 1 & nb.sp[i] >= 3) {
      if (Warning)
        thresh <- 4
      if (!Warning)
        thresh <- 3
      if (nb.sp[i] >= thresh) {
        # Option QJ is helpfull in case of planar hulls, Pp removes warning
        convhull <- geometry::convhulln(tr.FRic, c("QJ", "FA", "Pp"))
        FRic[i] <- convhull$vol
      }
    }
    if (ncol(tr.FRic) == 1) {
      tr.range <- range(tr.FRic[, 1])
      t.range <- tr.range[2] - tr.range[1]
      FRic[i] <- t.range
    }

    if (nb.sp[i] >= 3) {
      tr.dist <- dist(tr) # pair-wise distance of ordination coordinates
      linkmst <- ape::mst(tr.dist)
      mstvect <- as.dist(linkmst)
      #abund2 <- matrix(0, nrow = S, ncol = S)
      #for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] + abundrel[r]
      # the *apply family is faster than the original code
      # with more than three species
      # Move this outside of the loop 'cause its always the same:
      #abund2 <- sapply( c(abundrel), function(x) x + abundrel)
      #abund2vect <- as.dist(abund2)
      #EW <- rep(0, S - 1)
      #flag <- 1
      #for (m in 1: ((S - 1) * S/2) ) {
      #  if (mstvect[m] != 0) {
      #    EW[flag] <- tr.dist[m]/(abund2vect[m])
      #    flag <- flag + 1
      #  }
      #}
      # Faster:
      EW <- c((tr.dist * mstvect) / abund2vect)
      EW <- EW[EW > 0]

      minPEW <- rep(0, S - 1)
      OdSmO <- 1/(S - 1)
      for (l in 1:(S - 1)) minPEW[l] <- min( (EW[l]/sum(EW)), OdSmO)
      # Slower:
      #sapply( EW/sum(EW), function(x) min( x, OdSmO ))
      FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    }
    if (ncol(tr.FRic) > 1 & nb.sp[i] >= 3) {
      # Option QJ is helpfull in case of planar hulls, Pp removes warning
      vert0 <- geometry::convhulln(tr.FRic, c("Fx TO 'vert.txt'", "QJ", "Pp"))
      vert1 <- scan("vert.txt", quiet = TRUE)
      vert2 <- vert1 + 1
      vertices <- vert2[-1]
      trvertices <- tr.FRic[vertices, ]
      baryv <- colMeans(trvertices) #apply(trvertices, 2, mean)
      #Faster:
      #distbaryv <- rep(0, S)
      #for (j in 1:S)
      # distbaryv[j] <- ( sum( (tr.FRic[j, ] - baryv)^2) )^0.5
      distbaryv <- sqrt( rowSums( (tr.FRic - baryv)^2) )

      meandB <- mean(distbaryv)
      devdB <- distbaryv - meandB
      abdev2 <- abundrel * devdB
      ababsdev2 <- abundrel * abs(devdB)
      FDiv[i] <- (sum(abdev2) + meandB)/(sum(ababsdev2) + meandB)
    }
  }
  res <- list()
  res$FRic <- FRic
  res$FEve <- FEve
  res$FDiv <- FDiv
  return(res)
}
