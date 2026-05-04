# removal of FD from CRAN necessitates this blatant copy:
FD_gowdis_deprec <- function (x, w, asym.bin = NULL, ord = c("podani", "metric", 
                                                      "classic")) 
{
  if (length(dx <- dim(x)) != 2 || !(is.data.frame(x) || is.numeric(x))) 
    stop("x is not a dataframe or a numeric matrix\n")
  n <- dx[1]
  p <- dx[2]
  ord <- match.arg(ord)
  varnames <- dimnames(x)[[2]]
  if (!missing(w)) {
    if (length(w) != p | !is.numeric(w)) 
      stop("w needs to be a numeric vector of length = number of variables in x\n")
    if (all(w == 0)) 
      stop("Cannot have only 0's in 'w'\n")
    w <- w/sum(w)
  }
  else w <- rep(1, p)/sum(rep(1, p))
  if (is.data.frame(x)) {
    type <- sapply(x, data.class)
  }
  else {
    type <- rep("numeric", p)
    names(type) <- colnames(x)
  }
  if (any(type == "character")) 
    for (i in 1:p) if (type[i] == "character") 
      x[, i] <- as.factor(x[, i])
  is.bin <- function(k) all(k[!is.na(k)] %in% c(0, 1))
  bin.var <- rep(NA, p)
  names(bin.var) <- varnames
  for (i in 1:p) bin.var[i] <- is.bin(x[, i])
  if (any(type[bin.var] != "numeric")) 
    stop("Binary variables should be of class 'numeric'\n")
  type[type %in% c("numeric", "integer")] <- 1
  type[type == "ordered"] <- 2
  type[type %in% c("factor", "character")] <- 3
  type[bin.var] <- 4
  if (!is.null(asym.bin)) {
    if (!all(bin.var[asym.bin])) 
      stop("Asymetric binary variables must only contain 0 or 1\n")
    else type[asym.bin] <- 5
  }
  type <- as.numeric(type)
  x <- data.matrix(x)
  if (any(type == 2)) {
    if (ord != "classic") 
      for (i in 1:p) if (type[i] == 2) 
        x[, i] <- rank(x[, i], na.last = "keep")
    else for (i in 1:p) if (type[i] == 2) 
      x[, i] <- as.numeric(x[, i])
  }
  range.Data <- function(v) {
    r.Data <- range(v, na.rm = T)
    res <- r.Data[2] - r.Data[1]
    return(res)
  }
  range2 <- apply(x, 2, range.Data)
  comp.Timax <- function(v) {
    Ti.max <- max(v, na.rm = T)
    no.na <- v[!is.na(v)]
    res <- length(no.na[no.na == Ti.max])
    return(res)
  }
  Timax <- apply(x, 2, comp.Timax)
  comp.Timin <- function(v) {
    Ti.min <- min(v, na.rm = T)
    no.na <- v[!is.na(v)]
    res <- length(no.na[no.na == Ti.min])
    return(res)
  }
  Timin <- apply(x, 2, comp.Timin)
  if (ord == "podani") 
    pod <- 1
  else pod <- 2
  res <- .C("gowdis", as.double(x), as.double(w), as.integer(type), 
            as.integer(n), as.integer(p), as.double(range2), as.integer(pod), 
            as.double(Timax), as.double(Timin), res = double(n * (n - 1)/2), NAOK = T, PACKAGE = "FD")$res
  type[type == 1] <- "C"
  type[type == 2] <- "O"
  type[type == 3] <- "N"
  type[type == 4] <- "B"
  type[type == 5] <- "A"
  if (any(is.na(res))) 
    attr(res, "NA.message") <- "NA's in the dissimilarity matrix!"
  attr(res, "Labels") <- dimnames(x)[[1]]
  attr(res, "Size") <- n
  attr(res, "Metric") <- "Gower"
  attr(res, "Types") <- type
  class(res) <- "dist"
  return(res)
}

FD_functcomp <- function(x, a, CWM.type = c("dom", "all"), bin.num = NULL) 
{
  if (!is.matrix(x) & !is.data.frame(x)) 
    stop("'x' must be a matrix or a data frame.", "\n")
  else x <- data.frame(x)
  if (!is.matrix(a)) 
    stop("'a' must be a matrix.", "\n")
  if (is.null(row.names(x))) 
    stop("'x' must have row names.", "\n")
  else x.n <- row.names(x)
  if (is.null(colnames(a))) 
    stop("'a' must have column names.", "\n")
  else a.n <- colnames(a)
  s.x <- dim(x)[1]
  s.a <- dim(a)[2]
  if (s.x != s.a) 
    stop("Different number of species in 'x' and 'a'.", "\n")
  if (any(x.n != a.n)) 
    stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  com <- dim(a)[1]
  t <- dim(x)[2]
  com.names <- row.names(a)
  sp.names <- row.names(x)
  tr.names <- names(x)
  a[which(is.na(a))] <- 0
  CWM.type <- match.arg(CWM.type)
  is.bin <- function(k) all(k[!is.na(k)] %in% c(0, 1))
  bin.var <- rep(NA, t)
  names(bin.var) <- tr.names
  for (i in 1:t) bin.var[i] <- is.bin(x[, i])
  if (!all(bin.var[bin.num])) 
    stop("'bin.num' points to non-binary variables.\n")
  bin.var[bin.num] <- FALSE
  type <- sapply(x, data.class)
  type[type %in% c("numeric", "integer")] <- "C"
  type[type == "ordered"] <- "O"
  type[type == "factor"] <- "N"
  type[bin.var] <- "B"
  sum.a <- apply(a, 1, sum)
  a <- a / sum.a
  a <- t(a)
  a <- data.frame(a)
  temp <- list()
  for (i in 1:t) {
    if (type[i] == "C") {
      vec <- numeric(com)
      for (j in 1:com) vec[j] <- weighted.mean(x[, i], 
                                               a[, j], na.rm = T)
      temp[[i]] <- matrix(vec, com, 1, dimnames = list(com.names, 
                                                       tr.names[i]))
    }
    else {
      x[, i] <- as.factor(x[, i])
      fac <- data.frame()
      which.dom <- rep(NA, com)
      for (k in 1:com) {
        temp2 <- tapply(a[, k], x[, i], sum)
        fac <- rbind(fac, temp2)
        which.dom[k] <- sample(levels(x[, i])[which(fac[k, 
        ] == max(fac[k, ]))], size = 1)
      }
      colnames(fac) <- paste(tr.names[i], "_", levels(x[, 
                                                        i]), sep = "")
      rownames(fac) <- com.names
      which.dom <- data.frame(which.dom)
      colnames(which.dom) <- tr.names[i]
      rownames(which.dom) <- com.names
      if (CWM.type == "dom") 
        temp[[i]] <- which.dom
      if (CWM.type == "all") 
        temp[[i]] <- fac
    }
  }
  temp <- data.frame(temp)
  return(temp)
}

FD_fdisp <- function(d, a, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (!is.matrix(a)) 
    stop("'a' must be a matrix.")
  if (ncol(a) != n) 
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a)) 
    stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    stop("At least one community has zero-sum abundances (no species).", 
         "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    stop("At least one species does not occur in any community (zero total abundance across all communities).", 
         "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- ade4::bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}

