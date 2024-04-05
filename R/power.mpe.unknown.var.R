power.mpe.unknown.var <- function(K, n = NULL, delta = NULL, Sigma, SD, rho, sig.level = 0.05,
                              power = NULL, M = 10000, n.min = NULL, n.max = NULL, 
                              tol = .Machine$double.eps^0.25, use.uniroot = TRUE){
  ## check of input
  if(missing(K))
    stop("specify the number of co-primary endpoints")
  if(!is.numeric(K))
    stop("'K' must be a natural number > 1")
  K <- as.integer(K)

  if(is.null(n) & is.null(power))
    stop("either 'n' or 'power' must be specified")
  if(!is.null(n) & !is.null(power))

    stop("either 'n' or 'power' must be NULL")
  if(!is.null(n)){
    if(length(n) > 1){
      warning("length of 'n' is greater than 1, only the first entry is used")
      n <- n[1]
    }
    n <- as.integer(n)
  }
  if(!is.null(power)){
    if(power <= 0 | power >= 1)
      stop("power must be in (0, 1)")
  }
  if(is.null(delta))
    stop("expected effect size 'delta' is missing")
  if(length(delta) < 2)
    stop("length of 'delta' < 2: effect for at least two co-primary endpoints is required")
  if(!all(delta > 0))
    stop("all effect sizes need to be positive")
  if(!missing(Sigma)){
    if(nrow(Sigma) != ncol(Sigma))
      stop("covariance matrix 'Sigma' must be quadratic")
    if(nrow(Sigma) != K)
      stop("covariance matrix must have dimension 'K' x 'K'")
    if(max(abs(Sigma - t(Sigma))) > 1e-10)
      stop("matrix 'Sigma' must be symmetric")
  }
  if(missing(Sigma)){
    if(missing(SD) || missing(rho))
      stop("if 'Sigma' is missing 'SD' and 'rho' must be given.")
    if(length(SD) != K)
      stop("length of 'SD' must be equal to 'K'")
    if(length(rho) != 0.5*K*(K-1))
      stop("length of 'rho' must be equal to '0.5*K*(K-1)'")
    Sigma <- matrix(0, nrow = K, ncol = K)
    iter <- 0
    for(i in 1:(K-1)){
      for(j in (i+1):K){
        iter <- iter + 1
        Sigma[i,j] <- rho[iter]*SD[i]*SD[j]
      }
    }
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- SD^2
  }
  if(!all(eigen(Sigma)$values > 0))
    stop("matrix 'Sigma' must be positive definite")
  Sigma.cor <- cov2cor(Sigma)

  if(sig.level <= 0 | sig.level >= 1)
    stop("significance level must be in (0, 1)")

  if(length(M) > 1){
    warning("length of 'M' is greater than 1, only the first entry is used")
    M <- M[1]
  }
  M <- as.integer(M)
  if(is.null(n)){
    if(length(n.min) > 1){
      warning("length of 'n.min' is greater than 1, only the first entry is used")
      n.min <- n.min[1]
    }
    n.min <- as.integer(n.min)
    if(length(n.max) > 1){
      warning("length of 'n.max' is greater than 1, only the first entry is used")
      n.max <- n.max[1]
    }
    n.max <- as.integer(n.max)
    if(n.min < 4)
      stop("'n.min' must be larger or equal 4")
    if(n.min >= n.max)
      stop("'n.min' must be smaller than 'n.max'")
  }
  ## calculations
  if(is.null(power)){
    std.effect <- delta/sqrt(diag(Sigma))
    probs <- numeric(M)
    Ws <- rWishart(M, df = n-2, Sigma = Sigma.cor)
    for(i in 1:M){
      Wi <- diag(Ws[,,i])
      ci <- qt(1-sig.level, df = n-2)*sqrt(Wi/(n-2)) - sqrt(n/2)*std.effect
      probs[i] <- pmvnorm(upper = -ci, sigma = Sigma.cor)
    }
    power <- mean(probs)
  }

  if(is.null(n)){
    std.effect <- delta/sqrt(diag(Sigma))
    ssize.fct <- function(n, std.effect, Sigma.cor, power, M, verbose = FALSE){
      probs <- numeric(M)
      Ws <- rWishart(M, df = n-2, Sigma = Sigma.cor)
      for(i in 1:M){
        Wi <- diag(Ws[,,i])
        ci <- qt(1-sig.level, df = n-2)*sqrt(Wi/(n-2)) - sqrt(n/2)*std.effect
        probs[i] <- pmvnorm(upper = -ci, sigma = Sigma.cor) - power
      }
      if(verbose) cat("Current precision:\t", mean(probs), "\n")
      mean(probs)
    }
    if(use.uniroot){
      n <- uniroot(ssize.fct, c(n.min, n.max), tol = tol, extendInt = "yes",
                   std.effect = std.effect, Sigma.cor = Sigma.cor, power = power, M = M,
                   verbose = TRUE)$root
    }else{
      ns <- n.min:n.max
      res <- sapply(ns, ssize.fct, std.effect = std.effect, Sigma.cor = Sigma.cor,
                    power = power, M = M)
      ns.pos <- ns[res > 0]
      res.pos <- res[res > 0]
      cat("Precision:\t", min(res.pos), "\n")
      n <- ns.pos[which.min(res.pos)]
    }
  }
  NOTE <- "n is number in *each* group"
  METHOD <- "Power calculation for multiple co-primary endpoints (covariance unknown)"
  structure(list(n = n, delta = delta, SD = sqrt(diag(Sigma)),
                 rho = Sigma.cor[lower.tri(Sigma.cor)], Sigma = Sigma, sig.level = sig.level,
                 power = power, note = NOTE, method = METHOD),
            class = "power.mpe.test")
}
