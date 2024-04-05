power.mpe.atleast.one <- function(K, n = NULL, delta = NULL, Sigma, SD, rho, sig.level = 0.05/K,
                                  power = NULL, n.max = 1e5, tol = .Machine$double.eps^0.25){
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
  
  ## calculations
  if(is.null(power)){
    std.effect <- delta/sqrt(diag(Sigma))
    z.alpha <- qnorm(1-sig.level)
    crit.vals <- z.alpha - sqrt(n/2)*std.effect
   power <- 1-pmvnorm(lower = -crit.vals, sigma = Sigma.cor) 
  }
  if(is.null(n)){
    std.effect <- delta/sqrt(diag(Sigma))
    z.alpha <- qnorm(1-sig.level)
    ssize.fct <- function(n, std.effect, z.alpha, Sigma.cor, power){
      crit.vals <- z.alpha - sqrt(n/2)*std.effect
      
      pmvnorm(lower = -crit.vals, sigma = Sigma.cor) - (1- power)
    }
    n <- uniroot(ssize.fct, c(2, n.max), tol = tol, extendInt = "yes",
                 std.effect = std.effect, z.alpha = z.alpha, Sigma.cor = Sigma.cor,
                 power = power)$root
  }
  NOTE <- "n is number in *each* group"
  METHOD <- "Power calculation for multiple primary endpoints for at least one endpoint"
  structure(list(n = n, delta = delta, SD = sqrt(diag(Sigma)),
                 rho = Sigma.cor[lower.tri(Sigma.cor)], Sigma = Sigma,
                 sig.level = sig.level,
                 power = power, note = NOTE, method = METHOD),
            class = "power.mpe.test")
}
