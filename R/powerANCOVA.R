power.ancova <- function(n = NULL, mu = NULL, var = 1, nr.covs = 1L, 
                         group.ratio = NULL, contr.mat = NULL, 
                         sig.level = 0.05, power = NULL, n.max = 1000L,
                         rel.tol = .Machine$double.eps^0.25){
  if(is.null(mu)){
    stop("You have to specify the adjusted group means 'mu'.")
  }
  if(is.null(n) && is.null(power)){
    stop("Either 'n' or 'power' must be non NULL.")
  }
  if(!is.numeric(var) || length(var) != 1L || any(var <= 0)){
    stop("The standard error of the residuals 'var' must be a positive (real) number.")
  }
  if(!is.numeric(nr.covs) || length(nr.covs) != 1L || any(nr.covs <= 0)){
    stop("'nr.covs' must be a positive (integer) number.")
  }
  if(!is.integer(nr.covs)) nr.covs <- as.integer(nr.covs)
  stopifnot(is.numeric(mu))
  nr.groups <- length(mu)
  if(!is.null(n) && (!is.numeric(n) || length(n) != nr.groups || any(n <= 0))){
    stop("'n' must be a positive (integer) number.")
  }
  if(!is.null(n) && !is.integer(n)) n <- as.integer(n)
  if(is.null(group.ratio)){ 
    ## balanced design
    group.ratio <- rep(1, nr.groups)
  }
  if(!is.numeric(group.ratio) || length(group.ratio) != nr.groups)
    stop("Length of 'group.ratio' must be identical to length of mu.")
  if(is.null(contr.mat)){
    contr.mat <- cbind(rep(1,nr.groups-1), -diag(nr.groups-1))
  }
  if(!is.numeric(contr.mat) || !is.matrix(contr.mat)){
    stop("'contr.mat' must be a numeric matrix including the contrasts.")
  }
  if(nr.groups != ncol(contr.mat)){
    stop("The number of columns of 'contr.mat' must be identical to length of 'mu'.")
  }
  if(!is.numeric(sig.level) || length(sig.level) != 1L || any(0 >= sig.level || sig.level >= 1))
    stop("'sig.level' must be numeric in (0, 1)")
  if(!is.null(power) && (!is.numeric(power) || length(sig.level) != 1L || any(0 >= power | power >= 1)))
    stop("'power' must be numeric in (0, 1)")
  if(!is.numeric(rel.tol) || length(rel.tol) != 1L || any(rel.tol <= 0)){
    stop("'rel.tol' must be positive number.")
  }
  
  Cmat <- contr.mat %*% matrix(mu, nrow = nr.groups)
  df1 <- nrow(contr.mat)
  Qmat <- diag(sum(group.ratio)/group.ratio)
  gamma.sq <- as.vector(t(Cmat) %*% solve(contr.mat %*% Qmat %*% t(contr.mat)) %*% Cmat/var)
  emp.power.fun <- function(n, nr.groups, group.ratio, nr.covs, sig.level, df1, 
                            gamma.sq, rel.tol, power = NULL){
    NT <- sum(n*group.ratio)
    df2 <- NT - nr.groups - nr.covs
    F.crit <- qf(1-sig.level, df1 = df1, df2 = df2)
    pow.fun.t <- function(x, F.crit, NT, gamma.sq, df1, df2){
      ncp <- NT*gamma.sq/(1+1/(df2+1)*x^2)
      pf(F.crit, df1 = df1, df2 = df2, ncp = ncp, lower.tail = FALSE)*dt(x, df = df2+1)
    }
    pow.fun.beta <- function(x, F.crit, NT, gamma.sq, df1, df2, nr.covs){
      ncp <- NT*gamma.sq*x
      pf(F.crit, df1 = df1, df2 = df2, ncp = ncp, lower.tail = FALSE)*dbeta(x, shape1 = (df2+1)/2, shape2 = nr.covs/2)
    }
    if(nr.covs == 1){
      res <- integrate(pow.fun.t, lower = -Inf, upper = Inf, rel.tol = rel.tol,
                       F.crit = F.crit, NT = NT, gamma.sq = gamma.sq, df1 = df1, 
                       df2 = df2)$value
    }else{
      res <- integrate(pow.fun.beta, lower = 0, upper = 1, rel.tol = rel.tol,
                       F.crit = F.crit, NT = NT, gamma.sq = gamma.sq, df1 = df1, 
                       df2 = df2, nr.covs = nr.covs)$value
    }
    if(is.null(power)){
      return(res)
    }else{
      return(res - power)
    }
  }
  if(is.null(n)){
    N.start <- 2*nr.groups - 1
    N1 <- uniroot(f = emp.power.fun, lower = N.start, upper = n.max, 
                 tol=10*rel.tol, nr.groups = nr.groups, group.ratio = group.ratio,
                 nr.covs = nr.covs, sig.level = sig.level, df1 = df1,
                 gamma.sq = gamma.sq, rel.tol = rel.tol, power = power)$root
    Ns <- N1*group.ratio
    NOTE <- paste("Total sample size:", sum(ceiling(Ns)))
  }
  if(is.null(power)){
    N1 <- n[1]
    Ns <- n
    power <- emp.power.fun(n = N1, nr.groups = nr.groups, group.ratio = group.ratio,
                           nr.covs = nr.covs, sig.level = sig.level, df1 = df1,
                           gamma.sq = gamma.sq, rel.tol = rel.tol)
    NOTE <- paste("Total sample size:", sum(n))
  }
  
  METHOD <- "ANCOVA power calculation"
  structure(list(ns = Ns, mu = mu, var = var, nr.covs = nr.covs, 
                 sig.level = sig.level, power = power, note = NOTE, method = METHOD), 
            class = "power.htest")
}
