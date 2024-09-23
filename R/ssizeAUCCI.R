ssize.auc.ci <- function(AUC = NULL, delta = NULL, n = NULL, sig.level = 0.05,
                         power = NULL, prev = NULL, NMAX = 1e4){
  if(sum(sapply(list(n, delta, sig.level, power), is.null)) != 1) 
    stop("exactly one of 'n', 'delta', 'sig.level', and 'power' must be NULL")
  if(!is.null(AUC) && any(0.5 >= AUC | AUC > 1))
    stop("'AUC' must be numeric in (0.5, 1]")
  if(!is.null(delta) && !is.numeric(delta) || any(0 > delta | delta > 1))
    stop("'delta' must be numeric in [0, 1]")
  if(!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
    stop("'sig.level' must be numeric in [0, 1]")
  if(!is.null(power) && !is.numeric(power) || any(0 > power | power > 1))
    stop("'power' must be numeric in [0, 1]")
  if(!is.null(prev) && !is.numeric(prev) || any(0 > prev | prev > 1))
    stop("'prev' must be numeric in [0, 1]")

  if(is.null(prev)) prev <- 0.5

  if(!is.null(delta)){
    AUC.delta <- AUC - delta
    if(!is.null(AUC.delta) && !is.numeric(AUC.delta) || any(0 > AUC.delta | AUC.delta > 1))
      stop("'AUC'-'delta' must be numeric in [0, 1]")    
  }

  ## variance from Hanley and McNeil (1982)  
  Varfun <- function(N, AUC, prev){
    n <- prev*N
    m <- (1-prev)*N
    Q1 <- AUC/(2 - AUC)
    Q2 <- 2*AUC^2/(1+AUC)
    1/(n*m)*(AUC*(1-AUC) + n*(Q1-AUC^2) + m*(Q2-AUC^2) - (Q1 + Q2 - 2*AUC^2))
  }
  
  if(is.null(n)){
    zalpha <- qnorm(1-sig.level)
    zbeta <- qnorm(power)
    fun.N <- function(N, AUC, delta, zalpha, zbeta, prev){
      Varfun(N = N, AUC = AUC, prev = prev)*(zalpha + zbeta)^2 - delta^2
    }
    n <- uniroot(fun.N, lower = 1, upper = NMAX, AUC = AUC, delta = delta, 
                 zalpha = zalpha, zbeta = zbeta, prev = prev)$root
  }else if(is.null(delta)){
    zalpha <- qnorm(1-sig.level)
    zbeta <- qnorm(power)
    fun.delta <- function(delta, N, AUC, zalpha, zbeta, prev){
      Varfun(N = N, AUC = AUC, prev = prev)*(zalpha + zbeta)^2 - delta^2
    }
    delta <- uniroot(fun.delta, lower = 1e-10, upper = AUC-1e-10, N = n, 
                     AUC = AUC, zalpha = zalpha, zbeta = zbeta, prev = prev)$root
  }else if(is.null(power)){
    zalpha <- qnorm(1-sig.level)
    fun.power <- function(power, N, AUC, delta, zalpha, prev){
      zbeta <- qnorm(power)
      Varfun(N = N, AUC = AUC, prev = prev)*(zalpha + zbeta)^2 - delta^2
    }
    power <- uniroot(fun.power, lower = 1e-3, upper = 1-1e-10, N = n, AUC = AUC, 
                     delta = delta, zalpha = zalpha, prev = prev)$root
  }else if(is.null(sig.level)){
    zbeta <- qnorm(power)
    fun.sig.level <- function(sig.level, N, AUC, delta, zbeta, prev){
      zalpha <- qnorm(1-sig.level)
      Varfun(N = N, AUC = AUC, prev = prev)*(zalpha + zbeta)^2 - delta^2
    }
    sig.level <- uniroot(fun.sig.level, lower = 1e-3, upper = 1-1e-3, N = n, 
                         AUC = AUC, delta = delta, zbeta = zbeta, prev = prev)$root
  }else stop("internal error")
  
  METHOD <- "Sample size calculation for AUC confidence interval"
  NOTE <- "n is number of cases, n1 is number of controls"
  structure(list(AUC = AUC, n = prev*n, n1 = (1-prev)*n, delta = delta, 
                 sig.level = sig.level, power = power, prev = prev, 
                 note = NOTE, method = METHOD), 
            class = "power.htest")  
}
