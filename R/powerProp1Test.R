power.prop1.test <- function(n = NULL, p1 = NULL, p0 = 0.5, sig.level = 0.05, 
                             power = NULL, alternative = c("two.sided", "less", "greater"),
                             cont.corr = TRUE, tol = .Machine$double.eps^0.25){
  if (sum(sapply(list(n, p1, power, sig.level), is.null)) != 1)
    stop("exactly one of 'n', 'p1', 'p0', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop("'sig.level' must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1))
    stop("'power' must be numeric in [0, 1]")
  if (!is.null(p1) && !is.numeric(p1) || any(0 > p1 | p1 > 1))
    stop("'p1' must be numeric in [0, 1]")
  if (!is.null(p0) && !is.numeric(p0) || any(0 > p0 | p0 > 1))
    stop("'p0' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)
  tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
  if(tside == 1) stopifnot(p1 < p0)
  if(tside == 3) stopifnot(p1 > p0)
  if(tside == 2){
    if(!is.null(sig.level)) z.alpha <- qnorm(1-sig.level/2)
  }else{
    if(!is.null(sig.level)) z.alpha <- qnorm(1-sig.level)
  }
  
  if(!is.null(power)) z.beta <- qnorm(power)
  if(is.null(n)){
    n <- ((z.alpha*sqrt(p0*(1-p0)) + z.beta*sqrt(p1*(1-p1)))/(p1-p0))^2
    if(cont.corr) n <- 0.25*n*(1 + sqrt(1 + 2/(n*abs(p1-p0))))^2
  }
  if(is.null(power))
    if(cont.corr)
      power <- pnorm((sqrt(n)*abs(p1-p0) - 0.5/sqrt(n) - 
                        z.alpha*sqrt(p0*(1-p0)))/sqrt(p1*(1-p1)))
    else
      power <- pnorm((sqrt(n)*abs(p1-p0) - 
                        z.alpha*sqrt(p0*(1-p0)))/sqrt(p1*(1-p1)))
  if(is.null(sig.level)){
    if(cont.corr)
      sig.level <- pnorm((sqrt(n)*abs(p1-p0) - 0.5/sqrt(n) - 
                            z.beta*sqrt(p1*(1-p1)))/sqrt(p0*(1-p0)),
                         lower.tail = FALSE)
    else
      sig.level <- pnorm((sqrt(n)*abs(p1-p0) - 
                            z.beta*sqrt(p1*(1-p1)))/sqrt(p0*(1-p0)),
                         lower.tail = FALSE)
    if(tside == 2) sig.level <- 2*sig.level
  }
  if(is.null(p1)){
    if(cont.corr)
      p1.fun <- function(p1){
        sqrt(n)*abs(p1-p0) - 0.5/sqrt(n) - z.beta*sqrt(p1*(1-p1)) - 
          z.alpha*sqrt(p0*(1-p0))
      }
    else
      p1.fun <- function(p1){
        sqrt(n)*abs(p1-p0) - z.beta*sqrt(p1*(1-p1)) - z.alpha*sqrt(p0*(1-p0))
      }
    interval <- if(p0 <= 0.5) c(p0, 1) else c(0, p0)
    p1 <- uniroot(p1.fun, interval = interval, tol = tol)$root
    delta <- abs(p1 - p0)
    if(tside == 1) p1 <- p0 - delta
    if(tside == 3) p1 <- p0 + delta
  }else{
    delta <- abs(p1-p0)
  }
  
  if(tside == 2){
    z.alpha <- qnorm(1-sig.level/2)
  }else{
    z.alpha <- qnorm(1-sig.level)
  }
  if(tside == 2){
    if(cont.corr) 
      crit.val <- n*p0 - z.alpha*sqrt(n*p0*(1-p0)) - 0.5
    else
      crit.val <- n*p0 - z.alpha*sqrt(n*p0*(1-p0))
    if(p1 > p0)
      exact.power <- pbinom(ceiling(n) - crit.val, size = ceiling(n), 
                            prob = p1, lower.tail = FALSE)
    else
      exact.power <- pbinom(crit.val, size = ceiling(n), prob = p1)
    exact.sig.level <- 2*pbinom(crit.val, size = ceiling(n), prob = p0)
  }
  if(tside == 1){
    if(cont.corr) 
      crit.val <- n*p0 - z.alpha*sqrt(n*p0*(1-p0)) - 0.5
    else
      crit.val <- n*p0 - z.alpha*sqrt(n*p0*(1-p0))
    exact.power <- pbinom(crit.val, size = ceiling(n), prob = p1) 
    exact.sig.level <- pbinom(crit.val, size = ceiling(n), prob = p0)
  }
  if(tside == 3){
    if(cont.corr) 
      crit.val <- n*p0 + z.alpha*sqrt(n*p0*(1-p0)) + 0.5
    else
      crit.val <- n*p0 + z.alpha*sqrt(n*p0*(1-p0))
    exact.power <- pbinom(crit.val+1, size = ceiling(n), prob = p1, 
                          lower.tail = FALSE)
    exact.sig.level <- pbinom(crit.val+1, size = ceiling(n), prob = p0, 
                              lower.tail = FALSE)
  }

  if(cont.corr)
    METHOD <- "Power calculation for testing a given proportion (with continuity correction)"
  else
    METHOD <- "Power calculation for testing a given proportion"
  NOTE <- "n = total sample size"
  structure(list(n = n, delta = delta, p1 = p1, p0 = p0, sig.level = sig.level, 
                 exact.sig.level = exact.sig.level,
                 power = power, exact.power = exact.power, alternative = alternative,
                 note = NOTE, method = METHOD), class = "power.htest")
}
