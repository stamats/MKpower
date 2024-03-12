ssize.reference.range <- function (n = NULL, delta = NULL, ref.prob = 0.95, conf.prob = NULL,
                                   alternative = c("two.sided", "one.sided"),
                                   method = "parametric", exact = TRUE, 
                                   tol = .Machine$double.eps^0.5) {
  if(sum(sapply(list(n, delta, ref.prob, conf.prob), is.null)) != 1)
    stop("exactly one of 'n', 'delta', 'ref.prob', and 'conf.prob' must be NULL")
  if (!is.null(ref.prob) && !is.numeric(ref.prob) || any(0 > ref.prob | ref.prob > 1))
    stop("'ref.prob' must be numeric in [0, 1]")
  if (!is.null(conf.prob) && !is.numeric(conf.prob) || any(0 > conf.prob | conf.prob > 1))
    stop("'conf.prob' must be numeric in [0, 1]")
  
  alternative <- match.arg(alternative)
  if(!is.null(ref.prob)){
    if(alternative == "two.sided"){
      qstar <- (1-ref.prob)/2
    } else {
      qstar <- 1-ref.prob
    }    
  }
  if(!is.null(delta)){
    if(qstar - delta < 0){
      stop("(1-ref.prob)/2-delta (two-sided) resp. (1-ref.prob)-delta (one-sided) must be > 0")
    }
  }
  
  if(method == "parametric"){
    if(exact){
      p.body <- if (alternative == "two.sided")
        quote({
          nu <- n - 1
          qstar <- (1-ref.prob)/2
          z <- qnorm(1-qstar)
          pt(sqrt(n)*z, nu, ncp = -sqrt(n)*qnorm(qstar+delta)) - 
            pt(sqrt(n)*z, nu, ncp = -sqrt(n)*qnorm(qstar-delta))
        })
      else 
       quote({
        nu <- n - 1
        qstar <- 1-ref.prob
        z <- qnorm(1-qstar)
        pt(sqrt(n)*z, nu, ncp = -sqrt(n)*qnorm(qstar+delta)) - 
          pt(sqrt(n)*z, nu, ncp = -sqrt(n)*qnorm(qstar-delta))
       })
      NOTE <- "Exact method for normal distribution"
    }else{
      p.body <- if (alternative == "two.sided")
        quote({
          qstar <- (1-ref.prob)/2
          z <- qnorm(1-qstar)
          2*pnorm(delta*sqrt(n/(1+0.5*z^2))/dnorm(z)) - 1
        })
      else
        quote({
          qstar <- 1-ref.prob
          z <- qnorm(1-qstar)
          2*pnorm(delta*sqrt(n/(1+0.5*z^2))/dnorm(z)) - 1
        })
      NOTE <- "Approximate method for normal distribution"
    }
  }
  if(method == "nonparametric"){
    if(exact){
      p.body <- if (alternative == "two.sided")
        quote({
          qstar <- (1-ref.prob)/2
          k <- round((n+1)*qstar, digits = 0)
          pbeta(qstar + delta, shape1 = k, shape2 = n-k+1) - 
            pbeta(qstar-delta, shape1 = k, shape2 = n-k+1)
        })
      else
        quote({
          qstar <- 1-ref.prob
          k <- round((n+1)*qstar, digits = 0)
          pbeta(qstar + delta, shape1 = k, shape2 = n-k+1) - 
            pbeta(qstar-delta, shape1 = k, shape2 = n-k+1)
        })
      NOTE <- "Exact non-parametric method"
    }else{
      p.body <- if (alternative == "two.sided")
        quote({
          qstar <- (1-ref.prob)/2
          z <- qnorm(1-qstar)
          2*pnorm(delta*sqrt(n/(qstar*(1-qstar)))) - 1
        })
      else
        quote({
          qstar <- 1-ref.prob
          z <- qnorm(1-qstar)
          2*pnorm(delta*sqrt(n/(qstar*(1-qstar)))) - 1
        })
      NOTE <- "Approximate non-parametric method"
    }
  }
  
  
  if (is.null(conf.prob))
    conf.prob <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - conf.prob, c(2, 1e+07),
                 tol = tol, extendInt = "upX")$root
  else if (is.null(delta)){
    if(alternative == "two.sided"){
      qstar <- (1-ref.prob)/2
    } else {
      qstar <- 1-ref.prob
    }
    delta <- uniroot(function(delta) eval(p.body) - conf.prob,
                     c(1e-10, qstar-1e-10), tol = tol, extendInt = "yes")$root
  }
  else if (is.null(ref.prob))
    sig.level <- uniroot(function(ref.prob) eval(p.body) - conf.prob, 
                         c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
  else stop("internal error", domain = NA)
  
  METHOD <- "Reference range sample size calculation"
  structure(list(n = n, delta = delta, ref.prob = ref.prob,
                 conf.prob = conf.prob, alternative = alternative, note = NOTE,
                 method = METHOD), class = "power.htest")
}
