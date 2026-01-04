sim.power.xiao.t.test <- function(n1, n2, delta = 1, sd1 = 1, sd2 = 2, mu = 0, sig.level = 0.05, 
                                  alternative = c("two.sided", "less", "greater"),
                                  iter = 10000, parallel = FALSE, cl = NULL){
  stopifnot(is.numeric(n1), length(n1) == 1, n1 >= 1)
  n1 <- trunc(n1)
  stopifnot(is.numeric(n2), length(n2) == 1, n2 >= 1)
  n2 <- trunc(n2)
  
  if (!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop("'sig.level' must be numeric in [0, 1]")
  if (!is.numeric(sd1) || any(0 > sd1))
    stop("'sd1' must be a positive numeric")
  if (!is.numeric(sd2) || any(0 > sd2))
    stop("'sd2' must be a positive numeric")
  
  stopifnot(is.numeric(mu), length(mu) == 1)
  stopifnot(is.numeric(iter), length(iter) == 1, iter >= 1)
  
  alternative <- match.arg(alternative)
  iter <- trunc(iter)

  if(parallel){
    if(!("package:MKinfer" %in% search())){
      DETACH <- TRUE
      attachNamespace("MKinfer")
    }else{
      DETACH <- FALSE
    }
    if(is.null(cl)){
      STOP <- TRUE
      ncores <- detectCores()-1
      cl <- makePSOCKcluster(rep("localhost", ncores))
    }else{
      STOP <- FALSE
    }
  }

  rx <- function(n) rnorm(n, mean = 0, sd = sd1)
  ry <- function(n) rnorm(n, mean = delta, sd = sd2)
  
  data.x <- matrix(rx(n1*iter), nrow = iter)
  data.y <- matrix(ry(n2*iter), nrow = iter)
  res <- row_t_welch(data.x, data.y, 
                     alternative = alternative, null = mu, 
                     conf.level = 1-sig.level)
  vx <- res$var.x
  vy <- res$var.y
  df <- ifelse(n1*(n1-1)/vx <= n2*(n2-1)/vy, 
               (n2-1)*(1 + n2/n1*vx/vy), (n1-1)*(1 + n1/n2*vy/vx))
  res$df <- df
  if(alternative == "less"){
    res$pvalue <- pgt(res$statistic, res$statistic, n1 = n1, n2 = n2, v1tov2 = vx/vy, 
                      parallel = parallel, cl = cl)
  }
  if(alternative == "greater"){
    res$pvalue <- pgt(res$statistic, n1 = n1, n2 = n2, v1tov2 = vx/vy, 
                      lower.tail = FALSE, parallel = parallel, 
                      cl = cl)
  }
  if(alternative == "two.sided"){
    res$pvalue <- 2 * pgt(-abs(res$statistic), n1 = n1, n2 = n2, 
                          v1tov2 = vx/vy, parallel = parallel, cl = cl)
  }
  empPower <- sum(res$pvalue < sig.level)/iter

  if(parallel){
    if(STOP) stopCluster(cl)
    if(DETACH) detach("package:MKinfer")
  }
  
  METHOD <- "Xiao two-sample t-test"
  NOTE <- "minimum total sample size (sum of group sample sizes)"
  res <- structure(list(n1 = n1, n2 = n2, delta = delta, sd1 = sd1, sd2 = sd2, 
                        sig.level = sig.level, emp.power = empPower, 
                        alternative = alternative, note = NOTE, 
                        method = METHOD), class = "power.htest")
  res
}

