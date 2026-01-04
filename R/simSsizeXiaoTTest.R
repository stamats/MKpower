sim.ssize.xiao.t.test <- function(delta = 1, sd1 = 1, sd2 = 2, mu = 0, sig.level = 0.05, 
                                  power = 0.8, alternative = c("two.sided", "one.sided"),
                                  n.delta = 5, step.size = 1, iter = 10000,
                                  parallel = FALSE, cl = NULL){
  alternative <- match.arg(alternative)
  if (!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop("'sig.level' must be numeric in [0, 1]")
  if (!is.numeric(power) || any(0 > power | power > 1))
    stop("'power' must be numeric in [0, 1]")
  if (!is.numeric(sd1) || any(0 > sd1))
    stop("'sd1' must be a positive numeric")
  if (!is.numeric(sd2) || any(0 > sd2))
    stop("'sd2' must be a positive numeric")
  
  stopifnot(is.numeric(mu), length(mu) == 1)
  stopifnot(is.numeric(n.delta), length(n.delta) == 1, n.delta >= 1)
  stopifnot(is.numeric(step.size), length(step.size) == 1)
  stopifnot(is.numeric(iter), length(iter) == 1, iter >= 1)
  
  n.delta <- trunc(n.delta)
  step.size <- trunc(step.size)
  iter <- trunc(iter)
  
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 2 && !is.null(delta)) delta <- abs(delta)
  
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
  
  res.welch <- ssize.welch.t.test(delta = delta, sd1 = sd1, sd2 = sd2, 
                                  sig.level = sig.level, power = power, 
                                  alternative = alternative)
  N1 <- seq(from = res.welch$n1 - n.delta, to = res.welch$n1+n.delta, 
            by = step.size)
  N2 <- seq(from = res.welch$n2 - n.delta, to = res.welch$n2+n.delta, 
            by = step.size)
  
  empPower <- matrix(0, nrow = length(N1), ncol = length(N2))
  for(i in seq_len(length(N1))){
    nx <- N1[i]
    for(j in seq_len(length(N2))){
      ny <- N2[j]
      data.x <- matrix(rx(nx*iter), nrow = iter)
      data.y <- matrix(ry(ny*iter), nrow = iter)
      res <- row_t_welch(data.x, data.y, 
                         alternative = alternative, null = mu, 
                         conf.level = 1-sig.level)
      vx <- res$var.x
      vy <- res$var.y
      df <- ifelse(nx*(nx-1)/vx <= ny*(ny-1)/vy, 
                   (ny-1)*(1 + ny/nx*vx/vy), (nx-1)*(1 + nx/ny*vy/vx))
      res$df <- df
      if(alternative == "one.sided"){
        tmp <- try(pgt(res$statistic, n1 = nx, n2 = ny, v1tov2 = vx/vy, 
                       lower.tail = FALSE, parallel = parallel, cl = cl))
        if(inherits(tmp, "try-error")){
          res$pvalue <- NA
        }else{
          res$pvalue <- tmp
        }
      }
      if(alternative == "two.sided"){
        tmp <- try(2 * pgt(-abs(res$statistic), n1 = nx, n2 = ny, 
                           v1tov2 = vx/vy, parallel = parallel, cl = cl))
        if(inherits(tmp, "try-error")){
          res$pvalue <- NA
        }else{
          res$pvalue <- tmp
        }
      }
      empPower[i,j] <- sum(res$pvalue < sig.level)/iter
    }
  }
  empPower[empPower == 0] <- NA
  SUM <- outer(N1, N2, "+")
  SSIZE <- outer(N1, N2, function(n1, n2) paste0(n1, ",", n2))
  ind <- which.min(SUM[empPower > power])
  ssize <- SSIZE[empPower > power][ind]
  POW <- empPower[empPower > power][ind]
  n <- as.integer(unlist(strsplit(ssize, "\\,")))
  
  if(parallel){
    if(STOP) stopCluster(cl)
    if(DETACH) detach("package:MKinfer")
  }
  
  METHOD <- "Xiao two-sample t-test"
  NOTE <- "minimum total sample size (sum of group sample sizes)"
  res <- structure(list(n1 = n[1], n2 = n[2], delta = delta, sd1 = sd1, sd2 = sd2, 
                        sig.level = sig.level, emp.power = POW, 
                        power = power, 
                        alternative = alternative, note = NOTE, 
                        method = METHOD), class = "power.htest")
  res
}

