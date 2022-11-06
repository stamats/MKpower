sim.power.wilcox.test <- function(nx, rx, rx.H0 = NULL, ny, ry, ry.H0 = NULL, 
                                  alternative = c("two.sided", "less", "greater"), 
                                  sig.level = 0.05, conf.int = FALSE, approximate = FALSE,
                                  ties = FALSE, iter = 10000, nresample = 10000,
                                  parallel = "no", ncpus = 1L, cl = NULL){
  alternative <- match.arg(alternative)
  if(length(sig.level) != 1)
    stop("'sig.level' has to be of length 1 (type I error)")
  if(sig.level <= 0 | sig.level >= 1)
    stop("'sig.level' has to be in (0, 1)")
  stopifnot(is.function(rx), is.function(ry))
  stopifnot(is.numeric(nx), length(nx) == 1, nx >= 1)
  nx <- trunc(nx)
  stopifnot(is.numeric(ny), length(ny) == 1, ny >= 1)
  ny <- trunc(ny)
  stopifnot(is.numeric(iter), length(iter) == 1, iter >= 1)
  iter <- trunc(iter)
  nresample <- trunc(nresample)
  
  conf.level <- 1 - sig.level
  alpha <- sig.level
  data.x <- matrix(rx(nx*iter), nrow = iter)
  data.y <- matrix(ry(ny*iter), nrow = iter)
  if(!is.null(rx.H0) & !is.null(ry.H0)){
    data.x.H0 <- matrix(rx.H0(nx*iter), nrow = iter)
    data.y.H0 <- matrix(ry.H0(ny*iter), nrow = iter)
  }
  wfun <- function(d, nx, ny, alternative, conf.level, distribution, conf.int = TRUE){
    group <- factor(c(rep(1, nx), rep(2, ny)))
    DF <- data.frame(x = d, group = group)
    temp <- wilcox_test(x ~ group, data = DF, distribution = distribution, 
                        conf.int = conf.int, alternative = alternative, 
                        conf.level = conf.level)
    if(conf.int){
      loc.diff <- confint(temp)$estimate
      names(loc.diff) <- NULL
      confInt <- confint(temp)$conf.int
    }else{
      loc.diff <- numeric(1)
      confInt <- numeric(2)
    }
    loc.null <- temp@nullvalue
    res <- c("obs.x" = nx, "obs.y" = ny, "obs.tot" = nx+ny, 
             "loc.diff" = ifelse(conf.int, loc.diff, FALSE), 
             "statistic" = statistic(temp), "pvalue" = pvalue(temp)[1], 
             "conf.low" = ifelse(conf.int, confInt[1], FALSE),
             "conf.high" = ifelse(conf.int, confInt[2], FALSE), 
             "loc.null" = loc.null, "conf.level" = conf.level)
    res
  }
  ## exact
  data.xy <- cbind(data.x, data.y)
  if(conf.int){
    res <- data.frame(t(apply(data.xy, 1, wfun, nx, ny, alternative, 
                              conf.level, "exact")))
  }else{
    if(ties){
      res <- data.frame(t(apply(data.xy, 1, wfun, nx, ny, alternative, 
                                conf.level, "exact", conf.int = FALSE)))
    }else{
      res <- row_wilcoxon_twosample(data.x, data.y, alternative, mu = 0, 
                                    exact = TRUE)
    }
  }
  if(!is.null(rx.H0) & !is.null(ry.H0)){
    data.xy.H0 <- cbind(data.x.H0, data.y.H0)
    if(conf.int){
      res.H0 <- data.frame(t(apply(data.xy.H0, 1, wfun, nx, ny, alternative, 
                                   conf.level, "exact")))
    }else{
      if(ties){
        res.H0 <- data.frame(t(apply(data.xy.H0, 1, wfun, nx, ny, alternative, 
                                     conf.level, "exact", conf.int = FALSE)))
      }else{
        res.H0 <- row_wilcoxon_twosample(data.x.H0, data.y.H0, alternative, 
                                         mu = 0, exact = TRUE)
      }
      
    }
    EXACT <- list("H1" = res, "H0" = res.H0)
  }else{
    EXACT <- list("H1" = res)  
  }
  ## asymptotic
  if(conf.int){
    res <- data.frame(t(apply(data.xy, 1, wfun, nx, ny, alternative, 
                              conf.level, "asymptotic")))
  }else{
    res <- row_wilcoxon_twosample(data.x, data.y, alternative, mu = 0, 
                                  exact = FALSE)
  }
  if(!is.null(rx.H0) & !is.null(ry.H0)){
    if(conf.int){
      res.H0 <- data.frame(t(apply(data.xy.H0, 1, wfun, nx, ny, alternative, 
                                   conf.level, "asymptotic")))
    }else{
      res.H0 <- row_wilcoxon_twosample(data.x.H0, data.y.H0, alternative, 
                                       mu = 0, exact = FALSE)
    }
    ASYMPTOTIC <- list("H1" = res, "H0" = res.H0)
  }else{
    ASYMPTOTIC <- list("H1" = res)  
  }
  ## approximate
  if(approximate){
    if(conf.int){
      res <- data.frame(t(apply(data.xy, 1, wfun, nx, ny, alternative, conf.level, 
                                approximate(nresample, parallel, ncpus, cl))))
    }else{
      res <- data.frame(t(apply(data.xy, 1, wfun, nx, ny, alternative, conf.level, 
                                approximate(nresample, parallel, ncpus, cl), 
                                conf.int = FALSE)))
    }
    if(!is.null(rx.H0) & !is.null(ry.H0)){
      if(conf.int){
        res.H0 <- data.frame(t(apply(data.xy.H0, 1, wfun, nx, ny, alternative, 
                                     conf.level, 
                                     approximate(nresample, parallel, ncpus, cl))))
      }else{
        res.H0 <- data.frame(t(apply(data.xy.H0, 1, wfun, nx, ny, alternative, 
                                     conf.level, 
                                     approximate(nresample, parallel, ncpus, cl),
                                     conf.int = FALSE)))
      }
      APPROXIMATE <- list("H1" = res, "H0" = res.H0)
    }else{
      APPROXIMATE <- list("H1" = res)  
    }
  }
  ## Set-up
  SetUp <- c("nx" = nx, "rx" = rx, "rx.H0" = rx.H0, 
             "ny" = ny, "ry" = ry, "ry.H0" = ry.H0, 
             "sig.level" = sig.level, "mu" = 0,
             "alternative" = alternative, "iter" = iter,
             "conf.int" = conf.int, "approximate" = approximate, 
             "ties" = ties)
  if(approximate){
    res <- list(Exact = EXACT, Asymptotic = ASYMPTOTIC, Approximate = APPROXIMATE,
                SetUp = SetUp)
  }else{
    res <- list(Exact = EXACT, Asymptotic = ASYMPTOTIC, Approximate = NULL,
                SetUp = SetUp)
  }
  class(res) <- "sim.power.wtest"
  res
}

print.sim.power.wtest <- function(x, digits = getOption("digits"), ...){
  cat("\n    Simulation Set-up\n")
  alpha <- x$SetUp$sig.level
  y0 <- x$SetUp
  iter <- x$SetUp$iter
  cat(paste(format(names(y0), width = 15L, justify = "right"), 
            format(y0, digits = digits), sep = " = "), sep = "\n")
  cat("\n    Exact Wilcoxon-Mann-Whitney Test\n")
  y1 <- c("emp.power" = sum(x$Exact$H1$pvalue < alpha)/iter)
  if(!is.null(x$Exact$H0)){
    y1 <- c(y1, "emp.type.I.error" = sum(x$Exact$H0$pvalue < alpha)/iter)
  }
  cat(paste(format(names(y1), width = 15L, justify = "right"), 
            format(y1, digits = digits), sep = " = "), sep = "\n")
  cat("\n    Asymptotic Wilcoxon-Mann-Whitney Test\n")
  y2 <- c("emp.power" = sum(x$Asymptotic$H1$pvalue < alpha)/iter)
  if(!is.null(x$Asymptotic$H0)){
    y2 <- c(y2, "emp.type.I.error" = sum(x$Asymptotic$H0$pvalue < alpha)/iter)
  }
  cat(paste(format(names(y2), width = 15L, justify = "right"), 
            format(y2, digits = digits), sep = " = "), sep = "\n")
  if(x$SetUp$approximate){
    cat("\n    Approximate Wilcoxon-Mann-Whitney Test\n")
    y3 <- c("emp.power" = sum(x$Approximate$H1$pvalue < alpha)/iter)
    if(!is.null(x$Approximate$H0)){
      y3 <- c(y3, "emp.type.I.error" = sum(x$Approximate$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y3), width = 15L, justify = "right"), 
              format(y3, digits = digits), sep = " = "), sep = "\n")
  }else{
    y3 <- NULL
  }
  cat("\n")
  invisible(c(y1, y2, y3))
}
