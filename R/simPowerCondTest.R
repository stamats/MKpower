sim.power.cond.test <- function(nx, rx = NULL, rx.H0 = NULL, ny, ry = NULL, ry.H0 = NULL, 
                                sig.level = 0.05, conf.int = FALSE, mu = 0, 
                                alternative = c("two.sided", "less", "greater"), 
                                iter = 10000, normality = TRUE, variance = TRUE, 
                                var.equal = FALSE, exact = TRUE){
  if(is.null(rx) && is.null(rx.H0) && is.null(ry) && is.null(ry.H0))
    stop("No functions for random number generation provided.")
  if((is.null(rx) && !is.null(ry)) || (!is.null(rx) && is.null(ry)))
    stop("'rx' or 'ry' is missing.")
  if((is.null(rx.H0) && !is.null(ry.H0)) || (!is.null(rx.H0) && is.null(ry.H0)))
    stop("'rx.H0' or 'ry.H0' is missing.")
  alternative <- match.arg(alternative)
  if(length(sig.level) != 1)
    stop("'sig.level' has to be of length 1 (type I error)")
  if(sig.level <= 0 | sig.level >= 1)
    stop("'sig.level' has to be in (0, 1)")
  stopifnot(is.numeric(nx), length(nx) == 1, nx >= 1)
  nx <- trunc(nx)
  stopifnot(is.numeric(ny), length(ny) == 1, ny >= 1)
  ny <- trunc(ny)
  stopifnot(is.numeric(mu), length(mu) == 1)
  stopifnot(is.numeric(iter), length(iter) == 1, iter >= 1)
  iter <- trunc(iter)
  if(!normality && !exact)
    stop("Pre-test for normality or pre-test for equal variance or both must be selected")
  
  if(conf.int){
    conf.level <- 1 - sig.level
    alpha <- sig.level
  }else{
    conf.level <- NA
  }
  
  SIM.POW <- ifelse(!is.null(rx), TRUE, FALSE)
  SIM.ALPHA <- ifelse(!is.null(rx.H0), TRUE, FALSE)
  
  ## generate the data
  if(SIM.POW){
    stopifnot(is.function(rx), is.function(ry))
    data.x <- matrix(rx(nx*iter), nrow = iter)
    data.y <- matrix(ry(ny*iter), nrow = iter)
  }
  if(SIM.ALPHA){
    stopifnot(is.function(rx.H0), is.function(ry.H0))
    data.x.H0 <- matrix(rx.H0(nx*iter), nrow = iter)
    data.y.H0 <- matrix(ry.H0(ny*iter), nrow = iter)
  }
  
  RESULT <- list()
  row.shapiro.test <- function(x){ shapiro.test(x)$p.value }
  
  if(SIM.POW){
    res <- data.frame(statistic = numeric(iter),
                      pvalue = numeric(iter),
                      normality = NA,
                      variance = NA)
    ## pre-testing for normality
    if(normality){
      norm.pval1 <- apply(data.x, 1, row.shapiro.test)
      norm.pval2 <- apply(data.y, 1, row.shapiro.test)
      ind.normality <- norm.pval1 >= 0.05 & norm.pval2 >= 0.05
      res.wilcox <- row_wilcoxon_twosample(data.x[!ind.normality,], 
                                           data.y[!ind.normality,], alternative, 
                                           null = mu, exact = exact)
      res$statistic[!ind.normality] <- res.wilcox$statistic
      res$pvalue[!ind.normality] <- res.wilcox$pvalue
    }else{
      ind.normality <- !logical(iter)
    }
    res$normality <- ind.normality
    
    ## pre-testing for equal variance
    if(variance){
      data.xy <- cbind(data.x[ind.normality,], data.y[ind.normality,])
      group <- factor(c(rep("1", nx), rep("2", ny)))
      res.levene <- row_levene(data.xy, group)$pvalue
      ind.var.equal <- res.levene >= 0.05
      res$variance[ind.normality] <- ind.var.equal
    }

    if(variance){
      res.student <- row_t_equalvar(data.x[ind.normality,][ind.var.equal,], 
                                    data.y[ind.normality,][ind.var.equal,], 
                                    alternative = alternative, null = mu, 
                                    conf.level = conf.level)
      res$statistic[ind.normality][ind.var.equal] <- res.student$statistic
      res$pvalue[ind.normality][ind.var.equal] <- res.student$pvalue
      res.welch <- row_t_welch(data.x[ind.normality,][!ind.var.equal,], 
                               data.y[ind.normality,][!ind.var.equal,], 
                               alternative = alternative, null = mu, 
                               conf.level = conf.level)
      res$statistic[ind.normality][!ind.var.equal] <- res.welch$statistic
      res$pvalue[ind.normality][!ind.var.equal] <- res.welch$pvalue
    }else if (var.equal){
      res.student <- row_t_equalvar(data.x[ind.normality,], 
                                    data.y[ind.normality,], 
                                    alternative = alternative, null = mu, 
                                    conf.level = conf.level)
      res$statistic[ind.normality] <- res.student$statistic
      res$pvalue[ind.normality] <- res.student$pvalue
      res$variance[ind.normality] <- TRUE
    }else{
      res.welch <- row_t_welch(data.x[ind.normality,], 
                               data.y[ind.normality,], 
                               alternative = alternative, null = mu, 
                               conf.level = conf.level)
      res$statistic[ind.normality] <- res.welch$statistic
      res$pvalue[ind.normality] <- res.welch$pvalue
      res$variance[ind.normality] <- FALSE
    }
  }
  if(SIM.ALPHA){
    res.H0 <- data.frame(statistic = numeric(iter),
                         pvalue = numeric(iter),
                         normality = NA,
                         variance = NA)
    ## pre-testing for normality
    if(normality){
      norm.pval1 <- apply(data.x.H0, 1, row.shapiro.test)
      norm.pval2 <- apply(data.y.H0, 1, row.shapiro.test)
      ind.normality <- norm.pval1 >= 0.05 & norm.pval2 >= 0.05
      res.wilcox <- row_wilcoxon_twosample(data.x.H0[!ind.normality,], 
                                           data.y.H0[!ind.normality,], alternative, 
                                           null = mu, exact = exact)
      res.H0$statistic[!ind.normality] <- res.wilcox$statistic
      res.H0$pvalue[!ind.normality] <- res.wilcox$pvalue
    }else{
      ind.normality <- !logical(iter)
    }
    res.H0$normality <- ind.normality
    
    ## pre-testing for equal variance
    if(variance){
      data.xy <- cbind(data.x.H0[ind.normality,], data.y.H0[ind.normality,])
      group <- factor(c(rep("1", nx), rep("2", ny)))
      res.levene <- row_levene(data.xy, group)$pvalue
      ind.var.equal <- res.levene >= 0.05
      res.H0$variance[ind.normality] <- ind.var.equal
    }
    
    if(variance){
      res.student <- row_t_equalvar(data.x.H0[ind.normality,][ind.var.equal,], 
                                    data.y.H0[ind.normality,][ind.var.equal,], 
                                    alternative = alternative, null = mu, 
                                    conf.level = conf.level)
      res.H0$statistic[ind.normality][ind.var.equal] <- res.student$statistic
      res.H0$pvalue[ind.normality][ind.var.equal] <- res.student$pvalue
      res.welch <- row_t_welch(data.x.H0[ind.normality,][!ind.var.equal,], 
                               data.y.H0[ind.normality,][!ind.var.equal,], 
                               alternative = alternative, null = mu, 
                               conf.level = conf.level)
      res.H0$statistic[ind.normality][!ind.var.equal] <- res.welch$statistic
      res.H0$pvalue[ind.normality][!ind.var.equal] <- res.welch$pvalue
    }else if (var.equal){
      res.student <- row_t_equalvar(data.x.H0[ind.normality,], 
                                    data.y.H0[ind.normality,], 
                                    alternative = alternative, null = mu, 
                                    conf.level = conf.level)
      res.H0$statistic[ind.normality] <- res.student$statistic
      res.H0$pvalue[ind.normality] <- res.student$pvalue
      res.H0$variance[ind.normality] <- TRUE
    }else{
      res.welch <- row_t_welch(data.x.H0[ind.normality,], 
                               data.y.H0[ind.normality,], 
                               alternative = alternative, null = mu, 
                               conf.level = conf.level)
      res.H0$statistic[ind.normality] <- res.welch$statistic
      res.H0$pvalue[ind.normality] <- res.welch$pvalue
      res.H0$variance[ind.normality] <- FALSE
    }
  }
  if(SIM.POW && SIM.ALPHA) RESULT <- list("H1" = res, "H0" = res.H0)
  if(SIM.POW && !SIM.ALPHA) RESULT <- list("H1" = res)
  if(!SIM.POW && SIM.ALPHA) RESULT <- list("H0" = res.H0)
  
  SetUp <- c("nx" = nx, "rx" = rx, "rx.H0" = rx.H0, 
             "ny" = ny, "ry" = ry, "ry.H0" = ry.H0, 
             "sig.level" = sig.level, "mu" = mu,
             "alternative" = alternative, "iter" = iter)
  RESULT <- c(RESULT, list(SetUp = SetUp))
  class(RESULT) <- "sim.power.cond.test"
  RESULT
}

print.sim.power.cond.test <- function(x, digits = getOption("digits"), ...){
  cat("\n    Simulation Set-up\n")
  alpha <- x$SetUp$sig.level
  y0 <- x$SetUp
  iter <- x$SetUp$iter
  cat(paste(format(names(y0), width = 15L, justify = "right"), 
            format(y0, digits = digits), sep = " = "), sep = "\n")
  RES <- NULL
  cat("\n    Conditional Two-sample Test\n")
  y1 <- NULL
  if(!is.null(x$H1)){
    y1 <- c("emp.power" = sum(x$H1$pvalue < alpha)/iter)
  }
  if(!is.null(x$H0)){
    y1 <- c(y1, "emp.type.I.error" = sum(x$H0$pvalue < alpha)/iter)
  }
  cat(paste(format(names(y1), width = 15L, justify = "right"), 
            format(y1, digits = digits), sep = " = "), sep = "\n")
  RES <- c(RES, y1)
  cat("\n")
  invisible(RES)
}
